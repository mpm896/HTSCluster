"""
Clustering methods for chemical compounds from high throughput screening data

@author: Matthew Martinez
https://github.com/mpm896/HTSCluster

1. Need to check that _calc_silhouette method works
"""

from typing import Callable, Dict, List, Optional

import numpy as np
from rdkit import Chem  # type: ignore
from rdkit.Chem import Mol  # type: ignore
from rdkit.Chem import AllChem # type: ignore # type: ignore
from rdkit.DataStructs import BulkTanimotoSimilarity, ExplicitBitVect  # type: ignore
from rdkit.ML.Cluster import Butina  # type: ignore
import polars as pl
from polars import DataFrame
from sklearn.cluster import KMeans  # type: ignore
from sklearn.metrics import silhouette_score  # type: ignore
from sklearn.decomposition import PCA  # type: ignore
from tqdm.auto import tqdm

from .utils.utils import get_fps


class ChemicalCluster:
    
    def __init__(
            self,
            fptype: str = "rdkit",
            cluster_type: str = "Butina",
            reduce: bool = False,
            n_clusters: Optional[int] = None,
            do_sillhoute: Optional[bool] = False,
    ) -> None:
        self.fptype: str = fptype
        self.cluster_type: str = cluster_type
        self.reduce = reduce
        if cluster_type == "KMeans":
            self.n_clusters: Optional[int] = n_clusters if n_clusters is not None and n_clusters > 0 else None
        self.do_sillhoute: Optional[bool] = do_sillhoute


    def cluster_smiles(
        self, 
        df: pl.DataFrame, 
        sim_cutoff: float = 0.8, 
        column_name: str = "SMILES"
    ) -> List[int]:
        """
        Cluster compounds based on the SMILES strings

        :param df: DataFrame object containing your compounds
        :param  sim_cutoff: Similarity cuttoff metric ONLY FOR Butina clustering
        :param column_name: Name of column to cluster. Default = "SMILES"

        :returns: list containing the cluster number of each compound
        """
        assert column_name in df.columns

        mols: List[Mol] = [Chem.MolFromSmiles(s) for s in df[column_name]]
        if self.cluster_type == "Butina":
            self.clusters = self._cluster_mols(mols, sim_cutoff=sim_cutoff)
        else:
            self.clusters = self._cluster_mols(mols)
        return self.clusters
    

    def _cluster_mols(self, mols: List[Mol], *, sim_cutoff: float=0.8) -> List[int]:
        """
        Cluster Mol objects

        :param mols: list of Mol objects
        :param sim_cutoff: Similarity cuttoff metric ONLY FOR Butina clustering

        :returns: list containing the cluster number for each compound


        *** FACTOR IN PCA REDUCTION ***
        """
        self.fp_list = get_fps(self.fptype, mols)
        if self.reduce:
            self.reduced_fps = self._reduce_hits(self.fp_list)
            fps = self.reduced_fps
        else:
            fps = self.fp_list

        cluster_method: Dict[str, Callable] = {
            "Butina": self._cluster_butina,
            "KMeans": self._cluster_kmeans
        }
        if cluster_method[self.cluster_type] is None:
            raise ValueError(f"Clustering method {self.cluster_type} is not supported.")

        if self.cluster_type == "KMeans":
            clustered = cluster_method[self.cluster_type](fps)
        else:
            clustered = cluster_method[self.cluster_type](fps, sim_cutoff)
        return clustered
    

    def get_mols(self, df: DataFrame, column_name: str="SMILES") -> List[Mol]:
        assert column_name in df.columns
        return [Chem.MolFromSmiles(s) for s in df[column_name]]
    
    
    def get_fps(self, mols: List[Mol]) -> List[ExplicitBitVect]:
        """
        Get the chemical bit fingerprint from Mol objects for each compound

        :param mols: a list of Mol objects

        :returns: a list of fingerprints (ExplicitBitVects) for each compound
        """
        fp_dict = {
            "rdkit": [Chem.RDKFingerprint(x) for x in mols],
            "morgan": [AllChem.GetMorganFingerprintAsBitVect(x, 2) for x in mols]
        }

        if fp_dict[self.fptype] is None:
            raise ValueError(f"Fingerprint method {self.fptype} is not supported.")
        
        return fp_dict[self.fptype]
    
    
    def _cluster_butina(
            self, 
            fps: list[ExplicitBitVect], 
            sim_cutoff: float
    ) -> np.ndarray[int]:
        """
        Cluster compounds using the Butina clustering algorithm

        :param fps: a list of fingerprints (ExplicitBitVects) for each compound
        :param sim_cutoff: Similarity cuttoff metric Butina clustering

        :returns: a list of cluster numbers for each compound
        """
        dist_cutoff = 1 - sim_cutoff
        dists = []
        nfps = len(self.fp_list)
        for i in tqdm(range(1, nfps), desc="Calculating Bulk Tanimoto Similarity"):
            sims = BulkTanimotoSimilarity(self.fp_list[i], self.fp_list[:i])
            dists.extend([1 - x for x in sims])

        self._dist_matrix = dists

        # Cluster
        mol_clusters = tqdm(Butina.ClusterData(
            dists, 
            nfps, 
            dist_cutoff, 
            isDistData=True
        ), desc="Clustering")
        cluster_ids = [0] * nfps
        for idx, cluster in enumerate(mol_clusters):
            for mol_idx in cluster:
                cluster_ids[mol_idx] = idx

        self.butina_clusters_ = np.array(
            [x - 1 for x in tqdm(cluster_ids, desc="Assigning clusters")]
        )
        return self.butina_clusters_
    

    def _cluster_kmeans(self, fps) -> np.ndarray[int]:
        """
        Cluster compounds using the Butina clustering algorithm

        :param fps: a list of fingerprints (ExplicitBitVects) for each compound

        :returns: a list of cluster numbers for each compound

        ***PLACEHOLD SELF.N_CLUSTERS FOR NOW. CHANGE TO BE 10% OF TOTAL COMPOUNDS)
        """
        if self.n_clusters is None:
            self.n_clusters = len(self.fp_list) // 10

        print(f"Size of FPs: {np.array(fps).shape}")

        k_means = KMeans(n_clusters=self.n_clusters, random_state=42, verbose=True)
        self.k_clusters_ = k_means.fit_predict(self.fp_list)
        return self.k_clusters_
    

    def _calc_sillhouette(self, X:  np.stack, min: int, max: int) -> DataFrame:
        """
        Calculate sillhoute scores for KMeans clustering

        :param X: numpy stack of fingerprints
        :param min: minimum number of clusters
        :param max: maximum number of clusters
        
        :returns: DataFrame object with cluster number and sillhoutte score

        *** CALCULATE MIN AS 5% OF TOTAL COMPOUNDS, MAX AS 20% OF TOTAL COMPOUNDS ***
        """
        clusters = range(min, max)
        scores: List = []
        for k in tqdm(clusters, desc="Calculating scores for {min} to {max} " \
                                     "clusters. This could take a while"):
            km = KMeans(n_clusters=k, n_init='auto', random_state=42)
            cluster_labels = km.fit_predict(X)
            score = silhouette_score(X, cluster_labels)
            scores.append([k, score])

        return pl.DataFrame(scores, schema={
            "K": pl.UInt8, 
            "Silhouette Score": pl.Float32
        })
    

    def _reduce_hits(
            self, 
            fps: List[ExplicitBitVect], 
            n_components: int | float=150
    ) -> List[ExplicitBitVect]:
        """
        Do a PCA reduction on the chemical fingerprints

        :param fps: Chemical bit fingerprints
        :param n_components: Number of components to retain in PCA reduction

        :returns: PCA reduced list of chemical bit fingerprints
        """
        reduction = PCA(n_components=n_components)
        reduced = reduction.fit_transform(fps)
        return reduced

        

        





    





