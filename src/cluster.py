from typing import Optional, List

import numpy as np
from rdkit import Chem  # type: ignore
from rdkit.Chem import Mol
from rdkit.Chem import AllChem, Draw, Descriptors  # type: ignore
from rdkit.Chem import rdFMCS # type: ignore
from rdkit.Chem.Scaffolds.MurckoScaffold import MakeScaffoldGeneric, MurckoScaffoldSmiles  # type: ignore
from rdkit.DataStructs import BulkTanimotoSimilarity, ExplicitBitVect  # type: ignore
from rdkit.ML.Cluster import Butina  # type: ignore
import polars as pl
from polars import DataFrame
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.cluster import KMeans  # type: ignore
from sklearn.decomposition import PCA  # type: ignore
from sklearn.neighbors import NearestNeighbors  # type: ignore
from tqdm.auto import tqdm


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
            self.n_clusters: Optional[int] = n_clusters if n_clusters > 0 else None
        self.do_sillhoute: bool = do_sillhoute


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

        mols: List[Mol] = [Chem.MolFromSmiles(s) for s in df.SMILES]
        if self.cluster_type == "Butina":
            return self.cluster_mols(mols, sim_cutoff=sim_cutoff)
        self.clusters = self._cluster_mols(mols)
        return self.clusters
    

    def _cluster_mols(self, mols: List[Mol], sim_cutoff: float) -> List[int]:
        """
        Cluster Mol objects

        :param mols: list of Mol objects
        :param sim_cutoff: Similarity cuttoff metric ONLY FOR Butina clustering

        :returns: list containing the cluster number for each compound


        *** FACTOR IN PCA REDUCTION ***
        """
        self.fp_list = self.get_fps(mols)
        if self.reduce:
            self.reduced_fps = self._reduce_hits(self.fp_list)
            fps = self.reduced_fps
        else:
            fps = self.fp_list

        cluster_method = {
            "Butina": self._cluster_butina(fps, sim_cutoff),
            "KMeans": self._cluster_kmeans(fps)
        }
        if cluster_method[self.cluster_type] is None:
            raise ValueError(f"Clustering method {self.cluster_type} is not supported.")

        return cluster_method[self.cluster_type]
    
    
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
    ) -> List[int]:
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

        self.butina_clusters_ = [x - 1 for x in tqdm(cluster_ids, desc="Assigning clusters")]
        return self.butina_clusters_
    

    def _cluster_kmeans(self, fps) -> List[int]:
        """
        Cluster compounds using the Butina clustering algorithm

        :param fps: a list of fingerprints (ExplicitBitVects) for each compound

        :returns: a list of cluster numbers for each compound

        ***PLACEHOLD SELF.N_CLUSTERS FOR NOW. CHANGE TO BE 10% OF TOTAL COMPOUNDS)
        """
        if self.n_clusters is None:
            self.n_clusters = 20

        k_means = KMeans(n_clusters=self.n_clusters, random_state=42, verbose=True)
        self.k_clusters_ = k_means.fit_predict(self.fp_list)
        return self._k_clusters_
    

    def _calc_sillhoutte(self, X:  np.stack, min: int, max: int) -> DataFrame:
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
            n_components: int | float=0.95
    ) -> List[ExplicitBitVect]:
        """
        Do a PCA reduction on the chemical fingerprints

        :param fps: Chemical bit fingerprints
        :param n_components: Number of components to retain in PCA reduction

        :returns: PCA reduced list of chemical bit fingerprints
        """
        pca = PCA(n_components=n_components)
        reduced = pca.fit_transform(fps)
        return reduced

        

        





    





