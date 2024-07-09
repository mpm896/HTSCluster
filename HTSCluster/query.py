"""
Functions for querying nearest neighbors of a compound

For example, if you want to search for all hits that contain a dihydroquinolone backbone, 
input the SMILES of a dihydroquinolone and get the most related compounds back
"""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import List #, Optional

import numpy as np
import polars as pl
from polars import DataFrame
from rdkit import Chem
from rdkit.Chem import AllChem, Mol
from sklearn.neighbors import NearestNeighbors  # type: ignore
from tqdm.auto import tqdm

from .utils.utils import get_fps

# from HTSCluster.cluster import ChemicalCluster
"""from HTSCluster.utils.utils import (
    insert_clusters, 
    insert_mols,
    mols_to_img,
    write_csv, 
    write_xlsx
)"""

class Query:

    def __init__(
            self,
            SMILES: List[str] | DataFrame,
            fptype: str="rdkit"
    ) -> None:
        self.SMILES: List[str] | DataFrame = (DataFrame(
            {"SMILES": SMILES}) if isinstance(SMILES, list) else SMILES
        )
        self.fptype: str = fptype


    def __repr__(self):
        return f"{self.SMILES}"


    @classmethod
    def from_file(cls, filaname: str | Path, fptype: str="rdkit") -> Query:
        file = Path(filaname)
        ext = file.suffix
        if not file.exists():
            print(f"""{file} does not exist. Check to ensure you spelled the file
                   correctly and have entered the correct path to the file.""")
        assert ext in ['.xlsx', '.csv'], "File must be .xlsx or .csv format"

        if ext == '.xlsx':
            df = pl.read_excel(file)
        elif ext == '.csv':
            df = pl.read_csv(file, truncate_ragged_lines=True)

        # Grab the SMILES column
        if len(df.columns) > 1:
            if 'SMILES' in df.columns:
                smiles = df.select('SMILES')
                return cls(smiles, fptype)
            else:
                raise Exception(
                    "Unclear format of data. Either Give Heading name 'SMILES' or input file with only one column"
                )
        return cls(df, fptype)


    @property
    def _mols(self) -> List[Mol]:
        return [Chem.MolFromSmiles(s) for (s) in self.SMILES['SMILES']]
    

    @property
    def _fps(self):
        return get_fps(self.fptype, self._mols)
    

    def query_neighbors(
            self,
            df: DataFrame,
            n_neighbors: int=50
    ) -> dict[str, DataFrame]:
        """
        Get the closest compounds to your query compound
        
        :param df: Polars DataFrame
            List of SMILES from which you want to find the most similar ones
        :param n_neighbors: int
            n most similar compounds to return

        TODO: STILL NEED TO TEST THIS -> Seems to be working based on example.py
        """
        assert 'SMILES' in df.columns, "No SMILES column"
        if n_neighbors == -1:
            n_neighbors = len(df)

        _lib_mols = [Chem.MolFromSmiles(s) for (s) in tqdm(
            df['SMILES'], desc='Getting mols from SMILES'
        )]
        _lib_fps = get_fps(self.fptype, _lib_mols)

        # Initialize nearest neighbors and find k nearest neighbors
        neigh = NearestNeighbors(n_neighbors=n_neighbors).fit(_lib_fps)

        self.query_neighbors = {}
        for i in tqdm(range(len(self.SMILES)), desc="Querying for nearest neighbors..."):
            kneighbors = neigh.kneighbors(np.array(self._fps[i]).reshape(1,-1))
            self.query_neighbors[self.SMILES[i].item()] = DataFrame(df[kneighbors[1][0], :])
        return self.query_neighbors
    

    def assign_query_hits(self, df: DataFrame) -> dict[str, DataFrame]:
        """
        Assign each neighbor as a hit or not

        :param df: Polars DataFrame
            SMILES of the hit compounds

        TODO: Optimize this to use Polars queries instead of vanilla python
        TODO: STILL NEED TO TEST THIS
        """
        assert 'SMILES' in df.columns

        for c in self.query_neighbors.copy():
            df_assigned = DataFrame(self.query_neighbors[c]['SMILES']).with_columns(
                pl.lit(None, dtype=pl.String)
                .alias('IS_HIT')
            )

            # Is this the efficient way to query for rows using Polars?
            for compound in df_assigned['SMILES']:
                if compound in self.SMILES['SMILES']:
                    df_assigned = df_assigned.with_columns(
                        (pl.when(pl.col('SMILES') == compound)
                                .then(pl.lit("QUERY"))
                                .otherwise(pl.col('IS_HIT')))
                            .alias('IS_HIT')
                    )
                elif compound in df['SMILES']:
                    df_assigned = df_assigned.with_columns(
                        (pl.when(pl.col('SMILES') == compound)
                                .then(pl.lit("HIT"))
                                .otherwise(pl.col('IS_HIT')))
                            .alias('IS_HIT')
                    )
            
        return self.query_neighbors



