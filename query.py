"""
Functions for querying nearest neighbors of a compound

For example, if you want to search for all hits that contain a dihydroquinolone backbone, 
input the SMILES of a dihydroquinolone and get the most related compounds back
"""
from __future__ import annotations

from os import PathLike
from pathlib import Path
from typing import List, Optional

import polars as pl

from HTSCluster.cluster import ChemicalCluster
from HTSCluster.utils.utils import (
    insert_clusters, 
    insert_mols,
    mols_to_img,
    write_csv, 
    write_xlsx
)

class Query:

    def __init__(self, SMILES: List[str]) -> None:
        self.SMILES: List[str] = SMILES

    @classmethod
    def from_file(cls, filaname: PathLike) -> Query:
        file = Path(filaname)
        ext = file.suffix
        if not file.exists():
            print(f"{file} does not exist. Check to ensure you spelled the file
                   correctly and have entered the correct path to the file.")
        assert ext in ['.xlsx', '.csv'], "File must be .xlsx or .csv format"

        if ext == '.xlsx':
            df = pl.read_excel(file, truncate_ragged_lines=True)
        elif ext == '.csv':
            df = pl.read_csv(file, truncate_ragged_lines=True)

        # Grab the SMILES column
        if 'SMILES' in df.columns:
            df = list(pl.select(['SMILES'])
                        .drop_nulls())
        elif len(df.columns) == 1:
            df = list(df)
            if 'SMILES' in df:
                df.remove('SMILES')
        else:
            raise Exception("""Unclear format of data. Either Give Heading name
                            'SMILES' or input file with only one column""")
        
        return cls(df)

