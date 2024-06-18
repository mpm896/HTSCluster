# %%
from os import PathLike
from pathlib import Path
from typing import List
from PIL import Image
from PIL.Image import Image

import numpy as np
import pandas as pd
import polars as pl
from polars import DataFrame
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import Mol, Draw

from cluster import ChemicalCluster

# %%
def insert_clusters(clusters: ChemicalCluster, df: DataFrame) -> DataFrame:
    """ Insert column containing the clusters """
    cluster_type = clusters.cluster_type
    return df.with_columns(
        pl.lit(clusters.clusters)
        .alias(f"{cluster_type} Cluster")
    )


def insert_mols(mols: List[np.ndarray], df: DataFrame) -> pd.DataFrame:
    """ 
    Insert Mols into the DataFrame. This MUST be done with Pandas
    and NOT Polars because Polars doesn't seem to support this datatype
    """
    column_name: str = "Molecule"

    new_df = df.with_columns(
        pl.Series(mols, dtype=pl.Object)
        .alias(column_name)
    )
    return new_df


def mols_to_img(mols: List[Mol]) -> List[Image]:
    """ Convert MOLs to PNG Images """
    return [Draw.MolToImage(m) for m in mols]


# %%
def write_csv(
        df: DataFrame, 
        path: PathLike=Path.cwd(), 
        base_name: str="cluster_results.csv"
) -> None:
    df.write_csv(f"{path}/{base_name}")


def write_xlsx(
        df: DataFrame,
        path: PathLike=Path.cwd(),
        base_name: str="cluster_results.xlsx"
) -> None:
    """ 
    Write to an Excel file. If has Molecule images,
    use rdkit PandasTools to save the file
    """
    molCol = "Molecule"
    if molCol in df.columns:
        PandasTools.SaveXlsxFromFrame(
            df, f"{path}/{base_name}", molCol=molCol
        )
        return
    df.write_excel(f"{path}/{base_name}")

# %%
