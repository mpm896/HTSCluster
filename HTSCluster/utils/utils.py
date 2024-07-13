# %%
from __future__ import annotations
from os import PathLike
from pathlib import Path
from typing import List
from PIL import Image
from PIL.Image import Image
from typing import List, TYPE_CHECKING

import numpy as np
import polars as pl
from polars import DataFrame
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Mol, Draw
from tqdm.auto import tqdm

from .polars_xlsx import xlsx_from_polarsdf

if TYPE_CHECKING:
    from ..cluster import ChemicalCluster

# %%
def get_fps(fptype: str, mols: List | pl.Series) -> List[Mol]:
    """ Get chemical bit fingerprints """
    if fptype not in ['rdkit', 'morgan']:
        raise ValueError(f"Fingerprint method {fptype} is not supported.")
    
    if fptype == 'rdkit':
        return [Chem.RDKFingerprint(x) for x in tqdm(
            mols, desc="Calculating chemical fingerprints..."
        )]
    # If fptype is 'morgan'
    return [AllChem.GetMorganFingerprintAsBitVect(x, 2) for x in tqdm(
            mols, desc="Calculating chemical fingerprints..."
    )]


def get_mols(df: DataFrame, column_name: str="SMILES") -> List[Mol]:
    assert column_name in df.columns
    return [Chem.MolFromSmiles(s) for s in df[column_name]]


def insert_clusters(clusters: ChemicalCluster, df: DataFrame) -> DataFrame:
    """ Insert column containing the clusters """
    cluster_type = clusters.cluster_type
    return df.with_columns(
        pl.lit(clusters.clusters)
        .alias(f"{cluster_type} Cluster")
    )


def insert_mols(mols: List[np.ndarray], df: DataFrame) -> DataFrame:
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
    return [Draw.MolToImage(m) for m in tqdm(mols, desc='Converting to images...')]


# %%
def write_csv(
        df: DataFrame, 
        path: PathLike=Path.cwd(), 
        base_name: str="cluster_results.csv"
) -> None:
    df.write_csv(f"{path}/{base_name}")


def write_file(df: DataFrame, filename: str) -> None:
    """
    Write a file to .csv or .xlsx format
    """
    path = filename.split('/')[:-1]
    base_name = filename.split('/')[-1]
    if Path(base_name).suffix == '.csv':
        write_csv(df, path, base_name)
    elif Path(base_name).suffix == '.xlsx':
        xlsx_from_polarsdf(df, filename, molCol='Molecule')
    else:
        raise ValueError('Not a supported file format for saving.')
       


# %%
