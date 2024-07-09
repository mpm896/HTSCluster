import argparse
from pathlib import Path
from typing import Optional, TypedDict

import polars as pl
from polars import DataFrame
from polars import selectors as cs 

from .cluster import ChemicalCluster
from .query import Query


class DataFrameClusters(TypedDict):
    df: Optional[DataFrame]
    cluster: Optional[ChemicalCluster]


def check_input(args: argparse.Namespace) -> None:
    """ Ensure that files are properly input """
    # If clustering
    if not args.query:
        if (not args.filename) and not (args.hits or args.lib):
            raise ValueError('Must provide at least one file, either by itself or by specifying the hits and the library with --hits and --lib, respectively')

        if args.filename and (args.hits or args.lib):
            raise ValueError('Provide either one file, or specify both a hits and library file.')

        if args.filename and args.hits and args.lib:
            raise ValueError('Too many input files! Can only enter 2 at max')
    # If Querying
    else:
        # Should be filename and EITHER hits or lib
        if not args.filename or sum(map(bool, [args.hits, args.lib])) == 0:
            raise ValueError('Must pass in a filename containing the SMILES to search for neighbors. Must also provide a hits or library file (or both) with --hits and --lib, respectively')
    


def file_to_df(file: str) -> DataFrame:
    """
    Read in file as Polars DataFrame 
    
    :param file: name of a file
    :return DataFrame: Polars DataFrame from your file
    """
    if file is None:
        return

    if Path(file).suffix == '.csv':
        return fix_column_names(
            (pl.read_csv(file, truncate_ragged_lines=True)
                             .select(cs.contains('SMILES'))
                             .drop_nulls())
        )
    elif Path(file).suffix == '.xlsx':
        return fix_column_names(
            (pl.read_excel(file)
                             .select(cs.contains('SMILES'))
                             .drop_nulls())
        )
    else:
        raise Exception('Not a supported file type. File types accepted are .csv and .xlsx (Excel)')


def fix_column_names(df: DataFrame, name: str='SMILES') -> DataFrame:
    """
    Rename column to desired name if the column name has a potential error, like an extra space

    :param df: Polars DataFrame
    :param name: str: desired column name
    """
    cols = df.columns
    for col in cols:
        if name in col:
            selected = col
            break
    return df.rename({selected: name})


def prepare_clusters(args: argparse.Namespace) -> dict[str, DataFrameClusters]:
    """
    Prepare the dataset of input file DataFrames and instantiated ChemicaCluster objects

    :param args: input CLI args
    :returns dict[str, DataFrameClusters]: Dict with each entry being the file entered with associated DataFrame and ChemicalCluster
    """
    # Prepare the ChemicalClusters
    inputdf, hits, lib = (
        file_to_df(args.filename.name) if args.filename else None,
        file_to_df(args.hits.name) if args.hits else None,
        file_to_df(args.lib.name) if args.lib else None
    )
    df_clusters: dict[str, DataFrameClusters] = {
        'input': {
            'df': inputdf,
            'cluster': None
            },
        'hits': {
            'df': hits,
            'cluster': None
            },
        'lib': {
            'df': lib,
            'cluster': None
            }
    }

    for df in df_clusters:
        if df_clusters[df]['df'] is not None:
            df_clusters[df]['cluster'] = ChemicalCluster(
                fptype=args.fptype,
                cluster_type=args.cluster_type,
                reduce=args.reduce,
                n_clusters=args.n_clusters,
                do_sillhoute=args.do_silhouette
            )
    return df_clusters
    

def prepare_query(args: argparse.Namespace) -> tuple[Query, DataFrame, Optional[DataFrame]]:
    """
    Prepare the dataset of input file DataFrames and instantiated Query objects

    :param args: input CLI args
    """
    # Prepare the Query
    input_smiles = Query.from_file(args.filename.name, fptype=args.fptype)
    if sum(map(bool, [args.hits, args.lib])) == 1:
        lib = (
            file_to_df(args.lib.name) if args.lib 
            else file_to_df(args.hits.name) if args.hits else None
        )
        hits = None
    elif args.lib and args.hits:
        lib = file_to_df(args.lib.name)
        hits = file_to_df(args.hits.name)

    if lib is None:
        raise ValueError('Must provide either --hits or --lib')
    
    return input_smiles, lib, hits

