"""
HTSCluster: A Python program to cluster chemical compounds and search for similar compounds
within a larger library

Initial intended use was to cluster chemical compounds from a highthroughput screen, and then
to search for similar compounds from a set of input compounds

@author: Matt Martinez
https://github.com/mpm896/HTSCluster
"""
from pathlib import Path
# import time

# import numpy as np
import polars as pl
from polars import DataFrame
from polars import selectors as cs

from .cluster import ChemicalCluster
from .parser import parse
# from .query import Query
from .utils.utils import (
    fix_column_names,
    # insert_clusters, 
    # insert_mols,
    # mols_to_img,
    # write_csv, 
    # write_xlsx
)
# from .utils.polars_xlsx import xlsx_from_polarsdf

def file_to_df(file: str) -> DataFrame:
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
    file.close()


# Gather the CLI args
args = parse()

# Ensure that either a file name, or the hits or lib are specified
if (not args.filename) and not (args.hits or args.lib):
    raise Exception('Must provide at least one file, either by itself or by specifying the hits and the library with --hits and --lib, respectively')

if args.filename and (args.hits or args.lib):
    raise Exception('Provide either one file, or specify both a hits and library file.')

if args.filename and args.hits and args.lib:
    raise Exception('Too many input files! Can only enter 2 at max')

inputdf, hits, lib = (
    file_to_df(args.filename.name) if args.filename else None,
    file_to_df(args.hits.name) if args.hits else None,
    file_to_df(args.lib.name) if args.lib else None
)
df_clusters = {
    'file': {
        'df': inputdf,
        'cluster': None},
    'hits': {
        'df': hits,
        'cluster': None
        },
    'lib': {
        'df': lib,
        'cluster': None}
}

# Do clustering
if not args.query:
    for df in df_clusters:
        if df_clusters[df]['df'] is not None:
            clusters = ChemicalCluster(
                fptype=args.fptype,
                cluster_type=args.cluster_type,
                reduce=args.reduce,
                n_clusters=args.n_clusters,
                do_sillhoute=args.do_sillhouette
            )

# Do query for nearest neighbors
else:
    pass