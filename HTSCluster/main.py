"""
HTSCluster: A Python program to cluster chemical compounds and search for similar compounds
within a larger library

Initial intended use was to cluster chemical compounds from a highthroughput screen, and then
to search for similar compounds from a set of input compounds

@author: Matt Martinez
https://github.com/mpm896/HTSCluster
"""
import time

import numpy as np
import polars as pl

from .cluster import ChemicalCluster
from .parser import parse
from .query import Query
from .utils.utils import (
    insert_clusters, 
    insert_mols,
    mols_to_img,
    write_csv, 
    write_xlsx
)
from .utils.polars_xlsx import xlsx_from_polarsdf

args = parse()
if args.filename:
    print(args.filename)
else:
    print("no")

