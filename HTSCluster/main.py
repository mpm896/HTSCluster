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

from .cluster import ChemicalCluster
from .parser import parse
from .prepare import check_input, prepare_clusters, prepare_query
from .process import process_clusters
# from .query import Query
from .utils.utils import write_file


# Gather the CLI args and ensure that either a file name, or the hits or lib are specified
args = parse()
check_input(args)
if not args.query:
    df_clusters = prepare_clusters(args)
    cluster_choice = [args.cluster_choice]
    if args.filename: cluster_choice = ['input']
    elif args.cluster_choice == 'both':
        cluster_choice = ['hits', 'lib']
    for file in cluster_choice:
        if df_clusters[file]['df'] is not None:
            df_clusters[file]['df'] = process_clusters(df_clusters[file])
            out_name = f'{args.out_path}/{file}-output{args.out_format}'
            
            print('----------------')
            print(df_clusters)
            print('----------------')
            print(f'Writing {out_name}.....')
            write_file(df_clusters[file]['df'], out_name)
         