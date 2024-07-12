"""
HTSCluster: A Python program to cluster chemical compounds and search for similar compounds
within a larger library

Initial intended use was to cluster chemical compounds from a highthroughput screen, and then
to search for similar compounds from a set of input compounds

@author: Matt Martinez
https://github.com/mpm896/HTSCluster
"""
import time
from .parser import parse
from .prepare import check_input, prepare_clusters, prepare_query
from .process import process_clusters, process_query
from .utils.utils import get_mols, insert_mols, write_file

DEBUGGING = False

def main():
    # Gather the CLI args and ensure that either a file name, or the hits or lib are specified
    args = parse()
    check_input(args)

    # Do the clustering
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


            if DEBUGGING:
                import pandas as pd
                import polars as pl
                from rdkit.Chem import PandasTools
                from .process import process_clusters_deprecated
                from .utils.utils import write_file_deprecated
                if df_clusters[file]['df'] is not None:
                    df_clusters[file]['df'] = process_clusters_deprecated(df_clusters[file])
                    out_name = f'{args.out_path}/{file}-output{args.out_format}'
                    
                    print('----------------')
                    print(df_clusters)
                    print('----------------')
                    print(f'Writing {out_name}.....')

                    start = time.time()
                    PandasTools.SaveXlsxFromFrame(pl.DataFrame.to_pandas(df_clusters[file]['df']), out_name, molCol='Molecule', size=(300,300))
                    print(f"{time.time() - start} seconds to write xlsx with PandasTools")

                    write_file_deprecated(df_clusters[file]['df'], out_name)
    # Do the query
    else:
        query, lib, hits = prepare_query(args)
        queried_smiles = process_query(query, lib, args.n_neighbors, hits)

        # Get Mols, save the file
        for smiles in queried_smiles:
            mols = get_mols(queried_smiles[smiles])
            df_with_mols = insert_mols(mols, queried_smiles[smiles])
            out_name = f'{args.out_path}neighbors-{smiles}{args.out_format}'

            print('----------------')
            print(f'Writing {out_name}.....')
            print('----------------')
            write_file(df_with_mols, out_name)


if __name__ == "__main__":
    main()