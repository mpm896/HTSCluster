"""
Parse the CLI args
"""

import argparse
from pathlib import Path


def parse() -> argparse.Namespace:
    """
    Parse CLI args with argparser

    Necessary CLI Args:
        filename OR lib and hits : str
        fptype (rdkit or morgan) : str
        cluster_type (Butina or KMeans) : str
        reduce (i.e. for PCA) : bool
        N_clusters : int
        do_silhouette : bool
        query : str (can be file or just plain string)
        n_neighbors : int (number of neighbors for Query)
    """
    parser = argparse.ArgumentParser(
        prog="HTSCluster",
        description="A CLI program to cluster chemical compounds and query a library for similar compounds",
    )
    parser.add_argument(
        "filename",
        nargs="?",
        type=argparse.FileType("r"),
        help="Input file containing compounds. MUST contain a SMILES column. If you want to input both a hits file and library file, use --hits and --lib options, respectively.",
    )
    parser.add_argument(
        "--hits", type=argparse.FileType("r"), help="file containing the chemical hits"
    )
    parser.add_argument(
        "--lib",
        type=argparse.FileType("r"),
        help="file containing the chemical library",
    )
    parser.add_argument(
        "-t",
        "--to-cluster",
        default="both",
        choices=["hits", "lib", "both"],
        dest="cluster_choice",
        help="Which file to cluster: hits, library, or both.",
    )
    parser.add_argument(
        "-f",
        "--fptype",
        dest="fptype",
        default="rdkit",
        choices=["rdkit", "morgan"],
        help="Chemical fingerprint type. Default is 'rdkit'. Likely don't need to change this.",
    )
    parser.add_argument(
        "-c",
        "--clustertype",
        dest="cluster_type",
        default="Butina",
        choices=["Butina", "KMeans"],
        help="Clustering method. If the library is very high diversity, use KMeans. Otherwise Butina is probably fine.",
    )
    parser.add_argument(
        "-r",
        "--reduce",
        action="store_true",
        dest="reduce",
        help="Whether or not to do PCA reduction before clustering. Defult is True",
    )
    parser.add_argument(
        "-n",
        "--nclusters",
        dest="n_clusters",
        type=int,
        help="Number of clusters, ONLY used for KMeans clustering. Defauls it 10%% of the number of compounds being clustered.",
    )
    parser.add_argument(
        "-s",
        "--silhouette",
        dest="do_silhouette",
        action="store_false",
        help="Do silhoutte analysis to identify optimal number of clusters for KMeans clustering. Default is False.",
    )
    parser.add_argument(
        "-q",
        "--query",
        dest="query",
        action="store_true",
        help="Query sequences for which you want to find similar compounds. Can be just one SMILES, or a CSV or Excel file containing a column of SMILES",
    )
    parser.add_argument(
        "-qn",
        "--queryneighbors",
        dest="n_neighbors",
        type=int,
        default=50,
        help="Number of desired similar compounds to get back. Default is 50. For the entire list of the input hits or libraries (ranked by similarity to query compound), enter -1",
    )
    parser.add_argument(
        "-o",
        "--out-format",
        dest="out_format",
        choices=[".xlsx", ".csv"],
        default=".xlsx",
        help="Format of output file. Default is .xlsx (Excel file)",
    )
    parser.add_argument(
        "-p",
        "--out-path",
        dest="out_path",
        default=Path.cwd(),
        help="Path to output file. Default is current directory.",
    )

    return parser.parse_args()
