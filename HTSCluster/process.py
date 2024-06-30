import polars as pl
from polars import DataFrame

from .prepare import DataFrameClusters
from .utils.utils import insert_clusters, insert_mols

def process_clusters(clusters: DataFrameClusters) -> DataFrame:
    """
    Cluster the data

    :param clusters: DataFrameClusters dict containing the DataFrame of chemical compounds and the ChemicalCluster object
    :returns DataFrame: DataFrame of chemical compounds with assigned clusters and molecule images
    """

    df = clusters['df']
    chem_cluster = clusters['cluster']

    print("----- CLUSTERING -----")
    print(f"{df = }")
    
    clusters = chem_cluster.cluster_smiles(df)
    mols = chem_cluster.get_mols(df)
    return insert_mols(mols, insert_clusters(chem_cluster, df))
    
