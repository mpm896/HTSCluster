from typing import Optional

from polars import DataFrame

from .prepare import DataFrameClusters
from .query import Query
from .utils.utils import insert_clusters


def process_clusters(clusters: DataFrameClusters) -> DataFrame | None:
    """
    Cluster the data

    :param clusters: DataFrameClusters dict containing the DataFrame of chemical compounds and the ChemicalCluster object
    :returns DataFrame: DataFrame of chemical compounds with assigned clusters and molecule images
    """

    df = clusters["df"]
    chem_cluster = clusters["cluster"]

    print("----- CLUSTERING -----")
    print(f"{df = }")

    if df is not None and chem_cluster is not None:
        cluster_nums = chem_cluster.cluster_smiles(df)
        return insert_clusters(chem_cluster, df)

    return None


def process_query(
    query: Query, df: DataFrame, n_neighbors: int = -1, hits: Optional[DataFrame] = None
) -> dict[str, DataFrame]:
    """
    Query for similar compounds

    :param query: Query object containing the SMILES of interest
    :param df: Polars DataFrame of chemical compounds

    :returns DataFrame: DataFrame of chemical compounds with assigned clusters and molecule images
    """
    queried = query.query_neighbors(df=df, n_neighbors=n_neighbors)

    # If the hits were provided
    if hits is not None:
        assigned = query.assign_query_hits(hits)
        return assigned

    return queried
