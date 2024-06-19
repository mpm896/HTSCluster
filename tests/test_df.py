import pytest

import polars as pl

from HTSCluster import Query

def sample_df():
    pass


def test_read_queries(
        query_file_csv,
        query_file_xlsx,
        query_file_multicol,
):
    """ Assert Query is properly constructes"""
    QUERY = [
    "CC1=C(C=C2NC(CC(C2=C1)C3=CCC=C3)=O)O",
    "O=C(NC1=CC(O)=C(C=C21)C)CC2C3=CC=CC=C3"
    ]
    query_csv = Query.from_file(query_file_csv)
    query_xslx = Query.from_file(query_file_xlsx)
    query_multicol = Query.from_file(query_file_multicol)
    
    assert query_csv.SMILES == QUERY
    assert query_xslx.SMILES == QUERY
    assert query_multicol.SMILES == QUERY


def test_query_file_error(query_file_error):
    """ Test that an incorrect file is recognized """
    msg = ("""Unclear format of data. Either Give Heading name
                            'SMILES' or input file with only one column""")
    with pytest.raises(Exception, match=msg):
        Query.from_file(query_file_error)