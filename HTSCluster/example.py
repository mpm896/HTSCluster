# %%
import time

import numpy as np
import polars as pl

from cluster import ChemicalCluster
from query import Query
from utils.utils import (
    insert_clusters, 
    insert_mols,
    mols_to_img,
    write_csv, 
    write_xlsx
)
from utils.polars_xlsx import xlsx_from_polarsdf

PATH_TO_CSV = "/Users/U1036725/Documents/PersonalProjects/HTSCluster/tests/data/Chembrigde_Div.csv"
QUERY_CSV = "../tests/data/query.csv"
QUERY_XLSX = "../tests/data/query.xlsx"

# %%
start = time.time()
df = (pl.read_csv(PATH_TO_CSV, truncate_ragged_lines=True)
                .select(['SMILES'])
                .drop_nulls()[:200])
print(f"{time.time() - start} seconds to read CSV with Polars")

""" --------------------------------------------------------- """
""" Running to make sure the ChemicalCluster class is working """
""" --------------------------------------------------------- """

# %%
df.head
# %%
lib_cluster = ChemicalCluster(cluster_type="Butina")
clusters = lib_cluster.cluster_smiles(df)

# %%
new_df = insert_clusters(lib_cluster, df)
print(new_df)
# %%
# Test writing to CSV and Excel
"""write_csv(new_df)
write_xlsx(new_df)"""

# %%
# Test getting Mols and inserting into DataFrame
mols = lib_cluster.get_mols(new_df)
df_with_mols = insert_mols(mols, new_df)

# %%
xlsx_from_polarsdf(df_with_mols, "testFile.xlsx", "Molecule")

# %%
""" ----------------------------------------------- """
""" Running to make sure the Query class is working """
""" ----------------------------------------------- """
# %%
import time

import numpy as np
import polars as pl
from polars import selectors as cs

from cluster import ChemicalCluster
from query import Query
from utils.utils import (
    insert_clusters, 
    insert_mols,
    mols_to_img,
    write_csv, 
    write_xlsx
)
from utils.polars_xlsx import xlsx_from_polarsdf

PATH_TO_CSV = "/Users/U1036725/Documents/PersonalProjects/HTSCluster/tests/data/Chembrigde_Div.csv"
PATH_TO_HITS = "/Users/U1036725/Documents/PersonalProjects/HTSCluster/tests/data/HTS_SG.csv"
QUERY_CSV = "../tests/data/query.csv"
QUERY_XLSX = "../tests/data/query.xlsx"

# %%
start = time.time()
df = (pl.read_csv(PATH_TO_CSV, truncate_ragged_lines=True)
                .select(cs.contains('SMILES'))
                .drop_nulls()) # [:5000])
print(f"{time.time() - start} seconds to read CSV with Polars")

start = time.time()
hits = (pl.read_csv(PATH_TO_HITS, truncate_ragged_lines=True)
                .select(cs.contains('SMILES')))

SMILES = [
    "CC1=C(C=C2NC(CC(C2=C1)C3=CCC=C3)=O)O",
    "O=C(NC1=CC(O)=C(C=C21)C)CC2C3=CC=CC=C3"
]
query = Query(SMILES=SMILES)
# %%
query_csv = Query.from_file(QUERY_CSV)
query_xlsx = Query.from_file(QUERY_XLSX)
# %%
start = time.time()
query_csv.query_neighbors(df=df, n_neighbors=5000)
print(f"{time.time() - start} seconds to query for neighbors")
# %%
# Rename the column if it is not 'SMILES'
name = hits.columns[0]
if name != 'SMILES':
    hits = hits.rename({name: 'SMILES'})
hits

# %%
assigned = query_csv.assign_query_hits(hits)
# %%
