# %%
import time

import numpy as np
import polars as pl
import pandas as pd

from cluster import ChemicalCluster
from utils.utils import (
    insert_clusters, 
    insert_mols,
    mols_to_img,
    write_csv, 
    write_xlsx
)
from utils.polars_xlsx import xlsx_from_polarsdf

PATH_TO_CSV = "/Users/U1036725/Documents/PersonalProjects/HTSCluster/tests/data/Chembrigde_Div.csv"

# %%
start = time.time()
df = (pl.read_csv(PATH_TO_CSV, truncate_ragged_lines=True)
                .select(['SMILES'])
                .drop_nulls()[:200])
print(f"{time.time() - start} seconds to read CSV with Polars")

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
