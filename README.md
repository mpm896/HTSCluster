# HTSCluster
## A CLI program for unsupervised clustering of chemical compounds, particularly libraries and hits from high throughput screens.
### This program has two main components:
1. Cluster (`ChemicalCluster` object) - Cluster your compounds
2. Query (`Query` object) - Search for similar compounds/nearest neighbors of query compounds

### Input compounds must be in SMILES format. 
### Input files can be .csv or .xlsx format.

## TODO
### Priority:
1. Add testing for clustering and for remainder of Query
2. Save the dataframe of neighbors of a query, with Mol images
3. Fix Query to produce one dataframe per query SMILES
4. Write up main script for CLI program with argparse or typer

### Extras:
1. Ensure that Silhouette analysis works