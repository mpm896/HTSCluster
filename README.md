# HTSCluster
## A CLI program for unsupervised clustering of chemical compounds, particularly libraries and hits from high throughput screens.
### This program has two main components:
1. Cluster (`ChemicalCluster` object) - Cluster your compounds
2. Query (`Query` object) - Search for similar compounds/nearest neighbors of query compounds

### Input compounds must be in SMILES format. 
### Input files can be .csv or .xlsx format.

## TODO
### Priority:
1. Write up main script for CLI program with argparse or typer -> IN PROGRESS (argparse - Clustering done, need to implement Query)
2. Move bulk of main into separate scripts: prepare_data, process_clusters, process_queries, etc 
3. Add testing for clustering and for remainder of Query
4. Update docs on how to install and use

### Extras:
1. Ensure that Silhouette analysis works
