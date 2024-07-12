# HTSCluster
## A CLI program for unsupervised clustering of chemical compounds, particularly libraries and hits from high throughput screens.
## Install:
First, [install miniconda](https://docs.anaconda.com/miniconda/miniconda-install/) if not already installed

Next, open your terminal and create a conda environment:

`conda create -n htscluster python=3.11`

Activate the conda environment:

`conda activate htscluster`

Now install the package in your conda environment `htscluster` in one of two ways:
1. `pip install git+https://github.com/mpm896/HTSCluster.git`
2. Or clone the git repository, `cd` into it, and install from there:
```
git clone https://github.com/mpm896/HTSCluster.git
cd HTScluster
pip install .
```

## Usage: Cluster
Run from the terminal with: `cluster-hits`

For help and usage: `cluster-hits -h` or `cluster-hits --help`

To run, you must provide **at least** one file when running. The file must contain chemical compounds in **SMILES** format and the file must be either .csv or .xlsx (Excel) formates

`cluster-hits filename.xlsx`

You can specify which file contains hits, and which file contains the chemical library with option `--hits` and `--lib`, respectifully:

`cluster-hits --hits hits.csv --lib lib.csv` 

## Usage: Query for similar compounds
If you have the SMILES of a chemical backbone and you want to find similar compounds in your hits or your library, you can **query** for this.

To do so, you must provide two files: one with your query SMILES, and one for the hits or the lib:

`cluster hits --query query_smiles.csv --hits hits.csv` (can use `--lib` instead if desired)

## Usage: Options
### Clustering options:
When clustering, you can select the method used to convert a SMILES to a bitwise fingerprint, either `rdkit` or `morgan`:

`-f [choice]` or `--fpytype [choice]` (where [choice] is rdkit or morgan)

You can select the clustering algorithm, `KMeans` or `Butina` (default is `Butina`):

`-c [choice]` or `--clustertype [choice]`

If using `KMeans` clustering, you can modify the number of clusters to create. Default is 10% of total compounds being clustered:

`-n [n]` or `--nclusters [n]`

If you feel you have very high chemical diversity and it will create too many clusters (mainly with the Butina clustering algorithm), you can reduce dimensionality with principal component analysis (PCA). It by default keeps 95% variance:

`-r` or `--reduce`

If you want to use an optimal number of clusters with KMeans clustering, you can perform a silhouette analysis to identify the optimal cluster number **(Functionality not yet confirmed)**:

`-s` or `--silhouette`

### Querying options:
You can specify the number of compounds to get back, ranked by similarity (default is 50). Enter -1 to get back **all** compounds, ranked by similarity to query SMILES:

`-qn [n]` or `--queryneighbors [n]`

### File saving options:
By default, the output file will be .xlsx Excel format. You can tell the program to save a .csv file, if desired:

`-o .csv` or `--out-format .csv`

You should specify the output location of the file. By default, it is saved wherever you run the script:

`-p /path/to/directory` or `--out-path /path/to/directory`

<br>
<br>

### The API has two main components:
1. Cluster (`ChemicalCluster` object) - Cluster your compounds
2. Query (`Query` object) - Search for similar compounds/nearest neighbors of query compounds

If doing an exploratory analysis, everything is centered around `ChemicalCluster` and `Query` objects

### Input compounds must be in SMILES format. 
### Input files can be .csv or .xlsx format.

## TODO
### Priority:
1. Add testing for clustering and for remainder of Query
2. Change implementation of `xlsx_from_polarsdf`

### Extras:
1. Ensure that Silhouette analysis works

### Times to write CSV
1. 315 compounds with images
PandasTools - 3.1 s
xlsx_from_polarsdf (same algo as PandasTools) - 3.0 s
xlsx_from_polarsdf (utilizing Polars methods) - 1.9 s

2. 43,000 compounds with images 
PandasTools - 656.0 s
xlsx_from_polarsdf (same algo as PandasTools) - 636.9 s
xlsx_from_polarsdf (utilizing Polars methods) - 463.6 s
