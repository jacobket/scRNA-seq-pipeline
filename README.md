### Bioinformatics Pipeline For Single Cell RNA-Seq Data

This pipeline automates the preprocessing and analysis of single-cell RNA sequencing (scRNA-seq) data. It integrates quality control, normalization, dimensionality reduction, clustering, differential expression analysis (DEA), and SQL-based data storage to streamline large-scale single-cell analysis.

## Dataset Requirements
This pipeline is designed for scRNA-seq datasets that contain the following three files:
* `matrix.mtx.gz`: Matrix of raw gene expression counts
* `features.tsv.gz`: List of gene features
* `barcodes.tsv.gz`: List of cell barcodes

## Example Datasets
Suitable datasets can be found on platforms like GEO. To download from GEO, modify the following URL template with your dataset ID:
`ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSEXXXnnn/GSEXXXXXX/suppl/GSEXXXXXX_matrix.mtx.gz`

## Running the Pipeline
Once an appropriate dataset has been found (from GEO, etc.) setup the function call as follows:
`python preprocess.py --mtx_url <matrix_url> --barcodes_url <barcodes_url> --features_url <features_url>`

**An example function call is as follows:**
`python preprocess.py --mtx_url ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE268nnn/GSE268249/suppl/GSE268249_matrix.mtx.gz \
--barcodes_url ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE268nnn/GSE268249/suppl/GSE268249_barcodes.tsv.gz \
--features_url ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE268nnn/GSE268249/suppl/GSE268249_features.tsv.gz`

Optionally, you can append the flag `--clear_db` at the end of your function call to clear the backing SQL database and refresh the schema to reflect any changes. 

## Pipeline Stages
1. **Automated Quality Control (QC)**: Filters cells based on gene count thresholds, mitochondrial content, and total counts, ensuring that low-quality cells are excluded from the dataset.
2. **Normalization and Scaling**:  Normalizes total counts, applies log transformation, identifies highly variable genes, and scales data for dimensionality reduction.
3. **Dimensionality Reduction with PCA/UMAP**: Performs PCA to gather major sources of variance in the dataset and UMAP to visualize the data.
4. **Clustering**: Groups cells into distinct clusters using the Leiden algorithm to identify subpopulations based on gene expression profiles.
5. **Differential Expression Analysis (DEA)**: Identifies genes with significant expression differences between clusters, providing insights into cluster-specific markers and biological processes.
6. **SQL Database Storage**: Stores processed results (clusters and DEGs) in a **PostgreSQL** database using **SQLAlchemy**, enabling efficient querying and data management for larger datasets.
7. **Visualizations**: Outputs various plots for initial data exploration, including violin plots for QC metrics, PCA variance ratios, and UMAP embeddings colored by top marker genes.

## Output 
- A **preprocessed AnnData object** with embeddings, clusters, and QC metrics is saved as `data/results_file.h5ad`.
- Processed **clustering and DEG results** are stored in a **PostgreSQL database** for efficient querying and retrieval.

## Installation Information
In order to automate installation of necessary packages using setup.py, run the following at the command line: 
`python setup.py install`

This will install:
* Scanpy (for scRNA-seq analysis)
* AnnData (for managing single-cell data)
* Pandas, NumPy (for data manipulation)
* PostgreSQL dependencies (psycopg2-binary)
* SQLAlchemy (for database integration)

## Setting up PostgreSQL for the Pipeline
To use the SQL storage feature, ensure that PostgreSQL is installed and running on your local environment. The pipeline will automatically create the necessary database and tables on first run. 

# Steps to Setup PostgreSQL:
1. Ensure PostgreSQL is installed and running:
* On Linux: `sudo systemctl start postgresql`
* On macOS: `brew services start postgresql`

## Additional Pipeline Features In Progress
This pipeline will soon include the following features:

1. **Dockerization**: The entire pipeline will be containerized using Docker, allowing easy deployment and ensuring reproducibility across environments.
2. **Batch Effect Corrections**: Add support for removing batch effects in data.
... **More to be determined**
