### Bioinformatics Pipeline For Single Cell RNA-Seq Data

This pipeline automates the preprocessing, quality control, normalization, and dimensionality reduction of single-cell RNA-seq data, providing an efficient and reproducible workflow for initial data exploration and analysis.

## Prerequisites

To use this pipeline, install the necessary packages:
```bash
pip install scanpy anndata pandas numpy
```
Alternatively, use the `setup.py` file (see instructions below) for automated package installation.

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

## Pipeline Stages
1. **Automated Quality Control (QC)**: Filters cells based on gene count thresholds, mitochondrial content, and total counts, ensuring that low-quality cells are excluded from the dataset.
2. **Normalization and Scaling**:  Normalizes total counts, applies log transformation, identifies highly variable genes, and scales data for dimensionality reduction.
3. **Dimensionality Reduction with PCA/UMAP**: Performs PCA to gather major sources of variance in the dataset and UMAP to visualize the data.
4. **Visualizations**: Outputs violin plots for QC metrics, PCA variance ratios, and UMAP embeddings colored by top genes for an initial exploration.

## Output 
The preprocessed AnnData object with embeddings and QC metrics is saved in `data/results_file.h5ad` for easy access in subsequent analysis steps.

## Using setup.py 
In order to automate installation of necessary packages using setup.py, run the following at the command line: 
`python setup.py install`

## Additional Pipeline Features In Progress
This pipeline will soon include the following features:

1. **Clustering**: After dimensionality reduction, the pipeline will cluster cells based on gene expression profiles to identify distinct cell populations.

2. **Differential Expression Analysis (DEA)**: This step will identify genes with significant expression differences between clusters, providing insights into cluster-specific markers and biological processes.

3. **SQL Database Storage**: Processed results and analysis outputs will be stored in a SQL database for efficient querying and retrieval, enhancing data management for larger datasets.

4. **Dockerization**: The entire pipeline will be containerized using Docker, allowing easy deployment and ensuring reproducibility across environments.