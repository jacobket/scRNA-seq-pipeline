import scanpy as sc
import pandas as pd
import os
import gzip
import shutil

# Set paths for the dataset files
mtx_url = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSEXXnnnnn/GSEXXnnnnn/suppl/GSE268249_matrix.mtx.gz'
barcodes_url = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSEXXnnnnn/GSEXXnnnnn/suppl/GSE268249_barcodes.tsv.gz'
features_url = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSEXXnnnnn/GSEXXnnnnn/suppl/GSE268249_features.tsv.gz'

# Download function for the files
def download_file(url, save_path):
    with gzip.open(url, 'rb') as f_in:
        with open(save_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

# Create a directory to store files if it doesn't exist
if not os.path.exists('data'):
    os.makedirs('data')

# Download files
download_file(mtx_url, 'data/matrix.mtx')
download_file(barcodes_url, 'data/barcodes.tsv')
download_file(features_url, 'data/features.tsv')

# Load the files into Scanpy
adata = sc.read_10x_mtx(
    'data',  # Path where the matrix.mtx, barcodes.tsv, and features.tsv are stored
    var_names='gene_symbols',  # Use gene symbols from features.tsv
    cache=True  # Cache the result for faster access later
)

# Preprocessing: Quality control, normalization, etc.
# Filter cells and genes
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Calculate QC metrics and add them to the AnnData object
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # Mitochondrial genes start with 'MT-'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Filter out cells with too many mitochondrial genes or low gene count
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

# Normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)  # Logarithmize the data

# Identify highly-variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

# Scale the data
sc.pp.scale(adata, max_value=10)

# Save the preprocessed AnnData object
adata.write('data/preprocessed_data.h5ad')

# Plot QC metrics and check data
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
sc.pl.highest_expr_genes(adata, n_top=20)

print(adata)
