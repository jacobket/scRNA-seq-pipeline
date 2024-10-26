import scanpy as sc
import pandas as pd
import numpy as np
import os
import argparse

def get_param():
    """
    Sets up command-line arguments for running the scRNA-seq pipeline using argparse.
    
    Returns:
        argparse.Namespace: An object containing the parsed command-line args, including URLs for 
        matrix, barcodes, and feature files.
    
    """
    parser = argparse.ArgumentParser(description="Download GEO dataset files")
    parser.add_argument('--mtx_url', type=str, required=False, help="URL for the matrix file")
    parser.add_argument('--barcodes_url', type=str, required=False, help="URL for the barcodes file")
    parser.add_argument('--features_url', type = str, required=False, help="URL for the features file")
    return parser.parse_args()

def download_files(mtx_url=None, barcodes_url=None, features_url=None):
    """
    Downloads scRNA-seq dataset files (matrix, barcodes, and features) from provided URLs into a local folder.
    Checks for existing local files to avoid redundant downloads when re-running the pipeline.
    
    Parameters:
        mtx_url (str or None): URL for the matrix file; if None, this file is not downloaded.
        barcodes_url (str or None): URL for the cell barcodes file; if None, this file is not downloaded.
        features_url (str or None): URL for the gene features file; if None, this file is not downloaded.
    
    Returns:
        None
    """
    urls = {"Matrix": mtx_url, "Barcodes": barcodes_url, "Features": features_url} # Dictionary, key-value pairs
    filenames = {
        "Matrix": f"{prefix_path}matrix.mtx.gz",
        "Barcodes": f"{prefix_path}barcodes.tsv.gz",
        "Features": f"{prefix_path}features.tsv.gz"
    }

    for file_type, url in urls.items(): # urls.items() iterates over key-value pairs
        if url:
            file_path = os.path.join("data", filenames[file_type])
            print(f"File path is: {file_path}")
            if not os.path.exists(file_path): 
                print(f"Downloading {file_type} file...")
                os.system(f"wget -P data/ {url}") # Download file to 'data' folder
            else:
                print(f"{file_type} file already exists, skipping download.")
        else:
            print(f"No {file_type} URL provided, skipping...")

def load_and_plot_qc(adata):
    """
    Performs quality control (QC) on single cell RNA-seq data by identifying mitochondrial genes, calculating QC metrics, 
    and creating a violin plot for gene counts, total UMI counts, and % of mitochondrial content across cells.
    
    Parameters: 
        adata (AnnData): Annotated data matrix with scRNA-seq counts.
    Returns: 
        adata (AnnData): Updated data matrix with added QC metrics for downstream filtering.
    """
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # Identify mitochondrial genes (pct_counts_mt)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter=0.4, multi_panel=True)
    return adata

def filter_data_based_on_qc(adata):
    """
    Automatically filters cells and genes based on QC thresholds for gene count, total counts, and mitochondrial content %.
    
    Parameters:
        adata (AnnData): Annotated data matrix with scRNA-seq counts.
    Returns:
        adata (AnnData): Filtered data matric after applying automated QC thresholds.
    """
    # Initial baseline filtering
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # Automatically use %-based thresholds
    print("Running QC with automated thresholds...")

    lower_n_genes_threshold = adata.obs['n_genes_by_counts'].quantile(0.05)
    upper_n_genes_threshold = adata.obs['n_genes_by_counts'].quantile(0.95) # to filter doulets/multiplets (multiple cells counted as one)
    # Gene count value that 95% of cells are below, to filter abnormally high gene counts

    # total_counts generalized filtering
    min_total_counts = adata.obs['total_counts'].quantile(0.05) # Value below which 5% of total_counts values fall
    max_total_counts = adata.obs['total_counts'].quantile(0.95)

    # pct_counts_mt generalized filtering
    max_pct_mt = adata.obs['pct_counts_mt'].quantile(0.5)  # Using median

    # Applying filtering
    adata = adata[(adata.obs['n_genes_by_counts'] > lower_n_genes_threshold) & (adata.obs['n_genes_by_counts'] < upper_n_genes_threshold), :]
    adata = adata[(adata.obs['total_counts'] > min_total_counts) & (adata.obs['total_counts'] < max_total_counts), :]
    adata = adata[adata.obs['pct_counts_mt'] < max_pct_mt, :]
    return adata

def normalize_and_scale_data(adata):
    """
    Performs normalization, log-transformation, identifies highly variable genes, and scales the scRNA-seq dataset.
    
    Parameters:
        adata (AnnData): Annotated data matric with scRNA-seq counts.
    Returns:
        adata (AnnData): Preprocessed data matric, filtered to highly variable genes and scaled for downstream analysis.
    """
    sc.pp.normalize_total(adata, target_sum=1e4) #Total count normalize data matrix to 10000 reads / cell, so counts are comparable across cells
    sc.pp.log1p(adata) # Logarithmize the data to reduce influence of outlier cells
    
    #Identify highly variable genes (preferred for studying differences across cell populations)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pl.highly_variable_genes(adata)
    
    # Make raw copy for downstream use during Diff Expression analysis
    adata.raw = adata.copy()
    adata = adata[:, adata.var.highly_variable] # Filter to keep highly variable genes
    sc.pl.highest_expr_genes(adata, n_top=20)
    # Look for highly expressed genes driving variability in dataset, can be cell-type specific marker genes
    # Check for genes that may dominate dataset unusually due to batch effects
    
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
    sc.pp.scale(adata, max_value=10)
    return adata

def find_elbow_pt(explained_var, n_points=30):
    """
    Finds the elbow point in the variance curve using the perpendicular line method.
    
    Parameters:
        explained_var (array): Variance explained by each principal component.
        n_points (int): Number of PCs to consider for finding the elbow (default is 30).
    Returns:
        int: The index of the PC representing the elbow point.
    """
    # Use only the first 'n_points' PCs, default is 30
    explained_var = explained_var[:n_points]
    x_coords = np.arange(n_points)
    
    #Line endpoints: First and last points in variance plot
    p1 = np.array([x_coords[0], explained_var[0]])
    p2 = np.array([x_coords[-1], explained_var[-1]])
    
    #Calculate distances from each point to the line defined by p1 and p2
    distances=[]
    for i in range(n_points):
        p = np.array([x_coords[i], explained_var[i]])
        # Calculate the perpendicular distance to the line
        distance = np.abs(np.cross(p2 - p1, p1 - p)) / np.linalg.norm(p2-p1)
        distances.append(distance)
    # Find index of max distance from line (elbow pt)
    elbow_pt = np.argmax(distances) + 1
    
    print(f"Elbow point found at PC {elbow_pt} based on the perpendicular line method.")
    return elbow_pt

def dimensionality_reduction(adata):
    """
    Reduces the dimensionality of the scRNA-seq data using PCA and UMAP,
    storing the results in the AnnData object for downstream analysis and visualization.

    Parameters:
        adata (AnnData): Annotated data matrix with single-cell RNA-seq counts.

    Returns:
        AnnData: Modified data matrix with PCA and UMAP embeddings.
    """
    # Step 1: Run PCA
    sc.tl.pca(adata, svd_solver="arpack")
    # Step 2: Automate selection of color gene for visualization
    top_gene = adata.var.highly_variable.idxmax() #Captures genes driving main sources of biovariability (gene expression that distinguishes cell types)
    sc.pl.pca(adata, color=top_gene) # ! Can update later to include a marker gene for specific cell type
    # Step 3: Plot variance ratios and determine num of PCs to use
    #import pdb; pdb.set_trace()
    sc.pl.pca_variance_ratio(adata, log=True) # x-axis is the PCs, y-axis is the variance explained by each
    explained_var = adata.uns['pca']['variance_ratio'] #Access stored annotations regarding variance explained by each PC
    n_pcs = find_elbow_pt(explained_var, n_points=30)
    print(f"Using {n_pcs} PCs for downstream analysis")
    adata.write("data/results_file.h5ad") #! Save current AnnData object
    # Step 4: Compute neighborhood graph based on optimal n_pcs
    # Graph w cell-to-cell relationships by connecting each cell to its nearest neighbor
    sc.pp.neighbors(adata, n_neighbors = 15, n_pcs=n_pcs) #!Fine tune n_neighbors based on nature of clusters in dataset
    # Step 5: Visualize using UMAP
    sc.tl.umap(adata) # Uses neighborhood graph from sc.pp.neighbors to calc UMAP coords for subsequent plotting
    # Select the top 3 highly variable genes
    top_genes = adata.var[adata.var.highly_variable].sort_values(by="means", ascending=False).head(3).index.tolist()
    # Plot UMAP colored w 3 most variable genes to reveal patterns/clusters and ID diff cell types/ explain biovariability
    sc.pl.umap(adata, color=top_genes, use_raw=False) # use_raw=False bc we want to use normalized, scaled data
    return adata
    
if __name__ == "__main__":
    args= get_param()
    prefix_path = args.mtx_url.split("/")[-1].split('_')[0] + '_'
    download_files(args.mtx_url, args.barcodes_url, args.features_url)
    adata = sc.read_10x_mtx('data', var_names='gene_symbols', cache=True, prefix = prefix_path) # Load the files into Scanpy
    adata = load_and_plot_qc(adata)
    adata = filter_data_based_on_qc(adata)
    adata = normalize_and_scale_data(adata)
    adata = dimensionality_reduction(adata)
    # Next- clustering method

#First step- Preprocessing: Quality control, normalization, scaling, filtering
#Second step- Dimensionality Redution- PCA and UMAP
# Third step- Clustering the neighborhood graph to id distinct subpopulations of cancer cells

# sc.tl.leiden( # Leiden for clustering 
#     adata,
#     resolution=0.9,
#     random_state=0,
#     flavor="igraph",
#     n_iterations=2,
#     directed=False,
# )
# sc.pl.umap(adata, color=[top_genes])
# adata.write('data/preprocessed_data.h5ad') # ?

# Differential Expression analysis, visualizations
# Use t-test to rank highly differential genes in each cluster
# sc.tl.rank_genes_groups(adata, "leiden", method="t-test")
# sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# sc.settings.verbosity = 2  # reduce the verbosity

# Using Wilcoxon method
# sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
# sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False) - write 

# sc.tl.rank_genes_groups(adata, "leiden", method="logreg", max_iter=1000)
# sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

#DESEQ2- using R
# Marker gene identification
# Database Integration uding PostgreSQL
# Future: DE analysis, results storing using SQL db, dockerization