from setuptools import setup, find_packages

setup(
    name="scRNAseq_pipeline",
    version="1.0",
    packages=find_packages(),
    install_requires=[
        "scanpy",
        "anndata",
        "pandas",
        "numpy",
        "python-igraph",
        "psycopg2-binary"
    ],
    description="Automated pipeline for scRNA-seq preprocessing and analysis",
    author="Jacob Ketchum",
)