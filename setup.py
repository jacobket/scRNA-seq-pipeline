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
        "psycopg2-binary",
    ],
    entry_points={
        "console_scripts": [
            "run-pipeline=scRNAseq_pipeline.main:main"
        ],
    },
    description="Automated pipeline for scRNA-seq preprocessing and analysis",
    author="Jacob Ketchum",
)