FROM ubuntu:22.04

RUN apt-get update && apt-get install -y python3 python3-pip wget
RUN pip install scanpy anndata psycopg2-binary seaborn sqlalchemy 
# dependencies are ScanPy and Anndata

COPY preprocess.py /
COPY models.py /

# final configuration
CMD ["python3", "preprocess.py", "--mtx_url", "mtx_url", "--barcodes_url", "barcodes_url", "--features_url", "features_url", "--step", "qc"]
# each argument is a separate string