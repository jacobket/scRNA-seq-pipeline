FROM ubuntu:22.04.4

RUN apt-get update && apt-get install -y python3 python3-pip
RUN pip install scanpy anndata
# dependencies are ScanPy and Anndata

COPY preprocess.py /

# final configuration
CMD ["python3", "preprocess.py", "--mtx_url", "mtx_url", "--barcodes_url", "barcodes_url", "--features_url", "features_url", "--step", "qc"]
# each argument is a separate string