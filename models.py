from sqlalchemy import Column, Integer, String, Float, ForeignKey, create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker

Base = declarative_base()

# Table #1: Store cluster information
class Cluster(Base):
    __tablename__ = 'cluster'
    id = Column(Integer, primary_key=True)
    run_id = Column(Integer, nullable = False)
    cluster_id = Column(String, nullable=False) # ID for leiden cluster this entry corresponds to (0,1,2...)
    num_cells = Column(Integer) # total num of cells in the cluster
    avg_expression = Column(Float) # avg gene expression across cells in the cluster
    def __repr__(self):
        return f"<Cluster(id={self.id}, cluster_id={self.cluster_id}, cell_count={self.num_cells})>"

# Table #2: Store DEG results
class DEG(Base):
    __tablename__ = 'deg'
    id = Column(Integer, primary_key=True)
    run_id = Column(Integer, nullable = False)
    cluster_id = Column(Integer, ForeignKey('cluster.id'), nullable=False) # id col of cluster table
    gene_id = Column(String, nullable=False)
    logFC = Column(Float, nullable=False)
    p_value = Column(Float, nullable=False)
    adj_p_value = Column(Float, nullable=False)
    def __repr__(self):
        return f"<DEG(id={self.id}, gene_name={self.gene_id}, log_fold_change={self.logFC})>"

# Define relationships 
Cluster.degs = relationship("DEG", back_populates="cluster")
DEG.cluster = relationship("Cluster", back_populates="degs")
