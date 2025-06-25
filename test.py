from sqlalchemy import create_engine, inspect, func
from sqlalchemy.orm import sessionmaker
from models import Cluster, DEG # import table models

engine = create_engine('postgresql+psycopg2://db_user:db_pass@localhost:5432/scrna_seq_db')

inspector = inspect(engine)
tables = inspector.get_table_names()

print("Tables in the database:", tables)

Session = sessionmaker(bind=engine)
session = Session()

# Check data in Cluster table
clusters = session.query(Cluster).all()
print(f"Number of clusters: {len(clusters)}")

unique_clusters = session.query(Cluster.cluster_id, func.count()).group_by(Cluster.cluster_id).all()
print("Unique Clusters:", len(unique_clusters))

for cluster in unique_clusters:
    print(cluster)
    
unique_degs = session.query(DEG.cluster_id, func.count()).group_by(DEG.cluster_id).all()
print("DEGs per Cluster:", unique_degs)
# for deg in session.query(DEG).limit(10).all():
#     print(deg)
# Check data in DEG table
degs = session.query(DEG).all()
print(f"Number of DEGs: {len(degs)}")

if clusters:
    print("Sample cluster entry:", clusters[0])
if degs: 
    print("Sample DEG entry:", degs[0])