import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.cluster import AgglomerativeClustering
import warnings
warnings.filterwarnings('ignore')

DATA_PATH = (
    "/path/to/data/"
    "Cell2location/spatial_mapping/exploration_mean/"
    "spatial_compositional_space_mean.h5ad"
)

print("Loading data...")
adata = sc.read_h5ad(DATA_PATH)
print(f"Shape: {adata.shape}")

X = adata.X.toarray() if sp.issparse(adata.X) else adata.X
print(f"CLR matrix shape: {X.shape}")

print("\nRunning hierarchical clustering (Ward linkage)...")
for k in [7, 8, 9, 10, 11, 12]:
    print(f"  Fitting k={k}...")
    hc = AgglomerativeClustering(
        n_clusters=k,
        metric='euclidean',
        linkage='ward'
    )
    labels = hc.fit_predict(X)
    col_name = f'hc_ward_k{k}'
    adata.obs[col_name] = labels.astype(str)

    counts = pd.Series(labels).value_counts().sort_index()
    max_pct = counts.max() / len(adata) * 100
    min_pct = counts.min() / len(adata) * 100
    print(f"  k={k}: largest={max_pct:.1f}%, "
          f"smallest={min_pct:.1f}%")

print("\nSaving updated h5ad with hierarchical cluster labels...")
adata.write_h5ad(DATA_PATH)
print("Done. Columns saved:")
hc_cols = [c for c in adata.obs.columns if 'hc_ward' in c]
print(hc_cols)
