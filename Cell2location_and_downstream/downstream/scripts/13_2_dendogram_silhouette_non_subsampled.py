import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.metrics import silhouette_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
os.environ["OMP_NUM_THREADS"] = "16"
os.environ["MKL_NUM_THREADS"] = "16"
os.environ["OPENBLAS_NUM_THREADS"] = "16"

OUT_DIR = ("/path/to/data/"
           "Cell2location/spatial_mapping/exploration_mean/"
           "clustering_comparison_nt13")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Load data ────────────────────────────────────────────────────────
print("Loading adata...")
adata = sc.read_h5ad(
    "/path/to/data/"
    "Cell2location/spatial_mapping/exploration_mean/"
    "spatial_compositional_space_mean.h5ad"
)

X = adata.X.toarray() if sp.issparse(adata.X) else adata.X
print(f"Data shape: {X.shape}")

# ── PART 1: Full Ward linkage and dendrogram ─────────────────────────
print(f"\nComputing Ward linkage on full {X.shape[0]:,} spots...")
Z = linkage(X, method='ward', metric='euclidean')
print("Linkage computed.")

# Save linkage matrix
np.save(os.path.join(OUT_DIR, 'ward_linkage_full.npy'), Z)
print("Linkage matrix saved.")

# Plot dendrogram
fig, ax = plt.subplots(figsize=(14, 7))
dendrogram(
    Z,
    truncate_mode='lastp',
    p=30,
    leaf_rotation=90,
    leaf_font_size=9,
    ax=ax,
    color_threshold=0.7 * max(Z[:, 2]),
    above_threshold_color='grey'
)
ax.set_title(
    'Hierarchical Clustering Dendrogram (Ward linkage)\n'
    f'Full dataset: {X.shape[0]:,} spots | Truncated to top 30 merges\n'
    'Large vertical gaps = natural cluster boundaries',
    fontsize=11
)
ax.set_xlabel('Cluster (spot count in brackets)', fontsize=10)
ax.set_ylabel('Ward Distance', fontsize=10)

for k, color, ls in [(8, 'red', '--'), (10, 'blue', '--'), (12, 'green', '--')]:
    threshold = Z[-(k-1), 2]
    ax.axhline(y=threshold, color=color, linestyle=ls, alpha=0.7,
               label=f'k={k} cut')
ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'dendrogram_ward_full.png'),
            dpi=150, bbox_inches='tight')
print("Saved: dendrogram_ward_full.png")

# ── PART 2: Silhouette on ALL spots ─────────────────────────────────
print(f"\nComputing silhouette scores on full {X.shape[0]:,} spots...")
print("This may take 30-60 minutes...")

cluster_cols = {
    'Leiden res=0.3': 'leiden_res0.3',
    'Leiden res=0.5': 'spatial_niche',
    'Leiden res=0.7': 'leiden_res0.7',
    'Leiden res=0.8': 'leiden_res0.8',
    'Leiden res=0.9': 'leiden_res0.9',
    'HC Ward k=7':  'hc_ward_k7',
    'HC Ward k=8':  'hc_ward_k8',
    'HC Ward k=9':  'hc_ward_k9',
    'HC Ward k=10': 'hc_ward_k10',
    'HC Ward k=11': 'hc_ward_k11',
    'HC Ward k=12': 'hc_ward_k12',
}

sil_results = []

for name, col in cluster_cols.items():
    if col not in adata.obs.columns:
        print(f"  Skipping {col} — not found")
        continue

    labels = adata.obs[col].values
    n_unique = len(np.unique(labels))
    if n_unique < 2:
        print(f"  Skipping {name} — only {n_unique} cluster")
        continue

    print(f"  Computing silhouette for {name}...")
    sil = silhouette_score(X, labels, metric='euclidean', n_jobs = 16)
    n_clusters = adata.obs[col].nunique()
    counts = adata.obs[col].value_counts()
    max_pct = counts.max() / len(adata) * 100

    sil_results.append({
        'Method': name,
        'N clusters': n_clusters,
        'Silhouette': round(sil, 4),
        'Largest niche (%)': round(max_pct, 1)
    })
    print(f"    {name}: silhouette={sil:.4f}, n={n_clusters}")

sil_df = pd.DataFrame(sil_results)
print("\n── Silhouette Summary (full dataset) ──")
print(sil_df.sort_values('Silhouette', ascending=False).to_string(index=False))

# Save results table
sil_df.to_csv(os.path.join(OUT_DIR, 'silhouette_full_dataset.csv'), index=False)
print("Saved: silhouette_full_dataset.csv")

# Plot
fig, ax = plt.subplots(figsize=(10, 5))
colors = ['#2196F3' if 'Leiden' in m else '#FF9800'
          for m in sil_df['Method']]
bars = ax.barh(sil_df['Method'], sil_df['Silhouette'],
               color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
ax.axvline(x=0, color='black', linewidth=0.8)
ax.set_xlabel('Silhouette Score', fontsize=11)
ax.set_title(
    'Cluster Quality: Silhouette Score (full dataset — no subsampling)\n'
    'Blue = Leiden | Orange = Hierarchical',
    fontsize=11
)
for bar, val in zip(bars, sil_df['Silhouette']):
    ax.text(val + 0.002, bar.get_y() + bar.get_height()/2,
            f'{val:.4f}', va='center', fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'silhouette_full_dataset.png'),
            dpi=150, bbox_inches='tight')
print("Saved: silhouette_full_dataset.png")

print("\nAll done.")