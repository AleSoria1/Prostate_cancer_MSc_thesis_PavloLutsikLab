import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp
import os
import matplotlib.pyplot as plt
import cell2location

# --- 1. SETUP & PATHS ---
SPATIAL_FILE = SPATIAL_FILE = "/path/to/data/Analyses_Erickson_Junyi/Thesis_QC_C2L_Ready_Fixed.h5ad"
SIG_FILE = "/path/to/data/Cell2location/reference_model/inferred_reference_signatures.csv"
OUT_DIR = "/path/to/data/Cell2location/spatial_mapping"
os.makedirs(OUT_DIR, exist_ok=True)

print(f"Loading spatial atlas: {SPATIAL_FILE}")
adata_vis = sc.read_h5ad(SPATIAL_FILE)

print(f"Loading reference signatures: {SIG_FILE}")
inf_aver = pd.read_csv(SIG_FILE, index_col=0)

# --- 2. GENE INTERSECTION ---
print("Intersecting spatial genes with reference signatures...")
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()
print(f"Retained {len(intersect)} overlapping genes for mapping.")

# --- 3. MATRIX OPTIMIZATION ---
print("Formatting raw integer counts for GPU optimization...")
# Enforce strict integer casting for the Gamma-Poisson distribution
adata_vis.layers["counts"] = sp.csr_matrix(adata_vis.layers["counts"]).astype(np.int32)

# --- 4. MODEL SETUP ---
print("Registering Spatial AnnData...")
# Force the batch key to be a strict categorical string
adata_vis.obs["sample"] = adata_vis.obs["sample"].astype(str).astype("category")

cell2location.models.Cell2location.setup_anndata(
    adata=adata_vis, 
    batch_key="slide_id", 
    layer="counts"
)

# Initialize mapping model
print("Initializing Spatial Mapping Model...")
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    N_cells_per_location=8, 
    detection_alpha=20
)

# --- 5. TRAINING ---
print("Commencing spatial mapping with mini-batching (max_epochs=30000)...")
# CRITICAL FIX: Slicing the 100k+ spots into 2500-spot chunks to prevent A100 OOM
mod.train(max_epochs=1000, batch_size=2500)

# --- 6. POSTERIOR EXTRACTION & THRESHOLD JUSTIFICATION ---
print("Exporting spatial posteriors (num_samples=250 to prevent CPU OOM)...")
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 250, 'batch_size': 2500}
)

print("Generating empirical plots for threshold justification...")
# Plot 1: Reconstruction accuracy (Observed vs Expected counts)
mod.plot_QC()
# Use plt.savefig to capture the active canvas instead of trying to save a NoneType object
plt.savefig(os.path.join(OUT_DIR, "reconstruction_accuracy_QC.png"), bbox_inches='tight')
plt.close()

# --- 7. SAVE OUTPUTS ---
print("Saving trained model and annotated AnnData...")
mod.save(os.path.join(OUT_DIR, "trained_spatial_model"), overwrite=True)
adata_vis.write(os.path.join(OUT_DIR, "spatial_mapped_with_posterior.h5ad"))

print("✅ Spatial mapping execution complete.")