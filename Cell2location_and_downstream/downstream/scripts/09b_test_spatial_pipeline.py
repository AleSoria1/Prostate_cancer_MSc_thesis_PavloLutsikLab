import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import cell2location
import warnings
warnings.filterwarnings("ignore")

# --- 1. SETUP ---
SPATIAL_FILE = "/path/to/data/Analyses_Erickson_Junyi/Thesis_QC_C2L_Ready.h5ad" 
SIG_FILE = "/path/to/data/Cell2location/reference_model/inferred_reference_signatures.csv"

print("Loading spatial atlas and reference signatures...")
adata_vis = sc.read_h5ad(SPATIAL_FILE)
inf_aver = pd.read_csv(SIG_FILE, index_col=0)

# --- 2. SUBSAMPLING ---
print("Subsampling spatial spots to 100 for dry run...")
sc.pp.subsample(adata_vis, n_obs=100, random_state=42)

# --- 3. GENE INTERSECTION ---
print("Intersecting spatial genes with reference signatures...")
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()
print(f"Intersection complete. Retained {len(intersect)} overlapping genes.")

# --- 4. MATRIX OPTIMIZATION ---
print("Formatting raw integer counts for GPU...")
# Do NOT overwrite with adata_vis.X. 
# We pull the raw counts saved from script 06 and enforce integer type.
adata_vis.layers["counts"] = sp.csr_matrix(adata_vis.layers["counts"]).astype(np.int32)

# --- 5. MODEL SETUP ---
print("Registering Spatial AnnData...")

# CRITICAL FIX: Force the batch key to be a strict categorical string
adata_vis.obs["sample"] = adata_vis.obs["sample"].astype(str).astype("category")

cell2location.models.Cell2location.setup_anndata(
    adata=adata_vis, batch_key="sample", layer="counts"
)

# --- 6. INITIALIZATION & CPU TEST ---
print("Initializing Spatial Mapping Model...")
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    N_cells_per_location=8, detection_alpha=20
)

print("Commencing CPU dry run (max_epochs=2)...")
mod.train(max_epochs=2)

print("\n✅ SPATIAL DRY RUN SUCCESSFUL. The matrices intersect correctly.")