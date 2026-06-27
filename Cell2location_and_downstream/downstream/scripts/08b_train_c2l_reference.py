import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp
import os
import cell2location

# --- 1. SETUP ---
# CRITICAL FIX: Pointing to the newly created stencil file, NOT the raw file.
REF_FILE = "/path/to/data/Cell2location/reference_model/c2l_curated_reference.h5ad"
OUT_DIR = "/path/to/data/Cell2location/reference_model"
os.makedirs(OUT_DIR, exist_ok=True)

print(f"Loading curated reference atlas: {REF_FILE}")
adata_ref = sc.read_h5ad(REF_FILE)

# --- 2. MATRIX OPTIMIZATION FOR TENSOR CORES ---
# QC filtering is intentionally omitted here because the stencil already applied the curated thresholds.
print("Casting matrix to Compressed Sparse Row (CSR) format for GPU optimization...")
if sp.issparse(adata_ref.X):
    adata_ref.X = adata_ref.X.tocsr()

print(f"Matrix dimensions (Curated Stencil): {adata_ref.shape}")

# Verify the unnormalized integer counts layer exists
if "counts" not in adata_ref.layers:
    adata_ref.layers["counts"] = adata_ref.X.copy()

# --- 3. CELL2LOCATION MODEL SETUP ---
print("Registering AnnData with Cell2location covariates...")
cell2location.models.RegressionModel.setup_anndata(
    adata=adata_ref,
    layer="counts",
    labels_key="annotation_manual_fine",
    batch_key="sample_id"
)

# --- 4. MODEL TRAINING ---
print("Initializing Negative Binomial Regression Model...")
mod = cell2location.models.RegressionModel(adata_ref)

print("Commencing model training (max_epochs=250)...")
mod.train(max_epochs=250)

# --- 5. POSTERIOR EXTRACTION ---
print("Exporting posterior distributions (num_samples=250 to prevent WICE OOM)...")
adata_ref = mod.export_posterior(
    adata_ref, 
    sample_kwargs={'num_samples': 250, 'batch_size': 2500}
)

# Extract the inferred average expression (the reference signatures)
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg']
    inf_aver.columns = adata_ref.uns['mod']['factor_names']
    
    sig_path = os.path.join(OUT_DIR, "inferred_reference_signatures.csv")
    inf_aver.to_csv(sig_path)
    print(f"✅ Reference signatures successfully extracted and saved to: {sig_path}")
else:
    print("CRITICAL ERROR: Posterior extraction failed to compute means_per_cluster.")

# Save the trained model object and the updated AnnData for reproducibility
mod.save(os.path.join(OUT_DIR, "trained_reference_model"), overwrite=True)
adata_ref.write(os.path.join(OUT_DIR, "reference_processed_with_posterior.h5ad"))

print("Execution complete.")