import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
from cell2location.models import RegressionModel

# --- 1. SETUP & PATHS ---
REF_DIR = "/path/to/data/Cell2location/reference_model"
ADATA_FILE = f"{REF_DIR}/reference_processed_with_posterior.h5ad"
SIG_FILE = f"{REF_DIR}/inferred_reference_signatures.csv"
OUT_DIR = f"{REF_DIR}/reference_plots"

os.makedirs(OUT_DIR, exist_ok=True)

# --- 2. BIOLOGICAL SANITY CHECK (HEATMAP) ---
print("Loading inferred signatures for biological validation...")
inf_aver = pd.read_csv(SIG_FILE, index_col=0)

# We search for canonical prostate and structural markers to verify the signatures
expected_markers = ['KLK3', 'AR', 'PTPRC', 'VIM', 'CD68', 'CD3E', 'EPCAM', 'ACTA2']
valid_markers = [g for g in expected_markers if g in inf_aver.index]

if valid_markers:
    print(f"Plotting signature heatmap for identified markers: {valid_markers}")
    plt.figure(figsize=(14, 6))
    
    # IMPROVED HEATMAP SETTINGS
    sns.heatmap(
        np.log1p(inf_aver.loc[valid_markers]), 
        cmap='YlOrRd',           # Changed to a crisp light-to-dark palette
        robust=True, 
        linewidths=0.1,          # Adds a subtle grid line between cells
        linecolor='lightgrey'    # Keeps grid lines clean and unobtrusive
    )
    
    plt.title("Biological Sanity Check: Inferred Signatures of Key Markers (Log Scale)", fontsize=14, pad=15)
    plt.xlabel("Inferred Cell Types", fontsize=12)
    plt.ylabel("Marker Genes", fontsize=12)
    plt.savefig(f"{OUT_DIR}/biological_marker_heatmap.png", bbox_inches='tight', dpi=300)
    plt.close()

# --- 3. MATHEMATICAL RECONSTRUCTION CHECK ---
print(f"Loading reference AnnData (this takes a moment, ~7.2GB)...")
adata_ref = sc.read_h5ad(ADATA_FILE)

print("Loading saved neural network weights...")
mod = RegressionModel.load(f"{REF_DIR}/trained_reference_model", adata_ref)

print("Populating posterior samples for the QC plot (Full Atlas)...")
# Drawing 20 samples across the full atlas to ensure matrix alignment
mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 20, 'batch_size': 2500}
)

print("Generating Reconstruction QC Plot...")
mod.plot_QC()
# Use plt.savefig directly to capture the active canvas
plt.savefig(f"{OUT_DIR}/reconstruction_accuracy_QC.png", bbox_inches='tight')
plt.close()

print(f"✅ Diagnostic plots successfully generated in: {OUT_DIR}")