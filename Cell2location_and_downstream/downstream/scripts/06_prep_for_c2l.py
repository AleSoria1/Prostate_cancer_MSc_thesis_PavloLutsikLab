import scanpy as sc
import os
import matplotlib
matplotlib.use('Agg')

# --- 1. SETUP ---
sc.settings.figdir = './figures_qc_and_c2l'
if not os.path.exists(sc.settings.figdir): os.makedirs(sc.settings.figdir)

MASTER_FILE = "/path/to/data/Analyses_Erickson_Junyi/Thesis_Raw_Combined_Master.h5ad"
C2L_READY_FILE = "/path/to/data/Analyses_Erickson_Junyi/Thesis_QC_C2L_Ready.h5ad"

print("🚀 Loading Raw Master Object...")
adata = sc.read_h5ad(MASTER_FILE)
print(f"Initial shape: {adata.shape}")

# --- 2. QC METRICS & BEFORE PLOTS ---
print("⚙️ Calculating QC metrics...")
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save="_QC_before_filtering.png")

# --- 3. APPLY JUSTIFIED THRESHOLDS ---
print("✂️ Filtering (< 500 genes, > 20% MT)...")
sc.pp.filter_cells(adata, min_genes=500)
adata = adata[adata.obs['pct_counts_mt'] < 20].copy()
sc.pp.filter_genes(adata, min_cells=10)
print(f"✅ Post-QC Shape: {adata.shape}")

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save="_QC_after_filtering.png")

# --- 4. SECURE RAW COUNTS FOR CELL2LOCATION ---
print("🔒 Saving raw integer counts to adata.layers['counts']...")
# Cell2location strictly requires unnormalized integers. We hide them here so Scanpy doesn't alter them.
adata.layers["counts"] = adata.X.copy()

# --- 5. PREPROCESS FOR FILTERED BATCH EFFECT PLOTS ---
print("🧬 Normalizing main layer for UMAP generation...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# CRITICAL FIX: Removed subset=True. 
# This flags the top 2000 genes in adata.var['highly_variable'] but keeps all 34,180 genes in the matrix.
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')

print("🏗️ Computing PCA and UMAP on filtered real cells...")
# CRITICAL FIX: Added use_highly_variable=True. 
# PCA will only use the 2000 genes to build the UMAP, leaving the full counts layer intact for Cell2location.
sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=True)

# Plot PCA variance ratio to mathematically justify the n_pcs=30 parameter
sc.pl.pca_variance_ratio(adata, log=True, save="_variance_ratio.png")

sc.pp.neighbors(adata, n_pcs=30)
sc.tl.umap(adata)

print("🎨 Plotting Cleaned Batch Effect...")
sc.pl.umap(adata, color=['study', 'tissue_status', 'patient'], 
           title=['Filtered: Study', 'Filtered: Status', 'Filtered: Patient'],
           save="_filtered_batch_effect.png", ncols=3)

# --- 6. SAVE ---
print(f"💾 Saving C2L-Ready Atlas to {C2L_READY_FILE}...")
adata.write(C2L_READY_FILE)
print("🎉 Done! Ready to build the Single-Cell Reference.")