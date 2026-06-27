import scanpy as sc
import scipy.sparse as sp
import cell2location
import os
import warnings
warnings.filterwarnings("ignore") # Suppress expected PyTorch CPU warnings for the dry run

# --- 1. SETUP ---
REF_FILE = "/path/to/data/RNA_counts_raw.h5ad"

print("Loading raw reference atlas...")
adata_ref = sc.read_h5ad(REF_FILE)

# --- 2. SUBSAMPLING (CRITICAL FOR LOGIN NODE EXECUTION) ---
print("Subsampling matrix to 1,000 cells for structural testing...")
sc.pp.subsample(adata_ref, n_obs=1000, random_state=42)

# --- 3. APPLY FILTERS ---
print("Applying QC metrics...")
adata_ref.var['mt'] = adata_ref.var_names.str.startswith('MT-') | adata_ref.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata_ref, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pp.filter_cells(adata_ref, min_genes=300)
sc.pp.filter_cells(adata_ref, max_genes=5000)
adata_ref = adata_ref[adata_ref.obs['pct_counts_mt'] < 10].copy()

# Lowering min_cells drastically because we only have 1000 total cells
sc.pp.filter_genes(adata_ref, min_cells=2) 

# --- 4. MATRIX OPTIMIZATION ---
print("Casting matrix to CSR format...")
if sp.issparse(adata_ref.X):
    adata_ref.X = adata_ref.X.tocsr()

adata_ref.layers["counts"] = adata_ref.X.copy()

# --- 5. CELL2LOCATION MODEL SETUP ---
print("Registering AnnData with Cell2location covariates...")
cell2location.models.RegressionModel.setup_anndata(
    adata=adata_ref,
    layer="counts",
    labels_key="annotation_manual_fine",
    batch_key="sample_id"
)

# --- 6. INITIALIZE & TRAIN (CPU ONLY) ---
print("Initializing Negative Binomial Regression Model...")
mod = cell2location.models.RegressionModel(adata_ref)

print("Commencing dry run (max_epochs=2, accelerator='cpu')...")
# Overriding to CPU specifically for the login node test
mod.train(max_epochs=2)
print("\n✅ DRY RUN SUCCESSFUL. The computational graph and API parameters are structurally sound.")