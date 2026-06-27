import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import os

# --- 1. SETUP & PATHS ---
IN_DIR = "/path/to/data/Cell2location/spatial_mapping"
SPATIAL_MAPPED_FILE = os.path.join(IN_DIR, "spatial_mapped_with_posterior.h5ad")
OUT_DIR = os.path.join(IN_DIR, "plots")
os.makedirs(OUT_DIR, exist_ok=True)

print(f"Loading mapped spatial data: {SPATIAL_MAPPED_FILE}")
adata_vis = sc.read_h5ad(SPATIAL_MAPPED_FILE)

# --- 2. EXTRACT STRICT THRESHOLD (q05) ---
q05_df = adata_vis.obsm['q05_cell_abundance_w_sf'].copy()
prefix = 'q05cell_abundance_w_sf_' 
q05_df.columns = [col.replace(prefix, '') for col in q05_df.columns]

for col in q05_df.columns:
    adata_vis.obs[col] = q05_df[col]

# --- 3. BIOLOGICAL VALIDATION LINEAGES ---
cell_types_to_plot = [
    'epithelial_luminal',    
    'epithelial_basal',      
    'fibroblast',            
    'smooth_muscle',         
    'macrophage_group_1'     
]

valid_cell_types = [ct for ct in cell_types_to_plot if ct in adata_vis.obs.columns]

# --- 4. GENERATE & SAVE PLOTS (CORRECTED FOR MULTI-SAMPLE) ---
sc.settings.set_figure_params(dpi=200, color_map='magma')

# Get a list of unique sample/slide IDs
sample_ids = adata_vis.obs['sample'].unique()
print(f"Found {len(sample_ids)} individual Visium slides.")

# To avoid generating 100 plots (5 lineages * 20 slides), let's plot the first slide
# You can change target_sample to any specific slide ID you want to inspect
target_sample = sample_ids[0] 
print(f"Generating plots for sample: {target_sample}")

# Slice the data to strictly one slide
adata_slide = adata_vis[adata_vis.obs['sample'] == target_sample].copy()

for cell_type in valid_cell_types:
    fig, ax = plt.subplots(figsize=(6, 6))
    
    sc.pl.spatial(
        adata_slide, 
        color=cell_type,
        cmap='magma', 
        vmin=0,           
        vmax='p99',       
        spot_size=100,    # You may need to adjust this (e.g., 50 or 150) now that it's a single slide
        show=False, 
        ax=ax
    )
    
    out_path = os.path.join(OUT_DIR, f"spatial_mapping_q05_{cell_type}_{target_sample}.png")
    plt.savefig(out_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Saved: {out_path}")

print("✅ Spatial abundance plotting complete.")