# Multimodal Spatial and Epigenomic Analysis of Prostate Cancer

**MSc Bioinformatics Thesis — KU Leuven, 2026**

**Supervisors:** Prof. Pavlo Lutsik | MSc Guillaume Sacchetti  
**Department:** Computational Cancer Biology and Epigenomics

---

## Scientific Summary

Prostate cancer is characterised by profound intratumoural heterogeneity, yet the spatial organisation of its cellular and molecular components within intact tissue remains incompletely understood. Bulk RNA-sequencing loses cellular resolution by averaging across millions of cells, single-cell RNA-sequencing dissolves the spatial context essential for understanding tissue architecture, and DNA methylation profiling — while capturing stable epigenomic regulatory states independent of transient transcriptional activity — lacks the cellular and spatial resolution to localise epigenomic programmes to specific tissue compartments. This thesis addressed these complementary limitations by integrating Visium spatial transcriptomics, single-cell RNA-seq deconvolution, bulk transcriptomics and DNA methylation epigenomics in a unified multimodal framework.

Two Visium datasets were integrated — a morphologically normal prostate atlas and a prostate cancer dataset — yielding 101,490 quality-controlled spots across 38 slides. Cell2location deconvolution against a 232,255-cell reference atlas estimated absolute abundances of 22 cell types per spot. A joint compositional-spatial graph clustering approach outperformed all composition-only methods, identifying five biologically interpretable niches: fibromuscular stroma, luminal epithelium, periurethral immune-ductal epithelium, high-grade aggressive tumour, and an incidental cancer niche validated against histopathological annotations. Clinical analysis across 707 PPCG patients revealed a two-phase stromal remodelling model — early smooth muscle loss followed by late cancer-associated fibroblast and macrophage accumulation — with fibroblast proportion independently prognostic beyond Gleason grade. Epigenomic LMC analysis identified desmoplastic and neuro-stromal remodelling programmes as the dominant risk signatures. Most interestingly, the incidental cancer niche was protective in bulk RNA-seq survival yet carried the highest spatial activity for risk LMC programmes, demonstrating that spatial transcriptomics detects localised epigenomic reprogramming diluted and masked in bulk sequencing.

---

## Repository Structure
├── Cell2location_and_downstream/

│   └── downstream/

│       ├── notebooks/          # Jupyter notebooks (numbered by analysis step)

│       └── scripts/            # Python and R scripts for cluster submission

│

└── ST_Erickson_and_JunyiHu_analyses/

# Early exploratory analyses of spatial datasets
> **Note:** All data files (`.h5ad`, `.csv`, `.tsv`, `.qs`, `.RDS`) and figures (`.png`, `.pdf`) are excluded from this repository. Patient-level data from the PPCG cohort is not publicly shareable. Spatial transcriptomics data is available from the original publications (see Data Availability).

---

## Analysis Pipeline

| Step | Notebook/Script | Description |
|------|----------------|-------------|
| 01 | `11_explore_spatial_abundances_entropies_batch_effects.ipynb` | Data integration, QC, batch effect assessment |
| 02 | `08b_train_c2l_reference.py` | Cell2location reference model training |
| 03 | `09_run_c2l_spatial.py` | Cell2location spatial mapping |
| 04 | `11_2_joint_graph_and_clustermap.ipynb` | Compositional space construction, joint graph clustering |
| 05 | `13_cluster_comparison_leiden_hierarchical.ipynb` | Clustering benchmarking |
| 06 | `14_niche_biology.ipynb` | Niche characterisation: DEGs, identity scoring, z-scores |
| 07 | `14_2_niches_H&E.ipynb` | H&E overlay validation |
| 08 | `16_2_PPCG.ipynb` | PPCG clinical correlations: Kruskal-Wallis, Mann-Whitney |
| 09 | `16_4_Celltype_proportions_survivals.ipynb` | Cell type survival analysis: univariate and multivariate Cox |
| 10 | `16_5_LMC_proportions_survivals.ipynb` | LMC survival analysis |
| 11 | `16_3_SpatialDEGsvsLMCs.ipynb` | AUCell LMC scoring in spatial data |
| 12 | `16_8_Spatial_scores_survivals.ipynb` | Spatial niche signature scoring in bulk RNA-seq + survival |
| 13 | `14_gsea_lmc.py` + `14_GSEA.slurm` | GSEA pathway enrichment (submitted as cluster job) |
| 14 | `16_7_Naming_LMCs.ipynb` | LMC annotation: GSEA + cell type correlation clustermap |
| 15 | `16_6_Liana.ipynb` | LIANA ligand-receptor interaction analysis |

---

## Methods Summary

| Component | Details |
|-----------|---------|
| **Data** | Erickson et al. 2022 (prostate cancer, 2 patients, 22 Visium slides) + Hu et al. 2025 (normal prostate, 7 patients, 16 slides) |
| **Deconvolution** | cell2location v0.1.5 — Bayesian model, 232,255-cell scRNA-seq reference atlas, 22 cell types |
| **Clustering** | Joint graph: A_joint = α·A_comp + (1−α)·A_spat, benchmarked against Leiden and HC Ward |
| **Clinical cohort** | PPCG — 707 patients, bulk RNA-seq (BayesPrim, 14 cell types), methylation array (MeDeCom, 14 LMCs), 158 relapse events |
| **Survival analysis** | Rank-based inverse normal transformation + Cox PH (lifelines) + Kaplan-Meier |
| **Pathway enrichment** | Preranked GSEA (gseapy) — MSigDB Hallmarks, Reactome 2022, GO BP 2023, CellMarker 2024 |
| **Cell-cell interactions** | LIANA v1.7.1, Rank Aggregate consensus method |

---

## Environment

```bash
# Python (conda)
conda activate thesis
# Python 3.11 — key packages: scanpy 1.11.5, cell2location 0.1.5,
# squidpy 1.6.5, lifelines, gseapy, liana 1.7.1

# R (for LMC extraction — HPC module)
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2
# Key packages: MeDeCom, RnBeads, qs, SummarizedExperiment
```

All analyses were performed on a high-performance computing cluster using SLURM. GPU-intensive steps required an NVIDIA A100-80GB GPU.

---

## Data Availability

| Dataset | Access |
|---------|--------|
| Erickson et al. 2022 spatial data | [Mendeley Data](https://doi.org/10.17632/svw96g68dv.1) |
| Hu et al. 2025 spatial data | [figshare](https://doi.org/10.6084/m9.figshare.25965613) |
| PPCG cohort | Available through ICGC data access — not publicly shareable |

---

## Key Results

| Analysis | Key Finding |
|----------|------------|
| Joint graph clustering | α=0.9 outperforms composition-only; 5 biological niches identified |
| Cell type survival | Fibroblast independent of Gleason grade (HR=1.39, p=0.03) |
| LMC survival | LMC7 HR=1.72 (desmoplastic/neurogenic); LMC8 HR=0.67 (protective) |
| Niche signature survival | Niche 5 and 6 protective in bulk (HR=0.69, 0.63) |
| AUCell | Niche 5 highest risk LMC activity despite low-grade transcriptome |
| LIANA | Niche 2 immune evasion: NECTIN4→TIGIT, HLA-F→LILRB2 |

---

