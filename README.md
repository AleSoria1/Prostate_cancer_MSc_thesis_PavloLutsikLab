Multimodal Spatial and Epigenomic Analysis of Prostate Cancer
MSc Bioinformatics Thesis — KU Leuven, 2026

Supervisors: Prof. Pavlo Lutsik | MSc Guillaume Sacchetti

Department: Computational Cancer Biology and Epigenomics

Scientific Summary
Prostate cancer is characterised by profound intratumoural heterogeneity, yet the spatial organisation of its cellular and molecular components within intact tissue remains incompletely understood. Bulk RNA-sequencing loses cellular resolution by averaging across millions of cells, single-cell RNA-sequencing dissolves the spatial context essential for understanding tissue architecture, and DNA methylation profiling — while capturing stable epigenomic regulatory states independent of transient transcriptional activity — lacks the cellular and spatial resolution to localise epigenomic programmes to specific tissue compartments. This thesis addressed these complementary limitations by integrating Visium spatial transcriptomics, single-cell RNA-seq deconvolution, bulk transcriptomics and DNA methylation epigenomics in a unified multimodal framework.
Two Visium datasets were integrated — a morphologically normal prostate atlas and a prostate cancer dataset — yielding 101,490 quality-controlled spots across 38 slides. Cell2location deconvolution against a 232,255-cell reference atlas estimated absolute abundances of 22 cell types per spot. A joint compositional-spatial graph clustering approach outperformed all composition-only methods, identifying five biologically interpretable niches: fibromuscular stroma, luminal epithelium, periurethral immune-ductal epithelium, high-grade aggressive tumour, and an incidental cancer niche validated against histopathological annotations. Clinical analysis across 707 PPCG patients revealed a two-phase stromal remodelling model — early smooth muscle loss followed by late cancer-associated fibroblast and macrophage accumulation — with fibroblast proportion independently prognostic beyond Gleason grade. Epigenomic LMC analysis identified desmoplastic and neuro-stromal remodelling programmes as the dominant risk signatures. Most interestingly, the incidental cancer niche was protective in bulk RNA-seq survival yet carried the highest spatial activity for risk LMC programmes, demonstrating that spatial transcriptomics detects localised epigenomic reprogramming diluted and masked in bulk sequencing.

Repository Structure
├── Cell2location_and_downstream/
│   └── downstream/
│       ├── notebooks/          # Jupyter notebooks (numbered by analysis step)
│       └── scripts/            # Python and R scripts for cluster submission
│
└── ST_Erickson_and_JunyiHu_analyses/
                                # Early exploratory analyses of spatial datasets

Note: All data files (.h5ad, .csv, .tsv, .qs, .RDS) and figures (.png, .pdf) are excluded from this repository. Patient-level data from the PPCG cohort is not publicly shareable. Spatial transcriptomics data is available from the original publications (see Data Availability).


Analysis Pipeline
StepNotebook/ScriptDescription0111_explore_spatial_abundances_entropies_batch_effects.ipynbData integration, QC, batch effect assessment0208b_train_c2l_reference.pyCell2location reference model training0309_run_c2l_spatial.pyCell2location spatial mapping0411_2_joint_graph_and_clustermap.ipynbCompositional space construction, joint graph clustering0513_cluster_comparison_leiden_hierarchical.ipynbClustering benchmarking0614_niche_biology.ipynbNiche characterisation: DEGs, identity scoring, z-scores0714_2_niches_H&E.ipynbH&E overlay validation0816_2_PPCG.ipynbPPCG clinical correlations: Kruskal-Wallis, Mann-Whitney0916_4_Celltype_proportions_survivals.ipynbCell type survival analysis: univariate and multivariate Cox1016_5_LMC_proportions_survivals.ipynbLMC survival analysis1116_3_SpatialDEGsvsLMCs.ipynbAUCell LMC scoring in spatial data1216_8_Spatial_scores_survivals.ipynbSpatial niche signature scoring in bulk RNA-seq + survival1314_gsea_lmc.py + 14_GSEA.slurmGSEA pathway enrichment (submitted as cluster job)1416_7_Naming_LMCs.ipynbLMC annotation: GSEA + cell type correlation clustermap1516_6_Liana.ipynbLIANA ligand-receptor interaction analysis

Methods Summary
Data: Erickson et al. 2022 (prostate cancer, 2 patients, 22 Visium slides) + Hu et al. 2025 (normal prostate, 7 patients, 16 slides)
Deconvolution: cell2location v0.1.5 — Bayesian model estimating absolute cell type abundances per spot against a 232,255-cell scRNA-seq reference atlas (22 cell types)
Clustering: Joint compositional-spatial graph — A_joint = α·A_comp + (1−α)·A_spat — benchmarked against Leiden and HC Ward composition-only methods across silhouette, spatial homophily and cross-slide entropy
Clinical cohort: Pan-Prostate Cancer Group (PPCG) — 707 patients, bulk RNA-seq deconvolution (BayesPrim, 14 cell types), DNA methylation array (MeDeCom, 14 LMCs), relapse-free survival (158 events, median follow-up 1,081 days)
Survival analysis: Rank-based inverse normal transformation (Blom's formula) + univariate/multivariate Cox PH (lifelines v0.27) + Kaplan-Meier curves
Pathway enrichment: Preranked GSEA (gseapy) against MSigDB Hallmarks, Reactome 2022, GO Biological Process 2023 and CellMarker 2024
Cell-cell interactions: LIANA v1.7.1, Rank Aggregate consensus method, consensus ligand-receptor database

Environment
bash# Python environment
conda activate thesis
# Python 3.11 — key packages: scanpy 1.11.5, cell2location 0.1.5,
# squidpy 1.6.5, lifelines, gseapy, liana 1.7.1

# R environment (for LMC extraction)
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.4.2
# Key packages: MeDeCom, RnBeads, qs, SummarizedExperiment
All analyses were performed on a high-performance computing cluster using SLURM. GPU-intensive steps (cell2location) required an NVIDIA A100-80GB GPU.
