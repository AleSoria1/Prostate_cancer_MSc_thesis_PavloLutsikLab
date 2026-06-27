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
