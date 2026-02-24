# Junco_reseq

Scripts and lightweight workflows used to analyze resequencing (WGS) data in the *Junco* system.

> This repository is a personal, practical collection of analysis scripts (not a packaged software).  
> Expect a mix of one-off utilities, reproducible pipelines, and plotting helpers.

---

## Contents

This repo is organized by analysis theme:

- **Assembly/** — assembly-related utilities
- **D-statistics/Dsuite/** — ABBA–BABA / f-branch workflows and plotting
- **Demographic/** — SFS preparation, demographic inference helpers (e.g., fastsimcoal2/SMC-style inputs)
- **GO_enrichment/** — gene list formatting and enrichment utilities
- **GWAS/Bypass/** — GWAS workflows and post-processing
- **Genomic_Scans/** — genome scans (FST, Dxy, π, Tajima’s D, iHS, XP-EHH, CLR, etc.)
- **Phylogeny/** — phylogenomic / tree-building helpers
- **Population_structure/** — PCA/ADMIXTURE-like workflows and plotting
- **Structural_variation/** — inversion/SV utilities (e.g., local PCA, SV filtering)
- **Various/** — small general-purpose scripts

---

## Quick start

Clone the repository:

```bash
git clone https://github.com/<your-user>/Junco_reseq.git
cd Junco_reseq
