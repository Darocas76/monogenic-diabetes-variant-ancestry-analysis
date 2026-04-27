# Ancestry-stratified variant classification in monogenic diabetes genes
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![medRxiv](https://img.shields.io/badge/medRxiv-2026.04.06.26350230-red.svg)](https://doi.org/10.64898/2026.04.06.26350230)

Data and analysis scripts for:

**Dario P.** Ancestry-stratified variant classification in monogenic diabetes genes: annotation coverage and differential curation burden. *Genetics in Medicine* (submitted, 2026).

**Preprint:** [medRxiv DOI 10.64898/2026.04.06.26350230](https://doi.org/10.64898/2026.04.06.26350230) (CC BY 4.0)

## Author

**Paulo Dario, PhD**
Instituto Nacional de Saúde Doutor Ricardo Jorge (INSA), Lisboa, Portugal
Centro Cardiovascular da Universidade de Lisboa (CCUL), Faculdade de Medicina, Universidade de Lisboa
BioSystems & Integrative Sciences Institute (BioISI), Faculdade de Ciências, Universidade de Lisboa
ORCID: [0000-0002-4203-9179](https://orcid.org/0000-0002-4203-9179)
Correspondence: paulo.dario@insa.min-saude.pt

## Overview

This repository contains the data and code to reproduce the analysis presented in the manuscript. The study cross-references ClinVar variant classifications (GRCh38, April 2026; 4,421,188 variants) with gnomAD v4.0 genome allele frequency data for 17 monogenic diabetes genes (*HNF1A, HNF4A, HNF1B, GCK, KCNJ11, ABCC8, INS, PDX1, NEUROD1, PTF1A, CEL, PPARG, APPL1, BLK, KLF11, PAX4, WFS1*), stratified by genetic ancestry group.

**Key findings:** 70.3% annotation gap (10,325 of 14,691 gnomAD variants without ClinVar classification); divergent mechanisms producing apparent VUS-rate symmetry between EUR and non-EUR groups (NFE submission backlog vs AFR functional-evidence deficit); pattern inversion at *GCK* (non-EUR VUS 18.5% > EUR 15.0%) consistent with progressive European reclassification absent in non-European cohorts.

## Repository structure

```
scripts/          Analysis pipeline (Python)
data/             Processed data files (CSV)
figures/          Publication figures (PNG 300 DPI + SVG)
```

## Data sources

- **ClinVar:** variant_summary.txt.gz (GRCh38, accessed April 2026) — https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/
- **gnomAD:** v4.0 genomes via public GraphQL API — https://gnomad.broadinstitute.org/api (dataset: gnomad_r4)

## Reproducing the analysis

### Requirements

```bash
pip install -r requirements.txt
```

Python 3.10+, pandas ≥2.1, scipy ≥1.11, matplotlib ≥3.8, requests ≥2.31, adjustText ≥1.0

### Pipeline

1. **Download ClinVar data:** `python scripts/01_download_clinvar.py`
2. **Query gnomAD API:** `python scripts/02_query_gnomad_api.py`
3. **Cross-reference databases:** `python scripts/03_merge_clinvar_gnomad.py`
4. **Statistical analysis:** `python scripts/04_statistical_analysis.py`
5. **Generate figures:** `python scripts/05_generate_figures.py`

> **Note:** Step 2 queries the gnomAD GraphQL API for 17 genes and may take several minutes. Steps 3–5 use the pre-computed data files in `data/` and can be run independently.

## Key data files

| File | Description |
|------|-------------|
| `data/supplementary_table1.csv` | 4,366 ClinVar-annotated variants with population allele frequencies (matches Supplementary Table S3 in the manuscript) |
| `data/gene_by_gene_analysis.csv` | ClinVar/gnomAD coverage and VUS rates by gene and ancestry |
| `data/mody_vus_by_population.csv` | VUS rates by individual gnomAD ancestry group |
| `data/mody_vus_EUR_vs_nonEUR.csv` | Aggregated EUR vs non-EUR classification comparison |
| `data/gnomad_clinvar_merged.csv` | Full cross-referenced dataset (ClinVar × gnomAD by population) |
| `data/clinvar_global_summary.csv` | Global ClinVar classification distribution (GRCh38, all chromosomes) |
| `data/clingen_results.csv` | ClinGen gene-disease validity curations for the 17 target genes |

## Figures

Numbering matches the manuscript:

| Figure | Description |
|--------|-------------|
| Figure 1 | ClinVar annotation coverage vs. EUR–non-EUR VUS divergence (scatter plot, 17 genes; gene labels positioned with adjustText) |
| Figure 2 | VUS rate by genetic ancestry group (horizontal bar chart, 8 gnomAD groups) |
| Figure 3 | Gene × ancestry heatmap of VUS rates (17 genes × 8 populations) |

## Supplementary materials

Supplementary tables and the VariantValidator HGVS verification report are included with the journal submission and in the medRxiv preprint deposit.

## License

This work is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). You are free to share and adapt the material with appropriate attribution.

## Citation

> Dario P. Ancestry-stratified variant classification in monogenic diabetes genes: annotation coverage and differential curation burden. medRxiv 2026.04.06.26350230 (preprint). doi:10.64898/2026.04.06.26350230

(Citation will be updated with the published journal reference upon acceptance.)

## Issues

Please open an issue on this repository if you find a bug, have a question about reproducing the analysis, or want to discuss extensions to other gene panels.
