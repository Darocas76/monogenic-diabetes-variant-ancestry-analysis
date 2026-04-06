# Ancestry-stratified variant classification in monogenic diabetes genes
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

Data and analysis scripts for:

**Dario P.** Ancestry-stratified variant classification in monogenic diabetes genes : annotation coverage and differential curation burden. *European Journal of Human Genetics* (submitted, 2026).

**Preprint:** [medRxiv DOI pending]

## Author

**Paulo Dario, PhD**
Instituto Nacional de Saúde Doutor Ricardo Jorge (INSA), Lisboa, Portugal
ORCID: [0000-0002-4203-9179](https://orcid.org/0000-0002-4203-9179)

## Overview

This repository contains the data and code to reproduce the analysis presented in the manuscript. The study cross-references ClinVar variant classifications (GRCh38, April 2026; 4,421,188 variants) with gnomAD v4.0 genome allele frequency data for 17 monogenic diabetes genes, stratified by genetic ancestry group.

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

Python 3.10+, pandas ≥2.1, scipy ≥1.11, matplotlib ≥3.8, requests ≥2.31

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
| `data/supplementary_table1.csv` | 4,366 ClinVar-annotated variants with population allele frequencies (Supplementary Table 1) |
| `data/gene_by_gene_analysis.csv` | ClinVar/gnomAD coverage and VUS rates by gene and ancestry |
| `data/mody_vus_by_population.csv` | VUS rates by individual gnomAD ancestry group |
| `data/mody_vus_EUR_vs_nonEUR.csv` | Aggregated EUR vs non-EUR classification comparison |
| `data/gnomad_clinvar_merged.csv` | Full cross-referenced dataset (ClinVar × gnomAD by population) |
| `data/clinvar_global_summary.csv` | Global ClinVar classification distribution (GRCh38, all chromosomes) |
| `data/clingen_results.csv` | ClinGen gene-disease validity curations for the 17 target genes |

## Figures

| Figure | Description |
|--------|-------------|
| Figure 1 | VUS rate by genetic ancestry group (horizontal bar chart) |
| Figure 2 | Gene × ancestry heatmap of VUS rates (17 genes × 8 populations) |
| Figure 3 | ClinVar annotation coverage vs. EUR–non-EUR VUS divergence (scatter plot) |

## License

This work is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## Citation

> Dario P. Ancestry-stratified variant classification in monogenic diabetes genes : annotation coverage and differential curation burden. *Eur J Hum Genet*. 2026 (submitted).

[To be updated upon publication]
