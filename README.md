# Ancestry-stratified variant classification in monogenic diabetes genes

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

Analysis scripts for the paper: "Ancestry-stratified analysis of variant classification in monogenic diabetes genes reveals a 70% annotation gap and differential curation burden across genomic databases"

**Author:** Paulo Dario, PhD  
**Affiliation:** Instituto Nacional de Saúde Doutor Ricardo Jorge (INSA), Lisboa, Portugal  
**ORCID:** 0000-0002-4203-9179

## Contents

- `gene_by_gene_analysis.csv` — ClinVar/gnomAD coverage and VUS rates for 17 MODY genes by ancestry group
- `mody_vus_by_population.csv` — VUS rates by individual gnomAD ancestry group
- `supplementary_table1.csv` — 4,365 classified variants with population allele frequencies
- `figures/` — Publication figures (PNG and SVG, 300 DPI)
- Analysis scripts for ClinVar download, gnomAD API query, and statistical analysis

## Data sources

- ClinVar variant_summary GRCh38 (March 2026): https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/
- gnomAD v4.0 genomes GraphQL API: https://gnomad.broadinstitute.org/api (dataset: gnomad_r4)

## Requirements

Python 3.10+, pandas 2.1, scipy 1.11, matplotlib 3.8, requests

## Citation

[To be updated upon publication]
