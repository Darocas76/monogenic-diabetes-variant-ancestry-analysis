# Ancestry-stratified variant classification in monogenic diabetes genes

[![DOI](https://img.shields.io/badge/DOI-10.1111%2Fahg.70054-blue.svg)](https://doi.org/10.1111/ahg.70054)

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![medRxiv](https://img.shields.io/badge/medRxiv-2026.04.06.26350230-red.svg)](https://doi.org/10.64898/2026.04.06.26350230)

Data and analysis code for:

**Dario P.** Ancestry-Stratified Variant Classification in Monogenic Diabetes Genes: Annotation Coverage and Differential Curation Burden. *Annals of Human Genetics*, published online 23 August 2026. https://doi.org/10.1111/ahg.70054

**Preprint (earlier version):** [medRxiv 10.64898/2026.04.06.26350230](https://doi.org/10.64898/2026.04.06.26350230) — reports the earlier genomes-only analysis; please cite the published article above.

## Author

**Paulo Dario, PhD**
Departamento da Promocao da Saude e Prevencao de Doencas Nao Transmissiveis, Instituto Nacional de Saude Doutor Ricardo Jorge (INSA), Lisboa, Portugal
Centro Cardiovascular da Universidade de Lisboa (CCUL), Faculdade de Medicina, Universidade de Lisboa
BioSystems & Integrative Sciences Institute (BioISI), Faculdade de Ciencias, Universidade de Lisboa
ORCID: [0000-0002-4203-9179](https://orcid.org/0000-0002-4203-9179)
Correspondence: paulo.dario@insa.min-saude.pt

## Overview

This repository reproduces the analysis in the manuscript. It cross-references ClinVar clinical classifications (GRCh38 VCF, accessed 2026-06-29) with gnomAD v4.0 allele-count data (genomes and exomes combined; 807,162 individuals) for 16 primary monogenic diabetes (MODY) genes, stratified by genetic-ancestry group. Variants are classified as **population-private** when observed (allele count > 0) in exactly one ancestry macro-group — European {NFE, FIN, ASJ} or non-European {AFR, AMR, EAS, SAS, MID} — with variants seen in both reported as *shared*. *KLF11*, refuted as a MODY gene by the ClinGen Gene Curation Expert Panel, is excluded from the primary analysis (16 genes) and retained only in the supplement.

Gene set (17 total; 16 primary): *HNF1A, HNF4A, HNF1B, GCK, KCNJ11, ABCC8, INS, PDX1, NEUROD1, PTF1A, CEL, PPARG, APPL1, BLK, PAX4, WFS1* (primary) and *KLF11* (refuted; supplement only).

## Key findings

- **Annotation gap (near-universal, not ancestry-specific):** of 54,865 gnomAD v4 variants in the 16 primary genes, 86.7% have no ClinVar classification (89.4% European-private, 87.5% non-European-private, 62.8% shared). In absolute terms the unannotated set holds more European-private (17,608) than non-European-private (9,860) variants, reflecting gnomAD's European-weighted composition.
- **Actionability gap (the central result):** among classified population-private variants, the pathogenic/likely-pathogenic rate is 19.8% for European-private (n = 1,935) but only 8.6% for non-European-private variants (n = 1,274); Fisher exact OR 2.6 (95% CI 2.1-3.3), p = 1x10^-18, with a higher uncertain-significance rate in non-European-private variants (52.1% vs 45.8%).
- **Global ClinVar composition (GRCh38, release 2026-06-28):** 4,439,382 records; VUS 52.1%, likely benign 24.5%, pathogenic + likely pathogenic 6.8%.

## Repository structure

```
scripts/          Analysis pipeline (Python)
data/             Frozen inputs and merged/derived data (JSON + CSV)
figures/          Manuscript figures (PNG, 300 DPI)
supplementary/    Supplementary tables (XLSX)
RESULTS.md        Frozen headline numbers for this revision
requirements.txt  Python dependencies
```

## Data sources

- **ClinVar:** GRCh38 VCF (NCBI), accessed 2026-06-29 — https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
- **gnomAD:** v4.0 genomes + exomes via the public GraphQL API — https://gnomad.broadinstitute.org/api (dataset: `gnomad_r4`)

A single dated snapshot is used throughout. No filtering by ClinVar review status is applied; records with conflicting interpretations are mapped to *other* and excluded from classification-rate calculations.

## Reproducing the analysis

```bash
pip install -r requirements.txt
```

Pipeline (run from `data/`, which holds the frozen inputs):

0. **Build ClinVar lookup + global composition (Table 1):** `python ../scripts/00_clinvar_vcf_to_lookup.py clinvar.vcf.gz` -> `clinvar_lookup.json` + `table1_global.csv` (needs the ClinVar GRCh38 VCF; gzip only, no pysam)
1. **Fetch gnomAD v4 (genomes + exomes):** `python ../scripts/01_fetch_gnomad_v4.py` -> `gnomad_v4_all.json`
2. **Cross-reference and analyse:** `python ../scripts/02_merge_analyze.py` -> `merged_v4.json` + printed headline numbers
3. **Generate figures:** `python ../scripts/03_make_figures.py` -> `Figure_1.png`, `Figure_2.png`, `Figure_3.png`

Steps 2-3 run from the pre-computed JSON files in `data/`, so the analysis can be reproduced without re-querying the gnomAD API. Rebuilding `clinvar_lookup.json` from the ClinVar VCF requires `pysam`.

## Figures

- **Figure 1 - Annotation gap by ancestry.** Percentage of gnomAD v4 variants without any ClinVar classification, overall and by ancestry group.
- **Figure 2 - Classification of population-private variants by ancestry.** Grouped bars (VUS; P/LP; B/LB) for European-private vs non-European-private variants.
- **Figure 3 - Gene-level P/LP rate among classified population-private variants.** Heatmap by gene and ancestry.

## License

Content licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
