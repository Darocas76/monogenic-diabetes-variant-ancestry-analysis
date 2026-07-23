# AHG Major Revision — reanalysis results (frozen snapshot)

Data freeze: gnomAD v4 (genomes+exomes, dataset gnomad_r4) + ClinVar GRCh38 VCF accessed 2026-06-29.
Primary analysis = 16 genes (KLF11 excluded; refuted by ClinGen GCEP, kept in supplementary).
"Population-private" = present (AC>0) in exactly one macro-group. EUR={NFE,FIN,ASJ}; non-EUR={AFR,AMR,EAS,SAS,MID}.
Independently reproduced by a separate agent (exact match on classified counts; gap strata ±<0.05%).

## Headline numbers (16 genes)
- gnomAD v4 variants: 54,865 (vs 14,691 genomes-only; ~3.9x). 17 genes: 57,168.
- Annotation gap: 86.7% unannotated (was 70.3% genomes-only).
- Gap-rate by ancestry (near-universal, NOT disproportionately non-European):
  EUR-private 89.4% | non-EUR-private 87.5% | shared 62.8%.
  Absolute unannotated counts: EUR-private 17,608 > non-EUR-private 9,860 (gnomAD is European-heavy).
- Actionability among classified population-private variants (THE finding):
  EUR-private (n=1,935): VUS 45.8% | P/LP 19.8% | B/LB 34.3%
  non-EUR-private (n=1,274): VUS 52.1% | P/LP 8.6% | B/LB 39.2%
  Fisher exact (P/LP vs rest, EUR vs non-EUR): OR 2.6, p≈1e-18.

## Files
- data/table2_coverage.csv — per-gene gnomAD total, ClinVar-annotated, coverage% (KLF11 flagged).
- data/table_gene_ancestry.csv — per-gene gap-rate and classified VUS/P/LP by ancestry.
- data/{gnomad_v4_all,clinvar_lookup,merged_v4}.json — frozen inputs/merged.
- figures/Figure_1.png (annotation gap by ancestry), Figure_2.png (classification of population-private variants), Figure_3.png (per-gene P/LP heatmap).
- scripts/00_clinvar_vcf_to_lookup.py, 01_fetch_gnomad_v4.py, 02_merge_analyze.py, 03_make_figures.py.

## Reviewer references verified (PubMed)
- Sharp LN et al. 2026, J Clin Endocrinol Metab — genotype-first MODY in 454,275 UK Biobank exomes. PMID 41175096; DOI 10.1210/clinem/dgaf599.
- Elashi AA et al. 2022, Int J Mol Sci — MODY genetic spectrum in Qatar (QBB WGS 14,364). PMID 36613572; DOI 10.3390/ijms24010130.
