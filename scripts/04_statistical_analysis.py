#!/usr/bin/env python3
"""
04_statistical_analysis.py
--------------------------
Computes ancestry-stratified classification statistics from the cross-
referenced ClinVar × gnomAD dataset.

Inputs:
    data/supplementary_table1.csv   (from 03_merge_clinvar_gnomad.py)
    data/gnomad_clinvar_merged.csv  (from 03_merge_clinvar_gnomad.py)

Outputs:
    data/gene_by_gene_analysis.csv       Per-gene EUR vs non-EUR statistics
    data/mody_vus_by_population.csv      VUS/P+LP/B+LB% per gnomAD population
    data/mody_vus_EUR_vs_nonEUR.csv      Aggregated EUR vs non-EUR comparison
    data/clingen_results.csv             ClinGen gene-disease validity curations
    data/summary_table.csv              Overall summary statistics

Usage:
    python scripts/04_statistical_analysis.py

Methodology (EUR vs non-EUR aggregation):
    EUR  group: AF_NFE > 0  OR  AF_FIN > 0  OR  AF_ASJ > 0
    Non-EUR:    AF_AFR > 0  OR  AF_AMR > 0  OR  AF_EAS > 0  OR
                AF_SAS > 0  OR  AF_MID > 0

    This OR-logic means a variant can appear in both groups simultaneously.
    The ratio (sum of group totals) / (unique variants) ≈ 2.05 for MODY genes,
    indicating ~38% of variants are observed in both ancestry groups.
    This overlapping methodology must be declared in the methods section.
"""

import requests
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

MODY_GENES = [
    "HNF1A", "HNF4A", "HNF1B", "GCK", "KCNJ11", "ABCC8", "INS",
    "PDX1", "NEUROD1", "PTF1A", "CEL", "PPARG", "APPL1", "BLK",
    "KLF11", "PAX4", "WFS1",
]

# Population AF columns in supplementary_table1
POP_COLS = {
    "AFR": "AF_AFR",
    "AMR": "AF_AMR",
    "EAS": "AF_EAS",
    "SAS": "AF_SAS",
    "NFE": "AF_NFE",
    "FIN": "AF_FIN",
    "ASJ": "AF_ASJ",
    "MID": "AF_MID",
}

EUR_POPS = ["AF_NFE", "AF_FIN", "AF_ASJ"]
NONEUR_POPS = ["AF_AFR", "AF_AMR", "AF_EAS", "AF_SAS", "AF_MID"]


# ---------------------------------------------------------------------------
# Per-population analysis
# ---------------------------------------------------------------------------

def per_population_stats(df: pd.DataFrame) -> pd.DataFrame:
    """Compute VUS/P+LP/B+LB% per individual gnomAD population."""
    rows = []
    for pop_label, col in POP_COLS.items():
        if col not in df.columns:
            continue
        sub = df[df[col].fillna(0) > 0]
        n = len(sub)
        if n == 0:
            continue
        counts = sub["clinvar_category"].value_counts()
        rows.append(
            {
                "Population": pop_label,
                "N_variants": n,
                "VUS_n": counts.get("VUS", 0),
                "VUS_%": counts.get("VUS", 0) / n * 100,
                "PLP_n": counts.get("P", 0) + counts.get("LP", 0),
                "PLP_%": (counts.get("P", 0) + counts.get("LP", 0)) / n * 100,
                "BLB_n": counts.get("B", 0) + counts.get("LB", 0),
                "BLB_%": (counts.get("B", 0) + counts.get("LB", 0)) / n * 100,
            }
        )
    return pd.DataFrame(rows).sort_values("VUS_%", ascending=False)


# ---------------------------------------------------------------------------
# EUR vs non-EUR aggregation
# ---------------------------------------------------------------------------

def eur_noneur_stats(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregated EUR and non-EUR statistics (OR-logic, overlapping)."""
    eur_mask = df[EUR_POPS].fillna(0).gt(0).any(axis=1)
    noneur_mask = df[NONEUR_POPS].fillna(0).gt(0).any(axis=1)

    results = []
    for label, mask in [("EUR", eur_mask), ("Non-EUR", noneur_mask)]:
        sub = df[mask]
        n = len(sub)
        counts = sub["clinvar_category"].value_counts()
        results.append(
            {
                "Group": label,
                "N_total": n,
                "VUS_n": counts.get("VUS", 0),
                "VUS_%": counts.get("VUS", 0) / n * 100 if n else 0,
                "PLP_n": counts.get("P", 0) + counts.get("LP", 0),
                "PLP_%": (counts.get("P", 0) + counts.get("LP", 0)) / n * 100 if n else 0,
                "BLB_n": counts.get("B", 0) + counts.get("LB", 0),
                "BLB_%": (counts.get("B", 0) + counts.get("LB", 0)) / n * 100 if n else 0,
            }
        )

    # Chi-squared test on VUS vs non-VUS
    eur_sub = df[eur_mask]["clinvar_category"]
    noneur_sub = df[noneur_mask]["clinvar_category"]
    contingency = pd.DataFrame(
        {
            "EUR": [(eur_sub == "VUS").sum(), (eur_sub != "VUS").sum()],
            "Non-EUR": [(noneur_sub == "VUS").sum(), (noneur_sub != "VUS").sum()],
        },
        index=["VUS", "non-VUS"],
    )
    chi2, p_val, dof, _ = stats.chi2_contingency(contingency)
    print(
        f"[STATS] EUR vs non-EUR VUS rate: χ²={chi2:.3f}, p={p_val:.4f}, df={dof}"
    )

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Gene-by-gene analysis
# ---------------------------------------------------------------------------

def gene_by_gene(df: pd.DataFrame, merged_csv: Path) -> pd.DataFrame:
    """
    For each gene: ClinVar coverage %, EUR and non-EUR VUS/P+LP/B+LB rates.
    Total variant count comes from gnomad_clinvar_merged.csv.
    """
    merged = pd.read_csv(merged_csv, low_memory=False)
    gnomad_totals = merged.groupby("gene")["variant_id"].nunique().rename("gnomAD_total_vars")

    eur_mask = df[EUR_POPS].fillna(0).gt(0).any(axis=1)
    noneur_mask = df[NONEUR_POPS].fillna(0).gt(0).any(axis=1)

    rows = []
    for gene in MODY_GENES:
        gdf = df[df["gene"] == gene]
        g_total = gnomad_totals.get(gene, len(gdf))
        annotated = len(gdf)
        coverage = annotated / g_total * 100 if g_total > 0 else 0

        def stats_for(mask):
            sub = gdf[mask[gdf.index]]
            n = len(sub)
            if n == 0:
                return 0, 0.0, 0.0, 0.0
            counts = sub["clinvar_category"].value_counts()
            vus_pct = counts.get("VUS", 0) / n * 100
            plp_pct = (counts.get("P", 0) + counts.get("LP", 0)) / n * 100
            blb_pct = (counts.get("B", 0) + counts.get("LB", 0)) / n * 100
            return n, vus_pct, plp_pct, blb_pct

        e_n, e_vus, e_plp, e_blb = stats_for(eur_mask)
        ne_n, ne_vus, ne_plp, ne_blb = stats_for(noneur_mask)

        rows.append(
            {
                "Gene": gene,
                "gnomAD_total_vars": g_total,
                "ClinVar_annotated": annotated,
                "ClinVar_coverage_%": round(coverage, 1),
                "EUR_total_obs": e_n,
                "EUR_VUS_n": int(gdf[eur_mask[gdf.index]]["clinvar_category"].eq("VUS").sum()),
                "EUR_VUS_%": round(e_vus, 2),
                "EUR_PLP_%": round(e_plp, 2),
                "EUR_BLB_%": round(e_blb, 2),
                "nonEUR_total_obs": ne_n,
                "nonEUR_VUS_n": int(gdf[noneur_mask[gdf.index]]["clinvar_category"].eq("VUS").sum()),
                "nonEUR_VUS_%": round(ne_vus, 2),
                "nonEUR_PLP_%": round(ne_plp, 2),
                "nonEUR_BLB_%": round(ne_blb, 2),
                "Delta_VUS_%_EUR_minus_nonEUR": round(e_vus - ne_vus, 2),
                "Delta_PLP_%_EUR_minus_nonEUR": round(e_plp - ne_plp, 2),
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# ClinGen gene-disease validity
# ---------------------------------------------------------------------------

def fetch_clingen() -> pd.DataFrame:
    """Query ClinGen Evidence Repository for gene-disease validity curations."""
    url = "https://erepo.clinicalgenome.org/evrepo/api/classifications"
    genes_query = ",".join(MODY_GENES)
    params = {"genes": genes_query, "limit": 200}
    try:
        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        rows = []
        for item in data.get("data", []):
            gene = item.get("gene", {}).get("symbol", "")
            disease = item.get("disease", {}).get("label", "")
            guidelines = item.get("guidelines", [])
            classification = (
                guidelines[0].get("outcome", {}).get("label", "Unknown")
                if guidelines
                else "Unknown"
            )
            rows.append({"gene": gene, "disease": disease, "classification": classification})
        return pd.DataFrame(rows)
    except Exception as exc:
        print(f"[WARN] ClinGen API query failed: {exc}")
        return pd.DataFrame(columns=["gene", "disease", "classification"])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    base = Path(__file__).resolve().parent.parent
    data_dir = base / "data"

    supp_csv = data_dir / "supplementary_table1.csv"
    merged_csv = data_dir / "gnomad_clinvar_merged.csv"

    print("[INFO] Loading supplementary table …")
    df = pd.read_csv(supp_csv, low_memory=False)
    print(f"[INFO] {len(df):,} annotated variants loaded.")

    # Ensure pop columns are numeric
    for col in list(POP_COLS.values()):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Per-population
    print("[INFO] Computing per-population statistics …")
    pop_df = per_population_stats(df)
    pop_df.to_csv(data_dir / "mody_vus_by_population.csv", index=False)
    print(f"[INFO] Saved mody_vus_by_population.csv")

    # EUR vs non-EUR
    print("[INFO] Computing EUR vs non-EUR aggregated statistics …")
    eur_df = eur_noneur_stats(df)
    eur_df.to_csv(data_dir / "mody_vus_EUR_vs_nonEUR.csv", index=False)
    print(f"[INFO] Saved mody_vus_EUR_vs_nonEUR.csv")
    print(eur_df.to_string(index=False))

    # Gene-by-gene
    if merged_csv.exists():
        print("[INFO] Computing gene-by-gene analysis …")
        g2g = gene_by_gene(df, merged_csv)
        g2g.to_csv(data_dir / "gene_by_gene_analysis.csv", index=False)
        print(f"[INFO] Saved gene_by_gene_analysis.csv")

    # ClinGen
    print("[INFO] Querying ClinGen API …")
    clingen = fetch_clingen()
    clingen.to_csv(data_dir / "clingen_results.csv", index=False)
    print(f"[INFO] Saved clingen_results.csv ({len(clingen)} records).")

    # Summary table
    total = len(df)
    vus_n = (df["clinvar_category"] == "VUS").sum()
    plp_n = df["clinvar_category"].isin(["P", "LP"]).sum()
    blb_n = df["clinvar_category"].isin(["B", "LB"]).sum()
    summary = pd.DataFrame(
        [
            {
                "Total_annotated_variants": total,
                "VUS_n": vus_n,
                "VUS_%": round(vus_n / total * 100, 2),
                "PLP_n": plp_n,
                "PLP_%": round(plp_n / total * 100, 2),
                "BLB_n": blb_n,
                "BLB_%": round(blb_n / total * 100, 2),
                "Genes_analysed": df["gene"].nunique(),
            }
        ]
    )
    summary.to_csv(data_dir / "summary_table.csv", index=False)
    print(f"[INFO] Saved summary_table.csv")


if __name__ == "__main__":
    main()
