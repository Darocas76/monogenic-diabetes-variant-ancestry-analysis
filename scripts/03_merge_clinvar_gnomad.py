#!/usr/bin/env python3
"""
03_merge_clinvar_gnomad.py
--------------------------
Cross-references gnomAD v4.0 variants with ClinVar classifications using
GRCh38 VCF coordinates (chr-pos-ref-alt) as the shared key.

Inputs:
    data/clinvar_processed.csv      (from 01_download_clinvar.py)
    data/gnomad_mody_raw.csv        (from 02_query_gnomad_api.py)

Outputs:
    data/gnomad_clinvar_merged.csv  Full cross-referenced dataset
    data/supplementary_table1.csv   ClinVar-annotated variants with AF per
                                    population (wide format; Supp. Table 1)

Usage:
    python scripts/03_merge_clinvar_gnomad.py

Methodology note:
    Each variant is counted independently in each gnomAD population where
    AF > 0 (overlapping/OR-logic methodology). A variant observed in both
    European and non-European populations is counted in both groups. The
    total sum of group-specific observations therefore exceeds the number
    of unique variants; this must be declared in the methods section.
"""

import pickle
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

MODY_GENES = [
    "HNF1A", "HNF4A", "HNF1B", "GCK", "KCNJ11", "ABCC8", "INS",
    "PDX1", "NEUROD1", "PTF1A", "CEL", "PPARG", "APPL1", "BLK",
    "KLF11", "PAX4", "WFS1",
]

POPULATIONS = ["afr", "amr", "eas", "sas", "nfe", "fin", "asj", "mid"]

CLINSIG_MAP = {
    "Pathogenic": "P",
    "Likely pathogenic": "LP",
    "Uncertain significance": "VUS",
    "Likely benign": "LB",
    "Benign": "B",
    "Pathogenic/Likely pathogenic": "LP",
    "Benign/Likely benign": "LB",
}


def map_clinsig(raw: str) -> str:
    if not isinstance(raw, str):
        return "Other"
    for key, cat in CLINSIG_MAP.items():
        if key.lower() in raw.lower():
            return cat
    return "Other"


# ---------------------------------------------------------------------------
# Build ClinVar lookup dict
# ---------------------------------------------------------------------------

def build_clinvar_dict(clinvar_csv: Path) -> dict:
    """
    Read ClinVar processed CSV and build a lookup dict:
        "chr-pos-ref-alt" → ClinicalSignificance string

    Uses PositionVCF / ReferenceAlleleVCF / AlternateAlleleVCF columns
    (not ReferenceAllele/AlternateAllele, which contain 'na' for non-SNVs).
    """
    df = pd.read_csv(clinvar_csv, low_memory=False)

    # Accept either column-name variants
    pos_col = next(
        (c for c in df.columns if "positionvcf" in c.lower() or "PositionVCF" == c),
        None,
    )
    ref_col = next(
        (c for c in df.columns if "referenceallelevcf" in c.lower()),
        None,
    )
    alt_col = next(
        (c for c in df.columns if "alternateallelevcf" in c.lower()),
        None,
    )
    chr_col = next((c for c in df.columns if c.lower() == "chromosome"), None)
    sig_col = next(
        (c for c in df.columns if "clinicalsignificance" in c.lower()), None
    )

    for col, name in [
        (pos_col, "PositionVCF"),
        (ref_col, "ReferenceAlleleVCF"),
        (alt_col, "AlternateAlleleVCF"),
        (chr_col, "Chromosome"),
        (sig_col, "ClinicalSignificance"),
    ]:
        if col is None:
            raise ValueError(f"Column '{name}' not found in {clinvar_csv}.")

    lookup = {}
    for _, row in df.iterrows():
        chrom = str(row[chr_col]).replace("chr", "")
        pos = str(row[pos_col])
        ref = str(row[ref_col])
        alt = str(row[alt_col])
        key = f"{chrom}-{pos}-{ref}-{alt}"
        lookup[key] = str(row[sig_col])

    return lookup


# ---------------------------------------------------------------------------
# Cross-reference gnomAD with ClinVar
# ---------------------------------------------------------------------------

def parse_variant_id(variant_id: str) -> tuple[str, str, str, str]:
    """Parse gnomAD variant ID 'chr-pos-ref-alt' → (chrom, pos, ref, alt)."""
    parts = variant_id.split("-")
    if len(parts) != 4:
        return ("", "", "", "")
    chrom, pos, ref, alt = parts
    chrom = chrom.replace("chr", "")
    return chrom, pos, ref, alt


def merge_datasets(
    gnomad_csv: Path,
    clinvar_dict: dict,
) -> pd.DataFrame:
    """
    Pivot gnomAD long-format data to wide (one row per variant),
    then annotate with ClinVar classification.
    """
    raw = pd.read_csv(gnomad_csv, low_memory=False)

    # Pivot: variant_id → AF per population
    wide = raw.pivot_table(
        index=[
            "gene", "variant_id", "hgvsc", "hgvsp",
            "consequence", "lof", "global_ac", "global_an",
        ],
        columns="pop",
        values="af",
        aggfunc="first",
    ).reset_index()
    wide.columns.name = None

    # Ensure all population columns exist
    for pop in POPULATIONS:
        if pop not in wide.columns:
            wide[pop] = np.nan

    # Compute global AF
    wide["gnomAD_AF_global"] = wide["global_ac"] / wide["global_an"].replace(0, np.nan)

    # Look up ClinVar classification
    def get_clinvar(row):
        chrom, pos, ref, alt = parse_variant_id(row["variant_id"])
        key = f"{chrom}-{pos}-{ref}-{alt}"
        return clinvar_dict.get(key, "Not in ClinVar")

    wide["ClinicalSignificance"] = wide.apply(get_clinvar, axis=1)
    wide["clinvar_category"] = wide["ClinicalSignificance"].apply(map_clinsig)

    return wide


# ---------------------------------------------------------------------------
# Build supplementary_table1
# ---------------------------------------------------------------------------

def build_supplementary_table(merged: pd.DataFrame) -> pd.DataFrame:
    """
    Filter to ClinVar-annotated variants (not 'Not in ClinVar' / 'Other')
    and select columns matching Supplementary Table 1 structure.
    """
    annotated = merged[
        ~merged["clinvar_category"].isin(["Other"])
        & (merged["ClinicalSignificance"] != "Not in ClinVar")
    ].copy()

    col_map = {
        "gene": "gene",
        "variant_id": "variant_id",
        "hgvsc": "hgvsc",
        "hgvsp": "hgvsp",
        "consequence": "consequence",
        "lof": "lof",
        "ClinicalSignificance": "clinical_significance",
        "clinvar_category": "clinvar_category",
        "gnomAD_AF_global": "gnomAD_AF_global",
        "afr": "AF_AFR",
        "amr": "AF_AMR",
        "eas": "AF_EAS",
        "sas": "AF_SAS",
        "nfe": "AF_NFE",
        "fin": "AF_FIN",
        "asj": "AF_ASJ",
        "mid": "AF_MID",
    }
    out_cols = [c for c in col_map if c in annotated.columns]
    supp = annotated[out_cols].rename(columns=col_map)
    return supp.sort_values(["gene", "variant_id"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    base = Path(__file__).resolve().parent.parent
    data_dir = base / "data"
    data_dir.mkdir(exist_ok=True)

    clinvar_csv = data_dir / "clinvar_processed.csv"
    gnomad_csv = data_dir / "gnomad_mody_raw.csv"
    merged_out = data_dir / "gnomad_clinvar_merged.csv"
    supp_out = data_dir / "supplementary_table1.csv"

    print("[INFO] Building ClinVar lookup dictionary …")
    clinvar_dict = build_clinvar_dict(clinvar_csv)
    print(f"[INFO] {len(clinvar_dict):,} ClinVar entries indexed.")

    print("[INFO] Merging gnomAD and ClinVar data …")
    merged = merge_datasets(gnomad_csv, clinvar_dict)
    merged.to_csv(merged_out, index=False)
    print(
        f"[INFO] Saved {merged_out} "
        f"({len(merged):,} unique variants, "
        f"{merged['gnomAD_AF_global'].notna().sum():,} with global AF)."
    )

    print("[INFO] Building supplementary table …")
    supp = build_supplementary_table(merged)
    supp.to_csv(supp_out, index=False)
    match_rate = len(supp) / len(merged) * 100
    print(
        f"[INFO] Saved {supp_out} "
        f"({len(supp):,} annotated variants, {match_rate:.1f}% ClinVar match rate)."
    )


if __name__ == "__main__":
    main()
