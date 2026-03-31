#!/usr/bin/env python3
"""
01_download_clinvar.py
----------------------
Downloads ClinVar variant_summary.txt.gz (GRCh38) and extracts variants
for the 17 monogenic diabetes genes of interest.

Output: data/clinvar_processed.csv

Usage:
    python scripts/01_download_clinvar.py

Notes:
    - The full ClinVar file is ~400 MB; extraction is done in chunks to
      minimise memory usage.
    - Only GRCh38 (assembly == 'GRCh38') rows are retained.
    - Clinical significance is mapped to five categories:
        P (Pathogenic), LP (Likely pathogenic), VUS (Uncertain significance),
        LB (Likely benign), B (Benign).
"""

import gzip
import os
import urllib.request
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

CLINVAR_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/"
    "variant_summary.txt.gz"
)

MODY_GENES = [
    "HNF1A", "HNF4A", "HNF1B", "GCK", "KCNJ11", "ABCC8", "INS",
    "PDX1", "NEUROD1", "PTF1A", "CEL", "PPARG", "APPL1", "BLK",
    "KLF11", "PAX4", "WFS1",
]

CLINSIG_MAP = {
    "Pathogenic": "P",
    "Likely pathogenic": "LP",
    "Uncertain significance": "VUS",
    "Likely benign": "LB",
    "Benign": "B",
    "Pathogenic/Likely pathogenic": "LP",
    "Benign/Likely benign": "LB",
}

CHUNK_SIZE = 150_000

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def download_clinvar(dest: Path) -> None:
    """Download ClinVar variant_summary.txt.gz if not already present."""
    if dest.exists():
        print(f"[INFO] {dest.name} already present, skipping download.")
        return
    print(f"[INFO] Downloading ClinVar data to {dest} …")
    urllib.request.urlretrieve(CLINVAR_URL, dest)
    print(f"[INFO] Download complete ({dest.stat().st_size / 1e6:.1f} MB).")


def map_clinsig(raw: str) -> str:
    """Return simplified category for a raw ClinicalSignificance string."""
    for key, cat in CLINSIG_MAP.items():
        if key.lower() in raw.lower():
            return cat
    return "Other"


def extract_mody_variants(gz_path: Path, genes: list[str]) -> pd.DataFrame:
    """Iterate over ClinVar in chunks and keep GRCh38 MODY gene rows."""
    gene_set = set(genes)
    kept_chunks = []

    with gzip.open(gz_path, "rt", encoding="utf-8", errors="replace") as fh:
        for chunk in pd.read_csv(
            fh,
            sep="\t",
            chunksize=CHUNK_SIZE,
            low_memory=False,
        ):
            sub = chunk[
                (chunk["Assembly"] == "GRCh38")
                & (chunk["GeneSymbol"].isin(gene_set))
            ]
            if not sub.empty:
                kept_chunks.append(sub)

    if not kept_chunks:
        raise ValueError("No MODY gene variants found — check column names.")

    df = pd.concat(kept_chunks, ignore_index=True)
    return df


def build_processed_table(df: pd.DataFrame) -> pd.DataFrame:
    """Select relevant columns and add clinvar_category."""
    cols = [
        "GeneSymbol", "Name", "Type",
        "ClinicalSignificance", "ReviewStatus",
        "Chromosome", "PositionVCF",
        "ReferenceAlleleVCF", "AlternateAlleleVCF",
        "RS# (dbSNP)", "VariationID",
    ]
    # Keep only columns that exist
    cols = [c for c in cols if c in df.columns]
    out = df[cols].copy()
    out.columns = [c.replace(" ", "_").replace("#", "num") for c in out.columns]
    out["clinvar_category"] = out["ClinicalSignificance"].apply(map_clinsig)
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    base = Path(__file__).resolve().parent.parent
    data_dir = base / "data"
    data_dir.mkdir(exist_ok=True)

    gz_path = base / "variant_summary.txt.gz"
    out_path = data_dir / "clinvar_processed.csv"

    download_clinvar(gz_path)
    print("[INFO] Extracting MODY gene variants …")
    df = extract_mody_variants(gz_path, MODY_GENES)
    print(f"[INFO] Found {len(df):,} GRCh38 rows for {df['GeneSymbol'].nunique()} genes.")

    processed = build_processed_table(df)
    processed.to_csv(out_path, index=False)
    print(f"[INFO] Saved {out_path} ({len(processed):,} rows).")

    # Global ClinVar summary (chromosome breakdown)
    global_summary_path = data_dir / "clinvar_global_summary.csv"
    if not global_summary_path.exists():
        print("[INFO] Building global ClinVar chromosome summary …")
        gz_path2 = base / "variant_summary.txt.gz"
        chroms = []
        with gzip.open(gz_path2, "rt", encoding="utf-8", errors="replace") as fh:
            for chunk in pd.read_csv(fh, sep="\t", chunksize=CHUNK_SIZE, low_memory=False):
                sub = chunk[chunk["Assembly"] == "GRCh38"]
                if not sub.empty:
                    chroms.append(sub[["Chromosome", "ClinicalSignificance"]])
        all_data = pd.concat(chroms, ignore_index=True)
        all_data["cat"] = all_data["ClinicalSignificance"].apply(map_clinsig)
        summary = (
            all_data.groupby("Chromosome")["cat"]
            .value_counts()
            .unstack(fill_value=0)
            .reset_index()
        )
        summary.to_csv(global_summary_path, index=False)
        print(f"[INFO] Saved {global_summary_path}.")


if __name__ == "__main__":
    main()
