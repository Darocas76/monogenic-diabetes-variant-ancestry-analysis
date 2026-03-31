#!/usr/bin/env python3
"""
02_query_gnomad_api.py
----------------------
Queries the gnomAD v4.0 public GraphQL API
(https://gnomad.broadinstitute.org/api, dataset gnomad_r4)
for all variants in the 17 monogenic diabetes genes, stratified by
genetic ancestry group.

Output: data/gnomad_mody_raw.csv  (long format, one row per variant×population)

Usage:
    python scripts/02_query_gnomad_api.py

Notes:
    - Uses gnomad_r4 (gnomAD v4.0 genomes).
    - The VariantPopulation object exposes: id, ac, an, homozygote_count.
      There is NO 'af' field in this schema — allele frequency is computed
      as AC/AN where AN > 0.
    - WFS1 has ~5,238 variants and is queried with a longer timeout (300 s).
    - Rate-limiting: a 1-second pause between gene queries is included.
"""

import json
import time
from pathlib import Path

import pandas as pd
import requests

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
DATASET = "gnomad_r4"

MODY_GENES = [
    "HNF1A", "HNF4A", "HNF1B", "GCK", "KCNJ11", "ABCC8", "INS",
    "PDX1", "NEUROD1", "PTF1A", "CEL", "PPARG", "APPL1", "BLK",
    "KLF11", "PAX4", "WFS1",
]

POPULATIONS = ["afr", "amr", "eas", "sas", "nfe", "fin", "asj", "mid"]

QUERY_TEMPLATE = """
{
  gene(gene_symbol: "GENE_SYMBOL", reference_genome: GRCh38) {
    variants(dataset: DATASET_NAME) {
      variant_id
      hgvsc
      hgvsp
      consequence
      lof
      genome {
        ac
        an
        populations {
          id
          ac
          an
          homozygote_count
        }
      }
    }
  }
}
"""

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def query_gene(gene: str, timeout: int = 120) -> list[dict]:
    """Query gnomAD API for a single gene. Returns list of variant dicts."""
    query = QUERY_TEMPLATE.replace("GENE_SYMBOL", gene).replace(
        "DATASET_NAME", DATASET
    )
    payload = json.dumps({"query": query})
    headers = {"Content-Type": "application/json"}

    resp = requests.post(
        GNOMAD_API, data=payload, headers=headers, timeout=timeout
    )
    resp.raise_for_status()
    data = resp.json()

    if "errors" in data:
        raise RuntimeError(f"GraphQL errors for {gene}: {data['errors']}")

    variants = data.get("data", {}).get("gene", {}).get("variants", [])
    return variants


def flatten_variants(gene: str, variants: list[dict]) -> list[dict]:
    """
    Flatten nested variant×population structure to a list of flat dicts.
    Only populations listed in POPULATIONS are kept (avoids sex-specific
    pseudo-populations such as 'XX' / 'XY').
    """
    rows = []
    for v in variants:
        genome = v.get("genome") or {}
        global_ac = genome.get("ac", 0) or 0
        global_an = genome.get("an", 0) or 0
        populations = genome.get("populations") or []

        for pop in populations:
            pop_id = pop["id"].lower()
            if pop_id not in POPULATIONS:
                continue
            rows.append(
                {
                    "gene": gene,
                    "variant_id": v["variant_id"],
                    "hgvsc": v.get("hgvsc", ""),
                    "hgvsp": v.get("hgvsp", ""),
                    "consequence": v.get("consequence", ""),
                    "lof": v.get("lof", ""),
                    "pop": pop_id,
                    "ac": pop["ac"],
                    "an": pop["an"],
                    "global_ac": global_ac,
                    "global_an": global_an,
                }
            )
    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    base = Path(__file__).resolve().parent.parent
    data_dir = base / "data"
    data_dir.mkdir(exist_ok=True)
    out_path = data_dir / "gnomad_mody_raw.csv"

    all_rows: list[dict] = []
    failed: list[str] = []

    for gene in MODY_GENES:
        timeout = 300 if gene == "WFS1" else 120
        print(f"[INFO] Querying {gene} (timeout={timeout}s) …", end=" ", flush=True)
        try:
            variants = query_gene(gene, timeout=timeout)
            rows = flatten_variants(gene, variants)
            all_rows.extend(rows)
            print(f"{len(variants)} variants → {len(rows)} pop rows")
        except Exception as exc:
            print(f"FAILED: {exc}")
            failed.append(gene)
        time.sleep(1)

    if all_rows:
        df = pd.DataFrame(all_rows)
        # Compute per-population allele frequency
        df["af"] = df.apply(
            lambda r: r["ac"] / r["an"] if r["an"] > 0 else 0.0, axis=1
        )
        df.to_csv(out_path, index=False)
        print(
            f"\n[INFO] Saved {out_path} "
            f"({len(df):,} rows, {df['variant_id'].nunique():,} unique variants)."
        )
    else:
        print("[ERROR] No data retrieved.")

    if failed:
        print(f"[WARN] Failed genes (retry manually): {failed}")


if __name__ == "__main__":
    main()
