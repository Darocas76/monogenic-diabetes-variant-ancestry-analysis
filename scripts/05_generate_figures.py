#!/usr/bin/env python3
"""
05_generate_figures.py
----------------------
Generates the three publication figures for:
    "Ancestry-stratified variant classification in monogenic diabetes genes"

Inputs:
    data/mody_vus_by_population.csv   (from 04_statistical_analysis.py)
    data/supplementary_table1.csv     (from 03_merge_clinvar_gnomad.py)
    data/gene_by_gene_analysis.csv    (from 04_statistical_analysis.py)

Outputs (PNG 300 DPI + SVG):
    figures/Figure_1.png / .svg   — VUS rate by genetic ancestry group (bar chart)
    figures/Figure_2.png / .svg   — Gene × ancestry VUS heatmap
    figures/Figure_3.png / .svg   — ClinVar coverage vs EUR–non-EUR VUS divergence

Usage:
    python scripts/05_generate_figures.py
"""

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

matplotlib.use("Agg")  # non-interactive backend

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------

POP_FULL_NAMES = {
    "AFR": "African/African American",
    "AMR": "Latino/Admixed American",
    "EAS": "East Asian",
    "SAS": "South Asian",
    "MID": "Middle Eastern",
    "NFE": "Non-Finnish European",
    "FIN": "Finnish",
    "ASJ": "Ashkenazi Jewish",
}

EUR_POPS = {"NFE", "FIN", "ASJ"}
COLORS = {
    "EUR": "#2166ac",
    "Non-EUR": "#d73027",
    "neutral": "#4d4d4d",
}


def setup_style():
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "figure.dpi": 100,
            "savefig.dpi": 300,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


# ---------------------------------------------------------------------------
# Figure 1 — horizontal bar chart, VUS% by population
# ---------------------------------------------------------------------------

def figure1(pop_df: pd.DataFrame, out_dir: Path):
    """Horizontal bar chart: VUS rate by gnomAD ancestry group."""
    pop_df = pop_df.sort_values("VUS_%", ascending=True).copy()
    pop_df["label"] = pop_df["Population"].map(
        lambda p: POP_FULL_NAMES.get(p, p)
    )
    colors = [
        COLORS["EUR"] if p in EUR_POPS else COLORS["Non-EUR"]
        for p in pop_df["Population"]
    ]

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.barh(pop_df["label"], pop_df["VUS_%"], color=colors, edgecolor="white")

    # Value labels
    for bar, val in zip(bars, pop_df["VUS_%"]):
        ax.text(
            bar.get_width() + 0.3, bar.get_y() + bar.get_height() / 2,
            f"{val:.1f}%", va="center", ha="left", fontsize=9,
        )

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=COLORS["EUR"], label="European ancestry"),
        Patch(facecolor=COLORS["Non-EUR"], label="Non-European ancestry"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", frameon=False)

    ax.set_xlabel("Variants of Uncertain Significance (%)")
    ax.set_title(
        "VUS rate by genetic ancestry group\n"
        "17 monogenic diabetes genes, ClinVar × gnomAD v4.0",
        pad=10,
    )
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
    ax.set_xlim(0, max(pop_df["VUS_%"]) * 1.12)

    plt.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"Figure_1.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("[INFO] Figure_1 saved.")


# ---------------------------------------------------------------------------
# Figure 2 — heatmap, gene × ancestry VUS%
# ---------------------------------------------------------------------------

def figure2(supp_df: pd.DataFrame, out_dir: Path):
    """Heatmap of VUS rates: 17 genes × 8 gnomAD ancestry groups."""
    POP_COLS = {
        "AFR": "AF_AFR", "AMR": "AF_AMR", "EAS": "AF_EAS", "SAS": "AF_SAS",
        "NFE": "AF_NFE", "FIN": "AF_FIN", "ASJ": "AF_ASJ", "MID": "AF_MID",
    }
    MODY_GENES = [
        "HNF1A", "HNF4A", "HNF1B", "GCK", "KCNJ11", "ABCC8", "INS",
        "PDX1", "NEUROD1", "PTF1A", "CEL", "PPARG", "APPL1", "BLK",
        "KLF11", "PAX4", "WFS1",
    ]

    matrix = pd.DataFrame(index=MODY_GENES, columns=list(POP_COLS.keys()), dtype=float)

    for pop_label, af_col in POP_COLS.items():
        if af_col not in supp_df.columns:
            matrix[pop_label] = np.nan
            continue
        for gene in MODY_GENES:
            gdf = supp_df[supp_df["gene"] == gene]
            sub = gdf[gdf[af_col].fillna(0) > 0]
            n = len(sub)
            if n == 0:
                matrix.loc[gene, pop_label] = np.nan
            else:
                matrix.loc[gene, pop_label] = (sub["clinvar_category"] == "VUS").sum() / n * 100

    fig, ax = plt.subplots(figsize=(10, 7))
    data = matrix.values.astype(float)

    im = ax.imshow(data, cmap="YlOrRd", aspect="auto", vmin=0, vmax=100)

    # Axes labels
    ax.set_xticks(range(len(POP_COLS)))
    ax.set_xticklabels(list(POP_COLS.keys()), rotation=45, ha="right")
    ax.set_yticks(range(len(MODY_GENES)))
    ax.set_yticklabels(MODY_GENES)

    # Annotate cells
    for i in range(len(MODY_GENES)):
        for j in range(len(POP_COLS)):
            val = data[i, j]
            if not np.isnan(val):
                color = "white" if val > 65 else "black"
                ax.text(j, i, f"{val:.0f}", ha="center", va="center",
                        fontsize=7.5, color=color)
            else:
                ax.text(j, i, "—", ha="center", va="center", fontsize=7, color="#aaa")

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("VUS (%)")

    ax.set_title(
        "VUS rate by gene and genetic ancestry group\n"
        "ClinVar × gnomAD v4.0",
        pad=10,
    )

    plt.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"Figure_2.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("[INFO] Figure_2 saved.")


# ---------------------------------------------------------------------------
# Figure 3 — scatter, coverage vs EUR–non-EUR delta
# ---------------------------------------------------------------------------

def figure3(g2g_df: pd.DataFrame, out_dir: Path):
    """Scatter: ClinVar annotation coverage vs EUR–non-EUR VUS% divergence."""
    fig, ax = plt.subplots(figsize=(8, 6))

    x = g2g_df["ClinVar_coverage_%"]
    y = g2g_df["Delta_VUS_%_EUR_minus_nonEUR"]

    # Colour by delta direction
    colors = [COLORS["EUR"] if v >= 0 else COLORS["Non-EUR"] for v in y]

    ax.scatter(x, y, c=colors, s=80, edgecolors="white", linewidths=0.5, zorder=3)
    ax.axhline(0, color="#888", linewidth=0.8, linestyle="--", zorder=1)

    # Gene labels (offset to avoid overlap)
    for _, row in g2g_df.iterrows():
        ax.annotate(
            row["Gene"],
            (row["ClinVar_coverage_%"], row["Delta_VUS_%_EUR_minus_nonEUR"]),
            textcoords="offset points",
            xytext=(5, 3),
            fontsize=8,
            color="#333",
        )

    # Reference line (linear trend)
    if len(x) > 2:
        z = np.polyfit(x.fillna(0), y.fillna(0), 1)
        xfit = np.linspace(x.min(), x.max(), 100)
        ax.plot(xfit, np.polyval(z, xfit), color="#aaa", linewidth=1, linestyle="--", zorder=2)

    ax.set_xlabel("ClinVar annotation coverage (%)")
    ax.set_ylabel("ΔVUS%  (EUR − non-EUR)")
    ax.set_title(
        "ClinVar annotation coverage vs. EUR – non-EUR VUS divergence\n"
        "17 monogenic diabetes genes",
        pad=10,
    )

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=COLORS["EUR"], label="Higher VUS in EUR"),
        Patch(facecolor=COLORS["Non-EUR"], label="Higher VUS in non-EUR"),
    ]
    ax.legend(handles=legend_elements, frameon=False)

    plt.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"Figure_3.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("[INFO] Figure_3 saved.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    base = Path(__file__).resolve().parent.parent
    data_dir = base / "data"
    fig_dir = base / "figures"
    fig_dir.mkdir(exist_ok=True)

    setup_style()

    pop_csv = data_dir / "mody_vus_by_population.csv"
    supp_csv = data_dir / "supplementary_table1.csv"
    g2g_csv = data_dir / "gene_by_gene_analysis.csv"

    print("[INFO] Loading data …")
    pop_df = pd.read_csv(pop_csv)
    supp_df = pd.read_csv(supp_csv, low_memory=False)
    g2g_df = pd.read_csv(g2g_csv)

    figure1(pop_df, fig_dir)
    figure2(supp_df, fig_dir)
    figure3(g2g_df, fig_dir)

    print("[INFO] All figures saved to figures/ (PNG 300 DPI + SVG).")


if __name__ == "__main__":
    main()
