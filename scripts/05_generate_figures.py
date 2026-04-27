#!/usr/bin/env python3
"""
05_generate_figures.py
----------------------
Generates the three publication figures for:
    "Ancestry-stratified variant classification in monogenic diabetes genes:
     annotation coverage and differential curation burden"

Figure numbering matches the manuscript:
    Figure 1 — Scatter: ClinVar annotation coverage vs EUR-non-EUR VUS divergence
    Figure 2 — Bar chart: VUS rate by genetic ancestry group
    Figure 3 — Heatmap: gene x ancestry VUS rate

Inputs:
    data/mody_vus_by_population.csv   (from 04_statistical_analysis.py)
    data/supplementary_table1.csv     (from 03_merge_clinvar_gnomad.py)
    data/gene_by_gene_analysis.csv    (from 04_statistical_analysis.py)

Outputs (PNG 300 DPI + SVG):
    figures/Figure_1.png / .svg   - scatter
    figures/Figure_2.png / .svg   - bar chart
    figures/Figure_3.png / .svg   - heatmap

Usage:
    python scripts/05_generate_figures.py
"""

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

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
    "scatter_EUR": "#c0392b",   # red, used for delta>=0 in scatter
    "scatter_Non": "#2980b9",   # blue, used for delta<0 in scatter
    "neutral": "#4d4d4d",
    "ci_band": "#bbbbbb",
}


def setup_style():
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "figure.dpi": 100,
            "savefig.dpi": 300,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


# ---------------------------------------------------------------------------
# Figure 1 - scatter: coverage vs EUR-non-EUR delta (with adjustText for labels)
# ---------------------------------------------------------------------------

def figure1(g2g_df: pd.DataFrame, out_dir: Path):
    """Scatter: ClinVar annotation coverage vs EUR-non-EUR VUS% divergence.

    Uses adjustText to prevent gene labels from overlapping their markers.
    """
    try:
        from adjustText import adjust_text
    except ImportError:
        adjust_text = None
        print("[WARN] adjustText not available; labels may overlap markers. "
              "Install with: pip install adjustText")

    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    x = g2g_df["ClinVar_coverage_%"].values
    y = g2g_df["Delta_VUS_%_EUR_minus_nonEUR"].values
    labels = g2g_df["Gene"].values
    colors = [COLORS["scatter_EUR"] if v >= 0 else COLORS["scatter_Non"] for v in y]

    fig, ax = plt.subplots(figsize=(10, 7.2))

    # Bootstrap regression CI
    rng = np.random.default_rng(42)
    n_boot = 5000
    n = len(x)
    xs = np.linspace(x.min(), x.max(), 100)
    boot_lines = np.zeros((n_boot, len(xs)))
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        z = np.polyfit(x[idx], y[idx], 1)
        boot_lines[i] = np.polyval(z, xs)
    ci_lower = np.percentile(boot_lines, 2.5, axis=0)
    ci_upper = np.percentile(boot_lines, 97.5, axis=0)
    r, p = pearsonr(x, y)

    ax.fill_between(xs, ci_lower, ci_upper, color=COLORS["ci_band"], alpha=0.4, zorder=1)
    z = np.polyfit(x, y, 1)
    ax.plot(xs, np.polyval(z, xs), color=COLORS["neutral"], linewidth=1.5, zorder=2)
    ax.axhline(0, color="#888", linewidth=0.7, linestyle="--", zorder=1)

    ax.scatter(x, y, c=colors, s=120, edgecolors="white",
               linewidths=1.2, zorder=4, alpha=0.92)

    texts = [ax.text(xi, yi, lab, fontsize=10, fontstyle="italic",
                     fontweight="500", color="#222", zorder=5)
             for xi, yi, lab in zip(x, y, labels)]

    if adjust_text is not None:
        adjust_text(
            texts, x=x, y=y, ax=ax,
            arrowprops=dict(arrowstyle="-", color="#666", lw=0.7, alpha=0.85),
            expand=(1.4, 1.6),
            force_text=(0.6, 0.8),
            force_static=(0.4, 0.6),
            only_move={'text': 'xy', 'static': 'xy', 'explode': 'xy', 'pull': 'xy'},
            min_arrow_len=4,
            avoid_self=True,
        )

    ax.text(0.02, 0.97, "EUR VUS% higher", transform=ax.transAxes,
            fontsize=10, color=COLORS["scatter_EUR"], va="top", ha="left", alpha=0.85)
    ax.text(0.02, 0.03, "Non-EUR VUS% higher", transform=ax.transAxes,
            fontsize=10, color=COLORS["scatter_Non"], va="bottom", ha="left", alpha=0.85)

    ax.set_xlabel("ClinVar annotation coverage (%)\n[annotated variants / gnomAD variants x 100]")
    ax.set_ylabel(u"ΔVUS%  (European - Non-European)")
    ax.set_title("ClinVar annotation coverage vs. " + u"Δ" + "VUS% (EUR - Non-EUR)\nacross 17 monogenic diabetes genes",
                 pad=12, fontweight="bold")

    legend_elements = [
        Line2D([0], [0], marker='o', linestyle='', color='w',
               markerfacecolor=COLORS["scatter_EUR"], markeredgecolor='white', markersize=10,
               label=u'EUR VUS% > Non-EUR (δ > 0)'),
        Line2D([0], [0], marker='o', linestyle='', color='w',
               markerfacecolor=COLORS["scatter_Non"], markeredgecolor='white', markersize=10,
               label=u'Non-EUR VUS% > EUR (δ < 0, inversion)'),
        Line2D([0], [0], color=COLORS["neutral"], linewidth=1.5,
               label=f'Regression  r = {r:+.2f}, p = {p:.3f}'),
        Patch(facecolor=COLORS["ci_band"], alpha=0.4, label='95% CI (bootstrap, n=5,000)'),
    ]
    ax.legend(handles=legend_elements, frameon=True, loc='upper right',
              framealpha=0.95, edgecolor='#ccc')

    ax.set_xlim(x.min() - 4, x.max() + 4)
    y_pad = (y.max() - y.min()) * 0.12
    ax.set_ylim(y.min() - y_pad, y.max() + y_pad)
    plt.tight_layout()

    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"Figure_1.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("[INFO] Figure_1 saved (scatter).")


# ---------------------------------------------------------------------------
# Figure 2 - horizontal bar chart, VUS% by population
# ---------------------------------------------------------------------------

def figure2(pop_df: pd.DataFrame, out_dir: Path):
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

    for bar, val in zip(bars, pop_df["VUS_%"]):
        ax.text(
            bar.get_width() + 0.3, bar.get_y() + bar.get_height() / 2,
            f"{val:.1f}%", va="center", ha="left", fontsize=9,
        )

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=COLORS["EUR"], label="European ancestry"),
        Patch(facecolor=COLORS["Non-EUR"], label="Non-European ancestry"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", frameon=False)

    ax.set_xlabel("Variants of Uncertain Significance (%)")
    ax.set_title(
        "VUS rate by genetic ancestry group\n"
        "17 monogenic diabetes genes, ClinVar x gnomAD v4.0",
        pad=10,
    )
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
    ax.set_xlim(0, max(pop_df["VUS_%"]) * 1.12)

    plt.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"Figure_2.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("[INFO] Figure_2 saved (bar chart).")


# ---------------------------------------------------------------------------
# Figure 3 - heatmap, gene x ancestry VUS%
# ---------------------------------------------------------------------------

def figure3(supp_df: pd.DataFrame, out_dir: Path):
    """Heatmap of VUS rates: 17 genes x 8 gnomAD ancestry groups."""
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

    ax.set_xticks(range(len(POP_COLS)))
    ax.set_xticklabels(list(POP_COLS.keys()), rotation=45, ha="right")
    ax.set_yticks(range(len(MODY_GENES)))
    ax.set_yticklabels(MODY_GENES)

    for i in range(len(MODY_GENES)):
        for j in range(len(POP_COLS)):
            val = data[i, j]
            if not np.isnan(val):
                color = "white" if val > 65 else "black"
                ax.text(j, i, f"{val:.0f}", ha="center", va="center",
                        fontsize=7.5, color=color)
            else:
                ax.text(j, i, "-", ha="center", va="center", fontsize=7, color="#aaa")

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("VUS (%)")

    ax.set_title(
        "VUS rate by gene and genetic ancestry group\n"
        "ClinVar x gnomAD v4.0",
        pad=10,
    )

    plt.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"Figure_3.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("[INFO] Figure_3 saved (heatmap).")


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

    print("[INFO] Loading data ...")
    pop_df = pd.read_csv(pop_csv)
    supp_df = pd.read_csv(supp_csv, low_memory=False)
    g2g_df = pd.read_csv(g2g_csv)

    figure1(g2g_df, fig_dir)   # scatter
    figure2(pop_df, fig_dir)   # bar chart
    figure3(supp_df, fig_dir)  # heatmap

    print("[INFO] All figures saved to figures/ (PNG 300 DPI + SVG).")


if __name__ == "__main__":
    main()
