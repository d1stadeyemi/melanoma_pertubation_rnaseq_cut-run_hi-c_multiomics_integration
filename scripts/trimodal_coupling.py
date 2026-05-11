#!/usr/bin/env python3
"""
Q2: Trimodal coupling — H3K27me3, RNA, and Hi-C compartment changes

Usage:
    python scripts/Q2_trimodal_coupling.py \
        outputs/tables/merged_hic_rnaseq_cutrun_parent_bottom10.txt \
        --outdir outputs/figures/Q2_coupling

Outputs:
    Fig2a_upset_cooccurrence.pdf/png
    Fig2b_spearman_correlations.pdf/png
    Fig2c_concordance_bar.pdf/png
    Q2_coupling_classified.txt
"""

import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import spearmanr
from upsetplot import UpSet, from_memberships

plt.rcParams.update({
    "font.family":       "Liberation Sans",
    "font.size":         11,
    "axes.labelsize":    12,
    "axes.titlesize":    13,
    "xtick.labelsize":   10,
    "ytick.labelsize":   10,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.linewidth":    1.0,
    "xtick.direction":   "out",
    "ytick.direction":   "out",
    "pdf.fonttype":      42,
    "ps.fonttype":       42,
    "savefig.dpi":       600,
    "figure.facecolor":  "white",
    "axes.facecolor":    "white",
})


def parse_args():
    p = argparse.ArgumentParser(
        description="Q2: Trimodal coupling analysis"
    )
    p.add_argument("merged_file", help="merged multi-omics table")
    p.add_argument("--outdir",       default="outputs/figures/Q2_coupling")
    p.add_argument("--lfc_thresh",   type=float, default=1.0,
                   help="log2FC threshold for RNA change (default 1.0)")
    p.add_argument("--me3_thresh",   type=float, default=0.05,
                   help="|H3K27me3 delta| threshold (default 0.05)")
    p.add_argument("--pc1_thresh",   type=float, default=0.03,
                   help="|PC1 delta| threshold (default 0.03)")
    p.add_argument("--padj_thresh",  type=float, default=0.05,
                   help="padj threshold for RNA (default 0.05)")
    return p.parse_args()


def save_fig(fig, outdir, stem):
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        path = os.path.join(outdir, f"{stem}.{ext}")
        fig.savefig(path, bbox_inches="tight", pad_inches=0.05, dpi=600)
        print(f"  Saved -> {path}")
    plt.close(fig)


def classify_coupling(df, args):
    """
    Classify each gene by which modalities show a meaningful change.
    RNA changed:        |log2FC| > lfc_thresh AND padj < padj_thresh
    H3K27me3 changed:  |b10_minus_parental_h3k27me3| > me3_thresh
    Hi-C changed:      |B10-P_Orig| > pc1_thresh

    Also record the direction of each change for concordance analysis.
    """
    # Force numeric
    for col in ["log2FC", "padj", "b10_minus_parental_h3k27me3"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    pc1_col = None
    for c in ["B10-P_Orig", "B10_P_Orig"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
            pc1_col = c
            break

    # Binary change flags
    rna_ok = "log2FC" in df.columns and "padj" in df.columns
    me3_ok = "b10_minus_parental_h3k27me3" in df.columns

    df["rna_changed"] = False
    df["me3_changed"] = False
    df["pc1_changed"] = False

    if rna_ok:
        df["rna_changed"] = (
            (np.abs(df["log2FC"]) > args.lfc_thresh) &
            (df["padj"] < args.padj_thresh)
        ).fillna(False)

    if me3_ok:
        df["me3_changed"] = (
            df["b10_minus_parental_h3k27me3"].abs() > args.me3_thresh
        ).fillna(False)

    if pc1_col:
        df["pc1_changed"] = (
            df[pc1_col].abs() > args.pc1_thresh
        ).fillna(False)

    # Direction of each change
    if rna_ok:
        df["rna_dir"] = np.where(
            df["log2FC"] > 0, "Up",
            np.where(df["log2FC"] < 0, "Down", "None")
        )
    if me3_ok:
        df["me3_dir"] = np.where(
            df["b10_minus_parental_h3k27me3"] > 0, "Gain",
            np.where(df["b10_minus_parental_h3k27me3"] < 0, "Loss", "None")
        )
    if pc1_col:
        df["pc1_dir"] = np.where(
            df[pc1_col] > 0, "To_A",
            np.where(df[pc1_col] < 0, "To_B", "None")
        )

    # Coupling label
    def make_label(row):
        parts = []
        if row["rna_changed"]:  parts.append("RNA")
        if row["me3_changed"]:  parts.append("H3K27me3")
        if row["pc1_changed"]:  parts.append("Hi-C")
        return "+".join(parts) if parts else "None"

    df["coupling_state"] = df.apply(make_label, axis=1)

    print("[INFO] Coupling state distribution:")
    print(df["coupling_state"].value_counts().to_string())
    return df, pc1_col


def plot_upset(df, outdir):
    """
    Fig 2a: UpSet plot — standard for multi-set intersection visualization.
    Shows what fraction of changes are trimodal vs modality-specific.
    """
    memberships = []
    for _, row in df.iterrows():
        m = []
        if row.get("rna_changed", False): m.append("Transcription")
        if row.get("me3_changed", False): m.append("H3K27me3")
        if row.get("pc1_changed", False): m.append("Hi-C Compartment")
        memberships.append(m)

    upset_data = from_memberships(memberships)

    fig = plt.figure(figsize=(9, 5))
    upset = UpSet(
        upset_data,
        subset_size="count",
        show_counts=True,
        sort_by="cardinality",
        element_size=42
    )
    upset.style_subsets(
        present=["Transcription", "H3K27me3", "Hi-C Compartment"],
        facecolor="#C0392B",
        label="All 3 co-changed"
    )
    upset.plot(fig)
    plt.suptitle(
        "Co-occurrence of Multi-omics Changes\nParental \u2192 Bottom10 (A375 Melanoma)",
        fontweight="bold", fontsize=12, y=1.02
    )
    save_fig(fig, outdir, "Fig2a_upset_cooccurrence")


def plot_spearman(df, pc1_col, outdir):
    """
    Fig 2b: Pairwise Spearman correlations between the three modality
    delta values. Each point is a gene. Rho and p annotated.
    This tests whether the three layers are mechanistically linked
    or can dissociate.
    """
    pairs = []
    if pc1_col:
        pairs = [
            ("b10_minus_parental_h3k27me3", "log2FC",
             "\u0394H3K27me3 (B10 \u2212 Parental)",
             r"$\log_2 FC$ (RNA-seq)", "#1A5276"),
            ("b10_minus_parental_h3k27me3", pc1_col,
             "\u0394H3K27me3 (B10 \u2212 Parental)",
             "\u0394PC1 (Hi-C)", "#C0392B"),
            (pc1_col, "log2FC",
             "\u0394PC1 (Hi-C)",
             r"$\log_2 FC$ (RNA-seq)", "#27AE60"),
        ]
    elif "log2FC" in df.columns:
        pairs = [
            ("b10_minus_parental_h3k27me3", "log2FC",
             "\u0394H3K27me3 (B10 \u2212 Parental)",
             r"$\log_2 FC$ (RNA-seq)", "#1A5276"),
        ]

    if not pairs:
        print("[WARN] No valid column pairs for Spearman. Skipping.")
        return

    fig, axes = plt.subplots(1, len(pairs),
                              figsize=(4.5 * len(pairs), 4.5))
    if len(pairs) == 1:
        axes = [axes]

    for ax, (xcol, ycol, xlabel, ylabel, color) in zip(axes, pairs):
        if xcol not in df.columns or ycol not in df.columns:
            ax.set_visible(False)
            continue
        sub = df[[xcol, ycol]].dropna()
        if len(sub) < 10:
            ax.set_visible(False)
            continue

        rho, p = spearmanr(sub[xcol], sub[ycol])

        # Scatter with rasterized points (thousands of genes)
        ax.scatter(sub[xcol], sub[ycol],
                   c=color, s=5, alpha=0.25,
                   linewidths=0, rasterized=True)

        # Trend line
        z = np.polyfit(sub[xcol].values, sub[ycol].values, 1)
        xl = np.linspace(sub[xcol].min(), sub[xcol].max(), 100)
        ax.plot(xl, np.polyval(z, xl),
                color="black", lw=1.2, ls="--", alpha=0.7)

        p_str = f"{p:.2e}" if p >= 1e-4 else r"$< 10^{-4}$"
        ax.text(0.05, 0.95,
                f"\u03c1 = {rho:.3f}\np = {p_str}",
                transform=ax.transAxes, fontsize=9,
                va="top",
                bbox=dict(boxstyle="round,pad=0.25",
                          fc="white", ec="none", alpha=0.85))

        ax.axhline(0, color="#888888", lw=0.6, ls=":")
        ax.axvline(0, color="#888888", lw=0.6, ls=":")
        ax.set_xlabel(xlabel, labelpad=5)
        ax.set_ylabel(ylabel, labelpad=5)

    fig.suptitle(
        "Pairwise Correlations: H3K27me3 \u2194 RNA \u2194 Hi-C",
        fontweight="bold", fontsize=12, y=1.02
    )
    plt.tight_layout()
    save_fig(fig, outdir, "Fig2b_spearman_correlations")


def plot_concordance_bar(df, outdir):
    """
    Fig 2c: For genes with a compartment change call, what fraction
    also show concordant RNA and H3K27me3 changes?

    B->A switches should correlate with:
        RNA upregulation (gene activation)
        H3K27me3 loss (derepression)
    A->B switches should correlate with:
        RNA downregulation
        H3K27me3 gain

    This directly tests the mechanistic model proposed in the report.
    """
    def get_comp_dir(row):
        for col in ["P->B10_Both_original",
                    "P->B10_Switch_Orig",
                    "P->B10_Shift_Orig"]:
            v = str(row.get(col, "")).strip()
            if v not in ["-", "", "nan", "None"]:
                if "toA" in v or v == "BtoA" or v == "TowardA":
                    return "B\u2192A"
                if "toB" in v or v == "AtoB" or v == "TowardB":
                    return "A\u2192B"
        return "Stable"

    df["comp_dir"] = df.apply(get_comp_dir, axis=1)

    # For each compartment direction, compute fraction with concordant changes
    comp_dirs = ["B\u2192A", "A\u2192B", "Stable"]
    scenarios = [
        ("Concordant RNA", lambda d, g:
            ((g == "B\u2192A") & (d["rna_dir"] == "Up")) |
            ((g == "A\u2192B") & (d["rna_dir"] == "Down"))),
        ("Concordant H3K27me3", lambda d, g:
            ((g == "B\u2192A") & (d["me3_dir"] == "Loss")) |
            ((g == "A\u2192B") & (d["me3_dir"] == "Gain"))),
        ("Both concordant", lambda d, g:
            (((g == "B\u2192A") & (d["rna_dir"] == "Up") & (d["me3_dir"] == "Loss")) |
             ((g == "A\u2192B") & (d["rna_dir"] == "Down") & (d["me3_dir"] == "Gain")))),
    ]
    colors_bar = ["#C0392B", "#1A5276", "#8E44AD"]

    x     = np.arange(len(comp_dirs))
    width = 0.25
    offsets = [-width, 0, width]

    fig, ax = plt.subplots(figsize=(7, 5))

    for (label, fn), color, offset in zip(scenarios, colors_bar, offsets):
        fracs = []
        for gd in comp_dirs:
            sub = df[df["comp_dir"] == gd]
            n   = len(sub)
            if n == 0:
                fracs.append(0)
                continue
            try:
                n_yes = fn(sub, gd).sum()
            except Exception:
                n_yes = 0
            fracs.append(n_yes / n * 100)
        ax.bar(x + offset, fracs, width=width * 0.9,
               color=color, label=label,
               edgecolor="white", linewidth=0.4)

    # Sample size annotations
    for xi, gd in zip(x, comp_dirs):
        n = (df["comp_dir"] == gd).sum()
        ax.text(xi, -3, f"n={n:,}", ha="center", va="top",
                fontsize=8, color="#555555")

    ax.set_xticks(x)
    ax.set_xticklabels(comp_dirs, fontsize=11)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylim(bottom=-6)
    ax.set_ylabel("% of genes with concordant change", labelpad=6)
    ax.set_xlabel("Compartment Change Direction", labelpad=6)
    ax.set_title(
        "Directional Concordance Across Modalities\nParental \u2192 Bottom10",
        fontweight="bold", pad=10
    )
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    plt.tight_layout()
    save_fig(fig, outdir, "Fig2c_concordance_bar")


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print("[INFO] Loading merged table...")
    df = pd.read_csv(args.merged_file, sep="\t", low_memory=False)
    print(f"  {len(df):,} rows x {len(df.columns)} columns")

    df, pc1_col = classify_coupling(df, args)

    # Save coupling classification
    out_path = os.path.join(args.outdir, "Q2_coupling_classified.txt")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"[INFO] Saved -> {out_path}")

    print("\n[INFO] Generating Fig 2a: UpSet plot...")
    plot_upset(df, args.outdir)

    print("[INFO] Generating Fig 2b: Spearman correlations...")
    plot_spearman(df, pc1_col, args.outdir)

    print("[INFO] Generating Fig 2c: Concordance bar chart...")
    plot_concordance_bar(df, args.outdir)

    print(f"\n[INFO] Q2 complete -> {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()