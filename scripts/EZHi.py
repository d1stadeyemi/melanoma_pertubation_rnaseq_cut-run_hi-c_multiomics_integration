#!/usr/bin/env python3
"""
Q4: EZH1/2 inhibition reversibility — which loci are pharmacologically reversible?

Usage:
    python scripts/Q4_EZHi.py \
        outputs/tables/cut_run_signal_table.txt \
        outputs/tables/merged_hic_rnaseq_cutrun_parent_bottom10.txt \
        --outdir outputs/figures/Q4_EZHi

Outputs:
    Fig4a_EZHi_classification.pdf/png
    Fig4b_H3K27me3_trajectory.pdf/png
    Fig4c_RNA_by_class.pdf/png
    Fig4d_compartment_by_class.pdf/png
    Q4_EZHi_classes.txt              <- input for R GO enrichment
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
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests

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

CLASS_COLORS = {
    "Reversible":   "#27AE60",
    "Persistent":   "#C0392B",
    "Loss_B10":     "#1A5276",
    "Unresponsive": "#95A5A6",
}


def parse_args():
    p = argparse.ArgumentParser(description="Q4: EZHi reversibility")
    p.add_argument("cutrun_file", help="cut_run_signal_table.txt")
    p.add_argument("merged_file", help="merged multi-omics table")
    p.add_argument("--outdir",      default="outputs/figures/Q4_EZHi")
    p.add_argument("--me3_thresh",  type=float, default=0.03,
                   help="H3K27me3 delta threshold (default 0.03)")
    p.add_argument("--lfc_thresh",  type=float, default=1.0)
    p.add_argument("--padj_thresh", type=float, default=0.05)
    return p.parse_args()


def save_fig(fig, outdir, stem):
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        path = os.path.join(outdir, f"{stem}.{ext}")
        fig.savefig(path, bbox_inches="tight", pad_inches=0.05, dpi=600)
        print(f"  Saved -> {path}")
    plt.close(fig)


def find_col(df, must, must_not=None):
    must_not = must_not or []
    cands = []
    for col in df.columns:
        c = col.lower()
        if all(k in c for k in must) and not any(k in c for k in must_not):
            cands.append(col)
    pref = [c for c in cands
            if not any(f"_r{i}" in c.lower() for i in range(1, 10))]
    return pref[0] if pref else (cands[0] if cands else None)


def load_data(cutrun_file, merged_file):
    """Load both files and detect the four H3K27me3 condition columns."""
    cutrun = pd.read_csv(cutrun_file, sep="\t", low_memory=False)
    merged = pd.read_csv(merged_file, sep="\t", low_memory=False)

    cols = {
        "parental":     find_col(cutrun, ["parental", "h3k27me3", "mean"],
                                 ["minus", "unc", "ezh", "b10", "bottom"]),
        "bottom10":     find_col(cutrun, ["b10", "h3k27me3", "mean"],
                                 ["minus", "unc", "ezh", "parental", "control"]),
        "parental_ezhi":find_col(cutrun, ["parental", "unc", "h3k27me3", "mean"],
                                 ["minus", "b10", "bottom"]),
        "bottom10_ezhi":find_col(cutrun, ["b10", "unc", "h3k27me3", "mean"],
                                 ["minus", "parental", "control"]),
    }

    # Fix: bottom10 should not match parental
    if cols["bottom10"] and "parental" in cols["bottom10"].lower():
        for col in cutrun.columns:
            c = col.lower()
            if ("b10" in c or "bottom" in c) and "h3k27me3" in c and \
                    "mean" in c and "parental" not in c and "control" not in c and \
                    "unc" not in c and "minus" not in c:
                if not any(f"_r{i}" in c for i in range(1, 10)):
                    cols["bottom10"] = col
                    break

    print("[INFO] H3K27me3 columns:")
    for k, v in cols.items():
        print(f"  {k}: {v}")

    # Convert to numeric
    for col in cols.values():
        if col and col in cutrun.columns:
            cutrun[col] = pd.to_numeric(cutrun[col], errors="coerce")

    # CUT&RUN subset
    keep = ["gene_name"] + [v for v in cols.values() if v and v in cutrun.columns]
    cr_sub = cutrun[keep].drop_duplicates("gene_name").copy()

    # Merged subset
    out_cols = ["gene_name", "log2FC", "padj",
                "b10_minus_parental_h3k27me3"]
    for c in ["P->B10_Switch_Orig", "P->B10_Shift_Orig",
              "P->B10_Both_original"]:
        if c in merged.columns:
            out_cols.append(c)
    m_sub = merged[[c for c in out_cols if c in merged.columns]
                   ].drop_duplicates("gene_name").copy()

    df = cr_sub.merge(m_sub, on="gene_name", how="inner")
    for c in ["log2FC", "padj", "b10_minus_parental_h3k27me3"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    print(f"[INFO] Merged: {len(df):,} genes")
    return df, cols


def classify_reversibility(df, cols, me3_thresh):
    """
    Classify H3K27me3 trajectory per gene.

    Gained in Bottom10:
        bottom10 signal > parental signal + threshold
    Reversible:
        gained in Bottom10 AND bottom10_ezhi signal returned
        toward parental (i.e., < bottom10 - threshold)
    Persistent:
        gained in Bottom10 AND NOT reversed by EZHi
    Loss_B10:
        lost in Bottom10 (decreased vs parental)
    Unresponsive:
        no meaningful change
    """
    par  = cols.get("parental")
    bot  = cols.get("bottom10")
    pezi = cols.get("parental_ezhi")
    bezi = cols.get("bottom10_ezhi")

    has_all = all(c and c in df.columns for c in [par, bot, pezi, bezi])

    if has_all:
        df["delta_base"] = df[bot]  - df[par]   # Bottom10 - Parental
        df["delta_ezhi"] = df[bezi] - df[bot]   # EZHi Bottom10 - Bottom10

        def classify(row):
            db = row["delta_base"]
            de = row["delta_ezhi"]
            if pd.isna(db) or pd.isna(de):
                return "Unknown"
            if db > me3_thresh and de < -me3_thresh:
                return "Reversible"
            if db > me3_thresh and de >= -me3_thresh:
                return "Persistent"
            if db < -me3_thresh:
                return "Loss_B10"
            return "Unresponsive"

        df["EZHi_class"] = df.apply(classify, axis=1)
    else:
        # Fallback: use pre-computed delta
        print("[WARN] Not all four EZHi columns found. Using pre-computed delta.")
        if "b10_minus_parental_h3k27me3" in df.columns:
            df["delta_base"] = df["b10_minus_parental_h3k27me3"]
            df["EZHi_class"] = np.where(
                df["delta_base"] > me3_thresh,  "Gained_B10",
                np.where(df["delta_base"] < -me3_thresh,
                         "Loss_B10", "Unresponsive")
            )
        else:
            df["EZHi_class"] = "Unknown"

    print("[INFO] EZHi class distribution:")
    print(df["EZHi_class"].value_counts().to_string())
    return df, has_all


def plot_classification(df, outdir):
    """Fig 4a: Bar chart of EZHi response classes."""
    cls_order = ["Reversible", "Persistent", "Loss_B10",
                 "Unresponsive", "Gained_B10", "Unknown"]
    cls_order = [c for c in cls_order if c in df["EZHi_class"].unique()]
    counts    = [df["EZHi_class"].value_counts().get(c, 0) for c in cls_order]
    colors    = [CLASS_COLORS.get(c, "#95A5A6") for c in cls_order]

    labels = {
        "Reversible":   "Reversible\n(EZHi restores)",
        "Persistent":   "Persistent\n(EZHi-resistant)",
        "Loss_B10":     "Lost in\nBottom10",
        "Unresponsive": "Unresponsive",
        "Gained_B10":   "Gained in\nBottom10",
        "Unknown":      "Unknown",
    }

    fig, ax = plt.subplots(figsize=(max(5, len(cls_order) * 1.5), 4.5))
    bars = ax.bar(
        [labels.get(c, c) for c in cls_order], counts,
        color=colors, edgecolor="white", linewidth=0.5, width=0.6
    )
    for bar, val in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(counts) * 0.012,
                f"{val:,}", ha="center", va="bottom", fontsize=9)

    ax.set_ylabel("Number of genes", labelpad=6)
    ax.set_title(
        "H3K27me3 Response to EZH1/2 Inhibition (UNC1999)\nGene Classification",
        fontweight="bold", pad=10
    )
    plt.tight_layout()
    save_fig(fig, outdir, "Fig4a_EZHi_classification")


def plot_trajectory(df, cols, has_all, outdir):
    """
    Fig 4b: H3K27me3 signal trajectory.
    Thick line = class mean ± SEM. Individual gene lines faint.
    Shows the reversal pattern — the defining result of Q4.
    """
    if not has_all:
        print("[INFO] EZHi columns incomplete. Generating simplified trajectory.")

    par  = cols.get("parental")
    bot  = cols.get("bottom10")
    pezi = cols.get("parental_ezhi")
    bezi = cols.get("bottom10_ezhi")

    if has_all:
        cond_cols  = [par, bot, pezi, bezi]
        cond_labs  = ["Parental", "Bottom10",
                      "Parental\n+EZHi", "Bottom10\n+EZHi"]
        positions  = [0, 1, 2, 3]
    else:
        cond_cols  = [c for c in [par, bot] if c and c in df.columns]
        cond_labs  = ["Parental", "Bottom10"][:len(cond_cols)]
        positions  = list(range(len(cond_cols)))

    classes = [c for c in ["Reversible", "Persistent", "Unresponsive"]
               if c in df["EZHi_class"].unique()]
    if not classes:
        classes = df["EZHi_class"].unique().tolist()[:3]

    n_cls = len(classes)
    if n_cls == 0:
        return

    fig, axes = plt.subplots(1, n_cls, figsize=(4.0 * n_cls, 4.5),
                              sharey=True)
    if n_cls == 1:
        axes = [axes]

    for ax, cls in zip(axes, classes):
        sub   = df[df["EZHi_class"] == cls]
        color = CLASS_COLORS.get(cls, "#95A5A6")

        # Sample gene lines
        n_samp = min(200, len(sub))
        sub_s  = sub.sample(n_samp, random_state=42) if len(sub) > 200 else sub

        for _, row in sub_s.iterrows():
            vals = [pd.to_numeric(row.get(c, np.nan), errors="coerce")
                    for c in cond_cols if c in df.columns]
            if len(vals) != len(positions) or any(pd.isna(v) for v in vals):
                continue
            ax.plot(positions, vals, color=color, alpha=0.06, lw=0.8)

        # Class mean ± SEM
        means, sems = [], []
        for c in [c for c in cond_cols if c and c in df.columns]:
            vals = pd.to_numeric(sub[c], errors="coerce").dropna()
            means.append(vals.mean())
            sems.append(vals.sem())

        ax.errorbar(
            positions[:len(means)], means, yerr=sems,
            color=color, lw=2.5, marker="o",
            markersize=7, capsize=4, capthick=1.5, zorder=5
        )
        ax.set_xticks(positions[:len(cond_labs)])
        ax.set_xticklabels(cond_labs[:len(means)], fontsize=9)
        ax.set_title(f"{cls}\n(n={len(sub):,})",
                     fontweight="bold", fontsize=10, color=color)
        ax.set_xlabel("Condition", labelpad=5)
        ax.axhline(0, color="#888888", lw=0.6, ls=":")

    axes[0].set_ylabel("H3K27me3 Signal (CPM normalized)", labelpad=6)
    fig.suptitle(
        "H3K27me3 Signal Trajectory\nby EZH1/2 Inhibition Response Class",
        fontweight="bold", fontsize=12, y=1.02
    )
    plt.tight_layout()
    save_fig(fig, outdir, "Fig4b_H3K27me3_trajectory")


def plot_rna_by_class(df, outdir):
    """Fig 4c: log2FC by EZHi class. Tests transcriptional coupling."""
    if "log2FC" not in df.columns:
        return

    df["log2FC"] = pd.to_numeric(df["log2FC"], errors="coerce")
    classes = [c for c in ["Reversible", "Persistent",
                            "Loss_B10", "Unresponsive"]
               if c in df["EZHi_class"].unique()]
    data = [df[df["EZHi_class"] == c]["log2FC"].dropna().values
            for c in classes]
    valid = [(c, d) for c, d in zip(classes, data) if len(d) > 1]
    if not valid:
        return
    cls_v, data_v = zip(*valid)

    fig, ax = plt.subplots(figsize=(max(4.5, len(cls_v) * 1.5), 4.5))
    parts = ax.violinplot(data_v, showmedians=True,
                           showextrema=False, widths=0.7)
    for pc, cls in zip(parts["bodies"], cls_v):
        pc.set_facecolor(CLASS_COLORS.get(cls, "#95A5A6"))
        pc.set_alpha(0.75)
        pc.set_edgecolor("none")
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(1.5)

    if len(data_v) >= 2:
        _, p = kruskal(*data_v)
        ax.text(0.97, 0.97, f"Kruskal-Wallis p={p:.2e}",
                transform=ax.transAxes, ha="right", va="top", fontsize=8)

    ax.set_xticks(range(1, len(cls_v) + 1))
    ax.set_xticklabels(cls_v, fontsize=9, rotation=15, ha="right")
    ax.axhline(0, color="#888888", lw=0.8, ls="--")
    ax.set_ylabel("log\u2082FC (RNA-seq, Bottom10 / Parental)", labelpad=6)
    ax.set_title(
        "Transcriptional Change by\nH3K27me3 EZHi Response Class",
        fontweight="bold", pad=10
    )
    plt.tight_layout()
    save_fig(fig, outdir, "Fig4c_RNA_by_class")


def plot_compartment_by_class(df, outdir):
    """Fig 4d: Compartment change co-occurrence by EZHi class."""
    def any_switch(row):
        for col in ["P->B10_Both_original",
                    "P->B10_Switch_Orig",
                    "P->B10_Shift_Orig"]:
            v = str(row.get(col, "")).strip()
            if v not in ["-", "", "nan", "None"]:
                return True
        return False

    df["has_comp"] = df.apply(any_switch, axis=1)
    classes = [c for c in ["Reversible", "Persistent",
                            "Loss_B10", "Unresponsive"]
               if c in df["EZHi_class"].unique()]

    fracs  = [df[df["EZHi_class"] == c]["has_comp"].mean() * 100
              for c in classes]
    ns     = [len(df[df["EZHi_class"] == c]) for c in classes]
    colors = [CLASS_COLORS.get(c, "#95A5A6") for c in classes]

    fig, ax = plt.subplots(figsize=(max(4.5, len(classes) * 1.5), 4.5))
    bars = ax.bar(classes, fracs, color=colors,
                  edgecolor="white", linewidth=0.5, width=0.6)
    for bar, frac, n in zip(bars, fracs, ns):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.5,
                f"{frac:.1f}%\n(n={n:,})",
                ha="center", va="bottom", fontsize=8)

    ax.set_ylabel("% genes with compartment change", labelpad=6)
    ax.set_xlabel("EZHi Response Class", labelpad=6)
    ax.set_title(
        "Compartment Change Co-occurrence\nby H3K27me3 EZHi Response Class",
        fontweight="bold", pad=10
    )
    plt.tight_layout()
    save_fig(fig, outdir, "Fig4d_compartment_by_class")


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df, cols     = load_data(args.cutrun_file, args.merged_file)
    df, has_all  = classify_reversibility(df, cols, args.me3_thresh)

    # Save for GO enrichment
    out_cols = (["gene_name", "EZHi_class"] +
                [c for c in ["delta_base", "delta_ezhi", "log2FC", "padj",
                              "b10_minus_parental_h3k27me3",
                              "P->B10_Switch_Orig",
                              "P->B10_Both_original"]
                 if c in df.columns])
    out_path = os.path.join(args.outdir, "Q4_EZHi_classes.txt")
    df[[c for c in out_cols if c in df.columns]].to_csv(
        out_path, sep="\t", index=False
    )
    print(f"[INFO] Saved -> {out_path}")

    print("\n[INFO] Fig 4a: Classification bar...")
    plot_classification(df, args.outdir)

    print("[INFO] Fig 4b: Trajectory...")
    plot_trajectory(df, cols, has_all, args.outdir)

    print("[INFO] Fig 4c: RNA by class...")
    plot_rna_by_class(df, args.outdir)

    print("[INFO] Fig 4d: Compartment by class...")
    plot_compartment_by_class(df, args.outdir)

    print(f"\n[INFO] Q4 complete -> {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()