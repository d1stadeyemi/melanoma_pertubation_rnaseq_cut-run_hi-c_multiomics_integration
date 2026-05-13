#!/usr/bin/env python3
"""
Q3: Bivalent chromatin state primes migration genes for activation

Usage:
    python scripts/Q3_bivalency.py \
        outputs/tables/cut_run_signal_table.txt \
        outputs/tables/merged_hic_rnaseq_cutrun_parent_bottom10.txt \
        --outdir outputs/figures/Q3_bivalency

Outputs:
    Fig3a_chromatin_state_pie.pdf/png
    Fig3b_fate_by_state.pdf/png
    Fig3c_bivalency_vs_PC1.pdf/png
    Fig3d_bivalency_vs_RNA.pdf/png
    Q3_chromatin_states.txt           <- input for R GO enrichment
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
from scipy.stats import spearmanr, fisher_exact
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

STATE_COLORS = {
    "Active":   "#F4C542",
    "Bivalent": "#8E44AD",
    "Silenced": "#1A5276",
    "Inactive": "#D5D8DC",
}


def parse_args():
    p = argparse.ArgumentParser(description="Q3: Bivalency analysis")
    p.add_argument("cutrun_file",  help="cut_run_signal_table.txt")
    p.add_argument("merged_file",  help="merged multi-omics table")
    p.add_argument("--outdir",       default="outputs/figures/Q3_bivalency")
    p.add_argument("--me3_pctile",   type=float, default=60,
                   help="H3K27me3 high percentile threshold (default 60)")
    p.add_argument("--ac_pctile",    type=float, default=60,
                   help="H3K27Ac high percentile threshold (default 60)")
    p.add_argument("--me1_pctile",   type=float, default=40,
                   help="H3K4me1 high percentile threshold (default 40)")
    p.add_argument("--lfc_thresh",   type=float, default=1.0)
    p.add_argument("--padj_thresh",  type=float, default=0.05)
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


def load_and_merge(cutrun_file, merged_file):
    """Load and merge CUT&RUN features with compartment/RNA outcomes."""
    cutrun = pd.read_csv(cutrun_file, sep="\t", low_memory=False)
    merged = pd.read_csv(merged_file, sep="\t", low_memory=False)

    # Detect Parental mark columns
    cols = {
        "me3": find_col(cutrun, ["parental", "h3k27me3", "mean"],
                        ["minus", "unc", "ezh", "b10", "bottom"]),
        "ac":  find_col(cutrun, ["parental", "h3k27ac",  "mean"],
                        ["minus", "unc", "ezh", "b10", "bottom"]),
        "me1": find_col(cutrun, ["parental", "h3k4me1",  "mean"],
                        ["minus", "unc", "ezh", "b10", "bottom"]),
        "me3h":find_col(cutrun, ["parental", "h3k9me3",  "mean"],
                        ["minus", "unc", "ezh", "b10", "bottom"]),
    }
    print("[INFO] Parental mark columns:")
    for k, v in cols.items():
        print(f"  {k}: {v}")

    # Check minimum
    if cols["me3"] is None or cols["me1"] is None:
        sys.exit("[ERROR] Need at least H3K27me3 and H3K4me1 columns.")

    # Build CUT&RUN features
    keep = ["gene_name"] + [v for v in cols.values() if v]
    cr_sub = cutrun[keep].drop_duplicates("gene_name").copy()
    for v in cols.values():
        if v:
            cr_sub[v] = pd.to_numeric(cr_sub[v], errors="coerce")

    # Build merged outcomes
    out_cols = ["gene_name", "log2FC", "padj",
                "b10_minus_parental_h3k27me3"]
    for c in ["P->B10_Switch_Orig", "P->B10_Shift_Orig",
              "P->B10_Both_original", "B10-P_Orig",
              "Parental_orig", "Bottom10_orig"]:
        if c in merged.columns:
            out_cols.append(c)

    merged_sub = merged[[c for c in out_cols if c in merged.columns]
                        ].drop_duplicates("gene_name").copy()
    for c in ["log2FC", "padj", "b10_minus_parental_h3k27me3"]:
        if c in merged_sub.columns:
            merged_sub[c] = pd.to_numeric(merged_sub[c], errors="coerce")

    df = cr_sub.merge(merged_sub, on="gene_name", how="inner")
    print(f"[INFO] Merged: {len(df):,} genes")
    return df, cols


def classify_chromatin_state(df, cols, args):
    """
    Classify each gene into chromatin state using percentile thresholds.
    Data-driven thresholds are preferred over absolute cutoffs —
    this avoids bias from CPM scaling differences between marks.

    States:
        Active:   high H3K27Ac, low H3K27me3 (actively transcribed)
        Bivalent: high H3K27me3 AND high H3K4me1 (poised/bivalent)
        Silenced: high H3K27me3, low H3K4me1 (stably silenced)
        Inactive: low all marks (unmarked heterochromatin)
    """
    me3_col = cols.get("me3")
    me1_col = cols.get("me1")
    ac_col  = cols.get("ac")

    me3_vals = df[me3_col].dropna().values if me3_col else np.array([])
    me1_vals = df[me1_col].dropna().values if me1_col else np.array([])
    ac_vals  = df[ac_col].dropna().values  if ac_col  else np.array([])

    me3_thr = np.percentile(me3_vals, args.me3_pctile) if len(me3_vals) else None
    me1_thr = np.percentile(me1_vals, args.me1_pctile) if len(me1_vals) else None
    ac_thr  = np.percentile(ac_vals,  args.ac_pctile)  if len(ac_vals)  else None

    print(f"[INFO] State thresholds:")
    print(f"  H3K27me3 >= {me3_thr:.4f} (p{args.me3_pctile})")
    print(f"  H3K4me1  >= {me1_thr:.4f} (p{args.me1_pctile})")
    if ac_thr:
        print(f"  H3K27Ac  >= {ac_thr:.4f}  (p{args.ac_pctile})")

    def assign(row):
        m3 = row[me3_col] if me3_col else np.nan
        m1 = row[me1_col] if me1_col else np.nan
        ac = row[ac_col]  if ac_col  else np.nan

        hi_me3 = pd.notna(m3) and me3_thr and m3 >= me3_thr
        hi_me1 = pd.notna(m1) and me1_thr and m1 >= me1_thr
        hi_ac  = pd.notna(ac) and ac_thr  and ac >= ac_thr

        if hi_ac and not hi_me3:
            return "Active"
        if hi_me3 and hi_me1:
            return "Bivalent"
        if hi_me3 and not hi_me1:
            return "Silenced"
        return "Inactive"

    df["chromatin_state"] = df.apply(assign, axis=1)

    # Bivalency index
    if me3_col and me1_col:
        df["bivalency_index"] = (
            df[me1_col] / (df[me3_col] + 1e-9)
        )

    print("[INFO] Chromatin state distribution:")
    print(df["chromatin_state"].value_counts().to_string())
    return df


def classify_fate(df, args):
    """Classify fate of each gene in Bottom10."""
    # Compartment change
    def any_switch(row):
        for col in ["P->B10_Both_original",
                    "P->B10_Switch_Orig",
                    "P->B10_Shift_Orig"]:
            v = str(row.get(col, "")).strip()
            if v not in ["-", "", "nan", "None"]:
                return True
        return False

    df["fate_comp"] = df.apply(any_switch, axis=1)

    # B->A specifically
    def is_b_to_a(row):
        for col in ["P->B10_Both_original",
                    "P->B10_Switch_Orig",
                    "P->B10_Shift_Orig"]:
            v = str(row.get(col, "")).strip()
            if "toA" in v or v in ["BtoA", "TowardA"]:
                return True
        return False

    df["fate_B_to_A"]  = df.apply(is_b_to_a, axis=1)

    # RNA activated
    if "log2FC" in df.columns and "padj" in df.columns:
        df["fate_rna_up"] = (
            (df["log2FC"] > args.lfc_thresh) &
            (df["padj"] < args.padj_thresh)
        ).fillna(False)
        df["fate_rna_down"] = (
            (df["log2FC"] < -args.lfc_thresh) &
            (df["padj"] < args.padj_thresh)
        ).fillna(False)
    else:
        df["fate_rna_up"] = df["fate_rna_down"] = False

    df["fate_both"] = df["fate_comp"] & (df["fate_rna_up"] | df["fate_rna_down"])
    return df


def plot_state_pie(df, outdir):
    """Fig 3a: Chromatin state distribution in Parental cells."""
    counts = df["chromatin_state"].value_counts()
    cats   = [c for c in ["Active", "Bivalent", "Silenced", "Inactive"]
              if c in counts.index]
    vals   = [counts[c] for c in cats]
    colors = [STATE_COLORS[c] for c in cats]
    explode = [0.06 if c == "Bivalent" else 0 for c in cats]

    fig, ax = plt.subplots(figsize=(5, 5))
    wedges, texts, autotexts = ax.pie(
        vals, labels=cats, colors=colors,
        autopct="%1.1f%%", explode=explode,
        startangle=90, pctdistance=0.78,
        textprops={"fontsize": 10}
    )
    for at in autotexts:
        at.set_fontsize(9)
        at.set_fontweight("bold")
    ax.set_title(
        "Parental A375 Chromatin State Distribution\n"
        "(H3K27me3 + H3K4me1 + H3K27Ac)",
        fontweight="bold", pad=15
    )
    save_fig(fig, outdir, "Fig3a_chromatin_state_pie")


def plot_fate_by_state(df, outdir):
    """
    Fig 3b: Fate of each chromatin state in Bottom10.
    Fisher's exact test: Bivalent vs Silenced for B->A activation.
    Key test of the bivalency priming hypothesis.
    """
    states = [c for c in ["Active", "Bivalent", "Silenced", "Inactive"]
              if c in df["chromatin_state"].unique()]

    fates  = ["fate_B_to_A", "fate_rna_up", "fate_both"]
    flabs  = ["B\u2192A\nSwitch", "RNA\nActivation", "Both"]
    fcolors = ["#1A5276", "#C0392B", "#8E44AD"]

    x      = np.arange(len(states))
    width  = 0.25
    offs   = [-width, 0, width]

    fig, ax = plt.subplots(figsize=(7, 5))
    for fate_col, flab, fcol, off in zip(fates, flabs, fcolors, offs):
        fracs = []
        for state in states:
            sub = df[df["chromatin_state"] == state]
            n   = len(sub)
            fracs.append(sub[fate_col].sum() / n * 100 if n > 0 else 0)
        ax.bar(x + off, fracs, width=width * 0.9,
               color=fcol, label=flab,
               edgecolor="white", linewidth=0.4)

    # Sample sizes below x-axis
    for xi, state in zip(x, states):
        n = (df["chromatin_state"] == state).sum()
        ax.text(xi, -3.5, f"n={n:,}", ha="center", va="top",
                fontsize=8, color="#444444")

    # Fisher's exact: Bivalent vs Silenced for B->A activation
    if "Bivalent" in states and "Silenced" in states:
        biv = df[df["chromatin_state"] == "Bivalent"]["fate_B_to_A"]
        sil = df[df["chromatin_state"] == "Silenced"]["fate_B_to_A"]
        ct  = [[biv.sum(), len(biv) - biv.sum()],
               [sil.sum(), len(sil) - sil.sum()]]
        try:
            _, p = fisher_exact(ct)
            ax.text(0.98, 0.97,
                    f"Bivalent vs Silenced\n(B\u2192A): p={p:.3e}",
                    transform=ax.transAxes,
                    ha="right", va="top", fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.3",
                              fc="white", ec="#cccccc", alpha=0.9))
        except Exception:
            pass

    ax.set_xticks(x)
    ax.set_xticklabels(states, fontsize=11)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylim(bottom=-6)
    ax.set_ylabel("% of genes (within state)", labelpad=6)
    ax.set_xlabel("Parental Chromatin State", labelpad=6)
    ax.set_title(
        "Fate of Parental Chromatin States in Bottom10\n"
        "Bivalent loci preferentially switch B\u2192A and activate",
        fontweight="bold", pad=10
    )
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    plt.tight_layout()
    save_fig(fig, outdir, "Fig3b_fate_by_state")


def plot_bivalency_scatter(df, x_col, xlabel, stem, outdir):
    """Fig 3c/d: Bivalency index vs PC1 delta or log2FC."""
    if "bivalency_index" not in df.columns or x_col not in df.columns:
        print(f"[WARN] Missing bivalency_index or {x_col}. Skipping.")
        return

    df[x_col] = pd.to_numeric(df[x_col], errors="coerce")
    sub = df[["bivalency_index", x_col, "chromatin_state"]].dropna()
    sub["bivalency_index"] = np.clip(sub["bivalency_index"], 0, 20)

    rho, p = spearmanr(sub["bivalency_index"], sub[x_col])

    fig, ax = plt.subplots(figsize=(5.5, 5))
    draw_order = ["Inactive", "Silenced", "Active", "Bivalent"]
    for state in draw_order:
        mask = sub["chromatin_state"] == state
        ax.scatter(
            sub.loc[mask, "bivalency_index"],
            sub.loc[mask, x_col],
            c=STATE_COLORS.get(state, "#95A5A6"),
            s=8 if state in ["Inactive"] else 18,
            alpha=0.20 if state == "Inactive" else 0.65,
            linewidths=0, rasterized=True,
            label=state,
            zorder=4 if state == "Bivalent" else 1
        )

    z = np.polyfit(sub["bivalency_index"].values, sub[x_col].values, 1)
    xl = np.linspace(0, 20, 100)
    ax.plot(xl, np.polyval(z, xl), "k--", lw=1.5, alpha=0.7, zorder=5)

    p_str = f"{p:.2e}" if p >= 1e-4 else "< 10\u207b\u2074"
    ax.text(0.05, 0.95, f"\u03c1 = {rho:.3f}\np = {p_str}",
            transform=ax.transAxes, fontsize=9, va="top",
            bbox=dict(boxstyle="round,pad=0.25",
                      fc="white", ec="none", alpha=0.85))

    ax.axhline(0, color="#888888", lw=0.6, ls=":")
    ax.set_xlabel("Bivalency Index (H3K4me1 / H3K27me3, Parental)",
                  labelpad=6)
    ax.set_ylabel(xlabel, labelpad=6)
    ax.set_title(
        f"Parental Bivalency vs {xlabel.split('(')[0].strip()}",
        fontweight="bold", pad=10
    )
    ax.legend(frameon=False, fontsize=8, markerscale=2)
    save_fig(fig, outdir, stem)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df, cols = load_and_merge(args.cutrun_file, args.merged_file)
    df       = classify_chromatin_state(df, cols, args)
    df       = classify_fate(df, args)

    # Save for downstream R GO enrichment
    out_path = os.path.join(args.outdir, "Q3_chromatin_states.txt")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"[INFO] Saved -> {out_path}")

    print("\n[INFO] Generating Fig 3a: Chromatin state pie...")
    plot_state_pie(df, args.outdir)

    print("[INFO] Generating Fig 3b: Fate by state...")
    plot_fate_by_state(df, args.outdir)

    pc1_col = next((c for c in ["B10-P_Orig", "B10_P_Orig"]
                    if c in df.columns), None)
    if pc1_col:
        print("[INFO] Generating Fig 3c: Bivalency vs PC1...")
        plot_bivalency_scatter(
            df, pc1_col,
            "\u0394PC1 (Hi-C, Bottom10 \u2212 Parental)",
            "Fig3c_bivalency_vs_PC1", args.outdir
        )

    if "log2FC" in df.columns:
        print("[INFO] Generating Fig 3d: Bivalency vs RNA...")
        plot_bivalency_scatter(
            df, "log2FC",
            "log\u2082FC (RNA-seq, Bottom10 / Parental)",
            "Fig3d_bivalency_vs_RNA", args.outdir
        )

    print(f"\n[INFO] Q3 complete -> {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()