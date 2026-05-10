#!/lustre/isaac24/scratch/majagbe/envs/epp604_env/bin/python3

# =============================================================================
# Multi-Omics Integration: Publication-Quality Figures
# A375 Melanoma — CUT&RUN + RNA-seq + Hi-C Integration
#
# USAGE:
#   python multiomics_plots.py merged_data.txt
#   python multiomics_plots.py merged_data.txt --outdir ./figures
#   python multiomics_plots.py merged_data.csv --sep , --padj 0.05 --lfc 1.5
#
# POSITIONAL:
#   input_file    Path to merged data table (tab-delimited by default)
#
# OPTIONAL:
#   --outdir      Output directory (default: current working directory)
#   --sep         Column delimiter: "tab" or "," or ";"  (default: tab)
#   --padj        Adjusted p-value threshold (default: 0.05)
#   --lfc         log2 fold-change threshold (default: 2.0)
#   --me3         H3K27me3 |Bottom-Parent| threshold for Fig4 (default: 0.05)
# =============================================================================

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
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# ── Global publication aesthetics ─────────────────────────────────────────────
plt.rcParams.update({
    "font.family":        "DejaVu Sans",
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.linewidth":     1.2,
    "xtick.major.width":  1.2,
    "ytick.major.width":  1.2,
    "xtick.major.size":   4,
    "ytick.major.size":   4,
    "pdf.fonttype":       42,
    "ps.fonttype":        42,
    "savefig.dpi":        300,
    "figure.dpi":         150,
})


# =============================================================================
# CLI ARGUMENT PARSING
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate publication-quality multi-omics volcano and bar plots\n"
            "for A375 melanoma DEG + H3K27me3 CUT&RUN + Hi-C compartment data."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python multiomics_plots.py merged.txt\n"
            "  python multiomics_plots.py merged.txt --outdir ./figures\n"
            "  python multiomics_plots.py merged.csv --sep , --padj 0.05\n"
        )
    )

    parser.add_argument(
        "input_file",
        metavar="INPUT_FILE",
        help=(
            "Path to merged data table (tab-delimited by default). "
            "Required columns: gene_name, padj, log2FC, baseMean, "
            "b10_minus_parental_h3k27me3, P->B10_Switch_Orig, P->B10_Shift_Orig, "
            "P->B10_Both_original."
        )
    )
    parser.add_argument(
        "--outdir", "-o",
        metavar="DIR",
        default=os.getcwd(),
        help="Directory to save output figures (default: current directory)."
    )
    parser.add_argument(
        "--sep",
        metavar="DELIM",
        default="tab",
        help='Column delimiter: "tab" (default), "," or ";". '
             'Use "tab" for .txt / .tsv files.'
    )
    parser.add_argument(
        "--padj",
        type=float,
        default=0.01,
        help="Adjusted p-value threshold for DEG significance (default: 0.01)."
    )
    parser.add_argument(
        "--lfc",
        type=float,
        default=2.0,
        help="log2 fold-change threshold (default: 2.0)."
    )
    parser.add_argument(
        "--me3",
        type=float,
        default=0.05,
        help="Minimum |Bottom-Parent| H3K27me3 signal for Fig4 (default: 0.05)."
    )

    return parser.parse_args()


# =============================================================================
# DATA LOADING
# =============================================================================

def load_data(input_file, sep="tab"):
    delim = "\t" if sep == "tab" else sep

    if not os.path.isfile(input_file):
        sys.exit(f"[ERROR] File not found: {input_file}")

    try:
        df = pd.read_csv(input_file, sep=delim, low_memory=False)
    except Exception as e:
        sys.exit(f"[ERROR] Could not read {input_file}: {e}")

    required = [
        "gene_name", "padj", "log2FC", "baseMean", "b10_minus_parental_h3k27me3",
        "P->B10_Switch_Orig", "P->B10_Shift_Orig", "P->B10_Both_original"
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        sys.exit(
            f"[ERROR] Missing required columns: {missing}\n"
            f"Found columns: {df.columns.tolist()}"
        )

    for col in ["padj", "log2FC", "baseMean", "b10_minus_parental_h3k27me3"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    print(f"[INFO] Loaded {len(df):,} rows x {len(df.columns)} columns")
    return df


# =============================================================================
# SHARED UTILITIES
# =============================================================================

def get_sig_mask(df, padj_thresh, lfc_thresh):
    return (df["padj"] < padj_thresh) & (np.abs(df["log2FC"]) > lfc_thresh)


def get_top10_labels(df, sig_mask):
    sig_df = df[sig_mask].copy()
    sig_df["_lp"] = -np.log10(sig_df["padj"].replace(0, 1e-300))
    top_up   = sig_df[sig_df["log2FC"] > 0].nlargest(5, "_lp")
    top_down = sig_df[sig_df["log2FC"] < 0].nlargest(5, "_lp")
    return pd.concat([top_up, top_down])


def annotate_top_genes(ax, top_genes):
    for _, row in top_genes.iterrows():
        xi = row["log2FC"]
        yi = -np.log10(row["padj"]) if row["padj"] > 0 else 300
        x_off = 1.5 if xi > 0 else -1.5
        y_off = max(3, yi * 0.04)
        ax.annotate(
            row["gene_name"],
            xy=(xi, yi),
            xytext=(xi + x_off, yi + y_off),
            fontsize=8, fontweight="bold", color="black",
            arrowprops=dict(
                arrowstyle="-", color="#444444", lw=0.7,
                connectionstyle="arc3,rad=0.1"
            ),
            zorder=10,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8)
        )


def draw_threshold_lines(ax, padj_thresh, lfc_thresh):
    y_thresh = -np.log10(padj_thresh)
    ax.axhline(y_thresh,    color="#999999", lw=0.9, ls="--", zorder=0, alpha=0.7)
    ax.axvline( lfc_thresh, color="#999999", lw=0.9, ls="--", zorder=0, alpha=0.7)
    ax.axvline(-lfc_thresh, color="#999999", lw=0.9, ls="--", zorder=0, alpha=0.7)
    ax.text(lfc_thresh * 1.02, y_thresh * 0.97,
            f"padj={padj_thresh}", color="#999999", fontsize=7.5, va="top")


def classify_compartment(row):
    """Priority: Both_original > Switch > Shift > No change."""
    for col in ["P->B10_Both_original", "P->B10_Switch_Orig", "P->B10_Shift_Orig"]:
        val = row[col]
        if pd.notna(val) and str(val).strip() not in ["-", "", "nan"]:
            return str(val).strip()
    return "No change"


def broad_category(cat):
    if cat in ["BothtoA", "BtoA"]:   return "Both to A"
    if cat in ["TowardA"]:           return "Shift to A"
    if cat in ["BothtoB", "AtoB"]:   return "Both to B"
    if cat in ["TowardB"]:           return "Shift to B"
    return None


# ── Compartment color palettes ────────────────────────────────────────────────
COMP_COLORS = {
    "BothtoA":   "#C0392B",
    "BtoA":      "#E74C3C",
    "TowardA":   "#F1948A",
    "BothtoB":   "#1A5276",
    "AtoB":      "#2980B9",
    "TowardB":   "#85C1E9",
    "No change": "#CACFD2",
}

BROAD_COLORS = {
    "Both to A": "#C0392B",
    "Shift to A":  "#F1948A",
    "Both to B": "#1A5276",
    "Shift to B":  "#85C1E9",
}

CAT_ORDER = ["Both to A", "Shift to A", "Shift to B", "Both to B"]


# =============================================================================
# SAVE HELPER
# =============================================================================

def _save(fig, outdir, stem):
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        path = os.path.join(outdir, f"{stem}.{ext}")
        fig.savefig(path, bbox_inches="tight", dpi=300)
        print(f"  Saved -> {path}")
    plt.close(fig)


# =============================================================================
# FIGURE 1 — Classic Volcano: Up/Down DEGs + Top 10 Labels
# =============================================================================

def plot_fig1(df, args, outdir):
    print("\n" + "="*60)
    print("Figure 1: Volcano -- Up/Down DEGs")
    print("="*60)

    sig  = get_sig_mask(df, args.padj, args.lfc)
    up   = sig & (df["log2FC"] >  args.lfc)
    down = sig & (df["log2FC"] < -args.lfc)
    y_all = -np.log10(df["padj"].replace(0, 1e-300))

    fig, ax = plt.subplots(figsize=(7.5, 7))

    ax.scatter(df.loc[~sig, "log2FC"], y_all[~sig],
               c="#d0d0d0", s=10, alpha=0.20, linewidths=0,
               rasterized=True, zorder=1)
    ax.scatter(df.loc[up,   "log2FC"], y_all[up],
               c="#D62728", s=55, alpha=0.88, linewidths=0.3,
               edgecolors="white", rasterized=True, zorder=4,
               label=f"Upregulated (n={up.sum()})")
    ax.scatter(df.loc[down, "log2FC"], y_all[down],
               c="#1F77B4", s=55, alpha=0.88, linewidths=0.3,
               edgecolors="white", rasterized=True, zorder=4,
               label=f"Downregulated (n={down.sum()})")

    draw_threshold_lines(ax, args.padj, args.lfc)
    annotate_top_genes(ax, get_top10_labels(df, sig))

    #ax.text(0.97, 0.97, f"up {up.sum()} genes",
    #       transform=ax.transAxes, ha="right", va="top",
    #        color="#D62728", fontsize=11, fontweight="bold")
    #ax.text(0.03, 0.97, f"down {down.sum()} genes",
    #        transform=ax.transAxes, ha="left", va="top",
    #        color="#1F77B4", fontsize=11, fontweight="bold")

    legend_elements = [
        mpatches.Patch(fc="#D62728", label=f"Upregulated (n={up.sum()})"),
        mpatches.Patch(fc="#1F77B4", label=f"Downregulated (n={down.sum()})"),
        mpatches.Patch(fc="#d0d0d0", label="Not significant"),
    ]
    ax.legend(handles=legend_elements, frameon=False, fontsize=9,
              loc="center", bbox_to_anchor=(0.5, 0.85))

    ax.set_xlabel("log$_2$ Fold Change (Bottom10 / Parental)", fontsize=12, labelpad=6)
    ax.set_ylabel("-log$_{10}$(adjusted $p$-value)", fontsize=12, labelpad=6)
    ax.set_title(
        "Differentially Expressed Genes\nParental vs Bottom10 (A375 Melanoma)",
        fontsize=13, fontweight="bold", pad=12
    )
    _save(fig, outdir, "Fig1_Volcano_UpDown")


# =============================================================================
# FIGURE 2 — Volcano Colored by H3K27me3 Signal Change
# =============================================================================

def plot_fig2(df, args, outdir):
    print("\n" + "="*60)
    print("Figure 2: Volcano -- H3K27me3 Signal Change")
    print("="*60)

    sig    = get_sig_mask(df, args.padj, args.lfc)
    sig_df = df[sig].copy()
    y_all  = -np.log10(df["padj"].replace(0, 1e-300))

    signal_vals = sig_df["b10_minus_parental_h3k27me3"].values
    vmax = max(abs(np.nanmin(signal_vals)), abs(np.nanmax(signal_vals)))
    vmin = -vmax
    cmap = plt.cm.RdBu_r
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm   = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    fig, ax = plt.subplots(figsize=(8, 7))

    ax.scatter(df.loc[~sig, "log2FC"], y_all[~sig],
               c="#d8d8d8", s=10, alpha=0.20, linewidths=0,
               rasterized=True, zorder=1)
    draw_threshold_lines(ax, args.padj, args.lfc)
    ax.scatter(sig_df["log2FC"], y_all[sig],
               c=signal_vals, cmap=cmap, vmin=vmin, vmax=vmax,
               s=65, alpha=0.90, linewidths=0.4,
               edgecolors="white", zorder=4)

    cbar = plt.colorbar(sm, ax=ax, shrink=0.50, pad=0.02, aspect=22)
    cbar.set_label("H3K27me3 Signal\n(Bottom10 - Parental)", fontsize=9, labelpad=8)
    cbar.ax.tick_params(labelsize=8)
    cbar.ax.axhline(0, color="black", lw=1.2)

    q_labels = [
        (sig & (df["log2FC"] >  args.lfc) & (df["b10_minus_parental_h3k27me3"] >  args.me3),
         "#D62728", "Up + H3K27me3 gain"),
        (sig & (df["log2FC"] >  args.lfc) & (df["b10_minus_parental_h3k27me3"] < -args.me3),
         "#1F77B4", "Up + H3K27me3 loss"),
        (sig & (df["log2FC"] < -args.lfc) & (df["b10_minus_parental_h3k27me3"] >  args.me3),
         "#E6550D", "Down + H3K27me3 gain"),
        (sig & (df["log2FC"] < -args.lfc) & (df["b10_minus_parental_h3k27me3"] < -args.me3),
         "#9ECAE1", "Down + H3K27me3 loss"),
    ]
    for j, (mask, color, label) in enumerate(q_labels):
        ax.text(0.3, 0.9 - j * 0.03, f"{label}: n={mask.sum()}",
                transform=ax.transAxes, fontsize=8,
                color=color, fontweight="bold")

    annotate_top_genes(ax, get_top10_labels(df, sig))

    ax.set_xlabel("log$_2$ Fold Change (Bottom10 / Parental)", fontsize=12, labelpad=6)
    ax.set_ylabel("-log$_{10}$(adjusted $p$-value)", fontsize=12, labelpad=6)
    ax.set_title(
        "DEGs Colored by H3K27me3 Signal Change\n(Bottom10 - Parental)",
        fontsize=13, fontweight="bold", pad=12
    )
    _save(fig, outdir, "Fig2_Volcano_H3K27me3Signal")


# =============================================================================
# FIGURE 3 — Volcano Colored by Hi-C Compartment Change
# =============================================================================

def plot_fig3(df, args, outdir):
    print("\n" + "="*60)
    print("Figure 3: Volcano -- Hi-C Compartment Change")
    print("="*60)

    df = df.copy()
    df["_comp_cat"] = df.apply(classify_compartment, axis=1)

    LABEL_MAP = {
        "BothtoA":   "Both -> A (Switch + Shift)",
        "BtoA":      "Switch B->A",
        "TowardA":   "Shift Toward A",
        "BothtoB":   "Both -> B (Switch + Shift)",
        "AtoB":      "Switch A->B",
        "TowardB":   "Shift Toward B",
        "No change": "No compartment change",
    }
    DRAW_ORDER = ["No change","TowardB","AtoB","BothtoB","TowardA","BtoA","BothtoA"]

    sig    = get_sig_mask(df, args.padj, args.lfc)
    sig_df = df[sig].copy()
    sig_df["_y"] = -np.log10(sig_df["padj"].replace(0, 1e-300))
    y_all  = -np.log10(df["padj"].replace(0, 1e-300))

    fig, ax = plt.subplots(figsize=(8.5, 7))

    ax.scatter(df.loc[~sig, "log2FC"], y_all[~sig],
               c="#e0e0e0", s=9, alpha=0.18, linewidths=0,
               rasterized=True, zorder=1)
    draw_threshold_lines(ax, args.padj, args.lfc)

    for cat in DRAW_ORDER:
        sub = sig_df[sig_df["_comp_cat"] == cat]
        if len(sub) == 0:
            continue
        is_nochange = (cat == "No change")
        ax.scatter(
            sub["log2FC"], sub["_y"],
            c=COMP_COLORS.get(cat, "#CACFD2"),
            s=18 if is_nochange else 60,
            alpha=0.30 if is_nochange else 0.90,
            linewidths=0 if is_nochange else 0.4,
            edgecolors="white",
            zorder=2 if is_nochange else 5,
            label=f"{LABEL_MAP.get(cat, cat)} (n={len(sub)})"
        )

    annotate_top_genes(ax, get_top10_labels(df, sig))

    handles, labels = ax.get_legend_handles_labels()
    order_map = {
        f"{LABEL_MAP.get(k,k)} (n={len(sig_df[sig_df['_comp_cat']==k])})": i
        for i, k in enumerate(DRAW_ORDER)
    }
    pairs = sorted(zip(labels, handles), key=lambda x: order_map.get(x[0], 99))
    if pairs:
        labels_s, handles_s = zip(*pairs)
        ax.legend(handles_s, labels_s, frameon=False, fontsize=8,
                  loc="center", bbox_to_anchor=(0.5, 0.84),
                  title="Compartment Change\n(P->B10, Original)",
                  title_fontsize=9)

    ax.set_xlabel("log$_2$ Fold Change (Bottom10 / Parental)", fontsize=12, labelpad=6)
    ax.set_ylabel("-log$_{10}$(adjusted $p$-value)", fontsize=12, labelpad=6)
    ax.set_title(
        "DEGs Colored by Hi-C Compartment Change\n(Parental -> Bottom10)",
        fontsize=13, fontweight="bold", pad=12
    )
    _save(fig, outdir, "Fig3_Volcano_CompartmentChange")
    return df


# =============================================================================
# FIGURE 4 — Side Bar: DEGs + Compartment Change + H3K27me3 Change
# =============================================================================

def plot_fig4(df, args, outdir):
    print("\n" + "="*60)
    print("Figure 4: Side Bar -- DEGs + Compartment + H3K27me3")
    print("="*60)

    df = df.copy()
    df["_comp_cat"]   = df.apply(classify_compartment, axis=1)
    df["_broad_comp"] = df["_comp_cat"].apply(broad_category)

    sig      = get_sig_mask(df, args.padj, args.lfc)
    has_comp = df["_broad_comp"].notna()
    has_me3  = df["b10_minus_parental_h3k27me3"].abs() > args.me3

    subset = df[sig & has_comp & has_me3].copy()
    print(f"[INFO] Genes passing all filters "
          f"(padj<{args.padj}, |lfc|>{args.lfc}, |dH3K27me3|>{args.me3}): "
          f"{len(subset)}")

    if len(subset) == 0:
        print("[WARNING] No genes pass all filters. Try lowering --me3.")
        return

    # ── Composite importance score ─────────────────────────────────────────────
    subset["_lp"]    = -np.log10(subset["padj"].replace(0, 1e-300))
    subset["_score"] = subset["b10_minus_parental_h3k27me3"].abs() * subset["_lp"]

    def _top_n(data, n):
        """
        Select top-n rows per _broad_comp by _score.
        Uses rank() to avoid the pandas groupby column-drop bug.
        """
        data = data.copy()
        data["_rank"] = data.groupby("_broad_comp")["_score"].rank(
            method="first", ascending=False
        )
        return data[data["_rank"] <= n].drop(columns="_rank").reset_index(drop=True)

    # Up to 5 per category initially
    top_df = _top_n(subset, 7)

    # If total < 10, relax per-category cap to reach at least 10
    if len(top_df) < 10:
        n_cats  = max(1, subset["_broad_comp"].nunique())
        per_cat = -(-10 // n_cats)          # ceiling division
        print(f"[INFO] < 10 genes; relaxing cap to {per_cat} per category.")
        top_df  = _top_n(subset, per_cat)

    # Final sort and hard cap at 20
    top_df["_cat_order"] = top_df["_broad_comp"].map(
        {c: i for i, c in enumerate(CAT_ORDER)}
    )
    top_df = (top_df
              .sort_values(["_cat_order", "b10_minus_parental_h3k27me3"], ascending=[True, True])
              .drop_duplicates(subset="gene_name")
              )

    print(f"[INFO] Final genes selected: {len(top_df)}")
    print(top_df[["gene_name", "_broad_comp", "b10_minus_parental_h3k27me3",
                  "log2FC", "padj"]].to_string(index=False))

    n_genes  = len(top_df)
    y_pos    = np.arange(n_genes)
    signals  = top_df["b10_minus_parental_h3k27me3"].values
    colors   = [BROAD_COLORS[c] for c in top_df["_broad_comp"]]
    gene_lbl = top_df["gene_name"].values
    lfc_vals = top_df["log2FC"].values
    lp_vals  = top_df["_lp"].values

    # ── Normalise dot sizes (40–200 range) ────────────────────────────────────
    lp_min, lp_max = lp_vals.min(), lp_vals.max()
    dot_sizes = 40 + (lp_vals - lp_min) / (lp_max - lp_min + 1e-9) * 160

    # ── Figure ────────────────────────────────────────────────────────────────
    fig_h = max(5, n_genes * 0.4)
    fig, ax = plt.subplots(figsize=(9.5, fig_h))

    # Horizontal bars
    ax.barh(y_pos, signals,
            color=colors, height=0.65,
            edgecolor="white", linewidth=0.5, zorder=3)

    # ── x-range for positioning dots and triangles ────────────────────────────
    x_span   = max(abs(np.nanmin(signals)), abs(np.nanmax(signals)))
    if x_span == 0:
        x_span = 0.1

    dot_gap = x_span * 0.07          # gap between bar tip and dot centre
    tri_gap = dot_gap * 2.4          # triangle sits beyond the dot

    for i, (sig_val, ds, col, lfc) in enumerate(
            zip(signals, dot_sizes, colors, lfc_vals)):

        direction = 1 if sig_val >= 0 else -1

        # Dot — significance marker
        x_dot = sig_val + direction * dot_gap
        ax.scatter(x_dot, i, s=ds, color=col,
                   edgecolors="black", linewidths=0.5,
                   zorder=5, alpha=0.90, clip_on=True)

        # Triangle — RNA expression direction, placed just after dot
        x_tri   = sig_val + direction * tri_gap
        symbol  = "▲" if lfc > 0 else "▼"
        sym_col = "#8B0000" if lfc > 0 else "#00008B"
        ax.text(x_tri, i, symbol,
                fontsize=8, va="center",
                ha="left" if direction > 0 else "right",
                color=sym_col, fontweight="bold", zorder=6)

    # ── Gene name labels ──────────────────────────────────────────────────────
    ax.set_yticks(y_pos)
    ax.set_yticklabels(gene_lbl, fontsize=9, fontweight="bold")

    # Zero line
    ax.axvline(0, color="black", lw=1.3, zorder=4)

    # ── Auto x-limits with breathing room for triangles ──────────────────────
    max_abs = np.nanmax(np.abs(signals))
    pad     = max_abs * 0.35        # 35 % padding either side
    ax.set_xlim(-max_abs - pad, max_abs + pad)

    # ── Category bracket annotations ─────────────────────────────────────────
    prev_cat, cat_start = None, 0
    boundaries = []
    for i, cat in enumerate(top_df["_broad_comp"].values):
        if cat != prev_cat:
            if prev_cat is not None:
                boundaries.append((prev_cat, cat_start, i - 1))
            cat_start = i
            prev_cat  = cat
    boundaries.append((prev_cat, cat_start, n_genes - 1))

    xlim      = ax.get_xlim()
    x_bracket = xlim[1] + (xlim[1] - xlim[0]) * 0.04

    for cat, start, end in boundaries:
        mid = (start + end) / 2
        col = BROAD_COLORS[cat]
        ax.annotate("",
                    xy=(x_bracket, end + 0.38),
                    xytext=(x_bracket, start - 0.38),
                    arrowprops=dict(arrowstyle="-", color=col, lw=2.5),
                    annotation_clip=False)
        ax.text(x_bracket + (xlim[1] - xlim[0]) * 0.015, mid,
                cat, fontsize=8.5, fontweight="bold",
                color=col, va="center", ha="left",
                clip_on=False)

    # ── Directional headers ───────────────────────────────────────────────────
    ax.text(0.01, 0.99,
            "<-- H3K27me3 Loss in Bottom10",
            transform=ax.transAxes, fontsize=10,
            color="#C0392B", style="italic", ha="left")
    ax.text(0.99, 0.99,
            "H3K27me3 Gain in Bottom10 -->",
            transform=ax.transAxes, fontsize=10,
            color="#1A5276", style="italic", ha="right")

    # Grid
    ax.grid(axis="x", lw=0.5, color="#ebebeb", zorder=0)
    ax.set_axisbelow(True)

# ── Legend ────────────────────────────────────────────────────────────────
    present_cats = [c for c in CAT_ORDER if c in top_df["_broad_comp"].values]
    legend_patches = [
        mpatches.Patch(fc=BROAD_COLORS[c], label=c) for c in present_cats
    ]

# ── Legend ────────────────────────────────────────────────────────────────
    present_cats = [c for c in CAT_ORDER if c in top_df["_broad_comp"].values]
    legend_patches = [
        mpatches.Patch(fc=BROAD_COLORS[c], label=c) for c in present_cats
    ]

    # Fixed graduated dot sizes for padj significance levels in legend
    legend_patches += [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor="grey", markeredgecolor="black",
               markersize=ms, label=label)
        for label, ms in [
            ("padj < 0.01",   6),
            ("padj < 0.001",  8),
            ("padj < 0.0001", 10),
        ]
    ]

    legend_patches += [
        Line2D([0], [0], marker="^", color="w",
               markerfacecolor="#8B0000", markeredgecolor="#8B0000",
               markersize=6, label="Upregulated"),
        Line2D([0], [0], marker="v", color="w",
               markerfacecolor="#00008B", markeredgecolor="#00008B",
               markersize=6, label="Downregulated"),
    ]
    
    ax.legend(handles=legend_patches, frameon=True, framealpha=0.9,
              edgecolor="#cccccc", fontsize=8,
              loc="lower left", bbox_to_anchor=(0.1, 0.1),
              title="Compartment Category", title_fontsize=8.5)

    ax.set_xlabel("H3K27me3 Signal Change",
                  fontsize=11, labelpad=8)
    ax.set_title(
        "Top DEGs with H3K27me3 Change and Compartment Reorganization\n"
        "A375 Melanoma: Parental -> Bottom10",
        fontsize=12, fontweight="bold", pad=20
    )

    ax.invert_yaxis()
    _save(fig, outdir, "Fig4_SideBar_DEGs_Compartment_H3K27me3")


# =============================================================================
# MAIN
# =============================================================================

def main():
    args = parse_args()

    print(f"\n{'='*60}")
    print("  Multi-Omics Visualization -- A375 Melanoma")
    print(f"{'='*60}")
    print(f"  Input     : {args.input_file}")
    print(f"  Output    : {args.outdir}")
    print(f"  padj      : {args.padj}")
    print(f"  |lfc|     : {args.lfc}")
    print(f"  |dH3K27me3|: {args.me3}")
    print(f"{'='*60}\n")

    df = load_data(args.input_file, sep=args.sep)

    plot_fig1(df, args, args.outdir)
    plot_fig2(df, args, args.outdir)
    df = plot_fig3(df, args, args.outdir)
    plot_fig4(df, args, args.outdir)

    print(f"\n{'='*60}")
    print("  All figures saved successfully.")
    print(f"  Output: {os.path.abspath(args.outdir)}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()