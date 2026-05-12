#!/usr/bin/env python3
"""
Flexible trimodal (or multi-modal) coupling analysis

Hi-C compartment is always included by default.
You must specify at least one of --rna or --marks.

Usage examples:
    # RNA + H3K27me3 + Hi-C (default thresholds)
    python scripts/trimodal_coupling.py merged.txt \
        --rna --marks h3k27me3

    # RNA + all four marks + Hi-C
    python scripts/trimodal_coupling.py merged.txt \
        --rna --marks h3k27me3 h3k27ac h3k4me1 h3k9me3

    # H3K27me3 + H3K9me3 + Hi-C (no RNA)
    python scripts/trimodal_coupling.py merged.txt \
        --marks h3k27me3 h3k9me3

    # Disable Hi-C
    python scripts/trimodal_coupling.py merged.txt \
        --rna --marks h3k27me3 --no_hic

Available mark keys (case-insensitive):
    h3k27me3   -> b10_minus_parental_h3k27me3
    h3k27ac    -> b10_minus_parental_h3k27ac
    h3k4me1    -> b10_minus_parental_h3k4me1
    h3k9me3    -> b10_minus_parental_h3k9me3
    unc        -> b10_minus_parental_unc_h3k27me3
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
from itertools import combinations
from scipy.stats import spearmanr

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

# ── Mark registry: short key -> (delta column, display label) ─────────────────
MARK_REGISTRY = {
    "h3k27me3": ("b10_minus_parental_h3k27me3",    "H3K27me3"),
    "h3k27ac":  ("b10_minus_parental_h3k27ac",     "H3K27Ac"),
    "h3k4me1":  ("b10_minus_parental_h3k4me1",     "H3K4me1"),
    "h3k9me3":  ("b10_minus_parental_h3k9me3",     "H3K9me3"),
    "unc":      ("b10_minus_parental_unc_h3k27me3", "UNC-H3K27me3"),
}

# Modality colors — up to 7 modalities without repeating
MODALITY_COLORS = [
    "#1A5276",   # Hi-C (always first)
    "#C0392B",   # RNA
    "#8E44AD",   # H3K27me3
    "#16A085",   # H3K27Ac
    "#F39C12",   # H3K4me1
    "#27AE60",   # H3K9me3
    "#E74C3C",   # UNC
]


def parse_args():
    p = argparse.ArgumentParser(
        description="Flexible multi-modal coupling analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    p.add_argument("merged_file", help="merged multi-omics table")
    p.add_argument("--outdir",      default="outputs/figures/Q2_coupling")

    # ── Modality selection ────────────────────────────────────────────────────
    p.add_argument("--rna",    action="store_true",
                   help="Include RNA-seq transcription as a modality")
    p.add_argument("--marks",  nargs="+",
                   metavar="MARK",
                   help="Epigenetic marks to include. Choose from: "
                        "h3k27me3 h3k27ac h3k4me1 h3k9me3 unc")
    p.add_argument("--no_hic", action="store_true",
                   help="Exclude Hi-C compartment (included by default)")

    # ── Thresholds ────────────────────────────────────────────────────────────
    p.add_argument("--lfc_thresh",  type=float, default=2.0,
                   help="RNA: |log2FC| threshold (default 2.0)")
    p.add_argument("--padj_thresh", type=float, default=0.01,
                   help="RNA: padj threshold (default 0.01)")
    p.add_argument("--me3_thresh",  type=float, default=0.05,
                   help="Epigenetic mark: |delta| threshold (default 0.05)")

    args = p.parse_args()

    # Validate: must have at least one modality besides Hi-C
    if not args.rna and not args.marks:
        p.error("Must specify at least one of --rna or --marks. "
                "Hi-C is included automatically unless --no_hic is set.")

    # Validate mark keys
    if args.marks:
        bad = [m for m in args.marks if m.lower() not in MARK_REGISTRY]
        if bad:
            p.error(f"Unknown mark(s): {bad}. "
                    f"Choose from: {list(MARK_REGISTRY.keys())}")
        args.marks = [m.lower() for m in args.marks]

    return args


def save_fig(fig, outdir, stem):
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        path = os.path.join(outdir, f"{stem}.{ext}")
        fig.savefig(path, bbox_inches="tight", pad_inches=0.05, dpi=600)
        print(f"  Saved -> {path}")
    plt.close(fig)


# ── Categorical Hi-C label checker ────────────────────────────────────────────
def has_comp_label(row):
    for col in ["P->B10_Both_original",
                "P->B10_Switch_Orig",
                "P->B10_Shift_Orig"]:
        v = str(row.get(col, "")).strip()
        if v not in ["-", "", "nan", "None"]:
            return True
    return False


def get_comp_dir(row):
    for col in ["P->B10_Both_original",
                "P->B10_Switch_Orig",
                "P->B10_Shift_Orig"]:
        v = str(row.get(col, "")).strip()
        if v not in ["-", "", "nan", "None"]:
            if "toA" in v or v in ["BtoA", "TowardA"]:
                return "B\u2192A"
            if "toB" in v or v in ["AtoB", "TowardB"]:
                return "A\u2192B"
    return "Stable"


def build_modalities(df, args):
    """
    Build a dict of:
        modality_label -> {
            "changed_col":  boolean column name in df,
            "delta_col":    continuous delta column (for Spearman),
            "color":        plot color,
        }

    Order: Hi-C first (if enabled), then RNA, then marks in order given.
    """
    modalities = {}
    color_idx = 0

    # 1. Hi-C compartment (always first, always categorical)
    if not args.no_hic:
        df["_hic_changed"] = df.apply(has_comp_label, axis=1)
        # For Spearman, use the continuous PC1 delta
        pc1_col = next((c for c in ["B10-P_Orig", "B10_P_Orig"]
                        if c in df.columns), None)
        if pc1_col:
            df[pc1_col] = pd.to_numeric(df[pc1_col], errors="coerce")
        modalities["Hi-C Compartment"] = {
            "changed_col": "_hic_changed",
            "delta_col":   pc1_col,
            "color":       MODALITY_COLORS[color_idx],
        }
        color_idx += 1

    # 2. RNA
    if args.rna:
        if "log2FC" not in df.columns or "padj" not in df.columns:
            print("[WARN] --rna requested but log2FC/padj not found. Skipping RNA.")
        else:
            df["log2FC"] = pd.to_numeric(df["log2FC"], errors="coerce")
            df["padj"]   = pd.to_numeric(df["padj"],   errors="coerce")
            df["_rna_changed"] = (
                (df["log2FC"].abs() > args.lfc_thresh) &
                (df["padj"] < args.padj_thresh)
            ).fillna(False)
            modalities["Transcription"] = {
                "changed_col": "_rna_changed",
                "delta_col":   "log2FC",
                "color":       MODALITY_COLORS[color_idx % len(MODALITY_COLORS)],
            }
            color_idx += 1

    # 3. Epigenetic marks
    if args.marks:
        for mark_key in args.marks:
            delta_col, label = MARK_REGISTRY[mark_key]
            if delta_col not in df.columns:
                print(f"[WARN] Column '{delta_col}' not found for mark "
                      f"'{mark_key}'. Skipping.")
                continue
            df[delta_col] = pd.to_numeric(df[delta_col], errors="coerce")
            changed_col = f"_{mark_key}_changed"
            df[changed_col] = (
                df[delta_col].abs() > args.me3_thresh
            ).fillna(False)
            modalities[label] = {
                "changed_col": changed_col,
                "delta_col":   delta_col,
                "color":       MODALITY_COLORS[color_idx % len(MODALITY_COLORS)],
            }
            color_idx += 1

    if len(modalities) < 2:
        sys.exit("[ERROR] Need at least 2 modalities to run coupling analysis.")

    # Print counts
    print("\n[INFO] Modality change counts:")
    for name, info in modalities.items():
        n = df[info["changed_col"]].sum()
        print(f"  {name}: {n:,} genes changed")

    return modalities


def make_coupling_label(row, modalities):
    parts = [name for name, info in modalities.items()
             if row.get(info["changed_col"], False)]
    return "+".join(parts) if parts else "None"


# ══════════════════════════════════════════════════════════════════════════════
# PLOT A: UpSet — manual matplotlib implementation (version-safe)
# ══════════════════════════════════════════════════════════════════════════════
def plot_upset(df, modalities, outdir):
    """
    Manual UpSet plot. Avoids upsetplot library color/nan crashes.
    Works for any number of modalities.
    """
    mod_names = list(modalities.keys())
    n_mods    = len(mod_names)

    # Build all 2^n intersection counts
    intersections = {}
    for size in range(n_mods + 1):
        for combo in combinations(range(n_mods), size):
            label = tuple(sorted([mod_names[i] for i in combo]))
            mask  = pd.Series([True] * len(df), index=df.index)
            for j, info in enumerate(modalities.values()):
                col = info["changed_col"]
                if j in combo:
                    mask = mask & df[col].fillna(False)
                else:
                    mask = mask & ~df[col].fillna(False)
            intersections[label] = int(mask.sum())

    # Sort by count descending
    sorted_items = sorted(intersections.items(),
                          key=lambda x: x[1], reverse=True)
    labels = [item[0] for item in sorted_items]
    counts = [item[1] for item in sorted_items]
    n_bars = len(labels)

    all_mods_set = frozenset(mod_names)

    # Determine if a bar has ALL modalities (highlight red)
    def is_all(label):
        return set(label) == all_mods_set

    # Figure layout
    fig = plt.figure(figsize=(max(10, n_bars * 0.9), 6 + n_mods * 0.4))
    ax_bar = fig.add_axes([0.22, 0.45, 0.74, 0.48])
    ax_dot = fig.add_axes([0.22, 0.05, 0.74, 0.38])
    ax_set = fig.add_axes([0.02, 0.05, 0.18, 0.38])

    x_pos = np.arange(n_bars)

    # ── Bar chart ──────────────────────────────────────────────────────────────
    bar_colors = ["#C0392B" if is_all(lab) else "#222222" for lab in labels]
    bars = ax_bar.bar(x_pos, counts, color=bar_colors,
                      edgecolor="none", width=0.6)
    for bar, val in zip(bars, counts):
        if val > 0:
            ax_bar.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + max(counts) * 0.01,
                        f"{val:,}", ha="center", va="bottom", fontsize=7.5,
                        fontweight="bold")
    ax_bar.set_xlim(-0.5, n_bars - 0.5)
    ax_bar.set_xticks([])
    ax_bar.set_ylabel("Intersection size", fontsize=11)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    if is_all(labels[0]):
        ax_bar.add_patch(mpatches.Patch(
            fc="#C0392B", ec="none",
            label=f"All {n_mods} co-changed (n={counts[0]:,})"
        ))
        ax_bar.legend(frameon=True, fontsize=9, loc="upper right")

    # ── Dot matrix ────────────────────────────────────────────────────────────
    y_positions = {name: i for i, name in enumerate(mod_names)}

    ax_dot.set_xlim(-0.5, n_bars - 0.5)
    ax_dot.set_ylim(-0.5, n_mods - 0.5)
    ax_dot.set_xticks([])
    ax_dot.set_yticks(range(n_mods))
    ax_dot.set_yticklabels(mod_names, fontsize=10)
    ax_dot.spines["top"].set_visible(False)
    ax_dot.spines["right"].set_visible(False)
    ax_dot.spines["bottom"].set_visible(False)

    for i in range(n_mods):
        if i % 2 == 0:
            ax_dot.axhspan(i - 0.5, i + 0.5, color="#F5F5F5", zorder=0)

    for xi, lab in enumerate(labels):
        active     = set(lab)
        all_active = is_all(lab)
        dot_color  = "#C0392B" if all_active else "black"
        active_ys  = sorted([y_positions[m] for m in active])

        for mod in mod_names:
            yi = y_positions[mod]
            if mod in active:
                ax_dot.scatter(xi, yi, s=110, color=dot_color,
                               zorder=3, linewidths=0)
            else:
                ax_dot.scatter(xi, yi, s=110, color="#DDDDDD",
                               zorder=2, linewidths=0)

        if len(active_ys) > 1:
            ax_dot.plot([xi, xi],
                        [min(active_ys), max(active_ys)],
                        color=dot_color, lw=2.5, zorder=2)

    # ── Set size bars ──────────────────────────────────────────────────────────
    set_sizes = {
        name: int(df[info["changed_col"]].fillna(False).sum())
        for name, info in modalities.items()
    }
    max_ss = max(set_sizes.values()) if set_sizes else 1

    ax_set.barh(
        [y_positions[m] for m in mod_names],
        [set_sizes[m] for m in mod_names],
        color=[modalities[m]["color"] for m in mod_names],
        height=0.5
    )
    for mod in mod_names:
        yi = y_positions[mod]
        ax_set.text(set_sizes[mod] + max_ss * 0.02,
                    yi, f"{set_sizes[mod]:,}",
                    va="center", fontsize=9)

    ax_set.set_ylim(-0.5, n_mods - 0.5)
    ax_set.invert_xaxis()
    ax_set.set_xlabel("Set size", fontsize=10)
    ax_set.set_yticks([])
    ax_set.spines["top"].set_visible(False)
    ax_set.spines["right"].set_visible(False)

    mod_str = " + ".join(mod_names)
    fig.suptitle(
        f"Co-occurrence of Multi-omics Changes\n"
        f"Parental \u2192 Bottom10  ({mod_str})",
        fontweight="bold", fontsize=12, y=1.01
    )
    save_fig(fig, outdir, "Fig2a_upset_cooccurrence")


# ══════════════════════════════════════════════════════════════════════════════
# PLOT B: Pairwise Spearman correlations
# ══════════════════════════════════════════════════════════════════════════════
def plot_spearman(df, modalities, outdir):
    """
    All pairwise Spearman correlations between continuous delta columns.
    One scatter panel per pair. Skips pairs where delta column is missing.
    """
    # Collect (label, delta_col, color) for modalities that have a delta col
    mod_list = [
        (name, info["delta_col"], info["color"])
        for name, info in modalities.items()
        if info["delta_col"] and info["delta_col"] in df.columns
    ]

    pairs = list(combinations(mod_list, 2))

    if not pairs:
        print("[WARN] No valid column pairs for Spearman. Skipping.")
        return

    n_pairs = len(pairs)
    # Best grid layout: max 3 columns, compute rows needed
    n_cols = min(n_pairs, 3)
    n_rows = int(np.ceil(n_pairs / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(4.5 * n_cols, 4.2 * n_rows),
                              squeeze=False)
    axes = axes.flatten()
    # Hide any unused panels
    for idx in range(n_pairs, len(axes)):
        axes[idx].set_visible(False)

    for ax, ((xname, xcol, xcolor), (yname, ycol, ycolor)) in zip(axes, pairs):
        sub = df[[xcol, ycol]].dropna()
        if len(sub) < 10:
            ax.set_visible(False)
            continue

        rho, p = spearmanr(sub[xcol], sub[ycol])

        # Use blended color between the two modalities
        ax.scatter(sub[xcol], sub[ycol],
                   c=xcolor, s=5, alpha=0.25,
                   linewidths=0, rasterized=True)

        z  = np.polyfit(sub[xcol].values, sub[ycol].values, 1)
        xl = np.linspace(sub[xcol].min(), sub[xcol].max(), 100)
        ax.plot(xl, np.polyval(z, xl),
                color="black", lw=1.2, ls="--", alpha=0.7)

        p_str = f"{p:.2e}" if p >= 1e-4 else "< 10\u207b\u2074"
        ax.text(0.05, 0.95,
                f"\u03c1 = {rho:.3f}\np = {p_str}",
                transform=ax.transAxes, fontsize=9, va="top",
                bbox=dict(boxstyle="round,pad=0.25",
                          fc="white", ec="none", alpha=0.85))

        ax.axhline(0, color="#888888", lw=0.6, ls=":")
        ax.axvline(0, color="#888888", lw=0.6, ls=":")

        # Axis labels: use delta symbol + modality name
        def delta_label(name, col):
            if "log" in col.lower() or "fc" in col.lower():
                return f"log\u2082FC ({name})"
            return f"\u0394{name} (B10 \u2212 Parental)"

        ax.set_xlabel(delta_label(xname, xcol), labelpad=5)
        ax.set_ylabel(delta_label(yname, ycol), labelpad=5)

    fig.suptitle(
        "Pairwise Spearman Correlations: " +
        " \u2194 ".join([m[0] for m in mod_list]),
        fontweight="bold", fontsize=12, y=1.02
    )
    plt.tight_layout()
    save_fig(fig, outdir, "Fig2b_spearman_correlations")


# ══════════════════════════════════════════════════════════════════════════════
# PLOT C: Directional concordance bar
# ══════════════════════════════════════════════════════════════════════════════
def plot_concordance_bar(df, modalities, outdir):
    """
    For genes with a Hi-C compartment call (B->A, A->B, Stable),
    compute the fraction that show concordant changes in each other modality.

    Concordance definition:
        B->A: other modality shows activation/loss direction
        A->B: other modality shows repression/gain direction

    Only runs if Hi-C is in modalities.
    """
    if "Hi-C Compartment" not in modalities:
        print("[INFO] Hi-C not included — skipping concordance bar.")
        return

    df["comp_dir"] = df.apply(get_comp_dir, axis=1)

    # Other modalities (excluding Hi-C)
    other_mods = {
        name: info for name, info in modalities.items()
        if name != "Hi-C Compartment"
    }

    if not other_mods:
        print("[INFO] No other modalities to compare — skipping concordance bar.")
        return

    comp_dirs = ["B\u2192A", "A\u2192B", "Stable"]
    x         = np.arange(len(comp_dirs))
    n_others  = len(other_mods)
    width     = 0.7 / n_others
    offsets   = np.linspace(-(n_others - 1) * width / 2,
                             (n_others - 1) * width / 2,
                             n_others)

    fig, ax = plt.subplots(figsize=(max(7, (n_others + 1) * 2.5), 5))

    # Track per-gene concordance flags for each modality — needed for "all concordant" bar
    concordance_flags = {}

    for (mod_name, info), offset in zip(other_mods.items(), offsets):
        delta_col = info["delta_col"]
        color     = info["color"]
        flags_per_dir = {}

        fracs = []
        for gd in comp_dirs:
            sub = df[df["comp_dir"] == gd]
            n   = len(sub)
            if n == 0 or not delta_col or delta_col not in df.columns:
                fracs.append(0)
                flags_per_dir[gd] = pd.Series(
                    [False] * len(sub), index=sub.index)
                continue

            vals = pd.to_numeric(sub[delta_col], errors="coerce")
            repressive = any(k in delta_col.lower()
                             for k in ["me3", "h3k9", "unc"])

            if gd == "B\u2192A":
                concordant = vals < 0 if repressive else vals > 0
            elif gd == "A\u2192B":
                concordant = vals > 0 if repressive else vals < 0
            else:
                concordant = pd.Series([False] * len(sub), index=sub.index)

            concordant = concordant.fillna(False)
            flags_per_dir[gd] = concordant
            fracs.append(concordant.sum() / n * 100)

        concordance_flags[mod_name] = flags_per_dir

        ax.bar(x + offset, fracs, width=width * 0.9,
               color=color, label=mod_name,
               edgecolor="white", linewidth=0.4)

    # "All concordant" bar — gene must be concordant in ALL other modalities
    all_offset = offsets[-1] + width
    all_fracs  = []
    for gd in comp_dirs:
        sub = df[df["comp_dir"] == gd]
        n   = len(sub)
        if n == 0 or not concordance_flags:
            all_fracs.append(0)
            continue
        # Intersect flags across all modalities for this direction
        combined = pd.Series([True] * len(sub), index=sub.index)
        for flags_per_dir in concordance_flags.values():
            if gd in flags_per_dir:
                combined = combined & flags_per_dir[gd]
        all_fracs.append(combined.sum() / n * 100)

    ax.bar(x + all_offset, all_fracs, width=width * 0.9,
           color="#2C3E50", label="All concordant",
           edgecolor="white", linewidth=0.4)

    # Widen x limits to accommodate the extra bar
    n_others  = len(other_mods)
    width     = 0.7 / (n_others + 1)
    offsets   = np.linspace(-(n_others) * width / 2,
                             (n_others) * width / 2,
                             n_others + 1)
    
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
        "Directional Concordance with Hi-C Compartment Switch\n"
        "Parental \u2192 Bottom10",
        fontweight="bold", pad=10
    )
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    plt.tight_layout()
    save_fig(fig, outdir, "Fig2c_concordance_bar")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print("[INFO] Loading merged table...")
    df = pd.read_csv(args.merged_file, sep="\t", low_memory=False)
    print(f"  {len(df):,} rows x {len(df.columns)} columns")

    # Build modality definitions and add boolean changed columns to df
    modalities = build_modalities(df, args)

    # Coupling label for each gene
    df["coupling_state"] = df.apply(
        lambda row: make_coupling_label(row, modalities), axis=1
    )
    print("\n[INFO] Coupling state distribution:")
    print(df["coupling_state"].value_counts().to_string())

    # Save
    out_path = os.path.join(args.outdir, "Q2_coupling_classified.txt")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"\n[INFO] Saved -> {out_path}")

    print("\n[INFO] Generating Fig 2a: UpSet plot...")
    plot_upset(df, modalities, args.outdir)

    print("[INFO] Generating Fig 2b: Spearman correlations...")
    plot_spearman(df, modalities, args.outdir)

    print("[INFO] Generating Fig 2c: Concordance bar chart...")
    plot_concordance_bar(df, modalities, args.outdir)

    print(f"\n[INFO] Q2 complete -> {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()