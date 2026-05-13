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
    p.add_argument("--outdir",      default="results/exploratory")

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
    UpSet plot using the upsetplot library
    """
    try:
        from upsetplot import UpSet, from_memberships
    except ImportError:
        print("[WARN] `upsetplot` not installed. "
              "Install with `pip install upsetplot` to enable Fig 2a. Skipping.")
        return

    # ── 1. Build per-gene membership lists ───────────────────────────────────
    # For each gene, record which modalities it is "changed" in.
    mod_names    = list(modalities.keys())
    changed_cols = [modalities[m]["changed_col"] for m in mod_names]

    # NaN-safe: missing values are treated as "not changed".
    changed_bool = df[changed_cols].fillna(False).astype(bool)

    # List-of-lists: for each row, which modalities are True.
    memberships = [
        [name for name, flag in zip(mod_names, row) if flag]
        for row in changed_bool.itertuples(index=False, name=None)
    ]

    # ── 2. Convert to the Series format UpSet expects ────────────────────────
    # `from_memberships(memberships, data=ones)` produces ONE Series row per
    # gene. That means the resulting MultiIndex has duplicate combinations,
    # which trips UpSet's default subset_size="auto" with:
    #   ValueError: subset_size="auto" cannot be used for a Series with
    #               non-unique groups.
    # The fix is twofold:
    #   (a) pass subset_size="sum" to UpSet so it groups duplicate index rows
    #       by summing their values (== counting genes per intersection), OR
    #   (b) pre-aggregate here with groupby().sum().
    # We use option (b) — it's cleaner, makes the data smaller, and avoids
    # any quirks with how different upsetplot versions parse "sum".
    upset_data = from_memberships(
        memberships,
        data=np.ones(len(memberships), dtype=int),
    )
    # Collapse the per-gene rows into one row per unique intersection.
    # The MultiIndex of booleans is preserved; values become counts.
    upset_data = upset_data.groupby(level=list(range(len(mod_names)))).sum()

    # ── 3. Configure the UpSet object ────────────────────────────────────────
    # show_counts='%d'              -> draws integer labels just outside bars,
    #                                  both intersection bars (top) and totals
    #                                  bars (side). Library handles offset so
    #                                  they don't intersect the bars.
    # facecolor='#222222'           -> dark grey bars (matches reference).
    # element_size=42               -> per-column width of the dot matrix.
    # intersection_plot_elements=8  -> taller top bar panel so big-count
    #                                  labels don't run into the title.
    # totals_plot_elements=4        -> wider side totals panel so side-count
    #                                  labels have room to the right of bars.
    upset = UpSet(
        upset_data,
        sort_by='cardinality',
        sort_categories_by='cardinality',
        show_counts='%d',
        facecolor='#222222',
        other_dots_color=0.18,        # grey level 0-1 for "absent" dots
        element_size=42,
        intersection_plot_elements=8,
        totals_plot_elements=4,
    )

    # ── 4. Highlight the "all modalities co-changed" subset ──────────────────
    n_mods = len(mod_names)
    upset.style_subsets(
        present=mod_names,        # the bar where every modality is True
        facecolor='#C0392B',      # red highlight
        label=f"All {n_mods} co-changed",
    )

    # ── 5. Draw onto a sized figure ──────────────────────────────────────────
    # Width scales gently with modality count so 2 vs 5 modalities both look
    # balanced. Height grows with category rows.
    try:
        fig = plt.figure(figsize=(max(11, 1.2 * n_mods + 8), 6.5 + 0.4 * n_mods))
        axes = upset.plot(fig=fig)   # returns dict of internal axes

        # ── 6. Polish the axes ───────────────────────────────────────────────────
        # 6a. Top bar panel y-label.
        axes['intersections'].set_ylabel("Intersection size", fontsize=11)

        # 6b. Side totals panel: widen the data limit so the numeric labels
        #     drawn by show_counts don't get clipped at the edge. The totals
        #     axis is normally inverted (large values extend leftward), so we
        #     extend whichever side represents the larger magnitude.
        totals_ax = axes['totals']
        totals_ax.set_xlabel("Set size", fontsize=10)
        xmin, xmax = totals_ax.get_xlim()
        pad = 0.20 * abs(xmin - xmax)
        if xmin > xmax:           # inverted (typical)
            totals_ax.set_xlim(xmin + pad, xmax)
        else:                     # non-inverted (defensive)
            totals_ax.set_xlim(xmin, xmax + pad)

        # 6c. Category tick labels on the matrix axis: push them away from the
        #     dots so "Hi-C Compartment" / "Transcription" have breathing room.
        matrix_ax = axes['matrix']
        for tick in matrix_ax.get_yticklabels():
            tick.set_horizontalalignment('right')
            tick.set_fontsize(10.5)
        matrix_ax.tick_params(axis='y', pad=8)

        # 6d. Format the legend created by style_subsets.
        leg = axes['intersections'].get_legend()
        if leg is not None:
            leg.set_frame_on(True)
            leg.set_bbox_to_anchor((1.0, 1.0))
            leg.set_loc('upper right')

        # ── 7. Title ─────────────────────────────────────────────────────────────
        mod_str = " + ".join(mod_names)
        fig.suptitle(
            f"Co-occurrence of Multi-omics Changes\n"
            f"Parental \u2192 Bottom10  ({mod_str})",
            fontweight="bold", fontsize=12, y=1.02,
        )

        # ── 8. Save (PDF + PNG, 600 dpi, via the existing helper) ────────────────
        save_fig(fig, outdir, "Fig2a_upset_cooccurrence")
    except (ValueError, TypeError) as e:
        print(f"[WARN] UpSet plot failed due to upsetplot/matplotlib compatibility issue: {e}")
        print("[WARN] Skipping Fig 2a. This is a known issue with certain upsetplot versions.")
        plt.close('all')    

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

    # Auto grid: square-ish layout. 1->1x1, 2->1x2, 3->1x3, 4->2x2,
    # 5->2x3, 6->2x3, 7-9->3x3, etc. Extra cells are hidden.
    n_pairs = len(pairs)
    ncols   = int(np.ceil(np.sqrt(n_pairs)))
    nrows   = int(np.ceil(n_pairs / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(4.5 * ncols, 4.3 * nrows),
                              squeeze=False)
    axes = axes.flatten()
    # Hide unused cells in non-balanced grids (e.g. 5 panels in a 2x3 grid)
    for ax in axes[n_pairs:]:
        ax.set_visible(False)

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
    fig, ax = plt.subplots(figsize=(max(7, len(other_mods) * 2.5), 5))

    # Add an extra series: "All concordant" = genes concordant in EVERY
    # non-Hi-C modality simultaneously. Renamed from "Both concordant" so
    # it scales to >2 modalities.
    n_series = len(other_mods) + 1                # +1 for "All concordant"
    width    = 0.75 / n_series
    offsets  = np.linspace(-(n_series - 1) * width / 2,
                            (n_series - 1) * width / 2,
                            n_series)

    # Helper: per-modality concordance boolean (sign-aware, repressive-aware)
    def concordance_mask(sub, info, gd):
        delta_col = info["delta_col"]
        if not delta_col or delta_col not in df.columns:
            return pd.Series([False] * len(sub), index=sub.index)
        vals       = pd.to_numeric(sub[delta_col], errors="coerce")
        repressive = any(k in delta_col.lower() for k in ["me3", "h3k9", "unc"])
        if gd == "B\u2192A":
            return (vals < 0) if repressive else (vals > 0)
        if gd == "A\u2192B":
            return (vals > 0) if repressive else (vals < 0)
        return pd.Series([False] * len(sub), index=sub.index)

    # Individual modality bars
    for (mod_name, info), offset in zip(other_mods.items(), offsets[:-1]):
        color = info["color"]
        fracs = []
        for gd in comp_dirs:
            sub = df[df["comp_dir"] == gd]
            n   = len(sub)
            if n == 0:
                fracs.append(0); continue
            mask = concordance_mask(sub, info, gd).fillna(False)
            fracs.append(mask.sum() / n * 100)
        ax.bar(x + offset, fracs, width=width * 0.9,
               color=color, label=f"Concordant {mod_name}",
               edgecolor="white", linewidth=0.4)

    # "All concordant" bar: gene must be concordant in EVERY non-Hi-C modality
    fracs_all = []
    for gd in comp_dirs:
        sub = df[df["comp_dir"] == gd]
        n   = len(sub)
        if n == 0:
            fracs_all.append(0); continue
        combined = pd.Series([True] * len(sub), index=sub.index)
        for info in other_mods.values():
            combined &= concordance_mask(sub, info, gd).fillna(False)
        fracs_all.append(combined.sum() / n * 100)
    ax.bar(x + offsets[-1], fracs_all, width=width * 0.9,
           color="#0A0B2C", label=f"All {len(other_mods)} concordant",
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
    
    # Create trimodal_coupling subdirectory for all figures
    fig_outdir = os.path.join(args.outdir, "trimodal_coupling")
    os.makedirs(fig_outdir, exist_ok=True)

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
    out_path = os.path.join(fig_outdir, "Q2_coupling_classified.txt")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"\n[INFO] Saved -> {out_path}")

    print("\n[INFO] Generating Fig 2a: UpSet plot...")
    plot_upset(df, modalities, fig_outdir)

    print("[INFO] Generating Fig 2b: Spearman correlations...")
    plot_spearman(df, modalities, fig_outdir)

    print("[INFO] Generating Fig 2c: Concordance bar chart...")
    plot_concordance_bar(df, modalities, fig_outdir)

    print(f"\n[INFO] Q2 complete -> {os.path.abspath(fig_outdir)}")


if __name__ == "__main__":
    main()