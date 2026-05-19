"""
Flexible trimodal (or multi-modal) coupling analysis

Hi-C compartment is always included by default.
You must specify at least one of --rna or --marks.

Usage examples:
    # Default behavior — no exports
    python scripts/analysis/trimodal_coupling_flex.py \
        data/processed/merged_hic_rnaseq_cutrun_parent_bottom10.txt \
        --rna --marks h3k27me3 h3k27ac

    # Concordant exports, both directions (bare flag — most common)
    python scripts/analysis/trimodal_coupling_flex.py \
        data/processed/merged_hic_rnaseq_cutrun_parent_bottom10.txt \
        --rna --marks h3k27me3 h3k27ac h3k4me1 h3k9me3 \
        --export_concordant

    # All 6 files — both flags bare
    python scripts/analysis/trimodal_coupling_flex.py \
        data/processed/merged_hic_rnaseq_cutrun_parent_bottom10.txt \
        --rna --marks h3k27me3 h3k27ac h3k4me1 h3k9me3 \
        --export_concordant --export_nonconcordant

    # One direction only
    python scripts/analysis/trimodal_coupling_flex.py \
        data/processed/merged_hic_rnaseq_cutrun_parent_bottom10.txt \
        --rna --marks h3k27me3 h3k27ac h3k4me1 h3k9me3 \
        --export_concordant BtoA

    # Mix-and-match: concordant both directions, nonconcordant only A→B
    python scripts/analysis/trimodal_coupling_flex.py \
        data/processed/merged_hic_rnaseq_cutrun_parent_bottom10.txt \
        --rna --marks h3k27me3 h3k27ac h3k4me1 h3k9me3 \
        --export_concordant \
        --export_nonconcordant AtoB

    # Future bin-mode (will fail cleanly today)
    python scripts/analysis/trimodal_coupling_flex.py merged_bins.txt \
        --rna --marks h3k27me3 --mode bins

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
    
    # ── Output / classification mode ──────────────────────────────────────────
    p.add_argument("--mode",
                   choices=["genes", "bins", "coordinates"],
                   default="genes",
                   help="Row classification mode (default: genes). "
                        "'bins' and 'coordinates' require an input table keyed "
                        "by genomic intervals rather than gene names.")

    p.add_argument("--export_concordant",
                   nargs="*",
                   choices=["BtoA", "AtoB", "both"],
                   metavar="DIR",
                   help="Export concordance gene tables. Optional values: "
                        "BtoA, AtoB, both. If the flag is given with no value, "
                        "defaults to both directions. For each chosen direction, "
                        "writes 2 files: all_concordant_<DIR>.txt (concordant in "
                        "EVERY non-Hi-C modality) and concordant_<DIR>.txt "
                        "(concordant in AT LEAST ONE non-Hi-C modality). "
                        "Stable genes excluded.")

    p.add_argument("--export_nonconcordant",
                   nargs="*",
                   choices=["BtoA", "AtoB", "both"],
                   metavar="DIR",
                   help="Export non-concordance gene tables. Optional values: "
                        "BtoA, AtoB, both. If the flag is given with no value, "
                        "defaults to both directions. For each chosen direction, "
                        "writes nonconcordant_<DIR>.txt (compartment switched "
                        "but NO concordant change in any non-Hi-C modality).")

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

    # Normalize export direction lists:
    #   None         -> flag absent       -> []
    #   []           -> flag bare         -> ["AtoB", "BtoA"]   (default: both)
    #   ["both"]     -> expand            -> ["AtoB", "BtoA"]
    #   ["BtoA"] etc -> as given (sorted, deduped)
    def _expand_dirs(lst):
        if lst is None:
            return []
        if len(lst) == 0:
            return ["AtoB", "BtoA"]
        out = set()
        for d in lst:
            if d == "both":
                out.update(["BtoA", "AtoB"])
            else:
                out.add(d)
        return sorted(out)

    args.export_concordant    = _expand_dirs(args.export_concordant)
    args.export_nonconcordant = _expand_dirs(args.export_nonconcordant)

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
                return "B→A"
            if "toB" in v or v in ["AtoB", "TowardB"]:
                return "A→B"
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
            "color":       MODALITY_COLORS[color_idx % len(MODALITY_COLORS)],
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
# PLOT A: UpSet plot
# ══════════════════════════════════════════════════════════════════════════════
def plot_upset(df, modalities, outdir):
    """
    UpSet plot using the upsetplot library
    """
    from upsetplot import UpSet, from_memberships

    # ── 1. Build per-gene membership lists ───────────────────────────────────
    mod_names    = list(modalities.keys())
    changed_cols = [modalities[m]["changed_col"] for m in mod_names]

    # NaN-safe: missing values are treated as "not changed".
    changed_bool = df[changed_cols].fillna(False).astype(bool)

    memberships = [
        [name for name, flag in zip(mod_names, row) if flag]
        for row in changed_bool.itertuples(index=False, name=None)
    ]

    # ── 2. Aggregate into one count per unique intersection ───────────────────
    upset_data = from_memberships(
        memberships,
        data=np.ones(len(memberships), dtype=int),
    )
    upset_data = upset_data.groupby(level=list(range(len(mod_names)))).sum()

    # ── 3. Configure the UpSet object ────────────────────────────────────────
    upset = UpSet(
        upset_data,
        sort_by='cardinality',
        sort_categories_by='cardinality',
        show_counts='%d',
        facecolor='#222222',
        other_dots_color="#837F7F",
        element_size=42,
        intersection_plot_elements=8,
        totals_plot_elements=4,
    )

    # ── 4. Highlight the "all modalities co-changed" subset ──────────────────
    n_mods = len(mod_names)
    upset.style_subsets(
        present=mod_names,
        facecolor='#C0392B',
        label=f"All {n_mods} co-changed",
    )

    # ── 5. Draw ──────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(max(11, 1.2 * n_mods + 8), 6.5 + 0.4 * n_mods))
    axes = upset.plot(fig=fig)

    # ── 6. Polish the axes ───────────────────────────────────────────────────
    # 6a. Top bar panel: label + push first bar away from y-axis so count
    #     labels don't overlap the spine.
    inter_ax = axes['intersections']
    inter_ax.set_ylabel("Intersection size", fontsize=11)
    ix_min, ix_max = inter_ax.get_xlim()
    inter_ax.set_xlim(ix_min - 0.4, ix_max)

    # 6b. Side totals panel: widen so count labels don't get clipped.
    totals_ax = axes['totals']
    totals_ax.set_xlabel("Set size", fontsize=10)
    xmin, xmax = totals_ax.get_xlim()
    pad = 0.20 * abs(xmin - xmax)
    if xmin > xmax:
        totals_ax.set_xlim(xmin + pad, xmax)
    else:
        totals_ax.set_xlim(xmin, xmax + pad)

    # 6c. Category tick labels: breathing room from the dot matrix.
    matrix_ax = axes['matrix']
    for tick in matrix_ax.get_yticklabels():
        tick.set_horizontalalignment('right')
        tick.set_fontsize(10.5)
    matrix_ax.tick_params(axis='y', pad=8)

    # 6d. Legend from style_subsets.
    leg = inter_ax.get_legend()
    if leg is not None:
        leg.set_frame_on(True)
        leg.set_bbox_to_anchor((1.0, 1.0))
        leg.set_loc('upper right')

    # ── 7. Title ─────────────────────────────────────────────────────────────
    mod_str = " + ".join(mod_names)
    fig.suptitle(
        f"Co-occurrence of Multi-omics Changes\n"
        f"Parental → Bottom10  ({mod_str})",
        fontweight="bold", fontsize=12, y=1.02,
    )

    # ── 8. Save ──────────────────────────────────────────────────────────────
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

        p_str = f"{p:.2e}" if p >= 1e-4 else "$< 10^{-4}$"
        ax.text(0.05, 0.95,
                f"ρ = {rho:.3f}\np = {p_str}",
                transform=ax.transAxes, fontsize=9, va="top",
                bbox=dict(boxstyle="round,pad=0.25",
                          fc="white", ec="none", alpha=0.85))

        ax.axhline(0, color="#888888", lw=0.6, ls=":")
        ax.axvline(0, color="#888888", lw=0.6, ls=":")

        # Axis labels: use delta symbol + modality name
        def delta_label(name, col):
            if "log" in col.lower() or "fc" in col.lower():
                return f"$\log_2 FC$ ({name})"
            return f"Δ{name} (B10 − Parental)"

        ax.set_xlabel(delta_label(xname, xcol), labelpad=5)
        ax.set_ylabel(delta_label(yname, ycol), labelpad=5)

    fig.suptitle(
        "Pairwise Spearman Correlations: " +
        " ↔ ".join([m[0] for m in mod_list]),
        fontweight="bold", fontsize=12, y=1.02
    )
    plt.tight_layout()
    save_fig(fig, outdir, "Fig2b_spearman_correlations")


# ══════════════════════════════════════════════════════════════════════════════
# PLOT C: Directional concordance bar
# ══════════════════════════════════════════════════════════════════════════════
def plot_concordance_bar(df, modalities, outdir):
    """
    For genes with a directional Hi-C compartment call (B->A or A->B),
    compute the fraction that show concordant changes in each other modality.

    Concordance definition:
        B->A: other modality shows activation/loss direction
              (positive log2FC, positive Δactivating mark, negative Δrepressive mark)
        A->B: opposite

    Stable-compartment genes are excluded — concordance is undefined for them.
    Their counts are still reported to stdout for the record.
    """
    if "Hi-C Compartment" not in modalities:
        print("[INFO] Hi-C not included — skipping concordance bar.")
        return

    df["comp_dir"] = df.apply(get_comp_dir, axis=1)

    # Report Stable count for the record (PI requested it not be plotted)
    n_stable = (df["comp_dir"] == "Stable").sum()
    print(f"[INFO] Stable-compartment rows excluded from concordance bar: "
          f"{n_stable:,}")

    # Other modalities (excluding Hi-C)
    other_mods = {
        name: info for name, info in modalities.items()
        if name != "Hi-C Compartment"
    }

    if not other_mods:
        print("[INFO] No other modalities to compare — skipping concordance bar.")
        return

    comp_dirs = ["B→A", "A→B"]                       # Stable removed
    x         = np.arange(len(comp_dirs))
    fig, ax = plt.subplots(figsize=(max(6, len(other_mods) * 2.2), 5.2))

    # +1 series for "All concordant"
    n_series = len(other_mods) + 1
    width    = 0.78 / n_series
    offsets  = np.linspace(-(n_series - 1) * width / 2,
                            (n_series - 1) * width / 2,
                            n_series)

    def concordance_mask(sub, info, gd):
        delta_col = info["delta_col"]
        if not delta_col or delta_col not in df.columns:
            return pd.Series([False] * len(sub), index=sub.index)
        vals       = pd.to_numeric(sub[delta_col], errors="coerce")
        repressive = any(k in delta_col.lower() for k in ["me3", "h3k9", "unc"])
        if gd == "B→A":
            return (vals < 0) if repressive else (vals > 0)
        if gd == "A→B":
            return (vals > 0) if repressive else (vals < 0)
        return pd.Series([False] * len(sub), index=sub.index)

    # Individual modality bars
    for (mod_name, info), offset in zip(other_mods.items(), offsets[:-1]):
        fracs = []
        for gd in comp_dirs:
            sub = df[df["comp_dir"] == gd]
            n   = len(sub)
            if n == 0:
                fracs.append(0); continue
            mask = concordance_mask(sub, info, gd).fillna(False)
            fracs.append(mask.sum() / n * 100)
        bars = ax.bar(x + offset, fracs, width=width * 0.92,
                      color=info["color"],
                      label=f"Concordant {mod_name}",
                      edgecolor="white", linewidth=0.4)
        ax.bar_label(bars, fmt="%.1f%%", padding=2, fontsize=8.5)

    # "All concordant" bar
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
    bars = ax.bar(x + offsets[-1], fracs_all, width=width * 0.92,
                  color="#0A0B2C",
                  label=f"All {len(other_mods)} concordant",
                  edgecolor="white", linewidth=0.4)
    ax.bar_label(bars, fmt="%.1f%%", padding=2, fontsize=8.5,
                 fontweight="bold")

    # Sample size annotations under each x-category
    ymin_text = -3
    for xi, gd in zip(x, comp_dirs):
        n = (df["comp_dir"] == gd).sum()
        ax.text(xi, ymin_text, f"n={n:,}", ha="center", va="top",
                fontsize=9, color="#444444")

    ax.set_xticks(x)
    ax.set_xticklabels(comp_dirs, fontsize=11)
    ax.axhline(0, color="black", lw=0.8)

    # Pad top so bar labels don't get clipped
    top_max = max([*fracs_all,
                   *(b.get_height() for c in ax.containers for b in c)])
    ax.set_ylim(bottom=-6, top=top_max * 1.15 + 2)

    ax.set_ylabel("% of genes with concordant change", labelpad=6)
    ax.set_xlabel("Compartment change direction", labelpad=12)
    ax.set_title(
        "Directional Concordance with Hi-C Compartment Switch\n"
        "Parental → Bottom10",
        fontweight="bold", pad=10
    )
    ax.legend(frameon=False, fontsize=9,
              loc="upper left", bbox_to_anchor=(1.01, 1.0),
              borderaxespad=0)
    plt.tight_layout(rect=[0, 0, 0.88, 1])
    save_fig(fig, outdir, "Fig2c_concordance_bar")

# ══════════════════════════════════════════════════════════════════════════════
# EXPORT: concordance / non-concordance gene tables
# ══════════════════════════════════════════════════════════════════════════════
def export_concordance_tables(df, modalities, id_cols, outdir,
                              concordant_dirs=None,
                              nonconcordant_dirs=None):
    """
    Write per-direction gene tables for downstream GO / GSEA / manual review.

    Requires Hi-C in modalities; uses df["comp_dir"] which must already
    be populated by plot_concordance_bar (or computed here if missing).

    Schema for every output file:
        id_cols + standard stat/delta columns + comp_dir + per-modality
        concordance booleans.

    Output files (only those whose flag is set):
        all_concordant_BtoA.txt / all_concordant_AtoB.txt
            -> concordant in EVERY non-Hi-C modality
        concordant_BtoA.txt     / concordant_AtoB.txt
            -> concordant in AT LEAST ONE non-Hi-C modality
        nonconcordant_BtoA.txt  / nonconcordant_AtoB.txt
            -> compartment switched, but NO modality is concordant
    """
    if not (concordant_dirs or nonconcordant_dirs):
        return
    if "Hi-C Compartment" not in modalities:
        print("[WARN] Cannot export concordance tables without Hi-C. Skipping.")
        return

    other_mods = {n: i for n, i in modalities.items() if n != "Hi-C Compartment"}
    if not other_mods:
        print("[WARN] No non-Hi-C modalities — nothing to export.")
        return

    # Ensure comp_dir is present
    if "comp_dir" not in df.columns:
        df["comp_dir"] = df.apply(get_comp_dir, axis=1)

    # Shared concordance mask helper (sign-aware, repressive-aware)
    def concordance_mask_full(info, comp_dir_series):
        delta_col = info["delta_col"]
        if not delta_col or delta_col not in df.columns:
            return pd.Series(False, index=df.index)
        vals       = pd.to_numeric(df[delta_col], errors="coerce")
        repressive = any(k in delta_col.lower() for k in ["me3", "h3k9", "unc"])
        out = pd.Series(False, index=df.index)
        if repressive:
            out.loc[comp_dir_series == "B→A"] = (vals < 0).loc[comp_dir_series == "B→A"]
            out.loc[comp_dir_series == "A→B"] = (vals > 0).loc[comp_dir_series == "A→B"]
        else:
            out.loc[comp_dir_series == "B→A"] = (vals > 0).loc[comp_dir_series == "B→A"]
            out.loc[comp_dir_series == "A→B"] = (vals < 0).loc[comp_dir_series == "A→B"]
        return out.fillna(False)

    # Per-modality concordance booleans (added as columns to the export)
    per_mod_conc_cols = {}
    for name, info in other_mods.items():
        col = f"_conc_{name.replace(' ', '_').replace('-', '_')}"
        df[col] = concordance_mask_full(info, df["comp_dir"])
        per_mod_conc_cols[name] = col

    # Aggregate masks
    any_conc = pd.Series(False, index=df.index)
    all_conc = pd.Series(True,  index=df.index)
    for col in per_mod_conc_cols.values():
        any_conc |= df[col]
        all_conc &= df[col]
    df["_concordant_any"] = any_conc
    df["_concordant_all"] = all_conc

    # Build standard schema: id_cols + common stat/delta columns + flags
    standard_cols = [
        "padj", "pvalue", "log2FC", "baseMean", "strand",
        "b10_minus_parental_h3k27ac",
        "b10_minus_parental_h3k27me3",
        "b10_minus_parental_h3k4me1",
        "b10_minus_parental_h3k9me3",
        "b10_minus_parental_unc_h3k27me3",
    ]
    standard_cols = [c for c in standard_cols if c in df.columns]
    flag_cols     = ["comp_dir", "_concordant_any", "_concordant_all"] + \
                    list(per_mod_conc_cols.values())
    export_cols   = id_cols + standard_cols + flag_cols

    # Output dir
    export_dir = os.path.join(outdir, "concordance_genes")
    os.makedirs(export_dir, exist_ok=True)

    def write(subset_df, fname):
        path = os.path.join(export_dir, fname)
        subset_df[export_cols].to_csv(path, sep="\t", index=False)
        print(f"  Exported {len(subset_df):>6,} rows -> {path}")

    print(f"\n[INFO] Writing concordance exports -> {export_dir}")

    # Map from CLI direction strings to the symbol used in comp_dir
    dir_map = {"BtoA": "B→A", "AtoB": "A→B"}

    if concordant_dirs:
        for fdir in concordant_dirs:
            label_dir = dir_map[fdir]
            sub_all = df[(df["comp_dir"] == label_dir) & df["_concordant_all"]]
            sub_any = df[(df["comp_dir"] == label_dir) & df["_concordant_any"]]
            write(sub_all, f"all_concordant_{fdir}.txt")
            write(sub_any, f"concordant_{fdir}.txt")

    if nonconcordant_dirs:
        for fdir in nonconcordant_dirs:
            label_dir = dir_map[fdir]
            sub_none = df[(df["comp_dir"] == label_dir) & (~df["_concordant_any"])]
            write(sub_none, f"nonconcordant_{fdir}.txt")

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

    # ── Resolve row-identifier columns based on --mode ───────────────────────
    if args.mode == "genes":
        id_cols = [c for c in ["gene_name", "gene_id"] if c in df.columns]
        if not id_cols:
            sys.exit("[ERROR] --mode genes requires a 'gene_name' or "
                     "'gene_id' column in the merged table.")
    elif args.mode == "bins":
        id_cols = [c for c in ["bin_id", "chrom", "start", "end"]
                   if c in df.columns]
        if "bin_id" not in df.columns:
            sys.exit("[ERROR] --mode bins requires a 'bin_id' column. "
                     "Re-merge inputs at bin resolution before running.")
    else:  # coordinates
        id_cols = [c for c in ["chrom", "start", "end", "strand"]
                   if c in df.columns]
        if not {"chrom", "start", "end"}.issubset(df.columns):
            sys.exit("[ERROR] --mode coordinates requires 'chrom', 'start', "
                     "and 'end' columns in the merged table.")

    print(f"[INFO] Classification mode: {args.mode} "
          f"(ID columns: {id_cols})")

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

    if args.export_concordant or args.export_nonconcordant:
        print("\n[INFO] Exporting concordance gene tables...")
        export_concordance_tables(
            df, modalities, id_cols, fig_outdir,
            concordant_dirs=args.export_concordant,
            nonconcordant_dirs=args.export_nonconcordant,
        )

    print(f"\n[INFO] Q2 complete -> {os.path.abspath(fig_outdir)}")


if __name__ == "__main__":
    main()
