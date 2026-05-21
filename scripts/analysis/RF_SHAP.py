"""
Epigenetic priming predicts compartment switching — Random Forest + SHAP

Adds a --comparison flag to refine the analysis along compartment-of-origin:
    A     : Stable_A vs A→B (silencing from A compartment)
    B     : Stable_B vs B→A (activation from B compartment)
    both  : Run A and B separately, side-by-side panels (DEFAULT)
    all   : Switcher-vs-stable framing dropped; all 4 categories included in RF and SHAP. The,
            violins use 4-category Kruskal-Wallis + Dunn's post-hoc with
            significance brackets and a pairwise heatmap supplement.

Usage:
    python scripts/RF_SHAP.py \
        data/processed/cut_run_signal_table.txt \
        data/processed/merged_hic_rnaseq_cutrun_parent_bottom10.txt \
        --outdir results/exploratory/RF_SHAP \
        --comparison both \
        --include_ctcf \
        --shap_dependence_partners H3K27me3:H3K27Ac H3K4me1:H3K9me3

Outputs (suffix depends on --comparison mode):
    Fig1a_switch_distribution_<mode>.pdf/png
    Fig1b_violin_marks_<mode>.pdf/png
    Fig1b_violin_pairwise_heatmap_<mode>.pdf/png  (all mode only)
    Fig1c_RF_ROC_PR_<mode>.pdf/png
    Fig1d_SHAP_<mode>.pdf/png
    Fig1d_dependence_<pair>_<mode>.pdf/png
    Fig1d_dependence_<pair>_<mode>_zoomed.pdf/png   (99th percentile clipped)
    Fig1e_PCA_marks_<mode>.pdf/png
    Statistics_<mode>.txt
    Statistics_dunn_<mark>_all.txt                (all mode only, per-mark)
    Switch_classified.txt                          (always emitted)
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
from itertools import combinations
from scipy.stats import mannwhitneyu, kruskal
from statsmodels.stats.multitest import multipletests
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (
    roc_curve, auc,
    average_precision_score, precision_recall_curve,
)
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import shap

# scikit-posthocs is needed for Dunn's test in 'all' mode.
# Install: pip install scikit-posthocs
try:
    import scikit_posthocs as sp
except ImportError:
    sp = None


# ── Aesthetics ────────────────────────────────────────────────────
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

SWITCH_COLORS = {
    "B_to_A":   "#C0392B",
    "A_to_B":   "#1A5276",
    "Stable_A": "#85C1E9",
    "Stable_B": "#D5D8DC",
}
SWITCH_LABELS = {
    "B_to_A":   "B\u2192A (Activated)",
    "A_to_B":   "A\u2192B (Silenced)",
    "Stable_A": "Stable A",
    "Stable_B": "Stable B",
}

# Comparison-mode definitions:
#   negative_class -> reference (Stable)
#   positive_class -> switcher target
#   label          -> short human label used in titles/filenames
COMPARISON_SPEC = {
    "A": {
        "neg":   "Stable_A",
        "pos":   "A_to_B",
        "label": "Stable_A vs A\u2192B",
        "short": "compA",
    },
    "B": {
        "neg":   "Stable_B",
        "pos":   "B_to_A",
        "label": "Stable_B vs B\u2192A",
        "short": "compB",
    },
}


# =============================================================================
# CLI
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Epigenetic priming predicts compartment switching"
    )
    p.add_argument("cutrun_file",  help="cut_run_signal_table.txt")
    p.add_argument("merged_file",  help="merged multi-omics table")
    p.add_argument("--outdir",     default="results/exploratory/RF_SHAP",
                   help="Output directory for figures and tables")
    p.add_argument("--include_ctcf", action="store_true",
                   help="Include parental_ctcf_mean as a 5th feature if present.")
    p.add_argument(
        "--shap_dependence_partners",
        nargs="+",
        default=None,
        metavar="MARK1:MARK2",
        help=(
            "Pairs of marks for SHAP dependence plots, formatted MARK1:MARK2. "
            "For each pair, two panels (or 4 in 'both' mode) are generated. "
            "Example: --shap_dependence_partners H3K27me3:H3K27Ac H3K4me1:H3K9me3"
        ),
    )
    p.add_argument(
        "--comparison",
        choices=["A", "B", "both", "all"],
        default="both",
        help=(
            "What to compare:\n"
            "  A    -> Stable_A vs A->B (silencing from A)\n"
            "  B    -> Stable_B vs B->A (activation from B)\n"
            "  both -> Run A and B separately, panels side-by-side (default)\n"
            "  all  -> Switcher-vs-stable binary RF; "
            "violin uses 4-category Kruskal-Wallis + Dunn's post-hoc."
        ),
    )
    return p.parse_args()


def save_fig(fig, outdir, stem):
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        path = os.path.join(outdir, f"{stem}.{ext}")
        fig.savefig(path, bbox_inches="tight", pad_inches=0.05, dpi=600)
        print(f"  Saved -> {path}")
    plt.close(fig)


# =============================================================================
# Helpers
# =============================================================================

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


def cliffs_delta_fast(x, y):
    """
    Memory-efficient Cliff's delta via sort + searchsorted.
    Convention: delta(x, y) = P(x > y) - P(x < y) in {-1, +1}.
    Magnitude thresholds: Romano et al. 2006.
    """
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return np.nan, "NA"
    ys = np.sort(y)
    lt_per_x = np.searchsorted(ys, x, side="left")
    le_per_x = np.searchsorted(ys, x, side="right")
    gt_per_x = ny - le_per_x
    total_gt = int(np.sum(gt_per_x))
    total_lt = int(np.sum(lt_per_x))
    d  = (total_gt - total_lt) / (nx * ny)
    ad = abs(d)
    if   ad < 0.147: mag = "negligible"
    elif ad < 0.330: mag = "small"
    elif ad < 0.474: mag = "medium"
    else:            mag = "large"
    return d, mag


def classify_switch(merged_df):
    cats = []
    for _, row in merged_df.iterrows():
        cat = None
        v = str(row.get("P->B10_Both_original", "")).strip()
        if v not in ["-", "", "nan", "None"]:
            if "toA" in v or "BtoA" in v:
                cat = "B_to_A"
            elif "toB" in v or "AtoB" in v:
                cat = "A_to_B"
        if cat is None:
            v = str(row.get("P->B10_Switch_Orig", "")).strip()
            if v not in ["-", "", "nan", "None"]:
                if v == "BtoA":   cat = "B_to_A"
                elif v == "AtoB": cat = "A_to_B"
        if cat is None:
            v = str(row.get("P->B10_Shift_Orig", "")).strip()
            if v not in ["-", "", "nan", "None"]:
                if "TowardA" in v:   cat = "B_to_A"
                elif "TowardB" in v: cat = "A_to_B"
        if cat is None:
            p_orig = row.get("Parental_orig", np.nan)
            b_orig = row.get("Bottom10_orig", np.nan)
            try:
                p_orig = float(p_orig); b_orig = float(b_orig)
                if   p_orig > 0 and b_orig > 0: cat = "Stable_A"
                elif p_orig < 0 and b_orig < 0: cat = "Stable_B"
            except (ValueError, TypeError):
                pass
        cats.append(cat)
    return cats


def load_data(cutrun_file, merged_file, include_ctcf=False):
    print("[INFO] Loading CUT&RUN signal table...")
    cutrun = pd.read_csv(cutrun_file, sep="\t", low_memory=False)
    print(f"  {len(cutrun):,} rows x {len(cutrun.columns)} columns")

    print("[INFO] Loading merged multi-omics table...")
    merged = pd.read_csv(merged_file, sep="\t", low_memory=False)
    print(f"  {len(merged):,} rows x {len(merged.columns)} columns")

    mark_spec = [
        ("H3K27me3", ["parental", "h3k27me3", "mean"]),
        ("H3K27Ac",  ["parental", "h3k27ac",  "mean"]),
        ("H3K4me1",  ["parental", "h3k4me1",  "mean"]),
        ("H3K9me3",  ["parental", "h3k9me3",  "mean"]),
    ]
    if include_ctcf:
        mark_spec.append(("CTCF", ["parental", "ctcf", "mean"]))

    mark_map = {}
    for mark, keywords in mark_spec:
        col = find_col(cutrun, keywords, ["minus", "unc", "ezh", "b10", "bottom"])
        if col:
            mark_map[mark] = col
            print(f"  Detected {mark}: {col}")
        else:
            if mark == "CTCF" and include_ctcf:
                print(f"  [WARN] --include_ctcf set but no parental CTCF column "
                      f"found; continuing without CTCF.")
            else:
                print(f"  [WARN] Could not detect {mark} column")

    if len(mark_map) < 2:
        sys.exit("[ERROR] Need at least 2 Parental mark columns.")

    feat_cols  = ["gene_name"] + list(mark_map.values())
    cutrun_sub = cutrun[feat_cols].dropna(subset=["gene_name"])
    cutrun_sub = cutrun_sub.drop_duplicates(subset="gene_name")
    for col in mark_map.values():
        cutrun_sub[col] = pd.to_numeric(cutrun_sub[col], errors="coerce")
    rename = {v: k for k, v in mark_map.items()}
    cutrun_sub = cutrun_sub.rename(columns=rename)

    print("[INFO] Classifying compartment switches...")
    switch_cols = [
        "gene_name", "Parental_orig", "Bottom10_orig",
        "P->B10_Switch_Orig", "P->B10_Shift_Orig",
        "P->B10_Both_original", "B10-P_Orig"
    ]
    avail = [c for c in switch_cols if c in merged.columns]
    merged_sw = merged[avail].drop_duplicates(subset="gene_name").copy()
    merged_sw["switch_cat"] = classify_switch(merged_sw)
    print(f"[INFO] Switch category counts:")
    print(merged_sw["switch_cat"].value_counts(dropna=False).to_string())

    df = cutrun_sub.merge(
        merged_sw[["gene_name", "switch_cat"] +
                  [c for c in avail if c != "gene_name"]],
        on="gene_name", how="inner"
    )
    df = df[df["switch_cat"].notna()].copy()
    df = df[df["switch_cat"].isin(["B_to_A", "A_to_B",
                                    "Stable_A", "Stable_B"])]

    mark_names = list(mark_map.keys())
    df[mark_names] = df[mark_names].apply(pd.to_numeric, errors="coerce")
    df = df.dropna(subset=mark_names)

    print(f"[INFO] Final dataset: {len(df):,} genes x {len(mark_names)} marks")
    return df, mark_names


def subset_for_comparison(df, comp_key):
    """
    Filter df to the two categories relevant for comparison A or B.

    comp_key in {'A', 'B'} -> returns the binary-filtered df + spec dict.

    The negative class (Stable) is encoded as 0 in 'is_switcher';
    the positive class (the switcher) is encoded as 1.
    """
    spec = COMPARISON_SPEC[comp_key]
    sub  = df[df["switch_cat"].isin([spec["neg"], spec["pos"]])].copy()
    sub["is_switcher"] = (sub["switch_cat"] == spec["pos"]).astype(int)
    return sub, spec


# =============================================================================
# Fig 1a — switch distribution
# =============================================================================

def plot_switch_distribution(df, outdir, comparison):
    """
    Fig 1a — bar chart of switch category counts.
    - both / all: all 4 categories
    - A:    only Stable_A and A_to_B
    - B:    only Stable_B and B_to_A
    """
    if comparison in ("both", "all"):
        cats = ["B_to_A", "A_to_B", "Stable_A", "Stable_B"]
        title_sfx = ""
    elif comparison == "A":
        cats = ["Stable_A", "A_to_B"]
        title_sfx = "  (A compartment only)"
    elif comparison == "B":
        cats = ["Stable_B", "B_to_A"]
        title_sfx = "  (B compartment only)"
    else:
        cats = ["B_to_A", "A_to_B", "Stable_A", "Stable_B"]
        title_sfx = ""

    cats   = [c for c in cats if c in df["switch_cat"].unique()]
    counts = [df["switch_cat"].value_counts().get(c, 0) for c in cats]
    colors = [SWITCH_COLORS[c] for c in cats]
    labels = [SWITCH_LABELS[c] for c in cats]

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    bars = ax.bar(labels, counts, color=colors,
                  edgecolor="white", linewidth=0.5, width=0.6)
    for bar, val in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(counts) * 0.012,
                f"{val:,}", ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("Number of genes", labelpad=6)
    ax.set_title("Hi-C Compartment Switch Categories\nParental \u2192 Bottom10"
                 + title_sfx,
                 fontweight="bold", pad=10)
    ax.tick_params(axis="x", rotation=20)
    ax.set_ylim(0, max(counts) * 1.15)
    plt.tight_layout()
    save_fig(fig, outdir, f"Fig1a_switch_distribution_{comparison}")


# =============================================================================
# Fig 1b — violins (mode-aware)
# =============================================================================

def _draw_violin_pair(ax, df, mark, neg_cat, pos_cat, title_extra=""):
    """
    Helper: draw a 2-category log1p violin (Stable vs Switcher) on `ax`,
    with colored medians + Cliff's delta + MWU annotation.
    Used in A/B/both modes.
    """
    vals_neg = np.log1p(df.loc[df["switch_cat"] == neg_cat, mark].dropna().values)
    vals_pos = np.log1p(df.loc[df["switch_cat"] == pos_cat, mark].dropna().values)

    if len(vals_neg) < 2 or len(vals_pos) < 2:
        ax.set_visible(False); return

    parts = ax.violinplot([vals_neg, vals_pos],
                          showmedians=False, showextrema=False, widths=0.7)
    for pc, cat in zip(parts["bodies"], [neg_cat, pos_cat]):
        pc.set_facecolor(SWITCH_COLORS[cat])
        pc.set_alpha(0.75)
        pc.set_edgecolor("none")
    for i, (vals, cat) in enumerate(
            zip([vals_neg, vals_pos], [neg_cat, pos_cat]), start=1):
        med = np.median(vals)
        ax.hlines(med, i - 0.22, i + 0.22, colors="black",
                  linewidth=2.2, zorder=3)
        ax.hlines(med, i - 0.18, i + 0.18, colors=SWITCH_COLORS[cat],
                  linewidth=1.4, zorder=4)

    _, p_mwu = mannwhitneyu(vals_pos, vals_neg, alternative="two-sided")
    d, mag   = cliffs_delta_fast(vals_pos, vals_neg)
    p_str = "p < 1e-300" if (p_mwu == 0 or p_mwu < 1e-300) else f"p = {p_mwu:.1e}"
    ax.text(
        0.97, 0.97,
        f"{SWITCH_LABELS[pos_cat]} vs {SWITCH_LABELS[neg_cat]}\n"
        f"Cliff's \u03b4 = {d:+.3f} ({mag})\n{p_str}",
        transform=ax.transAxes, ha="right", va="top",
        fontsize=7.5,
        bbox=dict(boxstyle="round,pad=0.25",
                  facecolor="white", edgecolor="lightgrey",
                  linewidth=0.5, alpha=0.85),
    )
    ax.set_xticks([1, 2])
    ax.set_xticklabels([SWITCH_LABELS[neg_cat], SWITCH_LABELS[pos_cat]],
                       fontsize=8, rotation=20, ha="right")
    ax.set_ylabel("log1p(CPM signal)", labelpad=4)
    ax.set_title(mark + title_extra, fontweight="bold", fontsize=10)


def plot_violin_marks_AB(df, mark_names, outdir, comparison):
    """
    Mode A/B: one panel per mark, all in a grid (3-column).
    Mode both: 4-column figure -> for each ROW: [mark_i compA | mark_i compB |
                                                  mark_{i+1} compA | mark_{i+1} compB]
              i.e. 2 marks per row, 4 panels per row.
    """
    if comparison in ("A", "B"):
        spec = COMPARISON_SPEC[comparison]
        n_marks = len(mark_names)
        n_cols  = 3
        n_rows  = (n_marks + n_cols - 1) // n_cols
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
        axes = np.atleast_2d(axes).flatten()
        for ax, mark in zip(axes, mark_names):
            _draw_violin_pair(ax, df, mark, spec["neg"], spec["pos"])
        for ax in axes[len(mark_names):]:
            ax.set_visible(False)
        fig.suptitle(
            f"Parental Epigenetic Marks: {spec['label']}",
            fontweight="bold", fontsize=12, y=1.00,
        )
        plt.tight_layout()
        save_fig(fig, outdir, f"Fig1b_violin_marks_{comparison}")
        return

    # comparison == "both": 4 columns, 2 marks per row
    n_marks = len(mark_names)
    marks_per_row = 2
    n_rows = (n_marks + marks_per_row - 1) // marks_per_row
    n_cols = 4  # 2 marks * 2 sub-panels (A, B)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.8 * n_cols, 4 * n_rows))
    axes = np.atleast_2d(axes)

    for row in range(n_rows):
        for slot in range(marks_per_row):
            mark_idx = row * marks_per_row + slot
            if mark_idx >= n_marks:
                # hide both panels for this slot
                axes[row, slot * 2].set_visible(False)
                axes[row, slot * 2 + 1].set_visible(False)
                continue
            mark = mark_names[mark_idx]
            ax_a = axes[row, slot * 2]
            ax_b = axes[row, slot * 2 + 1]
            _draw_violin_pair(ax_a, df, mark,
                              COMPARISON_SPEC["A"]["neg"],
                              COMPARISON_SPEC["A"]["pos"],
                              title_extra="\n(A compartment)")
            _draw_violin_pair(ax_b, df, mark,
                              COMPARISON_SPEC["B"]["neg"],
                              COMPARISON_SPEC["B"]["pos"],
                              title_extra="\n(B compartment)")

    fig.suptitle(
        "Parental Epigenetic Marks: Compartment-Resolved Stable vs Switcher",
        fontweight="bold", fontsize=12, y=1.00,
    )
    plt.tight_layout()
    save_fig(fig, outdir, "Fig1b_violin_marks_both")


def _bracket(ax, x1, x2, y, label, lw=1.0, fs=8):
    """Draw a significance bracket between x1 and x2 at height y."""
    dy = (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.012
    ax.plot([x1, x1, x2, x2], [y, y + dy, y + dy, y], lw=lw, c="black")
    ax.text((x1 + x2) / 2, y + dy * 1.2, label,
            ha="center", va="bottom", fontsize=fs)


def _p_to_stars(p):
    if p < 1e-4:  return "****"
    if p < 1e-3:  return "***"
    if p < 1e-2:  return "**"
    if p < 5e-2:  return "*"
    return "ns"


def plot_violin_marks_all(df, mark_names, outdir):
    """
    'all' mode: 4-category violins per mark with Kruskal-Wallis omnibus +
    Dunn's pairwise post-hoc (BH correction) + significance brackets.
    Pairwise Cliff's delta + p heatmap saved as a separate figure.
    """
    if sp is None:
        print("[ERROR] scikit-posthocs not installed.")
        print("        Install with: pip install scikit-posthocs")
        print("        Falling back to no post-hoc; brackets will be omitted.")

    cat_order = ["Stable_A", "Stable_B", "A_to_B", "B_to_A"]
    cat_order = [c for c in cat_order if c in df["switch_cat"].unique()]
    n_marks   = len(mark_names)
    n_cols    = 3
    n_rows    = (n_marks + n_cols - 1) // n_cols

    # Storage for Dunn's tables (per mark) and the heatmap (per mark)
    dunn_tables = {}

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.2 * n_cols, 4.2 * n_rows))
    axes = np.atleast_2d(axes).flatten()

    for ax, mark in zip(axes, mark_names):
        data_per_cat = []
        keep_cats    = []
        for c in cat_order:
            vals = df.loc[df["switch_cat"] == c, mark].dropna().values
            vals = np.log1p(vals)
            if len(vals) > 1:
                data_per_cat.append(vals); keep_cats.append(c)
        if not data_per_cat:
            ax.set_visible(False); continue

        parts = ax.violinplot(data_per_cat, showmedians=False,
                              showextrema=False, widths=0.7)
        for pc, cat in zip(parts["bodies"], keep_cats):
            pc.set_facecolor(SWITCH_COLORS[cat]); pc.set_alpha(0.75)
            pc.set_edgecolor("none")
        for i, (vals, cat) in enumerate(zip(data_per_cat, keep_cats), start=1):
            med = np.median(vals)
            ax.hlines(med, i - 0.22, i + 0.22, colors="black",
                      linewidth=2.2, zorder=3)
            ax.hlines(med, i - 0.18, i + 0.18, colors=SWITCH_COLORS[cat],
                      linewidth=1.4, zorder=4)

        # Kruskal-Wallis omnibus
        if len(data_per_cat) >= 2:
            kw_stat, kw_p = kruskal(*data_per_cat)
        else:
            kw_stat, kw_p = np.nan, np.nan
        kw_str = ("p < 1e-300" if (kw_p == 0 or kw_p < 1e-300)
                  else f"p = {kw_p:.1e}")
        ax.text(
            0.03, 0.97,
            f"Kruskal-Wallis\nH = {kw_stat:.1f}\n{kw_str}",
            transform=ax.transAxes, ha="left", va="top", fontsize=7.5,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      edgecolor="lightgrey", linewidth=0.5, alpha=0.85),
        )

        # Dunn's post-hoc (BH) + Cliff's delta per pair
        if sp is not None and len(data_per_cat) >= 2:
            # Long form for Dunn
            long_df = pd.concat([
                pd.DataFrame({"val": d, "cat": c})
                for d, c in zip(data_per_cat, keep_cats)
            ], ignore_index=True)
            dunn = sp.posthoc_dunn(
                long_df, val_col="val", group_col="cat", p_adjust="fdr_bh"
            )

            # Collect per-pair stats (p + Cliff's delta) for export
            rows = []
            for c1, c2 in combinations(keep_cats, 2):
                p_pair = float(dunn.loc[c1, c2])
                d_pair, mag_pair = cliffs_delta_fast(
                    df.loc[df["switch_cat"] == c2, mark].dropna().values,
                    df.loc[df["switch_cat"] == c1, mark].dropna().values,
                )
                rows.append({
                    "group1": c1, "group2": c2,
                    "dunn_p_BH": p_pair,
                    "cliffs_delta_g2_vs_g1": d_pair,
                    "effect_magnitude": mag_pair,
                    "stars": _p_to_stars(p_pair),
                })
            dunn_tables[mark] = pd.DataFrame(rows)

            # Significance brackets — draw only pairs with p < 0.05
            # to keep the figure readable.
            sig_pairs = [r for r in rows if r["dunn_p_BH"] < 0.05]
            # Order brackets bottom-up so closer pairs are lower
            cat_to_pos = {c: i + 1 for i, c in enumerate(keep_cats)}
            sig_pairs.sort(key=lambda r: abs(cat_to_pos[r["group1"]]
                                             - cat_to_pos[r["group2"]]))
            top = max(np.max(d) for d in data_per_cat)
            tick = (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.07
            base_y = top + tick * 0.6
            for k, r in enumerate(sig_pairs):
                x1 = cat_to_pos[r["group1"]]
                x2 = cat_to_pos[r["group2"]]
                _bracket(ax, x1, x2, base_y + k * tick, r["stars"])

            # Expand y-range so brackets fit
            ax.set_ylim(ax.get_ylim()[0],
                        base_y + (len(sig_pairs) + 1) * tick)

        ax.set_xticks(range(1, len(keep_cats) + 1))
        ax.set_xticklabels([SWITCH_LABELS[c] for c in keep_cats],
                           fontsize=8, rotation=25, ha="right")
        ax.set_title(mark, fontweight="bold", fontsize=10)
        ax.set_ylabel("log1p(CPM signal)", labelpad=4)

    for ax in axes[len(mark_names):]:
        ax.set_visible(False)

    fig.suptitle(
        "Parental Epigenetic Marks across Compartment-Switch Categories"
        "\nKruskal-Wallis omnibus + Dunn's BH-adjusted post-hoc",
        fontweight="bold", fontsize=12, y=1.00,
    )
    plt.tight_layout()
    save_fig(fig, outdir, "Fig1b_violin_marks_all")

    # Save per-mark Dunn tables
    for mark, tab in dunn_tables.items():
        out_tbl = os.path.join(outdir, f"Statistics_dunn_{mark}_all.txt")
        tab.to_csv(out_tbl, sep="\t", index=False)
        print(f"  Dunn table -> {out_tbl}")

    # Supplementary heatmap of pairwise p-values per mark
    if dunn_tables:
        plot_pairwise_heatmap(dunn_tables, mark_names, outdir)


def plot_pairwise_heatmap(dunn_tables, mark_names, outdir):
    """
    Supplementary heatmap: one mini-heatmap per mark showing pairwise
    Dunn p-values (BH-adjusted), -log10 transformed for visibility.
    """
    marks_with_data = [m for m in mark_names if m in dunn_tables]
    if not marks_with_data:
        return

    n_marks = len(marks_with_data)
    n_cols  = 3
    n_rows  = (n_marks + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3.6 * n_rows))
    axes = np.atleast_2d(axes).flatten()

    # Establish a stable category order present across marks
    cat_order = ["Stable_A", "Stable_B", "A_to_B", "B_to_A"]

    for ax, mark in zip(axes, marks_with_data):
        tab = dunn_tables[mark]
        # Build a square symmetric matrix
        cats = [c for c in cat_order
                if c in set(tab["group1"]).union(tab["group2"])]
        mat = pd.DataFrame(np.nan, index=cats, columns=cats)
        for _, r in tab.iterrows():
            v = -np.log10(max(r["dunn_p_BH"], 1e-300))
            mat.loc[r["group1"], r["group2"]] = v
            mat.loc[r["group2"], r["group1"]] = v
        np.fill_diagonal(mat.values, 0)

        im = ax.imshow(mat.values, cmap="viridis", aspect="auto",
                       vmin=0, vmax=max(50, np.nanmax(mat.values)))
        ax.set_xticks(range(len(cats)))
        ax.set_yticks(range(len(cats)))
        ax.set_xticklabels([SWITCH_LABELS[c] for c in cats],
                           rotation=30, ha="right", fontsize=8)
        ax.set_yticklabels([SWITCH_LABELS[c] for c in cats], fontsize=8)
        ax.set_title(f"{mark}\n-log10(Dunn p_BH)",
                     fontweight="bold", fontsize=10, pad=6)
        for i in range(len(cats)):
            for j in range(len(cats)):
                if not np.isnan(mat.iloc[i, j]) and i != j:
                    txt = f"{mat.iloc[i,j]:.0f}"
                    ax.text(j, i, txt, ha="center", va="center",
                            fontsize=7,
                            color="white" if mat.iloc[i,j] > 25 else "black")
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    for ax in axes[len(marks_with_data):]:
        ax.set_visible(False)

    fig.suptitle("Pairwise Dunn's Test (BH-adjusted) — per mark",
                 fontweight="bold", fontsize=12, y=1.00)
    plt.tight_layout()
    save_fig(fig, outdir, "Fig1b_violin_pairwise_heatmap_all")


# =============================================================================
# Fig 1c — ROC + PR (mode-aware)
# =============================================================================

def _run_rf_cv(X, y, seed=42):
    """
    Run 5-fold CV RF on (X, y) and return everything needed to plot ROC/PR
    and to compute SHAP later.
    Returns dict with keys: aucs, aps, mean_tpr, std_tpr, mean_prec, std_prec,
                            mfpr, mrec, prevalence, mean_auc, mean_ap, std_auc,
                            std_ap, rf_refit (RF fit on the full data).
    """
    rf = RandomForestClassifier(
        n_estimators=1000, max_depth=None,
        class_weight="balanced_subsample", random_state=seed, n_jobs=-1
    )
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    aucs, aps = [], []
    tprs, prcs = [], []
    mfpr = np.linspace(0, 1, 100)
    mrec = np.linspace(0, 1, 100)
    for tr_idx, te_idx in cv.split(X, y):
        rf.fit(X[tr_idx], y[tr_idx])
        y_prob = rf.predict_proba(X[te_idx])[:, 1]
        fpr, tpr, _ = roc_curve(y[te_idx], y_prob)
        aucs.append(auc(fpr, tpr))
        tprs.append(np.interp(mfpr, fpr, tpr))
        precision, recall, _ = precision_recall_curve(y[te_idx], y_prob)
        aps.append(average_precision_score(y[te_idx], y_prob))
        prcs.append(np.interp(mrec, recall[::-1], precision[::-1]))

    # Refit on full dataset for SHAP
    rf.fit(X, y)

    return {
        "aucs": np.array(aucs),
        "aps":  np.array(aps),
        "tprs": np.array(tprs),
        "prcs": np.array(prcs),
        "mfpr": mfpr, "mrec": mrec,
        "prevalence": float(y.mean()),
        "mean_auc": float(np.mean(aucs)), "std_auc": float(np.std(aucs)),
        "mean_ap":  float(np.mean(aps)),  "std_ap":  float(np.std(aps)),
        "rf_refit": rf,
    }


def _draw_roc_panel(ax, res, title_extra=""):
    """ROC sub-panel with random line + per-fold AUC inset."""
    for tpr_fold in res["tprs"]:
        ax.plot(res["mfpr"], tpr_fold, color="#85C1E9", alpha=0.15, lw=1)
    mean_tpr = res["tprs"].mean(axis=0)
    std_tpr  = res["tprs"].std(axis=0)
    ax.plot(res["mfpr"], mean_tpr, color="#C0392B", lw=2,
            label=f"Mean ROC (AUC = {res['mean_auc']:.3f} "
                  f"\u00b1 {res['std_auc']:.3f})")
    ax.fill_between(res["mfpr"], mean_tpr - std_tpr, mean_tpr + std_tpr,
                    color="#E74C3C", alpha=0.12)
    ax.plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.5,
            label="Random (AUC = 0.5)")
    ax.set_xlabel("False Positive Rate", labelpad=6)
    ax.set_ylabel("True Positive Rate", labelpad=6)
    ax.set_title("ROC Curve (5-fold CV)" + title_extra,
                 fontweight="bold", pad=8)
    ax.legend(frameon=False, fontsize=8.5, loc="lower right")

    ax_in = ax.inset_axes([0.10, 0.62, 0.25, 0.22])
    ax_in.scatter(res["aucs"], [1] * len(res["aucs"]),
                  color="#C0392B", s=35, zorder=3,
                  edgecolor="white", linewidth=0.6)
    ax_in.axvline(res["mean_auc"], color="black", lw=1.0,
                  linestyle="--", zorder=2)
    ax_in.set_title("Per-fold AUC", fontsize=8, pad=3)
    ax_in.tick_params(axis="x", labelsize=7)
    ax_in.set_yticks([])
    ax_in.set_xlim(res["aucs"].min() - 0.005, res["aucs"].max() + 0.005)
    ax_in.set_ylim(0.5, 1.5)
    for s in ("top", "right", "left"):
        ax_in.spines[s].set_visible(False)


def _draw_pr_panel(ax, res, title_extra=""):
    """PR sub-panel with prevalence baseline."""
    for prc_fold in res["prcs"]:
        ax.plot(res["mrec"], prc_fold, color="#85C1E9", alpha=0.15, lw=1)
    mean_prec = res["prcs"].mean(axis=0)
    std_prec  = res["prcs"].std(axis=0)
    ax.plot(res["mrec"], mean_prec, color="#C0392B", lw=2,
            label=f"Mean PR (PR-AUC = {res['mean_ap']:.3f} "
                  f"\u00b1 {res['std_ap']:.3f})")
    ax.fill_between(res["mrec"], mean_prec - std_prec, mean_prec + std_prec,
                    color="#E74C3C", alpha=0.12)
    ax.axhline(res["prevalence"], color="black", lw=0.8, linestyle="--",
               alpha=0.5,
               label=f"Random (PR-AUC = {res['prevalence']:.3f})")
    ax.set_xlabel("Recall", labelpad=6)
    ax.set_ylabel("Precision", labelpad=6)
    ax.set_title("Precision-Recall Curve (5-fold CV)" + title_extra,
                 fontweight="bold", pad=8)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")


def plot_roc_pr(results_by_mode, outdir, comparison):
    """
    Plot ROC + PR.
        single-mode (A / B / all): 1x2 figure.
        both:                      2x2 figure (top=A, bottom=B).
    """
    if comparison in ("A", "B", "all"):
        if comparison == "all":
            mode_label = "Switchers vs Stable"
        else:
            mode_label = COMPARISON_SPEC[comparison]["label"]
        res = results_by_mode[comparison]
        fig, (ax_roc, ax_pr) = plt.subplots(1, 2, figsize=(11, 5))
        _draw_roc_panel(ax_roc, res)
        _draw_pr_panel(ax_pr, res)
        fig.suptitle(f"Random Forest: Parental Marks Predict "
                     f"Compartment Switching\n({mode_label})",
                     fontweight="bold", fontsize=12, y=1.02)
        plt.tight_layout()
        save_fig(fig, outdir, f"Fig1c_RF_ROC_PR_{comparison}")
        return

    # both mode -> 2x2
    fig, axes = plt.subplots(2, 2, figsize=(11, 10))
    for r_idx, comp in enumerate(["A", "B"]):
        res = results_by_mode[comp]
        spec = COMPARISON_SPEC[comp]
        _draw_roc_panel(axes[r_idx, 0], res,
                        title_extra=f"\n({spec['label']})")
        _draw_pr_panel(axes[r_idx, 1], res,
                       title_extra=f"\n({spec['label']})")
    fig.suptitle(
        "Random Forest: Parental Marks Predict Compartment Switching"
        "\nCompartment-resolved (A and B)",
        fontweight="bold", fontsize=12, y=1.00,
    )
    plt.tight_layout()
    save_fig(fig, outdir, "Fig1c_RF_ROC_PR_both")


# =============================================================================
# Fig 1d — SHAP (mode-aware)
# =============================================================================

def _shap_values(rf, X):
    explainer = shap.TreeExplainer(rf)
    shap_vals = explainer.shap_values(X)
    if isinstance(shap_vals, list):
        sv = shap_vals[1]
    elif np.ndim(shap_vals) == 3:
        sv = shap_vals[:, :, 1]
    else:
        sv = shap_vals
    return sv


def _draw_shap_panels(ax_bar, ax_bee, sv, X, mark_names, title_extra=""):
    """Draw a (bar, beeswarm) pair on the provided axes."""
    mean_shap = np.abs(sv).mean(axis=0).ravel()
    order     = np.argsort(mean_shap)[::-1]
    top_idx   = int(np.argmax(mean_shap))

    bar_colors = ["#C0392B" if int(i) == top_idx else "#1A5276"
                  for i in order]
    bars = ax_bar.barh(
        [mark_names[i] for i in order[::-1]],
        [mean_shap[i]  for i in order[::-1]],
        color=bar_colors[::-1], edgecolor="none", height=0.6,
    )
    for bar, val in zip(bars, [mean_shap[i] for i in order[::-1]]):
        ax_bar.text(bar.get_width() + max(mean_shap) * 0.012,
                    bar.get_y() + bar.get_height() / 2,
                    f"{val:.3f}", va="center", ha="left", fontsize=9)
    ax_bar.set_xlabel("Mean |SHAP value|\n(magnitude of contribution)",
                      labelpad=6)
    ax_bar.set_title("Feature Importance (magnitude)" + title_extra,
                     fontweight="bold", pad=8)
    ax_bar.set_xlim(0, max(mean_shap) * 1.18)

    # Beeswarm
    feat_order = [mark_names[i] for i in order][::-1]
    X_df  = pd.DataFrame(X,  columns=mark_names)[feat_order]
    sv_df = pd.DataFrame(sv, columns=mark_names)[feat_order]
    plt.sca(ax_bee)
    shap.summary_plot(
        sv_df.values, X_df.values,
        feature_names=feat_order,
        plot_type="dot", show=False, color_bar=True,
        sort=False, plot_size=None,
    )
    ax_bee.set_xlabel("SHAP value (impact on switcher probability)",
                      labelpad=6)
    ax_bee.set_title("Feature Importance (direction)" + title_extra,
                     fontweight="bold", pad=8)


def plot_shap(results_by_mode, mark_names, X_by_mode, outdir, comparison):
    """
    SHAP bar + beeswarm.
        single-mode (A/B/all): 1x2 figure  -> [bar | beeswarm]
        both:                  2x2 figure  -> A on top row, B on bottom row
    Returns dict {mode: (sv, X, mark_names)} for downstream dependence plots.
    """
    sv_store = {}

    if comparison in ("A", "B", "all"):
        res   = results_by_mode[comparison]
        X     = X_by_mode[comparison]
        sv    = _shap_values(res["rf_refit"], X)
        sv_store[comparison] = (sv, X)

        fig, (ax_bar, ax_bee) = plt.subplots(
            1, 2, figsize=(11, 4.5),
            gridspec_kw={"width_ratios": [1.0, 1.2]},
        )
        title_extra = ""
        if comparison == "all":
            mode_label = "Switchers vs Stable"
        else:
            mode_label = COMPARISON_SPEC[comparison]["label"]
        _draw_shap_panels(ax_bar, ax_bee, sv, X, mark_names)
        fig.suptitle(
            f"Parental Marks \u2192 Compartment Switch Prediction (SHAP)"
            f"\n({mode_label})",
            fontweight="bold", fontsize=12, y=1.02,
        )
        plt.tight_layout()
        save_fig(fig, outdir, f"Fig1d_SHAP_{comparison}")
        return sv_store

    # both -> 2x2
    fig, axes = plt.subplots(
        2, 2, figsize=(11, 9),
        gridspec_kw={"width_ratios": [1.0, 1.2]},
    )
    for r_idx, comp in enumerate(["A", "B"]):
        res   = results_by_mode[comp]
        X     = X_by_mode[comp]
        sv    = _shap_values(res["rf_refit"], X)
        sv_store[comp] = (sv, X)
        spec = COMPARISON_SPEC[comp]
        _draw_shap_panels(
            axes[r_idx, 0], axes[r_idx, 1],
            sv, X, mark_names,
            title_extra=f"\n({spec['label']})",
        )
    fig.suptitle(
        "Parental Marks \u2192 Compartment Switch Prediction (SHAP)"
        "\nCompartment-resolved (A on top, B on bottom)",
        fontweight="bold", fontsize=12, y=1.00,
    )
    plt.tight_layout()
    save_fig(fig, outdir, "Fig1d_SHAP_both")
    return sv_store


# =============================================================================
# SHAP dependence plots
# =============================================================================

def _draw_dependence_panel(ax, x_vals, sv_vals, c_vals, feat_x, feat_c,
                            xlim=None, ylim=None):
    """Helper: draw a single dependence scatter on `ax`."""
    sc = ax.scatter(
        x_vals, sv_vals, c=c_vals,
        cmap="plasma", s=8, alpha=0.55,
        edgecolors="none", rasterized=True,
    )
    ax.axhline(0, color="black", lw=0.6, alpha=0.5)
    ax.set_xlabel(f"log1p({feat_x}) (Parental CPM signal)", labelpad=6)
    ax.set_ylabel(
        f"SHAP value for {feat_x}\n(impact on switcher probability)",
        labelpad=6,
    )
    ax.set_title(
        f"SHAP Dependence: {feat_x}\n(colored by {feat_c})",
        fontweight="bold", pad=8,
    )
    if xlim is not None: ax.set_xlim(*xlim)
    if ylim is not None: ax.set_ylim(*ylim)
    cbar = plt.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label(f"log1p({feat_c})", fontsize=9)
    cbar.ax.tick_params(labelsize=8)


def _resolve_pair(pair, mark_names):
    """
    Case-insensitive mark resolution.
    Returns (feat_a, feat_b, idx_a, idx_b) or None if either is missing.
    """
    name_to_idx_ci = {m.lower(): i for i, m in enumerate(mark_names)}
    feat_a_in, feat_b_in = pair
    feat_a_idx = name_to_idx_ci.get(feat_a_in.lower())
    feat_b_idx = name_to_idx_ci.get(feat_b_in.lower())
    if feat_a_idx is None or feat_b_idx is None:
        print(f"  [WARN] '{feat_a_in}' or '{feat_b_in}' not in {mark_names}; "
              f"skipping.")
        return None
    return (mark_names[feat_a_idx], mark_names[feat_b_idx],
            feat_a_idx, feat_b_idx)


def plot_shap_dependence(sv_store, mark_names, pairs, outdir, comparison,
                          zoom_pct=99):
    """
    For each requested pair, produce TWO figures: full-range and zoomed
    (clipped to 1st-99th percentile by default).

    Layouts:
        A or B (single):   1x2  -> [SHAP(A) vs A col by B | SHAP(B) vs B col by A]
        both:              2x2  -> top row A mode, bottom row B mode
        all:               1x2  -> like single, run on full dataset
    """
    if not pairs:
        print("[INFO] No --shap_dependence_partners specified; skipping "
              "dependence plots.")
        return

    def _draw_for_mode(ax_l, ax_r, sv, X, feat_a, feat_b, idx_a, idx_b,
                        xlim_a=None, ylim_a=None, xlim_b=None, ylim_b=None):
        _draw_dependence_panel(
            ax_l, np.log1p(X[:, idx_a]), sv[:, idx_a],
            np.log1p(X[:, idx_b]), feat_a, feat_b, xlim_a, ylim_a,
        )
        _draw_dependence_panel(
            ax_r, np.log1p(X[:, idx_b]), sv[:, idx_b],
            np.log1p(X[:, idx_a]), feat_b, feat_a, xlim_b, ylim_b,
        )

    def _compute_pct(arr, pct):
        lo = np.percentile(arr, 100 - pct)
        hi = np.percentile(arr, pct)
        pad = 0.04 * (hi - lo) if hi > lo else 0
        return (lo - pad, hi + pad)

    for pair in pairs:
        resolved = _resolve_pair(pair, mark_names)
        if resolved is None:
            continue
        feat_a, feat_b, idx_a, idx_b = resolved
        print(f"[INFO] SHAP dependence pair ({comparison}): "
              f"{feat_a} <-> {feat_b}")

        # ── Full-range figure ────────────────────────────────────────────
        if comparison in ("A", "B", "all"):
            sv, X = sv_store[comparison]
            fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(12, 4.8))
            _draw_for_mode(ax_l, ax_r, sv, X, feat_a, feat_b, idx_a, idx_b)
            mode_label = ("Switchers vs Stable" if comparison == "all"
                          else COMPARISON_SPEC[comparison]["label"])
            fig.suptitle(f"({mode_label})", fontweight="bold",
                         fontsize=10, y=1.02)
            plt.tight_layout()
            save_fig(fig, outdir,
                     f"Fig1d_dependence_{feat_a}_vs_{feat_b}_{comparison}")
        else:  # both -> 2x2
            fig, axes = plt.subplots(2, 2, figsize=(12, 9.6))
            for r_idx, comp in enumerate(["A", "B"]):
                sv, X = sv_store[comp]
                _draw_for_mode(axes[r_idx, 0], axes[r_idx, 1],
                                sv, X, feat_a, feat_b, idx_a, idx_b)
                spec = COMPARISON_SPEC[comp]
                # Subtle row-level annotation via the panel titles' extra
                axes[r_idx, 0].set_title(
                    axes[r_idx, 0].get_title() + f"\n[{spec['label']}]",
                    fontweight="bold", pad=8,
                )
                axes[r_idx, 1].set_title(
                    axes[r_idx, 1].get_title() + f"\n[{spec['label']}]",
                    fontweight="bold", pad=8,
                )
            fig.suptitle("SHAP Dependence — A on top, B on bottom",
                         fontweight="bold", fontsize=12, y=1.00)
            plt.tight_layout()
            save_fig(fig, outdir,
                     f"Fig1d_dependence_{feat_a}_vs_{feat_b}_both")

        # ── Zoomed figure (clip to 1st – zoom_pct percentile) ────────────
        # Limits are computed per-mode so each panel zooms to its own bulk.
        if comparison in ("A", "B", "all"):
            sv, X = sv_store[comparison]
            xlim_a = _compute_pct(np.log1p(X[:, idx_a]), zoom_pct)
            ylim_a = _compute_pct(sv[:, idx_a],          zoom_pct)
            xlim_b = _compute_pct(np.log1p(X[:, idx_b]), zoom_pct)
            ylim_b = _compute_pct(sv[:, idx_b],          zoom_pct)
            fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(12, 4.8))
            _draw_for_mode(ax_l, ax_r, sv, X,
                            feat_a, feat_b, idx_a, idx_b,
                            xlim_a, ylim_a, xlim_b, ylim_b)
            mode_label = ("Switchers vs Stable" if comparison == "all"
                          else COMPARISON_SPEC[comparison]["label"])
            fig.suptitle(
                f"({mode_label})  —  zoomed to {zoom_pct}th percentile",
                fontweight="bold", fontsize=10, y=1.02,
            )
            plt.tight_layout()
            save_fig(fig, outdir,
                     f"Fig1d_dependence_{feat_a}_vs_{feat_b}_"
                     f"{comparison}_zoomed")
        else:
            fig, axes = plt.subplots(2, 2, figsize=(12, 9.6))
            for r_idx, comp in enumerate(["A", "B"]):
                sv, X = sv_store[comp]
                xlim_a = _compute_pct(np.log1p(X[:, idx_a]), zoom_pct)
                ylim_a = _compute_pct(sv[:, idx_a],          zoom_pct)
                xlim_b = _compute_pct(np.log1p(X[:, idx_b]), zoom_pct)
                ylim_b = _compute_pct(sv[:, idx_b],          zoom_pct)
                _draw_for_mode(axes[r_idx, 0], axes[r_idx, 1],
                                sv, X,
                                feat_a, feat_b, idx_a, idx_b,
                                xlim_a, ylim_a, xlim_b, ylim_b)
                spec = COMPARISON_SPEC[comp]
                axes[r_idx, 0].set_title(
                    axes[r_idx, 0].get_title() + f"\n[{spec['label']}]",
                    fontweight="bold", pad=8,
                )
                axes[r_idx, 1].set_title(
                    axes[r_idx, 1].get_title() + f"\n[{spec['label']}]",
                    fontweight="bold", pad=8,
                )
            fig.suptitle(
                f"SHAP Dependence (zoomed to {zoom_pct}th percentile)"
                "\nA on top, B on bottom",
                fontweight="bold", fontsize=12, y=1.00,
            )
            plt.tight_layout()
            save_fig(fig, outdir,
                     f"Fig1d_dependence_{feat_a}_vs_{feat_b}_both_zoomed")


def parse_dependence_pairs(raw_pairs):
    """Parse 'MARK1:MARK2' CLI strings into list of tuples."""
    if not raw_pairs:
        return []
    out = []
    for s in raw_pairs:
        if ":" not in s:
            print(f"[WARN] '{s}' missing ':'; skipping.")
            continue
        parts = s.split(":")
        if len(parts) != 2 or not parts[0].strip() or not parts[1].strip():
            print(f"[WARN] '{s}' malformed; skipping.")
            continue
        out.append((parts[0].strip(), parts[1].strip()))
    return out


# =============================================================================
# Fig 1e — PCA (mode-aware)
# =============================================================================

def _draw_pca_panel(ax, df_subset, mark_names, cats_fg, cats_bg):
    """
    Draw one PCA on `ax`. cats_fg are foreground (visible, larger);
    cats_bg are faded background.
    """
    if df_subset.empty:
        ax.set_visible(False); return None

    scaler = StandardScaler()
    Xsc    = scaler.fit_transform(np.log1p(df_subset[mark_names].values))
    pca    = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(Xsc)
    var    = pca.explained_variance_ratio_ * 100

    for cat in cats_bg:
        mask = df_subset["switch_cat"].values == cat
        if mask.sum() == 0: continue
        ax.scatter(coords[mask, 0], coords[mask, 1],
                   c=SWITCH_COLORS[cat], s=4, alpha=0.06,
                   linewidths=0, rasterized=True, label=SWITCH_LABELS[cat])
    for cat in cats_fg:
        mask = df_subset["switch_cat"].values == cat
        if mask.sum() == 0: continue
        ax.scatter(coords[mask, 0], coords[mask, 1],
                   c=SWITCH_COLORS[cat], s=18, alpha=0.55,
                   linewidths=0, rasterized=True, label=SWITCH_LABELS[cat])

    xlo, xhi = np.percentile(coords[:, 0], [1, 99])
    ylo, yhi = np.percentile(coords[:, 1], [1, 99])
    pad_x = 0.05 * (xhi - xlo); pad_y = 0.05 * (yhi - ylo)
    ax.set_xlim(xlo - pad_x, xhi + pad_x)
    ax.set_ylim(ylo - pad_y, yhi + pad_y)

    ax.set_xlabel(f"PC1 ({var[0]:.1f}% variance)", labelpad=6)
    ax.set_ylabel(f"PC2 ({var[1]:.1f}% variance)", labelpad=6)
    legend_handles = [
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=SWITCH_COLORS[c],
                   markersize=8, label=SWITCH_LABELS[c])
        for c in (cats_fg + cats_bg)
        if c in df_subset["switch_cat"].values
    ]
    ax.legend(handles=legend_handles, frameon=False, fontsize=8.5, loc="best")
    return pca


def plot_pca(df, mark_names, outdir, comparison):
    """
    PCA of log1p(Parental marks).
        A / B (single): 1 panel, on filtered subset
        both:           1x2 panel (A on left, B on right)
        all:            unchanged from old behaviour (full 4-category)
    """
    if comparison == "all":
        df_sub = df[df["switch_cat"].isin(
            ["Stable_A", "Stable_B", "A_to_B", "B_to_A"]
        )].dropna(subset=mark_names)
        fig, ax = plt.subplots(figsize=(6.2, 5.6))
        _draw_pca_panel(ax, df_sub, mark_names,
                         cats_fg=["A_to_B", "B_to_A"],
                         cats_bg=["Stable_A", "Stable_B"])
        ax.set_title(
            "PCA: Parental Epigenetic Mark Space (log1p)"
            "\nby Compartment Switch Outcome",
            fontweight="bold", pad=10,
        )
        plt.tight_layout()
        save_fig(fig, outdir, "Fig1e_PCA_marks_all")
        return

    if comparison in ("A", "B"):
        spec = COMPARISON_SPEC[comparison]
        df_sub = df[df["switch_cat"].isin([spec["neg"], spec["pos"]])] \
                    .dropna(subset=mark_names)
        fig, ax = plt.subplots(figsize=(6.2, 5.6))
        _draw_pca_panel(ax, df_sub, mark_names,
                         cats_fg=[spec["pos"]],
                         cats_bg=[spec["neg"]])
        ax.set_title(
            f"PCA: Parental Epigenetic Mark Space (log1p)"
            f"\n{spec['label']}",
            fontweight="bold", pad=10,
        )
        plt.tight_layout()
        save_fig(fig, outdir, f"Fig1e_PCA_marks_{comparison}")
        return

    # both -> 1x2
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.6))
    for ax, comp in zip(axes, ["A", "B"]):
        spec = COMPARISON_SPEC[comp]
        df_sub = df[df["switch_cat"].isin([spec["neg"], spec["pos"]])] \
                    .dropna(subset=mark_names)
        _draw_pca_panel(ax, df_sub, mark_names,
                         cats_fg=[spec["pos"]],
                         cats_bg=[spec["neg"]])
        ax.set_title(f"PCA: {spec['label']}", fontweight="bold", pad=10)
    fig.suptitle(
        "PCA: Parental Epigenetic Mark Space (log1p)"
        "\nCompartment-resolved (A on left, B on right)",
        fontweight="bold", fontsize=12, y=1.00,
    )
    plt.tight_layout()
    save_fig(fig, outdir, "Fig1e_PCA_marks_both")


# =============================================================================
# Statistics output
# =============================================================================

def save_statistics(df, mark_names, outdir, comparison):
    """
    Per-mark Mann-Whitney + Cliff's delta, with the comparison column
    indicating which contrast was used.
        A:    rows for Stable_A vs A->B per mark
        B:    rows for Stable_B vs B->A per mark
        both: A and B rows interleaved
        all:  original Switchers (B->A ∪ A->B) vs Stable (Stable_A ∪ Stable_B)
              (Dunn's tables are saved separately by plot_violin_marks_all.)
    """
    def _per_mark(neg_cats, pos_cats, label):
        out = []
        neg_mask = df["switch_cat"].isin(neg_cats)
        pos_mask = df["switch_cat"].isin(pos_cats)
        for mark in mark_names:
            neg_vals = df[neg_mask][mark].dropna().values
            pos_vals = df[pos_mask][mark].dropna().values
            if len(neg_vals) == 0 or len(pos_vals) == 0:
                continue
            stat, p = mannwhitneyu(pos_vals, neg_vals, alternative="two-sided")
            d, mag  = cliffs_delta_fast(pos_vals, neg_vals)
            out.append({
                "comparison":         label,
                "mark":               mark,
                "n_positive":         len(pos_vals),
                "n_negative":         len(neg_vals),
                "median_positive":    np.median(pos_vals),
                "median_negative":    np.median(neg_vals),
                "mannwhitney_U":      stat,
                "p_value":            p,
                "cliffs_delta":       d,   # positive vs negative
                "effect_magnitude":   mag,
            })
        return out

    rows = []
    if comparison == "A":
        rows += _per_mark(["Stable_A"], ["A_to_B"], "Stable_A_vs_AtoB")
    elif comparison == "B":
        rows += _per_mark(["Stable_B"], ["B_to_A"], "Stable_B_vs_BtoA")
    elif comparison == "both":
        rows += _per_mark(["Stable_A"], ["A_to_B"], "Stable_A_vs_AtoB")
        rows += _per_mark(["Stable_B"], ["B_to_A"], "Stable_B_vs_BtoA")
    elif comparison == "all":
        rows += _per_mark(["Stable_A", "Stable_B"],
                          ["B_to_A", "A_to_B"],
                          "Switchers_vs_Stable")

    if not rows:
        print("[WARN] No statistics rows generated.")
        return

    res_df = pd.DataFrame(rows)
    res_df["padj_BH"] = multipletests(res_df["p_value"], method="fdr_bh")[1]
    out = os.path.join(outdir, f"Statistics_{comparison}.txt")
    res_df.to_csv(out, sep="\t", index=False)
    print(f"[INFO] Statistics saved -> {out}")
    print(res_df[["comparison", "mark", "median_positive", "median_negative",
                  "cliffs_delta", "effect_magnitude",
                  "p_value", "padj_BH"]].to_string())


# =============================================================================
# Orchestrator
# =============================================================================

def run_ml_for_comparison(df, mark_names, comparison):
    """
    Build X, y for the requested mode and run the RF + CV.
    Returns (results_dict, X_dict) keyed by sub-comparison.
        single mode 'A'/'B'/'all' -> dict with 1 key.
        'both'                    -> dict with 2 keys ('A' and 'B').
    """
    results = {}
    Xs      = {}

    if comparison == "all":
        sub = df[df["switch_cat"].isin(
            ["B_to_A", "A_to_B", "Stable_A", "Stable_B"]
        )].copy()
        sub["is_switcher"] = sub["switch_cat"].isin(
            ["B_to_A", "A_to_B"]
        ).astype(int)
        X = sub[mark_names].values
        y = sub["is_switcher"].values
        if len(X) < 50 or y.sum() < 5:
            print("[WARN] Insufficient data for RF (all). Skipping.")
            return results, Xs
        print(f"[INFO] RF (all): N={len(X)}, switchers={y.sum()}, "
              f"stable={(y == 0).sum()}, prevalence={y.mean():.3f}")
        results["all"] = _run_rf_cv(X, y)
        Xs["all"]      = X
        return results, Xs

    comps = ["A", "B"] if comparison == "both" else [comparison]
    for comp in comps:
        sub, spec = subset_for_comparison(df, comp)
        if sub["is_switcher"].sum() < 5 or len(sub) < 50:
            print(f"[WARN] Insufficient data for RF ({comp}). Skipping.")
            continue
        X = sub[mark_names].values
        y = sub["is_switcher"].values
        print(f"[INFO] RF ({spec['label']}): N={len(X)}, "
              f"positive={y.sum()} ({spec['pos']}), "
              f"negative={(y == 0).sum()} ({spec['neg']}), "
              f"prevalence={y.mean():.3f}")
        results[comp] = _run_rf_cv(X, y)
        Xs[comp]      = X
    return results, Xs


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if args.comparison == "all" and sp is None:
        print("[WARN] 'all' mode requested but scikit-posthocs is not "
              "installed. Dunn's post-hoc and brackets will be skipped.")
        print("       Install with: pip install scikit-posthocs")

    df, mark_names = load_data(args.cutrun_file, args.merged_file,
                                include_ctcf=args.include_ctcf)

    # Always emit the classified gene table
    out_genes = os.path.join(args.outdir, "Switch_classified.txt")
    df.to_csv(out_genes, sep="\t", index=False)
    print(f"[INFO] Classified gene table saved -> {out_genes}")

    print(f"\n[INFO] Comparison mode: {args.comparison}")
    print("\n[INFO] Fig 1a: switch distribution...")
    plot_switch_distribution(df, args.outdir, args.comparison)

    print("\n[INFO] Fig 1b: violins...")
    if args.comparison == "all":
        plot_violin_marks_all(df, mark_names, args.outdir)
    else:
        plot_violin_marks_AB(df, mark_names, args.outdir, args.comparison)

    print("\n[INFO] ML pipeline (RF + CV)...")
    results_by_mode, X_by_mode = run_ml_for_comparison(
        df, mark_names, args.comparison
    )
    if not results_by_mode:
        print("[ERROR] No RF results produced; aborting downstream plots.")
        sys.exit(1)

    print("\n[INFO] Fig 1c: ROC + PR...")
    plot_roc_pr(results_by_mode, args.outdir, args.comparison)

    print("\n[INFO] Fig 1d: SHAP bar + beeswarm...")
    sv_store = plot_shap(results_by_mode, mark_names, X_by_mode,
                          args.outdir, args.comparison)

    pairs = parse_dependence_pairs(args.shap_dependence_partners)
    if pairs:
        print(f"\n[INFO] Fig 1d: SHAP dependence pairs ({pairs})...")
        plot_shap_dependence(sv_store, mark_names, pairs,
                              args.outdir, args.comparison,
                              zoom_pct=99)

    print("\n[INFO] Fig 1e: PCA...")
    plot_pca(df, mark_names, args.outdir, args.comparison)

    print("\n[INFO] Statistics...")
    save_statistics(df, mark_names, args.outdir, args.comparison)

    print(f"\n[INFO] Q1 complete -> {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()