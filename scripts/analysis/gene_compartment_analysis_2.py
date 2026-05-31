"""
=============================================================================
gene_compartment_analysis_2.py  — STEP 2: gene density across switch categories
=============================================================================

GOAL
  Take the per-bin gene-density table from Step 1 and test whether gene density
  differs across the four compartment-switch categories, with the SAME
  statistical framework used by RF_SHAP.py's violins (so the whole study speaks
  one statistical language):
      - Kruskal-Wallis omnibus (non-parametric ANOVA; gene counts are skewed
        counts, not normal, so KW is the correct omnibus, not one-way ANOVA).
      - Dunn's post-hoc with Benjamini-Hochberg FDR correction (which pairs of
        categories differ).
      - Cliff's delta effect size per pair (magnitude + direction, robust to
        the skew; same estimator as RF_SHAP.py).

WHY THIS IS A REAL RESULT (not redundant)
  Gene density is the DEFINING property that distinguishes compartment A from B
  (the standard A/B caller assigns 'A' to the higher-gene-density sign). Showing
  the density gradient across switch categories (expected: Stable_A >> others,
  Stable_B lowest) both (a) validates the compartment calls internally and
  (b) characterizes switching regions as an intermediate, gene-density class.
  This is the 'structure' half of the structure->function bridge; Step 3 adds
  expression, Step 4 adds histone-mark change.

INPUT
  Step1_bin_gene_density.txt  (from Step 1): columns
      coordinate, chr, start, end, switch_cat, gene_count, gene_names

OUTPUTS (into --outdir)
  Step2_density_by_switch.pdf/png    violin+box of gene_count by category
  Step2_KW_omnibus.txt               Kruskal-Wallis H, df, p
  Step2_dunn_posthoc.txt             pairwise Dunn BH-adjusted p + Cliff's delta
  Step2_density_summary.txt          n_bins, mean/median/total genes per category

USAGE
  python gene_compartment_step2_density_stats.py \
      --step1 results/exploratory/gene_compartment/Step1_bin_gene_density.txt \
      --outdir results/exploratory/gene_compartment
"""

import argparse
import os
import sys
from itertools import combinations

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import kruskal

try:
    import scikit_posthocs as sp
except ImportError:
    sp = None


# ── Aesthetics — match RF_SHAP.py palette/labels for a consistent study ──────
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         11,
    "axes.labelsize":    12,
    "axes.titlesize":    13,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "pdf.fonttype":      42,
    "ps.fonttype":       42,
    "savefig.dpi":       600,
    "figure.facecolor":  "white",
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
# Order: gene-rich (Stable_A) -> gene-poor (Stable_B), switchers in between.
CAT_ORDER = ["Stable_A", "A_to_B", "B_to_A", "Stable_B"]


# ── Cliff's delta (identical estimator/bins to RF_SHAP.py) ───────────────────
def cliffs_delta_fast(x, y):
    """Cliff's delta of x vs y. Positive => x tends to exceed y."""
    x = np.asarray(x); y = np.asarray(y)
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return np.nan, "NA"
    ys = np.sort(y)
    lt = np.searchsorted(ys, x, side="left")   # y-values strictly below each x
    le = np.searchsorted(ys, x, side="right")
    gt = ny - le                               # y-values strictly above each x
    # Standard Cliff's delta convention: positive => x tends to EXCEED y.
    # (Earlier code used (gt - lt), which inverts the sign — fixed here.)
    d = (int(np.sum(lt)) - int(np.sum(gt))) / (nx * ny)
    ad = abs(d)
    if   ad < 0.147: mag = "negligible"
    elif ad < 0.330: mag = "small"
    elif ad < 0.474: mag = "medium"
    else:            mag = "large"
    return d, mag


def _p_to_stars(p):
    if p < 1e-4: return "****"
    if p < 1e-3: return "***"
    if p < 1e-2: return "**"
    if p < 5e-2: return "*"
    return "ns"


def _bracket(ax, x1, x2, y, label, lw=1.0, fs=9):
    dy = (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.012
    ax.plot([x1, x1, x2, x2], [y, y + dy, y + dy, y], lw=lw, c="black")
    ax.text((x1 + x2) / 2, y + dy * 1.2, label, ha="center", va="bottom", fontsize=fs)


def parse_args():
    p = argparse.ArgumentParser(
        description="Step 2: gene density across switch categories (KW + Dunn).")
    p.add_argument("--step1", required=True,
                   help="Step1_bin_gene_density.txt from Step 1.")
    p.add_argument("--outdir", default="results/exploratory/gene_compartment")
    p.add_argument("--log1p", action="store_true",
                   help="Plot log1p(gene_count) on the y-axis (counts are "
                        "right-skewed; log1p matches the RF_SHAP violins). "
                        "Stats are ALWAYS run on raw counts regardless.")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print("[INFO] Loading Step 1 density table...")
    df = pd.read_csv(args.step1, sep="\t", low_memory=False)
    for col in ("switch_cat", "gene_count"):
        if col not in df.columns:
            sys.exit(f"[ERROR] '{col}' not in {args.step1}. Re-run Step 1.")
    df = df[df["switch_cat"].isin(CAT_ORDER)].copy()
    df["gene_count"] = pd.to_numeric(df["gene_count"], errors="coerce")
    df = df.dropna(subset=["gene_count"])
    cats = [c for c in CAT_ORDER if c in df["switch_cat"].unique()]
    print(f"[INFO] {len(df):,} bins across {len(cats)} categories: {cats}")

    # ── Per-category data (raw counts for stats) ─────────────────────────────
    data = {c: df.loc[df["switch_cat"] == c, "gene_count"].values for c in cats}

    # ── Summary table ────────────────────────────────────────────────────────
    summ = (df.groupby("switch_cat")["gene_count"]
              .agg(n_bins="size", mean_genes="mean", median_genes="median",
                   total_genes="sum")
              .reindex(cats))
    summ.to_csv(os.path.join(args.outdir, "Step2_density_summary.txt"), sep="\t")
    print("[INFO] Density summary:\n", summ.to_string())

    # ── Kruskal-Wallis omnibus (raw counts) ──────────────────────────────────
    kw_H, kw_p = kruskal(*[data[c] for c in cats])
    with open(os.path.join(args.outdir, "Step2_KW_omnibus.txt"), "w") as fh:
        fh.write("test\tH_statistic\tdf\tp_value\n")
        fh.write(f"Kruskal-Wallis\t{kw_H:.4f}\t{len(cats)-1}\t{kw_p:.3e}\n")
    print(f"[INFO] Kruskal-Wallis: H={kw_H:.2f}, df={len(cats)-1}, p={kw_p:.3e}")

    # ── Dunn's post-hoc (BH) + Cliff's delta per pair ────────────────────────
    dunn_rows = []
    if sp is not None:
        long_df = df[["switch_cat", "gene_count"]].copy()
        dunn = sp.posthoc_dunn(long_df, val_col="gene_count",
                               group_col="switch_cat", p_adjust="fdr_bh")
        for c1, c2 in combinations(cats, 2):
            p_pair = float(dunn.loc[c1, c2])
            # Cliff's delta of c1 vs c2 (positive => c1 denser than c2)
            d_pair, mag = cliffs_delta_fast(data[c1], data[c2])
            dunn_rows.append({
                "group1": c1, "group2": c2, "dunn_p_BH": p_pair,
                "cliffs_delta_g1_vs_g2": d_pair, "effect_magnitude": mag,
                "stars": _p_to_stars(p_pair),
            })
        pd.DataFrame(dunn_rows).to_csv(
            os.path.join(args.outdir, "Step2_dunn_posthoc.txt"),
            sep="\t", index=False)
        print("[INFO] Dunn post-hoc (BH):")
        print(pd.DataFrame(dunn_rows).to_string(index=False))
    else:
        print("[WARN] scikit-posthocs not installed; skipping Dunn post-hoc "
              "and significance brackets. Install with: pip install "
              "scikit-posthocs --break-system-packages")

    # ── Figure: violin + box of gene density by category ─────────────────────
    use_log = args.log1p
    def transform(v):
        return np.log1p(v) if use_log else v
    plot_data = [transform(data[c]) for c in cats]

    fig, ax = plt.subplots(figsize=(7, 5.5))
    parts = ax.violinplot(plot_data, showmedians=False, showextrema=False,
                          widths=0.75)
    for pc, c in zip(parts["bodies"], cats):
        pc.set_facecolor(SWITCH_COLORS[c]); pc.set_alpha(0.75)
        pc.set_edgecolor("none")
    # Box overlay + median tick
    bp = ax.boxplot(plot_data, widths=0.18, showfliers=False, patch_artist=True,
                    medianprops=dict(color="black", linewidth=1.8),
                    boxprops=dict(facecolor="white", alpha=0.9, linewidth=0.8),
                    whiskerprops=dict(linewidth=0.8),
                    capprops=dict(linewidth=0.8))

    ax.set_xticks(range(1, len(cats) + 1))
    ax.set_xticklabels([SWITCH_LABELS[c] for c in cats], rotation=20, ha="right")
    ax.set_ylabel("log1p(genes per 250kb bin)" if use_log
                  else "Genes per 250kb bin", labelpad=6)
    ax.set_title("Gene density across compartment-switch categories\n"
                 f"Kruskal-Wallis H = {kw_H:.1f}, "
                 + ("p < 1e-300" if kw_p == 0 or kw_p < 1e-300
                    else f"p = {kw_p:.1e}"),
                 fontweight="bold", pad=10)

    # Significance brackets for the significant pairs (closest pairs lowest).
    if dunn_rows:
        cat_pos = {c: i + 1 for i, c in enumerate(cats)}
        sig = [r for r in dunn_rows if r["dunn_p_BH"] < 0.05]
        sig.sort(key=lambda r: abs(cat_pos[r["group1"]] - cat_pos[r["group2"]]))
        top = max(np.max(d) for d in plot_data)
        tick = (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.07
        base_y = top + tick * 0.6
        for k, r in enumerate(sig):
            # Show significance AND effect size together — stars alone overstate
            # trivial-but-significant pairs at large n; delta gives honest magnitude.
            label = f"{r['stars']}  \u03b4={r['cliffs_delta_g1_vs_g2']:+.2f} ({r['effect_magnitude']})"
            _bracket(ax, cat_pos[r["group1"]], cat_pos[r["group2"]],
                     base_y + k * tick, label)
        ax.set_ylim(ax.get_ylim()[0], base_y + (len(sig) + 1) * tick)

    plt.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(args.outdir, f"Step2_density_by_switch.{ext}")
        fig.savefig(out, bbox_inches="tight", pad_inches=0.05, dpi=600)
        print(f"  Saved -> {out}")
    plt.close(fig)

    print("\n[INFO] Step 2 complete. If Stable_A is significantly denser than "
          "the others (expected), proceed to Step 3 (expression by switch "
          "direction).")


if __name__ == "__main__":
    main()