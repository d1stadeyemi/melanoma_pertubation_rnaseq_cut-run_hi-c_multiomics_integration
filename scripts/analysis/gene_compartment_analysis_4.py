"""
=============================================================================
gene_compartment_analysis.py  — STEP 4: histone-mark change overlay
=============================================================================

GOAL — the integrative payoff
  Tie together STRUCTURE (compartment switch), CHROMATIN (histone-mark change),
  and FUNCTION (expression change) at 250kb resolution:
    4a  Does histone-mark CHANGE (B10 vs Parental) differ across switch
        categories?  e.g. do A->B (silenced) bins gain repressive / lose active
        marks, and B->A (activated) bins the reverse?
    4b  In switching bins, do mark-change and expression-change move
        concordantly?  (the three-way association figure)

  Frame as ASSOCIATION, not causation. All of this is observational (both
  conditions measured post-migration); the causal arm is the UNC/EZH2i
  perturbation, not this comparison.

=== THE PSEUDOREPLICATION RULE (critical) ===
  Histone-mark change is a BIN-LEVEL quantity: every gene in a 250kb bin shares
  the SAME mark-change value. Therefore:
    * MARK-CHANGE stats are computed at the BIN level (one row per bin, each
      bin contributes once). Treating genes as independent here would inflate
      significance ~N-fold (N = genes per bin) — textbook pseudoreplication.
    * EXPRESSION stats stay at the GENE level (expression is genuinely per-gene;
      genes are independent units there — handled in Step 3, reused in 4b by
      AGGREGATING genes to a per-bin mean before relating to mark-change).
  This separation is the whole reason Step 4 is split into 4a (bin) and 4b
  (bin-aggregated). A reviewer will check exactly this.

INPUTS
  --step1      Step1_bin_gene_density.txt        (coordinate + switch_cat)
  --cutrun     250kb CUT&RUN table               (coordinate + log2ratio_p-b10_<mark>_mean)
  --step3      Step3_gene_expression_joined.txt   (per-gene lfc; OPTIONAL, enables 4b)

OUTPUTS (into --outdir)
  Step4a_markchange_by_switch.pdf/png     per-mark violin of bin mark-change by category
  Step4a_markchange_stats.txt             per-mark KW + Dunn(BH) + Cliff's delta (bin-level)
  Step4b_mark_vs_expression.pdf/png       (if --step3) per switching bin: mark-change
                                          vs mean expression-change, one dot per bin
  Step4b_bin_integrated.txt               (if --step3) per-bin: switch_cat, each mark
                                          change, mean gene expression change, n_genes

USAGE
  python gene_compartment_step4_markchange.py \
      --step1   results/exploratory/gene_compartment/Step1_bin_gene_density.txt \
      --cutrun  data/processed/cut_run_signal_table.txt \
      --step3   results/exploratory/gene_compartment/Step3_gene_expression_joined.txt \
      --outdir  results/exploratory/gene_compartment
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


# ── Aesthetics — consistent with Steps 1-3 / RF_SHAP ─────────────────────────
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
SWITCH_COLORS = {"B_to_A": "#C0392B", "A_to_B": "#1A5276",
                 "Stable_A": "#85C1E9", "Stable_B": "#D5D8DC"}
SWITCH_LABELS = {"B_to_A": "B\u2192A (Activated)", "A_to_B": "A\u2192B (Silenced)",
                 "Stable_A": "Stable A", "Stable_B": "Stable B"}
CAT_ORDER = ["Stable_A", "A_to_B", "B_to_A", "Stable_B"]
DEFAULT_MARKS = ("H3K27me3", "H3K27Ac", "H3K4me1", "H3K9me3")


def cliffs_delta_fast(x, y):
    """Cliff's delta of x vs y. Positive => x tends to EXCEED y. (Correct sign.)"""
    x = np.asarray(x); y = np.asarray(y)
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return np.nan, "NA"
    ys = np.sort(y)
    lt = np.searchsorted(ys, x, side="left")
    le = np.searchsorted(ys, x, side="right")
    gt = ny - le
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


def _bracket(ax, x1, x2, y, label, lw=1.0, fs=8):
    dy = (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.012
    ax.plot([x1, x1, x2, x2], [y, y + dy, y + dy, y], lw=lw, c="black")
    ax.text((x1 + x2) / 2, y + dy * 1.2, label, ha="center", va="bottom", fontsize=fs)


def parse_args():
    p = argparse.ArgumentParser(
        description="Step 4: histone-mark change overlay (bin-level + integrative).")
    p.add_argument("--step1", required=True,
                   help="Step1_bin_gene_density.txt (coordinate + switch_cat).")
    p.add_argument("--cutrun", required=True,
                   help="250kb CUT&RUN table with log2ratio_p-b10_<mark>_mean columns.")
    p.add_argument("--step3", default=None,
                   help="Step3_gene_expression_joined.txt (enables 4b integrative view).")
    p.add_argument("--outdir", default="results/exploratory/gene_compartment")
    p.add_argument("--coord_col", default="coordinate",
                   help="Shared coordinate column in step1 [default %default].")
    p.add_argument("--cutrun_coord_col", default=None,
                   help="Coordinate column in cutrun table. Default: auto-detect "
                        "'coordinate' else 'extra_1'.")
    p.add_argument("--marks", nargs="+", default=list(DEFAULT_MARKS),
                   help="Marks to test (default: H3K27me3 H3K27Ac H3K4me1 H3K9me3).")
    p.add_argument("--gene_name_col", default="gene_name")
    p.add_argument("--zscore", action="store_true",
                   help="Z-score each mark's change before plotting so all "
                        "panels share a comparable spread (use only if a shared "
                        "raw axis makes a small-magnitude mark like H3K27me3 too "
                        "flat to see). Stats are ALWAYS run on raw values "
                        "regardless of this flag.")
    return p.parse_args()


def detect_markchange_cols(cutrun_cols, marks):
    """Find migration-axis mark-change columns: log2ratio_p-b10_<mark>_mean.
    EXCLUDES drug-axis (unc) columns — switching is a migration phenotype, so
    only migration log2 ratios are biologically appropriate here."""
    cols_lc = {c: c.lower() for c in cutrun_cols}
    found = {}
    for mark in marks:
        want = f"log2ratio_p-b10_{mark}_mean".lower()
        hits = [c for c, lc in cols_lc.items() if lc == want]
        if hits:
            found[mark] = hits[0]
    return found


def auto_cutrun_coord(cutrun, override):
    if override:
        if override in cutrun.columns:
            return override
        sys.exit(f"[ERROR] --cutrun_coord_col '{override}' not in cutrun table.")
    for cand in ("coordinate", "coordinates", "coord"):
        if cand in cutrun.columns:
            return cand
    if "extra_1" in cutrun.columns:
        return "extra_1"
    sys.exit("[ERROR] No coordinate column found in cutrun table; pass "
             "--cutrun_coord_col explicitly.")


# =============================================================================
# 4a — mark-change by switch category, BIN LEVEL
# =============================================================================
def step4a(bins, mark_cols, outdir, zscore=False):
    """bins: DataFrame with one row per bin: switch_cat + <mark> change columns.
    Runs KW + Dunn + Cliff's delta per mark at the BIN level (no pseudoreplication)."""
    cats = [c for c in CAT_ORDER if c in bins["switch_cat"].unique()]
    marks = list(mark_cols.keys())
    n = len(marks)
    n_cols = min(2, n)
    n_rows = (n + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6.5 * n_cols, 5 * n_rows))
    axes = np.atleast_1d(axes).flatten()
    all_stats = []

    # --- PASS 1: shared y-limit across ALL marks (view only; data untouched) ---
    # Stats below are ALWAYS computed on RAW values; only the DRAWN values and
    # the axis window are affected by --zscore / clipping. This makes the four
    # panels visually comparable without deleting any data point.
    def _draw_vals(mark):
        """Per-category values to PLOT for a mark: raw, or z-scored across all
        bins of that mark (mean/std from the pooled mark column)."""
        col = mark_cols[mark]
        raw = {c: bins.loc[bins["switch_cat"] == c, col].dropna().values
               for c in cats}
        if not zscore:
            return raw
        pooled = bins[col].dropna().values
        mu, sd = float(np.mean(pooled)), float(np.std(pooled))
        if sd == 0:
            return raw
        return {c: (v - mu) / sd for c, v in raw.items()}

    pooled_plot = np.concatenate(
        [np.concatenate([_draw_vals(m)[c] for c in cats]) for m in marks])
    ylo, yhi = np.percentile(pooled_plot, [1, 99])   # clip VIEW to 1-99th pct
    span = (yhi - ylo) if yhi > ylo else 1.0
    n_pairs = len(cats) * (len(cats) - 1) // 2          # 6 brackets max
    SHARED_YLIM = (ylo - 0.40 * span,                   # small bottom pad
                   yhi + (0.60 + 0.085 * n_pairs) * span)  # top headroom scales with #brackets
    ylab = ("standardized mark change (z)" if zscore
            else "mark change  log2(B10/Parental)")

    for ax, mark in zip(axes, marks):
        col = mark_cols[mark]
        # RAW values -> statistics (never clipped, never z-scored)
        data = {c: bins.loc[bins["switch_cat"] == c, col].dropna().values
                for c in cats}
        # PLOT values -> drawing (z-scored if requested)
        pdata = _draw_vals(mark)
        plot_data = [pdata[c] for c in cats]
        if any(len(d) < 2 for d in plot_data):
            ax.set_visible(False); continue

        # KW omnibus on RAW data
        kw_H, kw_p = kruskal(*[data[c] for c in cats])

        # violin + box + zero line (on PLOT data)
        ax.axhline(0, color="#888888", lw=0.8, linestyle="--", zorder=0)
        parts = ax.violinplot(plot_data, showmedians=False, showextrema=False,
                              widths=0.75)
        for pc, c in zip(parts["bodies"], cats):
            pc.set_facecolor(SWITCH_COLORS[c]); pc.set_alpha(0.75)
            pc.set_edgecolor("none")
        ax.boxplot(plot_data, widths=0.16, showfliers=False, patch_artist=True,
                   medianprops=dict(color="black", linewidth=1.6),
                   boxprops=dict(facecolor="white", alpha=0.9, linewidth=0.7),
                   whiskerprops=dict(linewidth=0.7), capprops=dict(linewidth=0.7))

        # SHARED axis (identical for every mark)
        ax.set_ylim(*SHARED_YLIM)

        # Dunn + Cliff's delta (bin level, RAW data)
        pair_rows = []
        if sp is not None:
            long_df = pd.DataFrame({
                "val": np.concatenate([data[c] for c in cats]),
                "cat": np.concatenate([[c] * len(data[c]) for c in cats])})
            dunn = sp.posthoc_dunn(long_df, val_col="val", group_col="cat",
                                   p_adjust="fdr_bh")
            for c1, c2 in combinations(cats, 2):
                p_pair = float(dunn.loc[c1, c2])
                d_pair, mag = cliffs_delta_fast(data[c1], data[c2])
                pair_rows.append({"mark": mark, "group1": c1, "group2": c2,
                                  "dunn_p_BH": p_pair,
                                  "cliffs_delta_g1_vs_g2": d_pair,
                                  "effect_magnitude": mag,
                                  "stars": _p_to_stars(p_pair)})
            all_stats.extend(pair_rows)

            # brackets anchored to the TOP of the shared window (closest pairs lowest)
            cat_pos = {c: i + 1 for i, c in enumerate(cats)}
            sig = list(pair_rows)   # show ALL pairs; non-significant print as "ns"
            sig.sort(key=lambda r: abs(cat_pos[r["group1"]] - cat_pos[r["group2"]]))
            if sig:
                # farthest-apart pairs on top so longer brackets don't cross shorter ones
                sig.sort(key=lambda r: abs(cat_pos[r["group1"]] - cat_pos[r["group2"]]))
                tick = span * 0.12                       # span from Pass-1
                base_y = yhi + 0.35 * span                # start just above the data
                for k, r in enumerate(sig):
                    label = (f"{r['stars']} \u03b4={r['cliffs_delta_g1_vs_g2']:+.2f} "
                             f"({r['effect_magnitude']})")
                    _bracket(ax, cat_pos[r["group1"]], cat_pos[r["group2"]],
                             base_y + k * tick, label)

        ax.set_xticks(range(1, len(cats) + 1))
        ax.set_xticklabels([SWITCH_LABELS[c] for c in cats], rotation=20, ha="right",
                           fontsize=9)
        ax.set_ylabel(f"{mark}\n{ylab}", labelpad=4)
        ax.set_title(f"{mark}  (KW H={kw_H:.0f}, "
                     + ("p<1e-300)" if kw_p == 0 or kw_p < 1e-300 else f"p={kw_p:.1e})"),
                     fontweight="bold", fontsize=10)

    for ax in axes[len(marks):]:
        ax.set_visible(False)

    fig.suptitle("Histone-mark CHANGE across compartment-switch categories\n"
                 "(bin-level; one value per 250kb bin — no pseudoreplication)",
                 fontweight="bold", fontsize=12, y=1.00)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(outdir, f"Step4a_markchange_by_switch.{ext}")
        fig.savefig(out, bbox_inches="tight", pad_inches=0.05, dpi=600)
        print(f"  Saved -> {out}")
    plt.close(fig)

    if all_stats:
        pd.DataFrame(all_stats).to_csv(
            os.path.join(outdir, "Step4a_markchange_stats.txt"),
            sep="\t", index=False)
        print(f"  Saved -> {os.path.join(outdir, 'Step4a_markchange_stats.txt')}")


# =============================================================================
# 4b — three-way integrative: per switching bin, mark-change vs mean expression
# =============================================================================
def step4b(bins, mark_cols, step3_path, outdir, coord_col, gene_name_col):
    """For SWITCHING bins, relate each mark's change to the mean expression
    change of the genes in that bin. Aggregating genes to a per-bin mean keeps
    the unit consistent with the bin-level mark value (no pseudoreplication)."""
    print("[INFO] Step 4b: integrating expression (per-bin mean) with mark change...")
    g = pd.read_csv(step3_path, sep="\t", low_memory=False)
    if coord_col not in g.columns or "lfc" not in g.columns:
        print(f"  [WARN] Step3 file lacks '{coord_col}' or 'lfc'; skipping 4b.")
        return

    # Mean expression change per bin (aggregate genes -> one value per bin).
    expr_by_bin = (g.groupby(coord_col)
                     .agg(mean_expr_lfc=("lfc", "mean"),
                          n_genes=("lfc", "size"))
                     .reset_index())

    # Join to bins (switch_cat + mark changes).
    merged = bins.merge(expr_by_bin, on=coord_col, how="inner")
    switch_bins = merged[merged["switch_cat"].isin(["A_to_B", "B_to_A"])].copy()
    if len(switch_bins) < 10:
        print(f"  [WARN] Only {len(switch_bins)} switching bins with expression; "
              f"skipping 4b.")
        return

    # Save the integrated per-bin table.
    keep = [coord_col, "switch_cat", "mean_expr_lfc", "n_genes"] + list(mark_cols.values())
    keep = [c for c in keep if c in merged.columns]
    merged[keep].to_csv(os.path.join(outdir, "Step4b_bin_integrated.txt"),
                        sep="\t", index=False)

    # Figure: one panel per mark, x=mark change, y=mean expression change,
    # colored by switch direction. One dot per SWITCHING bin.
    marks = list(mark_cols.keys())
    n = len(marks); n_cols = min(2, n); n_rows = (n + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
    axes = np.atleast_1d(axes).flatten()

    from scipy.stats import spearmanr
    for ax, mark in zip(axes, marks):
        col = mark_cols[mark]
        sub = switch_bins.dropna(subset=[col, "mean_expr_lfc"])
        if len(sub) < 10:
            ax.set_visible(False); continue
        for cat in ("A_to_B", "B_to_A"):
            d = sub[sub["switch_cat"] == cat]
            ax.scatter(d[col], d["mean_expr_lfc"], s=14, alpha=0.55,
                       color=SWITCH_COLORS[cat], edgecolors="none",
                       label=SWITCH_LABELS[cat])
        rho, p = spearmanr(sub[col], sub["mean_expr_lfc"])
        ax.axhline(0, color="#bbbbbb", lw=0.6); ax.axvline(0, color="#bbbbbb", lw=0.6)
        ax.set_xlabel(f"{mark} change  log2(B10/Parental)", labelpad=5)
        ax.set_ylabel("Mean expression change\n(bin mean log2FC)", labelpad=5)
        ax.set_title(f"{mark}  (Spearman rho={rho:+.2f}, "
                     + ("p<1e-300" if p == 0 else f"p={p:.1e})"),
                     fontweight="bold", fontsize=10)
        ax.legend(frameon=False, fontsize=8, loc="best")

    for ax in axes[len(marks):]:
        ax.set_visible(False)
    fig.suptitle("Mark change vs expression change in switching bins\n"
                 "(one dot per switching 250kb bin; association, not causation)",
                 fontweight="bold", fontsize=12, y=1.00)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(outdir, f"Step4b_mark_vs_expression.{ext}")
        fig.savefig(out, bbox_inches="tight", pad_inches=0.05, dpi=600)
        print(f"  Saved -> {out}")
    plt.close(fig)
    print(f"  Saved -> {os.path.join(outdir, 'Step4b_bin_integrated.txt')}")


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    cc = args.coord_col

    print("[INFO] Loading Step 1 (switch categories)...")
    step1 = pd.read_csv(args.step1, sep="\t", low_memory=False)
    for c in (cc, "switch_cat"):
        if c not in step1.columns:
            sys.exit(f"[ERROR] '{c}' not in {args.step1}.")

    print("[INFO] Loading 250kb CUT&RUN table...")
    cutrun = pd.read_csv(args.cutrun, sep="\t", low_memory=False)
    ck = auto_cutrun_coord(cutrun, args.cutrun_coord_col)
    print(f"[INFO] cutrun coordinate column: {ck}")

    mark_cols = detect_markchange_cols(cutrun.columns, args.marks)
    if not mark_cols:
        sys.exit("[ERROR] No migration mark-change columns found "
                 "(expected log2ratio_p-b10_<mark>_mean). "
                 f"Available: {[c for c in cutrun.columns if 'log2ratio' in c.lower()]}")
    print(f"[INFO] Mark-change columns detected: {mark_cols}")

    # Coerce mark columns numeric.
    for col in mark_cols.values():
        cutrun[col] = pd.to_numeric(cutrun[col], errors="coerce")

    # Join switch_cat onto the cutrun bins via coordinate.
    cutrun[ck] = cutrun[ck].astype(str).str.strip()
    step1[cc] = step1[cc].astype(str).str.strip()
    bins = cutrun.merge(step1[[cc, "switch_cat"]],
                        left_on=ck, right_on=cc, how="inner")
    if ck != cc and cc not in bins.columns:
        bins[cc] = bins[ck]
    bins = bins[bins["switch_cat"].isin(CAT_ORDER)].drop_duplicates(subset=cc)
    print(f"[INFO] Bins with switch category + mark change: {len(bins):,}")

    print("\n[INFO] === Step 4a: mark change by switch category (bin level) ===")
    step4a(bins, mark_cols, args.outdir, zscore=args.zscore)

    if args.step3:
        print("\n[INFO] === Step 4b: mark change vs expression (integrative) ===")
        step4b(bins, mark_cols, args.step3, args.outdir, cc, args.gene_name_col)
    else:
        print("\n[INFO] --step3 not provided; skipping 4b integrative view. "
              "Provide Step3_gene_expression_joined.txt to enable it.")

    print("\n[INFO] Step 4 complete. This is the integrative figure set: 4a shows "
          "mark dynamics differ by switch category (bin-level, no "
          "pseudoreplication); 4b shows mark-change tracks expression-change in "
          "switching bins. Frame as ASSOCIATION; causal arm is UNC.")


if __name__ == "__main__":
    main()