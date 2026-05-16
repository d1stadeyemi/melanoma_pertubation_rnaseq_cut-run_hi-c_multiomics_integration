"""
=============================================================================
H3K27me3 Disambiguation Analysis at B→A Compartment-Switching Loci
=============================================================================

PURPOSE
-------
Reconcile the apparent paradox between two findings:

  (A) Q1 ML result — Parental H3K27me3 predicts B→A compartment switching.
      Naive reading: H3K27me3 must be REMOVED at these loci during
      Parental → Bottom10 transition to permit activation.

  (B) Wet-lab UNC1999 (EZH1/2 inhibitor) result — UNC slows Bottom10
      migration. If the mark is "already gone" at B→A loci, what is UNC
      actually inhibiting?

This script runs three direct tests against the existing data tables to
disambiguate among four mechanistic possibilities:

  1. Mark removed at B→A loci  (canonical Polycomb derepression model)
  2. Mark retained but flux is ongoing (dynamic-equilibrium model)
  3. Mark redeposited at B→A loci (paradoxical / non-canonical)
  4. UNC acts at different loci than B→A switchers (off-target population)

INPUT FILES
-----------
  cut_run_signal_table.txt
      Must contain per-gene columns (case-sensitive as in the file header):
        gene_name
        parental_h3k27me3_mean
        b10_h3k27me3_mean
        parental_unc_h3k27me3_mean         (Parental + UNC)
        b10_unc_h3k27me3_mean              (Bottom10 + UNC)
        b10_minus_parental_h3k27me3        (Bottom10 - Parental delta)
        b10_minus_parental_unc_h3k27me3    (Bottom10+UNC - Parental+UNC delta)

  merged_hic_rnaseq_cutrun_parent_bottom10.txt
      Must contain per-gene columns:
        gene_name
        P->B10_Switch_Orig      (categorical: BtoA, AtoB, -)
        P->B10_Both_original    (categorical, more stringent)
      Optionally also:
        P->B10_Shift_Orig
        Parental_orig, Bottom10_orig    (raw PC1 values, used for Stable_A/B
                                          fallback classification)

OUTPUTS (in --outdir)
---------------------
  Tables (tab-separated):
    Disambiguation_summary.txt              Test-level summary statistics
    Test1_BtoA_H3K27me3_P_vs_B10.txt        Per-gene values for Test 1
    Test2_BtoA_H3K27me3_B10_vs_B10UNC.txt   Per-gene values for Test 2
    Test3_UNC_responsive_loci.txt           Loci where UNC changes H3K27me3
    Test3_overlap_with_BtoA.txt             Genes in BOTH UNC-responsive
                                            AND B→A sets (direct effect)
    Test3_UNC_only_not_BtoA.txt             UNC-responsive but NOT B→A
                                            (UNC acts elsewhere)

  Figures (pdf + png):
    Fig_Test1_BtoA_P_vs_B10.pdf             Violin: H3K27me3 at B→A loci
                                            in Parental vs Bottom10
    Fig_Test2_BtoA_B10_vs_B10UNC.pdf        Violin: H3K27me3 at B→A loci
                                            in Bottom10 vs Bottom10+UNC
    Fig_Test3_overlap_venn.pdf              Two-set Venn: UNC-responsive
                                            vs B→A switcher loci
    Fig_summary_4panel.pdf                  4-panel summary combining the
                                            three tests + mechanism table

USAGE
-----
  python scripts/h3k27me3_disambiguation.py \
      data/processed/cut_run_signal_table.txt \
      data/processed/merged_hic_rnaseq_cutrun_parent_bottom10.txt \
      --outdir results/exploratory/UNC_disambiguation \
      --unc_threshold 0.02

  Optional flags:
    --unc_threshold FLOAT
        |b10_unc_h3k27me3_mean - b10_h3k27me3_mean| cutoff for calling a
        locus "UNC-responsive" in Test 3. Default: 0.02 (CPM units).

    --include_aTob
        Include A→B silencing loci alongside B→A in an additional
        comparison panel. Default: B→A only.

NOTES ON STATISTICS
-------------------
- Within-locus (paired) comparisons use Wilcoxon signed-rank test, because
  each gene contributes paired observations across conditions (same locus,
  different treatment). Wilcoxon is non-parametric and handles the skewed
  CPM signal distribution that violates t-test assumptions.

- Effect size: Cliff's delta (Romano 2006 magnitude thresholds), the
  non-parametric companion the rest of the Q1 pipeline uses for
  consistency. Magnitudes: |δ| < 0.147 negligible, < 0.33 small,
  < 0.474 medium, else large.

- Between-locus (independent) comparisons (Test 3 UNC-responsive vs not)
  use Mann-Whitney U.

- Multiple-testing: only three primary tests are run, so BH-FDR is reported
  but the raw p-values are usable directly given the small comparison count.

=============================================================================
"""

import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import wilcoxon, mannwhitneyu
from statsmodels.stats.multitest import multipletests


# =============================================================================
# Publication aesthetics 
# =============================================================================
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

# Consistent colour palette 
COLORS = {
    "Parental":    "#85C1E9",   # light blue
    "Bottom10":    "#C0392B",   # dark red
    "Bottom10+UNC":"#7D3C98",   # purple — distinguishes UNC condition
    "Parental+UNC":"#5DADE2",   # mid blue
    "BtoA":        "#C0392B",
    "AtoB":        "#1A5276",
    "UNC_resp":    "#F39C12",   # amber for UNC-responsive loci
    "neutral":     "#95A5A6",
}


# =============================================================================
# CLI
# =============================================================================
def parse_args():
    p = argparse.ArgumentParser(
        description="H3K27me3 disambiguation analysis: Q1 ML vs UNC wet-lab"
    )
    p.add_argument("cutrun_file",
                   help="cut_run_signal_table.txt with parental_/b10_/_unc_ columns")
    p.add_argument("merged_file",
                   help="merged_hic_rnaseq_cutrun_parent_bottom10.txt with switch labels")
    p.add_argument("--outdir",
                   default="results/exploratory/UNC_disambiguation",
                   help="Output directory (default: %(default)s)")
    p.add_argument("--unc_threshold", type=float, default=0.05,
                   help="|B10+UNC - B10| H3K27me3 delta to call locus "
                        "UNC-responsive. Default 0.05 CPM.")
    p.add_argument("--include_aTob", action="store_true",
                   help="Also generate a complementary panel for A→B silencing loci.")
    return p.parse_args()


# =============================================================================
# Effect size helper — same Cliff's delta used elsewhere in the pipeline
# =============================================================================
def cliffs_delta_fast(x, y):
    """
    Cliff's delta = P(x > y) - P(x < y), in [-1, 1].
    Memory-efficient via sort + searchsorted.
    Returns (delta, magnitude_string).
    """
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return np.nan, "NA"
    ys = np.sort(y)
    lt = np.searchsorted(ys, x, side="left")     # count of y<x per element
    le = np.searchsorted(ys, x, side="right")    # count of y<=x per element
    gt = ny - le                                  # count of y>x per element
    d  = (int(lt.sum()) - int(gt.sum())) / (nx * ny)
    # Sign convention: positive => x tends to exceed y
    # Above formula: lt - gt over nx*ny = (P(x>y) - P(x<y))
    ad = abs(d)
    if   ad < 0.147: mag = "negligible"
    elif ad < 0.330: mag = "small"
    elif ad < 0.474: mag = "medium"
    else:            mag = "large"
    return d, mag


# =============================================================================
# Switch classification — same priority cascade as RF_SHAP.py for consistency
# =============================================================================
def classify_switch(merged_df):
    """
    Assign each gene to B_to_A / A_to_B / Stable_A / Stable_B / None
    using the priority cascade: Both > Switch > Shift > infer from PC1.
    Returns a list aligned to merged_df rows.
    """
    cats = []
    for _, row in merged_df.iterrows():
        cat = None
        # 1. Most stringent: Both (switch + shift agree)
        v = str(row.get("P->B10_Both_original", "")).strip()
        if v not in ["-", "", "nan", "None"]:
            if   "toA" in v or "BtoA" in v: cat = "B_to_A"
            elif "toB" in v or "AtoB" in v: cat = "A_to_B"
        # 2. Sign switch
        if cat is None:
            v = str(row.get("P->B10_Switch_Orig", "")).strip()
            if v not in ["-", "", "nan", "None"]:
                if   v == "BtoA": cat = "B_to_A"
                elif v == "AtoB": cat = "A_to_B"
        # 3. Magnitude shift (no sign change but moves toward A/B)
        if cat is None:
            v = str(row.get("P->B10_Shift_Orig", "")).strip()
            if v not in ["-", "", "nan", "None"]:
                if   "TowardA" in v: cat = "B_to_A"
                elif "TowardB" in v: cat = "A_to_B"
        # 4. Fallback: stable categories from raw PC1 sign
        if cat is None:
            try:
                p_orig = float(row.get("Parental_orig", np.nan))
                b_orig = float(row.get("Bottom10_orig", np.nan))
                if   p_orig > 0 and b_orig > 0: cat = "Stable_A"
                elif p_orig < 0 and b_orig < 0: cat = "Stable_B"
            except (ValueError, TypeError):
                pass
        cats.append(cat)
    return cats


# =============================================================================
# Load + merge
# =============================================================================
def load_and_merge(cutrun_file, merged_file):
    """
    Load both tables, classify compartment switches, and inner-join on
    gene_name. Returns a DataFrame with:
      - gene_name
      - switch_cat
      - parental_h3k27me3_mean
      - b10_h3k27me3_mean
      - parental_unc_h3k27me3_mean
      - b10_unc_h3k27me3_mean
      - b10_minus_parental_h3k27me3
      - b10_minus_parental_unc_h3k27me3
    """
    print(f"[INFO] Loading cut-and-run table: {cutrun_file}")
    cr = pd.read_csv(cutrun_file, sep="\t", low_memory=False)
    print(f"  {len(cr):,} rows x {len(cr.columns)} columns")

    print(f"[INFO] Loading merged switch table: {merged_file}")
    mg = pd.read_csv(merged_file, sep="\t", low_memory=False)
    print(f"  {len(mg):,} rows x {len(mg.columns)} columns")

    required_cr = [
        "gene_name",
        "parental_h3k27me3_mean",
        "b10_h3k27me3_mean",
        "parental_unc_h3k27me3_mean",
        "b10_unc_h3k27me3_mean",
        "b10_minus_parental_h3k27me3",
        "b10_minus_parental_unc_h3k27me3",
    ]
    missing = [c for c in required_cr if c not in cr.columns]
    if missing:
        sys.exit(f"[ERROR] Missing columns in cut_run table: {missing}")

    cr_sub = cr[required_cr].drop_duplicates(subset="gene_name").copy()
    for c in required_cr[1:]:
        cr_sub[c] = pd.to_numeric(cr_sub[c], errors="coerce")

    # Classify switches on the merged file
    print("[INFO] Classifying compartment switches...")
    switch_cols_present = [
        c for c in ["gene_name", "Parental_orig", "Bottom10_orig",
                    "P->B10_Switch_Orig", "P->B10_Shift_Orig",
                    "P->B10_Both_original"]
        if c in mg.columns
    ]
    mg_sw = mg[switch_cols_present].drop_duplicates(subset="gene_name").copy()
    mg_sw["switch_cat"] = classify_switch(mg_sw)
    print("[INFO] Switch category counts:")
    print(mg_sw["switch_cat"].value_counts(dropna=False).to_string())

    df = cr_sub.merge(
        mg_sw[["gene_name", "switch_cat"]],
        on="gene_name", how="inner"
    )
    df = df[df["switch_cat"].notna()].copy()

    # Drop rows missing any of the four signal columns we need
    sig_cols = [
        "parental_h3k27me3_mean", "b10_h3k27me3_mean",
        "parental_unc_h3k27me3_mean", "b10_unc_h3k27me3_mean",
    ]
    df = df.dropna(subset=sig_cols)
    print(f"[INFO] After merge + drop-NA: {len(df):,} genes")
    print(df["switch_cat"].value_counts().to_string())
    return df


# =============================================================================
# Test 1 — H3K27me3 at B→A loci: Parental vs Bottom10
# =============================================================================
def test1_parental_vs_b10(df, outdir, switch_label="B_to_A"):
    """
    Paired Wilcoxon test comparing parental_h3k27me3_mean vs
    b10_h3k27me3_mean at B→A loci.

    Null hypothesis: H3K27me3 at B→A loci is unchanged Parental → Bottom10.
    Rejecting null with negative direction (B10 < Parental) supports the
    "mark removal" model.
    """
    sub = df[df["switch_cat"] == switch_label].copy()
    if len(sub) < 5:
        print(f"[WARN] Test 1: only {len(sub)} {switch_label} loci, skipping.")
        return None

    p_vals  = sub["parental_h3k27me3_mean"].values
    b_vals  = sub["b10_h3k27me3_mean"].values
    paired_diff = b_vals - p_vals    # negative => mark removed in B10

    # Wilcoxon signed-rank (paired, non-parametric)
    try:
        w_stat, w_p = wilcoxon(b_vals, p_vals, alternative="two-sided")
    except ValueError as e:
        # Happens if all pairs identical — extreme edge case
        print(f"[WARN] Wilcoxon failed: {e}")
        w_stat, w_p = np.nan, np.nan

    # Cliff's delta on the same paired vectors treated as independent samples
    # (Romano 2006 magnitude thresholds; direction = b10 vs parental)
    d, mag = cliffs_delta_fast(b_vals, p_vals)

    median_p = float(np.median(p_vals))
    median_b = float(np.median(b_vals))
    median_diff = float(np.median(paired_diff))

    summary = {
        "test":               "Test1_BtoA_Parental_vs_Bottom10",
        "n_loci":             len(sub),
        "median_parental":    median_p,
        "median_b10":         median_b,
        "median_paired_diff": median_diff,
        "wilcoxon_W":         float(w_stat) if w_stat == w_stat else np.nan,
        "wilcoxon_p":         float(w_p)    if w_p    == w_p    else np.nan,
        "cliffs_delta":       d,
        "effect_magnitude":   mag,
        "interpretation":     (
            "B10 < Parental (mark removed)" if median_diff < 0
            else "B10 > Parental (mark gained)" if median_diff > 0
            else "B10 ≈ Parental (mark retained)"
        ),
    }

    # Export per-gene table
    out_tbl = os.path.join(outdir, "Test1_BtoA_H3K27me3_P_vs_B10.txt")
    sub_out = sub[["gene_name", "switch_cat",
                   "parental_h3k27me3_mean", "b10_h3k27me3_mean",
                   "b10_minus_parental_h3k27me3"]].copy()
    sub_out["b10_minus_parental_calc"] = b_vals - p_vals
    sub_out.to_csv(out_tbl, sep="\t", index=False)
    print(f"[INFO] Test 1 per-gene table -> {out_tbl}")

    # Figure
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    parts = ax.violinplot(
        [np.log1p(p_vals), np.log1p(b_vals)],
        showmedians=False, showextrema=False, widths=0.7
    )
    for pc, c in zip(parts["bodies"], [COLORS["Parental"], COLORS["Bottom10"]]):
        pc.set_facecolor(c)
        pc.set_alpha(0.75)
        pc.set_edgecolor("none")
    # Colored median crossbars
    for i, (vals, c) in enumerate(
            zip([np.log1p(p_vals), np.log1p(b_vals)],
                [COLORS["Parental"], COLORS["Bottom10"]]), start=1):
        med = np.median(vals)
        ax.hlines(med, i - 0.22, i + 0.22, colors="black",
                  linewidth=2.2, zorder=3)
        ax.hlines(med, i - 0.18, i + 0.18, colors=c, linewidth=1.4, zorder=4)

    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Parental", "Bottom10"])
    ax.set_ylabel("log1p(H3K27me3 CPM signal)")
    ax.set_title(f"H3K27me3 at {switch_label.replace('_', '→')} loci\n"
                 f"Parental vs Bottom10",
                 fontweight="bold", pad=10)

    p_str = ("p < 1e-300" if (w_p == 0 or w_p < 1e-300)
             else f"p = {w_p:.1e}")
    ax.text(
        0.97, 0.97,
        f"n = {len(sub):,} loci\nMedian Δ = {median_diff:+.4f}\n"
        f"Cliff's δ = {d:+.3f} ({mag})\n{p_str}",
        transform=ax.transAxes, ha="right", va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.30",
                  facecolor="white", edgecolor="lightgrey",
                  linewidth=0.5, alpha=0.9)
    )
    plt.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(outdir, f"Fig_Test1_BtoA_P_vs_B10.{ext}")
        fig.savefig(out, bbox_inches="tight", pad_inches=0.05, dpi=600)
    plt.close(fig)
    print(f"[INFO] Test 1 figure -> Fig_Test1_BtoA_P_vs_B10.pdf/png")

    return summary


# =============================================================================
# Test 2 — H3K27me3 at B→A loci: Bottom10 vs Bottom10+UNC
# =============================================================================
def test2_b10_vs_b10unc(df, outdir, switch_label="B_to_A"):
    """
    Paired Wilcoxon test comparing b10_h3k27me3_mean vs
    b10_unc_h3k27me3_mean at B→A loci.

    If UNC's effect at these loci is meaningful for the migration phenotype,
    UNC should further reduce H3K27me3 here. If UNC acts elsewhere, signal
    should be unchanged.
    """
    sub = df[df["switch_cat"] == switch_label].copy()
    if len(sub) < 5:
        print(f"[WARN] Test 2: only {len(sub)} {switch_label} loci, skipping.")
        return None

    b_vals    = sub["b10_h3k27me3_mean"].values
    bunc_vals = sub["b10_unc_h3k27me3_mean"].values
    paired_diff = bunc_vals - b_vals

    try:
        w_stat, w_p = wilcoxon(bunc_vals, b_vals, alternative="two-sided")
    except ValueError as e:
        print(f"[WARN] Wilcoxon failed: {e}")
        w_stat, w_p = np.nan, np.nan

    d, mag = cliffs_delta_fast(bunc_vals, b_vals)

    median_b    = float(np.median(b_vals))
    median_bunc = float(np.median(bunc_vals))
    median_diff = float(np.median(paired_diff))

    summary = {
        "test":               "Test2_BtoA_Bottom10_vs_Bottom10UNC",
        "n_loci":             len(sub),
        "median_b10":         median_b,
        "median_b10_unc":     median_bunc,
        "median_paired_diff": median_diff,
        "wilcoxon_W":         float(w_stat) if w_stat == w_stat else np.nan,
        "wilcoxon_p":         float(w_p)    if w_p    == w_p    else np.nan,
        "cliffs_delta":       d,
        "effect_magnitude":   mag,
        "interpretation":     (
            "B10+UNC < B10 (UNC reduces mark further at B→A loci → direct effect)"
            if median_diff < 0 and abs(d) >= 0.147 else
            "B10+UNC ≈ B10 (UNC does not act at B→A loci → indirect mechanism)"
            if abs(d) < 0.147 else
            "B10+UNC > B10 (paradoxical: UNC increases mark — needs investigation)"
        ),
    }

    out_tbl = os.path.join(outdir, "Test2_BtoA_H3K27me3_B10_vs_B10UNC.txt")
    sub_out = sub[["gene_name", "switch_cat",
                   "b10_h3k27me3_mean", "b10_unc_h3k27me3_mean",
                   "b10_minus_parental_unc_h3k27me3"]].copy()
    sub_out["b10unc_minus_b10"] = bunc_vals - b_vals
    sub_out.to_csv(out_tbl, sep="\t", index=False)
    print(f"[INFO] Test 2 per-gene table -> {out_tbl}")

    # Figure
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    parts = ax.violinplot(
        [np.log1p(b_vals), np.log1p(bunc_vals)],
        showmedians=False, showextrema=False, widths=0.7
    )
    for pc, c in zip(parts["bodies"], [COLORS["Bottom10"], COLORS["Bottom10+UNC"]]):
        pc.set_facecolor(c)
        pc.set_alpha(0.75)
        pc.set_edgecolor("none")
    for i, (vals, c) in enumerate(
            zip([np.log1p(b_vals), np.log1p(bunc_vals)],
                [COLORS["Bottom10"], COLORS["Bottom10+UNC"]]), start=1):
        med = np.median(vals)
        ax.hlines(med, i - 0.22, i + 0.22, colors="black",
                  linewidth=2.2, zorder=3)
        ax.hlines(med, i - 0.18, i + 0.18, colors=c, linewidth=1.4, zorder=4)

    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Bottom10", "Bottom10+UNC"])
    ax.set_ylabel("log1p(H3K27me3 CPM signal)")
    ax.set_title(f"H3K27me3 at {switch_label.replace('_', '→')} loci\n"
                 f"Bottom10 vs Bottom10+UNC",
                 fontweight="bold", pad=10)

    p_str = ("p < 1e-300" if (w_p == 0 or w_p < 1e-300)
             else f"p = {w_p:.1e}")
    ax.text(
        0.97, 0.97,
        f"n = {len(sub):,} loci\nMedian Δ = {median_diff:+.4f}\n"
        f"Cliff's δ = {d:+.3f} ({mag})\n{p_str}",
        transform=ax.transAxes, ha="right", va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.30",
                  facecolor="white", edgecolor="lightgrey",
                  linewidth=0.5, alpha=0.9)
    )
    plt.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(outdir, f"Fig_Test2_BtoA_B10_vs_B10UNC.{ext}")
        fig.savefig(out, bbox_inches="tight", pad_inches=0.05, dpi=600)
    plt.close(fig)
    print(f"[INFO] Test 2 figure -> Fig_Test2_BtoA_B10_vs_B10UNC.pdf/png")

    return summary


# =============================================================================
# Test 3 — Differential analysis: UNC-responsive loci vs B→A overlap
# =============================================================================
def test3_overlap(df, outdir, unc_threshold=0.02):
    """
    Identify loci where UNC treatment changes H3K27me3 above threshold,
    then quantify overlap with the B→A switcher set.

    UNC-responsive locus  := |b10_unc - b10| H3K27me3 > unc_threshold
    B→A locus             := switch_cat == "B_to_A"

    High overlap => UNC acts AT the switched loci (direct effect).
    Low overlap  => UNC acts at DIFFERENT loci (indirect mechanism).
    """
    df_t = df.copy()
    df_t["b10unc_minus_b10"] = (df_t["b10_unc_h3k27me3_mean"]
                                - df_t["b10_h3k27me3_mean"])
    df_t["abs_unc_delta"]    = df_t["b10unc_minus_b10"].abs()
    df_t["unc_responsive"]   = df_t["abs_unc_delta"] > unc_threshold
    df_t["is_BtoA"]          = df_t["switch_cat"] == "B_to_A"

    # Counts for 2x2 contingency
    n_total          = len(df_t)
    n_unc_resp       = int(df_t["unc_responsive"].sum())
    n_BtoA           = int(df_t["is_BtoA"].sum())
    n_both           = int((df_t["unc_responsive"] & df_t["is_BtoA"]).sum())
    n_BtoA_not_UNC   = n_BtoA - n_both
    n_UNC_not_BtoA   = n_unc_resp - n_both
    n_neither        = n_total - n_unc_resp - n_BtoA + n_both

    # Fisher's exact test for independence of the two sets
    from scipy.stats import fisher_exact
    table = np.array([[n_both,         n_BtoA_not_UNC],
                      [n_UNC_not_BtoA, n_neither]])
    try:
        odds_ratio, fisher_p = fisher_exact(table, alternative="two-sided")
    except Exception:
        odds_ratio, fisher_p = np.nan, np.nan

    # Enrichment statistics
    p_overlap_obs = n_both / n_BtoA if n_BtoA > 0 else np.nan
    p_overlap_exp = n_unc_resp / n_total if n_total > 0 else np.nan
    enrichment    = (p_overlap_obs / p_overlap_exp
                     if (p_overlap_exp and p_overlap_exp > 0) else np.nan)

    summary = {
        "test":              "Test3_UNC_responsive_vs_BtoA_overlap",
        "unc_threshold":     unc_threshold,
        "n_total":           n_total,
        "n_unc_responsive":  n_unc_resp,
        "n_BtoA":            n_BtoA,
        "n_overlap":         n_both,
        "frac_BtoA_in_UNC":  round(p_overlap_obs, 4) if p_overlap_obs == p_overlap_obs else np.nan,
        "frac_UNC_in_BtoA":  round(n_both / n_unc_resp, 4) if n_unc_resp > 0 else np.nan,
        "expected_overlap_frac": round(p_overlap_exp, 4) if p_overlap_exp == p_overlap_exp else np.nan,
        "fold_enrichment":   round(enrichment, 3) if enrichment == enrichment else np.nan,
        "fisher_odds_ratio": float(odds_ratio) if odds_ratio == odds_ratio else np.nan,
        "fisher_p":          float(fisher_p)    if fisher_p   == fisher_p    else np.nan,
        "interpretation":    (
            "UNC acts AT B→A loci (direct mechanism)" if (enrichment is not np.nan
                and enrichment > 1.5 and fisher_p < 0.05) else
            "UNC acts at DIFFERENT loci (indirect mechanism)" if (enrichment is not np.nan
                and enrichment < 1.0) else
            "Ambiguous — overlap roughly matches expectation"
        ),
    }

    # Export per-set tables
    out_unc = os.path.join(outdir, "Test3_UNC_responsive_loci.txt")
    df_t[df_t["unc_responsive"]][
        ["gene_name", "switch_cat",
         "b10_h3k27me3_mean", "b10_unc_h3k27me3_mean",
         "b10unc_minus_b10", "is_BtoA"]
    ].to_csv(out_unc, sep="\t", index=False)
    print(f"[INFO] Test 3: UNC-responsive loci ({n_unc_resp:,}) -> {out_unc}")

    out_overlap = os.path.join(outdir, "Test3_overlap_with_BtoA.txt")
    df_t[df_t["unc_responsive"] & df_t["is_BtoA"]][
        ["gene_name", "switch_cat",
         "parental_h3k27me3_mean", "b10_h3k27me3_mean",
         "b10_unc_h3k27me3_mean", "b10unc_minus_b10"]
    ].to_csv(out_overlap, sep="\t", index=False)
    print(f"[INFO] Test 3: UNC-responsive AND B→A ({n_both:,}) -> {out_overlap}")

    out_unc_only = os.path.join(outdir, "Test3_UNC_only_not_BtoA.txt")
    df_t[df_t["unc_responsive"] & (~df_t["is_BtoA"])][
        ["gene_name", "switch_cat",
         "parental_h3k27me3_mean", "b10_h3k27me3_mean",
         "b10_unc_h3k27me3_mean", "b10unc_minus_b10"]
    ].to_csv(out_unc_only, sep="\t", index=False)
    print(f"[INFO] Test 3: UNC-responsive but NOT B→A ({n_UNC_not_BtoA:,}) "
          f"-> {out_unc_only}")

    # Two-set Venn diagram
    fig, ax = plt.subplots(figsize=(6.5, 5))
    try:
        from matplotlib_venn import venn2
        v = venn2(
            subsets=(n_BtoA_not_UNC, n_UNC_not_BtoA, n_both),
            set_labels=("B→A loci", "UNC-responsive loci"),
            ax=ax,
        )
        # Restyle if patches exist
        if v.get_patch_by_id("10"): v.get_patch_by_id("10").set_color(COLORS["BtoA"])
        if v.get_patch_by_id("01"): v.get_patch_by_id("01").set_color(COLORS["UNC_resp"])
        if v.get_patch_by_id("11"): v.get_patch_by_id("11").set_color("#8E44AD")
        for pid in ("10", "01", "11"):
            patch = v.get_patch_by_id(pid)
            if patch:
                patch.set_alpha(0.55)
                patch.set_edgecolor("white")
    except ImportError:
        # matplotlib_venn not installed → fall back to manual two-bar plot
        ax.bar(["B→A only", "Overlap", "UNC-resp only"],
               [n_BtoA_not_UNC, n_both, n_UNC_not_BtoA],
               color=[COLORS["BtoA"], "#8E44AD", COLORS["UNC_resp"]],
               edgecolor="white", linewidth=0.5)
        for i, v_count in enumerate([n_BtoA_not_UNC, n_both, n_UNC_not_BtoA]):
            ax.text(i, v_count + max([n_BtoA_not_UNC, n_both,
                                       n_UNC_not_BtoA]) * 0.02,
                    f"{v_count:,}", ha="center", va="bottom", fontsize=10)
        ax.set_ylabel("Number of genes")

    fold_str = (f"Fold enrichment = {enrichment:.2f}×"
                if enrichment == enrichment else "Fold enrichment = NA")
    p_str    = (f"Fisher p = {fisher_p:.2e}"
                if fisher_p == fisher_p else "Fisher p = NA")
    ax.set_title(
        f"UNC-responsive vs B→A Switcher Loci\n"
        f"{fold_str}, {p_str}",
        fontweight="bold", pad=10
    )
    fig.text(
        0.5, 0.02,
        f"UNC threshold: |Δ H3K27me3| > {unc_threshold}    "
        f"Total genes analyzed: {n_total:,}",
        ha="center", fontsize=9, style="italic", color="#555555"
    )
    plt.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(outdir, f"Fig_Test3_overlap_venn.{ext}")
        fig.savefig(out, bbox_inches="tight", pad_inches=0.05, dpi=600)
    plt.close(fig)
    print(f"[INFO] Test 3 figure -> Fig_Test3_overlap_venn.pdf/png")

    return summary


# =============================================================================
# Summary 4-panel figure: combine the three tests into a single story
# =============================================================================
def build_summary_figure(df, outdir, summaries, unc_threshold):
    """
    A single 2x2 figure summarising the three tests + mechanism interpretation.
    Useful for the lab meeting / one-slide briefing.
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    sub = df[df["switch_cat"] == "B_to_A"].copy()

    # Panel A — Test 1
    ax = axes[0, 0]
    p_vals = np.log1p(sub["parental_h3k27me3_mean"].values)
    b_vals = np.log1p(sub["b10_h3k27me3_mean"].values)
    parts = ax.violinplot([p_vals, b_vals], showmedians=False,
                          showextrema=False, widths=0.7)
    for pc, c in zip(parts["bodies"], [COLORS["Parental"], COLORS["Bottom10"]]):
        pc.set_facecolor(c); pc.set_alpha(0.75); pc.set_edgecolor("none")
    ax.set_xticks([1, 2]); ax.set_xticklabels(["Parental", "Bottom10"])
    ax.set_ylabel("log1p(H3K27me3 CPM)")
    ax.set_title("A. Parental vs Bottom10 at B→A loci",
                 fontweight="bold", loc="left")

    # Panel B — Test 2
    ax = axes[0, 1]
    b_vals    = np.log1p(sub["b10_h3k27me3_mean"].values)
    bunc_vals = np.log1p(sub["b10_unc_h3k27me3_mean"].values)
    parts = ax.violinplot([b_vals, bunc_vals], showmedians=False,
                          showextrema=False, widths=0.7)
    for pc, c in zip(parts["bodies"], [COLORS["Bottom10"], COLORS["Bottom10+UNC"]]):
        pc.set_facecolor(c); pc.set_alpha(0.75); pc.set_edgecolor("none")
    ax.set_xticks([1, 2]); ax.set_xticklabels(["Bottom10", "Bottom10+UNC"])
    ax.set_ylabel("log1p(H3K27me3 CPM)")
    ax.set_title("B. Bottom10 vs Bottom10+UNC at B→A loci",
                 fontweight="bold", loc="left")

    # Panel C — Test 3 (bar of three categories)
    ax = axes[1, 0]
    s3 = summaries.get("test3", {})
    counts = [
        s3.get("n_BtoA", 0) - s3.get("n_overlap", 0),
        s3.get("n_overlap", 0),
        s3.get("n_unc_responsive", 0) - s3.get("n_overlap", 0),
    ]
    bars = ax.bar(["B→A only", "Overlap", "UNC-resp only"], counts,
                  color=[COLORS["BtoA"], "#8E44AD", COLORS["UNC_resp"]],
                  edgecolor="white", linewidth=0.5)
    for bar, c in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(counts) * 0.02 if max(counts) else 1,
                f"{c:,}", ha="center", va="bottom", fontsize=10)
    ax.set_ylabel("Number of genes")
    ax.set_title("C. Set overlap: UNC-responsive vs B→A",
                 fontweight="bold", loc="left")

    # Panel D — Mechanism interpretation table
    ax = axes[1, 1]
    ax.axis("off")
    s1 = summaries.get("test1", {})
    s2 = summaries.get("test2", {})
    s3 = summaries.get("test3", {})
    interp_text = (
        f"D. Mechanism summary\n"
        f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
        f"Test 1 (P→B10):   {s1.get('interpretation', 'NA')}\n"
        f"  Cliff's δ = {s1.get('cliffs_delta', 'NA')} "
        f"({s1.get('effect_magnitude', 'NA')})\n\n"
        f"Test 2 (B10→B10+UNC): {s2.get('interpretation', 'NA')}\n"
        f"  Cliff's δ = {s2.get('cliffs_delta', 'NA')} "
        f"({s2.get('effect_magnitude', 'NA')})\n\n"
        f"Test 3 (overlap): {s3.get('interpretation', 'NA')}\n"
        f"  Fold enrichment = {s3.get('fold_enrichment', 'NA')}×, "
        f"Fisher p = {s3.get('fisher_p', 'NA')}\n"
        f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
        f"UNC threshold: |Δ H3K27me3| > {unc_threshold} CPM"
    )
    ax.text(0.02, 0.95, interp_text, transform=ax.transAxes,
            ha="left", va="top", fontsize=9.5, family="monospace",
            bbox=dict(boxstyle="round,pad=0.6", facecolor="#FAFAFA",
                      edgecolor="lightgrey", linewidth=0.6))

    fig.suptitle(
        "H3K27me3 Disambiguation at B→A Compartment-Switching Loci\n"
        "Reconciling Q1 ML prediction with UNC1999 wet-lab phenotype",
        fontweight="bold", fontsize=13, y=1.00
    )
    plt.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(outdir, f"Fig_summary_4panel.{ext}")
        fig.savefig(out, bbox_inches="tight", pad_inches=0.05, dpi=600)
    plt.close(fig)
    print(f"[INFO] Summary figure -> Fig_summary_4panel.pdf/png")


# =============================================================================
# Main
# =============================================================================
def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = load_and_merge(args.cutrun_file, args.merged_file)

    summaries = {}
    print("\n" + "=" * 70)
    print("Test 1 — H3K27me3 at B→A loci: Parental vs Bottom10")
    print("=" * 70)
    s1 = test1_parental_vs_b10(df, args.outdir, switch_label="B_to_A")
    if s1: summaries["test1"] = s1

    print("\n" + "=" * 70)
    print("Test 2 — H3K27me3 at B→A loci: Bottom10 vs Bottom10+UNC")
    print("=" * 70)
    s2 = test2_b10_vs_b10unc(df, args.outdir, switch_label="B_to_A")
    if s2: summaries["test2"] = s2

    print("\n" + "=" * 70)
    print("Test 3 — UNC-responsive loci vs B→A overlap")
    print("=" * 70)
    s3 = test3_overlap(df, args.outdir, unc_threshold=args.unc_threshold)
    if s3: summaries["test3"] = s3

    # Optional A→B complement
    if args.include_aTob:
        print("\n" + "=" * 70)
        print("Bonus — Repeating Tests 1 & 2 for A→B silencing loci")
        print("=" * 70)
        s1b = test1_parental_vs_b10(df, args.outdir, switch_label="A_to_B")
        s2b = test2_b10_vs_b10unc(df, args.outdir, switch_label="A_to_B")
        if s1b: summaries["test1_AtoB"] = s1b
        if s2b: summaries["test2_AtoB"] = s2b

    # Save text summary
    print("\n" + "=" * 70)
    print("Writing summary table & figure")
    print("=" * 70)
    summary_df = pd.DataFrame(list(summaries.values()))
    # BH correction across primary p-value columns
    p_cols = [c for c in ["wilcoxon_p", "fisher_p"] if c in summary_df.columns]
    for pc in p_cols:
        valid = summary_df[pc].notna()
        if valid.sum() > 1:
            adj = multipletests(summary_df.loc[valid, pc].values,
                                method="fdr_bh")[1]
            summary_df.loc[valid, pc + "_BH"] = adj

    out_summary = os.path.join(args.outdir, "Disambiguation_summary.txt")
    summary_df.to_csv(out_summary, sep="\t", index=False)
    print(f"[INFO] Summary table -> {out_summary}")
    print(summary_df.to_string(index=False))

    build_summary_figure(df, args.outdir, summaries, args.unc_threshold)

    print(f"\n[INFO] Complete. All outputs in: {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()