"""
=============================================================================
csaw_differential_summary.py
  One summary figure for the csaw differential-binding results across marks.

WHAT IT SHOWS
  A diverging horizontal bar chart: one row per mark (per contrast), with the
  number of significant differential regions GAINED in B10 (group1) drawn to
  the right and LOST in B10 drawn to the left of zero. Bar length = magnitude
  of differential response; the left/right split = directional bias. Each bar
  is annotated with totals and the % gained — so "which marks changed most?"
  and "in which direction?" are answered in one glance.

WHY RECOMPUTE (not hand-copy from stdout)
  The "N significant" and "% up" numbers csaw prints are just filters/counts on
  the FDR and rep.logFC columns of each *_all_regions.csv. Hand-copying breaks
  the audit trail and goes stale on re-run. This script GLOBS the all_regions
  CSVs and recomputes everything with the SAME thresholds as csaw.R
  (FDR<0.05, |log2FC|>=1, both exposed as flags), so the figure is reproducible
  and always matches your significance definition.

INTERPRETATION CAVEAT (state on the slide)
  Under TMM/library-size normalization, the direction balance (e.g. "76% gained
  in B10") is a RELATIVE-REDISTRIBUTION statement: among regions that changed,
  more gained than lost. It is NOT a genome-wide-increase claim (TMM forces the
  global average to ~0). A true global change requires spike-in (Orlando 2014).
  Phrase as "of differential regions, X% gained signal in B10", not "mark
  increased genome-wide". This is most defensible for the migration axis
  (no drug); flag the drug/UNC axis separately.

INPUT
  A directory of csaw '*_all_regions.csv' files (one per mark/contrast), each
  with seqnames/start/end + a logFC column (rep.logFC) + FDR.

USAGE
  python csaw_differential_summary.py \
      --csaw_dir results/exploratory/diff_analysis/csaw \
      --outdir   results/exploratory/diff_analysis/csaw/summary

  # restrict to a contrast, or change thresholds:
  python csaw_differential_summary.py --csaw_dir <dir> \
      --pattern "*B10-none_vs_Parental-none*all_regions.csv" \
      --fdr 0.05 --lfc 1.0

OUTPUTS
  csaw_differential_summary.pdf/png   diverging bar chart
  csaw_differential_summary.txt       tidy table: mark, contrast, n_tested,
                                       n_sig, n_up, n_down, pct_up, pct_down
"""

import argparse
import glob
import os
import re
import sys

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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

# Gained (up in B10/group1) = warm; Lost (down) = cool. Matches the study's
# B->A (red/activated) vs A->B (blue/silenced) intuition loosely, but here it
# is simply gain (right, warm) vs loss (left, cool).
COLOR_UP   = "#C0392B"   # gained in group1 (e.g. B10)
COLOR_DOWN = "#1A5276"   # lost in group1


def merged_coverage_bp(df):
    """Total genomic footprint (bp) covered by the tested/enriched regions,
    after merging any overlapping or book-ended intervals so genome is NOT
    double-counted. csaw windows are generated with overlapping spacing, so a
    naive sum(width) would overstate coverage; we collapse per chromosome first.
    Expects columns seqnames/start/end."""
    need = {"seqnames", "start", "end"}
    if not need.issubset(df.columns):
        return np.nan
    sub = df[["seqnames", "start", "end"]].dropna().copy()
    sub["start"] = pd.to_numeric(sub["start"], errors="coerce")
    sub["end"] = pd.to_numeric(sub["end"], errors="coerce")
    sub = sub.dropna()
    total = 0
    for _, chunk in sub.groupby("seqnames"):
        iv = chunk[["start", "end"]].sort_values("start").to_numpy()
        if len(iv) == 0:
            continue
        cur_s, cur_e = iv[0]
        for s, e in iv[1:]:
            if s <= cur_e:                 # overlap or adjacency -> extend
                cur_e = max(cur_e, e)
            else:                          # gap -> close current interval
                total += cur_e - cur_s
                cur_s, cur_e = s, e
        total += cur_e - cur_s
    return int(total)


def parse_args():
    p = argparse.ArgumentParser(description="csaw differential summary figure.")
    p.add_argument("--csaw_dir", required=True,
                   help="Directory containing csaw *_all_regions.csv files.")
    p.add_argument("--pattern", default="*_all_regions.csv",
                   help="Glob pattern for the region CSVs [default %default].")
    p.add_argument("--outdir", default=None,
                   help="Output dir [default: <csaw_dir>/summary].")
    p.add_argument("--fdr", type=float, default=0.05,
                   help="FDR threshold (match csaw.R) [default %default].")
    p.add_argument("--lfc", type=float, default=1.0,
                   help="Min |log2FC| (match csaw.R) [default %default].")
    p.add_argument("--lfc_col", default="rep.logFC",
                   help="logFC column in the CSVs [default %default].")
    p.add_argument("--fdr_col", default="FDR",
                   help="FDR column in the CSVs [default %default].")
    p.add_argument("--group1_label", default="B10",
                   help="Name of group1/numerator for axis labels [default %default].")
    p.add_argument("--genome_size", type=float, default=2.88e9,
                   help="Denominator for %% genome covered. Default 2.88e9 = "
                        "mappable autosomal hg19 (chr1-22), matching this "
                        "study's autosome-only scope. Adjust if your tested "
                        "regions span a different reference.")
    return p.parse_args()


def parse_mark_contrast(fname):
    """Pull a human label from a csaw filename like
    'H3K9me3_B10-none_vs_Parental-none_libsize_all_regions.csv'.
    Returns (mark, contrast_string)."""
    base = os.path.basename(fname)
    base = re.sub(r"_(libsize|spikein)?_?all_regions\.csv$", "", base)
    # mark is the leading token up to the first '_'
    parts = base.split("_", 1)
    mark = parts[0]
    contrast = parts[1] if len(parts) > 1 else ""
    contrast = contrast.replace("_vs_", " vs ").replace("-", "-")
    return mark, contrast


def main():
    args = parse_args()
    outdir = args.outdir or os.path.join(args.csaw_dir, "summary")
    os.makedirs(outdir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(args.csaw_dir, args.pattern)))
    if not files:
        sys.exit(f"[ERROR] No files matching '{args.pattern}' in {args.csaw_dir}.")
    print(f"[INFO] Found {len(files)} csaw region file(s):")
    for f in files:
        print("  ", os.path.basename(f))

    rows = []
    for f in files:
        df = pd.read_csv(f)
        if args.lfc_col not in df.columns or args.fdr_col not in df.columns:
            print(f"  [WARN] {os.path.basename(f)} lacks "
                  f"'{args.lfc_col}' or '{args.fdr_col}'; skipping.")
            continue
        lfc = pd.to_numeric(df[args.lfc_col], errors="coerce")
        fdr = pd.to_numeric(df[args.fdr_col], errors="coerce")
        n_tested = int(lfc.notna().sum())

        sig_mask = (fdr < args.fdr) & (lfc.abs() >= args.lfc)
        n_sig = int(sig_mask.sum())
        n_up = int((sig_mask & (lfc > 0)).sum())     # gained in group1
        n_down = int((sig_mask & (lfc < 0)).sum())   # lost in group1
        pct_up = 100 * n_up / n_sig if n_sig else np.nan
        pct_down = 100 * n_down / n_sig if n_sig else np.nan
        # Percent OF TESTED regions (comparable across marks regardless of how
        # many regions each mark tested): this is what the figure plots.
        pct_up_tested = 100 * n_up / n_tested if n_tested else np.nan
        pct_down_tested = 100 * n_down / n_tested if n_tested else np.nan
        pct_sig_tested = 100 * n_sig / n_tested if n_tested else np.nan

        # Genomic FOOTPRINT: bp covered by ALL tested/enriched regions (the
        # mark's distribution breadth), interval-merged so genome isn't
        # double-counted. Condition-AGNOSTIC: this is a between-marks breadth
        # comparison, NOT a between-condition quantity claim.
        cov_bp = merged_coverage_bp(df)
        cov_pct = (100 * cov_bp / args.genome_size
                   if (cov_bp == cov_bp and args.genome_size) else np.nan)

        mark, contrast = parse_mark_contrast(f)
        rows.append({"mark": mark, "contrast": contrast,
                     "n_tested": n_tested, "n_sig": n_sig,
                     "n_up": n_up, "n_down": n_down,
                     "pct_up": pct_up, "pct_down": pct_down,
                     "pct_up_tested": pct_up_tested,
                     "pct_down_tested": pct_down_tested,
                     "pct_sig_tested": pct_sig_tested,
                     "coverage_bp": cov_bp, "coverage_pct_genome": cov_pct,
                     "file": os.path.basename(f)})
        cov_str = (f", footprint={cov_bp/1e6:.1f} Mb ({cov_pct:.1f}% genome)"
                   if cov_bp == cov_bp else "")
        print(f"  [INFO] {mark}: {n_sig} sig / {n_tested} tested "
              f"(up={n_up} [{pct_up:.0f}%], down={n_down} [{pct_down:.0f}%])"
              f"{cov_str}")

    if not rows:
        sys.exit("[ERROR] No usable region files (missing logFC/FDR columns).")

    summ = pd.DataFrame(rows).sort_values("n_sig", ascending=True).reset_index(drop=True)
    summ.to_csv(os.path.join(outdir, "csaw_differential_summary.txt"),
                sep="\t", index=False)
    print(f"[INFO] Summary table -> "
          f"{os.path.join(outdir, 'csaw_differential_summary.txt')}")

    # ── Auto-drop zero-significance marks from the FIGURE (kept in the table) ─
    # A mark with 0 significant regions (e.g. H3K27Ac under its batch effect)
    # renders as an empty 'nan%' row. Drop it from the plot and record a note
    # naming which mark(s) were dropped, so the omission is explicit, not hidden.
    zero_marks = summ.loc[summ["n_sig"] == 0, "mark"].tolist()
    plot_df = summ[summ["n_sig"] > 0].copy().reset_index(drop=True)
    note = ""
    if zero_marks:
        note = ("Omitted (0 significant regions at this threshold): "
                + ", ".join(zero_marks))
        print(f"[INFO] {note}")
    if plot_df.empty:
        sys.exit("[ERROR] Every mark had 0 significant regions; nothing to plot.")

    # ── Diverging bar chart (x-axis = ABSOLUTE count of significant regions) ──
    # Bar length = number of significant regions gained (right) / lost (left)
    # in group1. Each bar is labelled with its count AND the tested denominator
    # (n = sig / tested), so magnitude and base are both visible at a glance.
    # Direction (left/right split) is the instant takeaway.
    multi_contrast = plot_df["contrast"].nunique() > 1
    labels = [(f"{r.mark}\n{r.contrast}" if multi_contrast else r.mark)
              for r in plot_df.itertuples()]
    y = np.arange(len(plot_df))

    fig_h = max(2.5, 0.85 * len(plot_df) + 1.6)
    fig, ax = plt.subplots(figsize=(9, fig_h))

    ax.barh(y, plot_df["n_up"], color=COLOR_UP, edgecolor="white",
            height=0.62, label=f"Gained in {args.group1_label}")
    ax.barh(y, -plot_df["n_down"], color=COLOR_DOWN, edgecolor="white",
            height=0.62, label=f"Lost in {args.group1_label}")

    xmax = max(plot_df["n_up"].max(), plot_df["n_down"].max())
    pad = xmax * 0.02
    for i, r in plot_df.iterrows():
        # Bar ends: absolute count of regions gained / lost.
        if r["n_up"] > 0:
            ax.text(r["n_up"] + pad, i, f"{r['n_up']:,} ({r['pct_up']:.0f}%)",
                    va="center", ha="left", fontsize=8, color=COLOR_UP)
        if r["n_down"] > 0:
            ax.text(-r["n_down"] - pad, i, f"{r['n_down']:,} ({r['pct_down']:.0f}%)",
                    va="center", ha="right", fontsize=8, color=COLOR_DOWN)
        # Top of bar: total significant over tested denominator (fast to read).
        ax.text(0, i + 0.33,
                f"n = {r['n_sig']:,} / {r['n_tested']:,} tested",
                va="bottom", ha="center", fontsize=7.5, color="#555555")

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel(f"Significant differential regions "
                  f"(FDR<{args.fdr}, |log2FC|\u2265{args.lfc})\n"
                  f"\u2190 lost in {args.group1_label}      "
                  f"gained in {args.group1_label} \u2192", labelpad=8)
    ax.set_title("csaw differential binding by mark\n"
                 f"(diverging: gain vs loss of signal in {args.group1_label})",
                 fontweight="bold", pad=10)
    ax.legend(frameon=False, fontsize=9, loc="lower right")

    # Symmetric x-limits so the zero line is the true center; absolute-value ticks.
    lim = xmax * 1.30
    ax.set_xlim(-lim, lim)
    ticks = ax.get_xticks()
    ax.set_xticklabels([f"{abs(int(t)):,}" for t in ticks])

    # NOTE: the zero-significance omission is intentionally NOT drawn on the
    # figure (console-only, printed above) to keep the slide clean.

    plt.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(outdir, f"csaw_differential_summary.{ext}")
        fig.savefig(out, bbox_inches="tight", pad_inches=0.05, dpi=600)
        print(f"  Saved -> {out}")
    plt.close(fig)

    # ── Coverage figure: genomic footprint (% genome) per mark ───────────────
    # Distribution BREADTH comparison BETWEEN marks (condition-agnostic). Uses
    # ALL tested/enriched regions, interval-merged. This is the safe footprint
    # claim ("H3K9me3 is broadly distributed"), NOT a between-condition quantity
    # claim. Includes every mark with a computable footprint (even those with 0
    # significant differential regions — they still occupy the genome).
    cov_df = summ[summ["coverage_bp"].notna() & (summ["coverage_bp"] > 0)].copy()
    cov_df = cov_df.sort_values("coverage_pct_genome", ascending=True).reset_index(drop=True)
    if not cov_df.empty:
        cov_labels = [(f"{r.mark}\n{r.contrast}"
                       if cov_df["contrast"].nunique() > 1 else r.mark)
                      for r in cov_df.itertuples()]
        yc = np.arange(len(cov_df))
        fig2_h = max(2.2, 0.7 * len(cov_df) + 1.3)
        fig2, ax2 = plt.subplots(figsize=(7.5, fig2_h))
        ax2.barh(yc, cov_df["coverage_pct_genome"], color="#5D6D7E",
                 edgecolor="white", height=0.6)
        cmax = cov_df["coverage_pct_genome"].max()
        for i, r in cov_df.iterrows():
            ax2.text(r["coverage_pct_genome"] + cmax * 0.02, i,
                     f"{r['coverage_pct_genome']:.1f}%  "
                     f"({r['coverage_bp']/1e6:,.0f} Mb)",
                     va="center", ha="left", fontsize=8, color="#34495E")
        ax2.set_yticks(yc); ax2.set_yticklabels(cov_labels, fontsize=9)
        ax2.set_xlabel("Genomic footprint: % of analyzed genome covered by "
                       "enriched regions", labelpad=8)
        ax2.set_title("Distribution breadth by mark\n"
                      "(all tested/enriched regions, interval-merged; "
                      "condition-agnostic)",
                      fontweight="bold", pad=10)
        ax2.set_xlim(0, cmax * 1.30)
        plt.tight_layout()
        for ext in ("pdf", "png"):
            out = os.path.join(outdir, f"csaw_genomic_footprint.{ext}")
            fig2.savefig(out, bbox_inches="tight", pad_inches=0.05, dpi=600)
            print(f"  Saved -> {out}")
        plt.close(fig2)

    print("\n[INFO] Done. NOTE for interpretation: under TMM normalization the "
          "gain/loss split is a RELATIVE redistribution among differential "
          "regions, not a genome-wide-increase claim (which needs spike-in). "
          "The footprint figure is a BETWEEN-MARKS breadth comparison "
          "(condition-agnostic) — do NOT read it as a between-condition "
          "quantity change.")


if __name__ == "__main__":
    main()