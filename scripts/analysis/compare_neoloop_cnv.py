#!/usr/bin/env python3
"""
compare_neoloop_cnv.py — independent check of our Hi-C CNV calls with NeoLoopFinder.

Takes the four NeoLoopFinder `calculate-cnv` bedGraphs and builds the SAME four
comparisons as hic_cnv_ratio.py (2 test + 2 null), plotted in the SAME layout and
colours, with our CNV segments shaded. Then quantifies agreement:

  * per our CNV segment: NeoLoopFinder B10-vs-Parental difference in each pair
  * genome-wide Spearman correlation between NeoLoopFinder and our ratio track

Scale handling: NeoLoopFinder values are converted to log2 relative to each
sample's genome median (log2(value/median) if all values are positive,
value - median otherwise), so pairs are compared on the same log2 scale as ours.

Outputs (only aggregated numbers, safe to share): neoloop_vs_ours_summary.txt,
neoloop_cnv_genome.png/.pdf, neoloop_segment_check.tsv

Example (multiomics-py env):
python scripts/analysis/compare_neoloop_cnv.py \
  --par-dmso results/exploratory/hic_cnv/neoloop/Parental_DMSO.cnv.bedGraph \
  --par-unc  results/exploratory/hic_cnv/neoloop/Parental_UNC.cnv.bedGraph \
  --b10-dmso results/exploratory/hic_cnv/neoloop/B10_DMSO.cnv.bedGraph \
  --b10-unc  results/exploratory/hic_cnv/neoloop/B10_UNC.cnv.bedGraph \
  --segments results/exploratory/hic_cnv/res_250000/cnv_segments.tsv \
  --our-bins results/exploratory/hic_cnv/res_250000/cnv_bins.tsv \
  --outdir   results/exploratory/hic_cnv/neoloop
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# shared visual language with hic_cnv_ratio.py / rnaseq_cnv_shift.py
SERIES = ["#2a78d6", "#eb6834"]
NULL_COLOR = "#8a8984"
INK, INK_2, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"
LOSS_TINT, LOSS_ALPHA = "#4a3aa7", 0.13
GAIN_TINT, GAIN_ALPHA = "#eda100", 0.22
HG19 = {"chr1": 249250621, "chr2": 243199373, "chr3": 198022430, "chr4": 191154276, "chr5": 180915260,
        "chr6": 171115067, "chr7": 159138663, "chr8": 146364022, "chr9": 141213431, "chr10": 135534747,
        "chr11": 135006516, "chr12": 133851895, "chr13": 115169878, "chr14": 107349540, "chr15": 102531392,
        "chr16": 90354753, "chr17": 81195210, "chr18": 78077248, "chr19": 59128983, "chr20": 63025520,
        "chr21": 48129895, "chr22": 51304566, "chrX": 155270560}
CHROMS = list(HG19)


def read_bg(path: str) -> pd.Series:
    if not os.path.exists(path):
        sys.exit(f"ERROR: not found: {path}")
    d = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[0, 1, 2, 3],
                    names=["chrom", "start", "end", "v"])
    d = d[pd.to_numeric(d["start"], errors="coerce").notna()]          # drop track/header lines
    d["start"] = d["start"].astype(int)
    d["v"] = pd.to_numeric(d["v"], errors="coerce")
    d["chrom"] = np.where(d["chrom"].astype(str).str.startswith("chr"), d["chrom"], "chr" + d["chrom"].astype(str))
    d = d[d["chrom"].isin(CHROMS)].dropna(subset=["v"])
    v = d["v"].to_numpy(float)
    if (v > 0).all():
        lv = np.log2(v / np.median(v))
        scale = "ratio -> log2(value/median)"
    else:
        lv = v - np.median(v)
        scale = "log-like -> value - median"
    s = pd.Series(lv, index=pd.MultiIndex.from_arrays([d["chrom"], d["start"]]))
    s = s[~s.index.duplicated()]
    s.attrs["scale"] = scale
    return s


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("par-dmso", "par-unc", "b10-dmso", "b10-unc"):
        ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--segments", required=True, help="cnv_segments.tsv (our calls)")
    ap.add_argument("--our-bins", default=None, help="cnv_bins.tsv (our per-bin ratios) for correlation")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    s = {"PD": read_bg(a.par_dmso), "PU": read_bg(a.par_unc), "BD": read_bg(a.b10_dmso), "BU": read_bg(a.b10_unc)}
    idx = s["PD"].index
    for k in s:
        idx = idx.intersection(s[k].index)
    if len(idx) < 100:
        sys.exit("ERROR: bedGraphs share too few bins - were they run at the same resolution?")
    pairs = {  # same names and order as our Hi-C figure
        "B10vsPar_DMSO": ("test", s["BD"][idx] - s["PD"][idx]),
        "B10vsPar_UNC": ("test", s["BU"][idx] - s["PU"][idx]),
        "Par_UNCvsDMSO": ("null", s["PU"][idx] - s["PD"][idx]),
        "B10_UNCvsDMSO": ("null", s["BU"][idx] - s["BD"][idx]),
    }
    bins = pd.DataFrame({"chrom": idx.get_level_values(0), "start": idx.get_level_values(1)})
    for name, (_, tr) in pairs.items():
        bins[name] = tr.to_numpy() - np.nanmedian(tr.to_numpy())   # centre like our ratio
    bins = bins.sort_values(["chrom", "start"], key=lambda c: c.map(CHROMS.index) if c.name == "chrom" else c)
    bins = bins.reset_index(drop=True)
    res = int(np.median(np.diff(bins.loc[bins["chrom"] == bins["chrom"].iloc[0], "start"])))

    seg = pd.read_csv(a.segments, sep="\t")
    alt = seg[seg["call"] != "neutral"].copy()

    # ---- per-segment check ----
    rows = []
    for r in alt.itertuples():
        m = (bins["chrom"] == r.chrom) & (bins["start"] >= r.start) & (bins["start"] < r.end)
        row = {"chrom": r.chrom, "start": r.start, "end": r.end, "call": r.call,
               "ours_log2": round(r.mean_all, 3), "n_bins": int(m.sum())}
        for name in pairs:
            row[f"neoloop_{name}"] = round(float(bins.loc[m, name].mean()), 3) if m.any() else np.nan
        t = [row["neoloop_B10vsPar_DMSO"], row["neoloop_B10vsPar_UNC"]]
        want = -1 if r.call == "loss" else 1
        agree = all(np.isfinite(t)) and all(np.sign(x) == want and abs(x) >= 0.15 for x in t)
        nulls_flat = all(abs(row[f"neoloop_{n}"]) < 0.15 for n in ("Par_UNCvsDMSO", "B10_UNCvsDMSO")
                         if np.isfinite(row[f"neoloop_{n}"]))
        row["neoloop_agrees"] = "YES" if agree and nulls_flat else ("direction only" if all(np.sign(x) == want for x in t if np.isfinite(x)) else "NO")
        rows.append(row)
    chk = pd.DataFrame(rows)
    chk.to_csv(os.path.join(a.outdir, "neoloop_segment_check.tsv"), sep="\t", index=False)

    # ---- genome-wide agreement with our ratio track ----
    corr_lines = []
    if a.our_bins and os.path.exists(a.our_bins):
        ob = pd.read_csv(a.our_bins, sep="\t")
        ob = ob[ob["valid"].astype(str) == "True"]
        mg = bins.merge(ob, on=["chrom", "start"], how="inner", suffixes=("", "_ours"))
        for name in pairs:
            col = f"log2_{name}"
            if col in mg.columns:
                ok = mg[[name, col]].dropna()
                rho = spearmanr(ok[name], ok[col]).correlation
                corr_lines.append(f"Spearman rho NeoLoopFinder vs ours, {name}: {rho:.2f}  (n bins = {len(ok):,})")

    with open(os.path.join(a.outdir, "neoloop_vs_ours_summary.txt"), "w") as fh:
        fh.write(f"NeoLoopFinder bins compared: {len(bins):,} (~{res // 1000} kb)\n")
        fh.write("scale: " + "; ".join(f"{k}: {v.attrs['scale']}" for k, v in s.items()) + "\n\n")
        fh.write("\n".join(corr_lines) + ("\n\n" if corr_lines else ""))
        fh.write("Our CNV segments re-measured with NeoLoopFinder (log2, B10 minus Parental):\n")
        fh.write(chk.to_string(index=False) + "\n\n")
        fh.write("YES = both test pairs |log2| >= 0.15 in the called direction AND both nulls |log2| < 0.15.\n"
                 "NeoLoopFinder corrects GC/mappability/restriction-site bias with a Poisson GAM, so agreement\n"
                 "means our calls are not an artefact of the simple coverage ratio.\n")
    print(open(os.path.join(a.outdir, "neoloop_vs_ours_summary.txt")).read())

    # ---- plot: identical layout to hic_cnv_ratio.py ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    present = [c for c in CHROMS if (bins["chrom"] == c).any()]
    offsets, mids, cum = {}, [], 0
    for c in present:
        offsets[c] = cum
        mids.append(cum + HG19[c] / 2)
        cum += HG19[c]
    gx = bins["chrom"].map(offsets).to_numpy(float) + bins["start"].to_numpy() + res / 2
    fig, axes = plt.subplots(4, 1, figsize=(15, 2.6 * 4), sharex=True, facecolor=SURFACE)
    n_test = 0
    for ax, (name, (kind, _)) in zip(axes, pairs.items()):
        ax.set_facecolor(SURFACE)
        for r in alt.itertuples():
            if r.chrom in offsets:
                ax.axvspan(offsets[r.chrom] + r.start, offsets[r.chrom] + r.end, lw=0, zorder=0,
                           color=GAIN_TINT if r.call == "gain" else LOSS_TINT,
                           alpha=GAIN_ALPHA if r.call == "gain" else LOSS_ALPHA)
        col = SERIES[n_test] if kind == "test" else NULL_COLOR
        n_test += kind == "test"
        ax.scatter(gx, np.clip(bins[name], -1.5, 1.5), s=2, color=col, alpha=0.5, linewidths=0, rasterized=True)
        for c in present:  # per-chromosome rolling median as the summary line
            g = bins[bins["chrom"] == c]
            roll = g[name].rolling(9, center=True, min_periods=3).median()
            ax.plot(offsets[c] + g["start"] + res / 2, roll, color=INK, lw=1)
            ax.axvline(offsets[c], color="#d9d8d3", lw=0.6, zorder=0)
        ax.axhline(0, color=INK_2, lw=0.8)
        ax.set_ylim(-1.5, 1.5)
        ax.set_ylabel("log2 (B10 - Par)" if kind == "test" else "log2 (UNC - DMSO)", color=INK_2, fontsize=9)
        what = "test: B10 vs Parental" if kind == "test" else "null: same cells, UNC vs DMSO - should be flat"
        ax.set_title(f"{name}  ({what})", loc="left", fontsize=10, color=INK)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(colors=INK_2, labelsize=8)
    axes[-1].set_xlim(0, cum)
    axes[-1].set_xticks(mids)
    axes[-1].set_xticklabels([c.replace("chr", "") for c in present], fontsize=8, color=INK_2)
    fig.suptitle("NeoLoopFinder CNV (independent method) - black: rolling median; shading: OUR CNV calls "
                 "(violet loss, yellow gain)", x=0.01, ha="left", fontsize=12, color=INK)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(a.outdir, f"neoloop_cnv_genome.{ext}"), dpi=200, facecolor=SURFACE)


if __name__ == "__main__":
    main()
