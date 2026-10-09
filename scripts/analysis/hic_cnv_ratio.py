#!/usr/bin/env python3
"""
hic_cnv_ratio.py — differential copy-number (CNV) detection from Hi-C coverage.

Idea
----
Every Hi-C read end lands on a genomic bin, so total reads per bin ~ copies of
that bin x local biases (GC, mappability, restriction sites, ...). When two
samples were made with the SAME protocol, those biases are identical and cancel
in the ratio:

        log2( coverage_B10 / coverage_Parental )  ~  log2( copies_B10 / copies_Parental )

Copy-number changes therefore appear as long, flat, step-like blocks in the
ratio track (x1.5 = +0.58, x0.67 = -0.58, x2 = +1, ...).

Design used here
----------------
* TEST pairs     (e.g. B10_DMSO/Parental_DMSO, B10_UNC/Parental_UNC): where CNV
                  differences are expected. A real CNV must replicate in ALL of them.
* NULL pairs     (e.g. Parental_UNC/Parental_DMSO, B10_UNC/B10_DMSO): same cells,
                  acute treatment -> copy number should NOT differ. Their segment
                  amplitudes set a data-driven noise threshold.

Steps
-----
1. Raw (unbalanced) read coverage per bin from each .mcool.
2. Library-size normalise, log2 ratio per pair, median-centre.
3. Segment the mean of the test-pair tracks per chromosome (binary segmentation,
   robust noise estimate).
4. Threshold = max(--min-log2, largest |segment mean| seen in the null pairs).
5. Call a segment ALTERED if |mean| >= threshold and every test pair agrees in sign
   with at least half the threshold.
6. Optional: test whether compartment-switch bins are enriched in ALTERED segments.

Outputs (in --outdir)
---------------------
  cnv_bins.tsv        per-bin log2 ratios (test + null pairs) and segment assignment
  cnv_segments.tsv    segments with per-pair means, call (gain/loss/neutral), approx ratio
  cnv_summary.txt     replication metrics, threshold, % genome altered, enrichment test
  cnv_genome.png/.pdf genome-wide plot (one panel per test pair + null panel)
  switch_bins_in_cnv.tsv  (only with --switch-bed)

Example
-------
python scripts/preprocessing/hic_cnv_ratio.py \
  --resolution 250000 \
  --test "B10vsPar_DMSO:data/processed/hic/mcool/B10_DMSO.hg19.mcool:data/processed/hic/mcool/Parental_DMSO.hg19.mcool" \
  --test "B10vsPar_UNC:data/processed/hic/mcool/B10_UNC.hg19.mcool:data/processed/hic/mcool/Parental_UNC.hg19.mcool" \
  --null "Par_UNCvsDMSO:data/processed/hic/mcool/Parental_UNC.hg19.mcool:data/processed/hic/mcool/Parental_DMSO.hg19.mcool" \
  --null "B10_UNCvsDMSO:data/processed/hic/mcool/B10_UNC.hg19.mcool:data/processed/hic/mcool/B10_DMSO.hg19.mcool" \
  --outdir results/exploratory/hic_cnv
"""
from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd

try:
    import cooler
except ImportError:  # pragma: no cover
    sys.exit("ERROR: cooler is not installed in this environment (conda activate multiomics-py)")

# --------------------------------------------------------------------------- #
# Plot colours (validated default categorical palette, light surface)
# --------------------------------------------------------------------------- #
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7"]
NULL_COLOR = "#8a8984"
INK = "#0b0b0b"
INK_2 = "#52514e"
SURFACE = "#fcfcfb"
BAND = "#f0efec"


def log(msg: str) -> None:
    print(f"[{datetime.now():%H:%M:%S}] {msg}", flush=True)


# --------------------------------------------------------------------------- #
# I/O
# --------------------------------------------------------------------------- #
def parse_pair(spec: str) -> tuple[str, str, str]:
    parts = spec.split(":")
    if len(parts) != 3:
        sys.exit(f"ERROR: pair must be LABEL:numerator.mcool:denominator.mcool, got '{spec}'")
    return parts[0], parts[1], parts[2]


def open_cooler(path: str, res: int) -> cooler.Cooler:
    if not os.path.exists(path.split("::")[0]):
        sys.exit(f"ERROR: file not found: {path}")
    uri = path if "::" in path else f"{path}::resolutions/{res}"
    try:
        return cooler.Cooler(uri)
    except Exception as e:  # noqa: BLE001
        sys.exit(f"ERROR: cannot open {uri}: {e}")


def raw_coverage(clr: cooler.Cooler, ignore_diags: int, chunksize: int = 10_000_000) -> np.ndarray:
    """Reads per bin from RAW counts (each contact adds one read end to bin1 and one to bin2)."""
    nbins = clr.info["nbins"]
    nnz = clr.info["nnz"]
    chrom_of_bin = clr.bins()["chrom"][:].astype(str).values
    cov = np.zeros(nbins, dtype=np.float64)
    for lo in range(0, nnz, chunksize):
        px = clr.pixels()[lo:min(lo + chunksize, nnz)]
        b1 = px["bin1_id"].to_numpy()
        b2 = px["bin2_id"].to_numpy()
        c = px["count"].to_numpy(dtype=np.float64)
        if ignore_diags > 0:
            near = (chrom_of_bin[b1] == chrom_of_bin[b2]) & ((b2 - b1) < ignore_diags)
            b1, b2, c = b1[~near], b2[~near], c[~near]
        cov += np.bincount(b1, weights=c, minlength=nbins)
        cov += np.bincount(b2, weights=c, minlength=nbins)
    return cov


# --------------------------------------------------------------------------- #
# Segmentation
# --------------------------------------------------------------------------- #
def robust_sd(x: np.ndarray) -> float:
    """Noise SD from first differences (MAD), insensitive to step changes."""
    d = np.diff(x[np.isfinite(x)])
    if d.size < 3:
        return np.nan
    mad = np.median(np.abs(d - np.median(d)))
    return 1.4826 * mad / np.sqrt(2.0)


def binary_segment(x: np.ndarray, sd: float, min_size: int, t_thresh: float) -> list[tuple[int, int]]:
    """Recursive binary segmentation on a gap-free vector. Returns [(lo, hi), ...]."""
    segs: list[tuple[int, int]] = []

    def rec(lo: int, hi: int) -> None:
        n = hi - lo
        if n < 2 * min_size or not np.isfinite(sd) or sd <= 0:
            segs.append((lo, hi))
            return
        y = x[lo:hi]
        cs = np.cumsum(y)
        tot = cs[-1]
        k = np.arange(min_size, n - min_size + 1)  # size of left part
        lm = cs[k - 1] / k
        rm = (tot - cs[k - 1]) / (n - k)
        t = np.abs(lm - rm) / (sd * np.sqrt(1.0 / k + 1.0 / (n - k)))
        i = int(np.argmax(t))
        if t[i] < t_thresh:
            segs.append((lo, hi))
            return
        cut = lo + int(k[i])
        rec(lo, cut)
        rec(cut, hi)

    rec(0, len(x))
    return sorted(segs)


def merge_similar(x: np.ndarray, segs: list[tuple[int, int]], delta: float) -> list[tuple[int, int]]:
    """Merge neighbouring segments whose means differ by < delta (log2 units)."""
    if not segs:
        return segs
    out = [list(segs[0])]
    for lo, hi in segs[1:]:
        plo, phi = out[-1]
        if abs(np.mean(x[plo:phi]) - np.mean(x[lo:hi])) < delta:
            out[-1][1] = hi
        else:
            out.append([lo, hi])
    return [tuple(s) for s in out]


def segment_track(bins: pd.DataFrame, values: np.ndarray, valid: np.ndarray,
                  min_size: int, t_thresh: float, merge_delta: float) -> pd.DataFrame:
    """Segment per chromosome using valid bins only; returns segment table with bin indices."""
    sd = robust_sd(values[valid])
    rows = []
    for chrom, idx in bins.groupby("chrom", sort=False).groups.items():
        idx = np.asarray(idx)
        idx = idx[valid[idx]]
        if idx.size == 0:
            continue
        x = values[idx]
        segs = merge_similar(x, binary_segment(x, sd, min_size, t_thresh), merge_delta)
        for lo, hi in segs:
            members = idx[lo:hi]
            rows.append({
                "chrom": chrom,
                "start": int(bins.at[members[0], "start"]),
                "end": int(bins.at[members[-1], "end"]),
                "n_bins": int(members.size),
                "first_bin": int(members[0]),
                "last_bin": int(members[-1]),
            })
    seg = pd.DataFrame(rows)
    seg.attrs["noise_sd"] = sd
    return seg


def seg_means(seg: pd.DataFrame, bins_valid_idx_map: dict, values: np.ndarray, valid: np.ndarray) -> np.ndarray:
    out = np.full(len(seg), np.nan)
    for i, r in enumerate(seg.itertuples()):
        sl = np.arange(r.first_bin, r.last_bin + 1)
        sl = sl[valid[sl]]
        if sl.size:
            out[i] = float(np.mean(values[sl]))
    return out


# --------------------------------------------------------------------------- #
# Plot
# --------------------------------------------------------------------------- #
def plot_genome(bins, valid, tracks, null_tracks, seg, outprefix, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    chroms = [c for c in dict.fromkeys(bins["chrom"]) if valid[(bins["chrom"] == c).to_numpy()].any()]
    offsets, mids, cum = {}, [], 0
    for c in chroms:
        sub = bins[bins["chrom"] == c]
        offsets[c] = cum
        L = int(sub["end"].max())
        mids.append(cum + L / 2)
        cum += L
    gx = bins["chrom"].map(offsets).to_numpy(dtype=float) + (bins["start"].to_numpy() + bins["end"].to_numpy()) / 2

    n_panels = len(tracks) + (1 if null_tracks else 0)
    fig, axes = plt.subplots(n_panels, 1, figsize=(15, 2.6 * n_panels), sharex=True,
                             facecolor=SURFACE)
    axes = np.atleast_1d(axes)
    ylim = 1.5

    def decorate(ax, label):
        ax.set_facecolor(SURFACE)
        ax.axhline(0, color=INK_2, lw=0.8)
        for c in chroms:
            ax.axvline(offsets[c], color="#d9d8d3", lw=0.6, zorder=0)
        ax.set_ylim(-ylim, ylim)
        ax.set_ylabel("log2 ratio", color=INK_2, fontsize=9)
        ax.set_title(label, loc="left", color=INK, fontsize=10)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.tick_params(colors=INK_2, labelsize=8)

    altered = seg[seg["call"] != "neutral"]
    for ax, (i, (label, vals)) in zip(axes, enumerate(tracks.items())):
        for r in altered.itertuples():
            ax.axvspan(offsets[r.chrom] + r.start, offsets[r.chrom] + r.end, color=BAND, zorder=0)
        y = np.clip(vals, -ylim, ylim)
        ax.scatter(gx[valid], y[valid], s=2, color=SERIES[i % len(SERIES)], alpha=0.5,
                   linewidths=0, rasterized=True)
        for r in seg.itertuples():
            m = getattr(r, f"mean_{label}")
            if np.isfinite(m):
                ax.plot([offsets[r.chrom] + r.start, offsets[r.chrom] + r.end], [m, m],
                        color=INK, lw=2, solid_capstyle="round")
        decorate(ax, f"{label}  (test pair; shaded = called CNV difference)")

    if null_tracks:
        ax = axes[-1]
        for label, vals in null_tracks.items():
            y = np.clip(vals, -ylim, ylim)
            ax.scatter(gx[valid], y[valid], s=2, color=NULL_COLOR, alpha=0.4, linewidths=0,
                       rasterized=True)
        decorate(ax, "Null pairs (same cells, acute treatment) — should be flat: "
                 + ", ".join(null_tracks))

    axes[-1].set_xticks(mids)
    axes[-1].set_xticklabels([c.replace("chr", "") for c in chroms], fontsize=8, color=INK_2)
    axes[-1].set_xlim(0, cum)
    fig.suptitle(title, x=0.01, ha="left", color=INK, fontsize=12)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{outprefix}.{ext}", dpi=200, facecolor=SURFACE)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--test", action="append", required=True,
                    help="LABEL:numerator.mcool:denominator.mcool (repeat; CNV must replicate in all)")
    ap.add_argument("--null", action="append", default=[],
                    help="LABEL:numerator.mcool:denominator.mcool for same-cell pairs (repeat)")
    ap.add_argument("--resolution", type=int, default=250000, help="bin size (default 250000)")
    ap.add_argument("--exclude-chroms", default="chrY,chrM,chrMT",
                    help="comma-separated chromosomes to drop (default chrY,chrM,chrMT)")
    ap.add_argument("--ignore-diags", type=int, default=0,
                    help="drop contacts within N bins of the diagonal (default 0 = keep all reads)")
    ap.add_argument("--min-cov-frac", type=float, default=0.2,
                    help="mask bins with coverage < frac x median in ANY sample (default 0.2)")
    ap.add_argument("--min-bins", type=int, default=8,
                    help="minimum segment length in bins (default 8 = 2 Mb at 250 kb)")
    ap.add_argument("--t-thresh", type=float, default=6.0,
                    help="breakpoint t-statistic threshold (default 6)")
    ap.add_argument("--merge-delta", type=float, default=0.1,
                    help="merge adjacent segments whose means differ < this (log2; default 0.1)")
    ap.add_argument("--min-log2", type=float, default=0.2,
                    help="floor for the calling threshold, |log2 ratio| (default 0.2 ~ x1.15)")
    ap.add_argument("--switch-bed", default=None,
                    help="optional BED/TSV (chrom start end ...) of compartment-switch bins to test")
    ap.add_argument("--n-perm", type=int, default=1000, help="permutations for the switch-overlap test")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    tests = [parse_pair(s) for s in args.test]
    nulls = [parse_pair(s) for s in args.null]

    # ---- load coverage (each file once) ----
    files = sorted({f for _, a, b in tests + nulls for f in (a, b)})
    cov, bins = {}, None
    for f in files:
        log(f"coverage: {f}")
        clr = open_cooler(f, args.resolution)
        b = clr.bins()[["chrom", "start", "end"]][:]
        b["chrom"] = b["chrom"].astype(str)
        if bins is None:
            bins = b.reset_index(drop=True)
        elif len(b) != len(bins) or not (b["start"].values == bins["start"].values).all():
            sys.exit(f"ERROR: bins in {f} differ from the first file — use coolers built on the same chrom sizes")
        cov[f] = raw_coverage(clr, args.ignore_diags)
        log(f"   total read ends: {cov[f].sum():,.0f}")

    excl = {c.strip() for c in args.exclude_chroms.split(",") if c.strip()}
    keep_chrom = ~bins["chrom"].isin(excl).to_numpy()
    valid = keep_chrom.copy()
    for f, v in cov.items():
        med = np.median(v[keep_chrom & (v > 0)])
        valid &= v >= args.min_cov_frac * med
    log(f"valid bins: {valid.sum():,} / {len(bins):,} (excluded chroms: {sorted(excl)})")

    def ratio(num, den):
        r = np.full(len(bins), np.nan)
        a = cov[num] / cov[num][valid].sum()
        b = cov[den] / cov[den][valid].sum()
        r[valid] = np.log2(a[valid] / b[valid])
        r[valid] -= np.median(r[valid])
        return r

    test_tracks = {lab: ratio(a, b) for lab, a, b in tests}
    null_tracks = {lab: ratio(a, b) for lab, a, b in nulls}

    # ---- segmentation on the mean of test tracks ----
    stack = np.vstack(list(test_tracks.values()))
    mean_track = np.full(len(bins), np.nan)
    mean_track[valid] = stack[:, valid].mean(axis=0)
    seg = segment_track(bins, mean_track, valid, args.min_bins, args.t_thresh, args.merge_delta)
    log(f"test segments: {len(seg)}  (noise SD of mean track = {seg.attrs['noise_sd']:.3f})")
    seg["mean_all"] = seg_means(seg, None, mean_track, valid)
    for lab, tr in test_tracks.items():
        seg[f"mean_{lab}"] = seg_means(seg, None, tr, valid)

    # ---- null-derived threshold ----
    null_max = 0.0
    null_report = []
    for lab, tr in null_tracks.items():
        nseg = segment_track(bins, tr, valid, args.min_bins, args.t_thresh, args.merge_delta)
        m = seg_means(nseg, None, tr, valid)
        big = np.abs(m[nseg["n_bins"].to_numpy() >= args.min_bins])
        mx = float(np.nanmax(big)) if big.size else 0.0
        null_max = max(null_max, mx)
        null_report.append((lab, len(nseg), mx))
    thr = max(args.min_log2, null_max)
    log(f"calling threshold |log2| >= {thr:.3f} (floor {args.min_log2}, null max {null_max:.3f})")

    # ---- calls ----
    pair_cols = [f"mean_{lab}" for lab in test_tracks]
    calls = []
    for r in seg.itertuples():
        m = r.mean_all
        per = np.array([getattr(r, c) for c in pair_cols])
        ok = (np.isfinite(m) and abs(m) >= thr and r.n_bins >= args.min_bins
              and np.all(np.sign(per) == np.sign(m)) and np.all(np.abs(per) >= thr / 2))
        calls.append(("gain" if m > 0 else "loss") if ok else "neutral")
    seg["call"] = calls
    seg["fold_change"] = np.round(2 ** seg["mean_all"], 3)
    simple = {"1/2": 0.5, "2/3": 2 / 3, "3/4": 0.75, "4/3": 4 / 3, "3/2": 1.5, "2": 2.0}
    seg["nearest_simple_ratio"] = [
        min(simple, key=lambda k: abs(np.log2(simple[k]) - m)) if c != "neutral" else ""
        for m, c in zip(seg["mean_all"], seg["call"])
    ]
    seg["length_mb"] = (seg["end"] - seg["start"]) / 1e6

    # ---- per-bin table ----
    bin_tab = bins.copy()
    bin_tab["valid"] = valid
    for lab, tr in {**test_tracks, **null_tracks}.items():
        bin_tab[f"log2_{lab}"] = np.round(tr, 4)
    bin_tab["log2_mean_test"] = np.round(mean_track, 4)
    bin_tab["segment_call"] = "neutral"
    for r in seg.itertuples():
        bin_tab.loc[r.first_bin:r.last_bin, "segment_call"] = r.call
    bin_tab.loc[~valid, "segment_call"] = "masked"

    # ---- replication metrics ----
    from scipy.stats import spearmanr, fisher_exact
    lines = []
    labs = list(test_tracks)
    if len(labs) >= 2:
        for i in range(len(labs)):
            for j in range(i + 1, len(labs)):
                rho = spearmanr(test_tracks[labs[i]][valid], test_tracks[labs[j]][valid]).correlation
                lines.append(f"Spearman rho, bins: {labs[i]} vs {labs[j]} = {rho:.3f}")
    nl = list(null_tracks)
    for lab in nl:
        for t in labs:
            rho = spearmanr(null_tracks[lab][valid], test_tracks[t][valid]).correlation
            lines.append(f"Spearman rho, bins: null {lab} vs test {t} = {rho:.3f}  (expected near 0; pairs sharing a sample carry some shared noise)")

    alt = seg[seg["call"] != "neutral"]
    genome_valid_mb = (bins.loc[valid, "end"] - bins.loc[valid, "start"]).sum() / 1e6
    alt_mb = alt["length_mb"].sum()

    # ---- optional compartment-switch enrichment ----
    enrich_lines = []
    if args.switch_bed:
        sw = pd.read_csv(args.switch_bed, sep="\t", header=None, comment="#", usecols=[0, 1, 2],
                         names=["chrom", "start", "end"])
        if not str(sw.iloc[0, 1]).isdigit():  # header row present
            sw = sw.iloc[1:]
        sw["start"] = sw["start"].astype(int)
        sw["end"] = sw["end"].astype(int)
        flag = np.zeros(len(bins), dtype=bool)
        for c, g in sw.groupby("chrom"):
            bidx = np.where(bins["chrom"].to_numpy() == c)[0]
            if bidx.size == 0:
                continue
            bs = bins.loc[bidx, "start"].to_numpy()
            be = bins.loc[bidx, "end"].to_numpy()
            for s, e in zip(g["start"], g["end"]):
                flag[bidx[(bs < e) & (be > s)]] = True
        is_alt = bin_tab["segment_call"].isin(["gain", "loss"]).to_numpy()
        v = valid
        frac_sw = flag[v].mean()
        if frac_sw > 0.5:
            enrich_lines.append(f"WARNING: {100 * frac_sw:.0f}% of valid bins are flagged as switch bins - "
                                "the BED probably contains ALL labelled bins, not only switching ones. "
                                "The test below is not interpretable.")
        chrom_arr = bins["chrom"].to_numpy()
        rng = np.random.default_rng(args.seed)

        def overlap_test(alt_mask, label):
            a = int(np.sum(flag & alt_mask & v)); b = int(np.sum(flag & ~alt_mask & v))
            c = int(np.sum(~flag & alt_mask & v)); d = int(np.sum(~flag & ~alt_mask & v))
            odds, p = fisher_exact([[a, b], [c, d]])
            # circular-shift permutation within each chromosome: keeps the clustering of switch
            # bins and of CNV segments, so it does not assume bins are independent (Fisher does)
            obs = a
            null = np.empty(args.n_perm)
            for k in range(args.n_perm):
                shifted = np.zeros_like(flag)
                for ch in np.unique(chrom_arr[v]):
                    idx = np.where((chrom_arr == ch) & v)[0]
                    shifted[idx] = np.roll(flag[idx], rng.integers(idx.size))
                null[k] = np.sum(shifted & alt_mask & v)
            p_enr = (1 + np.sum(null >= obs)) / (1 + args.n_perm)
            p_dep = (1 + np.sum(null <= obs)) / (1 + args.n_perm)
            return [f"[{label}] switch bins in segments: {a:,} of {a + b:,} ({100 * a / max(a + b, 1):.1f}%) | "
                    f"non-switch: {100 * c / max(c + d, 1):.1f}% | expected by chance: {null.mean():.0f}",
                    f"[{label}] Fisher OR = {odds:.2f}, p = {p:.3g} (assumes independent bins - optimistic)",
                    f"[{label}] circular-permutation p (enrichment) = {p_enr:.4f}, p (depletion) = {p_dep:.4f}"]

        enrich_lines += overlap_test(is_alt, "all CNV")
        for call in ("gain", "loss"):
            m = (bin_tab["segment_call"] == call).to_numpy()
            if m.any():
                enrich_lines += overlap_test(m, call)
        bin_tab["is_switch"] = flag
        bin_tab[flag & is_alt & v].to_csv(os.path.join(args.outdir, "switch_bins_in_cnv.tsv"),
                                          sep="\t", index=False)

    # ---- write ----
    bin_tab.to_csv(os.path.join(args.outdir, "cnv_bins.tsv"), sep="\t", index=False)
    seg_out = seg.drop(columns=["first_bin", "last_bin"])
    seg_out.round(4).to_csv(os.path.join(args.outdir, "cnv_segments.tsv"), sep="\t", index=False)

    with open(os.path.join(args.outdir, "cnv_summary.txt"), "w") as fh:
        fh.write(f"hic_cnv_ratio.py  {datetime.now():%Y-%m-%d %H:%M}\n")
        fh.write(f"resolution: {args.resolution:,} bp; valid bins: {valid.sum():,}\n")
        fh.write("test pairs:\n" + "".join(f"  {l}: {a} / {b}\n" for l, a, b in tests))
        fh.write("null pairs:\n" + "".join(f"  {l}: {a} / {b}\n" for l, a, b in nulls))
        fh.write(f"\nnoise SD (mean test track, from first differences): {seg.attrs['noise_sd']:.3f}\n")
        for lab, n, mx in null_report:
            fh.write(f"null {lab}: {n} segments, max |segment mean| = {mx:.3f}\n")
        fh.write(f"calling threshold: |log2 ratio| >= {thr:.3f} (x{2 ** thr:.2f})\n\n")
        fh.write("\n".join(lines) + "\n\n")
        fh.write(f"CNV-changed segments: {len(alt)} "
                 f"(gain {int((alt['call'] == 'gain').sum())}, loss {int((alt['call'] == 'loss').sum())})\n")
        fh.write(f"genome altered: {alt_mb:,.1f} Mb of {genome_valid_mb:,.1f} Mb valid "
                 f"({100 * alt_mb / genome_valid_mb:.1f}%)\n")
        if len(alt):
            fh.write("\nlargest CNV-changed segments:\n")
            top = alt.sort_values("length_mb", ascending=False).head(15)
            fh.write(top[["chrom", "start", "end", "length_mb", "call", "mean_all", "fold_change",
                          "nearest_simple_ratio"]].round(3).to_string(index=False) + "\n")
        if enrich_lines:
            fh.write("\ncompartment-switch enrichment in CNV-changed segments:\n")
            fh.write("\n".join(enrich_lines) + "\n")
        fh.write("\nInterpretation notes:\n"
                 "- Real CNVs: long flat blocks, near-simple ratios, replicated in every test pair,\n"
                 "  absent from null pairs. Small patchy shifts can come from compartment changes\n"
                 "  (A regions yield slightly more reads) - these are NOT copy-number calls.\n"
                 "- Ratios are relative (B10 vs Parental); absolute copy number needs ploidy.\n")

    log("plotting")
    plot_genome(bins, valid, test_tracks, null_tracks, seg,
                os.path.join(args.outdir, "cnv_genome"),
                f"Hi-C coverage ratio — differential copy number ({args.resolution // 1000} kb bins)")
    log(f"done -> {args.outdir}")
    print(open(os.path.join(args.outdir, "cnv_summary.txt")).read())


if __name__ == "__main__":
    main()
