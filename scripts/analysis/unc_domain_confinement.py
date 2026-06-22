#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unc_domain_confinement.py
===========================================

Does constricted migration CONFINE H3K27me3 domain architecture -- i.e. does the
mark coalesce into fewer, tighter, more contiguous domains -- and does UNC1999
(EZH2/EZH1 inhibition) RELAX that confinement back toward the parental layout?

Motivated by Krug et al. 2026 Nat Genet (https://doi.org/10.1038/s41588-026-02586-y): 
H3K27me3 confinement / cPRC1 under constriction.

The 2x2 it reads
----------------
From the Question-A master table, the four per-bin normalized H3K27me3 tracks:
    parental_none , parental_unc , b10_none , b10_unc
and asks, for each architecture metric M:
    migration confinement  = M(b10_none)      - M(parental_none)
    UNC reversal in B10     = M(b10_unc)       - M(b10_none)
    UNC effect in parental  = M(parental_unc)  - M(parental_none)   (control)
    reversal fraction (B10) = -[M(b10_unc)-M(b10_none)] / [M(b10_none)-M(parental_none)]
                              (1 = UNC fully undoes the migration change; 0 = no reversal;
                               <0 = UNC pushes further the same way migration did)

Why this is scale-free -- EXACTLY, not just in contrast (the whole point)
-------------------------------------------------------------------------
The unknown global H3K27me3 offset is a positive multiplicative constant on each
track (or additive in log).  Every metric here is invariant to that constant:

  * Gini concentration          : Gini(c*x) = Gini(x) for c>0           (exact)
  * rank spatial autocorrelation: built on within-track RANKS            (exact)
  * domain segmentation         : domains = top-q bins BY RANK, so the
                                  number of domain bins is identical in every
                                  condition; only their spatial CLUSTERING varies (exact)

So a difference in any metric between two tracks reflects a difference in spatial
ORGANIZATION, never in the (unrecoverable) global amount.  The self-test multiplies
a whole track by 10 and asserts every metric is unchanged to numerical precision.

Metrics (all "higher = more confined", except n_domains where lower = more confined)
------------------------------------------------------------------------------------
  gini          concentration of signal across bins (Lorenz/Gini)
  rank_autocorr lag-1 autocorrelation of within-track ranks along the genome
                (smoother/blockier track -> higher)
  n_domains     number of contiguous domains among the top-q bins (fewer = confined)
  mean_run_len  mean domain length in bins (longer = confined)
  frac_in_large fraction of domain bins sitting in "large" domains (>= --large-domain-bins)

Inference
---------
Block / chromosome bootstrap for CIs (default: resample whole chromosomes, which
preserves all within-chromosome spatial structure exactly), and a block-swap
permutation for each paired between-condition difference (swap the two tracks
within blocks under the no-effect null, preserving spatial structure).

Environment
-----------
conda env: multiomics-py   (numpy, scipy, pandas, matplotlib). No deeptools / BAMs.

------------------------------------------------------------------------------
USAGE EXAMPLES
------------------------------------------------------------------------------
# 0) Validate the scale-free metrics + inference synthetically (no data needed):
python scripts/analysis/unc_domain_confinement.py --self-test

# 1) Standard run against the Question-A master table:
python scripts/analysis/unc_domain_confinement.py \
    --master  results/exploratory/unc_selective_erasure/master_bins_H3K27me3_250kb.tsv \
    --outdir  results/exploratory/unc_domain_confinement

# 2) Tighter domain definition (top 15% bins, "large" >= 1.5 Mb = 6 bins) and a
#    block bootstrap instead of whole-chromosome:
python scripts/analysis/unc_domain_confinement.py \
    --master  .../master_bins_H3K27me3_250kb.tsv \
    --top-q 0.15 --large-domain-bins 6 --bootstrap block --block-len 20 \
    --outdir  results/exploratory/unc_domain_confinement_tight

NOTE on the bin universe: domain contiguity depends on which bins are present. If
the master was built with --min-baseline-quantile > 0 (occupancy filter), dropped
low-signal bins become genomic gaps that legitimately break runs (they are not
H3K27me3 domains anyway). For the definitive domain figure you may also run a
master built with --min-baseline-quantile 0 and compare; the metrics are robust
to this because domains are the TOP-q bins regardless.
------------------------------------------------------------------------------
"""

import argparse
import os
import sys
import subprocess
import datetime as _dt

import numpy as np
import pandas as pd
from scipy.stats import rankdata

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


COND_COLS = ["parental_none", "parental_unc", "b10_none", "b10_unc"]
COND_LABEL = {"parental_none": "Parental DMSO", "parental_unc": "Parental UNC",
              "b10_none": "B10 DMSO", "b10_unc": "B10 UNC"}
REQUIRED = ["chr", "bed_start", "bed_end"] + COND_COLS
METRICS = ["rank_autocorr", "n_domains", "mean_run_len", "frac_in_large"]
# direction in which each metric INCREASES with confinement (+1) or decreases (-1)
# (all four are SPATIAL: built on within-track ranks / top-q domain calling, so
#  they are exactly invariant to any positive global scaling of a track. A
#  concentration metric such as Gini is deliberately NOT used: it is permutation-
#  invariant and therefore blind to spatial domain architecture.)
CONFINE_SIGN = {"rank_autocorr": +1, "n_domains": -1,
                "mean_run_len": +1, "frac_in_large": +1}


# =============================================================================
# Metric primitives (all scale-invariant)
# =============================================================================
def _frame_arrays(frames, col):
    """Pool values across frames, compute genome-wide pct RANKS once (scale-free),
    and hand back per-frame (rank, start, end) so spatial metrics use ranks that
    are consistent across the whole analysed set."""
    vals = [np.asarray(f[col].values, float) for f in frames]
    starts = [np.asarray(f["bed_start"].values, np.int64) for f in frames]
    ends = [np.asarray(f["bed_end"].values, np.int64) for f in frames]
    allv = np.concatenate(vals) if vals else np.array([])
    if len(allv) == 0:
        return allv, []
    pr = rankdata(allv) / len(allv)        # pct ranks, average ties
    out, k = [], 0
    for v in vals:
        out.append(pr[k:k + len(v)]); k += len(v)
    return allv, list(zip(out, starts, ends))


def compute_metrics(frames, col, top_q, large_len):
    """All five architecture metrics for one track (one condition).

    frames : list of per-chromosome (or per-block) DataFrames, each sorted by
             start, with columns bed_start/bed_end and the value column `col`.
    """
    allv, per = _frame_arrays(frames, col)
    if len(allv) == 0:
        return {m: np.nan for m in METRICS}

    # ---- lag-1 rank autocorrelation along the genome (adjacent bins only) ----
    a_acc, b_acc = [], []
    # ---- domain runs among top-q bins (genomically contiguous) ----
    n_dom_bins = 0; n_runs = 0; run_lengths = []
    thr = 1.0 - top_q
    for r, s, e in per:
        if len(r) == 0:
            continue
        adj = np.r_[False, s[1:] == e[:-1]]          # adj[i]: bin i abuts bin i-1
        # autocorr pairs
        m = s[1:] == e[:-1]
        if m.any():
            a_acc.append(r[:-1][m]); b_acc.append(r[1:][m])
        # domains
        dom = r >= thr
        if dom.any():
            new = dom & ~(np.r_[False, dom[:-1]] & adj)   # run-start markers
            n_runs += int(new.sum())
            n_dom_bins += int(dom.sum())
            runid = np.where(dom, np.cumsum(new), 0)
            if runid.max() > 0:
                lens = np.bincount(runid)[1:]             # drop id 0 (non-domain)
                run_lengths.append(lens[lens > 0])

    if a_acc:
        A = np.concatenate(a_acc); B = np.concatenate(b_acc)
        ac = np.corrcoef(A, B)[0, 1] if (np.std(A) > 0 and np.std(B) > 0) else np.nan
    else:
        ac = np.nan

    if run_lengths:
        L = np.concatenate(run_lengths)
        mean_run = float(L.mean())
        frac_large = float(L[L >= large_len].sum()) / float(L.sum()) if L.sum() else np.nan
    else:
        mean_run = np.nan; frac_large = np.nan

    return {"rank_autocorr": ac, "n_domains": float(n_runs),
            "mean_run_len": mean_run, "frac_in_large": frac_large}


# =============================================================================
# Block construction + bootstrap / permutation
# =============================================================================
def make_blocks(df, mode, block_len):
    """Split the (chr-sorted) master into resampling units.
      chrom : one block per chromosome (preserves all within-chr structure; coarse)
      block : contiguous runs of `block_len` bins within each chromosome (finer)"""
    frames = []
    for _, g in df.groupby("chr", sort=False):
        g = g.sort_values("bed_start")
        if mode == "chrom":
            frames.append(g)
        else:
            for i in range(0, len(g), block_len):
                frames.append(g.iloc[i:i + block_len])
    return [f for f in frames if len(f)]


def bootstrap(frames, cols, top_q, large_len, n_boot, seed):
    """Paired block bootstrap: each iteration resamples the SAME block set for all
    conditions, so per-condition metrics and their differences are paired."""
    rng = np.random.default_rng(seed)
    nb = len(frames)
    store = {c: {m: np.empty(n_boot) for m in METRICS} for c in cols}
    for it in range(n_boot):
        sel = rng.integers(0, nb, nb)
        fr = [frames[i] for i in sel]
        for c in cols:
            mm = compute_metrics(fr, c, top_q, large_len)
            for m in METRICS:
                store[c][m][it] = mm[m]
    return store


def perm_pvalue(frames, colA, colB, top_q, large_len, n_perm, seed):
    """Block-swap permutation p for the paired difference M(colA)-M(colB).
    Under H0 the two tracks are exchangeable; swap them within whole blocks
    (preserving spatial structure) and rebuild the null difference."""
    rng = np.random.default_rng(seed)
    obs = {m: compute_metrics(frames, colA, top_q, large_len)[m]
           - compute_metrics(frames, colB, top_q, large_len)[m] for m in METRICS}
    cnt = {m: 0 for m in METRICS}; valid = {m: 0 for m in METRICS}
    for _ in range(n_perm):
        swapped = []
        for f in frames:
            if rng.random() < 0.5:
                f2 = f.copy()
                f2[colA], f2[colB] = f[colB].values, f[colA].values
                swapped.append(f2)
            else:
                swapped.append(f)
        for m in METRICS:
            da = compute_metrics(swapped, colA, top_q, large_len)[m]
            db = compute_metrics(swapped, colB, top_q, large_len)[m]
            d = da - db
            if np.isfinite(d) and np.isfinite(obs[m]):
                valid[m] += 1
                if abs(d) >= abs(obs[m]):
                    cnt[m] += 1
    p = {m: ((cnt[m] + 1) / (valid[m] + 1) if valid[m] else np.nan) for m in METRICS}
    return obs, p


def ci(arr):
    a = np.asarray(arr, float); a = a[np.isfinite(a)]
    if len(a) == 0:
        return (np.nan, np.nan)
    return tuple(np.percentile(a, [2.5, 97.5]))


# =============================================================================
# Plots
# =============================================================================
COND_COLOR = {"parental_none": "#f4a582", "parental_unc": "#b2182b",
              "b10_none": "#92c5de", "b10_unc": "#2166ac"}


def plot_metric_panels(metric_vals, metric_ci, effects, outpath):
    """One panel per metric: 4 condition bars with bootstrap-CI whiskers, plus a
    one-line readout of the migration effect and the B10 reversal fraction."""
    try:
        ms = METRICS
        fig, axes = plt.subplots(1, len(ms), figsize=(3.0 * len(ms), 4.2),
                                 constrained_layout=True)
        for ax, m in zip(np.atleast_1d(axes), ms):
            xs = np.arange(len(COND_COLS))
            vals = [metric_vals[c][m] for c in COND_COLS]
            los = [metric_vals[c][m] - metric_ci[c][m][0] for c in COND_COLS]
            his = [metric_ci[c][m][1] - metric_vals[c][m] for c in COND_COLS]
            ax.bar(xs, vals, color=[COND_COLOR[c] for c in COND_COLS],
                   edgecolor="black", linewidth=0.6,
                   yerr=[np.abs(los), np.abs(his)], capsize=3, error_kw=dict(lw=1))
            ax.set_xticks(xs)
            ax.set_xticklabels([COND_LABEL[c].replace(" ", "\n") for c in COND_COLS],
                               fontsize=7)
            ax.set_title(m, fontsize=9)
            arrow = "↑confined" if CONFINE_SIGN[m] > 0 else "↓confined"
            ax.set_xlabel(arrow, fontsize=7, color="0.4")
            row = effects.get(("migration", m), {})
            rev = effects.get(("reversal_b10", m), {})
            txt = ""
            if row:
                txt += f"mig Δ={row['diff']:+.3g}{_stars(row['p'])}"
            if rev and np.isfinite(rev.get("diff", np.nan)):
                txt += f"\nrev={rev['diff']:.2f}"
            if txt:
                ax.annotate(txt, xy=(0.5, 0.98), xycoords="axes fraction",
                            ha="center", va="top", fontsize=7)
        fig.suptitle("H3K27me3 domain architecture across the 2×2\n"
                     "(migration confinement and UNC reversal)", fontsize=11)
        fig.savefig(outpath, dpi=150)
        plt.close(fig)
    except Exception as e:
        print(f"[warn] metric panel plot failed: {e}")


def _stars(p):
    if not np.isfinite(p):
        return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def plot_domain_raster(df, chrom, top_q, outpath):
    """Domain occupancy (top-q bins) along one chromosome for the four conditions
    -- the visual of confinement (fewer/longer blocks) and its UNC reversal."""
    try:
        g = df[df["chr"] == chrom].sort_values("bed_start")
        if len(g) < 10:
            print(f"[warn] raster: chromosome {chrom} has <10 bins, skipping")
            return
        x0 = g["bed_start"].values / 1e6
        fig, ax = plt.subplots(figsize=(10, 3.2), constrained_layout=True)
        for row, c in enumerate(COND_COLS):
            r = rankdata(g[c].values) / len(g)
            dom = r >= (1 - top_q)
            ax.broken_barh(
                [(x0[i], (g["bed_end"].values[i] - g["bed_start"].values[i]) / 1e6)
                 for i in np.where(dom)[0]],
                (row - 0.4, 0.8), facecolors=COND_COLOR[c], edgecolor="none")
        ax.set_yticks(range(len(COND_COLS)))
        ax.set_yticklabels([COND_LABEL[c] for c in COND_COLS], fontsize=8)
        ax.set_xlabel(f"{chrom} position (Mb)")
        ax.set_title(f"H3K27me3 domains (top {top_q:.0%} bins) on {chrom}\n"
                     "fewer/longer blocks = more confined")
        ax.set_ylim(-0.6, len(COND_COLS) - 0.4)
        fig.savefig(outpath, dpi=150)
        plt.close(fig)
    except Exception as e:
        print(f"[warn] raster plot failed: {e}")


# =============================================================================
# Provenance
# =============================================================================
def git_hash():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"],
                                       stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "NA"


def write_manifest(outdir, args, n_bins, metric_vals):
    lines = ["# MANIFEST  --  unc_domain_confinement.py (Question C)",
             f"generated  : {_dt.datetime.now().isoformat(timespec='seconds')}",
             f"git_commit : {git_hash()}",
             f"python     : {sys.version.split()[0]}", "", "## parameters"]
    for k, v in sorted(vars(args).items()):
        lines.append(f"  {k}: {v}")
    lines += ["", f"## bins analysed : {n_bins}", "", "## metric values (per condition)"]
    for c in COND_COLS:
        vs = "  ".join(f"{m}={metric_vals[c][m]:.4g}" for m in METRICS)
        lines.append(f"  {COND_LABEL[c]:14s} {vs}")
    with open(os.path.join(outdir, "MANIFEST.md"), "w") as fh:
        fh.write("\n".join(lines))


# =============================================================================
# Main
# =============================================================================
def _load(path, min_baseline_quantile):
    m = pd.read_csv(path, sep="\t")
    miss = [c for c in REQUIRED if c not in m.columns]
    if miss:
        sys.exit(f"[ERROR] master missing columns {miss}; regenerate with "
                 f"unc_selective_erasure.py")
    m = m[np.all([np.isfinite(m[c]) for c in COND_COLS], axis=0)].copy()
    if min_baseline_quantile > 0 and "baseline_b10" in m.columns:
        fb = np.nanquantile(m["baseline_b10"], min_baseline_quantile)
        fp = np.nanquantile(m["baseline_parental"], min_baseline_quantile)
        m = m[(m["baseline_b10"] >= fb) | (m["baseline_parental"] >= fp)].copy()
    m["chr"] = m["chr"].astype(str)
    return m.sort_values(["chr", "bed_start"]).reset_index(drop=True)


def main(args):
    os.makedirs(args.outdir, exist_ok=True)
    df = _load(args.master, args.min_baseline_quantile)
    frames = make_blocks(df, args.bootstrap, args.block_len)
    print(f"[info] {len(df)} bins, {df['chr'].nunique()} chromosomes, "
          f"{len(frames)} {args.bootstrap} blocks")

    # ---- point estimates ----
    metric_vals = {c: compute_metrics(frames, c, args.top_q, args.large_domain_bins)
                   for c in COND_COLS}
    print("[info] architecture metrics (per condition):")
    for c in COND_COLS:
        print(f"    {COND_LABEL[c]:14s} " +
              "  ".join(f"{m}={metric_vals[c][m]:.4g}" for m in METRICS))

    # ---- bootstrap CIs (paired) ----
    boot = bootstrap(frames, COND_COLS, args.top_q, args.large_domain_bins,
                     args.n_boot, args.seed)
    metric_ci = {c: {m: ci(boot[c][m]) for m in METRICS} for c in COND_COLS}

    # ---- effects: migration, UNC-B10 (reversal), UNC-parental (control) ----
    contrasts = [("migration", "b10_none", "parental_none"),
                 ("unc_b10",   "b10_unc",  "b10_none"),
                 ("unc_parental", "parental_unc", "parental_none")]
    effects = {}
    rows = []
    for name, cA, cB in contrasts:
        obs, p = perm_pvalue(frames, cA, cB, args.top_q, args.large_domain_bins,
                             args.n_perm, args.seed + 7)
        for m in METRICS:
            d_boot = boot[cA][m] - boot[cB][m]
            lo, hi = ci(d_boot)
            effects[(name, m)] = {"diff": obs[m], "ci_low": lo, "ci_high": hi,
                                  "p": p[m]}
            rows.append({"effect": name, "metric": m, "diff": obs[m],
                         "ci_low": lo, "ci_high": hi, "perm_p": p[m]})

    # ---- reversal fraction in B10 (paired bootstrap of the ratio) ----
    print("[info] reversal fraction in B10 "
          "(1 = UNC fully reverses migration; 0 = none; <0 = same direction):")
    for m in METRICS:
        mig_b = boot["b10_none"][m] - boot["parental_none"][m]
        unc_b = boot["b10_unc"][m] - boot["b10_none"][m]
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = -unc_b / mig_b
        ratio = ratio[np.isfinite(ratio)]
        mig_obs = metric_vals["b10_none"][m] - metric_vals["parental_none"][m]
        unc_obs = metric_vals["b10_unc"][m] - metric_vals["b10_none"][m]
        rev_obs = (-unc_obs / mig_obs) if mig_obs != 0 else np.nan
        lo, hi = ci(ratio) if len(ratio) else (np.nan, np.nan)
        effects[("reversal_b10", m)] = {"diff": rev_obs, "ci_low": lo, "ci_high": hi,
                                        "p": np.nan}
        rows.append({"effect": "reversal_b10", "metric": m, "diff": rev_obs,
                     "ci_low": lo, "ci_high": hi, "perm_p": np.nan})
        flag = "" if abs(mig_obs) > 1e-9 else "  (migration effect ~0; ratio unstable)"
        print(f"    {m:14s} migΔ={mig_obs:+.3g}  uncΔ={unc_obs:+.3g}  "
              f"reversal={rev_obs:+.2f} [{lo:+.2f},{hi:+.2f}]{flag}")

    pd.DataFrame(rows).to_csv(os.path.join(args.outdir, "architecture_effects.tsv"),
                              sep="\t", index=False)
    mv = pd.DataFrame({c: metric_vals[c] for c in COND_COLS}).T
    mv.index.name = "condition"
    mv.to_csv(os.path.join(args.outdir, "architecture_metrics.tsv"), sep="\t")

    # ---- console summary of the headline contrasts ----
    print("[info] migration confinement & UNC reversal (perm p):")
    for m in METRICS:
        e = effects[("migration", m)]; u = effects[("unc_b10", m)]
        dir_mig = ("confines" if np.sign(e["diff"]) == CONFINE_SIGN[m] else "relaxes")
        dir_unc = ("relaxes" if np.sign(u["diff"]) == -CONFINE_SIGN[m] else "confines")
        print(f"    {m:14s} migration {dir_mig} (Δ={e['diff']:+.3g}, p={e['p']:.3f}); "
              f"UNC {dir_unc} (Δ={u['diff']:+.3g}, p={u['p']:.3f})")

    # ---- plots ----
    plot_metric_panels(metric_vals, metric_ci, effects,
                       os.path.join(args.outdir, "architecture_panels.png"))
    chrom = args.raster_chrom
    if chrom in (None, "auto"):
        chrom = df.groupby("chr").size().idxmax()        # densest chromosome
    plot_domain_raster(df, str(chrom), args.top_q,
                       os.path.join(args.outdir, f"domain_raster_{chrom}.png"))

    write_manifest(args.outdir, args, len(df), metric_vals)
    print(f"[done] effects -> {os.path.join(args.outdir, 'architecture_effects.tsv')}")


# =============================================================================
# Self-test
# =============================================================================
def self_test():
    """(a) injected confinement is detected in the right direction across all
    spatial metrics; (b) every metric is EXACTLY invariant to a x10 global offset;
    (c) a no-change control track gives ~0 effects; (d) an INTERMEDIATE-confinement
    UNC track yields a reversal fraction in (0,1)."""
    print("=== SELF-TEST: Question C confinement core ===")
    rng = np.random.default_rng(0)
    # build a genome: 6 chromosomes x 400 contiguous 250 kb bins
    recs = []
    for ch in range(1, 7):
        base = np.abs(rng.normal(0, 1, 400)) + 0.05
        for i in range(400):
            recs.append((f"chr{ch}", i * 250000, (i + 1) * 250000, base[i]))
    df = pd.DataFrame(recs, columns=["chr", "bed_start", "bed_end", "raw"])

    def smooth(col_in, k):
        """per-chromosome spatial smoothing (raises autocorr / contiguity), then
        sharpen so domains are well defined. Larger k = more confined."""
        out = []
        for _, g in df.groupby("chr", sort=False):
            x = g[col_in].values
            out.append(np.convolve(x, np.ones(k) / k, mode="same") ** 3)
        return np.concatenate(out)

    df["parental_none"] = df["raw"]
    df["b10_none"]      = smooth("raw", 11)     # strongly confined
    df["b10_unc"]       = smooth("raw", 5)      # intermediate -> partial reversal
    df["parental_unc"]  = df["raw"]             # UNC ~ no-op in parental (null control)

    frames = make_blocks(df, "chrom", 20)
    M = {c: compute_metrics(frames, c, 0.25, 4) for c in COND_COLS}

    # (a) migration confines: autocorr up, n_domains down, runs longer
    mig = {m: M["b10_none"][m] - M["parental_none"][m] for m in METRICS}
    print("(a) migration Δ:", {m: round(mig[m], 3) for m in METRICS})
    assert mig["rank_autocorr"] > 0 and mig["n_domains"] < 0 and mig["mean_run_len"] > 0, \
        "injected confinement not detected in the expected direction"

    # (b) EXACT scale invariance: x10 a track -> identical metrics
    df2 = df.copy(); df2["b10_none"] = df2["b10_none"] * 10.0
    M2 = compute_metrics(make_blocks(df2, "chrom", 20), "b10_none", 0.25, 4)
    diffs = {m: abs(M2[m] - M["b10_none"][m]) for m in METRICS}
    print("(b) x10 offset |Δmetric|:", {m: f"{diffs[m]:.2e}" for m in METRICS})
    assert max(diffs.values()) < 1e-9, "a metric was NOT scale-invariant"

    # (c) control: parental_unc == parental_none -> ~0 effect
    ctrl = {m: M["parental_unc"][m] - M["parental_none"][m] for m in METRICS}
    print("(c) parental UNC Δ (should be ~0):", {m: round(ctrl[m], 4) for m in METRICS})
    assert max(abs(v) for v in ctrl.values()) < 1e-6, "control track not null"

    # (d) reversal fraction in (0,1): UNC partly undoes the confinement
    for m in ("rank_autocorr", "mean_run_len"):
        rev = -(M["b10_unc"][m] - M["b10_none"][m]) / mig[m]
        print(f"(d) {m} reversal fraction = {rev:.2f} (expect 0<rev<1)")
        assert 0.0 < rev < 1.0, f"{m} reversal fraction outside partial-reversal range"
    print("ALL SELF-TESTS PASSED")


# =============================================================================
# CLI
# =============================================================================
def build_parser():
    p = argparse.ArgumentParser(
        description="Question C: does migration confine H3K27me3 domain "
                    "architecture and does UNC reverse it? (scale-free).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--self-test", action="store_true",
                   help="run synthetic validation (metrics + scale invariance) and exit")
    p.add_argument("--master",
                   help="master_bins_H3K27me3_250kb.tsv from unc_selective_erasure.py")
    p.add_argument("--outdir", default="results/exploratory/unc_domain_confinement")
    p.add_argument("--top-q", type=float, default=0.25,
                   help="top fraction of bins (by within-track rank) called 'domain'; "
                        "fixed per condition so 'amount' is matched and only spatial "
                        "organisation is compared")
    p.add_argument("--large-domain-bins", type=int, default=4,
                   help="a domain run of this many bins (x250 kb) or more counts as "
                        "'large' for frac_in_large (4 = 1 Mb)")
    p.add_argument("--bootstrap", choices=["chrom", "block"], default="chrom",
                   help="resampling unit: 'chrom' (whole chromosomes; preserves all "
                        "within-chr structure exactly) or 'block' (--block-len bins)")
    p.add_argument("--block-len", type=int, default=20,
                   help="block length in bins for --bootstrap block / permutation")
    p.add_argument("--min-baseline-quantile", type=float, default=0.0,
                   help="extra occupancy filter (master is usually already filtered; "
                        "0 = use as-is)")
    p.add_argument("--raster-chrom", default="auto",
                   help="chromosome for the domain raster figure ('auto' = densest)")
    p.add_argument("--n-boot", type=int, default=1000)
    p.add_argument("--n-perm", type=int, default=1000)
    p.add_argument("--seed", type=int, default=0)
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    if args.self_test:
        self_test(); sys.exit(0)
    if not args.master:
        sys.exit("[ERROR] --master is required (or use --self-test)")
    main(args)