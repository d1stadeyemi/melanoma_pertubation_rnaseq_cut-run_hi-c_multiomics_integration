#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unc_domain_confinement.py   --   Question C
===========================================

Does constricted migration CONFINE H3K27me3 domain architecture -- i.e. does the
mark coalesce into fewer, tighter, more contiguous domains -- and does UNC1999
(EZH2/EZH1 inhibition) RELAX that confinement back toward the parental layout?

(Motivated by Krug et al. 2026 Nat Genet: H3K27me3 confinement / cPRC1 under
constriction.)

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


def _stars(p):
    if not np.isfinite(p):
        return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def _sig(eff):
    """A bootstrap CI that excludes 0 = a 'real' effect (used to gate reversal)."""
    lo, hi = eff.get("ci_low", np.nan), eff.get("ci_high", np.nan)
    return np.isfinite(lo) and np.isfinite(hi) and (lo > 0 or hi < 0)


def plot_metric_panels(metric_vals, metric_ci, effects, outpath):
    """One panel per metric: 4 condition bars with bootstrap-CI whiskers.

    Significance is shown on the contrasts that matter: a bracket over each
    DMSO->UNC pair (Parental DMSO vs Parental UNC; B10 DMSO vs B10 UNC) carries
    the UNC-effect stars -- these are the real signal. The migration effect
    (B10 vs Parental DMSO) is reported in the x-label, and the B10 reversal
    fraction is shown ONLY when the migration effect is itself significant
    (otherwise the ratio is undefined and misleading)."""
    try:
        ms = METRICS
        fig, axes = plt.subplots(1, len(ms), figsize=(3.2 * len(ms), 4.6),
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

            # significance brackets over the DMSO->UNC pairs (the UNC effect)
            tops = [vals[i] + his[i] for i in range(len(COND_COLS))]
            ymax = max(tops) if tops else 1.0
            gap = 0.05 * (ymax if ymax > 0 else 1.0)

            def bracket(i, j, eff):
                y = max(tops[i], tops[j]) + gap
                ax.plot([i, i, j, j], [y, y + gap, y + gap, y], lw=1, color="black")
                ax.text((i + j) / 2.0, y + gap, _stars(eff.get("p", np.nan)),
                        ha="center", va="bottom", fontsize=9)
                return y + 2.2 * gap

            top_used = ymax
            # COND_COLS order = [parental_none, parental_unc, b10_none, b10_unc]
            top_used = max(top_used, bracket(0, 1, effects.get(("unc_parental", m), {})))
            top_used = max(top_used, bracket(2, 3, effects.get(("unc_b10", m), {})))
            ax.set_ylim(0, top_used + gap)

            # migration effect -> x-label; reversal only if migration is real
            arrow = "↑ confined" if CONFINE_SIGN[m] > 0 else "↓ confined"
            mig = effects.get(("migration", m), {})
            lab = f"{arrow}\nmig Δ={mig.get('diff', np.nan):+.3g} {_stars(mig.get('p', np.nan))}"
            if _sig(mig):
                rev = effects.get(("reversal_b10", m), {})
                if np.isfinite(rev.get("diff", np.nan)):
                    lab += f"  rev={rev['diff']:.2f}"
            ax.set_xlabel(lab, fontsize=7, color="0.35")

        fig.suptitle("H3K27me3 domain architecture across the 2×2\n"
                     "brackets = UNC effect (DMSO→UNC); * p<0.05, ** p<0.01, *** p<0.001, ns = not sig.",
                     fontsize=10)
        fig.savefig(outpath, dpi=150)
        plt.close(fig)
    except Exception as e:
        print(f"[warn] metric panel plot failed: {e}")


def plot_domain_raster(df, chrom, top_q, outpath, window=None, diff_cellline="b10"):
    """Domain occupancy (top-q bins) along one chromosome for the four conditions,
    plus a 'UNC Δ' track showing which domain bins are LOST vs GAINED under UNC in
    one cell line -- the visual of confinement / consolidation.

    Domain calls are GENOME-WIDE (rank over the whole analysed set), so they match
    the n_domains / mean_run_len metrics exactly. Rows are ordered Parental DMSO,
    Parental UNC, B10 DMSO, B10 UNC top-to-bottom. Use `window=(start_mb,end_mb)`
    to zoom in so individual domains and their lengths are visible."""
    try:
        # genome-wide top-q domain calls (consistent with compute_metrics)
        d = df.copy()
        for c in COND_COLS:
            d[c + "_dom"] = (rankdata(d[c].values) / len(d)) >= (1.0 - top_q)
        g = d[d["chr"] == chrom].sort_values("bed_start")
        if window is not None:
            s_mb, e_mb = window
            g = g[(g["bed_end"].values / 1e6 > s_mb) & (g["bed_start"].values / 1e6 < e_mb)]
        if len(g) < 5:
            print(f"[warn] raster: {chrom}{(' window' if window else '')} has <5 bins, skipping")
            return

        starts = g["bed_start"].values / 1e6
        widths = (g["bed_end"].values - g["bed_start"].values) / 1e6
        n_main = len(COND_COLS)
        fig, ax = plt.subplots(figsize=(11, 3.8), constrained_layout=True)

        # four condition rows (row index = COND_COLS order; y inverted below)
        for i, c in enumerate(COND_COLS):
            dom = g[c + "_dom"].values
            ax.broken_barh([(starts[k], widths[k]) for k in np.where(dom)[0]],
                           (i - 0.4, 0.8), facecolors=COND_COLOR[c], edgecolor="none")

        # fifth row: domain bins LOST (DMSO-only) vs GAINED (UNC-only) under UNC
        cn, cu = (("b10_none", "b10_unc") if diff_cellline == "b10"
                  else ("parental_none", "parental_unc"))
        lost = g[cn + "_dom"].values & ~g[cu + "_dom"].values
        gained = g[cu + "_dom"].values & ~g[cn + "_dom"].values
        drow = n_main
        ax.broken_barh([(starts[k], widths[k]) for k in np.where(lost)[0]],
                       (drow - 0.4, 0.8), facecolors="#d6604d", edgecolor="none")
        ax.broken_barh([(starts[k], widths[k]) for k in np.where(gained)[0]],
                       (drow - 0.4, 0.8), facecolors="#2166ac", edgecolor="none")

        labels = [COND_LABEL[c] for c in COND_COLS] + \
                 [f"{'B10' if diff_cellline == 'b10' else 'Parental'} UNC Δ"]
        ax.set_yticks(range(n_main + 1))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_ylim(-0.6, n_main + 0.6)
        ax.invert_yaxis()                      # Parental DMSO on top
        ax.axhline(n_main - 0.5, color="0.7", lw=0.8, ls=":")   # separate the Δ track
        ax.set_xlabel(f"{chrom} position (Mb)")
        ttl = f"H3K27me3 domains (top {top_q:.0%} bins) on {chrom}"
        if window is not None:
            ttl += f"  [{window[0]:g}-{window[1]:g} Mb]"
        ax.set_title(ttl + "\nfewer/longer blocks = more confined; "
                           "Δ track: red = lost under UNC, blue = gained")
        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(color="#d6604d", label="lost under UNC"),
                           Patch(color="#2166ac", label="gained under UNC")],
                  fontsize=7, frameon=False, loc="upper right", ncol=2)
        fig.savefig(outpath, dpi=150)
        plt.close(fig)
    except Exception as e:
        print(f"[warn] raster plot failed: {e}")


# =============================================================================
# Domain-set dynamics (Jaccard + lost/gained spatial signature)
# =============================================================================
def _runlen_and_adj(dom, retained, starts, ends):
    """Within one chromosome (sorted): per-bin run length in the `dom` track, and
    whether each bin is genomically adjacent to a `retained` domain bin."""
    n = len(dom)
    adjL = np.r_[False, starts[1:] == ends[:-1]]          # bin i abuts bin i-1
    new = dom & ~(np.r_[False, dom[:-1]] & adjL)
    runid = np.where(dom, np.cumsum(new), 0)
    lengths = np.bincount(runid)
    runlen = np.where(dom, lengths[runid], 0).astype(float)
    left_ret = np.r_[False, retained[:-1]] & adjL
    adjR = np.r_[starts[1:] == ends[:-1], False]          # bin i abuts bin i+1
    right_ret = np.r_[retained[1:], False] & adjR
    return runlen, (left_ret | right_ret)


def analyze_domain_dynamics(df, top_q):
    """For each cell line, compare the top-q domain SET between DMSO and UNC.

    Discriminates 'consolidation by attrition' (biologically real differential
    erasure) from a 'thresholding artifact':
      * Jaccard(DMSO domains, UNC domains): set stability. A PURE global-scaling
        artifact preserves rank order -> Jaccard=1 and no metric change; the
        observed metric change therefore REQUIRES a rank-structure change.
      * lost bins should sit in SHORT DMSO runs / be singletons (the dispersed
        inter-domain spread), retained bins in LONG runs (the cores).
      * gained bins should be ADJACENT to retained cores (they extend/merge
        domains) rather than seed new isolated ones.
    """
    d = df.sort_values(["chr", "bed_start"]).reset_index(drop=True)
    dom = {c: (rankdata(d[c].values) / len(d)) >= (1.0 - top_q) for c in COND_COLS}
    base_rank = {}
    if "baseline_b10" in d:
        base_rank["b10"] = rankdata(d["baseline_b10"].values) / len(d)
    if "baseline_parental" in d:
        base_rank["parental"] = rankdata(d["baseline_parental"].values) / len(d)

    pairs = {"b10": ("b10_none", "b10_unc"),
             "parental": ("parental_none", "parental_unc")}
    rows = []
    for cell, (cn, cu) in pairs.items():
        dn, du = dom[cn], dom[cu]
        union = (dn | du).sum()
        jac = (dn & du).sum() / union if union else np.nan
        lost = dn & ~du; gained = du & ~dn; retained = dn & du

        runlen = np.zeros(len(d)); adj_ret = np.zeros(len(d), bool)
        for _, g in d.groupby("chr", sort=False):
            idx = g.index.values
            rl, ar = _runlen_and_adj(dn[idx], retained[idx],
                                     g["bed_start"].values, g["bed_end"].values)
            runlen[idx] = rl; adj_ret[idx] = ar

        row = {"cell_line": cell, "jaccard": jac,
               "n_dmso_domains": int(dn.sum()), "n_unc_domains": int(du.sum()),
               "n_lost": int(lost.sum()), "n_gained": int(gained.sum()),
               "n_retained": int(retained.sum()),
               "lost_mean_dmso_runlen": float(runlen[lost].mean()) if lost.any() else np.nan,
               "retained_mean_dmso_runlen": float(runlen[retained].mean()) if retained.any() else np.nan,
               "lost_singleton_frac": float((runlen[lost] == 1).mean()) if lost.any() else np.nan,
               "retained_singleton_frac": float((runlen[retained] == 1).mean()) if retained.any() else np.nan,
               "gained_adj_to_core_frac": float(adj_ret[gained].mean()) if gained.any() else np.nan}
        if cell in base_rank:
            row["lost_mean_baseline_rank"] = float(base_rank[cell][lost].mean()) if lost.any() else np.nan
            row["retained_mean_baseline_rank"] = float(base_rank[cell][retained].mean()) if retained.any() else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


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

    # ---- domain-set dynamics: is UNC consolidation attrition or thresholding? ----
    if args.jaccard:
        dyn = analyze_domain_dynamics(df, args.top_q)
        dyn.to_csv(os.path.join(args.outdir, "domain_set_dynamics.tsv"),
                   sep="\t", index=False)
        print("[info] domain-set dynamics (DMSO vs UNC top-q domain set):")
        for _, r in dyn.iterrows():
            verdict = ("attrition (lost=dispersed singletons, gained extend cores)"
                       if (np.isfinite(r["lost_singleton_frac"])
                           and r["lost_singleton_frac"] > r["retained_singleton_frac"]
                           and np.isfinite(r["gained_adj_to_core_frac"])
                           and r["gained_adj_to_core_frac"] > 0.5)
                       else "spatially unstructured (caution: closer to thresholding)")
            print(f"    {r['cell_line']:8s} Jaccard={r['jaccard']:.3f}  "
                  f"lost={int(r['n_lost'])} (singleton {r['lost_singleton_frac']:.2f} "
                  f"vs retained {r['retained_singleton_frac']:.2f}; "
                  f"DMSO run-len {r['lost_mean_dmso_runlen']:.2f} vs "
                  f"{r['retained_mean_dmso_runlen']:.2f})  "
                  f"gained adj-to-core={r['gained_adj_to_core_frac']:.2f}")
            print(f"             -> {verdict}")

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
    window = tuple(args.raster_window) if args.raster_window else None
    wtag = f"_{window[0]:g}-{window[1]:g}Mb" if window else ""
    plot_domain_raster(df, str(chrom), args.top_q,
                       os.path.join(args.outdir, f"domain_raster_{chrom}{wtag}.png"),
                       window=window, diff_cellline=args.raster_cellline)

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
    p.add_argument("--raster-window", type=float, nargs=2, default=None,
                   metavar=("START_MB", "END_MB"),
                   help="zoom the raster to a genomic window in Mb (e.g. --raster-window "
                        "20 40); omit to show the whole chromosome")
    p.add_argument("--raster-cellline", choices=["b10", "parental"], default="b10",
                   help="cell line whose UNC lost/gained domain track is drawn under "
                        "the four condition rows")
    p.add_argument("--jaccard", action=argparse.BooleanOptionalAction, default=True,
                   help="report domain-set dynamics (Jaccard overlap of the top-q "
                        "domain set between DMSO and UNC, plus the lost/gained spatial "
                        "signature that separates real attrition from a thresholding "
                        "effect). Use --no-jaccard to skip.")
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