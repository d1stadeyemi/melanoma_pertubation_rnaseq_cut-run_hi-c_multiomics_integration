#!/usr/bin/env python
"""
=============================================================================
integrate_hic_cutrun_rnaseq.py
  Three-way integration: 3D genome (Hi-C PC1) x chromatin (CUT&RUN) x
  expression (DESeq2), at bin level.
=============================================================================

WHAT THIS ANSWERS
  For each histone mark: does the CHANGE in that mark track the CHANGE in
  expression, and is that association independent of the change in compartment
  (PC1)? The partial correlation is the point of the whole script -- a raw
  correlation between mark change and expression change can arise purely
  because both track the compartment switch, with no direct mark-expression
  relationship at all.

WHY THE STATISTICS ARE NOT JUST scipy.stats.spearmanr
  Two problems make a naive Spearman p-value badly anticonservative here.

  1. PSEUDOREPLICATION. Mark change and PC1 change are BIN-level quantities:
     every gene in a 250 kb bin shares the same value. Correlating at gene
     level inflates n by the number of genes per bin. Expression is therefore
     aggregated to a per-bin summary before anything is tested, and each bin
     contributes exactly one observation.

  2. SPATIAL AUTOCORRELATION. A/B compartments are megabase-scale domains, so
     adjacent 250 kb bins are strongly correlated -- PC1 autocorrelation is a
     well-documented property of Hi-C eigenvectors. Classical tests assume
     independent observations; under autocorrelation the variance of a
     correlation estimate is underestimated and p-values are far too small.
     This script therefore reports:
       * lag-1 autocorrelation and an EFFECTIVE SAMPLE SIZE (Clifford-Richardson
         style variance inflation), so the reader can see how much of the
         nominal n is real;
       * a MOVING BLOCK BOOTSTRAP confidence interval on rho, resampling
         contiguous runs of bins so within-block dependence is preserved;
       * a CIRCULAR BLOCK PERMUTATION p-value, shifting one variable relative
         to the other within each chromosome. This destroys the cross-variable
         association while leaving each variable's own autocorrelation exactly
         intact -- the correct null here.
     Block length is estimated from the autocorrelation decay of the data
     itself, and can be overridden.

  3. MULTIPLE TESTING across marks: Benjamini-Hochberg over the set of marks
     tested, for both the raw and partial correlations.

FRAMING
  Everything here is ASSOCIATION. Both conditions are measured after migration
  has already happened, so nothing in this script establishes direction of
  causation. The causal arm is the inhibitor perturbation.

INPUTS
  --hic       Hi-C table. Needs a coordinate column and a numeric delta-PC1
              column, both named on the command line.
  --cutrun    CUT&RUN signal table. Needs a coordinate column and one or more
              numeric mark-change columns, named on the command line.
  --deseq2    DESeq2 results (one row per gene) with log2FoldChange and padj.
  --gene_map  Optional gene -> bin coordinate table, if --deseq2 has no
              coordinate column of its own.

OUTPUT (into --outdir)
  integration_bin_table.txt      per-bin merged values actually analysed
  integration_stats.txt          per-mark rho, partial rho, CIs, p-values
  integration_scatter.pdf/png    per-mark scatter, coloured by delta PC1
  integration_diagnostics.txt    autocorrelation and effective-n report
=============================================================================
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import rankdata, spearmanr, norm

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 10, "axes.labelsize": 11,
    "axes.titlesize": 11, "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600,
    "figure.facecolor": "white",
})


# ---------------------------------------------------------------------------
# Autocorrelation diagnostics
# ---------------------------------------------------------------------------
def lag1_acf(x):
    """Lag-1 autocorrelation of a 1-D series (NaNs dropped pairwise)."""
    x = np.asarray(x, dtype=float)
    a, b = x[:-1], x[1:]
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return np.nan
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def acf_block_length(x, max_lag=40):
    """Block length = first lag at which the autocorrelation falls below 1/e.

    A block must be long enough to contain the dependence, so that resampled
    blocks are approximately independent of one another. 1/e is the standard
    'correlation length' convention; the value is reported and can be
    overridden with --block_len.
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = len(x)
    if n < 10:
        return 2
    x = x - x.mean()
    denom = np.dot(x, x)
    if denom == 0:
        return 2
    for lag in range(1, min(max_lag, n - 2)):
        r = np.dot(x[:-lag], x[lag:]) / denom
        if r < 1.0 / np.e:
            return max(2, lag)
    return max(2, min(max_lag, n // 10))


def effective_n(x, y):
    """Effective sample size under autocorrelation (Clifford-Richardson idea).

    n_eff = n / (1 + 2*sum_k rho_x(k)*rho_y(k)). When both series are
    autocorrelated the variance of the correlation estimate is inflated by
    roughly this factor, so the honest n is much smaller than the bin count.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok] - np.mean(x[ok]), y[ok] - np.mean(y[ok])
    n = len(x)
    if n < 10:
        return float(n)
    dx, dy = np.dot(x, x), np.dot(y, y)
    if dx == 0 or dy == 0:
        return float(n)
    s = 0.0
    for k in range(1, min(50, n // 4)):
        s += (np.dot(x[:-k], x[k:]) / dx) * (np.dot(y[:-k], y[k:]) / dy)
    return float(np.clip(n / max(1e-9, 1.0 + 2.0 * s), 2, n))


# ---------------------------------------------------------------------------
# Correlation with autocorrelation-aware inference
# ---------------------------------------------------------------------------
def partial_spearman(x, y, z):
    """Spearman partial correlation of x and y given z.

    Rank-transform all three, then correlate the residuals of x~z and y~z.
    This is what separates 'the mark tracks expression' from 'the mark and
    expression both track the compartment'.
    """
    x, y, z = map(lambda v: np.asarray(v, float), (x, y, z))
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if ok.sum() < 10:
        return np.nan
    rx, ry, rz = (rankdata(v[ok]) for v in (x, y, z))
    Z = np.column_stack([np.ones(ok.sum()), rz])
    ex = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
    ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    if np.std(ex) == 0 or np.std(ey) == 0:
        return np.nan
    return float(np.corrcoef(ex, ey)[0, 1])


def _blocks(idx, block_len, rng):
    """Moving-block resample of an index array to (about) its original length."""
    n = len(idx)
    starts = rng.integers(0, max(1, n - block_len + 1), size=int(np.ceil(n / block_len)))
    out = np.concatenate([idx[s:s + block_len] for s in starts])
    return out[:n]


def block_bootstrap_ci(df, xcol, ycol, zcol, chrom_col, block_len,
                       n_boot=2000, partial=False, seed=0, alpha=0.05):
    """Percentile CI for (partial) Spearman via moving block bootstrap.

    Blocks are drawn within chromosomes, so a resampled dataset preserves the
    local dependence structure that makes naive CIs too narrow.
    """
    rng = np.random.default_rng(seed)
    groups = [g.index.to_numpy() for _, g in df.groupby(chrom_col, sort=False)]
    stats = []
    for _ in range(n_boot):
        take = np.concatenate([_blocks(g, block_len, rng) for g in groups])
        sub = df.loc[take]
        if partial:
            r = partial_spearman(sub[xcol], sub[ycol], sub[zcol])
        else:
            xr, yr = sub[xcol].to_numpy(), sub[ycol].to_numpy()
            ok = np.isfinite(xr) & np.isfinite(yr)
            r = spearmanr(xr[ok], yr[ok]).statistic if ok.sum() > 5 else np.nan
        if np.isfinite(r):
            stats.append(r)
    if len(stats) < 50:
        return (np.nan, np.nan, np.nan)
    s = np.asarray(stats)
    return (float(np.percentile(s, 100 * alpha / 2)),
            float(np.percentile(s, 100 * (1 - alpha / 2))),
            float(np.std(s)))


def block_permutation_p(df, xcol, ycol, zcol, chrom_col, block_len,
                        n_perm=2000, partial=False, seed=0):
    """Two-sided p-value by circular block shift within each chromosome.

    Shifting y circularly relative to x leaves y's own autocorrelation exactly
    intact while destroying any real x-y association. This is the appropriate
    null for two autocorrelated genomic tracks: a plain permutation would
    destroy the autocorrelation too and give a null that is far too tight.
    """
    rng = np.random.default_rng(seed)
    if partial:
        obs = partial_spearman(df[xcol], df[ycol], df[zcol])
    else:
        obs = spearmanr(df[xcol], df[ycol], nan_policy="omit").statistic
    if not np.isfinite(obs):
        return np.nan, obs
    groups = [g.index.to_numpy() for _, g in df.groupby(chrom_col, sort=False)]
    y = df[ycol].to_numpy()
    ge = 0
    for _ in range(n_perm):
        ysh = y.copy()
        for g in groups:
            pos = df.index.get_indexer(g)
            k = int(rng.integers(1, max(2, len(g))))
            ysh[pos] = np.roll(y[pos], k)
        tmp = df.assign(**{"_y_": ysh})
        if partial:
            r = partial_spearman(tmp[xcol], tmp["_y_"], tmp[zcol])
        else:
            r = spearmanr(tmp[xcol], tmp["_y_"], nan_policy="omit").statistic
        if np.isfinite(r) and abs(r) >= abs(obs):
            ge += 1
    return (ge + 1) / (n_perm + 1), obs


def bh(pvals):
    p = np.asarray(pvals, float)
    ok = np.isfinite(p)
    out = np.full(p.shape, np.nan)
    q = p[ok]
    order = np.argsort(q)
    ranked = q[order]
    m = len(ranked)
    adj = ranked * m / (np.arange(m) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    res = np.empty(m); res[order] = np.clip(adj, 0, 1)
    out[ok] = res
    return out


# ---------------------------------------------------------------------------
# IO / assembly
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Three-way Hi-C x CUT&RUN x RNA-seq integration at bin level.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python integrate_hic_cutrun_rnaseq.py \\
    --hic data/processed/hic_pc1_table.txt \\
    --hic_coord_col Coord \\
    --hic_delta_col "Bottom10-Parental_DMSO" \\
    --hic_category_col "Bottom10-Parental_DMSO_category" \\
    --cutrun data/processed/cut_run_signal_table_250kb.eps_fixed.txt \\
    --cutrun_coord_col extra_1 \\
    --mark_cols log2ratio_p-b10_h3k27me3_mean log2ratio_p-b10_h3k27ac_mean \\
                log2ratio_p-b10_h3k4me1_mean log2ratio_p-b10_h3k9me3_mean \\
    --deseq2 results/deseq2/b10_vs_parental_results.txt \\
    --gene_map data/processed/hg19_gene_list_binned.txt \\
    --outdir results/exploratory/integration
""")
    p.add_argument("--hic", required=True)
    p.add_argument("--hic_coord_col", default="Coord")
    p.add_argument("--hic_delta_col", required=True,
                   help="Numeric delta-PC1 column, e.g. 'Bottom10-Parental_DMSO'.")
    p.add_argument("--hic_category_col", default=None,
                   help="Optional switch-category column, used only for plot symbols.")
    p.add_argument("--hic_chr_col", default="Chr")

    p.add_argument("--cutrun", required=True)
    p.add_argument("--cutrun_coord_col", default=None,
                   help="Default: auto-detect 'coordinate' then 'extra_1'.")
    p.add_argument("--mark_cols", nargs="+", required=True,
                   help="One or more numeric mark-CHANGE columns from --cutrun.")
    p.add_argument("--mark_names", nargs="+", default=None,
                   help="Display names, same order as --mark_cols.")

    p.add_argument("--deseq2", required=True)
    p.add_argument("--deseq2_lfc_col", default="log2FoldChange")
    p.add_argument("--deseq2_padj_col", default="padj")
    p.add_argument("--deseq2_gene_col", default=None,
                   help="Gene id/symbol column; default: first column.")
    p.add_argument("--deseq2_coord_col", default=None,
                   help="Bin coordinate column IF already present in --deseq2.")
    p.add_argument("--gene_map", default=None,
                   help="Gene -> bin table, needed if --deseq2 lacks coordinates.")
    p.add_argument("--gene_map_gene_col", default="gene_name")
    p.add_argument("--gene_map_coord_col", default="coordinate")

    p.add_argument("--expr_agg", default="mean", choices=["mean", "median"],
                   help="How gene LFCs are summarised to one value per bin.")
    p.add_argument("--min_genes", type=int, default=1,
                   help="Drop bins with fewer than this many genes (default 1). "
                        "Raising it trades bins for a more stable per-bin mean.")
    p.add_argument("--sig_only", action="store_true",
                   help="Use only genes passing --padj_cut when aggregating.")
    p.add_argument("--padj_cut", type=float, default=0.05)

    p.add_argument("--restrict_to", nargs="+", default=None,
                   help="Optional values of --hic_category_col to keep "
                        "(e.g. Switch_AtoB Switch_BtoA).")
    p.add_argument("--block_len", type=int, default=None,
                   help="Block length in bins for bootstrap/permutation. "
                        "Default: estimated from the delta-PC1 autocorrelation.")
    p.add_argument("--n_boot", type=int, default=2000)
    p.add_argument("--n_perm", type=int, default=2000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--outdir", default="results/exploratory/integration")
    return p.parse_args()


def main():
    a = parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    names = a.mark_names or a.mark_cols
    if len(names) != len(a.mark_cols):
        sys.exit("[ERROR] --mark_names must match --mark_cols in length.")

    # ---- Hi-C ------------------------------------------------------------
    print("[INFO] Loading Hi-C table...")
    hic = pd.read_csv(a.hic, sep="\t", low_memory=False)
    for c in (a.hic_coord_col, a.hic_delta_col):
        if c not in hic.columns:
            sys.exit(f"[ERROR] '{c}' not in {a.hic}. Columns: {list(hic.columns)[:25]}")
    hic[a.hic_coord_col] = hic[a.hic_coord_col].astype(str).str.strip()
    hic["_dpc1_"] = pd.to_numeric(
        hic[a.hic_delta_col].astype(str).str.replace(",", "", regex=False),
        errors="coerce")
    if a.hic_chr_col in hic.columns:
        hic["_chr_"] = hic[a.hic_chr_col].astype(str)
    else:
        hic["_chr_"] = hic[a.hic_coord_col].str.split(":").str[0]
    keep = [a.hic_coord_col, "_dpc1_", "_chr_"]
    if a.hic_category_col and a.hic_category_col in hic.columns:
        hic["_cat_"] = hic[a.hic_category_col].astype(str).str.strip()
        keep.append("_cat_")
    hic = hic[keep].rename(columns={a.hic_coord_col: "coord"})
    print(f"       {len(hic):,} Hi-C bins")

    # ---- CUT&RUN ---------------------------------------------------------
    print("[INFO] Loading CUT&RUN table...")
    cr = pd.read_csv(a.cutrun, sep="\t", low_memory=False)
    ck = a.cutrun_coord_col
    if ck is None:
        for cand in ("coordinate", "extra_1"):
            if cand in cr.columns:
                ck = cand; break
    if ck is None or ck not in cr.columns:
        sys.exit("[ERROR] no CUT&RUN coordinate column; pass --cutrun_coord_col.")
    missing = [c for c in a.mark_cols if c not in cr.columns]
    if missing:
        sys.exit(f"[ERROR] mark columns not found: {missing}")
    cr[ck] = cr[ck].astype(str).str.strip()
    cr = cr[[ck] + a.mark_cols].rename(columns={ck: "coord"})
    for c in a.mark_cols:
        cr[c] = pd.to_numeric(cr[c], errors="coerce")
    cr = cr.drop_duplicates(subset="coord")
    print(f"       {len(cr):,} CUT&RUN bins, marks: {a.mark_cols}")

    # ---- DESeq2 ----------------------------------------------------------
    print("[INFO] Loading DESeq2 results...")
    de = pd.read_csv(a.deseq2, sep="\t", low_memory=False)
    if a.deseq2_lfc_col not in de.columns:
        sys.exit(f"[ERROR] '{a.deseq2_lfc_col}' not in {a.deseq2}.")
    gcol = a.deseq2_gene_col or de.columns[0]
    de["_lfc_"] = pd.to_numeric(de[a.deseq2_lfc_col], errors="coerce")
    if a.deseq2_padj_col in de.columns:
        de["_padj_"] = pd.to_numeric(de[a.deseq2_padj_col], errors="coerce")
    else:
        de["_padj_"] = np.nan

    if a.deseq2_coord_col and a.deseq2_coord_col in de.columns:
        de["coord"] = de[a.deseq2_coord_col].astype(str).str.strip()
    elif a.gene_map:
        gm = pd.read_csv(a.gene_map, sep="\t", low_memory=False)
        for c in (a.gene_map_gene_col, a.gene_map_coord_col):
            if c not in gm.columns:
                sys.exit(f"[ERROR] '{c}' not in {a.gene_map}.")
        gm = gm[[a.gene_map_gene_col, a.gene_map_coord_col]].copy()
        gm.columns = ["_g_", "coord"]
        gm["_g_"] = gm["_g_"].astype(str).str.strip()
        gm["coord"] = gm["coord"].astype(str).str.strip()
        de["_g_"] = de[gcol].astype(str).str.strip()
        before = len(de)
        de = de.merge(gm.drop_duplicates("_g_"), on="_g_", how="inner")
        print(f"       {len(de):,}/{before:,} genes mapped to a bin")
    else:
        sys.exit("[ERROR] DESeq2 file has no coordinate column; supply "
                 "--deseq2_coord_col or --gene_map.")

    if a.sig_only:
        n0 = len(de)
        de = de[de["_padj_"] < a.padj_cut]
        print(f"       {len(de):,}/{n0:,} genes kept at padj < {a.padj_cut}")

    agg = (de.dropna(subset=["_lfc_"])
             .groupby("coord")
             .agg(expr=("_lfc_", a.expr_agg),
                  expr_sd=("_lfc_", "std"),
                  n_genes=("_lfc_", "size"))
             .reset_index())
    print(f"       {len(agg):,} bins carry at least one gene")

    # ---- Merge -----------------------------------------------------------
    df = (hic.merge(cr, on="coord", how="inner")
             .merge(agg, on="coord", how="inner"))
    df = df[df["n_genes"] >= a.min_genes]
    if a.restrict_to and "_cat_" in df.columns:
        df = df[df["_cat_"].isin(a.restrict_to)]
        print(f"[INFO] restricted to categories {a.restrict_to}")
    df = df.dropna(subset=["_dpc1_", "expr"])

    # Order by genomic position -- essential, the block methods assume the
    # rows are in chromosome order.
    df["_start_"] = (df["coord"].str.split(":").str[1]
                                .str.split("-").str[0].astype(float))
    df = df.sort_values(["_chr_", "_start_"]).reset_index(drop=True)
    if len(df) < 30:
        sys.exit(f"[ERROR] only {len(df)} bins survived the merge; check that the "
                 f"coordinate formats match across the three files.")
    print(f"[INFO] {len(df):,} bins with Hi-C + CUT&RUN + expression")

    # ---- Diagnostics -----------------------------------------------------
    blk = a.block_len or acf_block_length(df["_dpc1_"].to_numpy())
    diag = [f"bins analysed: {len(df)}",
            f"block length used: {blk} bins ({blk*250} kb, if 250 kb bins)",
            f"lag-1 autocorrelation, delta PC1: {lag1_acf(df['_dpc1_']):.3f}",
            f"lag-1 autocorrelation, bin expression: {lag1_acf(df['expr']):.3f}"]
    for c, nm in zip(a.mark_cols, names):
        diag.append(f"lag-1 autocorrelation, {nm}: {lag1_acf(df[c]):.3f}")
    print("\n".join("  " + d for d in diag))

    # ---- Per-mark statistics --------------------------------------------
    rows = []
    for col, nm in zip(a.mark_cols, names):
        sub = df.dropna(subset=[col]).reset_index(drop=True)
        n_eff = effective_n(sub[col], sub["expr"])
        rho_naive = spearmanr(sub[col], sub["expr"]).statistic
        p_naive = spearmanr(sub[col], sub["expr"]).pvalue
        lo, hi, se = block_bootstrap_ci(sub, col, "expr", "_dpc1_", "_chr_",
                                        blk, a.n_boot, False, a.seed)
        p_blk, _ = block_permutation_p(sub, col, "expr", "_dpc1_", "_chr_",
                                       blk, a.n_perm, False, a.seed)
        prho = partial_spearman(sub[col], sub["expr"], sub["_dpc1_"])
        plo, phi, pse = block_bootstrap_ci(sub, col, "expr", "_dpc1_", "_chr_",
                                           blk, a.n_boot, True, a.seed)
        p_pblk, _ = block_permutation_p(sub, col, "expr", "_dpc1_", "_chr_",
                                        blk, a.n_perm, True, a.seed)
        rows.append(dict(
            mark=nm, column=col, n_bins=len(sub), n_eff=round(n_eff, 1),
            spearman_rho=rho_naive, p_naive_INVALID=p_naive,
            rho_CI95_lo=lo, rho_CI95_hi=hi, rho_block_se=se, p_block_perm=p_blk,
            partial_rho_given_dPC1=prho, partial_CI95_lo=plo,
            partial_CI95_hi=phi, p_block_perm_partial=p_pblk))
        print(f"  {nm}: rho={rho_naive:+.3f} [{lo:+.3f},{hi:+.3f}] "
              f"p_block={p_blk:.4g} | partial={prho:+.3f} "
              f"[{plo:+.3f},{phi:+.3f}] p={p_pblk:.4g} | n={len(sub)} n_eff={n_eff:.0f}")

    st = pd.DataFrame(rows)
    st["q_block_perm"] = bh(st["p_block_perm"])
    st["q_block_perm_partial"] = bh(st["p_block_perm_partial"])
    st.to_csv(os.path.join(a.outdir, "integration_stats.txt"), sep="\t", index=False)

    # ---- Figure ----------------------------------------------------------
    n = len(a.mark_cols)
    nc = min(2, n); nr = (n + nc - 1) // nc
    fig, axes = plt.subplots(nr, nc, figsize=(5.6 * nc, 4.9 * nr))
    axes = np.atleast_1d(axes).flatten()
    vmax = float(np.nanpercentile(np.abs(df["_dpc1_"]), 98)) or 1.0
    sc = None
    for ax, (col, nm) in zip(axes, zip(a.mark_cols, names)):
        sub = df.dropna(subset=[col])
        sc = ax.scatter(sub[col], sub["expr"], c=sub["_dpc1_"], cmap="RdBu_r",
                        vmin=-vmax, vmax=vmax, s=16, alpha=0.8,
                        edgecolors="none", rasterized=True)
        ax.axhline(0, color="#BBBBBB", lw=0.6); ax.axvline(0, color="#BBBBBB", lw=0.6)
        r = st.loc[st["mark"] == nm].iloc[0]
        ax.set_title(f"{nm}\n"
                     f"rho={r.spearman_rho:+.2f} "
                     f"[{r.rho_CI95_lo:+.2f}, {r.rho_CI95_hi:+.2f}], "
                     f"p={r.p_block_perm:.3g}\n"
                     f"partial (|dPC1) = {r.partial_rho_given_dPC1:+.2f}, "
                     f"p={r.p_block_perm_partial:.3g}", fontsize=9)
        ax.set_xlabel(f"{nm} change  log2(B10/Parental)")
        ax.set_ylabel(f"bin {a.expr_agg} expression change (log2FC)")
    for ax in axes[n:]:
        ax.set_visible(False)
    if sc is not None:
        cb = fig.colorbar(sc, ax=axes[:n].tolist(), fraction=0.03, pad=0.02)
        cb.set_label("change in PC1 (B10 - Parental)")
    fig.suptitle("Chromatin change vs expression change, coloured by compartment change\n"
                 f"one dot per 250 kb bin; block-resampled statistics "
                 f"(block = {blk} bins); association, not causation",
                 fontweight="bold", fontsize=11)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(a.outdir, f"integration_scatter.{ext}"),
                    bbox_inches="tight", pad_inches=0.06, dpi=600)
    plt.close(fig)

    out = df.rename(columns={"_dpc1_": "delta_PC1", "_chr_": "chr"})
    out.drop(columns=["_start_"], errors="ignore").to_csv(
        os.path.join(a.outdir, "integration_bin_table.txt"), sep="\t", index=False)
    with open(os.path.join(a.outdir, "integration_diagnostics.txt"), "w") as f:
        f.write("\n".join(diag) + "\n")

    print(f"\n[INFO] Outputs -> {a.outdir}")
    print("[NOTE] 'p_naive_INVALID' is the ordinary Spearman p-value. It is "
          "reported only to show how much smaller it is than the block-"
          "permutation p-value; do NOT quote it. Use p_block_perm / "
          "q_block_perm, and read partial_rho_given_dPC1 to see whether the "
          "mark-expression link survives once compartment change is "
          "accounted for.")


if __name__ == "__main__":
    main()