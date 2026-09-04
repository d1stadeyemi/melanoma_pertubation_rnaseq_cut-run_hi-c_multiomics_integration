#!/usr/bin/env python
"""
=============================================================================
integrate_hic_cutrun_rnaseq_2.py
  Compartment-STRATIFIED mark-vs-expression correlation.
  Runs immediately after integrate_hic_cutrun_rnaseq_1.py and consumes its
  per-bin output, so the bins analysed are identical by construction.
=============================================================================

WHY THIS EXISTS
  integrate_hic_cutrun_rnaseq_1.py asks one question genome-wide: does mark
  change track expression change? A single pooled correlation can hide two
  different biologies, because a bin that starts in A and a bin that starts in
  B are not the same kind of locus. The same pooling problem already showed up
  in the compartment-resolved Random Forest analysis, where features that
  ranked first in the pooled model dropped to last once the model was forced to
  discriminate WITHIN a compartment -- the pooled ranking was partly reporting
  which compartment a locus lived in, not what predicted its behaviour.

  This script therefore splits bins by their STARTING (Parental) compartment
  and reports the correlation separately in each, alongside the pooled value.

WHAT IT ADDS BEYOND TWO SEPARATE NUMBERS
  Reporting rho_A and rho_B with their confidence intervals does NOT test
  whether they differ: two intervals can overlap while the difference is real,
  and can both exclude zero while being indistinguishable from each other. The
  script therefore bootstraps the DIFFERENCE rho_A - rho_B directly, resampling
  both strata within the same iteration, and reports a CI and p-value for the
  contrast itself. Only that number licenses a claim that the coupling is
  stronger in one compartment.

  All resampling is block-based (contiguous runs of bins within chromosomes),
  because adjacent 250 kb bins are not independent -- same reasoning as the
  upstream script.

STRATIFICATION -- two modes
  (a) BY COLUMN (default). --stratify_col names a categorical column already in
      the bin table. The natural choice is '_cat_', the switch category carried
      over from the upstream run, which splits the switching bins into their
      two directions (e.g. Switch AtoB vs Switch BtoA). Needs no extra input.
  (b) BY BASELINE COMPARTMENT. --pc1_col names the CONTROL-condition PC1 column
      in --hic; bins are split by its sign into A-resident and B-resident, i.e.
      by where each bin STARTED rather than by how it moved. --pc1_min drops
      bins near the A/B boundary, whose assignment is arbitrary.
  Exactly two strata are compared; use --strata to choose which two if the
  column holds more.

INPUTS
  --bin_table   integration_bin_table.txt written by integrate_hic_cutrun_rnaseq_1.py
  --hic, --pc1_col   only for mode (b)

OUTPUT (into --outdir)
  stratified_stats.txt        per mark x stratum, plus the A-vs-B difference
  stratified_scatter.pdf/png  one row per mark, A and B panels side by side
  stratified_diagnostics.txt  per-stratum autocorrelation and effective n
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
from scipy.stats import rankdata, spearmanr

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 10, "axes.labelsize": 11,
    "axes.titlesize": 11, "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600,
    "figure.facecolor": "white",
})

A_COL, B_COL = "#B54B26", "#1F6187"    # warm = toward A, cool = toward B


def _dir_colour(level):
    """Colour by DESTINATION compartment, taken from the last character of the
    level name ('Switch AtoB', 'Shift in A towardB' -> B). Colouring by panel
    position instead would let the same category change colour between runs,
    because the level order comes from whichever bin appears first in the file.
    """
    return B_COL if str(level).strip()[-1:].upper() == "B" else A_COL


def read_table(path, **kw):
    with open(path) as fh:
        head = fh.readline()
    sep = "," if head.count(",") > head.count("\t") else "\t"
    return pd.read_csv(path, sep=sep, low_memory=False, **kw)


# ---------------------------------------------------------------------------
# Autocorrelation-aware machinery (mirrors the upstream script)
# ---------------------------------------------------------------------------
def lag1_acf(x):
    x = np.asarray(x, float)
    a, b = x[:-1], x[1:]
    ok = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[ok], b[ok])[0, 1]) if ok.sum() >= 3 else np.nan


def acf_block_length(x, max_lag=40):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    n = len(x)
    if n < 10:
        return 2
    x = x - x.mean()
    denom = np.dot(x, x)
    if denom == 0:
        return 2
    for lag in range(1, min(max_lag, n - 2)):
        if np.dot(x[:-lag], x[lag:]) / denom < 1.0 / np.e:
            return max(2, lag)
    return max(2, min(max_lag, n // 10))


def effective_n(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok] - np.mean(x[ok]), y[ok] - np.mean(y[ok])
    n = len(x)
    if n < 10:
        return float(n)
    dx, dy = np.dot(x, x), np.dot(y, y)
    if dx == 0 or dy == 0:
        return float(n)
    s = sum((np.dot(x[:-k], x[k:]) / dx) * (np.dot(y[:-k], y[k:]) / dy)
            for k in range(1, min(50, n // 4)))
    return float(np.clip(n / max(1e-9, 1 + 2 * s), 2, n))


def partial_spearman(x, y, z):
    """Spearman of x,y after removing z, on ranks."""
    x, y, z = (np.asarray(v, float) for v in (x, y, z))
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


def _rho(sub, xcol, ycol, zcol, partial):
    if partial:
        return partial_spearman(sub[xcol], sub[ycol], sub[zcol])
    x, y = sub[xcol].to_numpy(), sub[ycol].to_numpy()
    ok = np.isfinite(x) & np.isfinite(y)
    return spearmanr(x[ok], y[ok]).statistic if ok.sum() > 5 else np.nan


def _block_idx(idx, blk, rng):
    n = len(idx)
    starts = rng.integers(0, max(1, n - blk + 1), size=int(np.ceil(n / blk)))
    return np.concatenate([idx[s:s + blk] for s in starts])[:n]


def block_resample(df, blk, rng, chrom_col="chr"):
    """One moving-block bootstrap replicate, blocks drawn within chromosome."""
    take = np.concatenate([_block_idx(g.index.to_numpy(), blk, rng)
                           for _, g in df.groupby(chrom_col, sort=False)])
    return df.loc[take]


def stratified_bootstrap(dfA, dfB, xcol, ycol, zcol, blk, n_boot, partial, seed):
    """Bootstrap rho_A, rho_B and their DIFFERENCE in the same iterations.

    Resampling both strata together per iteration is what makes the difference
    distribution valid: it propagates the uncertainty of both estimates into
    the contrast rather than comparing two intervals by eye.
    """
    rng = np.random.default_rng(seed)
    rA, rB, rD = [], [], []
    for _ in range(n_boot):
        a = _rho(block_resample(dfA, blk, rng), xcol, ycol, zcol, partial)
        b = _rho(block_resample(dfB, blk, rng), xcol, ycol, zcol, partial)
        if np.isfinite(a) and np.isfinite(b):
            rA.append(a); rB.append(b); rD.append(a - b)
    if len(rD) < 50:
        return {k: np.nan for k in
                ("A_lo","A_hi","B_lo","B_hi","D_lo","D_hi","p_diff")}
    rA, rB, rD = map(np.asarray, (rA, rB, rD))
    # Two-sided bootstrap p for the difference: how often the resampled
    # difference sits on the opposite side of zero from the observed one.
    frac = min((rD <= 0).mean(), (rD >= 0).mean())
    return dict(
        A_lo=float(np.percentile(rA, 2.5)),  A_hi=float(np.percentile(rA, 97.5)),
        B_lo=float(np.percentile(rB, 2.5)),  B_hi=float(np.percentile(rB, 97.5)),
        D_lo=float(np.percentile(rD, 2.5)),  D_hi=float(np.percentile(rD, 97.5)),
        p_diff=float(min(1.0, 2 * frac + 1.0 / (len(rD) + 1))))


def block_permutation_p(df, xcol, ycol, zcol, blk, n_perm, partial, seed,
                        chrom_col="chr"):
    """Circular block shift of y within chromosome; preserves autocorrelation."""
    rng = np.random.default_rng(seed)
    obs = _rho(df, xcol, ycol, zcol, partial)
    if not np.isfinite(obs):
        return np.nan
    groups = [g.index.to_numpy() for _, g in df.groupby(chrom_col, sort=False)]
    y = df[ycol].to_numpy()
    ge = 0
    for _ in range(n_perm):
        ysh = y.copy()
        for g in groups:
            pos = df.index.get_indexer(g)
            ysh[pos] = np.roll(y[pos], int(rng.integers(1, max(2, len(g)))))
        r = _rho(df.assign(**{ycol: ysh}), xcol, ycol, zcol, partial)
        if np.isfinite(r) and abs(r) >= abs(obs):
            ge += 1
    return (ge + 1) / (n_perm + 1)


def cliffs_delta(x, y):
    """Cliff's delta of x vs y. Positive => x tends to EXCEED y.

    delta = (#(x>y) - #(x<y)) / (nx * ny). The denominator counts only CROSS
    pairs -- every x against every y -- never x-vs-x or y-vs-y. Dividing by all
    pairwise comparisons instead caps the statistic near 0.5 and is the usual
    way this gets computed wrongly by hand. Range is -1 to +1; +1 means every x
    exceeds every y.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    x = x[np.isfinite(x)]; y = np.sort(y[np.isfinite(y)])
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return np.nan
    lt = np.searchsorted(y, x, side="left")        # y strictly below each x
    gt = ny - np.searchsorted(y, x, side="right")  # y strictly above each x
    return float((int(lt.sum()) - int(gt.sum())) / (nx * ny))


def delta_magnitude(d):
    """Conventional Cliff's delta labels (Romano et al.)."""
    if not np.isfinite(d):
        return "NA"
    ad = abs(d)
    return ("negligible" if ad < 0.147 else "small" if ad < 0.330
            else "medium" if ad < 0.474 else "large")


def group_delta_perm(sub, xcol, ycol, mask, n_perm, seed, chrom_col="chr"):
    """Cliff's delta for mask vs rest, with a circular-block-shift p-value.

    The GROUPS are defined by x (the mark change) and are held fixed; only the
    expression values are shifted along the chromosome. That keeps expression's
    own autocorrelation intact while breaking its correspondence to the groups
    -- the same null used for the correlations, so every p in this script comes
    from one framework rather than mixing a rank-sum test in alongside.
    """
    y = sub[ycol].to_numpy()
    m = np.asarray(mask, bool)
    if m.sum() < 5 or (~m).sum() < 5:
        return np.nan, np.nan
    obs = cliffs_delta(y[m], y[~m])
    if not np.isfinite(obs):
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    groups = [sub.index.get_indexer(g.index) for _, g in sub.groupby(chrom_col, sort=False)]
    ge = 0
    for _ in range(n_perm):
        ysh = y.copy()
        for pos in groups:
            ysh[pos] = np.roll(y[pos], int(rng.integers(1, max(2, len(pos)))))
        d = cliffs_delta(ysh[m], ysh[~m])
        if np.isfinite(d) and abs(d) >= abs(obs):
            ge += 1
    return obs, (ge + 1) / (n_perm + 1)


def bh(p):
    p = np.asarray(p, float); ok = np.isfinite(p)
    out = np.full(p.shape, np.nan)
    q = p[ok]; order = np.argsort(q); m = len(q)
    adj = np.minimum.accumulate((q[order] * m / (np.arange(m) + 1))[::-1])[::-1]
    res = np.empty(m); res[order] = np.clip(adj, 0, 1); out[ok] = res
    return out


# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Compartment-stratified mark-vs-expression correlation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python integrate_hic_cutrun_rnaseq_2.py \\
    --bin_table results/.../integration_bin_table.txt \\
    --mark_cols log2ratio_p-b10_h3k27me3_mean log2ratio_p-b10_h3k27ac_mean \\
                log2ratio_p-b10_h3k4me1_mean log2ratio_p-b10_h3k9me3_mean \\
    --mark_names H3K27me3 H3K27ac H3K4me1 H3K9me3 \\
    --outdir results/.../stratified
""")
    p.add_argument("--bin_table", required=True,
                   help="integration_bin_table.txt from integrate_hic_cutrun_rnaseq.py")
    p.add_argument("--stratify_col", default="_cat_",
                   help="Categorical column IN THE BIN TABLE to split on "
                        "(default '_cat_', the switch category from the "
                        "upstream run). Ignored when --pc1_col is given.")
    p.add_argument("--strata", nargs=2, default=None,
                   help="The two levels to compare, if the column holds more. "
                        "Order matters: the difference is level1 - level2.")
    p.add_argument("--hic", default=None, help="Hi-C table (only with --pc1_col).")
    p.add_argument("--hic_coord_col", default="Coord")
    p.add_argument("--pc1_col", default=None,
                   help="Switch to baseline-compartment mode: split by the sign "
                        "of this CONTROL-condition PC1 column in --hic, i.e. by "
                        "where each bin started rather than how it moved.")
    p.add_argument("--pc1_min", type=float, default=0.0,
                   help="Drop bins with |baseline PC1| <= this as compartment-"
                        "ambiguous (default 0 = keep all).")
    p.add_argument("--mark_cols", nargs="+", required=True)
    p.add_argument("--mark_names", nargs="+", default=None)
    p.add_argument("--expr_col", default="expr")
    p.add_argument("--dpc1_col", default="delta_PC1")
    p.add_argument("--block_len", type=int, default=None)
    p.add_argument("--n_boot", type=int, default=2000)
    p.add_argument("--n_perm", type=int, default=2000)
    p.add_argument("--partial", action="store_true",
                   help="Report PARTIAL correlations (conditioning on delta PC1) "
                        "instead of raw. Within a stratum the compartment is "
                        "already held roughly constant, so the raw correlation "
                        "is usually the quantity of interest; use this to check "
                        "that a stratum result is not driven by residual "
                        "within-stratum compartment movement.")
    p.add_argument("--effect_cut", type=float, default=0.5,
                   help="Threshold on the mark change (log2 units) defining the "
                        "'strong loss' (x < -cut) and 'strong gain' (x > +cut) "
                        "groups compared against all remaining bins of the same "
                        "stratum by Cliff's delta. Default 0.5 = a one-and-half-fold "
                        "change. Groups with fewer than 5 bins are reported as "
                        "NA rather than as an unstable number.")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def main():
    a = parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    names = a.mark_names or a.mark_cols
    if len(names) != len(a.mark_cols):
        sys.exit("[ERROR] --mark_names must match --mark_cols in length.")

    df = read_table(a.bin_table)
    for c in ["coord", "chr", a.expr_col, a.dpc1_col] + list(a.mark_cols):
        if c not in df.columns:
            sys.exit(f"[ERROR] '{c}' not in {a.bin_table}. "
                     f"Columns: {list(df.columns)}")
    print(f"[INFO] {len(df):,} bins from {a.bin_table}")

    if a.pc1_col:
        if not a.hic:
            sys.exit("[ERROR] --pc1_col requires --hic.")
        hic = read_table(a.hic)
        for c in (a.hic_coord_col, a.pc1_col):
            if c not in hic.columns:
                sys.exit(f"[ERROR] '{c}' not in {a.hic}. Columns: {list(hic.columns)[:25]}")
        hic = hic[[a.hic_coord_col, a.pc1_col]].copy()
        hic.columns = ["coord", "_pc1_"]
        hic["coord"] = hic["coord"].astype(str).str.strip()
        hic["_pc1_"] = pd.to_numeric(
            hic["_pc1_"].astype(str).str.replace(",", "", regex=False), errors="coerce")
        df["coord"] = df["coord"].astype(str).str.strip()
        df = df.merge(hic.drop_duplicates("coord"), on="coord", how="inner").dropna(subset=["_pc1_"])
        if a.pc1_min > 0:
            n0 = len(df); df = df[df["_pc1_"].abs() > a.pc1_min]
            print(f"[INFO] dropped {n0 - len(df):,} bins with |PC1| <= {a.pc1_min}")
        df["_stratum_"] = np.where(df["_pc1_"] > 0, "A-resident", "B-resident")
        strat_desc = f"sign of {a.pc1_col} (|PC1| > {a.pc1_min})"
    else:
        if a.stratify_col not in df.columns:
            sys.exit(f"[ERROR] --stratify_col '{a.stratify_col}' not in {a.bin_table}.\n"
                     f"        Columns: {list(df.columns)}\n"
                     f"        (Re-run the upstream script with --hic_category_col "
                     f"so the category is carried into the bin table, or pass "
                     f"--pc1_col to split by baseline compartment instead.)")
        df["_stratum_"] = df[a.stratify_col].astype(str).str.strip()
        strat_desc = f"column '{a.stratify_col}'"

    levels = list(pd.unique(df["_stratum_"].dropna()))
    if a.strata:
        missing = [x for x in a.strata if x not in levels]
        if missing:
            sys.exit(f"[ERROR] --strata {missing} not present. Available: {levels}")
        levels = list(a.strata)
        df = df[df["_stratum_"].isin(levels)]
    elif len(levels) != 2:
        sys.exit(f"[ERROR] need exactly 2 strata, found {len(levels)}: {levels}. "
                 f"Choose two with --strata.")
    L1, L2 = levels
    print(f"[INFO] stratifying on {strat_desc}: '{L1}' vs '{L2}'")

    blk = a.block_len or acf_block_length(df[a.dpc1_col].to_numpy())
    diag = [f"bins analysed: {len(df)}",
            f"stratifier: {strat_desc}  ({L1} vs {L2})",
            f"block length: {blk} bins",
            f"quantity: {'partial (given dPC1)' if a.partial else 'raw'} Spearman"]
    for st, sub in df.groupby("_stratum_"):
        diag.append(f"  {st}: n={len(sub)}, "
                    f"lag-1 ACF expression={lag1_acf(sub[a.expr_col]):.3f}")
    print("\n".join("  " + d for d in diag))

    dfA = df[df["_stratum_"] == L1].reset_index(drop=True)
    dfB = df[df["_stratum_"] == L2].reset_index(drop=True)
    if min(len(dfA), len(dfB)) < 100:
        print(f"[WARN] smallest stratum has {min(len(dfA), len(dfB))} bins. Block "
              f"resampling on a group this small gives wide intervals and the "
              f"difference test will have little power; treat the estimate as "
              f"indicative and report the CI, not the point value alone.")

    rows = []
    for col, nm in zip(a.mark_cols, names):
        A = dfA.dropna(subset=[col]).reset_index(drop=True)
        B = dfB.dropna(subset=[col]).reset_index(drop=True)
        P = df.dropna(subset=[col]).reset_index(drop=True)
        rA = _rho(A, col, a.expr_col, a.dpc1_col, a.partial)
        rB = _rho(B, col, a.expr_col, a.dpc1_col, a.partial)
        rP = _rho(P, col, a.expr_col, a.dpc1_col, a.partial)
        ci = stratified_bootstrap(A, B, col, a.expr_col, a.dpc1_col,
                                  blk, a.n_boot, a.partial, a.seed)
        pA = block_permutation_p(A, col, a.expr_col, a.dpc1_col, blk,
                                 a.n_perm, a.partial, a.seed)
        pB = block_permutation_p(B, col, a.expr_col, a.dpc1_col, blk,
                                 a.n_perm, a.partial, a.seed)
        # Group comparison: does STRONG loss (or gain) of the mark separate
        # bins on expression? This is the right summary when the scatter shows
        # a displaced group rather than a gradient, and unlike the correlation
        # it makes no monotonicity assumption.
        eff = {}
        for tag, sub_df in (("A", A), ("B", B)):
            for dirn, msk in (("loss", sub_df[col] < -a.effect_cut),
                              ("gain", sub_df[col] > a.effect_cut)):
                d, pd_ = group_delta_perm(sub_df, col, a.expr_col, msk,
                                          a.n_perm, a.seed)
                eff[f"delta_{dirn}_{tag}"] = d
                eff[f"p_delta_{dirn}_{tag}"] = pd_
                eff[f"n_{dirn}_{tag}"] = int(msk.sum())
                eff[f"mag_{dirn}_{tag}"] = delta_magnitude(d)

        rows.append(dict(
            mark=nm, column=col,
            rho_pooled=rP,
            stratum_1=L1, stratum_2=L2,
            n_A=len(A), rho_A=rA, A_CI_lo=ci["A_lo"], A_CI_hi=ci["A_hi"],
            p_A=pA, n_eff_A=round(effective_n(A[col], A[a.expr_col]), 1),
            n_B=len(B), rho_B=rB, B_CI_lo=ci["B_lo"], B_CI_hi=ci["B_hi"],
            p_B=pB, n_eff_B=round(effective_n(B[col], B[a.expr_col]), 1),
            rho_diff=(rA - rB),
            diff_CI_lo=ci["D_lo"], diff_CI_hi=ci["D_hi"], p_diff=ci["p_diff"],
            effect_cut=a.effect_cut, **eff))
        print(f"  {nm}: pooled={rP:+.3f} | A={rA:+.3f} [{ci['A_lo']:+.3f},{ci['A_hi']:+.3f}] "
              f"p={pA:.3g} (n={len(A)}) | B={rB:+.3f} [{ci['B_lo']:+.3f},{ci['B_hi']:+.3f}] "
              f"p={pB:.3g} (n={len(B)}) | diff={rA-rB:+.3f} "
              f"[{ci['D_lo']:+.3f},{ci['D_hi']:+.3f}] p={ci['p_diff']:.3g}")
        for tag, lab in (("A", L1), ("B", L2)):
            bits = [f"{d} n={eff[f'n_{d}_{tag}']} delta={eff[f'delta_{d}_{tag}']:+.3f} "
                    f"({eff[f'mag_{d}_{tag}']}) p={eff[f'p_delta_{d}_{tag}']:.3g}"
                    for d in ("loss", "gain")
                    if np.isfinite(eff[f"delta_{d}_{tag}"])]
            print(f"      [{lab}] " + ("  |  ".join(bits) if bits else
                  f"no group reached 5 bins at |x| > {a.effect_cut}"))

    st = pd.DataFrame(rows)
    for src, dst in (("p_A","q_A"), ("p_B","q_B"), ("p_diff","q_diff"),
                     ("p_delta_loss_A","q_delta_loss_A"), ("p_delta_gain_A","q_delta_gain_A"),
                     ("p_delta_loss_B","q_delta_loss_B"), ("p_delta_gain_B","q_delta_gain_B")):
        st[dst] = bh(st[src])
    st.to_csv(os.path.join(a.outdir, "stratified_stats.txt"), sep="\t", index=False)

    # ---- Figure: one ROW per mark, A and B side by side -------------------
    nm_n = len(a.mark_cols)
    fig, axes = plt.subplots(nm_n, 2, figsize=(9.2, 4.75 * nm_n),
                             constrained_layout=True, squeeze=False)
    # Extra vertical gap between rows: each row carries a caption underneath it
    # (the A-vs-B difference), which would otherwise collide with the titles of
    # the row below. API moved in matplotlib 3.6, so try both.
    try:
        fig.get_layout_engine().set(hspace=0.16)
    except AttributeError:
        fig.set_constrained_layout_pads(hspace=0.16)
    for i, (col, nm) in enumerate(zip(a.mark_cols, names)):
        r = st.iloc[i]
        for j, (sub, lab, colr, rho, lo, hi, pv, nn) in enumerate([
            (dfA, L1, _dir_colour(L1), r.rho_A, r.A_CI_lo, r.A_CI_hi, r.p_A, r.n_A),
            (dfB, L2, _dir_colour(L2), r.rho_B, r.B_CI_lo, r.B_CI_hi, r.p_B, r.n_B)]):
            ax = axes[i][j]
            s = sub.dropna(subset=[col])
            ss = float(np.clip(2200.0 / max(len(s), 1), 4, 22))
            al = float(np.clip(900.0 / max(len(s), 1), 0.30, 0.85))
            ax.scatter(s[col], s[a.expr_col], c=colr, s=ss, alpha=al,
                       linewidths=0, rasterized=True, zorder=2)
            ax.axhline(0, color="#BBBBBB", lw=0.6, zorder=1)
            ax.axvline(0, color="#BBBBBB", lw=0.6, zorder=1)
            # binned-median trend, same convention as the upstream figure
            x = s[col].to_numpy(); y = s[a.expr_col].to_numpy()
            qs = np.unique(np.nanquantile(x, np.linspace(0, 1, 11)))
            if len(qs) > 3:
                cx, cm, cl, ch = [], [], [], []
                for k in range(len(qs) - 1):
                    m = (x >= qs[k]) & (x <= qs[k+1] if k == len(qs)-2 else x < qs[k+1])
                    if m.sum() >= 10:
                        cx.append(np.median(x[m])); cm.append(np.median(y[m]))
                        cl.append(np.percentile(y[m], 25)); ch.append(np.percentile(y[m], 75))
                if len(cx) >= 3:
                    ax.fill_between(cx, cl, ch, color="#333333", alpha=0.12, lw=0, zorder=3)
                    ax.plot(cx, cm, color="#222222", lw=1.8, zorder=4)
                    ax.plot(cx, cm, "o", color="#222222", ms=3.2, zorder=5)
            pstr = ("< 5e-4" if pv <= 1.0 / (a.n_perm + 1) + 1e-12 else f"= {pv:.2g}")
            # One extra line: the strong-loss / strong-gain group comparison for
            # THIS panel. Only the direction that actually has >=5 bins is shown,
            # so the panel never carries a number computed on a handful of bins.
            tag = "A" if j == 0 else "B"
            extra = ""
            for dirn, sym in (("loss", "<"), ("gain", ">")):
                dv = r[f"delta_{dirn}_{tag}"]; dp = r[f"p_delta_{dirn}_{tag}"]
                if np.isfinite(dv):
                    dps = ("< 5e-4" if dp <= 1.0 / (a.n_perm + 1) + 1e-12
                           else f"= {dp:.2g}")
                    extra += (f"\n$\\delta$ (x {sym} {'-' if dirn=='loss' else '+'}"
                              f"{a.effect_cut:g}) = {dv:+.2f}, $P$ {dps}"
                              f"  [n = {int(r[f'n_{dirn}_{tag}'])}]")
            ax.text(0.03, 0.97,
                    f"$\\rho$ = {rho:+.2f} [{lo:+.2f}, {hi:+.2f}]\n$P$ {pstr}\nn = {int(nn):,}" + extra,
                    transform=ax.transAxes, ha="left", va="top", fontsize=7.6,
                    linespacing=1.45, zorder=6,
                    bbox=dict(boxstyle="round,pad=0.38", facecolor="white",
                              edgecolor="#CCCCCC", linewidth=0.6, alpha=0.92))
            ax.set_title(f"{nm} — {lab}", fontweight="bold", fontsize=10)
            ax.set_xlabel(f"{nm} change, log$_2$(B10/Parental)")
            if j == 0:
                ax.set_ylabel("bin mean expression change (log$_2$FC)")
            ax.tick_params(labelsize=9)
    # The A-vs-B difference describes the PAIR, not either panel, so it is
    # centred beneath the row rather than placed inside one of the two panels.
    # constrained_layout finalises axes positions only at draw time, hence the
    # explicit draw before reading them.
    fig.canvas.draw()
    for i in range(nm_n):
        r = st.iloc[i]
        dstr = ("< 5e-4" if r.p_diff <= 2.0 / (a.n_boot + 1) else f"= {r.p_diff:.2g}")
        # Use the TIGHT bbox (which includes tick labels and the x-axis label)
        # so the caption clears them instead of landing on top.
        rend = fig.canvas.get_renderer()
        t0 = axes[i][0].get_tightbbox(rend).transformed(fig.transFigure.inverted())
        t1 = axes[i][1].get_tightbbox(rend).transformed(fig.transFigure.inverted())
        b0, b1 = t0, t1
        fig.text((b0.x0 + b1.x1) / 2, b0.y0 - 0.012,
                 f"$\\Delta\\rho$ ({L1} $-$ {L2}) = {r.rho_diff:+.2f} "
                 f"[{r.diff_CI_lo:+.2f}, {r.diff_CI_hi:+.2f}],  $P$ {dstr}",
                 ha="center", va="top", fontsize=8.2,
                 bbox=dict(boxstyle="round,pad=0.34", facecolor="#F6F6F6",
                           edgecolor="#CCCCCC", linewidth=0.6))
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(a.outdir, f"stratified_scatter.{ext}"),
                    bbox_inches="tight", pad_inches=0.04, dpi=600)
    plt.close(fig)

    df.drop(columns=["_start_"], errors="ignore").to_csv(
        os.path.join(a.outdir, "stratified_bin_table.txt"), sep="\t", index=False)
    with open(os.path.join(a.outdir, "stratified_diagnostics.txt"), "w") as f:
        f.write("\n".join(diag) + "\n")

    print(f"\n[INFO] Outputs -> {a.outdir}")
    print("[NOTE] Read p_diff / q_diff, NOT the overlap of the two CIs, to decide "
          "whether the correlation differs between strata. Overlapping "
          "intervals do not imply the difference is absent, and two intervals "
          "that both exclude zero do not imply they differ from each other.")


if __name__ == "__main__":
    main()
