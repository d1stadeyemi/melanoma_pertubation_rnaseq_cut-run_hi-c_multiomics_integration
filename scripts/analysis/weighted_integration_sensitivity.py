#!/usr/bin/env python3
"""
weighted_integration_sensitivity.py

PURPOSE
  Sensitivity analysis for the Hi-C x CUT&RUN x RNA-seq bin-level integration.

  Each bin's expression value (`expr`) is the MEAN log2FC of the genes falling
  in that bin. Roughly half to two-thirds of 250-kb bins contain exactly ONE
  gene, so their `expr` is a single noisy gene-level estimate, while other bins
  average several genes and are far better measured. An unweighted Spearman
  treats all of these equally.

  This script recomputes the mark-vs-expression association three ways so you
  can see whether it is carried by, or diluted by, the poorly-measured bins:

    1. UNWEIGHTED           - all bins equal (what the pipeline reports)
    2. WEIGHTED by n_genes  - bins averaging more genes count more
    3. WEIGHTED by 1/sd^2   - inverse-variance of the bin mean (see caveat)

  It also reports the hard-threshold series (n_genes >= 1,2,3) for comparison,
  since thresholds discard data whereas weights retain all of it.

INTERPRETATION
  - Association STRENGTHENS under weighting -> real effect previously
    attenuated by noise in single-gene bins.
  - Association WEAKENS/VANISHES under weighting -> it was being carried by the
    noisiest bins; treat as not robust.
  - Association STABLE -> not sensitive to bin measurement quality.

CAVEAT ON INVERSE-VARIANCE WEIGHTS
  `expr_sd` is the SD of gene LFCs within a bin and is undefined (NaN) for
  1-gene bins and 0 for bins whose genes happen to agree exactly. The SE of the
  bin mean is sd/sqrt(n_genes). Single-gene bins are handled by --single_sd
  (default: median of the observed SDs), which is an ASSUMPTION, not a
  measurement -- it says "assume a lone gene is about as variable as a typical
  bin's genes". Results under this weighting should be read as a sensitivity
  check, not a definitive estimate. The n_genes weighting makes no such
  assumption and is the more conservative of the two.

WEIGHTED SPEARMAN
  Computed as the weighted Pearson correlation of the RANKS of x and y (the
  standard extension of Spearman to weighted data). CIs come from a weighted
  bootstrap; p-values from a permutation test that shuffles y against x while
  holding the weights with their bins.

USAGE
  python weighted_integration_sensitivity.py \
      --base results/exploratory/integration_hic_cutrun_rnaseq

  # or point at specific run directories
  python weighted_integration_sensitivity.py \
      --dirs results/.../b10U_vs_parU results/.../b10D_vs_parD
"""

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats


# -----------------------------------------------------------------------------
# Weighted statistics
# -----------------------------------------------------------------------------
def weighted_pearson(x, y, w):
    """Weighted Pearson correlation."""
    w = np.asarray(w, dtype=float)
    wsum = w.sum()
    if wsum <= 0:
        return np.nan
    mx = np.sum(w * x) / wsum
    my = np.sum(w * y) / wsum
    cov = np.sum(w * (x - mx) * (y - my)) / wsum
    vx = np.sum(w * (x - mx) ** 2) / wsum
    vy = np.sum(w * (y - my) ** 2) / wsum
    if vx <= 0 or vy <= 0:
        return np.nan
    return cov / np.sqrt(vx * vy)


def weighted_spearman(x, y, w):
    """Weighted Spearman = weighted Pearson on ranks."""
    return weighted_pearson(stats.rankdata(x), stats.rankdata(y), w)


def weighted_spearman_ci(x, y, w, n_boot=5000, seed=0):
    """Percentile bootstrap CI for the weighted Spearman."""
    rng = np.random.default_rng(seed)
    n = len(x)
    vals = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        vals[b] = weighted_spearman(x[idx], y[idx], w[idx])
    vals = vals[np.isfinite(vals)]
    if vals.size < 100:
        return (np.nan, np.nan)
    return tuple(np.percentile(vals, [2.5, 97.5]))


def weighted_spearman_perm_p(x, y, w, n_perm=10000, seed=0):
    """Two-sided permutation p: shuffle y, keep (x, w) paired with their bins."""
    rng = np.random.default_rng(seed)
    obs = weighted_spearman(x, y, w)
    if not np.isfinite(obs):
        return np.nan
    count = 0
    for _ in range(n_perm):
        yp = rng.permutation(y)
        v = weighted_spearman(x, yp, w)
        if np.isfinite(v) and abs(v) >= abs(obs):
            count += 1
    return (count + 1) / (n_perm + 1)


# -----------------------------------------------------------------------------
# Per-run analysis
# -----------------------------------------------------------------------------
def analyse_run(path, n_boot, n_perm, single_sd_mode, seed):
    name = os.path.basename(os.path.normpath(path))
    table = os.path.join(path, "integration_bin_table.txt")
    if not os.path.isfile(table):
        print(f"\n{'='*78}\n{name}\n{'='*78}")
        print(f"  (no integration_bin_table.txt in {path}) -- skipped")
        return None

    df = pd.read_csv(table, sep=None, engine="python")

    mark_cols = [c for c in df.columns if "log2ratio" in c.lower()]
    if not mark_cols:
        mark_cols = [c for c in df.columns if "h3k" in c.lower()]
    if not mark_cols or "expr" not in df.columns:
        print(f"\n{'='*78}\n{name}\n{'='*78}")
        print("  (could not find mark / expr columns) -- skipped")
        return None
    mark = mark_cols[0]

    need = [mark, "expr", "n_genes"]
    if "expr_sd" in df.columns:
        need.append("expr_sd")
    d = df[need].apply(pd.to_numeric, errors="coerce")
    d = d.dropna(subset=[mark, "expr", "n_genes"])

    x = d[mark].to_numpy(float)
    y = d["expr"].to_numpy(float)
    ng = d["n_genes"].to_numpy(float)
    n = len(d)

    print(f"\n{'='*78}\n{name}\n{'='*78}")
    print(f"  mark column : {mark}")
    print(f"  bins        : {n}")
    print(f"  n_genes     : median={np.median(ng):.0f}  "
          f"range={ng.min():.0f}-{ng.max():.0f}  "
          f"1-gene bins={int((ng==1).sum())} ({100*(ng==1).mean():.0f}%)")

    rows = []

    # --- 1. unweighted -------------------------------------------------------
    rho_u, p_u = stats.spearmanr(x, y)
    ci_u = weighted_spearman_ci(x, y, np.ones(n), n_boot, seed)
    p_perm_u = weighted_spearman_perm_p(x, y, np.ones(n), n_perm, seed)
    rows.append(("unweighted", n, rho_u, ci_u, p_perm_u, p_u))

    # --- 2. weighted by n_genes ---------------------------------------------
    w_ng = ng.copy()
    rho_g = weighted_spearman(x, y, w_ng)
    ci_g = weighted_spearman_ci(x, y, w_ng, n_boot, seed)
    p_g = weighted_spearman_perm_p(x, y, w_ng, n_perm, seed)
    rows.append(("weighted: n_genes", n, rho_g, ci_g, p_g, np.nan))

    # --- 3. weighted by inverse variance of the bin mean ---------------------
    if "expr_sd" in d.columns:
        sd = d["expr_sd"].to_numpy(float)
        obs_sd = sd[np.isfinite(sd) & (sd > 0)]
        if obs_sd.size:
            if single_sd_mode == "median":
                fill = float(np.median(obs_sd))
            elif single_sd_mode == "mean":
                fill = float(np.mean(obs_sd))
            else:
                fill = float(single_sd_mode)
            sd_filled = np.where(np.isfinite(sd) & (sd > 0), sd, fill)
            se = sd_filled / np.sqrt(np.maximum(ng, 1.0))
            w_iv = 1.0 / np.maximum(se, 1e-12) ** 2
            # cap extreme weights so one bin cannot dominate
            cap = np.percentile(w_iv, 99)
            w_iv = np.minimum(w_iv, cap)
            rho_v = weighted_spearman(x, y, w_iv)
            ci_v = weighted_spearman_ci(x, y, w_iv, n_boot, seed)
            p_v = weighted_spearman_perm_p(x, y, w_iv, n_perm, seed)
            rows.append((f"weighted: 1/SE^2 (single_sd={fill:.3g})",
                         n, rho_v, ci_v, p_v, np.nan))
        else:
            print("  [note] no usable expr_sd values; inverse-variance skipped")
    else:
        print("  [note] no expr_sd column; inverse-variance skipped")

    # --- print weighted comparison ------------------------------------------
    print(f"\n  {'method':<34} {'n':>5} {'rho':>8} {'95% CI':>20} "
          f"{'p_perm':>9} {'p_naive':>9}")
    print("  " + "-" * 88)
    for label, nn, rho, ci, pp, pn in rows:
        ci_s = (f"[{ci[0]:+.3f},{ci[1]:+.3f}]"
                if np.isfinite(ci[0]) else "        n/a       ")
        pn_s = f"{pn:.2e}" if np.isfinite(pn) else "     -"
        print(f"  {label:<34} {nn:>5} {rho:>+8.3f} {ci_s:>20} "
              f"{pp:>9.4f} {pn_s:>9}")

    # --- threshold series for comparison ------------------------------------
    print(f"\n  hard thresholds (discard data; shown for comparison):")
    for thr in (1, 2, 3):
        m = ng >= thr
        if m.sum() > 10:
            r, p = stats.spearmanr(x[m], y[m])
            print(f"    n_genes >= {thr}: n={int(m.sum()):4d}  "
                  f"rho={r:+.3f}  p_naive={p:.2e}")

    # --- verdict -------------------------------------------------------------
    print()
    if np.isfinite(rho_u) and np.isfinite(rho_g):
        change = abs(rho_g) - abs(rho_u)
        same_sign = np.sign(rho_g) == np.sign(rho_u)
        if not same_sign:
            verdict = "SIGN FLIPS under weighting -- not robust"
        elif change > 0.05:
            verdict = ("STRENGTHENS under weighting -- consistent with a real "
                       "effect diluted by noisy single-gene bins")
        elif change < -0.05:
            verdict = ("WEAKENS under weighting -- association was carried by "
                       "poorly-measured bins; treat as not robust")
        else:
            verdict = "STABLE under weighting -- not sensitive to bin quality"
        print(f"  VERDICT: {verdict}")
        print(f"           |rho| unweighted={abs(rho_u):.3f} -> "
              f"n_genes-weighted={abs(rho_g):.3f} (change {change:+.3f})")

    return dict(run=name, mark=mark, n=n, rho_unw=rho_u, rho_ng=rho_g)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--base", help="Parent dir; every subdirectory containing "
                                  "integration_bin_table.txt is analysed.")
    g.add_argument("--dirs", nargs="+", help="Explicit run directories.")
    ap.add_argument("--n_boot", type=int, default=5000,
                    help="Bootstrap replicates for CIs (default: %(default)s)")
    ap.add_argument("--n_perm", type=int, default=10000,
                    help="Permutations for p-values (default: %(default)s)")
    ap.add_argument("--single_sd", default="median",
                    help="SD assumed for 1-gene bins in the inverse-variance "
                         "weighting: 'median', 'mean', or a number "
                         "(default: %(default)s)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    if args.base:
        dirs = sorted(os.path.dirname(p) for p in
                      glob.glob(os.path.join(args.base, "*",
                                             "integration_bin_table.txt")))
        if not dirs:
            sys.exit(f"No integration_bin_table.txt found under {args.base}")
    else:
        dirs = args.dirs

    print(f"Analysing {len(dirs)} run(s)")
    print(f"bootstrap={args.n_boot}  permutations={args.n_perm}  "
          f"single_sd={args.single_sd}")

    summary = []
    for d in dirs:
        r = analyse_run(d, args.n_boot, args.n_perm, args.single_sd, args.seed)
        if r:
            summary.append(r)

    if summary:
        print(f"\n\n{'='*78}\nSUMMARY\n{'='*78}")
        print(f"  {'run':<24} {'n':>5} {'rho unw':>9} {'rho n_genes':>12} "
              f"{'change':>9}")
        print("  " + "-" * 64)
        for s in summary:
            ch = abs(s['rho_ng']) - abs(s['rho_unw'])
            print(f"  {s['run']:<24} {s['n']:>5} {s['rho_unw']:>+9.3f} "
                  f"{s['rho_ng']:>+12.3f} {ch:>+9.3f}")
        print("\n  'change' = |rho weighted| - |rho unweighted|. Positive means "
              "the\n  association sharpens when better-measured bins are given "
              "more weight.")


if __name__ == "__main__":
    main()