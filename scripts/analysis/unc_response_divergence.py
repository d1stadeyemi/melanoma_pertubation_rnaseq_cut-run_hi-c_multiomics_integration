#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unc_response_divergence.py 
============================================

Did constricted migration REPROGRAM the H3K27me3 response to EZH2 inhibition?

The idea
--------
unc_selective_erasure.py (Question A) already wrote a master table with, per
250 kb bin, the within-cell-line UNC effect for BOTH cell lines:
    dUNC_b10       = log2(B10-UNC / B10-none)
    dUNC_parental  = log2(Parental-UNC / Parental-none)
plus the migration compartment annotations (mig_dir, is_mig_switcher,
is_reverter_b10, is_gainer, is_stable, ...).

If migration did NOT reprogram the drug response, the loci that lose H3K27me3
under UNC are the SAME in both cell lines -> dUNC_b10 and dUNC_parental are
concordant, and any bins where they disagree are scattered at random.

If migration DID reprogram it, the response DIVERGES between the cell lines, and
-- the striking part -- that divergence is concentrated at exactly the loci
migration structurally reprogrammed (switchers / reverters).

Why this is scale-free (same logic as Question A)
-------------------------------------------------
Each cell line's dUNC carries its own unknown global offset:
    dUNC_b10_i      = trueB10_i      + c_b10
    dUNC_parental_i = trueParental_i + c_parental
We compare the two cell lines only through RANKS. Ranks are invariant to adding a
constant, so c_b10 and c_parental drop out: the Spearman concordance and the
per-bin rank-divergence are free of both global offsets. No absolute magnitude is
used anywhere.

What it does
------------
1.  Loads the master table from Question A (no re-counting).
2.  Concordance: Spearman rho between dUNC_b10 and dUNC_parental over all tested
    bins, and separately within reprogrammed vs stable bins (a decoupled response
    shows LOWER rho inside reprogrammed bins).
3.  Per-bin divergence: percentile-rank each cell line's dUNC, then
        signed_div = pctrank(dUNC_b10) - pctrank(dUNC_parental)
    ( <0 = MORE depleted in B10 ; >0 = more depleted in parental ),
    abs_div = |signed_div|, and flags the top --divergent-pct as 'divergent'.
4.  Enrichment (the Question-D test): is abs_div HIGHER at migration-reprogrammed
    loci (switchers / reverters) than at baseline-matched stable loci? Uses the
    same matched Cliff's-delta + within-stratum label permutation + bootstrap CI
    machinery as Question A. A positive, significant delta = migration reprogrammed
    the EZH2i response specifically where it reprogrammed the genome.
    A categorical Fisher odds-ratio (divergent x reprogrammed) is reported too.
5.  Writes an annotated per-bin table, a results table, diagnostic plots, and a
    provenance MANIFEST.

Sign conventions (documented to avoid the historical Cliff's-delta sign bug)
----------------------------------------------------------------------------
  * Cliff's delta here uses value = abs_div, target = reprogrammed, stable = control.
    delta > 0  <=>  reprogrammed bins are MORE divergent (hypothesis-consistent).
    (Note this is the opposite expected sign from Question A, where the value was
     dUNC and depletion made delta negative.)

Environment
-----------
conda env: multiomics-py  (numpy, scipy, pandas, matplotlib). No deeptools / BAMs.

------------------------------------------------------------------------------
USAGE EXAMPLES
------------------------------------------------------------------------------
# 0) Validate the stats core synthetically (no data needed):
python scripts/analysis/unc_response_divergence.py --self-test

# 1) Standard run against the Question-A master table:
python scripts/analysis/unc_response_divergence.py \
    --master  results/exploratory/unc_selective_erasure/master_bins_H3K27me3_250kb.tsv \
    --outdir  results/exploratory/unc_response_divergence

# 2) Stricter divergence cut (top 5%) and a hard-switch reverter target:
python scripts/analysis/unc_response_divergence.py \
    --master  results/exploratory/unc_selective_erasure/master_bins_H3K27me3_250kb.tsv \
    --divergent-pct 5 --reverter-strict \
    --outdir  results/exploratory/unc_response_divergence_strict

# 3) Emit a divergent-bin gene list for pathway analysis (same hg19 gene->bin table
#    used by Question A), one gene per bin:
python scripts/analysis/unc_response_divergence.py \
    --master  .../master_bins_H3K27me3_250kb.tsv \
    --gene-table data/processed/hg19_gene_bins.tsv \
    --outdir  results/exploratory/unc_response_divergence
------------------------------------------------------------------------------
"""

import argparse
import os
import sys
import subprocess
import datetime as _dt

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu, fisher_exact

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


REQUIRED_COLS = [
    "coordinate", "dUNC_b10", "dUNC_parental", "baseline_b10",
    "mig_dir", "is_mig_switcher", "is_reverter_b10", "is_reverter_b10_strict",
    "is_gainer", "is_stable",
]


# =============================================================================
# Stats core  (identical to unc_selective_erasure.py -- copied so this script is
# standalone; offset-cancellation here comes from RANKS, the matched test reuses
# the same abundance-matched permutation design)
# =============================================================================
def cliffs_delta(x, y):
    """Cliff's delta via the Mann-Whitney U identity: delta = 2*U1/(n1*n2) - 1.
    delta > 0  <=>  x tends to be HIGHER than y."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return np.nan
    U1 = mannwhitneyu(x, y, alternative="two-sided").statistic
    return 2.0 * U1 / (n1 * n2) - 1.0


def stratified_matched_test(values, is_target, baseline, n_strata=10,
                            n_perm=1000, seed=0, n_boot=0):
    """Baseline-matched Cliff's delta (target vs stable) with within-stratum
    label-permutation p and optional within-stratum bootstrap CI. See Question A
    for the full rationale."""
    values   = np.asarray(values, float)
    is_target = np.asarray(is_target, bool)
    baseline = np.asarray(baseline, float)

    keep = np.isfinite(values) & np.isfinite(baseline)
    v, t, b = values[keep], is_target[keep], baseline[keep]

    edges = np.quantile(b, np.linspace(0, 1, n_strata + 1))
    edges[0] -= 1e-9; edges[-1] += 1e-9
    strata = np.clip(np.digitize(b, edges[1:-1]), 0, n_strata - 1)

    def pooled_delta(labels):
        num = den = 0.0
        for s in range(n_strata):
            idx = strata == s
            xs = v[idx & labels]; ys = v[idx & ~labels]
            if len(xs) == 0 or len(ys) == 0:
                continue
            d = cliffs_delta(xs, ys)
            if not np.isfinite(d):
                continue
            w = len(xs) * len(ys)
            num += d * w; den += w
        return num / den if den > 0 else np.nan

    obs = pooled_delta(t)

    rng = np.random.default_rng(seed)
    perm = np.empty(n_perm)
    for i in range(n_perm):
        lab = np.zeros_like(t)
        for s in range(n_strata):
            idx = np.where(strata == s)[0]
            k = int(t[idx].sum())
            if k == 0 or k == len(idx):
                lab[idx] = t[idx]; continue
            lab[rng.choice(idx, size=k, replace=False)] = True
        perm[i] = pooled_delta(lab)
    perm = perm[np.isfinite(perm)]
    p_two = (np.sum(np.abs(perm) >= abs(obs)) + 1) / (len(perm) + 1) if len(perm) else np.nan

    ci_low = ci_high = np.nan
    if n_boot > 0:
        tgt_idx = [np.where((strata == s) & t)[0] for s in range(n_strata)]
        stb_idx = [np.where((strata == s) & ~t)[0] for s in range(n_strata)]
        boots = np.empty(n_boot)
        for j in range(n_boot):
            num = den = 0.0
            for s in range(n_strata):
                ti, si = tgt_idx[s], stb_idx[s]
                if len(ti) == 0 or len(si) == 0:
                    continue
                xs = v[rng.choice(ti, size=len(ti), replace=True)]
                ys = v[rng.choice(si, size=len(si), replace=True)]
                d = cliffs_delta(xs, ys)
                if not np.isfinite(d):
                    continue
                w = len(ti) * len(si)
                num += d * w; den += w
            boots[j] = num / den if den > 0 else np.nan
        boots = boots[np.isfinite(boots)]
        if len(boots):
            ci_low, ci_high = np.percentile(boots, [2.5, 97.5])

    return {"cliffs_delta": obs, "perm_p_two_sided": p_two,
            "ci_low": ci_low, "ci_high": ci_high,
            "n_target": int(t.sum()), "n_stable": int((~t).sum()),
            "median_target": float(np.nanmedian(v[t])) if t.any() else np.nan,
            "median_stable": float(np.nanmedian(v[~t])) if (~t).any() else np.nan}


# =============================================================================
# Load + divergence
# =============================================================================
def load_master(path, min_baseline_quantile):
    m = pd.read_csv(path, sep="\t")
    miss = [c for c in REQUIRED_COLS if c not in m.columns]
    if miss:
        sys.exit(f"[ERROR] master table missing columns: {miss}. "
                 f"Re-run unc_selective_erasure.py to regenerate it.")
    # keep bins with a finite response in BOTH cell lines (needed for ranks)
    m = m[np.isfinite(m["dUNC_b10"]) & np.isfinite(m["dUNC_parental"])].copy()
    # optional extra occupancy filter (usually already applied upstream)
    if min_baseline_quantile > 0 and "baseline_parental" in m.columns:
        fb = np.nanquantile(m["baseline_b10"], min_baseline_quantile)
        fp = np.nanquantile(m["baseline_parental"], min_baseline_quantile)
        m = m[(m["baseline_b10"] >= fb) | (m["baseline_parental"] >= fp)].copy()
    return m.reset_index(drop=True)


def compute_divergence(m, divergent_pct):
    """Percentile-rank each cell line's dUNC and derive the divergence score."""
    m = m.copy()
    m["pr_b10"] = m["dUNC_b10"].rank(pct=True)          # low rank = most depleted
    m["pr_parental"] = m["dUNC_parental"].rank(pct=True)
    m["signed_div"] = m["pr_b10"] - m["pr_parental"]    # <0: B10 more depleted
    m["abs_div"] = m["signed_div"].abs()
    thr = m["abs_div"].quantile(1 - divergent_pct / 100.0)
    m["is_divergent"] = m["abs_div"] >= thr
    m["div_direction"] = np.where(
        m["signed_div"] < 0, "b10_more_depleted",
        np.where(m["signed_div"] > 0, "parental_more_depleted", "tie"))
    return m, thr


# =============================================================================
# Concordance
# =============================================================================
def spearman_on(m, mask=None):
    d = m if mask is None else m[mask]
    if len(d) < 10:
        return {"rho": np.nan, "p": np.nan, "n": int(len(d))}
    rho, p = spearmanr(d["dUNC_b10"], d["dUNC_parental"])
    return {"rho": float(rho), "p": float(p), "n": int(len(d))}


# =============================================================================
# Plots
# =============================================================================
def plot_rank_scatter(m, rho_all, outpath):
    try:
        fig, ax = plt.subplots(figsize=(8, 8), constrained_layout=True)
        stab = m["is_stable"].values
        sw = m["is_mig_switcher"].values
        ax.scatter(m["pr_parental"][stab], m["pr_b10"][stab], s=8, alpha=0.25,
                   color="0.6", label=f"stable (n={int(stab.sum())})")
        ax.scatter(m["pr_parental"][sw], m["pr_b10"][sw], s=10, alpha=0.5,
                   color="#080cee", label=f"switchers (n={int(sw.sum())})")
        dv = m["is_divergent"].values & sw
        ax.scatter(m["pr_parental"][dv], m["pr_b10"][dv], s=14, alpha=0.8,
                   color="#d62728", label=f"divergent switchers (n={int(dv.sum())})")
        ax.plot([0, 1], [0, 1], color="k", ls=":", lw=1)
        ax.set_xlabel("pctrank  dUNC  Parental   (left = more depleted)", fontsize=12)
        ax.set_ylabel("pctrank  dUNC  B10        (down = more depleted)", fontsize=12)
        ax.tick_params(axis="both", which="major", labelsize=12)
        ax.set_title(f"UNC response: B10 vs Parental\nSpearman rho = {rho_all:.3f}", fontsize=12)
        ax.legend(fontsize=12, frameon=False, loc="upper left", 
                  prop={"weight": "bold"}, markerscale=2.0)
        fig.savefig(outpath, dpi=150); plt.close(fig)
    except Exception as e:
        print(f"[warn] rank scatter failed: {e}")


def plot_divergence_hist(m, outpath):
    try:
        fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
        ax.hist(m.loc[m["is_stable"], "signed_div"].dropna(), bins=40, density=True,
                alpha=0.5, label="stable", color="0.6")
        ax.hist(m.loc[m["is_mig_switcher"], "signed_div"].dropna(), bins=40,
                density=True, alpha=0.5, label="switchers", color="#1f77b4")
        ax.axvline(0, color="k", ls=":", lw=1)
        ax.set_xlabel("signed divergence  =  pctrank(B10) - pctrank(Parental)\n"
                      "(<0: B10 more depleted ; >0: parental more depleted)", fontsize=12)
        ax.set_ylabel("density", fontsize=12)
        ax.set_title("Response-divergence by set", fontsize=12)
        ax.legend(fontsize=12, frameon=False)
        fig.savefig(outpath, dpi=150); plt.close(fig)
    except Exception as e:
        print(f"[warn] divergence hist failed: {e}")


# =============================================================================
# Optional gene mapping for the divergent set
# =============================================================================
def map_divergent_genes(m, gene_table_path, outdir, genes_per_bin="one"):
    gt = pd.read_csv(gene_table_path, sep="\t")
    if not {"gene_name", "coordinate"}.issubset(gt.columns):
        sys.exit("[ERROR] gene table needs columns gene_name, coordinate")
    gt = gt[gt["coordinate"].isin(set(m["coordinate"]))].copy()
    if genes_per_bin == "one" and {"TSS", "bin_start", "bin_end"}.issubset(gt.columns):
        ctr = (gt["bin_start"].astype(float) + gt["bin_end"].astype(float)) / 2
        gt["_d"] = (gt["TSS"].astype(float) - ctr).abs()
        gt = gt.sort_values("_d").drop_duplicates("coordinate").drop(columns="_d")
    elif genes_per_bin == "one":
        gt = gt.drop_duplicates("coordinate")
    gldir = os.path.join(outdir, "genelists"); os.makedirs(gldir, exist_ok=True)
    with open(os.path.join(gldir, "background_universe.txt"), "w") as fh:
        fh.write("\n".join(sorted(set(gt["gene_name"]))) + "\n")
    sets = {
        "divergent_all":        m["is_divergent"].values,
        "divergent_b10_more":   (m["is_divergent"] & (m["signed_div"] < 0)).values,
        "divergent_parental_more": (m["is_divergent"] & (m["signed_div"] > 0)).values,
        "divergent_switchers":  (m["is_divergent"] & m["is_mig_switcher"]).values,
        "divergent_reverters":  (m["is_divergent"] & m["is_reverter_b10"]).values,
    }
    for name, mask in sets.items():
        coords = set(m.loc[mask, "coordinate"])
        gl = sorted(set(gt.loc[gt["coordinate"].isin(coords), "gene_name"]))
        with open(os.path.join(gldir, f"{name}.txt"), "w") as fh:
            fh.write("\n".join(gl) + "\n")
        print(f"    {name:26s} {int(mask.sum()):>5} bins -> {len(gl):>5} genes")
    gt[["gene_name", "coordinate"]].to_csv(
        os.path.join(gldir, "gene_assignment.tsv"), sep="\t", index=False)
    print(f"[info] divergent-set gene lists -> {gldir}")


# =============================================================================
# Provenance
# =============================================================================
def git_hash():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"],
                                       stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "NA"


def write_manifest(outdir, args, n_bins, concord, results_path):
    lines = ["# MANIFEST  --  unc_response_divergence.py (Question D)",
             f"generated  : {_dt.datetime.now().isoformat(timespec='seconds')}",
             f"git_commit : {git_hash()}",
             f"python     : {sys.version.split()[0]}",
             "", "## parameters"]
    for k, v in sorted(vars(args).items()):
        lines.append(f"  {k}: {v}")
    lines += ["", "## concordance (Spearman rho of dUNC_b10 vs dUNC_parental)"]
    for k, r in concord.items():
        lines.append(f"  {k:14s} rho={r['rho']:+.3f}  p={r['p']:.3g}  n={r['n']}")
    lines += ["", f"## bins analysed : {n_bins}", f"## results : {results_path}", ""]
    with open(os.path.join(outdir, "MANIFEST.md"), "w") as fh:
        fh.write("\n".join(lines))


# =============================================================================
# Main
# =============================================================================
def main(args):
    os.makedirs(args.outdir, exist_ok=True)
    m = load_master(args.master, args.min_baseline_quantile)
    m, thr = compute_divergence(m, args.divergent_pct)
    print(f"[info] {len(m)} bins; divergence threshold (top {args.divergent_pct}%) "
          f"abs_div >= {thr:.3f}; {int(m['is_divergent'].sum())} divergent bins")

    # ---- concordance ----
    concord = {
        "overall":    spearman_on(m),
        "switchers":  spearman_on(m, m["is_mig_switcher"].values),
        "reverters":  spearman_on(m, m["is_reverter_b10"].values),
        "stable":     spearman_on(m, m["is_stable"].values),
    }
    print("[info] Spearman concordance dUNC_b10 vs dUNC_parental:")
    for k, r in concord.items():
        print(f"    {k:10s} rho={r['rho']:+.3f}  p={r['p']:.3g}  n={r['n']}")
    print("    (lower rho inside reprogrammed sets than stable = response decoupled "
          "by migration)")

    # ---- enrichment: is divergence higher at reprogrammed than matched stable? ----
    reverter_col = "is_reverter_b10_strict" if args.reverter_strict else "is_reverter_b10"
    targets = {
        "migration_switchers": m["is_mig_switcher"].values,
        "switchers_A_to_B":    (m["mig_dir"] == "A_to_B").values,
        "switchers_B_to_A":    (m["mig_dir"] == "B_to_A").values,
        "UNC_reverters":       m[reverter_col].values,
    }
    stab = m["is_stable"].values
    base = m["baseline_b10"].values

    rows = []
    for i, (name, tmask) in enumerate(targets.items()):
        sel = tmask | stab
        res = stratified_matched_test(
            values=m["abs_div"].values[sel], is_target=tmask[sel],
            baseline=base[sel], n_strata=args.n_strata, n_perm=args.n_perm,
            seed=args.seed + i, n_boot=args.n_boot)
        # categorical companion: divergent x reprogrammed odds ratio (vs stable)
        a = int((tmask & m["is_divergent"].values).sum())          # reprog & divergent
        b_ = int((tmask & ~m["is_divergent"].values).sum())        # reprog & concordant
        c = int((stab & m["is_divergent"].values).sum())           # stable & divergent
        d = int((stab & ~m["is_divergent"].values).sum())          # stable & concordant
        OR, fisher_p = (np.nan, np.nan)
        if min(a + b_, c + d) > 0:
            OR, fisher_p = fisher_exact([[a, b_], [c, d]], alternative="greater")
        interp = ("MORE divergent at reprogrammed (migration reprogrammed the response)"
                  if (np.isfinite(res["cliffs_delta"]) and res["cliffs_delta"] > 0
                      and res["perm_p_two_sided"] < 0.05) else "no excess divergence")
        rows.append({"target_set": name, "n_target": res["n_target"],
                     "n_stable": res["n_stable"], "cliffs_delta": res["cliffs_delta"],
                     "ci_low": res["ci_low"], "ci_high": res["ci_high"],
                     "perm_p_two_sided": res["perm_p_two_sided"],
                     "fisher_OR_divergent": OR, "fisher_p": fisher_p,
                     "interpretation": interp})
        print(f"  {name:20s} delta={res['cliffs_delta']:+.3f} "
              f"p={res['perm_p_two_sided']:.4f}  OR(divergent)={OR:.2f} "
              f"-> {interp}")

    results = pd.DataFrame(rows)
    results_path = os.path.join(args.outdir, "divergence_enrichment.tsv")
    results.to_csv(results_path, sep="\t", index=False)

    # divergence-direction breakdown among reprogrammed-divergent bins
    dv_sw = m[m["is_divergent"] & m["is_mig_switcher"]]
    print(f"[info] of {len(dv_sw)} divergent switchers: "
          f"{int((dv_sw['signed_div']<0).sum())} B10-more-depleted, "
          f"{int((dv_sw['signed_div']>0).sum())} parental-more-depleted")

    # ---- outputs ----
    keep_cols = (["coordinate", "chr", "dUNC_b10", "dUNC_parental", "pr_b10",
                  "pr_parental", "signed_div", "abs_div", "is_divergent",
                  "div_direction", "mig_dir", "is_mig_switcher", "is_reverter_b10",
                  "is_gainer", "is_stable"])
    m[[c for c in keep_cols if c in m.columns]].to_csv(
        os.path.join(args.outdir, "master_divergence.tsv"), sep="\t", index=False)
    plot_rank_scatter(m, concord["overall"]["rho"],
                      os.path.join(args.outdir, "rank_scatter_b10_vs_parental.png"))
    plot_divergence_hist(m, os.path.join(args.outdir, "signed_divergence_hist.png"))

    if args.gene_table:
        map_divergent_genes(m, args.gene_table, args.outdir, args.genes_per_bin)

    write_manifest(args.outdir, args, len(m), concord, results_path)
    print(f"[done] results -> {results_path}")


# =============================================================================
# Self-test
# =============================================================================
def self_test():
    """(a) when the B10 response is scrambled at 'reprogrammed' bins, those bins
    are more divergent (delta>0, sig) and overall concordance drops;
    (b) adding large global offsets to either dUNC does NOT change the answer
    (rank-based -> offset-free); (c) with no scrambling, reprogrammed bins are not
    more divergent."""
    print("=== SELF-TEST: Question D divergence core ===")
    rng = np.random.default_rng(0)
    N = 5000
    base = rng.normal(0, 1, N)
    par = rng.normal(0, 1, N)
    b10 = par + rng.normal(0, 0.3, N)              # concordant baseline
    is_t = np.zeros(N, bool); is_t[rng.choice(N, 600, replace=False)] = True

    def run(b10v, parv, c_b10=0.0, c_par=0.0):
        d = pd.DataFrame({"coordinate": np.arange(N).astype(str),
                          "dUNC_b10": b10v + c_b10, "dUNC_parental": parv + c_par,
                          "baseline_b10": base, "baseline_parental": base,
                          "mig_dir": np.where(is_t, "A_to_B", "Stable_A"),
                          "is_mig_switcher": is_t,
                          "is_reverter_b10": is_t, "is_reverter_b10_strict": is_t,
                          "is_gainer": False, "is_stable": ~is_t})
        d, _ = compute_divergence(d, 10)
        sel = d["is_mig_switcher"].values | d["is_stable"].values
        r = stratified_matched_test(d["abs_div"].values[sel],
                                    d["is_mig_switcher"].values[sel],
                                    d["baseline_b10"].values[sel],
                                    n_strata=10, n_perm=300, n_boot=0, seed=1)
        rho = spearmanr(d["dUNC_b10"], d["dUNC_parental"]).statistic
        return r, rho

    # (a) scramble B10 response at target bins -> divergence injected there
    b10_scr = b10.copy(); b10_scr[is_t] = rng.normal(0, 1, is_t.sum())
    r, rho = run(b10_scr, par)
    print(f"(a) scrambled targets   delta={r['cliffs_delta']:+.3f} "
          f"p={r['perm_p_two_sided']:.4f}  rho_all={rho:+.3f}")
    assert r["cliffs_delta"] > 0.1 and r["perm_p_two_sided"] < 0.05, "should detect divergence"

    # (b) same truth, huge offsets on both -> unchanged (rank/Spearman are offset-free)
    r2, rho2 = run(b10_scr, par, c_b10=-7.0, c_par=+9.0)
    print(f"(b) + big offsets       delta={r2['cliffs_delta']:+.3f} "
          f"p={r2['perm_p_two_sided']:.4f}  rho_all={rho2:+.3f}  (offset-free)")
    assert abs(r2["cliffs_delta"] - r["cliffs_delta"]) < 1e-9, "offsets changed ranks!"

    # (c) no scrambling -> reprogrammed not more divergent
    r3, rho3 = run(b10, par)
    print(f"(c) no scramble         delta={r3['cliffs_delta']:+.3f} "
          f"p={r3['perm_p_two_sided']:.4f}  rho_all={rho3:+.3f}  (should be ~null)")
    assert abs(r3["cliffs_delta"]) < 0.12 or r3["perm_p_two_sided"] >= 0.05, "false positive"
    print("ALL SELF-TESTS PASSED")


# =============================================================================
# CLI
# =============================================================================
def build_parser():
    p = argparse.ArgumentParser(
        description="Question D: did migration reprogram the H3K27me3 response to "
                    "EZH2i? (rank concordance + divergence enrichment, scale-free).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--self-test", action="store_true",
                   help="run synthetic validation of the stats core and exit")
    p.add_argument("--master",
                   help="master_bins_H3K27me3_250kb.tsv from unc_selective_erasure.py")
    p.add_argument("--outdir", default="results/exploratory/unc_response_divergence")
    p.add_argument("--divergent-pct", type=float, default=10.0,
                   help="top percentile of |rank divergence| flagged as 'divergent'")
    p.add_argument("--reverter-strict", action="store_true",
                   help="use the hard-switch reverter set as a target")
    p.add_argument("--min-baseline-quantile", type=float, default=0.0,
                   help="extra occupancy filter (usually already applied upstream; "
                        "0 = use the master as-is)")
    p.add_argument("--gene-table", default=None,
                   help="hg19 gene->bin TSV; if given, emit divergent-set gene lists "
                        "+ background for the pathway-enrichment script")
    p.add_argument("--genes-per-bin", choices=["one", "all"], default="one",
                   help="one representative gene per bin (recommended) or all genes")
    p.add_argument("--n-strata", type=int, default=10)
    p.add_argument("--n-perm", type=int, default=1000)
    p.add_argument("--n-boot", type=int, default=1000,
                   help="bootstrap resamples for the Cliff's-delta CI; 0 disables")
    p.add_argument("--seed", type=int, default=0)
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    if args.self_test:
        self_test(); sys.exit(0)
    if not args.master:
        sys.exit("[ERROR] --master is required (or use --self-test)")
    main(args)