#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
test_ctcf_anchor_stability.py

Tests whether the constitutive-CTCF reference regions used as Ryder's anchor are
remodeled across the constricted-migration selection. PRIMARY test = migration
axis on UNTREATED tracks (Parental_none vs B10_none) -- clean, no UNC global-loss
confound. Also reports the within-condition UNC axis as DESCRIPTIVE ONLY (RPM
tracks confound uniform global loss; absolute CTCF invariance under UNC needs
IF/western, not this script).

Outputs: <out>_perRegion.tsv, _summary.txt, _ecdf.pdf, _pairedScatter.pdf,
         _deltaHist.pdf, _metaprofile.pdf, _MANIFEST.txt

PRE-SPECIFIED DECISION RULE (migration axis; fix BEFORE looking):
  ANCHOR STABLE      if |median Δlog2| < --stable-med  AND frac(|Δlog2|>1) < --stable-frac
  ANCHOR COMPROMISED otherwise
  defaults: --stable-med 0.20  --stable-frac 0.10
"""
import argparse, os, sys, subprocess, datetime
import numpy as np
import pandas as pd

CANON = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ref", required=True, help="CTCF reference BED (chr-prefixed)")
    p.add_argument("--parental-none", required=True, help="Parental untreated avg bigWig")
    p.add_argument("--b10-none", required=True, help="B10 untreated avg bigWig")
    p.add_argument("--parental-unc", default=None, help="Parental UNC avg bigWig (optional)")
    p.add_argument("--b10-unc", default=None, help="B10 UNC avg bigWig (optional)")
    p.add_argument("--out", required=True, help="output prefix")
    p.add_argument("--pseudocount", type=float, default=1e-3, help="added before log2 (default 1e-3)")
    p.add_argument("--ext", type=int, default=5000, help="flank bp each side of center for metaprofile")
    p.add_argument("--meta-bins", type=int, default=100, help="bins across 2*ext window")
    p.add_argument("--exact", action="store_true", help="exact per-region mean (recommended for small regions)")
    p.add_argument("--canonical-only", action="store_true", default=True, help="restrict to chr1..22,X,Y")
    p.add_argument("--stable-med", type=float, default=0.20, help="|median Δlog2| threshold for STABLE")
    p.add_argument("--stable-frac", type=float, default=0.10, help="frac(|Δlog2|>1) threshold for STABLE")
    return p.parse_args()


def git_hash():
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"],
                                       stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "NA"


def load_regions(bed, canonical_only):
    rows = []
    with open(bed) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith(("#", "track", "browser")):
                continue
            f = ln.rstrip("\n").split("\t")
            c, s, e = f[0], int(f[1]), int(f[2])
            if canonical_only and c not in CANON:
                continue
            if e <= s:
                continue
            rows.append((c, s, e))
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


def chrom_resolver(bw):
    have = set(bw.chroms().keys())
    def resolve(c):
        if c in have:
            return c
        alt = c[3:] if c.startswith("chr") else "chr" + c
        return alt if alt in have else None
    return resolve


def region_means(bwpath, regions, exact):
    import pyBigWig
    bw = pyBigWig.open(bwpath); res = chrom_resolver(bw)
    out = np.full(len(regions), np.nan)
    for i, (c, s, e) in enumerate(zip(regions.chrom.values, regions.start.values, regions.end.values)):
        cc = res(c)
        if cc is None:
            continue
        L = bw.chroms()[cc]; s2, e2 = max(0, s), min(e, L)
        if e2 <= s2:
            continue
        try:
            v = bw.stats(cc, s2, e2, type="mean", exact=exact)[0]
        except (RuntimeError, OverflowError):
            v = None
        if v is not None:
            out[i] = v
    bw.close(); return out


def metaprofile(bwpath, regions, ext, nbins):
    import pyBigWig
    bw = pyBigWig.open(bwpath); res = chrom_resolver(bw)
    acc = np.zeros(nbins); cnt = np.zeros(nbins)
    for c, s, e in zip(regions.chrom.values, regions.start.values, regions.end.values):
        cc = res(c)
        if cc is None:
            continue
        L = bw.chroms()[cc]; ctr = (s + e) // 2; a, b = ctr - ext, ctr + ext
        if a < 0 or b > L:
            continue
        try:
            vals = bw.stats(cc, a, b, type="mean", nBins=nbins, exact=False)
        except (RuntimeError, OverflowError):
            continue
        v = np.array([np.nan if x is None else x for x in vals], dtype=float)
        ok = ~np.isnan(v); acc[ok] += v[ok]; cnt[ok] += 1
    bw.close(); cnt[cnt == 0] = np.nan
    return acc / cnt


def cliffs_delta(a, b):
    a = np.asarray(a); b = np.asarray(b)
    if len(a) * len(b) > 4_000_000:
        rng = np.random.default_rng(0)
        a = rng.choice(a, 2000); b = rng.choice(b, 2000)
    gt = (a[:, None] > b[None, :]).sum()
    lt = (a[:, None] < b[None, :]).sum()
    return (gt - lt) / (len(a) * len(b))


def main():
    args = parse_args()
    from scipy import stats
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    reg = load_regions(args.ref, args.canonical_only)
    if len(reg) == 0:
        sys.exit("ERROR: no reference regions loaded.")

    tracks = {"Parental_none": args.parental_none, "B10_none": args.b10_none}
    if args.parental_unc:
        tracks["Parental_UNC"] = args.parental_unc
    if args.b10_unc:
        tracks["B10_UNC"] = args.b10_unc

    sig = pd.DataFrame({"chrom": reg.chrom, "start": reg.start, "end": reg.end})
    for name, path in tracks.items():
        sig[name] = region_means(path, reg, args.exact)

    keep = sig[["Parental_none", "B10_none"]].notna().all(axis=1)
    n_drop = int((~keep).sum()); sig = sig[keep].reset_index(drop=True)

    pc = args.pseudocount
    log2 = {n: np.log2(sig[n].values + pc) for n in tracks}

    d_mig = log2["B10_none"] - log2["Parental_none"]
    med = float(np.median(d_mig)); q1, q3 = np.percentile(d_mig, [25, 75])
    frac2 = float(np.mean(np.abs(d_mig) > 1))
    try:
        w_p = stats.wilcoxon(log2["B10_none"], log2["Parental_none"]).pvalue
    except ValueError:
        w_p = float("nan")
    cd = cliffs_delta(log2["B10_none"], log2["Parental_none"])
    stable = (abs(med) < args.stable_med) and (frac2 < args.stable_frac)
    verdict = "ANCHOR STABLE" if stable else "ANCHOR COMPROMISED"

    desc = []
    for cond, none_k, unc_k in [("Parental", "Parental_none", "Parental_UNC"),
                                ("B10", "B10_none", "B10_UNC")]:
        if unc_k in log2:
            m = log2[unc_k] - log2[none_k]; v = ~np.isnan(m)
            desc.append((cond, float(np.median(m[v])), float(np.mean(np.abs(m[v]) > 1))))

    for n in tracks:
        sig["log2_" + n] = log2[n]
    sig["delta_mig_B10_minus_Parental_log2"] = d_mig
    sig.to_csv(args.out + "_perRegion.tsv", sep="\t", index=False)

    with open(args.out + "_summary.txt", "w") as fh:
        w = fh.write
        w("CTCF ANCHOR STABILITY TEST\n" + "=" * 60 + "\n")
        w(f"reference regions used: {len(sig)} (dropped {n_drop} for missing signal)\n\n")
        w("PRIMARY -- migration axis (untreated), B10_none vs Parental_none\n")
        w("  (clean: no UNC global-loss confound)\n")
        w(f"  median Δlog2 (B10 - Parental) : {med:+.4f}  ({2**med:.3f}x)\n")
        w(f"  IQR of Δlog2                  : [{q1:+.3f}, {q3:+.3f}]\n")
        w(f"  fraction |Δlog2| > 1 (2-fold) : {frac2:.4f}\n")
        w(f"  Cliff's delta (B10 vs Par)    : {cd:+.4f}\n")
        w(f"  Wilcoxon signed-rank p        : {w_p:.2e}  (n large -> read effect size, not p)\n\n")
        w(f"  DECISION RULE: |median|<{args.stable_med} AND frac2fold<{args.stable_frac}\n")
        w(f"  >>> {verdict} <<<\n\n")
        if desc:
            w("DESCRIPTIVE -- UNC axis within condition (none -> UNC)\n")
            w("  CONFOUNDED by uniform global loss on RPM tracks; NOT an invariance test.\n")
            w("  Absolute CTCF invariance under UNC requires IF/western, not this.\n")
            for cond, m, f2 in desc:
                w(f"  {cond:9s}: median Δlog2 = {m:+.4f} ({2**m:.3f}x), frac|Δ|>1 = {f2:.4f}\n")
            w("\n")
        w("INTERPRETATION\n")
        if stable:
            w("  Anchor region set effectively stable across migration. Ryder's differing\n")
            w("  B10 behavior is NOT explained by anchor remodeling; the B10 -0.29 residual\n")
            w("  is either real uniform UNC loss or a technical track effect -> IF/western.\n")
        else:
            w("  Anchor region set is itself remodeled across migration, so 'constitutive'\n")
            w("  is weaker on the B10 axis -- a principled reason to lean on library-size /\n")
            w("  greenlist (csaw) for B10 and treat Ryder's B10 anchor with caution.\n")
            w("  Absolute magnitude still requires IF/western.\n")

    order = list(tracks)
    plt.figure(figsize=(5, 4))
    for n in order:
        x = np.sort(log2[n][~np.isnan(log2[n])])
        plt.plot(x, np.linspace(0, 1, len(x)), label=n)
    plt.xlabel("log2(signal) at CTCF reference"); plt.ylabel("ECDF")
    plt.legend(frameon=False, fontsize=8); plt.tight_layout()
    plt.savefig(args.out + "_ecdf.pdf"); plt.close()

    plt.figure(figsize=(4.2, 4.2))
    plt.scatter(log2["Parental_none"], log2["B10_none"], s=1, alpha=.2)
    lo = np.nanmin([log2["Parental_none"].min(), log2["B10_none"].min()])
    hi = np.nanmax([log2["Parental_none"].max(), log2["B10_none"].max()])
    plt.plot([lo, hi], [lo, hi], "--", color="gray")
    plt.xlabel("log2 Parental_none"); plt.ylabel("log2 B10_none")
    plt.title(f"migration axis  (median Δ={med:+.3f})", fontsize=9)
    plt.tight_layout(); plt.savefig(args.out + "_pairedScatter.pdf"); plt.close()

    plt.figure(figsize=(5, 3.5))
    plt.hist(d_mig, bins=80); plt.axvline(0, color="gray", ls="--")
    plt.axvline(med, color="red", ls="-")
    plt.xlabel("Δlog2 (B10_none - Parental_none)"); plt.ylabel("CTCF regions")
    plt.title(f"frac |Δ|>1 = {frac2:.3f}", fontsize=9)
    plt.tight_layout(); plt.savefig(args.out + "_deltaHist.pdf"); plt.close()

    xs = np.linspace(-args.ext, args.ext, args.meta_bins)
    plt.figure(figsize=(5, 3.6))
    for n, col in [("Parental_none", "C3"), ("B10_none", "C0")]:
        mp = metaprofile(tracks[n], sig, args.ext, args.meta_bins)
        plt.plot(xs, mp, label=n, color=col)
    plt.xlabel("distance from CTCF center (bp)"); plt.ylabel("avg signal")
    plt.title("CTCF reference shape (untreated)", fontsize=9)
    plt.legend(frameon=False, fontsize=8); plt.tight_layout()
    plt.savefig(args.out + "_metaprofile.pdf"); plt.close()

    with open(args.out + "_MANIFEST.txt", "w") as fh:
        fh.write("test_ctcf_anchor_stability.py\n")
        fh.write(f"timestamp: {datetime.datetime.now().isoformat()}\n")
        fh.write(f"git_commit: {git_hash()}\n")
        fh.write(f"ref_bed: {args.ref}\n")
        for n, p in tracks.items():
            fh.write(f"track[{n}]: {p}\n")
        fh.write(f"pseudocount: {args.pseudocount}  exact: {args.exact}\n")
        fh.write(f"stable_med: {args.stable_med}  stable_frac: {args.stable_frac}\n")
        fh.write(f"n_regions_used: {len(sig)}  n_dropped: {n_drop}\n")

    print(f"{verdict}: median Δlog2 = {med:+.4f} ({2**med:.3f}x), "
          f"frac|Δ|>1 = {frac2:.4f}, Cliff's δ = {cd:+.4f}")


if __name__ == "__main__":
    main()