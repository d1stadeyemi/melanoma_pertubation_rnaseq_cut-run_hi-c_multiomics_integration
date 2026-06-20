#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
greenlist_deseq2_sizefactors.py

de Mello et al. 2024 (Briefings in Bioinformatics) EXACT greenlist protocol,
adapted to this project's sample sheet:
    1. quantify read pileup in greenlist regions with deeptools multiBamSummary
       (BED-file mode) -- the paper's exact counting tool;
    2. compute DESeq2 size factors = median-of-ratios (estimateSizeFactors,
       type="ratio") -- the paper's exact normalization;
    3. report size factors + outlier diagnostics, and OPTIONALLY cross-check
       against the csaw TMM-on-greenlist factors (normfactors_used.csv) and emit
       a csaw --normfactors injection file so the paper's factors can be run
       through your edgeR pipeline as a sensitivity arm.

Independent second implementation of greenlist normalization (DESeq2
median-of-ratios) to confirm the csaw greenlist run (TMM) from a separate
codebase and to let you cite the paper's exact method.

Size-factor math is bit-faithful to DESeq2::estimateSizeFactorsForMatrix:
  loggeomeans_i = mean_j log(count_ij)            # region i geometric mean (log)
  region used   <=> loggeomeans_i is finite        # i.e. NO zero in ANY sample
  sf_j          = exp( median_i [ log(count_ij) - loggeomeans_i ] )  over used i

CONVENTION NOTE (edgeR bridge): a DESeq2 size factor s_j is a standalone divisor
(normalized = raw / s_j). edgeR uses effective.lib = lib.size * norm.factor as the
divisor. So to inject de Mello factors into csaw --normfactors you must pass
  nf_j = s_j / lib.size_j   (csaw then centers it; effective.lib becomes ∝ s_j).
Passing the raw s_j would double-count library size. This script emits nf_j for you.
"""

import argparse, os, sys, subprocess, datetime, tempfile, shutil
import numpy as np
import pandas as pd

CANON = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sheet", required=True, help="sample sheet CSV (SampleID,Factor,Condition,Treatment,Replicate,Batch,bamReads,...)")
    p.add_argument("--mark", required=True, help="Factor to select, case-insensitive (e.g. H3K27me3)")
    p.add_argument("--greenlist", required=True, help="greenlist BED in your genome build (hg19)")
    p.add_argument("--out", required=True, help="output prefix")
    p.add_argument("--samples", default=None,
                   help="optional comma-separated SampleIDs to restrict to (e.g. drop an outlier replicate). Default: all rows of --mark.")
    p.add_argument("--counts", default=None,
                   help="optional precomputed multiBamSummary --outRawCounts table; skips running multiBamSummary.")
    p.add_argument("--csaw-normfactors", default=None,
                   help="optional csaw <tag>_normfactors_used.csv for cross-check + exact lib.size (needs nf_greenlist + lib.size columns).")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--min-mapq", type=int, default=20, help="multiBamSummary --minMappingQuality (match csaw minq=20). [20]")
    p.add_argument("--keep-duplicates", action="store_true",
                   help="do NOT pass --ignoreDuplicates (default ignores dups, matching csaw dedup=TRUE).")
    p.add_argument("--extend-reads", action="store_true", help="pass --extendReads to multiBamSummary.")
    p.add_argument("--canonical-only", action="store_true", default=True)
    p.add_argument("--tmp", default=None, help="scratch dir for multiBamSummary npz (default: system tmp).")
    return p.parse_args()


def git_hash():
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"],
                                       stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "NA"


def load_sheet(path, mark, samples):
    ss = pd.read_csv(path)
    req = {"SampleID", "Factor", "Condition", "Treatment", "bamReads"}
    miss = req - set(ss.columns)
    if miss:
        sys.exit(f"sheet missing columns: {miss}")
    ss = ss[ss.Factor.str.strip().str.lower() == mark.strip().lower()].copy()
    if samples:
        keep = [s.strip() for s in samples.split(",")]
        ss = ss[ss.SampleID.isin(keep)]
        missing = set(keep) - set(ss.SampleID)
        if missing:
            sys.exit(f"--samples not found for mark {mark}: {missing}")
    if ss.empty:
        sys.exit(f"no samples for mark {mark}")
    bad = ss[~ss.bamReads.apply(os.path.exists)]
    if len(bad):
        sys.exit("missing BAM(s):\n" + "\n".join(bad.bamReads))
    return ss.reset_index(drop=True)


def run_multibamsummary(ss, greenlist, args, tmpdir):
    npz = os.path.join(tmpdir, "mbs.npz")
    raw = os.path.join(tmpdir, "mbs_counts.tab")
    cmd = ["multiBamSummary", "BED-file", "--BED", greenlist,
           "--bamfiles", *ss.bamReads.tolist(),
           "--labels", *ss.SampleID.tolist(),
           "-o", npz, "--outRawCounts", raw,
           "-p", str(args.threads),
           "--minMappingQuality", str(args.min_mapq)]
    if not args.keep_duplicates:
        cmd.append("--ignoreDuplicates")
    if args.extend_reads:
        cmd.append("--extendReads")
    print("[multiBamSummary] " + " ".join(cmd[:6]) + " ... ("
          + str(len(ss)) + " bams)", flush=True)
    subprocess.run(cmd, check=True)
    return raw


def load_counts(raw, sample_order):
    # multiBamSummary --outRawCounts: header like  #'chr'  'start'  'end'  'lab1' ...
    df = pd.read_csv(raw, sep="\t")
    df.columns = [c.strip().strip("#").strip().strip("'").strip('"') for c in df.columns]
    coord = df.columns[:3].tolist()
    mat = df.drop(columns=coord)
    mat.columns = [c.strip().strip("'").strip('"') for c in mat.columns]
    missing = set(sample_order) - set(mat.columns)
    if missing:
        sys.exit(f"counts table missing sample columns: {missing}\nhave: {list(mat.columns)}")
    mat = mat[sample_order]
    region_id = df[coord[0]].astype(str) + ":" + df[coord[1]].astype(str) + "-" + df[coord[2]].astype(str)
    return mat.to_numpy(dtype=float), region_id.to_numpy(), df[coord[0]].astype(str).to_numpy()


def deseq2_size_factors(counts):
    """Exact DESeq2 estimateSizeFactorsForMatrix (type='ratio'). counts: regions x samples."""
    with np.errstate(divide="ignore"):
        logc = np.log(counts)               # log(0) -> -inf
    loggeo = logc.mean(axis=1)              # any zero in a row -> -inf
    used = np.isfinite(loggeo)
    n_used = int(used.sum())
    if n_used == 0:
        sys.exit("No greenlist region is non-zero in ALL samples -> cannot compute size factors. "
                 "Greenlist coverage too sparse for this sample set.")
    logratio = logc[used] - loggeo[used][:, None]   # (n_used x samples)
    sf = np.exp(np.median(logratio, axis=0))
    ratios = np.exp(logratio)                        # per-region ratio to geomean, per sample
    return sf, n_used, used, ratios


def main():
    a = parse_args()
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    ss = load_sheet(a.sheet, a.mark, a.samples)
    order = ss.SampleID.tolist()

    own_tmp = a.tmp is None
    tmpdir = a.tmp or tempfile.mkdtemp(prefix="greenlist_mbs_")
    os.makedirs(tmpdir, exist_ok=True)
    try:
        raw = a.counts if a.counts else run_multibamsummary(ss, a.greenlist, a, tmpdir)
        counts, region_id, chroms = load_counts(raw, order)
    finally:
        if own_tmp and not a.counts and os.path.isdir(tmpdir):
            shutil.rmtree(tmpdir, ignore_errors=True)

    if a.canonical_only:
        keep = np.isin(chroms, list(CANON))
        counts, region_id = counts[keep], region_id[keep]

    sf, n_used, used, ratios = deseq2_size_factors(counts)
    sfac = pd.Series(sf, index=order)

    # ---- per-sample outlier diagnostics (the manual "check outliers" step, automated)
    diag = []
    for j, s in enumerate(order):
        r = ratios[:, j]
        lr2 = np.log2(r / sf[j])                     # log2 deviation of each region from this sample's SF
        diag.append(dict(SampleID=s, size_factor=sf[j],
                         ratio_median=np.median(r), ratio_MAD=np.median(np.abs(r - np.median(r))),
                         frac_2fold_outlier=float(np.mean(np.abs(lr2) > 1))))
    diag = pd.DataFrame(diag)

    out = pd.DataFrame({"SampleID": order})
    for c in ["Condition", "Treatment", "Replicate", "Batch"]:
        if c in ss.columns:
            out[c] = ss[c].values
    out["deseq2_size_factor"] = sf
    out["n_regions_used"] = n_used

    # ---- optional cross-check vs csaw TMM greenlist + edgeR-injection column
    xcheck_txt = []
    if a.csaw_normfactors:
        cz = pd.read_csv(a.csaw_normfactors)
        if "lib.size" not in cz.columns or "nf_greenlist" not in cz.columns:
            xcheck_txt.append("csaw file lacks lib.size/nf_greenlist; cross-check skipped.")
        else:
            cz = cz.set_index("SampleID")
            common = [s for s in order if s in cz.index]
            if len(common) < len(order):
                xcheck_txt.append(f"WARNING: only {len(common)}/{len(order)} samples matched csaw file by SampleID.")
            lib = cz.loc[common, "lib.size"].astype(float).values
            csaw_nf = cz.loc[common, "nf_greenlist"].astype(float).values
            s_common = sfac.loc[common].values
            # de Mello -> edgeR norm.factor: nf = (s/lib) centered to geomean 1
            nf_raw = s_common / lib
            nf_dm = nf_raw / np.exp(np.mean(np.log(nf_raw)))
            csaw_nf_c = csaw_nf / np.exp(np.mean(np.log(csaw_nf)))
            pct = 100 * np.abs(nf_dm - csaw_nf_c) / csaw_nf_c
            rho = np.corrcoef(nf_dm, csaw_nf_c)[0, 1] if len(common) > 1 else float("nan")
            inj = pd.DataFrame({"SampleID": common, "SF": nf_dm, "QC": "pass"})
            inj.to_csv(a.out + "_csaw_normfactors.tsv", sep="\t", index=False)
            cmp = pd.DataFrame({"SampleID": common, "csaw_nf_greenlist": csaw_nf_c,
                                "deMello_nf": nf_dm, "pct_diff": pct})
            cmp.to_csv(a.out + "_crosscheck.tsv", sep="\t", index=False)
            out = out.merge(inj[["SampleID", "SF"]].rename(columns={"SF": "edger_normfactor"}),
                            on="SampleID", how="left")
            xcheck_txt += [
                f"cross-check vs csaw TMM greenlist (both centered to geomean 1):",
                f"  Pearson r           = {rho:.4f}",
                f"  max |%diff|         = {pct.max():.2f}%   median |%diff| = {np.median(pct):.2f}%",
                "  -> AGREE (independent DESeq2 and csaw-TMM greenlist factors concur)"
                if pct.max() < 10 else
                "  -> DIVERGE >10% on some sample; inspect counts/coverage before trusting.",
                f"  injection file -> {os.path.basename(a.out)}_csaw_normfactors.tsv  (use: csaw.R --normfactors <this> --norm sf)",
            ]

    out.to_csv(a.out + "_size_factors.tsv", sep="\t", index=False)
    diag.to_csv(a.out + "_diagnostics.tsv", sep="\t", index=False)

    with open(a.out + "_diagnostics.txt", "w") as fh:
        fh.write("de Mello greenlist DESeq2 size factors\n" + "=" * 50 + "\n")
        fh.write(f"mark: {a.mark}   samples: {len(order)}   greenlist regions used (non-zero in all): {n_used}\n\n")
        fh.write(diag.to_string(index=False) + "\n\n")
        if n_used < 100:
            fh.write(f"CAUTION: only {n_used} usable greenlist regions -> median-of-ratios is noisy. "
                     "Check liftOver yield / coverage.\n\n")
        if xcheck_txt:
            fh.write("\n".join(xcheck_txt) + "\n")

    # ---- plots
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].bar(range(len(order)), sf, color="#2A9D8F")
    ax[0].set_xticks(range(len(order))); ax[0].set_xticklabels(order, rotation=90, fontsize=7)
    ax[0].axhline(1, color="gray", ls="--"); ax[0].set_ylabel("DESeq2 size factor")
    ax[0].set_title(f"greenlist size factors (n_regions={n_used})", fontsize=9)
    ax[1].boxplot([np.log2(ratios[:, j] / sf[j]) for j in range(len(order))],
                  showfliers=False)
    ax[1].axhline(0, color="gray", ls="--")
    ax[1].set_xticks(range(1, len(order) + 1))
    ax[1].set_xticklabels(order, rotation=90, fontsize=7)
    ax[1].set_ylabel("log2(region ratio / size factor)")
    ax[1].set_title("per-region spread around SF (outlier check)", fontsize=9)
    plt.tight_layout(); plt.savefig(a.out + "_sizefactors.pdf"); plt.close()

    with open(a.out + "_MANIFEST.txt", "w") as fh:
        fh.write("greenlist_deseq2_sizefactors.py (de Mello protocol)\n")
        fh.write(f"timestamp: {datetime.datetime.now().isoformat()}\n")
        fh.write(f"git_commit: {git_hash()}\n")
        fh.write(f"sheet: {a.sheet}\nmark: {a.mark}\ngreenlist: {a.greenlist}\n")
        fh.write(f"samples: {','.join(order)}\n")
        fh.write(f"min_mapq: {a.min_mapq}  ignore_duplicates: {not a.keep_duplicates}  extend_reads: {a.extend_reads}\n")
        fh.write(f"n_regions_used: {n_used}\n")
        if a.csaw_normfactors:
            fh.write(f"crosscheck_against: {a.csaw_normfactors}\n")

    print("\nsize factors:")
    print(sfac.to_string())
    print(f"\nregions used: {n_used}")
    if xcheck_txt:
        print("\n" + "\n".join(xcheck_txt))
    print(f"\nwrote: {a.out}_size_factors.tsv (+ diagnostics, pdf, MANIFEST"
          + (", crosscheck, csaw_normfactors" if a.csaw_normfactors else "") + ")")


if __name__ == "__main__":
    main()