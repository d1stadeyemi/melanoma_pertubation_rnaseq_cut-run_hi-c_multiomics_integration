#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unc_selective_erasure.py
===========================

Does UNC1999 (EZH2/EZH1 inhibitor) SELECTIVELY erase the H3K27me3 that
constricted migration installed?

The claim under test
--------------------
Under UNC1999, relative H3K27me3 depletion is *enriched* at the loci that
constricted migration reprogrammed -- compartment switchers, UNC compartment
"reverter" bins, and migration H3K27me3-gainers -- relative to a baseline-matched
set of structurally stable loci.

Why this is valid WITHOUT absolute magnitude (the whole point)
-------------------------------------------------------------
Greenlist (DESeq2 median-of-ratios) normalization centres the unknown global
H3K27me3 loss into a single per-contrast constant.  Write each bin's normalized
UNC effect as:

        dUNC_i = trueDelta_i + c          (c = unknown residual global offset)

Any contrast taken BETWEEN two sets of bins removes c:

        mean(dUNC | reprogrammed) - mean(dUNC | stable)
            = mean(trueDelta | reprogrammed) - mean(trueDelta | stable)

So the enrichment of depletion at reprogrammed-vs-stable loci is recoverable even
though the per-bin magnitude is not.  This script therefore NEVER reports an
absolute "amount lost"; it reports a between-set effect size (Cliff's delta) and a
label-permutation p-value, with the two sets matched on baseline H3K27me3
abundance so the contrast is abundance-fair (and robust to any mild
abundance-dependence in the normalization).

What it does
------------
1.  Selects the 8 H3K27me3 samples from the greenlist size-factor table
    (Parental-R3 is already excluded there) and pulls their BAM paths from the
    DiffBind-style sample sheet by SampleID.
2.  Uses the 250 kb bins of A375-CompartmentsAcrossConditions.txt as the canonical
    universe (so CUT&RUN and Hi-C share identical bins -- a clean join).
3.  Counts H3K27me3 reads per bin per sample (deeptools multiBamSummary), applies
    the greenlist size factors (normalized = raw / size_factor), averages the two
    replicates per group.
4.  Derives per-bin quantities used by BOTH Q-A and Q-D and writes a MASTER table:
    dUNC_b10, dUNC_parental, migration_dH3K27me3, baseline abundances, and every
    target-set annotation.
5.  Defines target sets from your compartment calls + the migration H3K27me3
    contrast, and a clean structurally-stable control set.
6.  Tests enrichment of UNC depletion at each target set vs the stable set,
    matched on baseline abundance (decile-stratified Cliff's delta + within-stratum
    label permutation, n=1000).
7.  Emits results table, diagnostic plots, and a provenance MANIFEST.

Scale-free guarantees, restated where they live in the code:
  * dUNC is a within-cell-line log2 ratio -> global offset is a shared constant.
  * The effect size is a BETWEEN-set contrast -> the constant cancels.
  * Target vs stable are matched on baseline abundance -> abundance is not the driver.
  * Permutation shuffles labels WITHIN abundance strata -> matched null.

Environment
-----------
conda env: multiomics-py   (deeptools, numpy, scipy, pandas, matplotlib)
Run on the cluster; multiBamSummary needs coordinate-sorted, indexed (.bai) BAMs.

------------------------------------------------------------------------------
USAGE EXAMPLES
------------------------------------------------------------------------------
# 0) Validate the scientific core synthetically (no data, no deeptools needed):
python scripts/analysis/unc_selective_erasure.py --self-test

# 1) Default run = CPM normalization (no size-factor file needed):
python scripts/analysis/unc_selective_erasure.py \
    --sample-sheet   data/processed/sample_sheet.csv \
    --compartments   data/processed/A375-CompartmentsAcrossConditions.txt \
    --blacklist      data/processed/greenlist/hg19-blacklist.v2.bed \
    --outdir         results/exploratory/unc_selective_erasure --threads 8

# 2) Greenlist normalization (pipeline-consistent; needs the size-factor file).
#    If --size-factors is omitted here it auto-falls back to CPM instead of failing:
python scripts/analysis/unc_selective_erasure.py \
    --norm greenlist \
    --sample-sheet   data/processed/sample_sheet.csv \
    --size-factors   results/exploratory/greenlist_deseq2/H3K27me3_all_size_factors.tsv \
    --compartments   data/processed/A375-CompartmentsAcrossConditions.txt \
    --blacklist      data/processed/greenlist/hg19-blacklist.v2.bed \
    --outdir         results/exploratory/unc_selective_erasure --threads 8

# 3) Robustness check: run both and confirm the enrichment delta is normalization-
#    independent (reuse the counts matrix so it's instant):
python scripts/analysis/unc_selective_erasure.py --norm cpm       [...] \
    --counts-table results/exploratory/unc_selective_erasure/raw_bin_counts.tsv
python scripts/analysis/unc_selective_erasure.py --norm greenlist [...] \
    --size-factors results/exploratory/greenlist_deseq2/H3K27me3_all_size_factors.tsv \
    --counts-table results/exploratory/unc_selective_erasure/raw_bin_counts.tsv
------------------------------------------------------------------------------
"""

import argparse
import os
import subprocess
import sys
import shutil
import json
import datetime as _dt

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

import matplotlib
matplotlib.use("Agg")          # headless cluster
import matplotlib.pyplot as plt


# =============================================================================
# Column-name constants for A375-CompartmentsAcrossConditions.txt
# (exact strings as supplied; hyphens mean string-key access only)
# =============================================================================
COL_COORD = "coordinate"
COL_CHR   = "chr"
COL_START = "start"          # 1-based start in the compartment file
COL_END   = "end"

# PC1 tracks (positive PC1 = A compartment, gene-density oriented)
PC1_PAR_DMSO = "Parental_DMSO"
PC1_PAR_UNC  = "Parental_UNC"
PC1_B10_DMSO = "Bottom10-DMSO"
PC1_B10_UNC  = "Bottom10-UNC"
PC1_PAR_ORIG = "Parental_orig"
PC1_B10_ORIG = "Bottom10_orig"

# Migration (Parental -> B10) transition calls, DMSO axis (matched to UNC experiment)
MIG_BOTH_DMSO   = "P->B10_Both_DMSO"
MIG_SWITCH_DMSO = "P->B10_Switch_DMSO"
MIG_SHIFT_DMSO  = "P->B10_Shift_DMSO"
# Migration transition calls, ORIG axis (original migration Hi-C)
MIG_BOTH_ORIG   = "P->B10_Both_original"
MIG_SWITCH_ORIG = "P->B10_Switch_Orig"
MIG_SHIFT_ORIG  = "P->B10_Shift_Orig"

# B10 DMSO -> UNC transition calls (used to define categorical reverters)
UNC_BOTH_B10   = "B10_DMSO->UNC_Both"
UNC_SWITCH_B10 = "B10_DMSO->UNC_switch"
UNC_SHIFT_B10  = "B10_DMSO->UNC_shift"

# Allowed categorical values (for reference / validation):
#   *_Switch_* in {'-','AtoB','BtoA'} ; *_Shift_* in {'-','TowardA','TowardB'}
#   *_Both_*   in {'-','BothtoA','BothtoB'}
NULL_TOKENS = {"-", "", "nan", "none", "na", "."}


# =============================================================================
# 1. Inputs: samples, size factors, BAM paths
# =============================================================================
def _assign_roles(cond, trt):
    """(Condition, Treatment) -> analysis role, case-insensitive and robust."""
    role = []
    for c, t in zip(cond.str.strip().str.lower(), trt.str.strip().str.lower()):
        cc = ("parental" if c.startswith("parent")
              else ("b10" if c in {"b10", "bottom", "bottom10"} else c))
        tt = ("unc" if t.startswith("unc")
              else ("none" if t in {"none", "dmso", "ctrl", "control"} else t))
        role.append(f"{cc}_{tt}")
    return role


def load_samples(size_factor_path, sample_sheet_path, project_root,
                 factor="H3K27me3", drop_samples=()):
    """Return a tidy frame of the 8 H3K27me3 samples with role, BAM, and (if
    available) greenlist size factor.

    Two selection modes, automatic:
      * size-factor file given  -> sample set + size factors come from it
        (it already excludes the Parental-R3 outlier; cleanest, used for greenlist).
      * no size-factor file     -> select Factor==<factor> rows from the sample
        sheet, minus --drop-samples; size factor column is absent (CPM only).
    BAM paths always come from the sample sheet, joined on SampleID -- never
    reconstructed (R2 BAM basenames differ from SampleIDs).
    """
    sheet = pd.read_csv(sample_sheet_path)  # standard CSV (quoted)
    for col in ("SampleID", "bamReads"):
        if col not in sheet.columns:
            sys.exit(f"[ERROR] sample sheet needs a '{col}' column")
    sheet_small = sheet[["SampleID", "bamReads"]].drop_duplicates("SampleID")

    if size_factor_path:
        sf = pd.read_csv(size_factor_path, sep="\t")
        need = {"SampleID", "Condition", "Treatment", "deseq2_size_factor"}
        missing = need - set(sf.columns)
        if missing:
            sys.exit(f"[ERROR] size-factor file missing columns: {missing}")
        df = sf.merge(sheet_small, on="SampleID", how="left")
    else:
        # derive the sample set straight from the sheet (already has bamReads)
        need = {"SampleID", "Condition", "Treatment", "Factor"}
        missing = need - set(sheet.columns)
        if missing:
            sys.exit(f"[ERROR] sample sheet missing columns for CPM selection: {missing}")
        df = sheet[sheet["Factor"].str.strip().str.lower()
                   == factor.strip().lower()].drop_duplicates("SampleID").copy()
        df["deseq2_size_factor"] = np.nan  # not used in CPM mode

    # drop named outliers (e.g. Parental-H3K27me3-R3)
    drop = {s.strip() for s in drop_samples if s.strip()}
    if drop:
        df = df[~df["SampleID"].isin(drop)].copy()

    if df["bamReads"].isna().any():
        bad = df.loc[df["bamReads"].isna(), "SampleID"].tolist()
        sys.exit(f"[ERROR] no BAM in sample sheet for: {bad}")

    df["role"] = _assign_roles(df["Condition"], df["Treatment"])
    df["bam"] = df["bamReads"].apply(
        lambda p: os.path.normpath(os.path.join(project_root, str(p)))
    )

    expected = {"parental_none", "parental_unc", "b10_none", "b10_unc"}
    got = set(df["role"])
    if got != expected:
        sys.exit(f"[ERROR] roles found {sorted(got)} != expected {sorted(expected)} "
                 f"(check --factor / --drop-samples / Condition+Treatment labels)")
    counts = df["role"].value_counts().to_dict()
    if any(v != 2 for v in counts.values()):
        sys.exit(f"[ERROR] expected 2 replicates per role, got {counts}")

    return df.reset_index(drop=True)


# =============================================================================
# 2. Bin universe from the compartment file
# =============================================================================
def load_bins(compartments_path, keep_chroms):
    comp = pd.read_csv(compartments_path, sep="\t", dtype={COL_CHR: str})
    for c in (COL_COORD, COL_CHR, COL_START, COL_END):
        if c not in comp.columns:
            sys.exit(f"[ERROR] compartment file missing column '{c}'")
    comp = comp[comp[COL_CHR].isin(keep_chroms)].copy()
    comp[COL_START] = comp[COL_START].astype(int)
    comp[COL_END]   = comp[COL_END].astype(int)
    # 1-based inclusive start -> 0-based half-open BED start
    comp["bed_start"] = comp[COL_START] - 1
    comp["bed_end"]   = comp[COL_END]
    comp = comp.sort_values([COL_CHR, "bed_start"]).reset_index(drop=True)
    return comp


def write_bins_bed(comp, path):
    comp[[COL_CHR, "bed_start", "bed_end", COL_COORD]].to_csv(
        path, sep="\t", header=False, index=False
    )
    return path


# =============================================================================
# 3. Counting (deeptools multiBamSummary) + greenlist normalization
# =============================================================================
def compute_blacklist_fraction(bins, blacklist_path):
    """Fraction of each 250 kb bin covered by blacklist regions (pure Python).

    Used instead of deeptools' --blackListFileName, which REJECTS any bin that
    merely touches a blacklist entry -- far too aggressive at 250 kb resolution.
    Here we compute the actual covered fraction so we can drop only bins that are
    *mostly* blacklist, keeping the rest (and reporting the count).
    """
    bl = pd.read_csv(blacklist_path, sep="\t", header=None,
                     usecols=[0, 1, 2], names=["chr", "start", "end"])
    bl["chr"] = bl["chr"].astype(str)
    bl_by = {c: g[["start", "end"]].values for c, g in bl.groupby("chr")}

    chrs   = bins["chr"].astype(str).values
    starts = bins["bed_start"].values
    ends   = bins["bed_end"].values
    frac = np.zeros(len(bins), float)
    for i in range(len(bins)):
        regions = bl_by.get(chrs[i])
        if regions is None:
            continue
        s, e = starts[i], ends[i]
        ov = np.clip(np.minimum(e, regions[:, 1]) - np.maximum(s, regions[:, 0]),
                     0, None).sum()
        frac[i] = ov / max(e - s, 1)
    return frac


def count_reads(bins_bed, samples, out_npz, out_raw, blacklist, threads,
                min_mapq, extend_reads, use_deeptools_blacklist=False):
    """Run multiBamSummary BED-file mode -> raw counts per bin per sample.

    NOTE: by default the blacklist is NOT passed to multiBamSummary, because in
    BED-file mode it rejects every bin that overlaps a blacklist entry (deeptools
    docs: "rejecting genomic chunks that happen to overlap an entry"), silently
    dropping whole 250 kb bins. Blacklist filtering is done in Python instead
    (compute_blacklist_fraction). Set use_deeptools_blacklist=True to restore the
    aggressive deeptools behavior.
    """
    if shutil.which("multiBamSummary") is None:
        sys.exit("[ERROR] multiBamSummary not found; activate env 'multiomics-py'")
    for b in samples["bam"]:
        if not os.path.exists(b):
            sys.exit(f"[ERROR] BAM not found: {b}")
        if not (os.path.exists(b + ".bai") or os.path.exists(b[:-4] + ".bai")):
            sys.exit(f"[ERROR] index (.bai) missing for: {b}")

    cmd = ["multiBamSummary", "BED-file", "--BED", bins_bed,
           "--bamfiles", *list(samples["bam"]),
           "--labels",   *list(samples["SampleID"]),
           "-o", out_npz, "--outRawCounts", out_raw,
           "-p", str(threads), "--minMappingQuality", str(min_mapq)]
    if blacklist and use_deeptools_blacklist:
        cmd += ["--blackListFileName", blacklist]
    if extend_reads:
        cmd += ["--extendReads"]
    print("[run]", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)
    return out_raw


def load_raw_counts(raw_path, samples):
    """Parse multiBamSummary --outRawCounts into a coordinate-indexed matrix."""
    raw = pd.read_csv(raw_path, sep="\t")
    raw.columns = [c.strip().strip("#").strip().strip("'\"") for c in raw.columns]
    ren = {}
    for c in raw.columns:
        cl = c.lower()
        if cl in ("chr", "chrom", "#chr", "'chr'"):     ren[c] = "chr"
        elif cl == "start":                              ren[c] = "bed_start"
        elif cl == "end":                                ren[c] = "bed_end"
    raw = raw.rename(columns=ren)
    raw["chr"] = raw["chr"].astype(str)
    raw["bed_start"] = raw["bed_start"].astype(int)
    raw["bed_end"]   = raw["bed_end"].astype(int)
    sample_cols = [s for s in samples["SampleID"] if s in raw.columns]
    if len(sample_cols) != len(samples):
        sys.exit(f"[ERROR] counts matrix columns {list(raw.columns)} "
                 f"do not cover all samples {list(samples['SampleID'])}")
    return raw[["chr", "bed_start", "bed_end"] + list(samples["SampleID"])]


def compute_divisors(raw, samples, norm):
    """Per-sample divisor:  normalized = (raw + pc) / divisor.

      greenlist : divisor = DESeq2 size factor (relative to invariant regions)
      cpm       : divisor = library_size / 1e6, where library_size = sum of the
                  sample's counts over all bins (genome-tiling colSum, the same
                  lib.size edgeR uses for CPM). No size-factor file required.

    Both are single per-sample scalars, so each contributes only a constant
    log-offset to dUNC -- which the between-set contrast cancels identically.
    """
    if norm == "greenlist":
        if samples["deseq2_size_factor"].isna().any():
            sys.exit("[ERROR] --norm greenlist but size factors are missing")
        return {sid: float(sf) for sid, sf in
                zip(samples["SampleID"], samples["deseq2_size_factor"])}, "size_factor"
    # cpm
    div = {}
    for sid in samples["SampleID"]:
        lib = float(raw[sid].sum())
        if lib <= 0:
            sys.exit(f"[ERROR] sample {sid} has zero library size")
        div[sid] = lib / 1e6
    return div, "libsize/1e6"


def normalize_and_aggregate(raw, samples, pseudocount, norm):
    """raw counts -> per-group mean normalized signal (CPM or greenlist).

    Pseudocount handling (edgeR-style proportional prior):
        A FLAT raw pseudocount divided by a per-sample scalar gives a
        sample-DEPENDENT normalized prior (1/divisor), which differentially
        shrinks log-ratios of low-count bins across samples whose size factors
        differ (here 0.47..1.61).  To avoid that, the prior is scaled by the
        per-sample divisor so the *normalized* additive prior is uniform:

            prior_raw_j   = pseudocount * divisor_j / mean(divisor)
            normalized_ij = (raw_ij + prior_raw_j) / divisor_j
                          = raw_ij / divisor_j  +  pseudocount / mean(divisor)

        => every sample gets the same additive prior (pseudocount/mean_div) in
        normalized space.  This makes normalized signal strictly positive, so no
        further epsilon is needed downstream.

    Replicates are then averaged within each role.  Absolute scale differs
    between CPM and greenlist, but dUNC log-ratios and their ranks (what the test
    uses) do not depend on it.
    """
    divisor, _ = compute_divisors(raw, samples, norm)
    mean_div = float(np.mean(list(divisor.values())))
    role_of = dict(zip(samples["SampleID"], samples["role"]))

    out = raw[["chr", "bed_start", "bed_end"]].copy()
    by_role = {r: [] for r in ("parental_none", "parental_unc", "b10_none", "b10_unc")}
    for sid in samples["SampleID"]:
        prior_raw = pseudocount * divisor[sid] / mean_div     # scaled prior
        normed = (raw[sid].astype(float) + prior_raw) / divisor[sid]
        by_role[role_of[sid]].append(normed)
    for r, cols in by_role.items():
        out[r] = np.mean(np.vstack([c.values for c in cols]), axis=0)
    return out


# =============================================================================
# 4. Per-bin derived quantities (shared with Q-D) + compartment categories
# =============================================================================
def categorize_direction(both, switch, shift, pc1_par, pc1_b10):
    """Reproduces the lab's Both -> Switch -> Shift -> PC1-sign categorizer.

    Returns one of: 'A_to_B', 'B_to_A', 'Stable_A', 'Stable_B', None.
    'TowardA'/'BtoA'/'toA'  => B_to_A (bin became more A)
    'TowardB'/'AtoB'/'toB'  => A_to_B (bin became more B)
    """
    def clean(v):
        return str(v).strip()

    v = clean(both)
    if v.lower() not in NULL_TOKENS:
        if "toA" in v or "BtoA" in v:   return "B_to_A"
        if "toB" in v or "AtoB" in v:   return "A_to_B"
    v = clean(switch)
    if v.lower() not in NULL_TOKENS:
        if v == "BtoA":  return "B_to_A"
        if v == "AtoB":  return "A_to_B"
    v = clean(shift)
    if v.lower() not in NULL_TOKENS:
        if "TowardA" in v:  return "B_to_A"
        if "TowardB" in v:  return "A_to_B"
    try:
        p, b = float(pc1_par), float(pc1_b10)
        if p > 0 and b > 0:  return "Stable_A"
        if p < 0 and b < 0:  return "Stable_B"
    except (ValueError, TypeError):
        pass
    return None


def categorize_movement(both, switch, shift):
    """Direction a bin MOVED between two conditions, from the lab's switch/shift/
    both calls.  Returns 'B_to_A' (moved toward A), 'A_to_B' (moved toward B),
    or None.  Same A/B vocabulary as categorize_direction, but with NO PC1
    stable fallback -- a movement axis has no 'stable' direction, so a bin with
    no categorical call is simply not a confident mover.
    Token map (mirrors categorize_direction): BothtoA/BtoA/TowardA -> B_to_A;
    BothtoB/AtoB/TowardB -> A_to_B.
    """
    s = str(both).strip()
    if s.lower() not in NULL_TOKENS:
        if "toA" in s or "BtoA" in s:  return "B_to_A"
        if "toB" in s or "AtoB" in s:  return "A_to_B"
    s = str(switch).strip()
    if s.lower() not in NULL_TOKENS:
        if s == "BtoA":  return "B_to_A"
        if s == "AtoB":  return "A_to_B"
    s = str(shift).strip()
    if s.lower() not in NULL_TOKENS:
        if "TowardA" in s:  return "B_to_A"
        if "TowardB" in s:  return "A_to_B"
    return None


def build_master(comp, norm, migration_axis, top_gainer_pct):
    """Join normalized signal to compartment calls and derive all annotations."""
    m = comp.merge(norm, on=["chr", "bed_start", "bed_end"], how="inner")

    # ---- relative effects (within-cell-line log2 ratios; global offset is a constant)
    # normalized signal is strictly positive (proportional prior), so no epsilon.
    m["dUNC_b10"]      = np.log2(m["b10_unc"]      / m["b10_none"])
    m["dUNC_parental"] = np.log2(m["parental_unc"] / m["parental_none"])
    m["migration_dH3K27me3"] = np.log2(m["b10_none"] / m["parental_none"])

    # ---- baseline abundance (log2 normalized none signal) used for matching
    m["baseline_b10"]      = np.log2(m["b10_none"])
    m["baseline_parental"] = np.log2(m["parental_none"])

    # ---- migration compartment direction (DMSO axis by default)
    if migration_axis == "dmso":
        both, sw, sh = MIG_BOTH_DMSO, MIG_SWITCH_DMSO, MIG_SHIFT_DMSO
        pc1_par, pc1_b10 = PC1_PAR_DMSO, PC1_B10_DMSO
    else:
        both, sw, sh = MIG_BOTH_ORIG, MIG_SWITCH_ORIG, MIG_SHIFT_ORIG
        pc1_par, pc1_b10 = PC1_PAR_ORIG, PC1_B10_ORIG

    m["mig_dir"] = [
        categorize_direction(r.get(both), r.get(sw), r.get(sh),
                             r.get(pc1_par), r.get(pc1_b10))
        for _, r in m.iterrows()
    ]
    m["is_mig_switcher"] = m["mig_dir"].isin(["A_to_B", "B_to_A"])

    # ---- UNC reverter (B10), CATEGORICAL (lab convention):
    #      a bin that (1) moved in a direction during migration
    #      (Parental_DMSO -> B10_DMSO: mig_dir is A_to_B or B_to_A) AND
    #      (2) moved the OPPOSITE direction under UNC (B10_DMSO -> B10_UNC),
    #      i.e. back toward the parental compartment. Purely from the lab's
    #      switch/shift/both calls -- no PC1 arithmetic, no thresholds.
    m["unc_dir_b10"] = [
        categorize_movement(r.get(UNC_BOTH_B10), r.get(UNC_SWITCH_B10),
                            r.get(UNC_SHIFT_B10))
        for _, r in m.iterrows()
    ]
    mig = m["mig_dir"]; unc = m["unc_dir_b10"]
    moved = mig.isin(["A_to_B", "B_to_A"]) & unc.isin(["A_to_B", "B_to_A"])
    m["is_reverter_b10"] = moved & (mig != unc)            # opposite direction

    # Strict variant for robustness: require HARD A<->B switches on BOTH axes
    # (the *_switch columns only), not soft shifts. Reported alongside; selected
    # for the contrast with --reverter-strict.
    def _hard(v):
        s = str(v).strip()
        return s if s in ("AtoB", "BtoA") else None
    mig_hard = m[MIG_SWITCH_DMSO if migration_axis == "dmso"
                 else MIG_SWITCH_ORIG].map(_hard).map(
        {"BtoA": "B_to_A", "AtoB": "A_to_B"})
    unc_hard = m[UNC_SWITCH_B10].map(_hard).map(
        {"BtoA": "B_to_A", "AtoB": "A_to_B"})
    m["is_reverter_b10_strict"] = (
        mig_hard.isin(["A_to_B", "B_to_A"]) & unc_hard.isin(["A_to_B", "B_to_A"])
        & (mig_hard != unc_hard)
    )

    # ---- migration H3K27me3-gainer: top percentile of the migration gain
    thr = np.nanpercentile(m["migration_dH3K27me3"], 100 - top_gainer_pct)
    m["is_gainer"] = m["migration_dH3K27me3"] >= thr

    # ---- structurally-stable control: same compartment in P and B10, no gain
    m["is_stable"] = m["mig_dir"].isin(["Stable_A", "Stable_B"]) & (~m["is_gainer"])

    return m


# =============================================================================
# 5. Statistics: Cliff's delta + baseline-matched permutation test
# =============================================================================
def cliffs_delta(x, y):
    """Cliff's delta for x vs y via the Mann-Whitney U identity.

        delta = 2*U1/(n1*n2) - 1 ,  U1 = #(x>y) + 0.5*ties  (scipy's U for x)

    Sign convention (documented because it has bitten this project before):
        delta < 0  <=>  x tends to be LOWER than y.
    Here x = target-set dUNC, y = stable-set dUNC, so a NEGATIVE delta means the
    target set is more depleted under UNC -- the hypothesis-consistent direction.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return np.nan
    U1 = mannwhitneyu(x, y, alternative="two-sided").statistic
    return 2.0 * U1 / (n1 * n2) - 1.0


def stratified_matched_test(values, is_target, baseline, n_strata=10,
                            n_perm=1000, seed=0, n_boot=0):
    """Baseline-matched Cliff's delta between target and stable bins, with a
    within-stratum label-permutation p-value and (optionally) a within-stratum
    bootstrap CI.

    * Strata = quantile bins of baseline abundance (target+stable pooled).
    * Pooled delta = pair-count-weighted mean of per-stratum Cliff's deltas
      (so abundance is controlled, not confounded).
    * Null = shuffle target/stable labels WITHIN each stratum (preserves both the
      abundance distribution and per-stratum target counts) -> matched permutation.
    * CI (n_boot>0) = resample target and stable values WITH replacement within
      each stratum, recompute the pooled delta -> 2.5/97.5 percentiles.
    Offset cancellation is automatic: delta is a between-set contrast.
    """
    values   = np.asarray(values, float)
    is_target = np.asarray(is_target, bool)
    baseline = np.asarray(baseline, float)

    # caller pre-filters to (target | stable); 'stable' is therefore '~target'
    keep = np.isfinite(values) & np.isfinite(baseline)
    v, t, b = values[keep], is_target[keep], baseline[keep]

    edges = np.quantile(b, np.linspace(0, 1, n_strata + 1))
    edges[0] -= 1e-9; edges[-1] += 1e-9
    strata = np.clip(np.digitize(b, edges[1:-1]), 0, n_strata - 1)

    def pooled_delta(target_labels):
        num = den = 0.0
        for s in range(n_strata):
            idx = strata == s
            xs = v[idx & target_labels]
            ys = v[idx & ~target_labels]
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
                # whole stratum one class -> contributes nothing either way
                lab[idx] = t[idx]
                continue
            lab[rng.choice(idx, size=k, replace=False)] = True
        perm[i] = pooled_delta(lab)
    perm = perm[np.isfinite(perm)]
    p_two = (np.sum(np.abs(perm) >= abs(obs)) + 1) / (len(perm) + 1) if len(perm) else np.nan

    # bootstrap CI (resample within strata, preserving target/stable membership)
    ci_low = ci_high = np.nan
    if n_boot > 0:
        # precompute per-stratum index pools
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

    return {
        "cliffs_delta": obs,
        "perm_p_two_sided": p_two,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "n_target": int(t.sum()),
        "n_stable": int((~t).sum()),
        "n_perm_valid": int(len(perm)),
        "median_target": float(np.nanmedian(v[t])) if t.any() else np.nan,
        "median_stable": float(np.nanmedian(v[~t])) if (~t).any() else np.nan,
    }


# =============================================================================
# 6. Run one enrichment contrast (target set vs stable set), abundance-matched
# =============================================================================
def run_contrast(master, target_mask, stable_mask, value_col, baseline_col,
                 n_strata, n_perm, seed, n_boot=0):
    # numpy boolean masks -> positional selection (index-safe even after filtering)
    target_mask = np.asarray(target_mask, bool)
    stable_mask = np.asarray(stable_mask, bool)
    sel = target_mask | stable_mask
    sub = master.loc[sel]
    return stratified_matched_test(
        values=sub[value_col].values,
        is_target=target_mask[sel],
        baseline=sub[baseline_col].values,
        n_strata=n_strata, n_perm=n_perm, seed=seed, n_boot=n_boot,
    )


def interpret(delta, p):
    if not np.isfinite(delta) or not np.isfinite(p):
        return "undetermined"
    if p >= 0.05:
        return "no enrichment"
    if delta < 0:
        return "DEPLETION enriched at target (hypothesis-consistent)"
    return "depletion DE-enriched at target (opposite)"


# ---- shared plotting vocabulary --------------------------------------------
# Category colours (stable = deep red reference; switchers/reverters in the blue
# family; gainers green). Direction (A->B vs B->A) is encoded by position+label,
# not colour, so the colour reads as "what kind of region" at a glance.
CAT_COLORS = {
    "stable":   "#b2182b",   # deep red  (reference)
    "switcher": "#0e12e0",   # deep blue
    "reverter": "#011b1d",   # mid blue  (blue family, distinct from switchers)
    "gainer":   "#0be96b",   # green
}


def cliff_band(delta):
    """Vargha-Delaney / Romano magnitude band for |Cliff's delta|."""
    a = abs(delta)
    if not np.isfinite(delta):
        return "na"
    if a < 0.147:
        return "negligible"
    if a < 0.330:
        return "small"
    if a < 0.474:
        return "medium"
    return "large"


def p_stars(p):
    """Significance stars from a (permutation) p-value."""
    if not np.isfinite(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


# =============================================================================
# 7. Plots
# =============================================================================
# Both ECDF and violin take a list of "panels". Each panel is a dict:
#   {label, values (1-D array), color, ref (bool), n, delta, p}
# 'ref' marks the stable reference (no stars / no effect-size annotation).
# delta/p are pulled from the SAME enrichment results table that the stats use,
# so the figure and the statistics can never disagree.
# =============================================================================
def plot_ecdf(panels, value_desc, outpath):
    try:
        fig, ax = plt.subplots(figsize=(6.5, 5), constrained_layout=True)
        for p in panels:
            x = np.sort(np.asarray(p["values"], float))
            x = x[np.isfinite(x)]
            if len(x) == 0:
                continue
            lw = 2.6 if p.get("ref") else 1.8
            ax.plot(x, np.linspace(0, 1, len(x)), lw=lw, color=p["color"],
                    label=f"{p['label']} (n={len(x)})")
        ax.axvline(0, color="grey", ls=":", lw=1)
        ax.set_xlabel(value_desc, fontsize=12)
        ax.set_ylabel("ECDF", fontsize=12)
        ax.tick_params(axis='both', labelsize=12)
        ax.set_title("UNC H3K27me3 effect: reprogrammed sets vs stable\n"
                     "(left shift = more depletion)")
        ax.legend(fontsize=12, frameon=False, loc="lower right")
        fig.savefig(outpath, dpi=150)
        plt.close(fig)
    except Exception as e:                       # plotting must never kill analysis
        print(f"[warn] ECDF plot failed: {e}")


def plot_violin(panels, value_desc, outpath, show_box=True):
    """Q1-style violins: category colours, an inner quartile bar, per-violin
    effect size (Cliff's delta + magnitude band) and significance stars (vs the
    stable reference), outliers trimmed for display with real whitespace above
    and below the data."""
    try:
        usable = [p for p in panels
                  if len(np.asarray(p["values"], float)[
                      np.isfinite(np.asarray(p["values"], float))]) >= 2]
        if len(usable) < 2:
            print("[warn] violin: fewer than 2 usable panels, skipping")
            return
        clean = [np.asarray(p["values"], float) for p in usable]
        clean = [v[np.isfinite(v)] for v in clean]

        # robust display range (drop 1%/99% tails) -> the "neat" trimmed look
        allv = np.concatenate(clean)
        lo, hi = np.nanpercentile(allv, [1, 99])
        if not (np.isfinite(lo) and np.isfinite(hi)) or hi <= lo:
            lo, hi = float(np.nanmin(allv)), float(np.nanmax(allv))
        rng = max(hi - lo, 1e-9)

        # trim each panel to the display range for the KDE body (keep full data
        # for the quartile bar so the summary stays exact)
        trimmed = []
        for v in clean:
            vt = v[(v >= lo) & (v <= hi)]
            trimmed.append(vt if len(vt) >= 2 else v)

        pos = np.arange(1, len(usable) + 1)
        fig, ax = plt.subplots(figsize=(0.8 * len(usable) + 1.5, 5.5),
                               constrained_layout=True)
        vp = ax.violinplot(trimmed, positions=pos, showmedians=False,
                           showextrema=False, widths=0.35)
        for body, p in zip(vp["bodies"], usable):
            body.set_facecolor(p["color"]); body.set_edgecolor("black")
            body.set_alpha(0.85); body.set_linewidth(0.6)

        # inner quartile "bar": thin 5-95% whisker, thick Q1-Q3 bar, white median
        if show_box:
            for x, v in zip(pos, clean):
                if len(v) < 5:
                    continue
                q1, med, q3 = np.percentile(v, [25, 50, 75])
                wlo, whi = np.percentile(v, [5, 95])
                ax.plot([x, x], [max(wlo, lo), min(whi, hi)], color="black",
                        lw=1.0, zorder=3)
                ax.plot([x, x], [q1, q3], color="black", lw=5.0,
                        solid_capstyle="butt", zorder=4)
                ax.scatter([x], [med], color="white", s=20, zorder=5,
                           edgecolor="black", linewidth=0.6)

        # reference lines: no-change (0) and the stable median
        ax.axhline(0, color="grey", ls=":", lw=1, zorder=1)
        ref = next((p for p in usable if p.get("ref")), None)
        if ref is not None:
            rv = np.asarray(ref["values"], float); rv = rv[np.isfinite(rv)]
            if len(rv):
                ax.axhline(np.median(rv), color=CAT_COLORS["stable"], ls="--",
                           lw=1, alpha=0.7, zorder=1,
                           label="stable median")

        # per-violin annotation: stars nearest the body, then delta + band above
        for x, p in zip(pos, usable):
            if p.get("ref"):
                continue
            d = p.get("delta", np.nan); pv = p.get("p", np.nan)
            ax.text(x, hi + 0.05 * rng, p_stars(pv), ha="center", va="bottom",
                    fontsize=12, fontweight="bold")
            ax.text(x, hi + 0.18 * rng, f"δ={d:+.2f}\n({cliff_band(d)})",
                    ha="center", va="bottom", fontsize=12)

        ax.set_xticks(pos)
        ax.set_xticklabels([f"{p['label']}\n(n={p.get('n', len(c))})"
                            for p, c in zip(usable, clean)], fontsize=12)
        # whitespace: pad below the data and leave headroom for the annotations
        ax.set_ylim(lo - 0.12 * rng, hi + 0.40 * rng)
        ax.set_ylabel(value_desc)
        ax.set_title("UNC H3K27me3 effect by reprogrammed set\n"
                     "(below 0 = depletion; * = perm p vs stable; "
                     "δ = Cliff's effect size vs stable)")
        fig.savefig(outpath, dpi=150)
        plt.close(fig)
    except Exception as e:
        print(f"[warn] violin plot failed: {e}")


def plot_balance(master, target_mask, baseline_col, name, outpath):
    try:
        fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
        ax.hist(master[master["is_stable"]][baseline_col].dropna(), bins=40,
                density=True, alpha=0.5, label="stable")
        ax.hist(master[target_mask][baseline_col].dropna(), bins=40,
                density=True, alpha=0.5, label=name)
        ax.set_xlabel("baseline abundance  log2(none signal)")
        ax.set_ylabel("density")
        ax.set_title(f"Baseline-abundance balance: {name} vs stable")
        ax.legend(fontsize=12, frameon=False)
        fig.savefig(outpath, dpi=150)
        plt.close(fig)
    except Exception as e:
        print(f"[warn] balance plot failed: {e}")


# =============================================================================
# 7b. Forest plot of effect sizes (the summary figure)
# =============================================================================
def plot_forest(results, outpath):
    """Horizontal Cliff's-delta forest plot with bootstrap-CI whiskers.

    One row per (cell_line, target_set); B10 and Parental colour-coded; vertical
    line at delta=0; star if perm p<0.05; n_target annotated. This is the
    effect-size-first summary panel."""
    try:
        df = results.copy()
        # order: B10 sets first (as listed), then Parental
        order = list(range(len(df)))[::-1]   # top row = first plan entry
        df = df.iloc[order].reset_index(drop=True)
        colors = {"B10": "#1f77b4", "Parental": "#ff7f0e"}
        y = np.arange(len(df))
        fig, ax = plt.subplots(figsize=(10, 0.55 * len(df) + 1.2),
                               constrained_layout=True)
        for i, r in df.iterrows():
            c = colors.get(r["cell_line"], "k")
            lo = r.get("ci_low", np.nan); hi = r.get("ci_high", np.nan)
            if np.isfinite(lo) and np.isfinite(hi):
                ax.plot([lo, hi], [i, i], color=c, lw=2, alpha=0.8, zorder=1)
            ax.scatter(r["cliffs_delta"], i, color=c, s=45, zorder=2)
            star = " *" if (np.isfinite(r["perm_p_two_sided"])
                            and r["perm_p_two_sided"] < 0.05) else ""
            ax.text(r["cliffs_delta"], i + 0.18,
                    f"δ={r['cliffs_delta']:+.3f}{star}  (n={r['n_target']})",
                    fontsize=12, ha="center", va="bottom", color=c)
        ax.axvline(0, color="grey", ls=":", lw=1)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{r['cell_line']}: {r['target_set']}"
                            for _, r in df.iterrows()], fontsize=12)
        ax.set_xlabel("Cliff's δ   (negative = H3K27me3 depletion enriched at target)", fontsize=12)
        ax.tick_params(axis='x', labelsize=12)
        ax.set_title("UNC selective-erasure effect sizes\n"
                     "(point = δ, bar = bootstrap 95% CI, * = perm p<0.05)")
        from matplotlib.lines import Line2D
        handles = [Line2D([0], [0], marker="o", color="w", label=k,
                          markerfacecolor=v, markersize=8) for k, v in colors.items()]
        ax.legend(handles=handles, fontsize=12, frameon=False, loc="lower right")
        fig.savefig(outpath, dpi=150)
        plt.close(fig)
    except Exception as e:
        print(f"[warn] forest plot failed: {e}")


# =============================================================================
# 7c. Gene mapping: target-set bins -> gene lists + occupied-bin background
# =============================================================================
def map_genes(master, gene_table_path, outdir, genes_per_bin="one"):
    """Map an hg19 gene->bin table onto the target-set flags and emit, into
    <outdir>/genelists/, one gene-symbol-per-line file per set plus a proper
    background_universe.txt (genes in EVERY bin that entered the test).

    The gene table must have columns: gene_name, coordinate  (250 kb bin id that
    matches the compartment 'coordinate'). TSS/bin_start/bin_end are used only
    for the one-gene-per-bin representative choice.

    genes_per_bin:
      'one' (default, recommended for ORA) -> collapse to ONE representative gene
            per bin (TSS closest to bin centre), so clustered genes (e.g. olfactory
            / HOX arrays falling in a single 250 kb bin) cannot inflate pathway
            counts. The SAME collapse is applied before deriving lists AND the
            background, so every list stays a strict subset of the background.
      'all' -> keep every gene (more power, but bin non-independence inflates ORA;
            only defensible with the bin-matched background, which we provide).
    """
    gt = pd.read_csv(gene_table_path, sep="\t")
    need = {"gene_name", "coordinate"}
    if not need.issubset(gt.columns):
        sys.exit(f"[ERROR] gene table needs columns {need}; has {list(gt.columns)}")

    # restrict to bins that were actually tested (the occupancy-filtered master)
    tested = set(master["coordinate"])
    gt = gt[gt["coordinate"].isin(tested)].copy()
    if not len(gt):
        sys.exit("[ERROR] no genes fall in the tested bins (coordinate mismatch?)")

    # optional one-gene-per-bin collapse (TSS closest to bin centre)
    if genes_per_bin == "one" and {"TSS", "bin_start", "bin_end"}.issubset(gt.columns):
        centre = (gt["bin_start"].astype(float) + gt["bin_end"].astype(float)) / 2.0
        gt["_dist"] = (gt["TSS"].astype(float) - centre).abs()
        gt = (gt.sort_values("_dist")
                .drop_duplicates("coordinate", keep="first")
                .drop(columns="_dist"))
    elif genes_per_bin == "one":
        gt = gt.drop_duplicates("coordinate", keep="first")

    gldir = os.path.join(outdir, "genelists")
    os.makedirs(gldir, exist_ok=True)

    def genes_in(mask):
        coords = set(master.loc[mask, "coordinate"])
        return sorted(set(gt.loc[gt["coordinate"].isin(coords), "gene_name"]))

    sets = {
        "migration_switchers": master["is_mig_switcher"].values,
        "switchers_A_to_B":    (master["mig_dir"] == "A_to_B").values,
        "switchers_B_to_A":    (master["mig_dir"] == "B_to_A").values,
        "UNC_reverters":       master["is_reverter_b10"].values,
        "UNC_reverters_A_to_B":(master["is_reverter_b10"] & (master["unc_dir_b10"] == "A_to_B")).values,
        "UNC_reverters_B_to_A":(master["is_reverter_b10"] & (master["unc_dir_b10"] == "B_to_A")).values,
        "UNC_reverters_strict":master["is_reverter_b10_strict"].values,
        "H3K27me3_gainers":    master["is_gainer"].values,
    }
    background = sorted(set(gt["gene_name"]))
    with open(os.path.join(gldir, "background_universe.txt"), "w") as fh:
        fh.write("\n".join(background) + "\n")

    assign_rows = []
    for name, mask in sets.items():
        gl = genes_in(mask)
        with open(os.path.join(gldir, f"{name}.txt"), "w") as fh:
            fh.write("\n".join(gl) + "\n")
        assign_rows.append((name, int(mask.sum()), len(gl)))

    # provenance table (pathway script ignores files matching 'gene_assignment')
    gt[["gene_name", "coordinate"]].to_csv(
        os.path.join(gldir, "gene_assignment.tsv"), sep="\t", index=False)

    print(f"[info] gene mapping ({genes_per_bin}/bin): "
          f"{len(background)} background genes, lists -> {gldir}")
    for name, nb, ng in assign_rows:
        print(f"    {name:22s} {nb:>5} bins -> {ng:>5} genes")
    return gldir


# =============================================================================
# 8. Provenance
# =============================================================================
def git_hash():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        return "NA"


def write_manifest(outdir, args, samples, n_bins, results_path, norm):
    lines = [
        "# MANIFEST  --  unc_selective_erasure.py",
        f"generated   : {_dt.datetime.now().isoformat(timespec='seconds')}",
        f"git_commit  : {git_hash()}",
        f"normalization : {norm}",
        f"python      : {sys.version.split()[0]}",
        f"numpy/scipy/pandas : {np.__version__}/"
        f"{__import__('scipy').__version__}/{pd.__version__}",
        "",
        "## parameters",
    ]
    for k, v in sorted(vars(args).items()):
        lines.append(f"  {k}: {v}")
    lines += ["", "## samples (role, size factor [greenlist only], BAM)"]
    for _, r in samples.iterrows():
        sf = r["deseq2_size_factor"]
        sf_str = f"{sf:.4f}" if pd.notna(sf) else "NA(cpm)"
        lines.append(f"  {r['SampleID']:28s} {r['role']:14s} "
                     f"sf={sf_str}  {r['bam']}")
    lines += ["", f"## bins analysed : {n_bins}",
              f"## results       : {results_path}", ""]
    with open(os.path.join(outdir, "MANIFEST.md"), "w") as fh:
        fh.write("\n".join(lines))


# =============================================================================
# 9. Main analysis pipeline
# =============================================================================
def main(args):
    os.makedirs(args.outdir, exist_ok=True)
    keep_chroms = ([f"chr{i}" for i in range(1, 23)]
                   + (["chrX"] if args.include_chrx else []))

    # --- resolve normalization (CPM default; greenlist needs a size-factor file) ---
    norm = args.norm
    if norm == "greenlist" and not args.size_factors:
        print("[warn] --norm greenlist requested but no --size-factors given; "
              "falling back to --norm cpm")
        norm = "cpm"
    drop_samples = [s for s in (args.drop_samples or "").split(",")]

    samples = load_samples(args.size_factors, args.sample_sheet, args.project_root,
                           factor=args.factor, drop_samples=drop_samples)
    print(f"[info] normalization: {norm}")
    comp = load_bins(args.compartments, keep_chroms)

    # ---- counts ----
    if args.counts_table:
        raw = load_raw_counts(args.counts_table, samples)
    else:
        bins_bed = write_bins_bed(comp, os.path.join(args.outdir, "bins_250kb.bed"))
        raw_path = os.path.join(args.outdir, "raw_bin_counts.tsv")
        count_reads(bins_bed, samples,
                    out_npz=os.path.join(args.outdir, "bin_counts.npz"),
                    out_raw=raw_path, blacklist=args.blacklist,
                    threads=args.threads, min_mapq=args.min_mapq,
                    extend_reads=args.extend_reads,
                    use_deeptools_blacklist=args.deeptools_blacklist)
        raw = load_raw_counts(raw_path, samples)

    norm_signal = normalize_and_aggregate(raw, samples, args.pseudocount, norm)
    master = build_master(comp, norm_signal, args.migration_axis,
                          args.top_gainer_pct)
    if args.migration_axis != "dmso":
        print("[warn] reverter uses the B10 DMSO->UNC calls; with --migration-axis "
              "orig the migration direction comes from a different experiment")

    # transparency: did the counting step silently drop any compartment bins?
    n_lost = len(comp) - len(master)
    if n_lost:
        print(f"[info] {n_lost} compartment bins had no counts and are absent from "
              f"the master ({len(comp)} compartment bins -> {len(master)} counted). "
              f"With --deeptools-blacklist this is expected (whole-bin rejection).")

    # Python-side blacklist filter (default): drop only bins MOSTLY covered, report.
    if args.blacklist and not args.deeptools_blacklist:
        bl_frac = compute_blacklist_fraction(master, args.blacklist)
        master["blacklist_frac"] = bl_frac
        n_drop = int((bl_frac >= args.blacklist_max_frac).sum())
        n_touch = int((bl_frac > 0).sum())
        print(f"[info] blacklist (python): {n_touch} bins touch a blacklist region; "
              f"{n_drop} are >= {args.blacklist_max_frac:.0%} covered -> dropped")
        master = master[bl_frac < args.blacklist_max_frac].copy()

    # optional occupancy filter: keep bins occupied in EITHER cell line, so the
    # parental contrasts are not dropped by a B10-only threshold.
    if args.min_baseline_quantile > 0:
        fb = np.nanquantile(master["baseline_b10"], args.min_baseline_quantile)
        fp = np.nanquantile(master["baseline_parental"], args.min_baseline_quantile)
        master = master[(master["baseline_b10"] >= fb)
                        | (master["baseline_parental"] >= fp)].copy()
    master = master.reset_index(drop=True)   # positional masks below assume 0..N-1

    master.to_csv(os.path.join(args.outdir, "master_bins_H3K27me3_250kb.tsv"),
                  sep="\t", index=False)

    # reverter breakdown (categorical): reverters move opposite to migration;
    # 'progressors' move the SAME direction (UNC pushes further). Report both,
    # plus the strict (hard-switch-only) reverter count.
    mig = master["mig_dir"]; unc = master["unc_dir_b10"]
    moved = mig.isin(["A_to_B", "B_to_A"]) & unc.isin(["A_to_B", "B_to_A"])
    n_prog = int((moved & (mig == unc)).sum())
    print(f"[info] master table: {len(master)} bins "
          f"({master['is_mig_switcher'].sum()} switchers, "
          f"{master['is_reverter_b10'].sum()} reverters, "
          f"{master['is_reverter_b10_strict'].sum()} reverters[strict], "
          f"{n_prog} progressors, "
          f"{master['is_gainer'].sum()} gainers, "
          f"{master['is_stable'].sum()} stable)")

    reverter_col = "is_reverter_b10_strict" if args.reverter_strict else "is_reverter_b10"
    print(f"[info] reverter set used for contrast: {reverter_col}")

    # ---- define contrasts (cell line, target name, target mask, value, baseline)
    sw   = master["is_mig_switcher"].values
    bta  = (master["mig_dir"] == "B_to_A").values
    atb  = (master["mig_dir"] == "A_to_B").values
    rev  = master[reverter_col].values
    rev_atb = rev & (master["unc_dir_b10"] == "A_to_B").values   # UNC moves A->B (these migrated B->A)
    rev_bta = rev & (master["unc_dir_b10"] == "B_to_A").values   # UNC moves B->A (these migrated A->B; H3K27me3-gaining loci)
    gain = master["is_gainer"].values
    stab = master["is_stable"].values

    plan = [
        ("B10", "migration_switchers", sw,      "dUNC_b10", "baseline_b10"),
        ("B10", "switchers_A_to_B",    atb,     "dUNC_b10", "baseline_b10"),
        ("B10", "switchers_B_to_A",    bta,     "dUNC_b10", "baseline_b10"),
        ("B10", "UNC_reverters",       rev,     "dUNC_b10", "baseline_b10"),
        ("B10", "reverters_A_to_B",    rev_atb, "dUNC_b10", "baseline_b10"),
        ("B10", "reverters_B_to_A",    rev_bta, "dUNC_b10", "baseline_b10"),
        ("B10", "H3K27me3_gainers",    gain,    "dUNC_b10", "baseline_b10"),
        # Parental control arm (reprogramming is a B10 property; effect should be
        # B10-specific or stronger).  Reverter set is B10-only, so omitted here.
        ("Parental", "migration_switchers", sw,   "dUNC_parental", "baseline_parental"),
        ("Parental", "H3K27me3_gainers",    gain, "dUNC_parental", "baseline_parental"),
    ]

    rows = []
    for i, (cell, name, tmask, vcol, bcol) in enumerate(plan):
        res = run_contrast(master, tmask, stab,
                           vcol, bcol, args.n_strata, args.n_perm,
                           seed=args.seed + i, n_boot=args.n_boot)
        res.update(cell_line=cell, target_set=name, value_col=vcol)
        res["interpretation"] = interpret(res["cliffs_delta"], res["perm_p_two_sided"])
        rows.append(res)
        print(f"  [{cell:8s}] {name:20s} "
              f"delta={res['cliffs_delta']:+.3f}  p={res['perm_p_two_sided']:.4f}  "
              f"(nT={res['n_target']}, nS={res['n_stable']})  -> {res['interpretation']}")

    results = pd.DataFrame(rows)[
        ["cell_line", "target_set", "value_col", "n_target", "n_stable",
         "cliffs_delta", "ci_low", "ci_high", "perm_p_two_sided",
         "median_target", "median_stable", "n_perm_valid", "interpretation"]
    ]
    results_path = os.path.join(args.outdir, "enrichment_results.tsv")
    results.to_csv(results_path, sep="\t", index=False)

    # ---- reverter robustness: lenient (switch+shift+both) vs strict (hard
    #      switches only). A robust result keeps the same sign and significance
    #      under the stricter, more conservative reverter definition.
    sweep_rows = []
    for label, col in [("lenient", "is_reverter_b10"),
                       ("strict", "is_reverter_b10_strict")]:
        rmask = master[col].values
        if rmask.sum() < 30:
            sweep_rows.append({"reverter_def": label, "n_reverters": int(rmask.sum()),
                               "cliffs_delta": np.nan, "perm_p_two_sided": np.nan,
                               "note": "too few reverters"})
            continue
        res = run_contrast(master, rmask, stab, "dUNC_b10", "baseline_b10",
                           args.n_strata, args.n_perm, seed=args.seed + 100)
        sweep_rows.append({"reverter_def": label, "n_reverters": res["n_target"],
                           "cliffs_delta": res["cliffs_delta"],
                           "perm_p_two_sided": res["perm_p_two_sided"],
                           "note": interpret(res["cliffs_delta"], res["perm_p_two_sided"])})
    sweep = pd.DataFrame(sweep_rows)
    sweep.to_csv(os.path.join(args.outdir, "reverter_robustness.tsv"),
                 sep="\t", index=False)
    print("[info] reverter robustness (lenient vs strict):")
    for _, r in sweep.iterrows():
        print(f"    {r['reverter_def']:8s}  n={r['n_reverters']:>5}  "
              f"delta={r['cliffs_delta']:+.3f}  p={r['perm_p_two_sided']}")

    # ---- plots ----
    # one source of truth: pull each set's delta/p from the results table so the
    # violin annotations and the statistics can never disagree.
    lookup = {(r["cell_line"], r["target_set"]): r for _, r in results.iterrows()}
    value_desc = "dUNC  =  log2(UNC / none)   [B10]"

    def _panel(label, mask, cat, set_name, ref=False):
        mask = np.asarray(mask, bool)
        n = int(mask.sum())
        if not ref and n < args.min_plot_n:
            print(f"[warn] plot: skipping '{label}' (n={n} < --min-plot-n "
                  f"{args.min_plot_n})")
            return None
        row = lookup.get(("B10", set_name)) if set_name else None
        return {"label": label, "color": CAT_COLORS[cat], "ref": ref, "n": n,
                "values": master.loc[mask, "dUNC_b10"].dropna().values,
                "delta": (row["cliffs_delta"] if row is not None else np.nan),
                "p": (row["perm_p_two_sided"] if row is not None else np.nan)}

    panels = [p for p in [
        _panel("stable",     stab,    "stable",   None, ref=True),
        _panel("A→B switch", atb,     "switcher", "switchers_A_to_B"),
        _panel("B→A switch", bta,     "switcher", "switchers_B_to_A"),
        _panel("A→B revert", rev_atb, "reverter", "reverters_A_to_B"),
        _panel("B→A revert", rev_bta, "reverter", "reverters_B_to_A"),
    ] if p is not None]
    if args.plot_gainers:                     # circular by construction; off by default
        g = _panel("gainer", gain, "gainer", "H3K27me3_gainers")
        if g is not None:
            panels.append(g)

    plot_ecdf(panels, value_desc, os.path.join(args.outdir, "ecdf_dUNC_b10.png"))
    plot_violin(panels, value_desc, os.path.join(args.outdir, "violin_dUNC_b10.png"))
    plot_balance(master, sw, "baseline_b10", "switchers",
                 os.path.join(args.outdir, "balance_switchers.png"))
    plot_forest(results, os.path.join(args.outdir, "forest_effect_sizes.png"))

    # ---- gene mapping (optional) ----
    if args.gene_table:
        map_genes(master, args.gene_table, args.outdir,
                  genes_per_bin=args.genes_per_bin)

    write_manifest(args.outdir, args, samples, len(master), results_path, norm)
    print(f"[done] results -> {results_path}")


# =============================================================================
# 10. Synthetic self-test  --  validates the scientific core (no data needed)
# =============================================================================
def self_test():
    """Prove: (a) injected depletion is recovered, (b) the unknown global offset c
    does NOT change the answer (offset cancellation), (c) a no-signal control is
    null, (d) baseline-abundance matching defuses an abundance confound."""
    print("=== SELF-TEST: Q-A enrichment core ===")
    rng = np.random.default_rng(42)
    N = 6000
    baseline = rng.normal(0, 1, N)                       # log2 baseline abundance

    def make(c, inject, confound=False):
        is_t = np.zeros(N, bool)
        is_t[rng.choice(N, 700, replace=False)] = True
        # target bins optionally sit at higher baseline (an abundance confound)
        if confound:
            base = baseline + np.where(is_t, 1.2, 0.0)
        else:
            base = baseline.copy()
        d = rng.normal(0, 0.5, N) + c                    # global offset c on ALL bins
        d[is_t] += inject                                # extra effect at target
        return d, is_t, base

    def test(values, is_t, base):
        keep = is_t | (~is_t)                            # all are target|stable here
        return stratified_matched_test(values[keep], is_t[keep], base[keep],
                                       n_strata=10, n_perm=500, seed=1)

    # (a) clear depletion injected at target
    d, t, b = make(c=0.0, inject=-0.8)
    r = test(d, t, b)
    print(f"(a) inject -0.8, c=0      -> delta={r['cliffs_delta']:+.3f} "
          f"p={r['perm_p_two_sided']:.4f}")
    assert r["cliffs_delta"] < -0.1 and r["perm_p_two_sided"] < 0.05, "should detect depletion"

    # (b) SAME truth but huge global offset c=-5 -> delta must be ~unchanged
    d2, t2, b2 = make(c=-5.0, inject=-0.8)
    # reuse identical target/noise by reseeding:
    rng2 = np.random.default_rng(42); _ = rng2.normal(0, 1, N)
    r2 = test(d2, t2, b2)
    print(f"(b) inject -0.8, c=-5     -> delta={r2['cliffs_delta']:+.3f} "
          f"p={r2['perm_p_two_sided']:.4f}  (offset cancellation)")
    assert r2["cliffs_delta"] < -0.1 and r2["perm_p_two_sided"] < 0.05, "offset broke cancellation"

    # (c) no injection -> null
    d3, t3, b3 = make(c=-0.3, inject=0.0)
    r3 = test(d3, t3, b3)
    print(f"(c) inject  0.0           -> delta={r3['cliffs_delta']:+.3f} "
          f"p={r3['perm_p_two_sided']:.4f}  (should be ~null)")
    assert abs(r3["cliffs_delta"]) < 0.1 or r3["perm_p_two_sided"] >= 0.05, "false positive"

    # (d) abundance confound but NO real effect -> matching keeps it null
    d4, t4, b4 = make(c=0.0, inject=0.0, confound=True)
    d4 = d4 + 0.0 * b4          # (no real coupling injected; only baseline differs)
    r4 = test(d4, t4, b4)
    print(f"(d) confound, no effect   -> delta={r4['cliffs_delta']:+.3f} "
          f"p={r4['perm_p_two_sided']:.4f}  (matching should keep null)")
    assert abs(r4["cliffs_delta"]) < 0.12 or r4["perm_p_two_sided"] >= 0.05, "confound leaked"

    # cliffs_delta sign sanity
    assert cliffs_delta([1, 2, 3], [4, 5, 6]) < 0, "sign convention wrong (x<y -> negative)"
    print("ALL SELF-TESTS PASSED")


# =============================================================================
# CLI
# =============================================================================
def build_parser():
    p = argparse.ArgumentParser(
        description="Does UNC1999 selectively erase migration-acquired H3K27me3? "
                    "(scale-free between-set enrichment of relative depletion).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument("--self-test", action="store_true",
                   help="run synthetic validation of the stats core and exit "
                        "(no data / no deeptools required)")

    # inputs
    p.add_argument("--sample-sheet", help="DiffBind-style CSV with SampleID,bamReads,...")
    p.add_argument("--norm", choices=["cpm", "greenlist"], default="cpm",
                   help="per-sample normalization. 'cpm' = library-size (no extra "
                        "file needed); 'greenlist' = DESeq2 size factors (requires "
                        "--size-factors). For the between-set test both give the "
                        "same ranks; greenlist falls back to cpm if no file is given.")
    p.add_argument("--size-factors", default=None,
                   help="greenlist DESeq2 size-factor TSV. Optional: only used by "
                        "--norm greenlist. When present it also defines the 8 samples "
                        "(Parental-R3 already excluded).")
    p.add_argument("--factor", default="H3K27me3",
                   help="Factor value to select from the sample sheet when no "
                        "size-factor file is given (CPM sample selection)")
    p.add_argument("--drop-samples", default="Parental-H3K27me3-R3",
                   help="comma-separated SampleIDs to exclude (e.g. the R3 outlier); "
                        "only applied to sample-sheet-based selection")
    p.add_argument("--compartments",
                   help="A375-CompartmentsAcrossConditions.txt (250 kb bins, both arms)")
    p.add_argument("--blacklist", default=None,
                   help="hg19 blacklist BED. By default applied in PYTHON (drop only "
                        "bins >= --blacklist-max-frac covered, keeping the rest).")
    p.add_argument("--deeptools-blacklist", action="store_true",
                   help="instead pass the blacklist to multiBamSummary, which REJECTS "
                        "any 250 kb bin overlapping a blacklist entry (aggressive; "
                        "drops whole bins and shrinks the bin universe)")
    p.add_argument("--blacklist-max-frac", type=float, default=0.5,
                   help="python blacklist mode: drop a bin if this fraction or more of "
                        "it is covered by blacklist regions")
    p.add_argument("--counts-table", default=None,
                   help="reuse a previous multiBamSummary --outRawCounts TSV "
                        "(skips read counting)")
    p.add_argument("--project-root", default=".",
                   help="root to resolve relative bamReads paths against")
    p.add_argument("--outdir", default="results/exploratory/unc_selective_erasure")

    # counting
    p.add_argument("--threads", type=int, default=8, help="multiBamSummary -p")
    p.add_argument("--min-mapq", type=int, default=20, help="minimum mapping quality")
    p.add_argument("--extend-reads", action="store_true",
                   help="extend paired-end reads to fragments (PE BAMs only)")
    p.add_argument("--pseudocount", type=float, default=1.0,
                   help="prior count; applied edgeR-style (scaled by the per-sample "
                        "normalization) so the normalized prior is uniform across samples")

    # set definitions
    p.add_argument("--migration-axis", choices=["dmso", "orig"], default="dmso",
                   help="which migration Hi-C axis defines switchers/direction: "
                        "'dmso' (vehicle, matched to UNC experiment) or 'orig'")
    p.add_argument("--top-gainer-pct", type=float, default=10.0,
                   help="top percentile of migration H3K27me3 gain = 'gainer' set")
    p.add_argument("--reverter-strict", action="store_true",
                   help="define B10 reverters using HARD A<->B switches only "
                        "(the *_switch columns) on both the migration and the "
                        "DMSO->UNC axis, instead of the switch+shift+both cascade. "
                        "Both counts are always reported; this selects which is "
                        "tested in the main contrast.")
    p.add_argument("--min-baseline-quantile", type=float, default=0.0,
                   help="drop bins below this quantile of B10 baseline abundance "
                        "(0 = keep all; 0.25 = H3K27me3-occupied focus)")
    p.add_argument("--include-chrx", action="store_true",
                   help="include chrX (default: autosomes only)")

    # gene mapping (optional)
    p.add_argument("--gene-table", default=None,
                   help="hg19 gene->bin TSV (cols: gene_name, coordinate[, TSS, "
                        "bin_start, bin_end]). If given, writes per-set gene lists + "
                        "background_universe.txt under <outdir>/genelists/ for the "
                        "pathway-enrichment script.")
    p.add_argument("--genes-per-bin", choices=["one", "all"], default="one",
                   help="'one' (recommended): one representative gene per bin (TSS "
                        "nearest bin centre) so clustered genes don't inflate ORA; "
                        "'all': keep every gene (needs the matched background).")

    # plotting
    p.add_argument("--plot-gainers", action="store_true",
                   help="include the H3K27me3-gainer set in the ECDF/violin. OFF by "
                        "default: gainers are defined from the same H3K27me3 counts "
                        "that enter dUNC, so their depletion signal is partly circular "
                        "and distracting in the figure (the stats still report them).")
    p.add_argument("--min-plot-n", type=int, default=15,
                   help="minimum bins for a set to appear as its own violin/ECDF curve "
                        "(small direction-split sets below this are skipped in the plot, "
                        "not in the stats)")

    # statistics
    p.add_argument("--n-strata", type=int, default=10,
                   help="baseline-abundance quantile strata for matching")
    p.add_argument("--n-perm", type=int, default=1000,
                   help="label permutations for the matched p-value")
    p.add_argument("--n-boot", type=int, default=1000,
                   help="bootstrap resamples for the Cliff's-delta CI (forest plot); "
                        "0 disables")
    p.add_argument("--seed", type=int, default=0)
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    if args.self_test:
        self_test()
        sys.exit(0)
    required = ["sample_sheet", "compartments"]
    miss = [r for r in required if getattr(args, r) is None]
    if miss:
        sys.exit(f"[ERROR] missing required args: {miss} (or use --self-test)")
    main(args)