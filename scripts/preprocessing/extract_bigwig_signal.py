# =============================================================================
# extract_bigwig_signal.py
#
# PURPOSE
#   Build a per-region signal table from CUT&RUN BigWig tracks for the A375
#   parent-vs-bottom10 constricted-migration multi-omics project, for use as
#   ML FEATURES and for figures.
#
#   For every region in a BED file, extract one or more SUMMARY STATISTICS from
#   each input BigWig (default: all of mean, max, std, coverage; user may
#   restrict via --summary). Then, OPTIONALLY:
#     - if replicate bigWigs are supplied (names contain R1/R2/R3), average
#       replicates within each condition+mark group (per stat), and
#     - if --log2 is requested, compute migration + drug log2 contrasts (per stat).
#
# MULTI-STAT OUTPUT  --------------------------------------------------------
#   One column per (sample, stat) pair, named '<sample>_<stat>':
#       parental_h3k27me3_mean   parental_h3k27me3_max
#       parental_h3k27me3_std    parental_h3k27me3_coverage
#       b10_h3k27me3_mean   ... etc.
#
#   This follows the published precedent for tree-model compartment-prediction
#   (Zheng et al., iScience 2024, "CoRNN"): their RF baseline used mean+std per
#   region, and they reported that feeding the raw sub-bin signal HURT the RF
#   relative to summary stats. We extend their idea to {mean, max, std,
#   coverage}: mean=bulk signal, max=peak height, std=peakiness vs flatness,
#   coverage=fraction of region with any signal. Four complementary stats.
#
# WHY NOT SUM/MIN
#   pyBigWig also supports 'sum' and 'min'. For fixed-width bins, sum is mean
#   times a constant (redundant); min over a wide bin is near-zero almost
#   everywhere (no discriminative power). Both excluded.
#
# WHAT THE OUTPUT IS / IS NOT  ------------------------------------------------
#   These are DESCRIPTIVE FEATURES, not a statistical test — no p-value, no
#   replicate-variance model. Any claim that a region "changed" must come from
#   csaw / DiffBind. (Inputs are assumed already CPM-normalized — do NOT scale.)
#
# BED INPUT FLEXIBILITY  ------------------------------------------------------
#   Cols 1-3 (chr,start,end) are the computed interval. Cols 4+ are carried
#   through verbatim as extra_1.. (gene BED and bin BED both work).
#   region_id = col4 if present else 'chr:start-end' (stable join handle).
# =============================================================================

import pandas as pd
import numpy as np
import pyBigWig
import argparse
import os
import re


# Supported summary statistics (what pyBigWig's stats() accepts and what we
# consider biologically informative for tree-model features).
ALLOWED_STATS = ("mean", "max", "std", "coverage")


# -----------------------------------------------------------------------------
# 1. BED PARSING  (col1-3 computed on; col4 -> region_id; col5+ pass-through)
# -----------------------------------------------------------------------------
def extract_region_info(bed_file):
    """
    Read BED into a DataFrame.
      col1 = chr, col2 = start, col3 = end  (always parsed numerically)
      col4 (if present) = region_id (stable join handle, e.g. underscore form
            'chr1_750001_1000000' or a gene name)
      col5+ (if present) = extra_1, extra_2, ... carried through verbatim.
            This is where a second coordinate string (e.g. the colon form
            'chr1:750001-1000000', matching a Hi-C 'coordinate' column) lives
            for joining downstream.
    If col4 is absent, region_id is synthesised as 'chrom:start-end' and there
    are no extras.

    CHANGE FROM AN EARLIER VERSION: col4 used to be emitted as BOTH region_id
    AND extra_1, which duplicated the column and pushed the actually-useful
    join-key (e.g. the colon-form coordinate) into extra_2. That confused
    downstream auto-detection of the join key. Now col4 -> region_id only,
    and pass-through starts at col5.
    """
    rows = []
    max_extra = 0
    with open(bed_file) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            # col4 (parts[3]) -> region_id if present, else synthesised.
            # col5+ (parts[4:]) -> extras.
            if len(parts) >= 4:
                region_id = parts[3]
                extras = parts[4:]
            else:
                region_id = f"{chrom}:{start}-{end}"
                extras = []
            max_extra = max(max_extra, len(extras))
            rows.append((chrom, start, end, region_id, extras))

    records = []
    for chrom, start, end, region_id, extras in rows:
        rec = {"chr": chrom, "start": start, "end": end, "region_id": region_id}
        for i in range(max_extra):
            rec[f"extra_{i + 1}"] = extras[i] if i < len(extras) else "."
        records.append(rec)

    df = pd.DataFrame(records)
    extra_cols = [f"extra_{i + 1}" for i in range(max_extra)]
    region_cols = ["chr", "start", "end", "region_id"] + extra_cols
    return df, region_cols


# -----------------------------------------------------------------------------
# 2. SAMPLE-NAME HANDLING
# -----------------------------------------------------------------------------
def sanitize_name(name):
    """Hyphens->underscores; lowercase. 'B10-H3K27Me3-R1'->'b10_h3k27me3_r1'."""
    return name.replace("-", "_").lower()


def get_bw_name(bw_file):
    """Sample name from a bigWig path: prefix before first dot, sanitized.
    'B10_H3K27me3.mean.bw' -> 'b10_h3k27me3'."""
    basename = os.path.basename(bw_file)
    raw_name = basename.split(".")[0]
    return sanitize_name(raw_name)


# Replicate tag detector: trailing _r1/_r2/_r3.
REP_TAG = re.compile(r"_r[123]$")


def has_replicate_tag(sample_name):
    return bool(REP_TAG.search(sample_name))


def get_replicate_base(sample_name):
    return REP_TAG.sub("", sample_name)


def _stat_suffix_of(col):
    """Return the stat suffix of a column or None.
    'b10_h3k27me3_mean' -> 'mean'; 'b10_h3k27me3' -> None."""
    for s in ALLOWED_STATS:
        if col.endswith(f"_{s}"):
            return s
    return None


def _strip_stat_suffix(col):
    """Strip a trailing _<stat>; return col unchanged if no stat suffix."""
    s = _stat_suffix_of(col)
    return col[: -(len(s) + 1)] if s else col


# -----------------------------------------------------------------------------
# 3. SIGNAL EXTRACTION  (one summary value per region per BigWig)
# -----------------------------------------------------------------------------
def extract_signals_from_bw(bw_file, regions_df, summary="mean"):
    """
    One summary stat over each region [start,end) for a single bigWig.
    Coordinates clamped to chrom bounds; missing chroms / errors -> NaN.
    """
    bw = pyBigWig.open(bw_file)
    chroms = bw.chroms()
    values = []
    for _, row in regions_df.iterrows():
        chrom, start, end = row["chr"], int(row["start"]), int(row["end"])
        if chrom not in chroms:
            values.append(np.nan)
            continue
        clen = chroms[chrom]
        s, e = max(0, start), min(end, clen)
        if s >= e:
            values.append(np.nan)
            continue
        try:
            v = bw.stats(chrom, s, e, type=summary)[0]
            values.append(v if v is not None else np.nan)
        except Exception:
            values.append(np.nan)
    bw.close()
    return values


# -----------------------------------------------------------------------------
# 4. REPLICATE MEANS  (CONDITIONAL: only if replicate-tagged columns exist)
# -----------------------------------------------------------------------------
def add_replicate_means(df, replicate_cols):
    """
    Group replicate-tagged columns by (base_without_rep_tag, stat) and average
    them per stat. So _r1's mean averages with _r2's mean, _r1's max with
    _r2's max, etc. A group with only one replicate still gets averaged into
    a clean '<base>_<stat>' name for uniform downstream handling.
    Returns (df, [(group_base, mean_col), ...]).
    """
    # A column qualifies as "replicate-tagged" only if its non-stat base has _r<n>.
    tagged = [c for c in replicate_cols if has_replicate_tag(_strip_stat_suffix(c))]
    if not tagged:
        print("  [MEAN] no replicate-tagged bigWigs (R1/R2/R3) — skipping means.")
        return df, []

    groups = {}
    for col in tagged:
        stat = _stat_suffix_of(col)               # may be None if user opted out
        base_with_rep = _strip_stat_suffix(col)
        rep_base = get_replicate_base(base_with_rep)
        groups.setdefault((rep_base, stat), []).append(col)

    group_mean_cols = []
    for (rep_base, stat), cols in groups.items():
        mean_col = f"{rep_base}_{stat}" if stat else rep_base
        df.insert(max(df.columns.get_loc(c) for c in cols) + 1,
                  mean_col, df[cols].mean(axis=1))
        group_mean_cols.append((rep_base, mean_col))
        print(f"  [MEAN] {mean_col} = mean{tuple(cols)}")
    return df, group_mean_cols


# -----------------------------------------------------------------------------
# 5. LOG2 RATIO — AUTO contrasts on the biologically-correct axis  (only --log2)
# -----------------------------------------------------------------------------
# Builds migration (B10 vs Parental untreated) + drug (each treated vs its
# OWN-condition untreated) log2 contrasts, PER STAT. Column-name convention:
# log2ratio_<den>-<num>_<mark>_<stat>, e.g. log2ratio_p-b10_h3k27me3_mean.
# -----------------------------------------------------------------------------

def parse_track_name(col, conditions, treatments):
    """
    Split '<cond>_<mark>_<stat>' or '<cond>_<treat>_<mark>_<stat>' into
    (condition, treatment_or_None, mark, stat_or_None).
    """
    stat = _stat_suffix_of(col)
    base = _strip_stat_suffix(col) if stat else col

    for cond in conditions:
        if base == cond or base.startswith(cond + "_"):
            rest = base[len(cond):].lstrip("_")
            for tr in treatments:
                if rest == tr or rest.startswith(tr + "_"):
                    mark = rest[len(tr):].lstrip("_")
                    return (cond, tr, mark, stat)
            return (cond, None, rest, stat)
    return None


def add_log2_ratio(df, signal_cols, conditions=("b10", "parental"),
                   treatments=("unc", "epz"), migr_num="b10", migr_den="parental",
                   eps=None, eps_frac=0.01, abbrev=None):
    """
    Auto-build migration + drug log2 contrasts, per stat.
    """
    conditions = [c.lower() for c in conditions]
    treatments = [t.lower() for t in treatments]
    migr_num, migr_den = migr_num.lower(), migr_den.lower()
    if abbrev is None:
        abbrev = {"parental": "p"}
    ab = lambda s: abbrev.get(s, s)

    idx = {}
    for col in signal_cols:
        p = parse_track_name(col, conditions, treatments)
        if p is None:
            print(f"  [LOG2] (skip) '{col}' matches no known condition.")
            continue
        idx[p] = col

    log2_cols = []

    def resolve_eps(a, b):
        """Pseudocount for THIS contrast.

        If --eps was given, use it verbatim. Otherwise derive it from the data:
        eps_frac (default 1%) of the median signal pooled across the two tracks
        being compared. Computed per contrast, so each mark AND each stat gets a
        pseudocount on its own scale (a 'mean' track sits near 0.02; a 'coverage'
        track near 0.5 -- one global eps cannot suit both).

        WHY THIS MATTERS: a fixed eps=1.0 is only appropriate for raw counts.
        Against CPM-normalised bigWig means (~0.02 per 250kb bin) it swamps the
        signal, collapsing log2((a+1)/(b+1)) into a scaled DIFFERENCE and
        compressing every ratio toward zero. It is also not monotone in the
        ratio, so rank-based tests downstream see a reordered variable.
        """
        if eps is not None:
            return float(eps)
        pooled = np.concatenate([
            pd.to_numeric(a, errors="coerce").dropna().to_numpy(),
            pd.to_numeric(b, errors="coerce").dropna().to_numpy()])
        med = float(np.nanmedian(pooled)) if pooled.size else 0.0
        e = eps_frac * med
        if not np.isfinite(e) or e <= 0:
            e = 1e-6
            print(f"  [LOG2] (warn) median signal non-positive; eps fallback {e:g}")
        return e

    def emit(num_key, den_key, num_label, den_label, mark, stat):
        if num_key in idx and den_key in idx:
            a, b = df[idx[num_key]], df[idx[den_key]]
            stat_sfx = f"_{stat}" if stat else ""
            out_name = f"log2ratio_{den_label}-{num_label}_{mark}{stat_sfx}"
            e = resolve_eps(a, b)
            df[out_name] = np.log2((a + e) / (b + e))
            log2_cols.append(out_name)
            print(f"  [LOG2] {out_name} = log2(({idx[num_key]}+{e:.6g})/({idx[den_key]}+{e:.6g}))")
            return True
        return False

    # Migration axis (per stat, untreated only).
    untreated = sorted({(m, s) for (c, t, m, s) in idx if t is None})
    for mark, stat in untreated:
        emit((migr_num, None, mark, stat), (migr_den, None, mark, stat),
             num_label=ab(migr_num), den_label=ab(migr_den), mark=mark, stat=stat)

    # Drug axis (per stat, each treated vs its own untreated).
    treated = sorted({(c, t, m, s) for (c, t, m, s) in idx if t is not None})
    for (cond, tr, mark, stat) in treated:
        if not emit((cond, tr, mark, stat), (cond, None, mark, stat),
                    num_label=f"{ab(cond)}_{tr}", den_label=ab(cond),
                    mark=mark, stat=stat):
            print(f"  [LOG2] (skip) {cond}_{tr}_{mark}_{stat}: no untreated "
                  f"baseline.")

    if not log2_cols:
        print("[WARN] No log2 contrasts could be formed (check track names).")
    return df, log2_cols


# -----------------------------------------------------------------------------
# 6. OUTPUT
# -----------------------------------------------------------------------------
def save_output(df, output_dir, output_file):
    os.makedirs(output_dir, exist_ok=True)
    full_path = os.path.join(output_dir, output_file)
    df.to_csv(full_path, sep="\t", index=False)
    print(f"[INFO] Signal table -> {full_path}  ({df.shape[0]} x {df.shape[1]})")


# -----------------------------------------------------------------------------
# 7. PIPELINE
# -----------------------------------------------------------------------------
def main(bed_file, output_dir, bw_files, summary_stats=None,
         output_file="cut_run_signal_table.txt",
         conditions=("b10", "parental"), treatments=("unc", "epz"),
         migr_num="b10", migr_den="parental", eps=None, eps_frac=0.01,
         want_log2=False):
    """
    summary_stats : iterable of stat names from ALLOWED_STATS. Default (None)
                    expands to ALL four ('mean','max','std','coverage').
    """
    if summary_stats is None or len(summary_stats) == 0:
        summary_stats = list(ALLOWED_STATS)
    # Validate and dedup while preserving order.
    seen, validated = set(), []
    for s in summary_stats:
        s = s.lower()
        if s not in ALLOWED_STATS:
            raise ValueError(f"Unknown stat '{s}'. Allowed: {ALLOWED_STATS}")
        if s not in seen:
            seen.add(s); validated.append(s)
    summary_stats = validated

    if not bw_files:
        raise ValueError("No BigWig files provided")

    df, region_cols = extract_region_info(bed_file)
    print(f"[INFO] Loaded {len(df)} regions from {bed_file}")
    print(f"       Region/pass-through columns: {region_cols}")
    print(f"[INFO] Extracting these stats per bigWig: {summary_stats}")

    # Per-bigWig signal for each stat. Column = '<sample>_<stat>'.
    signal_cols = []
    for stat in summary_stats:
        print(f"[INFO] Stat = {stat}:")
        for bw_file in bw_files:
            name = get_bw_name(bw_file)
            col_name = f"{name}_{stat}"
            if col_name in df.columns:
                print(f"  [WARN] Duplicate column '{col_name}' "
                      f"({os.path.basename(bw_file)}) — skipping.")
                continue
            df[col_name] = extract_signals_from_bw(bw_file, df, summary=stat)
            signal_cols.append(col_name)
            print(f"  {os.path.basename(bw_file)} -> '{col_name}'")

    # Replicate means (per stat) — only if any column is replicate-tagged.
    print("[INFO] Replicate means (conditional):")
    df, group_mean_cols = add_replicate_means(df, signal_cols)

    compare_cols = ([m for _, m in group_mean_cols] if group_mean_cols
                    else signal_cols)

    if want_log2:
        print("[INFO] Log2 ratios (auto: migration + drug axes), per stat:")
        df, _ = add_log2_ratio(df, compare_cols, conditions=conditions,
                               treatments=treatments, migr_num=migr_num,
                               migr_den=migr_den, eps=eps, eps_frac=eps_frac)
    else:
        print("[INFO] --log2 not set; extracting signal only.")

    save_output(df, output_dir, output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract per-region BigWig signal as ML features. Computes "
                    "one or more SUMMARY STATS per region per bigWig (default: "
                    "all of mean,max,std,coverage). Replicate means and log2 "
                    "ratios computed automatically when applicable.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default: all four stats (mean,max,std,coverage), no log2:
  python extract_bigwig_signal.py bins_250kb.bed \\
      --bw data/processed/mean_bigwigs/*.mean.bw --output_dir results/

  # Restrict to just mean and std:
  python extract_bigwig_signal.py bins_250kb.bed \\
      --bw data/processed/mean_bigwigs/*.mean.bw --output_dir results/ \\
      --summary mean std

  # All stats + auto log2 contrasts:
  python extract_bigwig_signal.py bins_250kb.bed \\
      --bw data/processed/mean_bigwigs/*.mean.bw --output_dir results/ --log2
        """)
    parser.add_argument("bed_file", help="BED (bin or gene BED)")
    parser.add_argument("--bw", nargs="+", required=True, dest="bw_files",
                        help="One or more BigWig files (averaged or replicate)")
    parser.add_argument("--summary", nargs="+", default=None,
                        choices=list(ALLOWED_STATS),
                        help="Summary statistic(s) per region. Accepts one or "
                             "more of: mean max std coverage. Default: all four.")
    parser.add_argument("--output_dir", default=".", help="Output directory")
    parser.add_argument("--output_file", default="cut_run_signal_table.txt",
                        help="Filename for the table")
    parser.add_argument("--conditions", default="B10,Parental",
                        help="Comma-sep condition names (default: B10,Parental)")
    parser.add_argument("--treatments", default="UNC,EPZ",
                        help="Comma-sep treatment tokens (default: UNC,EPZ)")
    parser.add_argument("--migr_num", default="B10",
                        help="Migration numerator condition (default: B10)")
    parser.add_argument("--migr_den", default="Parental",
                        help="Migration denominator condition (default: Parental)")
    parser.add_argument("--eps", type=float, default=None,
                        help="FIXED pseudocount for the log2 ratio. Default: not set, "
                             "in which case eps is computed from the data per contrast "
                             "as --eps_frac of the median pooled signal (recommended). "
                             "Only pass this if you need a specific value. NOTE: 1.0 is "
                             "appropriate for raw counts but far too large for "
                             "CPM-normalised bigWig means (~0.02 per 250kb bin), where "
                             "it compresses every log2 ratio toward zero.")
    parser.add_argument("--eps_frac", type=float, default=0.01,
                        help="Data-driven pseudocount, as a fraction of the median "
                             "pooled signal of the two tracks in each contrast "
                             "(default: 0.01 = 1%%). Scales automatically per mark and "
                             "per stat. Ignored when --eps is given.")
    parser.add_argument("--log2", dest="want_log2", action="store_true",
                        help="Auto-compute migration + drug log2 contrasts")
    parser.set_defaults(want_log2=False)

    args = parser.parse_args()
    conditions = tuple(c.strip() for c in args.conditions.split(",") if c.strip())
    treatments = tuple(t.strip() for t in args.treatments.split(",") if t.strip())
    main(args.bed_file, args.output_dir, args.bw_files,
         summary_stats=args.summary,
         output_file=args.output_file,
         conditions=conditions, treatments=treatments,
         migr_num=args.migr_num, migr_den=args.migr_den,
         eps=args.eps, eps_frac=args.eps_frac, want_log2=args.want_log2)