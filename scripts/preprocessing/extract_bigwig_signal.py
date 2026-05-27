# =============================================================================
# extract_bigwig_signal.py
#
# PURPOSE
#   Build a per-region signal table from CUT&RUN BigWig tracks for the A375
#   parent-vs-bottom10 constricted-migration multi-omics project, for use as
#   ML FEATURES and for figures.
#
#   For every region in a BED file, extract a summary signal (default: mean)
#   from each input BigWig. Then, OPTIONALLY:
#     - if replicate bigWigs are supplied (names contain R1/R2/R3), average
#       replicates within each condition+mark group, and
#     - if --log2 is requested, compute the log2 ratio between two conditions
#       per mark.
#   If neither applies, the script simply extracts and writes the per-bigWig
#   mean signal (the common case when you feed it ALREADY-AVERAGED tracks from
#   average_replicate_bigwigs.sh — one track per mark, no replicates to merge).
#
# WHAT CHANGED FROM THE ORIGINAL (and why) ------------------------------------
#   (1) Replicate-mean is now CONDITIONAL: it runs only if bigWig names contain
#       a replicate tag (_r1/_r2/_r3). If you pass already-averaged tracks (no
#       replicate tag), there is nothing to average, so the step is skipped.
#   (2) RAW DIFFERENCE removed entirely. A raw difference of chromatin means is
#       not biologically meaningful here (dominated by the highest-signal bins,
#       no scale invariance). Only the LOG2 RATIO is offered, behind --log2.
#   (3) No separate "change-only" table is written (there is no diff to split
#       out). One table is produced: regions + signal (+ means/log2 if computed).
#
# WHAT THE OUTPUT IS / IS NOT -------------------------------------------------
#   These are DESCRIPTIVE FEATURES (mean signal; optional log2 ratio). They are
#   NOT a statistical test — no p-value, no replicate-variance model. Any claim
#   that a region "changed" must come from the replicate-aware csaw differential
#   pipeline, not from this table. (If inputs are already CPM-normalized
#   *.cpm.norm.bw or averaged tracks thereof, do NOT apply further scaling.)
#
# BED INPUT FLEXIBILITY
#   Cols 1-3 (chr,start,end) are the computed interval. Cols 4+ are carried
#   through verbatim as extra_1.. (a gene BED and a 5 kb-bin BED both work).
#   A region_id is synthesised (col4 if present, else chr:start-end) so every
#   row is joinable to RNA-seq / Hi-C downstream.
# =============================================================================

import pandas as pd
import numpy as np
import pyBigWig
import argparse
import os
import re


# -----------------------------------------------------------------------------
# 1. BED PARSING  (cols 1-3 computed on; cols 4+ opaque pass-through)
# -----------------------------------------------------------------------------
def extract_region_info(bed_file):
    """
    Read BED into a DataFrame. Fixed: col1=chr, col2=start, col3=end. Cols 4+
    are kept verbatim as extra_1.. (uninterpreted). region_id = col4 if present
    else 'chr:start-end', guaranteeing a stable joinable handle per row.
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
            extras = parts[3:]
            max_extra = max(max_extra, len(extras))
            region_id = extras[0] if extras else f"{chrom}:{start}-{end}"
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
    """Hyphens->underscores (valid pandas col), then lowercase so 'H3K27Me3'
    and 'H3K27me3' collapse to one group. 'B10-H3K27Ac-R1'->'b10_h3k27ac_r1'."""
    return name.replace("-", "_").lower()


def get_bw_name(bw_file):
    """Sample name from a bigWig path: prefix before the first dot, sanitized.
    'B10_none_H3K27me3.mean.bw' -> 'b10_none_h3k27me3'."""
    basename = os.path.basename(bw_file)
    raw_name = basename.split(".")[0]
    return sanitize_name(raw_name)


# Replicate tag detector: matches a trailing _r1/_r2/_r3 (case-insensitive,
# already lowercased by sanitize). Used both to DECIDE whether to average and
# to derive the group base by stripping the tag.
REP_TAG = re.compile(r"_r[123]$")


def has_replicate_tag(sample_name):
    """True if the sample name ends in _r1/_r2/_r3 (i.e. it is a replicate)."""
    return bool(REP_TAG.search(sample_name))


def get_replicate_base(sample_name):
    """Strip a trailing _r1/_r2/_r3 to get the condition+mark group base.
    'b10_h3k27ac_r1' -> 'b10_h3k27ac'. If no tag, returns the name unchanged."""
    return REP_TAG.sub("", sample_name)


# -----------------------------------------------------------------------------
# 3. SIGNAL EXTRACTION  (one summary value per region per BigWig)
# -----------------------------------------------------------------------------
def extract_signals_from_bw(bw_file, regions_df, summary="mean"):
    """
    One summary stat over each region [start,end) for a single bigWig.
    Coordinates clamped to chrom bounds; missing chroms / errors -> NaN.
    Replicate column is named by sample (no stat suffix).
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
    Group only the replicate-TAGGED columns by base name (strip _r<n>) and add a
    row-wise '<group>_mean' after each group's last replicate. Columns WITHOUT a
    replicate tag are ignored here (they are already single tracks). A group with
    only one replicate still gets a *_mean (equal to that replicate) so naming is
    uniform downstream. Returns (df, [(group, mean_col), ...]).

    NOTE: this is summary-then-average; deepTools bigwigAverage is average-then-
    summary. For descriptive features either is fine, but do not mix them in one
    table. (If you pre-averaged with average_replicate_bigwigs.sh, you will have
    no replicate tags here and this step is skipped — by design.)
    """
    tagged = [c for c in replicate_cols if has_replicate_tag(c)]
    if not tagged:
        print("  [MEAN] no replicate-tagged bigWigs (R1/R2/R3) — skipping means.")
        return df, []

    groups = {}
    for col in tagged:
        groups.setdefault(get_replicate_base(col), []).append(col)

    group_mean_cols = []
    for group, cols in groups.items():
        mean_col = f"{group}_mean"
        df.insert(max(df.columns.get_loc(c) for c in cols) + 1,
                  mean_col, df[cols].mean(axis=1))
        group_mean_cols.append((group, mean_col))
        print(f"  [MEAN] {mean_col} = mean{tuple(cols)}")
    return df, group_mean_cols


# -----------------------------------------------------------------------------
# 5. LOG2 RATIO between two conditions per mark  (descriptive; only if --log2)
# -----------------------------------------------------------------------------
def add_log2_ratio(df, signal_cols, cond_a="b10", cond_b="parental", eps=1.0):
    """
    Per mark present in BOTH conditions, compute:
        log2ratio_<cond_a>_<cond_b>_<mark> = log2((a+eps)/(b+eps))
    where a,b are that mark's signal columns for cond_a / cond_b.

    signal_cols here are the per-group signal columns to compare. When inputs are
    already-averaged tracks named '<cond>_<treat>_<mark>', the group base is the
    whole name; we match a leading '<cond>_' prefix to assign condition and take
    the remainder (treat_mark) as the mark key, so 'none' treatment is preserved
    in the key (b10_none_h3k27me3 vs parental_none_h3k27me3 compare correctly).

    eps (pseudocount, default 1.0, deepTools convention) avoids log(0) and damps
    near-zero noise. Positive value => HIGHER in cond_a (b10), matching the
    project's parent->bottom10 direction. This is a FEATURE, not a stat test.
    """
    map_a, map_b = {}, {}
    for col in signal_cols:
        if col.startswith(f"{cond_a}_"):
            map_a[col[len(cond_a) + 1:]] = col      # key = treat_mark
        elif col.startswith(f"{cond_b}_"):
            map_b[col[len(cond_b) + 1:]] = col

    keys = sorted(set(map_a) & set(map_b))
    log2_cols = []
    if not keys:
        print(f"[WARN] No mark shared between '{cond_a}' and '{cond_b}' for log2.")
        return df, log2_cols

    for k in keys:
        a, b = df[map_a[k]], df[map_b[k]]
        c = f"log2ratio_{cond_a}_{cond_b}_{k}"
        df[c] = np.log2((a + eps) / (b + eps))
        log2_cols.append(c)
        print(f"  [LOG2] {c} = log2(({map_a[k]}+{eps})/({map_b[k]}+{eps}))")
    return df, log2_cols


# -----------------------------------------------------------------------------
# 6. OUTPUT  (single table — no separate change-only table)
# -----------------------------------------------------------------------------
def save_output(df, output_dir, output_file):
    os.makedirs(output_dir, exist_ok=True)
    full_path = os.path.join(output_dir, output_file)
    df.to_csv(full_path, sep="\t", index=False)
    print(f"[INFO] Signal table -> {full_path}  ({df.shape[0]} x {df.shape[1]})")


# -----------------------------------------------------------------------------
# 7. PIPELINE
# -----------------------------------------------------------------------------
def main(bed_file, output_dir, bw_files, summary="mean",
         output_file="cut_run_signal_table.txt",
         cond_a="b10", cond_b="parental", eps=1.0, want_log2=False):
    if not bw_files:
        raise ValueError("No BigWig files provided")
    cond_a, cond_b = cond_a.lower(), cond_b.lower()

    # Step 1: regions
    df, region_cols = extract_region_info(bed_file)
    print(f"[INFO] Loaded {len(df)} regions from {bed_file}")
    print(f"       Region/pass-through columns: {region_cols}")

    # Step 2: per-bigWig signal extraction (dedup identical sanitized names)
    signal_cols = []
    for bw_file in bw_files:
        name = get_bw_name(bw_file)
        if name in df.columns:
            print(f"[WARN] Duplicate sample name '{name}' "
                  f"({os.path.basename(bw_file)}) — skipping.")
            continue
        df[name] = extract_signals_from_bw(bw_file, df, summary=summary)
        signal_cols.append(name)
        print(f"[INFO] {os.path.basename(bw_file)} -> '{name}'")

    # Step 3: replicate means — ONLY if any column is replicate-tagged (R1/2/3)
    print("[INFO] Replicate means (conditional):")
    df, group_mean_cols = add_replicate_means(df, signal_cols)

    # The columns to compare for log2: prefer group means if we made them,
    # otherwise the raw per-bigWig signal columns (already-averaged tracks).
    compare_cols = ([m for _, m in group_mean_cols] if group_mean_cols
                    else signal_cols)

    # Step 4: log2 ratio — ONLY if requested
    if want_log2:
        print(f"[INFO] Log2 ratio ({cond_a} vs {cond_b}, eps={eps}):")
        df, _ = add_log2_ratio(df, compare_cols, cond_a=cond_a, cond_b=cond_b, eps=eps)
    else:
        print("[INFO] --log2 not set; extracting mean signal only.")

    # Step 5: save (single table)
    save_output(df, output_dir, output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract per-region BigWig signal as ML features. Computes "
                    "replicate means ONLY if bigWig names are replicate-tagged "
                    "(R1/R2/R3); computes a log2 ratio between two conditions "
                    "ONLY with --log2. No raw-difference, single output table. "
                    "Cols 1-3 of the BED are parsed; cols 4+ pass through.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Already-averaged tracks (one per mark) -> just extract mean signal:
  python extract_bigwig_signal.py coordinates_5kb_bins.bed \\
      --bw data/processed/mean_bigwigs/*.mean.bw --output_dir results/

  # Same, but also add log2(B10/Parental) per mark:
  python extract_bigwig_signal.py coordinates_5kb_bins.bed \\
      --bw data/processed/mean_bigwigs/*.mean.bw --output_dir results/ \\
      --log2 --cond_a B10 --cond_b Parental

  # Replicate-tagged bigWigs -> replicate means computed automatically:
  python extract_bigwig_signal.py coordinates_5kb_bins.bed \\
      --bw bw/*-R*.bw --output_dir results/ --log2
        """)
    parser.add_argument("bed_file", help="BED (gene BED or coordinate-bin BED)")
    parser.add_argument("--bw", nargs="+", required=True, dest="bw_files",
                        help="One or more BigWig files (averaged or replicate)")
    parser.add_argument("--summary", default="mean",
                        choices=["mean", "max", "min", "sum"],
                        help="Per-region summary statistic (default: mean)")
    parser.add_argument("--output_dir", default=".", help="Output directory")
    parser.add_argument("--output_file", default="cut_run_signal_table.txt",
                        help="Filename for the table")
    parser.add_argument("--cond_a", default="B10",
                        help="Numerator condition for log2 (default: B10)")
    parser.add_argument("--cond_b", default="Parental",
                        help="Denominator condition for log2 (default: Parental)")
    parser.add_argument("--eps", type=float, default=1.0,
                        help="Pseudocount for log2 ratio (default: 1.0)")
    parser.add_argument("--log2", dest="want_log2", action="store_true",
                        help="Also compute log2 ratio between cond_a and cond_b")
    parser.set_defaults(want_log2=False)

    args = parser.parse_args()
    main(args.bed_file, args.output_dir, args.bw_files,
         summary=args.summary, output_file=args.output_file,
         cond_a=args.cond_a, cond_b=args.cond_b, eps=args.eps,
         want_log2=args.want_log2)