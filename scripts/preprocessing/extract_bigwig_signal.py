#!/lustre/isaac24/scratch/majagbe/envs/epp604_env/bin/python3

import pandas as pd
import numpy as np
import pyBigWig
import argparse
import os
import re

def extract_region_info(bed_file):
    """
    Read a BED file and return a DataFrame with one row per genomic region.
    Expects at least 3 tab-separated columns: chr, start, end.
    Optional columns: gene_id (col4), gene_name (col5), strand (col6).
    """
    results = []
    with open(bed_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom     = parts[0]
            start     = int(parts[1])
            end       = int(parts[2])
            gene_id   = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
            gene_name = parts[4] if len(parts) > 4 else gene_id
            strand    = parts[5] if len(parts) > 5 else "."
            results.append({
                "chr": chrom, "start": start, "end": end,
                "gene_id": gene_id, "gene_name": gene_name, "strand": strand
            })
    return pd.DataFrame(results)


def sanitize_name(name):
    """
    Clean a sample name derived from a BigWig filename:
      1. Replace hyphens with underscores (safe for pandas column names)
      2. Convert to lowercase so capitalisation differences in filenames
         (e.g. H3K27Me3 vs H3K27me3) do not create spurious duplicate groups.
    e.g. 'B10-UNC-H3K27Me3-R1' -> 'b10_unc_h3k27me3_r1'
    """
    return name.replace("-", "_").lower()


def get_bw_name(bw_file):
    """
    Derive a clean sample name from a BigWig file path.
    Takes the filename prefix (everything before the first dot),
    then sanitizes it (hyphens -> underscores, lowercase).
    e.g. 'B10-H3K27Ac-R1.cpm.norm.bw' -> 'b10_h3k27ac_r1'
    """
    basename = os.path.basename(bw_file)   # strip directory path
    raw_name = basename.split('.')[0]      # take prefix before first dot
    return sanitize_name(raw_name)         # clean the name


def extract_signals_from_bw(bw_file, regions_df, summaries=None):
    """
    Extract one or more summary statistics from a single BigWig file
    for every region in regions_df.

    Column names for individual replicates do NOT include the stat suffix —
    only the sample name is used.
    e.g. sample 'b10_h3k27ac_r1' -> column 'b10_h3k27ac_r1'
    The stat suffix is only used internally for the replicate-mean columns.

    Parameters:
        bw_file    : path to the .bw file
        regions_df : DataFrame with columns chr, start, end
        summaries  : list of stat strings e.g. ['mean', 'max']

    Returns:
        copy of regions_df with new signal columns appended,
        named simply <sample> (no stat suffix on individual replicates)
        e.g. 'b10_h3k27ac_r1'
    """
    if summaries is None:
        summaries = ["mean"]

    bw_name = get_bw_name(bw_file)
    bw      = pyBigWig.open(bw_file)

    # For individual replicate columns we only store the primary stat (first in list).
    # Multiple stats are still computed but stored under the same column with stat label
    # only for the mean-of-replicates columns, not the raw replicate columns.
    # Store one list per stat internally, but expose only the primary stat as the
    # clean replicate column name.
    primary_stat = summaries[0]   # e.g. 'mean' — used for the clean replicate column
    signal_list  = []             # values for the primary stat across all regions

    for _, row in regions_df.iterrows():
        chrom = row["chr"]
        start = int(row["start"])
        end   = int(row["end"])

        if chrom not in bw.chroms():
            signal_list.append(np.nan)
            continue

        chrom_len = bw.chroms(chrom)
        start = max(0, start)
        end   = min(end, chrom_len)

        try:
            val = bw.stats(chrom, start, end, type=primary_stat)[0]
            signal_list.append(val if val is not None else np.nan)
        except Exception:
            signal_list.append(np.nan)

    bw.close()

    df_result = regions_df.copy()
    # Column name is just the sample name — no stat suffix on individual replicates
    df_result[bw_name] = signal_list

    return df_result


def get_replicate_base(sample_name):
    """
    Strip the replicate suffix (_r1, _r2, _r3, ...) from a sample name
    to obtain the condition+mark group name.
    Works on lowercased names so _R1 and _r1 both match.
    e.g. 'b10_h3k27ac_r1' -> 'b10_h3k27ac'
         'parental_h3k9me3_r2' -> 'parental_h3k9me3'
    """
    return re.sub(r'_r\d+$', '', sample_name)


def add_replicate_means(df, replicate_cols):
    """
    For each unique condition+mark group (e.g. b10_h3k27ac, parental_h3k27ac),
    compute the row-wise mean across all replicate columns of that group
    and insert it as a new column immediately after the last replicate column.

    New column name: <group>_mean
    e.g. mean of [b10_h3k27ac_r1, b10_h3k27ac_r2] -> b10_h3k27ac_mean

    Parameters:
        df             : DataFrame with individual replicate signal columns
        replicate_cols : list of column names that are individual replicates
                         (all columns added by extract_signals_from_bw)

    Returns:
        df with replicate-mean columns inserted after each group's replicates
        Also returns list of (group, mean_col_name) tuples for downstream use
    """
    # Group replicate columns by their base name (strip _r<n>)
    groups = {}
    for col in replicate_cols:
        group = get_replicate_base(col)   # e.g. 'b10_h3k27ac'
        groups.setdefault(group, []).append(col)

    group_mean_cols = []   # track (group, mean_col) for difference computation

    for group, cols in groups.items():
        mean_col  = f"{group}_mean"              # e.g. 'b10_h3k27ac_mean'
        mean_vals = df[cols].mean(axis=1)        # row-wise mean across replicates

        # Insert immediately after the last replicate column of this group
        last_idx  = max(df.columns.get_loc(c) for c in cols)
        df.insert(last_idx + 1, mean_col, mean_vals)
        group_mean_cols.append((group, mean_col))
        print(f"  [MEAN] {mean_col} = mean of {cols}")

    return df, group_mean_cols


def add_automatic_differences(df, group_mean_cols, cond_a="b10", cond_b="parental"):
    """
    Automatically find all mark groups that have a replicate-mean column for
    both cond_a and cond_b, then compute cond_a - cond_b for each matched mark.

    Matching logic:
      - For each group mean column starting with cond_a (e.g. 'b10_h3k27ac_mean'),
        extract the mark name by stripping the condition prefix
        e.g. 'b10_h3k27ac_mean' -> mark = 'h3k27ac'
      - Find the corresponding cond_b column with the same mark
        e.g. 'parental_h3k27ac_mean'
      - Compute difference and name it: <cond_a>_minus_<cond_b>_<mark>
        e.g. 'b10_minus_parental_h3k27ac'

    Parameters:
        df              : DataFrame with replicate-mean columns
        group_mean_cols : list of (group, mean_col) tuples from add_replicate_means
        cond_a          : first condition (minuend)     default 'b10'
        cond_b          : second condition (subtrahend) default 'parental'

    Returns:
        df with difference columns appended
        list of difference column names created
    """
    # Build lookup: mark -> mean_col for each condition
    # Strip condition prefix and trailing _mean to isolate the mark
    map_a, map_b = {}, {}
    for group, mean_col in group_mean_cols:
        if group.startswith(f"{cond_a}_"):
            mark = group[len(f"{cond_a}_"):]        # e.g. 'h3k27ac'
            map_a[mark] = mean_col
        elif group.startswith(f"{cond_b}_"):
            mark = group[len(f"{cond_b}_"):]        # e.g. 'h3k27ac'
            map_b[mark] = mean_col

    common_marks  = sorted(set(map_a) & set(map_b))
    diff_cols_created = []

    if not common_marks:
        print(f"[WARN] No matching marks found between {cond_a} and {cond_b}.")
        return df, diff_cols_created

    for mark in common_marks:
        diff_col = f"{cond_a}_minus_{cond_b}_{mark}"   # e.g. b10_minus_parental_h3k27ac
        df[diff_col] = df[map_a[mark]] - df[map_b[mark]]
        diff_cols_created.append(diff_col)
        print(f"  [DIFF] {diff_col} = {map_a[mark]} - {map_b[mark]}")

    return df, diff_cols_created


def save_outputs(df, diff_cols, base_cols, output_dir, output_file):
    """
    Save two output files:

    File 1 — full table (<output_file>):
        All columns: base region info + all replicate signals +
        replicate means + difference columns.

    File 2 — difference-only table (<stem>_diff_only.txt):
        Base region info columns only + difference columns.
        Compact table for direct downstream integration with RNA-seq / Hi-C.

    Parameters:
        df          : complete DataFrame
        diff_cols   : list of difference column names
        base_cols   : list of base region column names (chr, start, end, etc.)
        output_dir  : directory to write files into
        output_file : filename for the full table
    """
    os.makedirs(output_dir, exist_ok=True)

    # ── File 1: full table ────────────────────────────────────────────────────
    full_path = os.path.join(output_dir, output_file)
    df.to_csv(full_path, sep="\t", index=False)
    print(f"[INFO] Full table saved      -> {full_path}")
    print(f"       {df.shape[0]} rows x {df.shape[1]} columns")

    # ── File 2: base + difference columns only ────────────────────────────────
    stem      = output_file.rsplit(".", 1)[0] if "." in output_file else output_file
    diff_file = f"{stem}_diff_only.txt"
    diff_path = os.path.join(output_dir, diff_file)

    existing_base = [c for c in base_cols if c in df.columns]
    diff_subset   = df[existing_base + diff_cols]
    diff_subset.to_csv(diff_path, sep="\t", index=False)
    print(f"[INFO] Diff-only table saved -> {diff_path}")
    print(f"       {diff_subset.shape[0]} rows x {diff_subset.shape[1]} columns")
    print(f"       Columns: {list(diff_subset.columns)}")


def main(bed_file, output_dir, bw_files, summaries=None,
         output_file="cut_run_signal_table.txt",
         cond_a="b10", cond_b="parental"):
    """
    Full pipeline:
      1. Read BED regions
      2. Extract BigWig signal for each file (clean column names, no stat suffix)
      3. Compute replicate means per condition+mark group
      4. Automatically compute condition differences (cond_a - cond_b)
      5. Save full table and difference-only table
    """
    if summaries is None:
        summaries = ["mean"]
    if not bw_files:
        raise ValueError("No BigWig files provided")

    # Lowercase condition names to match sanitized column names
    cond_a = cond_a.lower()
    cond_b = cond_b.lower()

    # Base region columns included in both output files
    base_cols = ["chr", "start", "end", "gene_id", "gene_name", "strand"]

    # ── Step 1: Load genomic regions ─────────────────────────────────────────
    df_regions = extract_region_info(bed_file)
    print(f"[INFO] Loaded {len(df_regions)} regions from {bed_file}")

    # ── Step 2: Extract signal from all BigWig files ─────────────────────────
    replicate_cols = []   # track the clean replicate column names in order

    df_result = extract_signals_from_bw(bw_files[0], df_regions, summaries=summaries)
    bw_name_0 = get_bw_name(bw_files[0])
    replicate_cols.append(bw_name_0)
    print(f"[INFO] Extracted: {os.path.basename(bw_files[0])} -> column '{bw_name_0}'")

    for bw_file in bw_files[1:]:
        df_bw   = extract_signals_from_bw(bw_file, df_regions, summaries=summaries)
        bw_name = get_bw_name(bw_file)

        # Skip duplicates — if two files map to the same sanitized name, warn and skip
        if bw_name in df_result.columns:
            print(f"[WARN] Duplicate sanitized name '{bw_name}' for "
                  f"'{os.path.basename(bw_file)}' — skipping to avoid duplication.")
            continue

        df_result[bw_name] = df_bw[bw_name]
        replicate_cols.append(bw_name)
        print(f"[INFO] Extracted: {os.path.basename(bw_file)} -> column '{bw_name}'")

    # ── Step 3: Compute replicate means ──────────────────────────────────────
    print("[INFO] Computing replicate means...")
    df_result, group_mean_cols = add_replicate_means(df_result, replicate_cols)

    # ── Step 4: Compute automatic condition differences ───────────────────────
    print(f"[INFO] Computing differences: {cond_a} - {cond_b}...")
    df_result, diff_cols = add_automatic_differences(
        df_result, group_mean_cols, cond_a=cond_a, cond_b=cond_b
    )

    # ── Step 5: Save outputs ──────────────────────────────────────────────────
    save_outputs(df_result, diff_cols, base_cols, output_dir, output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Extract CUT&RUN BigWig signals, compute replicate means, "
            "automatically compute condition differences, "
            "and save a full table plus a difference-only table."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python extract_bigwig_signal.py regions.bed \\
      --bw raw/data/* --output_dir results/

  python extract_bigwig_signal.py regions.bed \\
      --bw raw/data/* --output_dir results/ \\
      --cond_a B10 --cond_b Parental
        """
    )
    parser.add_argument("bed_file",      help="BED file with genomic regions")
    parser.add_argument("--bw",          nargs="+", required=True, dest="bw_files",
                        help="One or more BigWig files")
    parser.add_argument("--summary",     nargs="+", default=["mean"],
                        choices=["mean", "max", "min", "sum"],
                        help="Summary statistic (default: mean)")
    parser.add_argument("--output_dir",  default=".", help="Output directory")
    parser.add_argument("--output_file", default="cut_run_signal_table.txt",
                        help="Output filename for full table")
    parser.add_argument("--cond_a",      default="B10",
                        help="First condition — minuend (default: B10)")
    parser.add_argument("--cond_b",      default="Parental",
                        help="Second condition — subtrahend (default: Parental)")

    args = parser.parse_args()

    main(args.bed_file, args.output_dir, args.bw_files,
         summaries=args.summary,
         output_file=args.output_file,
         cond_a=args.cond_a,
         cond_b=args.cond_b)