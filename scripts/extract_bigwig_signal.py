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
    Replace all hyphens with underscores so column names are safe to use
    in pandas and downstream tools.
    e.g. 'B10-H3K27Ac-R1' -> 'B10_H3K27Ac_R1'
    """
    return name.replace("-", "_")


def get_bw_name(bw_file):
    """
    Derive a clean sample name from a BigWig file path.
    Takes the filename prefix (everything before the first dot),
    then replaces hyphens with underscores.
    e.g. 'B10-H3K27Ac-R1.cpm.norm.bw' -> 'B10_H3K27Ac_R1'
    """
    basename = os.path.basename(bw_file)   # strip directory path
    raw_name = basename.split('.')[0]      # take prefix before first dot
    return sanitize_name(raw_name)         # replace hyphens with underscores


def extract_signals_from_bw(bw_file, regions_df, summaries=None):
    """
    Extract one or more summary statistics from a single BigWig file
    for every region in regions_df.

    For each region, queries pyBigWig for each requested stat type.
    Regions on chromosomes absent from the BigWig get NaN.
    Coordinates are clamped to chromosome boundaries before querying.

    Parameters:
        bw_file    : path to the .bw file
        regions_df : DataFrame with columns chr, start, end
        summaries  : list of stat strings e.g. ['mean', 'max']

    Returns:
        copy of regions_df with new columns appended, named <sample>_<stat>
        e.g. 'B10_H3K27Ac_R1_mean'
    """
    if summaries is None:
        summaries = ["mean"]

    bw_name = get_bw_name(bw_file)         # clean sample name
    bw      = pyBigWig.open(bw_file)

    # Pre-allocate one empty list per requested stat
    signal_data = {f"{bw_name}_{stat}": [] for stat in summaries}

    for _, row in regions_df.iterrows():
        chrom = row["chr"]
        start = int(row["start"])
        end   = int(row["end"])

        # If chromosome not in this BigWig, record NaN for all stats
        if chrom not in bw.chroms():
            for stat in summaries:
                signal_data[f"{bw_name}_{stat}"].append(np.nan)
            continue

        # Clamp start/end to valid chromosome range
        chrom_len = bw.chroms(chrom)
        start = max(0, start)
        end   = min(end, chrom_len)

        for stat in summaries:
            try:
                val = bw.stats(chrom, start, end, type=stat)[0]
                signal_data[f"{bw_name}_{stat}"].append(
                    val if val is not None else np.nan
                )
            except Exception:
                signal_data[f"{bw_name}_{stat}"].append(np.nan)

    bw.close()

    # Attach all signal columns to a copy of the regions DataFrame
    df_result = regions_df.copy()
    for stat in summaries:
        col = f"{bw_name}_{stat}"
        df_result[col] = signal_data[col]

    return df_result


def get_replicate_base(sample_name):
    """
    Strip the replicate suffix (R1, R2, R3, ...) from a sample name
    to obtain the condition+mark group name.
    Uses a regex that removes a trailing _R followed by one or more digits.
    e.g. 'B10_H3K27Ac_R1' -> 'B10_H3K27Ac'
         'Parental_H3K9me3_R2' -> 'Parental_H3K9me3'
    """
    return re.sub(r'_R\d+$', '', sample_name)


def add_replicate_means(df, summaries):
    """
    For each unique condition+mark group (e.g. B10_H3K27Ac, Parental_H3K27Ac),
    compute the row-wise mean across all replicate columns of that group
    and insert it as a new column immediately after the last replicate column.

    New column is named: <group>_<stat>_mean
    e.g. mean of [B10_H3K27Ac_R1_mean, B10_H3K27Ac_R2_mean] -> B10_H3K27Ac_mean

    Parameters:
        df        : DataFrame containing individual replicate signal columns
        summaries : list of stat types that were extracted e.g. ['mean']

    Returns:
        df with replicate-mean columns inserted after each group's replicates
    """
    for stat in summaries:
        # Find all individual replicate signal columns for this stat type
        stat_cols = [c for c in df.columns if c.endswith(f"_{stat}")]

        # Group columns by their condition+mark base name
        # e.g. 'B10_H3K27Ac_R1_mean' -> sample='B10_H3K27Ac_R1' -> group='B10_H3K27Ac'
        groups = {}
        for col in stat_cols:
            sample = col[: col.rfind(f"_{stat}")]   # strip trailing _<stat>
            group  = get_replicate_base(sample)      # strip trailing _R<n>
            groups.setdefault(group, []).append(col)

        # Compute and insert mean column for each group
        for group, cols in groups.items():
            mean_col  = f"{group}_{stat}_mean"           # new column name
            mean_vals = df[cols].mean(axis=1)            # row-wise mean across replicates

            # Insert immediately after the last replicate column of this group
            last_idx  = max(df.columns.get_loc(c) for c in cols)
            df.insert(last_idx + 1, mean_col, mean_vals)

    return df


def add_automatic_differences(df, summaries, cond_a="B10", cond_b="Parental"):
    """
    Automatically find all mark groups that have a replicate-mean column for
    both cond_a (e.g. B10) and cond_b (e.g. Parental), then compute
    cond_a - cond_b for each matched mark and append as a new column.

    How matching works:
      - Scans all columns ending in '_<stat>_mean'
      - Strips the condition prefix and stat suffix to isolate the mark name
        e.g. 'B10_H3K27Ac_mean_mean' -> mark = 'H3K27Ac'
      - Pairs columns sharing the same mark but differing by condition prefix
      - Difference column: <cond_a>_minus_<cond_b>_<mark>
        e.g. B10_minus_Parental_H3K27Ac

    Parameters:
        df        : DataFrame already containing replicate-mean columns
        summaries : list of stat types e.g. ['mean']
        cond_a    : name of the first condition (minuend)     default 'B10'
        cond_b    : name of the second condition (subtrahend) default 'Parental'

    Returns:
        df with difference columns appended at the end
        Also returns a list of the difference column names created,
        so the caller can use them to build the summary-only output file.
    """
    diff_cols_created = []   # track names of all difference columns created

    for stat in summaries:

        def build_mark_map(cond):
            """
            For a given condition prefix, scan all replicate-mean columns
            and return a dict mapping mark name -> column name.
            Strips '<cond>_' from the front and '_<stat>_mean' from the end.
            e.g. 'B10_H3K27Ac_mean_mean' -> {'H3K27Ac': 'B10_H3K27Ac_mean_mean'}
            """
            mark_map = {}
            suffix = f"_{stat}_mean"
            for col in df.columns:
                if col.startswith(f"{cond}_") and col.endswith(suffix):
                    # Remove condition prefix and stat suffix
                    inner = col[len(f"{cond}_") : col.rfind(suffix)]
                    # Strip any trailing _mean that crept in from the stat name
                    mark = re.sub(r'_mean$', '', inner)
                    mark_map[mark] = col
            return mark_map

        map_a = build_mark_map(cond_a)   # e.g. {'H3K27Ac': 'B10_H3K27Ac_mean_mean'}
        map_b = build_mark_map(cond_b)   # e.g. {'H3K27Ac': 'Parental_H3K27Ac_mean_mean'}

        # Compute difference only for marks present in both conditions
        common_marks = sorted(set(map_a) & set(map_b))

        if not common_marks:
            print(f"[WARN] No matching marks found between {cond_a} and {cond_b}.")
            continue

        for mark in common_marks:
            diff_col = f"{cond_a}_minus_{cond_b}_{mark}"   # e.g. B10_minus_Parental_H3K27Ac
            df[diff_col] = df[map_a[mark]] - df[map_b[mark]]
            diff_cols_created.append(diff_col)              # record for output filtering
            print(f"  [DIFF] {diff_col} = {map_a[mark]} - {map_b[mark]}")

    return df, diff_cols_created


def save_outputs(df, diff_cols, base_cols, output_dir, output_file):
    """
    Save two output files:

    File 1 — full table (<output_file>):
        All columns: base region info + all signal columns + all difference columns.

    File 2 — difference-only table (<output_file stem>_diff_only.txt):
        Base region info columns (chr, start, end, gene_id, gene_name, strand)
        plus only the difference columns (e.g. B10_minus_Parental_H3K27Ac).
        This is the compact table most useful for downstream integration.

    Parameters:
        df          : complete DataFrame with all columns
        diff_cols   : list of difference column names created by add_automatic_differences()
        base_cols   : list of base region column names to include in both files
        output_dir  : directory to write files into
        output_file : filename for the full table (stem used for diff-only name)
    """
    os.makedirs(output_dir, exist_ok=True)

    # ── File 1: full table ────────────────────────────────────────────────────
    full_path = os.path.join(output_dir, output_file)
    df.to_csv(full_path, sep="\t", index=False)
    print(f"[INFO] Full table saved    -> {full_path}")
    print(f"       Shape: {df.shape[0]} rows x {df.shape[1]} columns")

    # ── File 2: base columns + difference columns only ────────────────────────
    # Build the stem of the output filename (strip extension) for the diff file
    stem      = output_file.rsplit(".", 1)[0] if "." in output_file else output_file
    diff_file = f"{stem}_diff_only.txt"
    diff_path = os.path.join(output_dir, diff_file)

    # Select base region columns that actually exist in df, then add diff cols
    existing_base = [c for c in base_cols if c in df.columns]
    diff_subset   = df[existing_base + diff_cols]
    diff_subset.to_csv(diff_path, sep="\t", index=False)
    print(f"[INFO] Diff-only table saved -> {diff_path}")
    print(f"       Shape: {diff_subset.shape[0]} rows x {diff_subset.shape[1]} columns")
    print(f"       Columns: {list(diff_subset.columns)}")


def main(bed_file, output_dir, bw_files, summaries=None,
         output_file="cut_run_signal_table.txt",
         cond_a="B10", cond_b="Parental"):
    """
    Full pipeline:
      1. Read BED regions
      2. Extract BigWig signal for each file
      3. Compute replicate means per condition+mark group
      4. Automatically compute condition differences (cond_a - cond_b)
      5. Save full table and difference-only table
    """
    if summaries is None:
        summaries = ["mean"]
    if not bw_files:
        raise ValueError("No BigWig files provided")

    # Base region columns — always included in both output files
    base_cols = ["chr", "start", "end", "gene_id", "gene_name", "strand"]

    # ── Step 1: Load genomic regions from BED ────────────────────────────────
    df_regions = extract_region_info(bed_file)
    print(f"[INFO] Loaded {len(df_regions)} regions from {bed_file}")

    # ── Step 2: Extract signal from all BigWig files ─────────────────────────
    # Start with the first file, then attach new signal columns from the rest
    df_result = extract_signals_from_bw(bw_files[0], df_regions, summaries=summaries)
    print(f"[INFO] Extracted: {os.path.basename(bw_files[0])}")

    for bw_file in bw_files[1:]:
        df_bw   = extract_signals_from_bw(bw_file, df_regions, summaries=summaries)
        bw_name = get_bw_name(bw_file)
        # Attach only the new signal columns, not the duplicate region columns
        new_cols = [c for c in df_bw.columns if c.startswith(f"{bw_name}_")]
        for col in new_cols:
            df_result[col] = df_bw[col]
        print(f"[INFO] Extracted: {os.path.basename(bw_file)}")

    # ── Step 3: Compute replicate means per condition+mark group ─────────────
    print("[INFO] Computing replicate means...")
    df_result = add_replicate_means(df_result, summaries)

    # ── Step 4: Automatically compute condition differences ──────────────────
    print(f"[INFO] Computing automatic differences: {cond_a} - {cond_b}...")
    df_result, diff_cols = add_automatic_differences(
        df_result, summaries, cond_a=cond_a, cond_b=cond_b
    )

    # ── Step 5: Save both output files ───────────────────────────────────────
    save_outputs(df_result, diff_cols, base_cols, output_dir, output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Extract CUT&RUN BigWig signals, compute replicate means, "
            "automatically subtract Parental from B10 per mark, "
            "and save a full table plus a difference-only table."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Two marks, two replicates each, auto-difference B10 vs Parental:
  python extract_bigwig_signal.py regions.bed \\
      --bw B10-H3K27Ac-R1.cpm.norm.bw  B10-H3K27Ac-R2.cpm.norm.bw \\
           Parental-H3K27Ac-R1.cpm.norm.bw Parental-H3K27Ac-R2.cpm.norm.bw \\
           B10-H3K4me1-R1.cpm.norm.bw   B10-H3K4me1-R2.cpm.norm.bw \\
           Parental-H3K4me1-R1.cpm.norm.bw Parental-H3K4me1-R2.cpm.norm.bw \\
      --output_dir results/

  # Custom condition names:
  python extract_bigwig_signal.py regions.bed \\
      --bw Treated-H3K27me3-R1.cpm.norm.bw Control-H3K27me3-R1.cpm.norm.bw \\
      --output_dir results/ --cond_a Treated --cond_b Control
        """
    )
    parser.add_argument("bed_file",
                        help="BED file with genomic regions of interest")
    parser.add_argument("--bw", nargs="+", required=True, dest="bw_files",
                        help="One or more BigWig files")
    parser.add_argument("--summary", nargs="+", default=["mean"],
                        choices=["mean", "max", "min", "sum"],
                        help="Summary statistic(s) to extract (default: mean)")
    parser.add_argument("--output_dir",  default=".",
                        help="Output directory (default: current directory)")
    parser.add_argument("--output_file", default="cut_run_signal_table.txt",
                        help="Output filename for full table (default: cut_run_signal_table.txt)")
    parser.add_argument("--cond_a", default="B10",
                        help="First condition — minuend (default: B10)")
    parser.add_argument("--cond_b", default="Parental",
                        help="Second condition — subtrahend (default: Parental)")

    args = parser.parse_args()

    main(args.bed_file, args.output_dir, args.bw_files,
         summaries=args.summary,
         output_file=args.output_file,
         cond_a=args.cond_a,
         cond_b=args.cond_b)