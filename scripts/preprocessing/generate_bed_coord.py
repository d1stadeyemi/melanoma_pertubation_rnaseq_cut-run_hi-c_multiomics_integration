import pandas as pd
import argparse
import os
import re

STEP = 5000
COORD_CANDIDATES = ["coordinate", "coordinates", "coord", "coords", "coor", "coors"]


class GridConvention:
    """Encapsulates the bin-coordinate arithmetic that differs between
    1-based (inclusive) and 0-based (half-open) systems.

    For a grid of size `step`:
      1-based: valid starts {1, step+1, ...}, ends {step, 2*step, ...},
               bins emitted inclusive [start, end].
      0-based: valid starts {0, step, ...},   ends {step, 2*step, ...},
               bins emitted half-open [start, end).
    """
    def __init__(self, step=STEP, zero_based=False):
        self.step = step
        self.zero_based = zero_based
        self.offset = 0 if zero_based else 1  # subtracted from start before grid math

    @property
    def name(self):
        return "0-based half-open" if self.zero_based else "1-based inclusive"

    def start_on_grid(self, span_start):
        return (span_start - self.offset) % self.step == 0

    def end_on_grid(self, span_end):
        # End boundary is the same physical position in both conventions.
        return span_end % self.step == 0

    def snap_start(self, span_start):
        return ((span_start - self.offset) // self.step) * self.step + self.offset

    def example_starts(self):
        b = self.offset
        return f"{b}, {b + self.step}, {b + 2 * self.step}, ..."


def find_coord_col(df, user_col=None):
    if user_col:
        if user_col not in df.columns:
            raise ValueError(f"Column '{user_col}' not found. Available: {list(df.columns)}")
        return user_col
    lower_map = {c.lower(): c for c in df.columns}
    for cand in COORD_CANDIDATES:
        if cand in lower_map:
            return lower_map[cand]
    raise ValueError(
        f"No coordinate column auto-detected (looked for {COORD_CANDIDATES}). "
        f"Re-run with --coord-col <name>. Available columns: {list(df.columns)}"
    )


def parse_coord(coord):
    # Handles 'chr1:1-250000 - chr1:135250001-135500000' and single 'chr1:1-250000'
    parts = [p.strip() for p in coord.split(" - ")]
    first, last = parts[0], parts[-1]
    m_first = re.match(r"^(chr[\w]+):(\d+)-(\d+)$", first)
    m_last = re.match(r"^(chr[\w]+):(\d+)-(\d+)$", last)
    if not m_first or not m_last:
        raise ValueError(f"Could not parse coordinate: '{coord}'")
    chrom = m_first.group(1)
    span_start = int(m_first.group(2))
    span_end = int(m_last.group(3))
    return chrom, span_start, span_end


def check_consistency(records, conv):
    """records: list of (coord, chrom, span_start, span_end).
    Returns (start_offenders, end_offenders)."""
    start_offenders, end_offenders = [], []
    for coord, chrom, span_start, span_end in records:
        if not conv.start_on_grid(span_start):
            start_offenders.append((coord, span_start, conv.snap_start(span_start)))
        if not conv.end_on_grid(span_end):
            clean_end = (span_end // conv.step) * conv.step
            end_offenders.append((coord, span_end, clean_end))
    return start_offenders, end_offenders


def tile_coord(chrom, span_start, span_end, coord, conv, snap_grid=False):
    """Tile [span_start, span_end] into bins under the given convention."""
    step = conv.step
    start = conv.snap_start(span_start) if snap_grid else span_start
    rows = []
    while start <= span_end:
        if conv.zero_based:
            end = min(start + step, span_end)        # half-open [start, end)
        else:
            end = min(start + step - 1, span_end)    # inclusive [start, end]
        rows.append((chrom, start, end, coord))
        start += step
    return rows


def load_table(input_file):
    ext = os.path.splitext(input_file)[1].lower()
    if ext in (".xlsx", ".xls"):
        return pd.read_excel(input_file)
    if ext == ".csv":
        return pd.read_csv(input_file)
    return pd.read_csv(input_file, sep="\t")


def _format_offenders(offenders, label, fix_label, flag_hint, conv):
    expected = conv.example_starts() if label == "start" else f"{conv.step}, {2 * conv.step}, ..."
    lines = [
        f"ERROR: Inconsistent coordinate {label} boundaries detected "
        f"(convention: {conv.name}, step={conv.step}).",
        f"Expected {label}s on the grid: {expected}",
        f"Found {len(offenders)} off-grid {label}(s). First few:",
    ]
    for coord, val, fixed in offenders[:10]:
        lines.append(f"  - '{coord}': {label}={val} ({fix_label} {fixed})")
    if len(offenders) > 10:
        lines.append(f"  ... and {len(offenders) - 10} more.")
    lines.append("")
    lines.append(flag_hint)
    return lines


def generate_bed(input_file, coord_col=None, step=STEP, zero_based=False,
                 snap_grid=False, allow_partial_end=False):
    conv = GridConvention(step=step, zero_based=zero_based)
    df = load_table(input_file)
    col = find_coord_col(df, coord_col)

    records = []
    for coord in df[col].dropna().astype(str):
        chrom, span_start, span_end = parse_coord(coord)
        records.append((coord, chrom, span_start, span_end))

    start_offenders, end_offenders = check_consistency(records, conv)

    error_blocks = []
    if start_offenders and not snap_grid:
        error_blocks.append(_format_offenders(
            start_offenders, "start", "would snap to",
            "Re-run with --snap-grid to align every bin start to the grid and continue.",
            conv,
        ))
    if end_offenders and not allow_partial_end:
        error_blocks.append(_format_offenders(
            end_offenders, "end", "last clean grid end is",
            "Re-run with --allow-partial-end to permit a short final bin (clipped to the coordinate end) and continue.",
            conv,
        ))
    if error_blocks:
        raise SystemExit("\n\n".join("\n".join(block) for block in error_blocks))

    all_rows = []
    for coord, chrom, span_start, span_end in records:
        all_rows.extend(
            tile_coord(chrom, span_start, span_end, coord, conv, snap_grid=snap_grid)
        )

    bed = pd.DataFrame(all_rows, columns=["chr", "start", "end", "coord"])
    bed["start"] = bed["start"].astype(int)
    bed["end"] = bed["end"].astype(int)
    bed.drop_duplicates(subset=["chr", "start", "end"], inplace=True)
    bed.sort_values(["chr", "start"], inplace=True)
    bed.reset_index(drop=True, inplace=True)

    # Positional unique ID (delimiter-safe), placed in BED column 4 (name field).
    bed["bin_id"] = (
        bed["chr"] + "_" + bed["start"].astype(str) + "_" + bed["end"].astype(str)
    )
    bed = bed[["chr", "start", "end", "bin_id", "coord"]]
    return bed


def save_bed(bed, output_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    bed.to_csv(os.path.join(output_dir, output_file), sep="\t", header=False, index=False)


def save_txt(bed, output_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    bed.to_csv(os.path.join(output_dir, output_file), sep="\t", header=True, index=False)


def main(input_file, output_dir, coord_col=None, step=STEP, zero_based=False,
         snap_grid=False, allow_partial_end=False):
    bed = generate_bed(input_file, coord_col=coord_col, step=step, zero_based=zero_based,
                       snap_grid=snap_grid, allow_partial_end=allow_partial_end)
    save_bed(bed, "coordinates_5kb_bins.bed", output_dir)
    save_txt(bed, "coordinates_5kb_bins.txt", output_dir)
    print(f"Files saved to {output_dir}")
    print(f"Convention: {conv_name(zero_based)} | step: {step}")
    print(f"Columns: {list(bed.columns)}")
    print(f"Rows: {len(bed)}")
    if snap_grid:
        print("Note: --snap-grid was active; bin starts were aligned to the grid.")
    if allow_partial_end:
        print("Note: --allow-partial-end was active; some final bins may be shorter than the step size.")


def conv_name(zero_based):
    return "0-based half-open" if zero_based else "1-based inclusive"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Tile genomic coordinates into fixed-size-bin BED + txt files"
    )
    parser.add_argument("input_file", help="Path to input file (xlsx/csv/tsv/txt) with a coordinate column")
    parser.add_argument("--coord-col", type=str, default=None,
                        help="Name of the coordinate column (if not auto-detected)")
    parser.add_argument("--step", type=int, default=STEP, help=f"Bin size in bp (default: {STEP})")
    parser.add_argument("--output_dir", type=str, default=".", help="Output directory (default: .)")
    parser.add_argument("--zero-based", dest="zero_based", action="store_true",
                        help="Treat coordinates as 0-based and emit true half-open [start, end) BED "
                             "(default: 1-based inclusive)")
    parser.add_argument("--snap-grid", dest="snap_grid", action="store_true",
                        help="Align bin starts to the grid instead of erroring on off-grid starts")
    parser.add_argument("--allow-partial-end", dest="allow_partial_end", action="store_true",
                        help="Permit a short final bin instead of erroring on off-grid ends")

    args = parser.parse_args()
    main(args.input_file, args.output_dir, coord_col=args.coord_col, step=args.step,
         zero_based=args.zero_based, snap_grid=args.snap_grid, allow_partial_end=args.allow_partial_end)