import pandas as pd
import numpy as np
import pyBigWig
import argparse
import os

def extract_region_info(bed_file):
    """Extract base region information from BED file"""
    results = []
    
    with open(bed_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            gene_id = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
            gene_name = parts[4] if len(parts) > 4 else gene_id
            strand = parts[5] if len(parts) > 5 else "."
            
            results.append({
                "chr": chrom,
                "start": start,
                "end": end,
                "gene_id": gene_id,
                "gene_name": gene_name,
                "strand": strand
            })
    
    return pd.DataFrame(results)

def extract_signals_from_bw(bw_file, regions_df, summaries=None, bw_name=None):
    """Extract signals from a single BigWig file for multiple summary statistics"""
    if summaries is None:
        summaries = ["mean"]
    
    if bw_name is None:
        # Extract name up to first dot (e.g., "H3k27me-R1" from "H3k27me-R1.cpm.norm.bw")
        basename = os.path.basename(bw_file)
        bw_name = basename.split('.')[0]
    
    bw = pyBigWig.open(bw_file)
    
    # Create a dictionary to store signal data for each summary type
    signal_data = {f"{bw_name}_{stat}": [] for stat in summaries}
    
    for idx, row in regions_df.iterrows():
        chrom, start, end = row["chr"], int(row["start"]), int(row["end"])
        
        if chrom not in bw.chroms():
            # Add NaN for this region across all summary types
            for stat in summaries:
                signal_data[f"{bw_name}_{stat}"].append(np.nan)
            continue
        
        chrom_len = bw.chroms(chrom)
        start = max(0, start)
        end = min(end, chrom_len)
        
        # Extract all requested summary statistics for this region
        for stat in summaries:
            try:
                signal = bw.stats(chrom, start, end, type=stat)[0]
                signal = signal if signal is not None else np.nan
            except Exception:
                signal = np.nan
            
            signal_data[f"{bw_name}_{stat}"].append(signal)
    
    bw.close()
    
    # Add signal columns to result dataframe
    df_result = regions_df.copy()
    for stat in summaries:
        df_result[f"{bw_name}_{stat}"] = signal_data[f"{bw_name}_{stat}"]
    
    return df_result

def main(bed_file, output_dir, bw_files, summaries=None, output_file="signal_table.txt"):
    """Extract signals from multiple BigWig files with multiple summary statistics"""
    if summaries is None:
        summaries = ["mean"]
    
    if not bw_files:
        raise ValueError("No BigWig files provided")
    
    # Read BED file once
    df_regions = extract_region_info(bed_file)
    
    # Start with the first BigWig file
    df_result = extract_signals_from_bw(bw_files[0], df_regions, summaries=summaries)
    
    # Merge signals from additional BigWig files
    for bw_file in bw_files[1:]:
        df_bw = extract_signals_from_bw(bw_file, df_regions, summaries=summaries)
        # Get only the signal columns (not the base region info)
        basename = os.path.basename(bw_file)
        bw_name = basename.split('.')[0]
        signal_cols = [col for col in df_bw.columns if col.startswith(f"{bw_name}_")]
        for col in signal_cols:
            df_result[col] = df_bw[col]
    
    # Save to output
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    output_path = os.path.join(output_dir, output_file)
    df_result.to_csv(output_path, sep="\t", header=True, index=False)
    print(f"Signal table saved to {output_path}")
    print(f"Columns: {list(df_result.columns)}")
    print(f"Rows: {len(df_result)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract signals from multiple BigWig files for regions defined in a BED file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single BigWig, mean signal only:
  python extract_bigwig_signal.py regions.bed --bw file1.bw --output_dir results/
  
  # Multiple BigWigs, mean signal:
  python extract_bigwig_signal.py regions.bed --bw file1.bw file2.bw file3.bw --output_dir results/
  
  # Multiple BigWigs, multiple summary statistics:
  python extract_bigwig_signal.py regions.bed --bw file1.bw file2.bw --summary mean max min --output_dir results/
        """
    )
    parser.add_argument("bed_file", help="Path to the BED file with regions of interest")
    parser.add_argument("--bw", nargs="+", required=True, dest="bw_files", help="Path(s) to BigWig file(s) (can specify multiple)")
    parser.add_argument(
        "--summary", 
        nargs="+", 
        default=["mean"], 
        choices=["mean", "max", "min", "sum"],
        help="Summary statistic(s) to extract (default: mean). Can specify multiple: --summary mean max sum"
    )
    parser.add_argument("--output_dir", type=str, default=".", help="Output directory for the signal table (default: current directory)")
    parser.add_argument("--output_file", type=str, default="cut_run_signal_table.txt", help="Output file name (default: signal_table.txt)")
    
    args = parser.parse_args()
    
    main(args.bed_file, args.output_dir, args.bw_files, summaries=args.summary, output_file=args.output_file)