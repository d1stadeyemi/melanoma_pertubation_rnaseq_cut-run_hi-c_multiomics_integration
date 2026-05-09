import pandas as pd
import numpy as np
import os 
import argparse

def merge_all(rnaseq, hic, cut_run):
    rnaseq = pd.read_csv(rnaseq, sep="\t")
    hic = pd.read_csv(hic, sep="\t")
    cut_run = pd.read_csv(cut_run, sep="\t")

    cut_run["Bottom-Parent"] = cut_run["Bottom-H3K27me3-R1_mean"] - cut_run["Control-H3K27me3-R1_mean"]
    hic_clean = hic.drop(columns=["chr", "start", "end"])

    merged_df = (
        hic_clean
        .merge(rnaseq[["gene_name", "padj", "pvalue", "log2FC", "baseMean", "description"]], on="gene_name", how="outer") 
        .merge(cut_run[["gene_name", "Control-H3K27me3-R1_mean", "Bottom-H3K27me3-R1_mean", "Bottom-Parent"]], 
                on="gene_name", how="outer"))
    
    cols_to_move = ["gene_name", "coordinate", "padj", "pvalue", "log2FC", "baseMean", 
                "Control-H3K27me3-R1_mean", "Bottom-H3K27me3-R1_mean", "Bottom-Parent"]

    cols = cols_to_move + [c for c in merged_df.columns if c not in cols_to_move]
    merged_arr = merged_df[cols]
  
    return merged_arr

def save_txt(merged_arr, output_file, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    merged_arr.to_csv(os.path.join(output_dir, output_file), sep="\t", index=False)

def main(rnaseq, hic, cut_run, output_dir):
    merged_arr = merge_all(rnaseq, hic, cut_run)
    save_txt(merged_arr, "merged_hic_rnaseq_cutrun_parent_bottom10.txt", output_dir)

    print(f"Output txt file saved to {output_dir} as merged_hic_rnaseq_cutrun_parent_bottom10.txt")
    print(f"Rows: {len(merged_arr)}")
    print(f"Columns: {merged_arr.shape[1]}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Merge DESeq, Hi-C PC1 and cut and run extracted signal files")
    parser.add_argument("rnaseq", help="Path to DESeq output file in txt")
    parser.add_argument("hic", help="Path to Hi-C PC1 txt file")
    parser.add_argument("cut_run", help="cut and run extracted signal file in txt")
    parser.add_argument("--output_dir", type=str, default=".",
                        help="Output directory for the txt file (default: current directory)")
    
    args = parser.parse_args()

    main(args.rnaseq, args.hic, args.cut_run, args.output_dir)