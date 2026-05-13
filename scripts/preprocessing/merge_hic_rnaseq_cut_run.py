#!/lustre/isaac24/scratch/majagbe/envs/epp604_env/bin/python3

import pandas as pd
import numpy as np
import os
import argparse


def load_hic(hic_path):
    """
    Load the Hi-C compartment file and keep only the columns we need.
    Excluded: chr, start, end  (positional columns not needed after merging on gene_name).
    All remaining columns are kept: gene_name, coordinate, and all PC1/switch/shift columns.
    """
    hic  = pd.read_csv(hic_path, sep="\t")
    keep = [c for c in hic.columns if c not in ["chr", "start", "end"]]
    return hic[keep]


def load_rnaseq(rnaseq_path):
    """
    Load the DESeq2 output file and keep only the columns we need.
    Excluded: gene_id (redundant with gene_name), lfcSE, stat.
    Kept: gene_name, padj, pvalue, log2FC, baseMean, description.
    """
    rna  = pd.read_csv(rnaseq_path, sep="\t")
    keep = [c for c in rna.columns if c not in ["gene_id", "lfcSE", "stat"]]
    return rna[keep]


def load_cutrun(cutrun_path):
    """
    Load the CUT&RUN difference-only file and keep only the columns we need.
    Excluded: chr, start, end, gene_id  (positional/redundant columns).
    Kept: gene_name, strand, and all b10_minus_parental_* difference columns.
    """
    cr   = pd.read_csv(cutrun_path, sep="\t")
    keep = [c for c in cr.columns if c not in ["chr", "start", "end", "gene_id"]]
    return cr[keep]


def merge_all(hic, rna, cr):
    """
    Merge the three cleaned DataFrames on gene_name using outer joins
    so that all genes from all sources are retained even if absent in one source.

    Merge order:
      1. Hi-C (anchor — widest gene set, includes coordinate/compartment info)
      2. RNA-seq joined onto Hi-C
      3. CUT&RUN joined onto the result

    Final column order (explicit priority list):
      gene_name, coordinate
      -> RNA-seq stats: padj, pvalue, log2FC, baseMean
      -> strand
      -> all CUT&RUN b10_minus_parental_* difference columns (in their natural order)
      -> all Hi-C compartment columns (Parental_DMSO ... P->B10_Both_original)
      -> description (pushed to end)
      -> any remaining columns not explicitly listed
    """
    merged = (
        hic
        .merge(rna, on="gene_name", how="outer")   # attach RNA-seq columns
        .merge(cr,  on="gene_name", how="outer")   # attach CUT&RUN columns
    )

    # ── Dynamically collect column groups ─────────────────────────────────────

    # All CUT&RUN difference columns — detected by their naming pattern
    cutrun_diff_cols = [c for c in merged.columns
                        if c.startswith("b10_minus_parental_")]

    # Hi-C compartment value and switch/shift columns — everything from the
    # Hi-C file except gene_name and coordinate (which are in the front group)
    hic_compartment_cols = [c for c in merged.columns
                            if c not in (
                                ["gene_name", "coordinate",
                                 "padj", "pvalue", "log2FC", "baseMean",
                                 "description", "strand"]
                                + cutrun_diff_cols
                            )]

    # ── Build the explicit ordered column list ────────────────────────────────
    priority = (
        ["gene_name", "coordinate"]          # identifiers
        + ["padj", "pvalue", "log2FC", "baseMean"]   # RNA-seq statistics
        + ["strand"]                          # strand from CUT&RUN
        + cutrun_diff_cols                    # CUT&RUN epigenetic differences
        + hic_compartment_cols               # Hi-C PC1 and switch/shift columns
        + ["description"]                    # gene description at the end
    )

    # Keep only columns that actually exist in the merged DataFrame
    # (guards against any source file missing an expected column)
    front = [c for c in priority if c in merged.columns]

    # Append any remaining columns not captured above (future-proofing)
    rest  = [c for c in merged.columns if c not in front]

    return merged[front + rest]


def save_output(df, output_dir, output_file):
    """
    Save the merged DataFrame to a tab-delimited text file.
    Creates the output directory if it does not already exist.
    """
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, output_file)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"[INFO] Saved -> {out_path}")
    print(f"       Rows   : {df.shape[0]}")
    print(f"       Columns: {df.shape[1]}")
    print(f"       Column names: {list(df.columns)}")


def main(rnaseq_path, hic_path, cutrun_path, output_dir,
         output_file="merged_hic_rnaseq_cutrun_parent_bottom10.txt"):

    print("[INFO] Loading Hi-C compartment file...")
    hic = load_hic(hic_path)

    print("[INFO] Loading RNA-seq DESeq2 file...")
    rna = load_rnaseq(rnaseq_path)

    print("[INFO] Loading CUT&RUN difference file...")
    cr  = load_cutrun(cutrun_path)

    print("[INFO] Merging on gene_name (outer join)...")
    merged = merge_all(hic, rna, cr)

    save_output(merged, output_dir, output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Merge Hi-C compartment, RNA-seq DESeq2, and CUT&RUN difference "
            "files into a single tab-delimited multi-omics table."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python merge_multiomics.py \\
      outputs/tables/GeneExpressionComparisons_parentvsBottom10_orig.txt \\
      outputs/tables/A375-gene_CompartmentsAcrossConditions.txt \\
      outputs/tables/cut_run_signal_table_diff_only.txt \\
      --output_dir outputs/tables/
        """
    )
    parser.add_argument("rnaseq",
                        help="Path to DESeq2 output file (.txt, tab-delimited)")
    parser.add_argument("hic",
                        help="Path to Hi-C compartment file (.txt, tab-delimited)")
    parser.add_argument("cutrun",
                        help="Path to CUT&RUN diff-only file (.txt, tab-delimited)")
    parser.add_argument("--output_dir", default=".",
                        help="Output directory (default: current directory)")
    parser.add_argument("--output_file",
                        default="merged_hic_rnaseq_cutrun_parent_bottom10.txt",
                        help="Output filename")

    args = parser.parse_args()

    main(args.rnaseq, args.hic, args.cutrun,
         args.output_dir, args.output_file)