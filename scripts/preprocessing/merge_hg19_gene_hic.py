import pandas as pd
import numpy as np
import os 
import argparse

def add_gene_name(hg19, hiC_PC1):
    hg19 = pd.read_csv(hg19, sep="\t", header=0)
    hiC_PC1 = pd.read_csv(hiC_PC1, sep="\t", header=0)

    hic_gene = hiC_PC1.merge(
        hg19[["gene_name", "coordinate"]], 
        on="coordinate", how="left"
    )

    hic_gene.insert(0, "gene_name", hic_gene.pop("gene_name"))
    hic_gene["gene_name"] = hic_gene["gene_name"].fillna(hic_gene["coordinate"])
    hic_gene.drop_duplicates(inplace=True)

    return hic_gene

def save_txt(hic_gene, output_file, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    hic_gene.to_csv(os.path.join(output_dir, output_file), sep="\t", index=False)

def main(hg19, hiC_PC1, output_dir):
    hic_gene = add_gene_name(hg19, hiC_PC1)
    save_txt(hic_gene, "A375-gene_CompartmentsAcrossConditions.txt", output_dir)

    print(f"Output txt file saved to {output_dir} as A375-gene_CompartmentsAcrossConditions.txt")    
    print(f"Rows: {len(hic_gene)}")
    print(f"Columns: {hic_gene.shape[1]}")
    print(f"Column list: {list(hic_gene.columns)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Adds gene column to Hi-C PC1 data")
    parser.add_argument("hg19", help="Path to the hg19 data txt file which contains gene list")
    parser.add_argument("hiC_PC1", help="Path to the input txt file containing Hi-C PC1 data")
    parser.add_argument("--output_dir", type=str, default=".",
                        help="Output directory for the txt file (default: current directory)")
    
    args = parser.parse_args()

    main(args.hg19, args.hiC_PC1, args.output_dir)