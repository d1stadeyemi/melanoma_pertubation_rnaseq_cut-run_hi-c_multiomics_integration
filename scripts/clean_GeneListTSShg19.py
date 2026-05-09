import pandas as pd
import argparse
import os 


def clean_hg19(hg19):
    hg19 = pd.read_excel(hg19)
    df = pd.DataFrame()
    df["gene_name"] = hg19["Gene"]
    df["chr"] = hg19["chr"]
    df["TSS"] = hg19["TSS"]
    df["bin"] = hg19["Bin"]
    df["bin_start"] = hg19["BinStart"]
    df["bin_end"] = hg19["BinEnd"]
    df["coordinate"] = df["chr"] + ":" + df["bin_start"].astype(str) + "-" + df["bin_end"].astype(str)
    df["index"] = hg19["Index"]
    
    df["gene_name"] = df["gene_name"].fillna(df["coordinate"])

    df.sort_values("chr", inplace=True)
    
    return df

def save_txt(df, output_file, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    df.to_csv(os.path.join(output_dir, output_file), sep="\t", header=True, index=False)

def main(input_file, output_dir):
    new_hg19 = clean_hg19(input_file)
    save_txt(new_hg19, "GeneListTSShg19_cleaned.txt", output_dir)

    print(f"Txt file saved to {output_dir}")
    print(f"Columns: {list(new_hg19.columns)}")
    print(f"Rows: {len(new_hg19)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate BED and txt files from hg19 TSS data")
    parser.add_argument("hg19", help="Path to hg19 TSS Excel file")
    parser.add_argument("--output_dir", type=str, default=".", 
                        help="Output directory for the txt file (default: current directory)")

    args = parser.parse_args()
    
    main(args.hg19, args.output_dir)