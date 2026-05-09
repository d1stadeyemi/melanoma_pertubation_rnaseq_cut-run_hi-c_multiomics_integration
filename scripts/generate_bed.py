import pandas as pd
import argparse
import os

def generate_bed(hg19, upstream=2000, downstream=5000):
    hg19 = pd.read_excel(hg19)
    bed = pd.DataFrame()
    bed["chr"] = hg19["chr"]
    bed["start"] = (hg19["TSS"] - upstream).clip(lower=0).astype(int)
    bed["end"] = (hg19["TSS"] + downstream).astype(int)
    bed["gene_id"] = hg19["Gene"]
    bed["gene_name"] = hg19["Gene"]
    bed["strand"] = "."
    
    bed["gene_id"] = bed["gene_id"].fillna(
            bed["chr"] + ":" + bed["start"].astype(str) + "-" + bed["end"].astype(str))

    bed["gene_name"] = bed["gene_name"].fillna(bed["gene_id"])
    bed.drop_duplicates(subset="gene_id", inplace=True)
    bed.sort_values(["chr", "start"], inplace=True)
    
    return bed

def save_bed(bed, output_file, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    bed.to_csv(os.path.join(output_dir, output_file), sep="\t", header=False, index=False)

def save_txt(bed, output_file, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    bed.to_csv(os.path.join(output_dir, output_file), sep="\t", header=True, index=False)

def main(input_file, output_dir, upstream=2000, downstream=5000):
    bed = generate_bed(input_file, upstream=upstream, downstream=downstream)
    save_bed(bed, "GeneListTSShg19_2kb_up_5kb_down.bed", output_dir)
    save_txt(bed, "GeneListTSShg19_2kb_up_5kb_down.txt", output_dir)

    print(f"Bed file saved to {output_dir}")
    print(f"Columns: {list(bed.columns)}")
    print(f"Rows: {len(bed)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate BED and txt files from hg19 TSS data")
    parser.add_argument("hg19", help="Path to hg19 TSS Excel file")
    parser.add_argument("--upstream", type=int, default=2000, help="Upstream distance from TSS (default: 2000)")
    parser.add_argument("--downstream", type=int, default=5000, help="Downstream distance from TSS (default: 5000)")
    parser.add_argument("--output_dir", type=str, default=".", help="Output directory for the BED file (default: .)")

    args = parser.parse_args()
    
    main(args.hg19, args.output_dir, upstream=args.upstream, downstream=args.downstream)