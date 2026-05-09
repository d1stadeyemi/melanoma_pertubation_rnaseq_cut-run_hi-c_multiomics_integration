import pandas as pd
import numpy as np
import argparse
import os

def clean_hic_pc1(input_file):
    hic = pd.read_excel(input_file, header=1)
    
    # Keep only the first 8 columns and the specified columns for Control-R1 and Bottom10-R1
    columns_to_keep = list(hic.columns[:8]) + list(hic.columns[10:12]) + list(hic.columns[26:30]) + list(hic.columns[38:42])
    hic_cleaned = hic[columns_to_keep]

    # Clean column names
    hic_cleaned.columns = hic_cleaned.columns.str.strip().str.replace(" ", "_").str.replace("(", "").str.replace(")", "")
    
    # Rename columns for clarity
    hic_cleaned = hic_cleaned.rename(columns={"Chr": "chr", "Start": "start", "End": "end", "Coord" : "coordinate",
                                              "Control-R1": "Parental_orig", "Bottom10-R1": "Bottom10_orig"})
    
    # Save the cleaned DataFrame to a new file
    return hic_cleaned

def save_txt(hic_cleaned, output_file, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    hic_cleaned.to_csv(os.path.join(output_dir, output_file), sep="\t", header=True, index=False)

def main(input_file, output_dir):
    hic_cleaned = clean_hic_pc1(input_file)
    save_txt(hic_cleaned, "A375-CompartmentsAcrossConditions.txt", output_dir)

    print(f"Csv file saved to {output_dir} as A375-CompartmentsAcrossConditions.csv")
    print(f"Rows: {len(hic_cleaned)}")
    print(f"Columns: {hic_cleaned.shape[1]}")
    print(f"Columns list: {list(hic_cleaned.columns)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clean Hi-C PC1 data by keeping only relevant columns")
    parser.add_argument("input_file", help="Path to the input Excel file containing Hi-C PC1 data")
    parser.add_argument("--output_dir", type=str, default=".", 
                        help="Output directory for the txt file (default: current directory)")

    args = parser.parse_args()
    
    main(args.input_file, args.output_dir) 