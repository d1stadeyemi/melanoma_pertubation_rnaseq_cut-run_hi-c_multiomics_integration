"""
=============================================================================
clean_hic_pc1.py — extract & tidy Hi-C PC1 compartment calls from the raw xlsx
=============================================================================

WHAT THIS DOES
  The raw Excel sheet (A375 compartments across conditions) has a two-row
  header and many condition blocks. This script pulls the columns we actually
  use, gives them clean, machine-safe names (no spaces/parentheses/arrows that
  break downstream string parsing), and writes a tab-separated table.

WHY SELECT BY NAME, NOT POSITION
  The previous version used positional slices (e.g. cols[26:30]). Positional
  slicing silently breaks the moment a column is inserted/removed or a merged-
  header cell shifts the count — you get a plausible-looking file with the
  WRONG data under each name. Here we select by the column's real (cleaned)
  header text and HALT if any requested column is missing, so a layout change
  fails loudly instead of corrupting the output.

COLUMN BLOCKS KEPT (raw header text -> output name)
  Core coordinates + raw PC1 signal:
    Coord                              -> coordinate
    Chr, Start, End                    -> chr, start, end
    Parental DMSO                      -> Parental_DMSO
    Parental UNC                       -> Parental_UNC
    Bottom10-DMSO (fixed inverted)     -> Bottom10-DMed
    Bottom10-UNC (fixed inverted)      -> Bottom10-Ued
    Control-R1                         -> Parental_orig     (baseline parental PC1)
    Bottom10-R1                        -> Bottom10_orig     (baseline B10 PC1)

  Parental DMSO -> UNC (drug switch in parental):
    UNC-DMSO Parental                  -> UNC-DMSO_Parental
    P DMSO->UNC switch / shift / Both  -> P_DMSO->UNC_switch / _shift / _Both

  Bottom10 DMSO -> UNC (drug switch in B10):
    UNC-DMSO Bottom10                  -> UNC-DMSO_Bottom10
    B10 DMSO->UNC switch / shift / Both-> B10_DMSO->UNC_switch / _shift / _Both

  Parental -> Bottom10, DMSO axis (migration switch, DMSO):
    B10-P (DMSO)                       -> B10-P_DMSO
    P->B10 Switch/Shift/Both (DMSO)    -> P->B10_Switch_DMSO / _Shift_DMSO / _Both_DMSO

  Parental -> Bottom10, original axis (migration switch, the one we classify on):
    B10-P (Orig)                       -> B10-P_Orig
    P->B10 Switch/Shift (Orig)         -> P->B10_Switch_Orig / _Shift_Orig
    P->B10 Both (original)             -> P->B10_Both_original

USAGE
  python clean_hic_pc1.py raw_compartments.xlsx --output_dir data/processed
  # optional: --sheet 0  --header_row 1  --output_name A375-CompartmentsAcrossConditions.txt
============================================================================
"""

import argparse
import os
import sys
import re
import pandas as pd


def _clean_name(name):
    """Make a header machine-safe WITHOUT destroying meaning:
    strip, collapse internal whitespace to '_', drop parentheses, keep the
    '->' arrow and hyphens (downstream switch-classifiers look for 'toA'/'toB'
    etc., so we preserve the arrow rather than mangling it)."""
    s = str(name).strip()
    s = s.replace("(", "").replace(")", "")     # drop parens
    s = re.sub(r"\s+", "_", s)                   # any whitespace run -> single _
    return s


def clean_hic_pc1(input_file, sheet, header_row):
    # The real header lives on the 2nd row of the sheet (header=1 by default),
    # because row 0 holds the merged block titles ("Parental DMSO->UNC", etc.).
    hic = pd.read_excel(input_file, sheet_name=sheet, header=header_row)

    # --- Map RAW header text -> desired OUTPUT name. We match on raw text so a
    #     column moving left/right does not matter; only its name must exist. ---
    # NOTE: raw keys must match the spreadsheet's header EXACTLY (including the
    # "(fixed inverted)" / "(DMSO)" / "(original)" suffixes and the '->' arrow).
    rename_map = {
        # core
        "Coord":                            "coordinate",
        "Chr":                              "chr",
        "Start":                            "start",
        "End":                              "end",
        "Parental DMSO":                    "Parental_DMSO",
        "Parental UNC":                     "Parental_UNC",
        "Bottom10-DMSO (fixed inverted)":   "Bottom10-DMSO",
        "Bottom10-UNC (fixed inverted)":    "Bottom10-UNC",
        "Control-R1":                       "Parental_orig",
        "Bottom10-R1":                      "Bottom10_orig",
        # Parental DMSO -> UNC (drug axis, parental)
        "UNC-DMSO Parental":                "UNC-DMSO_Parental",
        "P DMSO->UNC switch":               "P_DMSO->UNC_switch",
        "P DMSO->UNC shift":                "P_DMSO->UNC_shift",
        "P DMSO->UNC Both":                 "P_DMSO->UNC_Both",
        # Bottom10 DMSO -> UNC (drug axis, B10)
        "UNC-DMSO Bottom10":                "UNC-DMSO_Bottom10",
        "B10 DMSO->UNC switch":             "B10_DMSO->UNC_switch",
        "B10 DMSO->UNC shift":              "B10_DMSO->UNC_shift",
        "B10 DMSO->UNC Both":               "B10_DMSO->UNC_Both",
        # Parental -> B10, DMSO axis (migration, DMSO)
        "B10-P (DMSO)":                     "B10-P_DMSO",
        "P->B10 Switch (DMSO)":             "P->B10_Switch_DMSO",
        "P->B10 Shift (DMSO)":              "P->B10_Shift_DMSO",
        "P->B10 Both (DMSO)":               "P->B10_Both_DMSO",
        # Parental -> B10, original axis (migration, the classification axis)
        "B10-P (Orig)":                     "B10-P_Orig",
        "P->B10 Switch (Orig)":             "P->B10_Switch_Orig",
        "P->B10 Shift (Orig)":              "P->B10_Shift_Orig",
        "P->B10 Both (original)":           "P->B10_Both_original",
    }

    # --- Resolve raw keys robustly: build a lookup of cleaned-actual -> actual
    #     so trailing spaces / double spaces in the sheet don't cause misses. --
    actual = {col: col for col in hic.columns}
    norm_actual = {re.sub(r"\s+", " ", str(c)).strip(): c for c in hic.columns}

    selected, out_names, missing = [], [], []
    for raw_key, out_name in rename_map.items():
        norm_key = re.sub(r"\s+", " ", raw_key).strip()
        if raw_key in actual:
            src = raw_key
        elif norm_key in norm_actual:                 # tolerant match
            src = norm_actual[norm_key]
        else:
            missing.append(raw_key)
            continue
        selected.append(src)
        out_names.append(out_name)

    # --- Fail LOUDLY if any requested column was not found ---
    if missing:
        sys.exit(
            "[ERROR] These expected columns were not found in the sheet:\n  "
            + "\n  ".join(missing)
            + "\n\nAvailable columns are:\n  "
            + "\n  ".join(map(str, hic.columns))
            + "\n\nFix: check the header text/suffixes above, or adjust "
              "--header_row (the real header may be on a different row)."
        )

    out = hic[selected].copy()
    out.columns = out_names      # apply clean output names in the same order
    return out


def save_txt(df, output_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, output_file)
    df.to_csv(path, sep="\t", header=True, index=False)
    return path


def main():
    ap = argparse.ArgumentParser(
        description="Clean Hi-C PC1 compartment calls from the raw xlsx, "
                    "keeping core PC1 + drug-axis + migration-axis switch columns.")
    ap.add_argument("input_file", help="Raw Excel file with Hi-C PC1 / compartment calls.")
    ap.add_argument("--output_dir", default=".", help="Output directory [default: .].")
    ap.add_argument("--output_name", default="A375-CompartmentsAcrossConditions.txt",
                    help="Output filename [default: %(default)s].")
    ap.add_argument("--sheet", default=0,
                    help="Sheet name or index [default: 0].")
    ap.add_argument("--header_row", type=int, default=1,
                    help="0-indexed row holding the real column headers "
                         "(row 0 is the merged block titles) [default: 1].")
    args = ap.parse_args()

    # allow --sheet to be an int index or a string name
    try:
        sheet = int(args.sheet)
    except (ValueError, TypeError):
        sheet = args.sheet

    df = clean_hic_pc1(args.input_file, sheet, args.header_row)
    path = save_txt(df, args.output_name, args.output_dir)

    print(f"[INFO] Saved -> {path}")
    print(f"[INFO] Rows: {len(df)}")
    print(f"[INFO] Columns ({df.shape[1]}): {list(df.columns)}")


if __name__ == "__main__":
    main()