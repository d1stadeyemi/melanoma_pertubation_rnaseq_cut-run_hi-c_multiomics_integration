"""
=============================================================================
gene_compartment_analysis.py  — STEP 1: gene density per 250kb bin
=============================================================================

GOAL OF THE FULL ANALYSIS (built step by step)
  Connect 3D-genome compartment switching (Parental -> Bottom10) to its
  FUNCTIONAL and EPIGENETIC context at 250kb resolution:
    Step 1 (THIS FILE) : gene density per bin (genes mapped to PC1 bins).
    Step 2             : gene density across switch categories (KW + Dunn).
    Step 3             : expression change (DESeq2) by switch direction.
    Step 4             : histone-mark change overlaid on switching genes.

WHY COORDINATE-MERGE (not interval overlap) FOR STEP 1
  The hg19 gene list is ALREADY binned to the SAME 250kb grid as the Hi-C PC1
  table, and both carry an identical 'coordinate' string (e.g.
  'chr12:6500001-6750000'). When both sides live on the same grid, a
  coordinate-EQUALITY merge is simpler and MORE correct than an interval
  overlap (which would just reproduce the equality join while risking
  boundary off-by-one errors). Gene density is then simply: how many gene-list
  rows share each bin coordinate (many genes -> one bin).

BUILD PROVENANCE (verified before writing this)
  PC1 chr1 max = 249,250,000 ; chr8 = 146,364,022  -> hg19 (confirmed).
  Gene list GAPDH TSS = 6,643,585 ; TP53 = 7,590,868 -> hg19 (confirmed).
  PC1 covers chr1-22 (autosomes). Gene list adds chrX. chrX genes will NOT
  match any PC1 bin and are dropped by the inner merge -- this is CORRECT
  (no compartment calls exist for chrX), not data loss. The script reports
  exactly how many genes drop and why.

MERGE GUARD
  Like the earlier coordinate-format bug, a silent format mismatch (stray
  space, 0- vs 1-based start) would yield ~0% overlap and all-empty density.
  Step 1 MEASURES the overlap fraction and halts if it is implausibly low,
  rather than producing a plausible-looking but empty result.

USAGE
  python gene_compartment_analysis.py \
      --pc1       data/processed/merged_hic_rnaseq_cutrun_parent_bottom10.txt \
      --gene_list data/processed/hg19_gene_list_binned.txt \
      --outdir    results/exploratory/gene_compartment

OUTPUT (Step 1)
  <outdir>/Step1_bin_gene_density.txt   one row per PC1 bin:
       coordinate, chr, start, end, switch_cat, gene_count, gene_names
  Console: overlap diagnostics, chrX drop count, density summary by category.
"""

import argparse
import os
import sys

import pandas as pd
import numpy as np


# -----------------------------------------------------------------------------
# Switch classification — IDENTICAL logic to RF_SHAP.py so categories match
# across every analysis in the study (consistency for the talk/paper).
# -----------------------------------------------------------------------------
def classify_switch(df):
    """Classify each bin into B_to_A / A_to_B / Stable_A / Stable_B / None.
    Reads switch columns in priority order, falling back to the sign of the
    original Parental/Bottom10 eigenvalues for the stable categories. Operates
    row-wise; identical to the RF_SHAP classifier."""
    cats = []
    for _, row in df.iterrows():
        cat = None
        v = str(row.get("P->B10_Both_original", "")).strip()
        if v not in ["-", "", "nan", "None"]:
            if "toA" in v or "BtoA" in v:
                cat = "B_to_A"
            elif "toB" in v or "AtoB" in v:
                cat = "A_to_B"
        if cat is None:
            v = str(row.get("P->B10_Switch_Orig", "")).strip()
            if v not in ["-", "", "nan", "None"]:
                if v == "BtoA":
                    cat = "B_to_A"
                elif v == "AtoB":
                    cat = "A_to_B"
        if cat is None:
            v = str(row.get("P->B10_Shift_Orig", "")).strip()
            if v not in ["-", "", "nan", "None"]:
                if "TowardA" in v:
                    cat = "B_to_A"
                elif "TowardB" in v:
                    cat = "A_to_B"
        if cat is None:
            p_orig = row.get("Parental_orig", np.nan)
            b_orig = row.get("Bottom10_orig", np.nan)
            try:
                p_orig = float(p_orig)
                b_orig = float(b_orig)
                if p_orig > 0 and b_orig > 0:
                    cat = "Stable_A"
                elif p_orig < 0 and b_orig < 0:
                    cat = "Stable_B"
            except (ValueError, TypeError):
                pass
        cats.append(cat)
    return cats


def parse_args():
    p = argparse.ArgumentParser(
        description="Step 1: gene density per 250kb compartment bin.")
    p.add_argument("--pc1", required=True,
                   help="Hi-C PC1 / merged table with a 'coordinate' column "
                        "and the P->B10_*_original switch columns.")
    p.add_argument("--gene_list", required=True,
                   help="hg19 gene list, 250kb-binned, with 'coordinate' and "
                        "'gene_name' columns.")
    p.add_argument("--outdir", default="results/exploratory/gene_compartment")
    p.add_argument("--coord_col", default="coordinate",
                   help="Name of the shared coordinate column [default %default].")
    p.add_argument("--gene_name_col", default="gene_name",
                   help="Gene-symbol column in the gene list [default %default].")
    p.add_argument("--min_overlap", type=float, default=0.80,
                   help="Min fraction of AUTOSOMAL gene-list bins that must "
                        "match a PC1 bin, else halt (probable format mismatch). "
                        "Default 0.80.")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # --- Load -----------------------------------------------------------------
    print("[INFO] Loading PC1 / merged table...")
    pc1 = pd.read_csv(args.pc1, sep="\t", low_memory=False)
    print(f"  {len(pc1):,} bins x {len(pc1.columns)} columns")

    print("[INFO] Loading hg19 gene list...")
    genes = pd.read_csv(args.gene_list, sep="\t", low_memory=False)
    print(f"  {len(genes):,} genes x {len(genes.columns)} columns")

    cc = args.coord_col
    gn = args.gene_name_col
    for label, df in (("PC1", pc1), ("gene_list", genes)):
        if cc not in df.columns:
            sys.exit(f"[ERROR] '{cc}' column not in {label}. "
                     f"Columns: {list(df.columns)}")
    if gn not in genes.columns:
        sys.exit(f"[ERROR] '{gn}' column not in gene_list. "
                 f"Columns: {list(genes.columns)}")

    # --- Normalize coordinate strings defensively (whitespace only; we do NOT
    #     reformat, to avoid masking a genuine grid mismatch) -----------------
    pc1[cc] = pc1[cc].astype(str).str.strip()
    genes[cc] = genes[cc].astype(str).str.strip()

    # --- chrX (and any non-autosome) reporting --------------------------------
    # PC1 is autosomes only; the gene list adds chrX. Quantify the expected drop.
    def chrom_of_coord(s):
        return s.split(":")[0] if isinstance(s, str) and ":" in s else None
    genes["_chrom"] = genes[cc].map(chrom_of_coord)
    pc1_chroms = set(pc1[cc].map(chrom_of_coord).dropna().unique())
    gene_chroms = set(genes["_chrom"].dropna().unique())
    extra_chroms = sorted(gene_chroms - pc1_chroms)
    if extra_chroms:
        n_extra = genes["_chrom"].isin(extra_chroms).sum()
        print(f"[INFO] Gene list has chromosomes absent from PC1: {extra_chroms}")
        print(f"  {n_extra:,} genes on these (e.g. chrX) will drop on merge "
              f"(no compartment calls exist there — expected, not data loss).")

    # --- Overlap guard (on AUTOSOMAL bins, where a match is expected) ---------
    pc1_coords = set(pc1[cc])
    autosomal_gene_coords = set(genes.loc[~genes["_chrom"].isin(extra_chroms), cc])
    if len(autosomal_gene_coords) == 0:
        sys.exit("[ERROR] No autosomal gene coordinates found; check --coord_col.")
    matched = autosomal_gene_coords & pc1_coords
    frac = len(matched) / len(autosomal_gene_coords)
    print(f"[INFO] Autosomal gene-bin coordinates matching a PC1 bin: "
          f"{len(matched):,}/{len(autosomal_gene_coords):,} ({frac:.1%})")
    if frac < args.min_overlap:
        sys.exit(
            f"[ERROR] Only {frac:.1%} of autosomal gene-bin coordinates match a "
            f"PC1 bin (below --min_overlap {args.min_overlap:.0%}). Likely a "
            f"coordinate-FORMAT mismatch (spacing, or 0- vs 1-based starts), "
            f"NOT biology.\n"
            f"  Example gene coords: {list(autosomal_gene_coords)[:3]}\n"
            f"  Example PC1 coords:  {list(pc1_coords)[:3]}\n"
            f"  Fix the formatting so both use identical 'chrN:start-end' "
            f"strings, then re-run."
        )

    # --- Classify switches on the PC1 table -----------------------------------
    print("[INFO] Classifying compartment switches (RF_SHAP logic)...")
    pc1 = pc1.copy()
    pc1["switch_cat"] = classify_switch(pc1)
    n_classified = pc1["switch_cat"].notna().sum()
    print(f"  {n_classified:,}/{len(pc1):,} bins classified")
    print(pc1["switch_cat"].value_counts(dropna=False).to_string())

    # --- Gene density per bin: count gene-list rows per coordinate ------------
    # Many genes -> one bin. Restrict to genes whose bin exists in PC1 (inner).
    genes_in_pc1 = genes[genes[cc].isin(pc1_coords)].copy()
    n_dropped = len(genes) - len(genes_in_pc1)
    print(f"[INFO] Genes mapped to a PC1 bin: {len(genes_in_pc1):,} "
          f"(dropped {n_dropped:,}: chrX + any unmatched).")

    density = (genes_in_pc1
               .groupby(cc)
               .agg(gene_count=(gn, "size"),
                    gene_names=(gn, lambda s: ";".join(map(str, s))))
               .reset_index())

    # --- Attach density to EVERY PC1 bin (bins with no gene -> count 0) -------
    keep_pc1_cols = [cc, "chr", "start", "end", "switch_cat"]
    keep_pc1_cols = [c for c in keep_pc1_cols if c in pc1.columns]
    out = pc1[keep_pc1_cols].merge(density, on=cc, how="left")
    out["gene_count"] = out["gene_count"].fillna(0).astype(int)
    out["gene_names"] = out["gene_names"].fillna("")

    # --- Report density by switch category (preview of Step 2) ----------------
    print("\n[INFO] Gene density (genes per 250kb bin) by switch category:")
    summ = (out[out["switch_cat"].notna()]
            .groupby("switch_cat")["gene_count"]
            .agg(n_bins="size", mean_genes="mean", median_genes="median",
                 total_genes="sum"))
    print(summ.to_string())

    # --- Save -----------------------------------------------------------------
    out_path = os.path.join(args.outdir, "Step1_bin_gene_density.txt")
    out.to_csv(out_path, sep="\t", index=False)
    print(f"\n[INFO] Step 1 output -> {out_path}")
    print("[INFO] Step 1 complete. Inspect the density-by-category summary "
          "above; if it looks sensible (expect Stable_A > Stable_B), proceed "
          "to Step 2 (Kruskal-Wallis + Dunn).")


if __name__ == "__main__":
    main()