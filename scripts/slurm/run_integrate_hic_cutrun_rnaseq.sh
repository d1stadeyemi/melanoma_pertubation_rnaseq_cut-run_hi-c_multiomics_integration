#!/bin/bash
#SBATCH --job-name=integrate
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --time=23:00:00
#SBATCH --output=logs/integrate_%j.out
#SBATCH --error=logs/integrate_%j.err
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

module load anaconda3/2024.06
source $ANACONDA3_SH

conda activate multiomics-py

python scripts/analysis/integrate_hic_cutrun_rnaseq_1.py \
  --hic data/processed/compartment_change_18AUG2026.txt \
  --hic_coord_col Coord --hic_chr_col Chr \
  --hic_delta_col "Bottom10-Parental_DMSO" \
  --hic_category_col "Bottom10-Parental_DMSO_category" \
  --restrict_to "Shift in A towardA" "Shift in B towardB" \
  --cutrun data/processed/cut_run_signal_table_250kb.eps_fixed.txt \
  --cutrun_coord_col extra_1 \
  --mark_cols log2ratio_p-b10_h3k9me3_mean \
  --mark_names H3K9me3 \
  --deseq2 /lustre/isaac24/proj/UTK0009/Taiwo/2026JUL-Taiwo-A375UNC-RNA/New_Deseq2_Output/BD_vs_PD/tables/20_results_Bottom-10-DMSO_vs_Parental-DMSO_full.csv \
  --deseq2_lfc_col log2FoldChange \
  --gene_map data/processed/GeneListTSShg19_2kb_up_5kb_down.txt \
  --deseq2_gene_col gene_name \
  --outdir results/exploratory/integration_hic_cutrun_rnaseq/H3K9me3/b10D_vs_parD_shiftin

python scripts/analysis/integrate_hic_cutrun_rnaseq_2.py \
  --bin_table results/exploratory/integration_hic_cutrun_rnaseq/H3K9me3/b10D_vs_parD_shiftin/integration_bin_table.txt \
  --mark_cols log2ratio_p-b10_h3k9me3_mean \
  --mark_names H3K9me3 \
  --outdir results/exploratory/integration_hic_cutrun_rnaseq/H3K9me3/b10D_vs_parD_shiftin/stratified
