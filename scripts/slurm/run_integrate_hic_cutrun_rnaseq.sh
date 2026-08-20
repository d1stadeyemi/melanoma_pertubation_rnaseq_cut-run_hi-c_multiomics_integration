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
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

module load anaconda3/2024.06
source $ANACONDA3_SH

conda activate multiomics-py

python scripts/analysis/integrate_hic_cutrun_rnaseq.py \
  --hic data/processed/compartment_change_18AUG2026.txt \
  --hic_coord_col Coord --hic_chr_col Chr \
  --hic_delta_col "Bottom10-Parental_DMSO" \
  --hic_category_col "Bottom10-Parental_DMSO_category" \
  --restrict_to "Switch AtoB" "Switch BtoA"
  --cutrun data/processed/cut_run_signal_table_250kb.eps_fixed.txt \
  --cutrun_coord_col extra_1 \
  --mark_cols log2ratio_p-b10_h3k27me3_mean \
  --mark_names H3K27me3 \
  --deseq2 <your DESeq2 results.txt> \
  --gene_map data/processed/GeneListTSShg19_2kb_up_5kb_down.txt \
  --outdir results/exploratory/integration_hic_cutrun_rnaseq