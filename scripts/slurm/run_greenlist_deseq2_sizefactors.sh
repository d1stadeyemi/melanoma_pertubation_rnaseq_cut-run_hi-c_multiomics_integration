#!/bin/bash
#SBATCH --job-name=greenlist_deseq2_diff
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/greenlist_deseq2_diff_%j.out
#SBATCH --error=logs/greenlist_deseq2_diff_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=all
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

module load anaconda3/2024.06
source $ANACONDA3_SH

#conda activate multiomics-py

#python scripts/analysis/greenlist_deseq2_sizefactors.py \
#  --sheet data/processed/sample_sheet.csv \
#  --mark H3K27me3 \
#  --greenlist data/processed/greenlist/greenlist_hg19.bed \
#  --out results/exploratory/greenlist_deseq2/H3K27me3_all \
#  --samples Control-H3K27me3-R1,Control-H3K27me3-R2,Bottom-H3K27me3-R1,Bottom-H3K27me3-R2,Parental-UNC-H3K27Me3-R1,Parental-UNC-H3K27Me3-R2,B10-UNC-H3K27Me3-R1,B10-UNC-H3K27Me3-R2

#conda deactivate

conda activate multiomics-r

# Parental UNC reversal, both units -> peaks/ and windows/ subdirs
Rscript scripts/analysis/greenlist_deseq2_diffbind.R \
  --sheet data/processed/sample_sheet.csv --mark H3K27me3 \
  --greenlist data/processed/greenlist/greenlist_hg19.bed \
  --group1 "Parental|UNC" --group2 "Parental|none" \
  --drop Parental-H3K27me3-R3 \
  --unit peaks,windows --peaks sheet \
  --blacklist data/processed/greenlist/hg19-blacklist.v2.bed \
  --outdir results/exploratory/diff_deseq2/H3K27me3_ParentalUNC_vs_Parentalnone
#--size_factors results/exploratory/greenlist_deseq2/H3K27me3_all_size_factors.tsv \

# B10 UNC
Rscript scripts/analysis/greenlist_deseq2_diffbind.R \
  --sheet data/processed/sample_sheet.csv --mark H3K27me3 \
  --greenlist data/processed/greenlist/greenlist_hg19.bed \
  --group1 "B10|UNC" --group2 "B10|none" \
  --unit peaks,windows --peaks sheet \
  --blacklist data/processed/greenlist/hg19-blacklist.v2.bed \
  --outdir results/exploratory/diff_deseq2/H3K27me3_B10UNC_vs_B10none
#--size_factors results/exploratory/greenlist_deseq2/H3K27me3_all_size_factors.tsv \

conda deactivate