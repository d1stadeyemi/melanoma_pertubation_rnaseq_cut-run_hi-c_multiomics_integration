#!/bin/bash
#SBATCH --job-name=selective_erasure
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/selective_erasure_%j.out
#SBATCH --error=logs/selective_erasure_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=all
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

module load anaconda3/2024.06
source $ANACONDA3_SH

conda activate multiomics-py

python scripts/analysis/unc_selective_erasure.py \
    --sample-sheet  data/processed/sample_sheet.csv \
    --size-factors results/exploratory/greenlist_deseq2/H3K27me3_all_size_factors.tsv \
    --compartments  data/processed/A375-CompartmentsAcrossConditions.txt \
    --blacklist     data/processed/greenlist/hg19-blacklist.v2.bed \
    --norm  greenlist --min-baseline-quantile 0.25 \
    --outdir   results/exploratory/unc_selective_erasure \
    --deeptools-blacklist 

conda deactivate