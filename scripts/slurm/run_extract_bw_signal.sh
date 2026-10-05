#!/bin/bash
#SBATCH --job-name=extract_bw_signal
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/extract_bw_signal_%j.out
#SBATCH --error=logs/extract_bw_signal_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=all
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

module load anaconda3/2024.06
source $ANACONDA3_SH

conda activate multiomics-py

python scripts/preprocessing/extract_bigwig_signal.py \
    data/processed/coordinates_250kb_bins.bed \
    --bw data/processed/mean_bigwigs/*.mean.bw \
    --output_dir data/processed \
    --output_file cut_run_signal_table_250kb_temp.txt 
#    --log2

conda deactivate
