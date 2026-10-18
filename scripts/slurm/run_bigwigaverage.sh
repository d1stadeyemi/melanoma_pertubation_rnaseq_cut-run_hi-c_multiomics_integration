#!/bin/bash
#SBATCH --job-name=bigwigaverage
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/bigwigaverage_%j.out
#SBATCH --error=logs/bigwigaverage_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

set -euo pipefail

module load anaconda3/2024.06
source "$ANACONDA3_SH"
conda activate multiomics-py

bash scripts/preprocessing/average_replicate_bigwigs.sh \
    --sheet data/processed/sample_sheet_temp.csv \
    --outdir data/processed/mean_bigwigs \
    --condition B10 \
    --treatment none

conda deactivate
