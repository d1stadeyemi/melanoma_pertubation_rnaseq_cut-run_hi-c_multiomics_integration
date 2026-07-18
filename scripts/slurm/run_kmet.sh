#!/bin/bash
#SBATCH -A ACF-UTK0009
#SBATCH -p condo-rmccord
#SBATCH -q condo
#SBATCH -J kmet
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 23:00:00
#SBATCH -o logs/kmet_%j.out
#SBATCH -e logs/kmet_%j.err

set -euo pipefail
module purge
unset PYTHONPATH
module load anaconda3/2024.06
source $ANACONDA3_SH
conda activate multiomics-py

cd "$SLURM_SUBMIT_DIR"
scripts/preprocessing/kmet_spikein_count.sh \
    --manifest-csv data/raw/manifest.csv \
    --root data/raw/raw_dirs \
    --trimdir trimmed \
    --factor H3K27me3,MIGRATED \
    --csv-target H3K27me3 \
    --outdir data/processed/kmet_spikein \
    --threads 8

conda deactivate