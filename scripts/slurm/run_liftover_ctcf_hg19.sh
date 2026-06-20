#!/bin/bash
#SBATCH --job-name=ryder
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/ryder_%j.out
#SBATCH --error=logs/ryder_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

set -euo pipefail

module load anaconda3/2024.06
source "$ANACONDA3_SH"
conda activate multiomics-r

Rscript scripts/preprocessing/liftover_ctcf_hg19.R \
    --bed ryder/refData/hg38_CTCF_constitutive.bed \
    --outdir data/processed/ryder \
    --out data/processed/ryder/CTCF_constitutive_hg19.bed \
    --chain data/processed/greenlist/hg38ToHg19.over.chain.gz \
    --blacklist data/processed/greenlist/hg19-blacklist.v2.bed \
    --min_width 200

conda deactivate