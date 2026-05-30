#!/bin/bash
#SBATCH --job-name=QC_Samples
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/qc_samples_%j.out
#SBATCH --error=logs/qc_samples_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

module load anaconda3/2024.06
source $ANACONDA3_SH

conda activate multiomics-r

Rscript scripts/analysis/qc_samples.R \
        --sheet data/processed/sample_sheet.csv \
        --outdir results/exploratory/qc_samples \
        --plot

Rscript scripts/analysis/qc_sample_rep_global.R \
        --sheet data/processed/sample_sheet.csv \
        --outdir results/exploratory/qc_sample_rep_global

conda deactivate