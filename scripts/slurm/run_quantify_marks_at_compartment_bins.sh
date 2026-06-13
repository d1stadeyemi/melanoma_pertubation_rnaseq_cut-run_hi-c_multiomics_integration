#!/bin/bash
#SBATCH --job-name=quantify_marks
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=logs/quantify_marks_%j.out
#SBATCH --error=logs/quantify_marks_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

module load anaconda3/2024.06
source $ANACONDA3_SH

conda activate multiomics-r

Rscript scripts/preprocessing/quantify_marks_at_compartment_bins.R \
    --sheet data/processed/sample_sheet.csv \
    --compartments data/processed/A375-CompartmentsAcrossConditions.txt \
    --coord_col coordinate \
    --expect_bins 10576 \
    --reversal_check results/exploratory/reversal/Reversal_classified.txt \
    --outdir data/processed/mark_signal_250kb

conda deactivate