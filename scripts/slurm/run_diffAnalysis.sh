#!/bin/bash
#SBATCH --job-name=DiffAnalysis
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:50:00
#SBATCH --output=logs/diff_analysis_%j.out
#SBATCH --error=logs/diff_analysis_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

module load anaconda3/2024.06
source $ANACONDA3_SH

conda activate multiomics-r

# Run differential analysis per mark/contrast:
Rscript scripts/analysis/csaw.R \
    --sheet data/processed/sample_sheet.csv \
    --outdir results/exploratory/diff_analysis/csaw \
    --mark H3K4me1 \
    --group1 "B10|none" --group2 "Parental|none"

Rscript scripts/analysis/csaw.R \
    --sheet data/processed/sample_sheet.csv \
    --outdir results/exploratory/diff_analysis/csaw \
    --mark H3K9me3 \
    --group1 "B10|none" --group2 "Parental|none"

Rscript scripts/analysis/csaw.R \
    --sheet data/processed/sample_sheet.csv \
    --outdir results/exploratory/diff_analysis/csaw \
    --mark H3K27ac \
    --group1 "B10|none" --group2 "Parental|none" \
#    --window 500 --spacing 250 --min_lfc_bg 1

conda deactivate