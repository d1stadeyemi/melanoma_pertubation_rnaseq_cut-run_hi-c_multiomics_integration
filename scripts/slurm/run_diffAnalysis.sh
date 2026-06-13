#!/bin/bash
#SBATCH --job-name=DiffAnalysis
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=logs/diff_analysis_int_%j.out
#SBATCH --error=logs/diff_analysis_int_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

module load anaconda3/2024.06
source $ANACONDA3_SH

conda activate multiomics-r

# Run differential analysis per mark/contrast:
Rscript scripts/analysis/csaw.R \
    --sheet data/processed/sample_sheet.csv \
    --outdir results/exploratory/diff_analysis/csaw_no_batch/H3K27me3_B10unc \
    --mark H3K27me3 \
    --group1 "B10|UNC" --group2 "B10|none" \
    --peaks sheet \
    --baseline_group "Parental|none" \
    --spikein sheet
# Default: --test contrast

# Run differential analysis per mark/contrast:
Rscript scripts/analysis/csaw.R \
    --sheet data/processed/sample_sheet.csv \
    --outdir results/exploratory/diff_analysis/csaw_no_batch/H3K27me3_Parunc \
    --mark H3K27me3 \
    --group1 "Parental|UNC" --group2 "Parental|none" \
    --peaks sheet \
    --baseline_group "B10|none" \
    --spikein sheet
# Default: --test contrast

# Run differential analysis per mark/interaction:
Rscript scripts/analysis/csaw.R \
    --sheet data/processed/sample_sheet.csv \
    --outdir results/exploratory/diff_analysis/csaw_no_batch/H3K27me3_interaction \
    --mark H3K27me3 \
    --test interaction \
    --peaks sheet \
    --spikein sheet
# Default reference levels:
# --cond_ref Parental --cond_alt B10 --treat_ref none --treat_alt UNC 

# Rscript scripts/analysis/csaw.R \
#    --sheet data/processed/sample_sheet.csv \
#    --outdir results/exploratory/diff_analysis/csaw \
#    --mark H3K27ac \
#    --group1 "B10|none" --group2 "Parental|none" \
#    --window 500 --spacing 250 --min_lfc_bg 1

conda deactivate