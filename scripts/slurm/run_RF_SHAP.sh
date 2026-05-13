#!/bin/bash
#SBATCH --job-name=RF_SHAP
#SBATCH -A ACF-UTK0009
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/rf_shap_%j.out
#SBATCH --error=logs/rf_shap_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

module load anaconda3/2024.06
source $ANACONDA3_SH

MERGED="data/processed/merged_hic_rnaseq_cutrun_parent_bottom10.txt"
CUTRUN="data/processed/cut_run_signal_table.txt"
FIGS="results/exploratory"

#conda activate multiomics-py

#echo "=========================================="
#echo " Q1: Epigenetic priming"
#echo "=========================================="
#./scripts/analysis/RF_SHAP.py \
#    $CUTRUN $MERGED \
#    --outdir $FIGS/Q1_RF_SHAP

#conda deactivate

conda activate multiomics-r

./scripts/analysis/PLSDA.R \
    $FIGS/Q1_RF_SHAP/Q1_switch_classified.txt \
    $FIGS/Q1_PLSDA

conda deactivate

echo "=========================================="
echo " All analyses complete."
echo " Figures in: $FIGS"
echo "=========================================="