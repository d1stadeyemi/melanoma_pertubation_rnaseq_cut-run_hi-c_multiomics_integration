#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH -A ACF-UTK0011
#SBATCH --qos=campus
#SBATCH --partition=campus
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --job-name=multiomics_RF_SHAP
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=ALL

module load anaconda3/2024.06
source $ANACONDA3_SH

conda activate shap_env

MERGED="outputs/tables/merged_hic_rnaseq_cutrun_parent_bottom10.txt"
CUTRUN="outputs/tables/cut_run_signal_table.txt"
FIGS="raw/tempo/"

echo "=========================================="
echo " Q1: Epigenetic priming"
echo "=========================================="
./scripts/RF_SHAP.py \
    $CUTRUN $MERGED \
    --outdir $FIGS/Q1_RF_SHAP

conda deactivate

conda activate r_analysis

./scripts/PLSDA.R \
    $FIGS/Q1_RF_SHAP/Q1_switch_classified.txt \
    $FIGS/Q1_PLSDA

conda deactivate

echo "=========================================="
echo " All analyses complete."
echo " Figures in: $FIGS"
echo "=========================================="