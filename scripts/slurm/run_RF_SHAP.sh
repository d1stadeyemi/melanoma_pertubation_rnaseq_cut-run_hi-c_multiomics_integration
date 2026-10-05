#!/bin/bash
#SBATCH --job-name=RF_SHAP_controls
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/rf_shap_controls_%j.out
#SBATCH --error=logs/rf_shap_controls_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=all
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

module load anaconda3/2024.06
source $ANACONDA3_SH

HIC="data/processed/A375-CompartmentsAcrossConditions.txt"
#CUTRUN="data/processed/cut_run_signal_table_250kb.txt"
CUTRUN="data/processed/cut_run_signal_table_250kb_temp.txt"
FIGS="results/exploratory"

conda activate multiomics-ml

echo "=========================================="
echo " Epigenetic priming"
echo "=========================================="
python scripts/analysis/RF_SHAP.py \
    $CUTRUN $HIC \
    --outdir $FIGS/RF_SHAP_both_temp \
    --comparison both \
    --stats mean \
    --model rf
# --shap_dependence_partners H3K27me3_mean:H3K27Ac_mean H3K9me3_mean:H3K4me1_mean H3K27me3_mean:H3K4me1_mean H3K27me3_mean:H3K9me3_mean \

# Permutation test for RF and SHAP significance
#python scripts/analysis/RF_permutation_test.py \
#    $FIGS/RF_SHAP_both/baseline/Switch_classified.txt \
#    --outdir $FIGS/RF_SHAP_both/baseline_permutation_2000_test \
#    --comparison both \
#    --model rf \
#    --n_perm 2000

# Control analyses
# python scripts/analysis/RF_SHAP_controls.py \
#    $CUTRUN $HIC \
#    --outdir $FIGS/RF_SHAP_both/controls \
#    --comparison both \
#    --model rf

conda deactivate

#conda activate multiomics-r

#Rscript scripts/analysis/PLSDA.R \
#    $FIGS/RF_SHAP_both/baseline/Switch_classified.txt \
#    $FIGS/PLSDA

#conda deactivate

echo "=========================================="
echo " all analyses complete."
echo " Figures in: $FIGS"
echo "=========================================="
