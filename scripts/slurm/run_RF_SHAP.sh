#!/bin/bash
#SBATCH --job-name=RF_SHAP
#SBATCH -A ACF-UTK0009
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/rf_shap_%j.out
#SBATCH --error=logs/rf_shap_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

module load anaconda3/2024.06
source $ANACONDA3_SH

MERGED="data/processed/merged_hic_rnaseq_cutrun_parent_bottom10.txt"
CUTRUN="data/processed/cut_run_signal_table.txt"
FIGS="results/exploratory_2"

conda activate multiomics-ml

echo "=========================================="
echo " Epigenetic priming"
echo "=========================================="
python scripts/analysis/RF_SHAP.py \
    $CUTRUN $MERGED \
    --outdir $FIGS/RF_SHAP \
    --include_ctcf \
    --shap_dependence_partners H3K27Ac:H3K27me3 H3K27me3:H3K4me1 CTCF:H3K27Ac CTCF:H3K27me3


conda deactivate

#conda activate multiomics-r

#Rscript scripts/analysis/PLSDA.R \
#    $FIGS/RF_SHAP/Switch_classified.txt \
#    $FIGS/PLSDA

#conda deactivate

echo "=========================================="
echo " All analyses complete."
echo " Figures in: $FIGS"
echo "=========================================="