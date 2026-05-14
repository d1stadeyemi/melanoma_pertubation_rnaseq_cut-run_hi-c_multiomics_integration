#!/bin/bash
#SBATCH --job-name=RF_permutation_test
#SBATCH -A ACF-UTK0009
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=logs/rf_permutation_test_%j.out
#SBATCH --error=logs/rf_permutation_test_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

module load anaconda3/2024.06
source $ANACONDA3_SH

conda activate multiomics-ml

echo "=========================================="
echo " RF Permutation Test"
echo "=========================================="

python scripts/RF_permutation_test.py \
    results/exploratory/RF_SHAP/Switch_classified.txt 
# Other parameters (e.g. --n_perm, --n_jobs) have been set to their defaults in the script, 
# but can be added here if desired.

conda deactivate

echo "=========================================="
echo " All analyses complete."
echo " Figures in: results/exploratory/RF_permutation_test"
echo "=========================================="
