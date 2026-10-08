#!/bin/bash
#SBATCH --job-name=gene_compartment
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=logs/gene_compartment_%j.out
#SBATCH --error=logs/gene_compartment_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

module load anaconda3/2024.06
source $ANACONDA3_SH

conda activate multiomics-py

python scripts/analysis/gene_compartment_analysis_4.py \
        --step1   results/exploratory/gene_compartment/Step1_bin_gene_density.txt \
        --cutrun  data/processed/cut_run_signal_table_250kb_temp.txt \
        --step3   results/exploratory/gene_compartment/Step3_gene_expression_joined.txt \
        --outdir  results/exploratory/gene_compartment_temp

conda deactivate
