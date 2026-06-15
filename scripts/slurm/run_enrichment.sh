#!/bin/bash
#SBATCH --job-name=enrichment
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/enrichment_%j.out
#SBATCH --error=logs/enrichment_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

module load anaconda3/2024.06
source $ANACONDA3_SH

conda activate multiomics-r

Rscript scripts/analysis/pathway_enrichment.R \
    --genelists results/exploratory/reversal/genelists \
    --background results/exploratory/reversal/genelists/background_universe.txt \
    --databases PRC2 --focal migration_specific_reversal \
    --max_gs 2000 --pvalue 1 --qvalue 1 \
    --outdir results/exploratory/enrichment_prc2_big

conda deactivate