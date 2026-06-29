#!/bin/bash
#SBATCH --job-name=enrichment
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
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
    --genelists results/exploratory/diff_analysis/csaw/Migrated_Parunc/spikein/genelists \
    --background results/exploratory/diff_analysis/csaw/Migrated_Parunc/spikein/genelists/background_universe.txt \
    --databases GO_BP,Reactome,Hallmark,PRC2 --max_gs 2000 \
    --outdir results/exploratory/diff_analysis/csaw/Migrated_Parunc/spikein/genelists_enrichment

Rscript scripts/analysis/pathway_enrichment.R \
    --genelists results/exploratory/diff_analysis/csaw/Migrated_B10unc/background/genelists \
    --background results/exploratory/diff_analysis/csaw/Migrated_B10unc/background/genelists/background_universe.txt \
    --databases GO_BP,Reactome,Hallmark,PRC2 --max_gs 2000 \
    --outdir results/exploratory/diff_analysis/csaw/Migrated_B10unc/background/genelists_enrichment

conda deactivate