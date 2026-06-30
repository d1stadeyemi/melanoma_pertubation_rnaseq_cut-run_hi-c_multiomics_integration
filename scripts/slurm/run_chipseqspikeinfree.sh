#!/bin/bash
#SBATCH --job-name=spikefree
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=logs/spikefree_%j.out
#SBATCH --error=logs/spikefree_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

module load anaconda3/2024.06
source $ANACONDA3_SH

conda activate multiomics-r

Rscript scripts/preprocessing/run_chipseqspikeinfree.R \
    --sheet  data/processed/sample_sheet.csv \
    --mark migrated \
    --treatment dmso,unc \
    --outdir data/processed/spikefree/H3K27Me3_migrated \
    --blacklist data/processed/greenlist/hg19-blacklist.v2.bed \
    --reference Parental-DMSO-Migrated-R1 --qc 

conda deactivate