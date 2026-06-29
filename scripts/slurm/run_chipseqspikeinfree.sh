#!/bin/bash
#SBATCH --job-name=spikefree_epz
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=logs/spikefree_epz_%j.out
#SBATCH --error=logs/spikefree_epz_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

module load anaconda3/2024.06
source $ANACONDA3_SH

conda activate multiomics-r

Rscript scripts/preprocessing/run_chipseqspikeinfree.R \
    --sheet  data/processed/sample_sheet_temp.csv \
    --mark H3K27me3 \
    --treatment dmso,epz \
    --outdir data/processed/spikefree/epz \
    --blacklist data/processed/greenlist/hg19-blacklist.v2.bed \
    --reference Parental-DMSO-H3K27me3-R1 --qc 

conda deactivate