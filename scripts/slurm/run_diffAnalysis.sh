#!/bin/bash
#SBATCH --job-name=csaw_migrated
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=logs/csaw_migrated_%j.out
#SBATCH --error=logs/csaw_migrated_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

module load anaconda3/2024.06
source $ANACONDA3_SH

conda activate multiomics-r

#Rscript scripts/preprocessing/prepare_greenlist_hg19.R \
#    --excel data/processed/greenlist/greenlist_hg38.xlsx \
#    --outdir data/processed/greenlist \
#    --out data/processed/greenlist/greenlist_hg19.bed \
#    --igg ../../Taiwo/20231024_Cut_Run/IgG-R1/aligned/dedup/IgG-R1.bam \
#    --blacklist data/processed/greenlist/hg19-blacklist.v2.bed \
#    --chain data/processed/greenlist/hg38ToHg19.over.chain.gz

# Run differential analysis per mark/contrast:
Rscript scripts/analysis/csaw_v6.R \
    --sheet data/processed/sample_sheet.csv \
    --outdir results/exploratory/diff_analysis/csaw_v6/Migrated_B10unc \
    --mark migrated \
    --group1 "B10|UNC" --group2 "B10|DMSO" \
    --norm libsize,background,spikein \
    --spikein sheet 
#    --normfactors data/processed/spikefree/H3K27me3_normfactors.csv \
#    --peaks sheet \
#    --baseline_group "Parental|none" \
#    --greenlist data/processed/greenlist/greenlist_hg19.bed
# Default: --test contrast

# Run differential analysis per mark/contrast:
Rscript scripts/analysis/csaw_v6.R \
    --sheet data/processed/sample_sheet.csv \
    --outdir results/exploratory/diff_analysis/csaw_v6/Migrated_Parunc \
    --mark migrated \
    --group1 "Parental|UNC" --group2 "Parental|DMSO" \
    --norm libsize,background,spikein \
    --spikein sheet 
#    --normfactors data/processed/spikefree/H3K27me3_normfactors.csv \
#    --peaks sheet \
#    --baseline_group "B10|none" \
#    --greenlist data/processed/greenlist/greenlist_hg19.bed
# Default: --test contrast

# Run differential analysis per mark/interaction:
Rscript scripts/analysis/csaw_v6.R \
    --sheet data/processed/sample_sheet.csv \
    --outdir results/exploratory/diff_analysis/csaw_v6/Migrated_H3K27me3_interaction \
    --mark migrated \
    --test interaction \
    --norm libsize,background,spikein \
    --spikein sheet 
#    --normfactors data/processed/spikefree/H3K27me3_normfactors.csv \
#    --peaks sheet \
#    --greenlist data/processed/greenlist/greenlist_hg19.bed
# Default reference levels:
# --cond_ref Parental --cond_alt B10 --treat_ref none --treat_alt UNC 

# Rscript scripts/analysis/csaw.R \
#    --sheet data/processed/sample_sheet.csv \
#    --outdir results/exploratory/diff_analysis/csaw \
#    --mark H3K27ac \
#    --group1 "B10|none" --group2 "Parental|none" \
#    --window 500 --spacing 250 --min_lfc_bg 1

conda deactivate