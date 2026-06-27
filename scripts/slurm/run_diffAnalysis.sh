#!/bin/bash
#SBATCH --job-name=diffAnalysis
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=logs/diffAnalysis_%j.out
#SBATCH --error=logs/diffAnalysis_%j.err
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
#Rscript scripts/analysis/csaw.R \
#    --sheet data/processed/sample_sheet.csv \
#    --outdir results/exploratory/diff_analysis/csaw/Migrated_B10unc \
#    --mark migrated \
#    --group1 "B10|UNC" --group2 "B10|DMSO" \
#    --norm background,spikein,sf \
#    --spikein sheet \
#    --normfactors data/processed/spikefree/H3K27Me3_migrated/Migrated_normfactors.csv \
#    --blacklist data/processed/greenlist/hg19-blacklist.v2.bed \
#    --peaks sheet \
#    --baseline_group "Parental|DMSO" \
#    --genes data/processed/GeneListTSShg19_cleaned.txt 
#    --greenlist data/processed/greenlist/greenlist_hg19.bed
# Default: --test contrast

# Run differential analysis per mark/contrast:
Rscript scripts/analysis/csaw_v7.R \
    --sheet data/processed/sample_sheet.csv \
    --outdir results/exploratory/diff_analysis/csaw_v7/Migrated_Parunc_3FC \
    --mark migrated \
    --group1 "Parental|UNC" --group2 "Parental|DMSO" \
    --norm background,spikein,sf --window 1000,2000,5000 \
    --spikein sheet --filter_mode group_aware \
    --blacklist data/processed/greenlist/hg19-blacklist.v2.bed \
    --normfactors data/processed/spikefree/H3K27Me3_migrated/Migrated_normfactors.csv \
    --peaks sheet \
    --baseline_group "B10|DMSO" 
#    --genes data/processed/GeneListTSShg19_cleaned.txt
#    --greenlist data/processed/greenlist/greenlist_hg19.bed
# Default: --test contrast

#Rscript scripts/analysis/csaw.R \
#    --sheet data/processed/sample_sheet.csv \
#    --outdir results/exploratory/diff_analysis/csaw/Chronic_B8unc \
#    --mark chronic \
#    --group1 "B8|UNC" --group2 "B8|DMSO" \
#    --norm background,spikein,sf \
#    --spikein sheet \
#    --blacklist data/processed/greenlist/hg19-blacklist.v2.bed \
#    --normfactors data/processed/spikefree/H3K27Me3_chronic/Chronic_normfactors.csv \
#    --genes data/processed/GeneListTSShg19_cleaned.txt

# Run differential analysis per mark/interaction:
#Rscript scripts/analysis/csaw.R \
#    --sheet data/processed/sample_sheet.csv \
#    --outdir results/exploratory/diff_analysis/csaw/Migrated_H3K27me3_interaction \
#    --mark migrated \
#    --test interaction \
#    --cond_ref Parental --cond_alt B10 --treat_ref DMSO --treat_alt UNC \
#    --norm background,spikein,sf \
#    --spikein sheet \
#    --blacklist data/processed/greenlist/hg19-blacklist.v2.bed \
#    --normfactors data/processed/spikefree/H3K27Me3_migrated/Migrated_normfactors.csv \
#    --peaks sheet \
#    --genes data/processed/GeneListTSShg19_cleaned.txt
#    --greenlist data/processed/greenlist/greenlist_hg19.bed

# Rscript scripts/analysis/csaw.R \
#    --sheet data/processed/sample_sheet.csv \
#    --outdir results/exploratory/diff_analysis/csaw \
#    --mark H3K27ac \
#    --group1 "B10|none" --group2 "Parental|none" \
#    --window 500 --spacing 250 --min_lfc_bg 1

#Rscript scripts/analysis/csaw_overlap_venn.R \
#       --a results/exploratory/diff_analysis/csaw/Migrated_Parunc/background/migrated_ParentalUNC_vs_ParentalDMSO_background_all_regions.csv \
#       --b results/exploratory/diff_analysis/csaw/Migrated_B10unc/background/migrated_B10UNC_vs_B10DMSO_background_all_regions.csv \
#       --a-name Parental --b-name B10 \
#       --up-label gained --down-label lost \
#       --genes data/processed/GeneListTSShg19_cleaned.txt \
#       --outdir results/exploratory/diff_analysis/csaw/overlap_migrated_UNC_Parental_vs_B10

conda deactivate