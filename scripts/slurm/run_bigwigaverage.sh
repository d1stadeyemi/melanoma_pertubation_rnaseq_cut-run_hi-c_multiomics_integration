#!/bin/bash
#SBATCH --job-name=bigwigaverage
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/bigwigaverage_%j.out
#SBATCH --error=logs/bigwigaverage_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

set -euo pipefail

module load anaconda3/2024.06
source "$ANACONDA3_SH"
conda activate multiomics-py

#BWDIR="data/raw/bigwigs"                       
#OUT="data/processed/ryder_output/avg_tracks"   

#mkdir -p "$OUT"

#bigwigAverage \
#    -b "$BWDIR"/Parental-H3K27me3-R1.cpm.norm.bw "$BWDIR"/Parental-H3K27me3-R2.cpm.norm.bw \
#    -o "$OUT"/Parental-H3K27Me3_avg_rep.rpm.bw -bs 50 -p "${SLURM_CPUS_PER_TASK:-8}"

#bigwigAverage \
#    -b "$BWDIR"/Parental-UNC-H3K27Me3-R1.cpm.norm.bw "$BWDIR"/Parental-UNC-H3K27me3-R2.cpm.norm.bw \
#    -o "$OUT"/Parental-UNC-H3K27Me3_avg_rep.rpm.bw -bs 50 -p "${SLURM_CPUS_PER_TASK:-8}"

#bigwigAverage \
#    -b "$BWDIR"/B10-H3K27me3-R1.cpm.norm.bw "$BWDIR"/B10-H3K27me3-R2.cpm.norm.bw  \
#    -o "$OUT"/B10-H3K27Me3_avg_rep.rpm.bw -bs 50 -p "${SLURM_CPUS_PER_TASK:-8}"

#bigwigAverage \
#    -b "$BWDIR"/B10-UNC-H3K27Me3-R1.cpm.norm.bw "$BWDIR"/B10-UNC-H3K27me3-R2.cpm.norm.bw \
#    -o "$OUT"/B10-UNC-H3K27Me3_avg_rep.rpm.bw -bs 50 -p "${SLURM_CPUS_PER_TASK:-8}"

#echo "Done. Tracks in $OUT:"
#ls -1 "$OUT"/*_avg_rep.rpm.bw

bigwigAverage \
    -b /lustre/isaac24/proj/UTK0009/Muhammad/RNAseq/rnaseq_output/A375_Control_R1_S1/A375_Control_R1_S1.cpm.norm.bw \
       /lustre/isaac24/proj/UTK0009/Muhammad/RNAseq/rnaseq_output/A375_Control_R3_S2/A375_Control_R3_S2.cpm.norm.bw \
    -o /lustre/isaac24/proj/UTK0009/Muhammad/RNAseq/rnaseq_output/avg_tracks/A375_Control_avg_reps.cpm.norm.bw \
    -bs 10 -p 8

bigwigAverage \
    -b /lustre/isaac24/proj/UTK0009/Muhammad/RNAseq/rnaseq_output/Bottom10_R2_S9/Bottom10_R2_S9.cpm.norm.bw \
       /lustre/isaac24/proj/UTK0009/Muhammad/RNAseq/rnaseq_output/Bottom10_R3_S10/Bottom10_R3_S10.cpm.norm.bw \
    -o /lustre/isaac24/proj/UTK0009/Muhammad/RNAseq/rnaseq_output/avg_tracks/Bottom10_avg_reps.cpm.norm.bw \
    -bs 10 -p 8

conda deactivate