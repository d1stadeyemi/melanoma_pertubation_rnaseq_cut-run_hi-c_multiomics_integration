#!/bin/bash
###############################################################################
# Batch submitter for run_rnaseq.sh
#
# Finds every R1/R2 paired-end sample in a fastq directory and submits one
# independent SLURM job per sample (they run in parallel, each in its own
# per-sample output subdirectory created by run_rnaseq.sh).
#
# Usage:
#   ./run_rnaseq_batch.sh <fastq_dir> <output_dir> [--rm_trimmed_fastqs] [--rm_sam]
#
# Example:
#   ./run_rnaseq_batch.sh \
#       /lustre/isaac24/proj/UTK0009/Muhammad/RNAseq/temp \
#       /lustre/isaac24/proj/UTK0009/Muhammad/RNAseq/rnaseq_output \
#       --rm_trimmed_fastqs --rm_sam
###############################################################################

set -euo pipefail

FastqDir=$1
OutDir=$2
shift 2
FLAGS="$@"   # forward any cleanup flags to each per-sample job

# Path to the main pipeline script (same directory as this submitter)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN_SCRIPT="${SCRIPT_DIR}/run_hg19_rnaseq_pipeline.sh"

if [[ ! -f "${MAIN_SCRIPT}" ]]; then
    echo "ERROR: cannot find run_hg19_rnaseq_pipeline.sh next to this submitter (${MAIN_SCRIPT})"
    exit 1
fi

submitted=0
for R1 in "${FastqDir}"/*_R1_001.fastq.gz; do
    # guard against the glob matching nothing
    [[ -e "${R1}" ]] || { echo "No *_R1_001.fastq.gz files found in ${FastqDir}"; exit 1; }

    R1name=$(basename "${R1}")
    R2name=${R1name/_R1_001/_R2_001}

    if [[ ! -f "${FastqDir}/${R2name}" ]]; then
        echo "WARNING: no R2 mate for ${R1name} — skipping"
        continue
    fi

    Sample=$(echo "${R1name}" | sed 's/_R1_001\.fastq\.gz//')
    echo "Submitting: ${Sample}"

    sbatch -J "rna_${Sample}" "${MAIN_SCRIPT}" \
        "${FastqDir}" "${OutDir}" "${R1name}" "${R2name}" ${FLAGS}

    submitted=$((submitted + 1))
done

echo "------------------------------------------------------------"
echo "Submitted ${submitted} sample(s)."
echo "Track with:  squeue -u \$USER"
echo "------------------------------------------------------------"