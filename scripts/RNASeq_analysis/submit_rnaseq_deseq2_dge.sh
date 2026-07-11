#!/bin/bash
#SBATCH -J deseq2_dge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00-02:00:00
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --error=deseq2_dge-%j-err.txt
#SBATCH --output=deseq2_dge-%j-out.txt

# Fail fast: exit on error, undefined variable, or any failed command in a pipe
set -euo pipefail

###############################################################################
# submit_rnaseq_deseq2_dge.sh
# -----------------------------------------------------------------------------
# Thin SLURM + conda wrapper around rnaseq_deseq2_dge.R.
#
# This is NOT a batch fan-out (differential expression is per-CONTRAST, not
# per-sample, so there is nothing to loop over). Its ONLY jobs are:
#   1. run on a compute node (the #SBATCH header above), and
#   2. activate the multiomics-r conda env with the project's env hygiene,
# then hand EVERY flag straight through to the R script unchanged.
#
# You therefore use exactly the R script's flags (run --help to see them all).
# The key point you asked about: --counts-dir is simply the OUTPUT DIRECTORY of
# your rnaseq_fastq_to_counts.sh run.
#
# USAGE (inline two-group mode):
#   sbatch submit_rnaseq_deseq2_dge.sh \
#       --ref grch38 \
#       --counts-dir /lustre/.../rnaseq_output \
#       --group1 Bottom10 --group1-samples B10_r1,B10_r2,B10_r3 \
#       --group2 Parental --group2-samples Par_r1,Par_r2,Par_r3 \
#       --control Parental \
#       --outdir results/exploratory/deseq2_B10_vs_Parental
#
# USAGE (sample-sheet mode, e.g. with a batch covariate):
#   sbatch submit_rnaseq_deseq2_dge.sh \
#       --ref grch38 --counts-dir /lustre/.../rnaseq_output \
#       --samplesheet samples.csv --control Parental --batch batch \
#       --outdir results/exploratory/deseq2_B10_vs_Parental
#
# Prefer to run interactively instead? You don't need this wrapper at all:
#   salloc -A ACF-UTK0009 -p condo-rmccord -q condo -c 4 --mem 32G -t 2:00:00
#   conda activate multiomics-r
#   Rscript rnaseq_deseq2_dge.R --ref grch38 --counts-dir ... [flags]
###############################################################################

# Locate the R script (must sit next to this wrapper) ------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT="${SCRIPT_DIR}/rnaseq_deseq2_dge.R"

if [[ $# -eq 0 || " $* " == *" --help "* || " $* " == *" -h "* ]]; then
    echo "Wrapper for rnaseq_deseq2_dge.R — all flags pass straight through."
    echo "Showing the R script's help:"; echo
    if [[ -f "${R_SCRIPT}" ]]; then Rscript "${R_SCRIPT}" --help 2>/dev/null || cat "${R_SCRIPT}" | sed -n '1,60p'; fi
    exit 0
fi
if [[ ! -f "${R_SCRIPT}" ]]; then
    echo "ERROR: cannot find rnaseq_deseq2_dge.R next to this wrapper (${R_SCRIPT})"; exit 1
fi

# --- Conda env holding R 4.4.3 + DESeq2 (project 'multiomics-r' env) ---------
CONDA_ENV=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration/conda-envs/multiomics-r

echo "============================================================"
echo "submit_rnaseq_deseq2_dge.sh"
echo "Started : $(date)"
echo "Node    : $(hostname)"
echo "CPUs    : ${SLURM_CPUS_PER_TASK:-NA}   Mem: ${SLURM_MEM_PER_NODE:-NA}"
echo "R script: ${R_SCRIPT}"
echo "Args    : $*"
echo "============================================================"

# --- Environment hygiene (matches the rest of the project) ------------------
module purge
module load anaconda3/2024.06
source "$ANACONDA3_SH"
unset PYTHONPATH                    # avoid stray PYTHONPATH leaking in
export PYTHONNOUSERSITE=1           # prevent ~/.local shadowing
conda activate "${CONDA_ENV}"

# Redirect caches to Lustre (never $HOME) ------------------------------------
export XDG_CACHE_HOME="${SLURM_SUBMIT_DIR:-$PWD}/.cache"
mkdir -p "${XDG_CACHE_HOME}"

# --- Run: forward every argument unchanged ----------------------------------
# --threads defaults inside the R script to $SLURM_CPUS_PER_TASK, so you don't
# need to pass it explicitly.
Rscript "${R_SCRIPT}" "$@"

echo "============================================================"
echo "Finished: $(date)"
echo "============================================================"