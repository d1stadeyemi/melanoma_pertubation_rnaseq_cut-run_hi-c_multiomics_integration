#!/bin/bash
#SBATCH -J deseq2_dge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=00-23:00:00
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --error=/lustre/isaac24/proj/UTK0009/RNA-Seq/logs/errors/%x-%j.err
#SBATCH --output=/lustre/isaac24/proj/UTK0009/RNA-Seq/logs/jobs/%x-%j.out

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
# The key point: --counts-dir is simply the OUTPUT DIRECTORY of your
# rnaseq_fastq_to_counts.sh run.
#
# USAGE (keyword-scanning mode — group membership auto-scanned from --counts-dir):
#   sbatch submit_rnaseq_deseq2_dge.sh \
#       --ref grch38 \
#       --counts-dir /lustre/.../rnaseq_output \
#       --group1 Bottom10 --group2 Parental \
#       --control Parental \
#       --outdir results/exploratory/deseq2_B10_vs_Parental
#   # optional replicate-count assertion + clean labels:
#   #   --group1-reps rep1,rep2,rep3 --group2-reps rep1,rep2,rep3
#
# NOTE on --control: it must be one of your two --group labels (it is the LFC
# denominator / reference level), NOT some other word. e.g. if the groups are
# Bottom10 and Control, then --control Control (or --control Bottom10).
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

# ---------------------------------------------------------------------------
# Locate the R script (it sits next to this wrapper).
# IMPORTANT: under `sbatch`, SLURM executes a SPOOLED COPY of this script, so
# ${BASH_SOURCE[0]} / $0 point at the spool dir (or are unset under `set -u`),
# NOT the real location. So: when running as a SLURM job, ask scontrol for the
# real submitted path; otherwise use BASH_SOURCE/$0. An env override
# (DESEQ2_DGE_RSCRIPT=/path/to/rnaseq_deseq2_dge.R) always wins.
# ---------------------------------------------------------------------------
SELF=""
if [[ -n "${SLURM_JOB_ID:-}" ]] && command -v scontrol >/dev/null 2>&1; then
    SELF="$(scontrol show job "${SLURM_JOB_ID}" 2>/dev/null \
            | awk -F'Command=' 'NF>1{print $2; exit}' | awk '{print $1}')"
fi
if [[ -z "${SELF:-}" || ! -f "${SELF:-}" ]]; then
    SELF="${BASH_SOURCE[0]:-$0}"
fi
SCRIPT_DIR="$(cd "$(dirname "${SELF}")" >/dev/null 2>&1 && pwd)"
R_SCRIPT="${DESEQ2_DGE_RSCRIPT:-${SCRIPT_DIR}/rnaseq_deseq2_dge.R}"

if [[ $# -eq 0 || " $* " == *" --help "* || " $* " == *" -h "* ]]; then
    echo "Wrapper for rnaseq_deseq2_dge.R — all flags pass straight through."
    echo "Showing the R script's help:"; echo
    if [[ -f "${R_SCRIPT}" ]]; then Rscript "${R_SCRIPT}" --help 2>/dev/null || sed -n '1,60p' "${R_SCRIPT}"; fi
    exit 0
fi
if [[ ! -f "${R_SCRIPT}" ]]; then
    echo "ERROR: cannot find rnaseq_deseq2_dge.R."
    echo "       looked at: ${R_SCRIPT}"
    echo "       resolved wrapper dir: ${SCRIPT_DIR}"
    echo "       Fix: keep rnaseq_deseq2_dge.R next to this wrapper, or set"
    echo "            DESEQ2_DGE_RSCRIPT=/full/path/to/rnaseq_deseq2_dge.R"
    exit 1
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