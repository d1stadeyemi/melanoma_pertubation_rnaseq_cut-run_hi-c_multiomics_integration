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

set -euo pipefail

###############################################################################
# submit_rnaseq_deseq2_dge.sh   (v3)  --  runs the 4-stage DESeq2 pipeline
#
# The analysis is split into four R scripts. This wrapper runs them IN ORDER,
# in ONE job, forwarding your flags unchanged to each. All four APPEND to a
# single <outdir>/run.log, tagged [1/setup] [2/counts] [3/model] [4/results].
#
#   deseq2_01_setup.R    samples, groups, design, batch, full-rank check
#   deseq2_02_counts.R   count matrix, annotation, library-composition QC
#   deseq2_03_model.R    prefilter, DESeq() fit, transforms, sample QC
#   deseq2_04_results.R  contrast, shrinkage, thresholds, DEG tables, figures
#
# Each R script documents its OWN parameters in its header, so to judge whether
# a threshold is right you read one header, not 1200 lines. Nothing about the
# split changes the statistics: stage N reads stage N-1's .rds and refuses to
# run if that file was produced with different arguments.
#
#  ==========================================================================
#  USAGE
#  ==========================================================================
#  # Validate the design without touching count data (ALWAYS DO THIS FIRST):
#    sbatch submit_rnaseq_deseq2_dge.sh --ref grch38 \
#        --counts-dir /lustre/.../rnaseq_output \
#        --samplesheet /lustre/.../rnaseq_output/samplesheet.csv \
#        --control Parental --treatment Bottom10 --batch batch \
#        --outdir results/exploratory/deseq2_B10_vs_Parental --dry-run
#
#  # Full run (defaults: padj<0.05 AND >=1.5-fold):
#    sbatch submit_rnaseq_deseq2_dge.sh --ref grch38 \
#        --counts-dir /lustre/.../rnaseq_output \
#        --samplesheet .../samplesheet.csv \
#        --control Parental --treatment Bottom10 --batch batch \
#        --outdir results/exploratory/deseq2_B10_vs_Parental
#
#  # Keyword scan instead of a sample sheet, batches with no file at all:
#    sbatch submit_rnaseq_deseq2_dge.sh --ref grch38 --counts-dir DIR \
#        --group1 Bottom10 --group2 Parental --control Parental \
#        --batch-map "rep1=A,rep2=B,rep3=C" --outdir OUT
#
#  # Change only the thresholds -- re-runs stage 4 in seconds, not the fit:
#    sbatch submit_rnaseq_deseq2_dge.sh <same flags> --stage 4 \
#        --alpha 0.01 --lfc-report 1
#
#  # Strict threshold TEST instead of post-hoc filtering:
#    ... --lfc-threshold 0.585        (p-values then test |LFC| > 0.585)
#
#  # Equivalence test ("which genes did NOT change?"):
#    ... --alt-hypothesis lessAbs --lfc-threshold 0.585
#
#  # Three conditions:
#    ... --control DMSO --treatment UNC1999          # UNC vs DMSO, all groups fitted
#    ... --control DMSO --treatment UNC1999 --subset-to-contrast   # drop the rest
#    ... --test lrt                                   # does condition matter at all
#
#  # Same library across two lanes (technical, NOT biological, replicates):
#    ... --collapse-technical biological_sample
#
#  ==========================================================================
#  SUBMITTER-ONLY FLAGS
#  ==========================================================================
#  --stage <1|2|3|4|all|N-|N-M>   Which stages to run. Default: all
#                                 4      only stage 4
#                                 3-     stages 3 and 4
#                                 2-3    stages 2 and 3
#  --force                        Reuse upstream state even if it was produced
#                                 with different arguments. Forwarded to R.
#  Everything else is passed to the R scripts unchanged; run --help for the list.
#
#  Prefer interactive? The wrapper is optional:
#    salloc -A ACF-UTK0009 -p condo-rmccord -q condo -c 4 --mem 32G -t 2:00:00
#    conda activate multiomics-r
#    Rscript deseq2_01_setup.R --ref grch38 --counts-dir ... [flags]
#    Rscript deseq2_02_counts.R  <same flags>     # etc.
###############################################################################

# ---------------------------------------------------------------------------
# Locate the R scripts (they sit next to this wrapper).
# Under `sbatch`, SLURM runs a SPOOLED COPY, so ${BASH_SOURCE[0]} points at the
# spool dir. Ask scontrol for the real submitted path when in a job.
# ---------------------------------------------------------------------------
SELF=""
if [[ -n "${SLURM_JOB_ID:-}" ]] && command -v scontrol >/dev/null 2>&1; then
    SELF="$(scontrol show job "${SLURM_JOB_ID}" 2>/dev/null \
            | awk -F'Command=' 'NF>1{print $2; exit}' | awk '{print $1}')"
fi
[[ -z "${SELF:-}" || ! -f "${SELF:-}" ]] && SELF="${BASH_SOURCE[0]:-$0}"
SCRIPT_DIR="${DESEQ2_DGE_DIR:-$(cd "$(dirname "${SELF}")" >/dev/null 2>&1 && pwd)}"

STAGE_FILES=( "deseq2_01_setup.R" "deseq2_02_counts.R" "deseq2_03_model.R" "deseq2_04_results.R" )
STAGE_NAMES=( "1/setup" "2/counts" "3/model" "4/results" )

# ---------------------------------------------------------------------------
# Split --stage off; everything else is forwarded verbatim
# ---------------------------------------------------------------------------
STAGE_SPEC="all"
PASS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --stage)   [[ $# -ge 2 ]] || { echo "ERROR: --stage needs a value."; exit 1; }
                   STAGE_SPEC="$2"; shift 2 ;;
        --stage=*) STAGE_SPEC="${1#--stage=}"; shift ;;
        *)         PASS+=( "$1" ); shift ;;
    esac
done

show_help() {
    echo "Wrapper for the 4-stage DESeq2 pipeline. All flags pass through to the R scripts."
    echo "Submitter-only: --stage <1|2|3|4|all|N-|N-M>"
    echo
    if [[ -f "${SCRIPT_DIR}/${STAGE_FILES[0]}" ]]; then
        Rscript "${SCRIPT_DIR}/${STAGE_FILES[0]}" --help 2>/dev/null \
          || sed -n '1,80p' "${SCRIPT_DIR}/${STAGE_FILES[0]}"
    fi
}
if [[ ${#PASS[@]} -eq 0 ]] || printf '%s\n' "${PASS[@]}" | grep -qxE -- '--help|-h'; then
    show_help; exit 0
fi

for f in "${STAGE_FILES[@]}" "deseq2_common.R"; do
    if [[ ! -f "${SCRIPT_DIR}/${f}" ]]; then
        echo "ERROR: cannot find ${f}"
        echo "       looked in: ${SCRIPT_DIR}"
        echo "       Keep all five R files together, or set DESEQ2_DGE_DIR=/path/to/them"
        exit 1
    fi
done

# --stage -> FIRST..LAST
case "${STAGE_SPEC}" in
    all)              FIRST=1; LAST=4 ;;
    [1-4])            FIRST="${STAGE_SPEC}"; LAST="${STAGE_SPEC}" ;;
    [1-4]-)           FIRST="${STAGE_SPEC:0:1}"; LAST=4 ;;
    [1-4]-[1-4])      FIRST="${STAGE_SPEC:0:1}"; LAST="${STAGE_SPEC:2:1}" ;;
    *) echo "ERROR: --stage must be 1|2|3|4|all|N-|N-M (got '${STAGE_SPEC}')"; exit 1 ;;
esac
(( FIRST <= LAST )) || { echo "ERROR: --stage range is backwards ('${STAGE_SPEC}')"; exit 1; }

# --dry-run stops inside stage 1 anyway; do not run later stages against no state.
if printf '%s\n' "${PASS[@]}" | grep -qx -- '--dry-run'; then
    (( FIRST == 1 )) || { echo "ERROR: --dry-run only makes sense from stage 1."; exit 1; }
    LAST=1
    echo "NOTE: --dry-run set -> running stage 1 only (design validation)."
fi

echo "============================================================"
echo "submit_rnaseq_deseq2_dge.sh (v3, 4 stages)"
echo "Started  : $(date)"
echo "Node     : $(hostname)"
echo "CPUs     : ${SLURM_CPUS_PER_TASK:-NA}   Mem: ${SLURM_MEM_PER_NODE:-NA}"
echo "Scripts  : ${SCRIPT_DIR}"
echo "Stages   : ${FIRST}..${LAST}"
echo "Args     : ${PASS[*]}"
echo "============================================================"

# --- Environment hygiene (matches the rest of the project) ------------------
CONDA_ENV=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration/conda-envs/multiomics-r
if [[ -z "${DESEQ2_SKIP_ENV:-}" ]]; then
    module purge
    module load anaconda3/2024.06
    unset PYTHONPATH                  # BEFORE activation
    export PYTHONNOUSERSITE=1
    # shellcheck disable=SC1090  # path provided by the module system at runtime
    source "$ANACONDA3_SH"
    conda activate "${CONDA_ENV}"
fi
export XDG_CACHE_HOME="${SLURM_SUBMIT_DIR:-$PWD}/.cache"
mkdir -p "${XDG_CACHE_HOME}"

# --- Run the stages in order ------------------------------------------------
# Stage 1 truncates run.log; stages 2-4 append. A failure stops the chain and
# names the stage, so you know whether the problem is the sample sheet, the
# annotation, the fit, or the contrast.
START_ALL=$SECONDS
for (( s=FIRST; s<=LAST; s++ )); do
    idx=$(( s - 1 ))
    echo
    echo "------------------------------------------------------------"
    echo ">> STAGE ${STAGE_NAMES[$idx]}  (${STAGE_FILES[$idx]})"
    echo "------------------------------------------------------------"
    t0=$SECONDS
    rc=0
    Rscript "${SCRIPT_DIR}/${STAGE_FILES[$idx]}" ${PASS[@]+"${PASS[@]}"} || rc=$?
    if (( rc != 0 )); then
        echo
        echo "============================================================"
        echo "FAILED at stage ${STAGE_NAMES[$idx]} (exit ${rc}) after $(( SECONDS - t0 ))s"
        echo "  The error above names the cause. See <outdir>/run.log for the"
        echo "  full record; earlier stages' output is intact, so once fixed you"
        echo "  can resume with:  --stage ${s}-"
        echo "============================================================"
        exit "${rc}"
    fi
    echo ">> stage ${STAGE_NAMES[$idx]} finished in $(( SECONDS - t0 ))s"
done

echo
echo "============================================================"
echo "All stages complete in $(( SECONDS - START_ALL ))s"
echo "Finished : $(date)"
echo "============================================================"