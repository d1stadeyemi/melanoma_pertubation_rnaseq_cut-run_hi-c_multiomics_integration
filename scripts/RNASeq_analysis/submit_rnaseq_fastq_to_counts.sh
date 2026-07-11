#!/bin/bash
###############################################################################
# Batch submitter for rnaseq_fastq_to_counts.sh
#
# This ONE script replaces the old hg19-only and grch38-only submitters.
# It finds every R1/R2 paired-end sample in a fastq folder and submits one
# independent SLURM job per sample (they run in parallel, each in its own
# per-sample output subdirectory created by rnaseq_fastq_to_counts.sh).
###############################################################################

set -euo pipefail

# Folder you are running this submitter from (used if --outdir is not given).
SUBMIT_DIR="$PWD"

# ---------------------------------------------------------------------------
# Help text
# ---------------------------------------------------------------------------
usage() {
cat <<'EOF'
------------------------------------------------------------------------------
RNA-seq BATCH submitter  —  submit one job per sample in a folder
------------------------------------------------------------------------------

USAGE:
  ./rnaseq_fastq_to_counts.sh --ref <hg19|grch38> \
         --fastq-dir <FOLDER> \
         [--outdir <FOLDER>] [--rm-trimmed-fastqs] [--rm-sam]

REQUIRED:
  --ref <hg19|grch38>   Reference genome. Must be exactly hg19 OR grch38.
  --fastq-dir <FOLDER>  Folder containing your *_R1_001.fastq.gz /
                        *_R2_001.fastq.gz files. Every R1/R2 pair found here
                        is submitted as its own job.

OPTIONAL:
  --outdir <FOLDER>     Where to write results.
                        If you leave this out, results are written into the
                        folder you are in right now.
  --rm-trimmed-fastqs   Delete trimmed fastqs when finished (saves disk space).
  --rm-sam              Delete the large SAM file when finished (saves disk space).
  -h, --help            Show this help and exit.

EXAMPLES:
  # hg19; results land in the current folder
  ./rnaseq_fastq_to_counts.sh --ref hg19 --fastq-dir /lustre/.../fastqs

  # grch38; results in a chosen folder; clean up intermediates afterward
  ./rnaseq_fastq_to_counts.sh --ref grch38 \
        --fastq-dir /lustre/.../fastqs \
        --outdir /lustre/.../rnaseq_output \
        --rm-trimmed-fastqs --rm-sam

Track your jobs with:  squeue -u $USER
------------------------------------------------------------------------------
EOF
}

# ---------------------------------------------------------------------------
# Parse flags (order does not matter; every value is clearly labeled)
# ---------------------------------------------------------------------------
REF=""
FASTQ_DIR=""
OUTDIR=""
RM_TRIMMED=false
RM_SAM=false

need_val() { [[ $# -ge 2 && "${2:0:1}" != "-" ]] || { echo "ERROR: option '$1' needs a value."; echo; usage; exit 1; }; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref)                    need_val "$@"; REF="$2";       shift 2 ;;
        --fastq-dir|--fastqdir)   need_val "$@"; FASTQ_DIR="$2"; shift 2 ;;
        --outdir|--output-dir)    need_val "$@"; OUTDIR="$2";    shift 2 ;;
        --rm-trimmed-fastqs|--rm_trimmed_fastqs) RM_TRIMMED=true; shift ;;
        --rm-sam|--rm_sam)                       RM_SAM=true;     shift ;;
        -h|--help)                usage; exit 0 ;;
        *) echo "ERROR: unknown option '$1'"; echo; usage; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Validate required flags — refuse to run (and explain) if anything is missing
# ---------------------------------------------------------------------------
missing=()
[[ -n "$REF" ]]       || missing+=("--ref")
[[ -n "$FASTQ_DIR" ]] || missing+=("--fastq-dir")
if [[ ${#missing[@]} -gt 0 ]]; then
    echo "ERROR: you did not provide the required option(s): ${missing[*]}"
    echo
    usage
    exit 1
fi

# Validate --ref value early (so the whole batch fails fast, not each job)
REF_LC=$(echo "$REF" | tr '[:upper:]' '[:lower:]')
case "$REF_LC" in
    hg19|grch38) : ;;   # ok
    *)
        echo "ERROR: --ref must be exactly 'hg19' or 'grch38' (you gave '$REF')."
        echo
        usage
        exit 1
        ;;
esac

# ---------------------------------------------------------------------------
# Output directory: use --outdir if given, otherwise the current folder.
# Resolve to an absolute path so every sample job writes to the same place.
# ---------------------------------------------------------------------------
if [[ -z "$OUTDIR" ]]; then
    OUTDIR="$SUBMIT_DIR"
    echo "NOTE: no --outdir given -> results will go into the current folder:"
    echo "      $OUTDIR"
fi
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

[[ -d "$FASTQ_DIR" ]] || { echo "ERROR: --fastq-dir folder not found: $FASTQ_DIR"; exit 1; }

# ---------------------------------------------------------------------------
# Locate the main pipeline script (must sit next to this submitter)
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN_SCRIPT="${SCRIPT_DIR}/rnaseq_fastq_to_counts.sh"

if [[ ! -f "${MAIN_SCRIPT}" ]]; then
    echo "ERROR: cannot find rnaseq_fastq_to_counts.sh next to this submitter."
    echo "       expected here: ${MAIN_SCRIPT}"
    exit 1
fi

# Build the list of cleanup flags to forward to each job (may be empty)
EXTRA=()
$RM_TRIMMED && EXTRA+=(--rm-trimmed-fastqs)
$RM_SAM     && EXTRA+=(--rm-sam)

# ---------------------------------------------------------------------------
# Find R1/R2 pairs and submit one job each
# ---------------------------------------------------------------------------
submitted=0
for R1 in "${FASTQ_DIR}"/*_R1_001.fastq.gz; do
    # guard against the glob matching nothing
    [[ -e "${R1}" ]] || { echo "No *_R1_001.fastq.gz files found in ${FASTQ_DIR}"; exit 1; }

    R1name=$(basename "${R1}")
    R2name=${R1name/_R1_001/_R2_001}

    if [[ ! -f "${FASTQ_DIR}/${R2name}" ]]; then
        echo "WARNING: no R2 mate for ${R1name} — skipping"
        continue
    fi

    Sample=$(echo "${R1name}" | sed 's/_R1_001\.fastq\.gz//')
    echo "Submitting: ${Sample}  (ref=${REF_LC})"

    sbatch -J "rna_${Sample}" "${MAIN_SCRIPT}" \
        --ref "${REF_LC}" \
        --fastq-dir "${FASTQ_DIR}" \
        --outdir "${OUTDIR}" \
        --r1 "${R1name}" --r2 "${R2name}" \
        ${EXTRA[@]+"${EXTRA[@]}"}

    submitted=$((submitted + 1))
done

echo "------------------------------------------------------------"
echo "Submitted ${submitted} sample(s)."
echo "Reference : ${REF_LC}"
echo "Output    : ${OUTDIR}"
echo "Track with:  squeue -u \$USER"
echo "------------------------------------------------------------"