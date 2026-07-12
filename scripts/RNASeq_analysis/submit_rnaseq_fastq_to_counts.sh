#!/bin/bash
###############################################################################
# submit_rnaseq_fastq_to_counts.sh   (batch submitter)
#
# Finds every R1/R2 paired-end sample in a fastq folder and submits one
# independent SLURM job per sample (they run in parallel, each in its own
# per-replicate output subdirectory created by rnaseq_fastq_to_counts.sh).
#
# It also writes  <outdir>/samples_manifest.tsv  mapping every replicate to its
# biological sample and R1/R2 files — handy for building the DESeq2 sample sheet.
#
# File-name flexibility (for fastqs from any sequencing vendor):
#   --fwd-keyword / --rev-keyword   tokens marking forward / reverse reads
#   --rep-keyword                   token marking the replicate number
# Defaults match Illumina/standard naming: R1_001 / R2_001 / rep.
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
  ./submit_rnaseq_fastq_to_counts.sh --ref <hg19|grch38> \
         --fastq-dir <FOLDER> \
         [--outdir <FOLDER>] \
         [--fwd-keyword R1_001] [--rev-keyword R2_001] [--rep-keyword rep] \
         [--rm-trimmed-fastqs] [--rm-sam]

REQUIRED:
  --ref <hg19|grch38>   Reference genome. Must be exactly hg19 OR grch38.
  --fastq-dir <FOLDER>  Folder containing your fastq.gz files. Every forward
                        read found (name contains the --fwd-keyword) is paired
                        with its reverse mate and submitted as its own job.

OPTIONAL:
  --outdir <FOLDER>     Where to write results (absolute path is resolved so all
                        jobs share it). Default: the folder you are in right now.

  --fwd-keyword <STR>   Token marking the FORWARD read.  Default: R1_001
  --rev-keyword <STR>   Token marking the REVERSE read.  Default: R2_001
                        (the reverse file name = forward name with fwd->rev).
  --rep-keyword <STR>   Token marking the REPLICATE number, used only to derive
                        the biological-sample name in samples_manifest.tsv.
                        Default: rep   (matches rep1/rep2/... case-insensitively)
                        Use distinctive keywords so they are not confused with
                        text elsewhere in the file name.

  --rm-trimmed-fastqs   Forwarded to each job: delete trimmed fastqs when done.
  --rm-sam              Forwarded to each job: delete the SAM when done.
  -h, --help            Show this help and exit.

EXAMPLES:
  # Standard Illumina naming (Sample_rep1_R1_001.fastq.gz), grch38
  ./submit_rnaseq_fastq_to_counts.sh --ref grch38 \
        --fastq-dir /lustre/.../fastqs \
        --outdir /lustre/.../rnaseq_output --rm-trimmed-fastqs --rm-sam

  # A different vendor naming, e.g. Ctrl_Rep1_R1.fastq.gz
  ./submit_rnaseq_fastq_to_counts.sh --ref hg19 \
        --fastq-dir /lustre/.../fastqs \
        --fwd-keyword R1 --rev-keyword R2 --rep-keyword Rep

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
FWD_KEYWORD="R1_001"
REV_KEYWORD="R2_001"
REP_KEYWORD="rep"
RM_TRIMMED=false
RM_SAM=false

need_val() { [[ $# -ge 2 && "${2:0:1}" != "-" ]] || { echo "ERROR: option '$1' needs a value."; echo; usage; exit 1; }; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref)                    need_val "$@"; REF="$2";         shift 2 ;;
        --fastq-dir|--fastqdir)   need_val "$@"; FASTQ_DIR="$2";   shift 2 ;;
        --outdir|--output-dir)    need_val "$@"; OUTDIR="$2";      shift 2 ;;
        --fwd-keyword|--fwd)      need_val "$@"; FWD_KEYWORD="$2"; shift 2 ;;
        --rev-keyword|--rev)      need_val "$@"; REV_KEYWORD="$2"; shift 2 ;;
        --rep-keyword|--rep)      need_val "$@"; REP_KEYWORD="$2"; shift 2 ;;
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

# ---------------------------------------------------------------------------
# Name-parsing helpers (keyword-driven; tested for Illumina + vendor naming)
#   rep_name  : R1 file  -> replicate/output-dir name (drop fwd keyword + ext)
#   mate_name : R1 file  -> R2 file    (replace last fwd keyword with rev)
#   samp_name : replicate-> biological sample (drop last <sep>rep<digits>)
# ---------------------------------------------------------------------------
rep_name() {
    local n; n=$(echo "$1" | sed -E 's/\.(fastq|fq)\.gz$//; s/\.(fastq|fq)$//')
    n=$(echo "$n" | sed -E "s/(.*)[._-]?${FWD_KEYWORD}.*/\1/" | sed -E 's/[._-]+$//')
    echo "$n"
}
mate_name() { echo "$1" | sed -E "s/(.*)${FWD_KEYWORD}/\1${REV_KEYWORD}/"; }
samp_name() { echo "$1" | sed -E "s/(.*)[._-]${REP_KEYWORD}[0-9]+.*/\1/I" | sed -E 's/[._-]+$//'; }

# Build the list of cleanup flags to forward to each job (may be empty)
EXTRA=()
$RM_TRIMMED && EXTRA+=(--rm-trimmed-fastqs)
$RM_SAM     && EXTRA+=(--rm-sam)

# Start the manifest (overwrites any previous one for this outdir)
MANIFEST="${OUTDIR}/samples_manifest.tsv"
printf "replicate\tsample\tR1\tR2\n" > "${MANIFEST}"

# ---------------------------------------------------------------------------
# Find forward-read files and submit one job per R1/R2 pair
# ---------------------------------------------------------------------------
shopt -s nullglob
fastqs=( "${FASTQ_DIR}"/*"${FWD_KEYWORD}"*.fastq.gz "${FASTQ_DIR}"/*"${FWD_KEYWORD}"*.fq.gz )
shopt -u nullglob
if [[ ${#fastqs[@]} -eq 0 ]]; then
    echo "ERROR: no forward-read fastqs matching *${FWD_KEYWORD}*.fastq.gz (or .fq.gz) in ${FASTQ_DIR}"
    echo "       if your files use a different token, set --fwd-keyword accordingly."
    exit 1
fi

submitted=0
for R1path in "${fastqs[@]}"; do
    R1name=$(basename "${R1path}")
    R2name=$(mate_name "${R1name}")

    if [[ "${R2name}" == "${R1name}" ]]; then
        echo "WARNING: could not derive R2 from ${R1name} (fwd keyword '${FWD_KEYWORD}' not found) — skipping"
        continue
    fi
    if [[ ! -f "${FASTQ_DIR}/${R2name}" ]]; then
        echo "WARNING: no R2 mate for ${R1name} (expected ${R2name}) — skipping"
        continue
    fi

    Rep=$(rep_name "${R1name}")
    Sample=$(samp_name "${Rep}")
    [[ -n "${Sample}" ]] || Sample="${Rep}"

    printf "%s\t%s\t%s\t%s\n" "${Rep}" "${Sample}" "${R1name}" "${R2name}" >> "${MANIFEST}"
    echo "Submitting: ${Rep}  (sample=${Sample}, ref=${REF_LC})"

    sbatch -J "rna_${Rep}" "${MAIN_SCRIPT}" \
        --ref "${REF_LC}" \
        --fastq-dir "${FASTQ_DIR}" \
        --outdir "${OUTDIR}" \
        --fwd-keyword "${FWD_KEYWORD}" \
        --r1 "${R1name}" --r2 "${R2name}" \
        ${EXTRA[@]+"${EXTRA[@]}"}

    submitted=$((submitted + 1))
done

echo "------------------------------------------------------------"
echo "Submitted ${submitted} sample(s)."
echo "Reference : ${REF_LC}"
echo "Output    : ${OUTDIR}"
echo "Manifest  : ${MANIFEST}"
echo "Track with:  squeue -u \$USER"
echo "------------------------------------------------------------"