#!/bin/bash
#SBATCH -J RNA-Seq_pipeline
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=00-23:00:00
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --error=/lustre/isaac24/proj/UTK0009/RNA-Seq/logs/errors/%x-%j-err.txt
#SBATCH --output=/lustre/isaac24/proj/UTK0009/RNA-Seq/logs/jobs/%x-%j-out.txt

# Fail fast: exit on error, undefined variable, or any failed command in a pipe
set -euo pipefail

###############################################################################
# RNA-seq pipeline  (one sample = one paired-end R1/R2 set)
#
# This script replaces the old hg19-only and grch38-only scripts.
# You choose the reference genome with the required --ref flag.
#
# Stages: adapter trim -> quality trim -> STAR align -> sort/index
#         -> htseq count -> CPM bigwig
#
# Output layout (each sample gets its own clean directory):
#   <outdir>/<sample>/
#       <sample>_genecounts.tsv     raw counts (DESeq2 input; NOT normalized)
#       <sample>_sorted.bam (+.bai) coordinate-sorted alignments
#       <sample>.cpm.norm.bw        CPM-normalized track for IGV
#       <sample>_Log.final.out      STAR mapping summary (+ other STAR logs)
#       intermediates/              trimmed fastqs + SAM (deletable via flags)
###############################################################################

# Remember the folder the job was submitted from, BEFORE we cd anywhere.
# Under SLURM this is where you ran `sbatch`; outside SLURM it is your current dir.
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

# ---------------------------------------------------------------------------
# Help text
# ---------------------------------------------------------------------------
usage() {
cat <<'EOF'
------------------------------------------------------------------------------
RNA-seq pipeline  —  run ONE paired-end sample (R1 + R2)
------------------------------------------------------------------------------

USAGE:
  sbatch rnaseq_fastq_to_counts.sh --ref <hg19|grch38> \
         --fastq-dir <FOLDER> \
         --r1 <SAMPLE_R1_001.fastq.gz> --r2 <SAMPLE_R2_001.fastq.gz> \
         [--outdir <FOLDER>] [--rm-trimmed-fastqs] [--rm-sam]

REQUIRED:
  --ref <hg19|grch38>   Reference genome. Must be exactly hg19 OR grch38.
  --fastq-dir <FOLDER>  Folder that contains your R1 and R2 fastq.gz files.
  --r1 <FILE>           R1 file NAME (just the name; it must live in --fastq-dir).
  --r2 <FILE>           R2 file NAME (just the name; it must live in --fastq-dir).

OPTIONAL:
  --outdir <FOLDER>     Where to write results.
                        If you leave this out, results are written into the
                        folder you were in when you submitted the job.
  --rm-trimmed-fastqs   Delete trimmed fastqs when finished (saves disk space).
  --rm-sam              Delete the large SAM file when finished (saves disk space).
  -h, --help            Show this help and exit.

EXAMPLES:
  # hg19; results land in the current folder (no --outdir given)
  sbatch run_rnaseq_pipeline.sh --ref hg19 \
         --fastq-dir /lustre/.../fastqs \
         --r1 Sample1_R1_001.fastq.gz --r2 Sample1_R2_001.fastq.gz

  # grch38; results in a chosen folder; clean up intermediates afterward
  sbatch run_rnaseq_pipeline.sh --ref grch38 \
         --fastq-dir /lustre/.../fastqs \
         --r1 Sample1_R1_001.fastq.gz --r2 Sample1_R2_001.fastq.gz \
         --outdir /lustre/.../rnaseq_output --rm-trimmed-fastqs --rm-sam
------------------------------------------------------------------------------
EOF
}

# ---------------------------------------------------------------------------
# Parse flags (order does not matter; every value is clearly labeled)
# ---------------------------------------------------------------------------
REF=""
FASTQ_DIR=""
R1=""
R2=""
OUTDIR=""
RM_TRIMMED=false
RM_SAM=false

# small helper so a flag given without a value fails with a clear message
need_val() { [[ $# -ge 2 && "${2:0:1}" != "-" ]] || { echo "ERROR: option '$1' needs a value."; echo; usage; exit 1; }; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref)                    need_val "$@"; REF="$2";       shift 2 ;;
        --fastq-dir|--fastqdir)   need_val "$@"; FASTQ_DIR="$2"; shift 2 ;;
        --r1|--R1)                need_val "$@"; R1="$2";        shift 2 ;;
        --r2|--R2)                need_val "$@"; R2="$2";        shift 2 ;;
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
[[ -n "$REF" ]]        || missing+=("--ref")
[[ -n "$FASTQ_DIR" ]]  || missing+=("--fastq-dir")
[[ -n "$R1" ]]         || missing+=("--r1")
[[ -n "$R2" ]]         || missing+=("--r2")
if [[ ${#missing[@]} -gt 0 ]]; then
    echo "ERROR: you did not provide the required option(s): ${missing[*]}"
    echo
    usage
    exit 1
fi

# ---------------------------------------------------------------------------
# Pick the reference genome from --ref (this is the ONLY thing that used to
# differ between the two old scripts)
# ---------------------------------------------------------------------------
REF_LC=$(echo "$REF" | tr '[:upper:]' '[:lower:]')
case "$REF_LC" in
    hg19)
        GENOME_DIR=/lustre/isaac24/proj/UTK0009/RNA-Seq/genome/hg19
        GTF=/lustre/isaac24/proj/UTK0009/RNA-Seq/reference_genomes/hg19/HomoSapiens_hg19.gtf
        ;;
    grch38)
        GENOME_DIR=/lustre/isaac24/proj/UTK0009/RNA-Seq/genome/grch38
        GTF=/lustre/isaac24/proj/UTK0009/RNA-Seq/reference_genomes/grch38/gencode.v49.annotation.gtf
        ;;
    *)
        echo "ERROR: --ref must be exactly 'hg19' or 'grch38' (you gave '$REF')."
        echo
        usage
        exit 1
        ;;
esac

# ---------------------------------------------------------------------------
# Output directory: use --outdir if given, otherwise the submission folder
# ---------------------------------------------------------------------------
if [[ -z "$OUTDIR" ]]; then
    OUTDIR="$SUBMIT_DIR"
    echo "NOTE: no --outdir given -> writing results into the current folder:"
    echo "      $OUTDIR"
fi
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"   # make it an absolute path

# ---------------------------------------------------------------------------
# Check the fastq inputs actually exist before doing any heavy work
# ---------------------------------------------------------------------------
[[ -d "$FASTQ_DIR" ]]        || { echo "ERROR: --fastq-dir folder not found: $FASTQ_DIR"; exit 1; }
[[ -f "$FASTQ_DIR/$R1" ]]    || { echo "ERROR: R1 file not found: $FASTQ_DIR/$R1"; exit 1; }
[[ -f "$FASTQ_DIR/$R2" ]]    || { echo "ERROR: R2 file not found: $FASTQ_DIR/$R2"; exit 1; }

# ---------------------------------------------------------------------------
# Tools / conda env (same for both references)
# ---------------------------------------------------------------------------
ADAPTERS=/lustre/isaac24/proj/UTK0009/RNA-Seq/bbmap/resources/adapters.fa
BBDUK=/lustre/isaac24/proj/UTK0009/RNA-Seq/conda_env/opt/bbmap-38.18/bbduk.sh
# Conda env holding deepTools (bamCoverage).
CONDA_ENV=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration/conda-envs/multiomics-py

# Threads = the CPUs SLURM allocated to this single task (fallback 20)
THREADS=${SLURM_CPUS_PER_TASK:-20}

# ---------------------------------------------------------------------------
# Derive sample name and set up a per-sample output directory
# ---------------------------------------------------------------------------
Sample=$(basename "${R1}" | sed 's/_R1_001\.fastq\.gz//')
SampleDir="${OUTDIR}/${Sample}"
InterDir="${SampleDir}/intermediates"
mkdir -p "${SampleDir}" "${InterDir}"
cd "${SampleDir}"

echo "============================================================"
echo "Sample    : ${Sample}"
echo "Reference : ${REF_LC}"
echo "Started   : $(date)"
echo "Node      : $(hostname)"
echo "Threads   : ${THREADS}   Mem: ${SLURM_MEM_PER_NODE:-NA}"
echo "Fastq dir : ${FASTQ_DIR}"
echo "Output    : ${SampleDir}"
echo "Genome    : ${GENOME_DIR}"
echo "GTF       : ${GTF}"
echo "Flags     : rm_trimmed=${RM_TRIMMED}  rm_sam=${RM_SAM}"
echo "============================================================"

# ---------------------------------------------------------------------------
# Filenames (intermediates in subdir; final outputs in the sample dir)
# ---------------------------------------------------------------------------
AdaptTrim1=${InterDir}/${Sample}_R1_adapter.fastq.gz
AdaptTrim2=${InterDir}/${Sample}_R2_adapter.fastq.gz
QualTrim1=${InterDir}/${Sample}_R1_qc.fastq.gz
QualTrim2=${InterDir}/${Sample}_R2_qc.fastq.gz
StarPrefix=${Sample}_
SamFile=${InterDir}/${Sample}_Aligned.out.sam
SortedBam=${Sample}_sorted.bam
Counts=${Sample}_genecounts.tsv
BigWig=${Sample}.cpm.norm.bw

# ---------------------------------------------------------------------------
# 1. Adapter trimming
# ---------------------------------------------------------------------------
echo ">> [1/6] adapter trimming"
${BBDUK} -Xmx8g threads=${THREADS} \
    in="${FASTQ_DIR}/${R1}" in2="${FASTQ_DIR}/${R2}" \
    ref="${ADAPTERS}" out1="${AdaptTrim1}" out2="${AdaptTrim2}" \
    ktrim=r k=23 mink=11 hdist=1

# ---------------------------------------------------------------------------
# 2. Quality trimming (3' end, Q28)
# ---------------------------------------------------------------------------
echo ">> [2/6] quality trimming"
${BBDUK} -Xmx8g threads=${THREADS} \
    in="${AdaptTrim1}" in2="${AdaptTrim2}" \
    out1="${QualTrim1}" out2="${QualTrim2}" \
    qtrim=r trimq=28

# ---------------------------------------------------------------------------
# 3. STAR alignment
# ---------------------------------------------------------------------------
echo ">> [3/6] STAR alignment"
module load star/2.7.6a
STAR --runThreadN ${THREADS} \
     --genomeDir "${GENOME_DIR}" \
     --readFilesIn "${QualTrim1}" "${QualTrim2}" \
     --readFilesCommand zcat \
     --outFilterScoreMinOverLread 0.2 \
     --outFilterMatchNminOverLread 0.2 \
     --outFileNamePrefix "${StarPrefix}"
# STAR writes <prefix>Aligned.out.sam into the cwd; move it into intermediates/
mv "${StarPrefix}Aligned.out.sam" "${SamFile}"

# ---------------------------------------------------------------------------
# 4. Sort + index (needed for htseq pos-mode and for bigwig)
# ---------------------------------------------------------------------------
echo ">> [4/6] sort + index"
module load samtools/1.16.1-gcc
samtools sort -@ ${THREADS} -o "${SortedBam}" "${SamFile}"
samtools index "${SortedBam}"

# ---------------------------------------------------------------------------
# 5. htseq gene counts (raw integers — DESeq2 normalizes internally)
# ---------------------------------------------------------------------------
echo ">> [5/6] htseq counts"
module load htseq/0.11.2-gcc
htseq-count --format bam --order pos --mode intersection-strict --stranded no \
    "${SortedBam}" "${GTF}" > "${Counts}"

# ---------------------------------------------------------------------------
# 6. bigwig (CPM-normalized, for IGV)
# ---------------------------------------------------------------------------
echo ">> [6/6] bigwig generation"
module purge
module load anaconda3/2024.06
source $ANACONDA3_SH
conda activate "${CONDA_ENV}"
unset PYTHONPATH
bamCoverage --bam "${SortedBam}" --outFileName "${BigWig}" \
    --binSize 10 --normalizeUsing CPM --numberOfProcessors ${THREADS}

# ---------------------------------------------------------------------------
# 7. Optional cleanup (only when flags are passed)
# ---------------------------------------------------------------------------
if [[ "${RM_TRIMMED}" == true ]]; then
    echo ">> removing trimmed fastqs"
    rm -f "${AdaptTrim1}" "${AdaptTrim2}" "${QualTrim1}" "${QualTrim2}"
fi
if [[ "${RM_SAM}" == true ]]; then
    echo ">> removing SAM"
    rm -f "${SamFile}"
fi
# remove intermediates/ if it ended up empty
rmdir "${InterDir}" 2>/dev/null || true

echo "============================================================"
echo "Finished: ${Sample} at $(date)"
echo "Kept    : ${SortedBam}, ${Counts}, ${BigWig}"
echo "============================================================"