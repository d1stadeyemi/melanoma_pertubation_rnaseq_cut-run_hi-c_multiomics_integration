#!/bin/bash
#SBATCH -J rnaseq
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
# rnaseq_fastq_to_counts.sh   (one sample = one paired-end R1/R2 set)
#
# This script replaces the old hg19-only and grch38-only scripts.
# You choose the reference genome with the required --ref flag.
#
# Stages: adapter trim -> quality trim -> STAR align -> sort/index
#         -> htseq count -> CPM bigwig
#
# Output layout (each fastq PAIR gets its own clean directory):
#   <outdir>/<replicate>/
#       <replicate>_genecounts.tsv     raw counts (htseq-style; DESeq2 input)
#       <replicate>_featurecounts.txt(.summary)  native featureCounts output
#       <replicate>_sorted.bam (+.bai) coordinate-sorted alignments
#       <replicate>.cpm.norm.bw        CPM-normalized track for IGV
#       <replicate>_Log.final.out      STAR mapping summary (+ other STAR logs)
#       trimmed_fastqs/                trimmed fastqs   (removable via --rm-trimmed-fastqs)
#       sam/                           STAR SAM         (removable via --rm-sam)
#
# The <replicate> directory name is derived from the R1 file name by removing
# the forward-read keyword (default R1_001) and the extension. Change the
# keyword with --fwd-keyword if your fastqs are named differently.
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
rnaseq_fastq_to_counts.sh  —  run ONE paired-end sample (R1 + R2)
------------------------------------------------------------------------------

USAGE:
  sbatch rnaseq_fastq_to_counts.sh --ref <hg19|grch38> \
         --fastq-dir <FOLDER> \
         --r1 <SAMPLE_R1_001.fastq.gz> --r2 <SAMPLE_R2_001.fastq.gz> \
         [--outdir <FOLDER>] [--fwd-keyword R1_001] \
         [--rm-trimmed-fastqs] [--rm-sam]

REQUIRED:
  --ref <hg19|grch38>   Reference genome. Must be exactly hg19 OR grch38.
  --fastq-dir <FOLDER>  Folder that contains your R1 and R2 fastq.gz files.
  --r1 <FILE>           R1 file NAME (just the name; it must live in --fastq-dir).
  --r2 <FILE>           R2 file NAME (just the name; it must live in --fastq-dir).

OPTIONAL:
  --outdir <FOLDER>     Where to write results.
                        If you leave this out, results are written into the
                        folder you were in when you submitted the job.
  --fwd-keyword <STR>   Token marking the forward read in the file name; it is
                        stripped to name the output directory. Default: R1_001
                        (e.g. use "R1" if your files are Sample_rep1_R1.fastq.gz).
  --stranded <no|yes|reverse>
                        Library strandedness for featureCounts (-s 0|1|2).
                        Default: no. VERIFY with RSeQC infer_experiment.py before
                        trusting counts; unstranded ("no") is the safe choice when
                        samples are mixed or unknown.
  --rm-trimmed-fastqs   Delete the trimmed_fastqs/ contents when finished.
  --rm-sam              Delete the sam/ contents when finished (saves disk space).
  -h, --help            Show this help and exit.

EXAMPLES:
  # hg19; results land in the current folder (no --outdir given)
  sbatch rnaseq_fastq_to_counts.sh --ref hg19 \
         --fastq-dir /lustre/.../fastqs \
         --r1 Sample1_R1_001.fastq.gz --r2 Sample1_R2_001.fastq.gz

  # grch38; results in a chosen folder; clean up intermediates afterward
  sbatch rnaseq_fastq_to_counts.sh --ref grch38 \
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
FWD_KEYWORD="R1_001"
STRANDED="no"
RM_TRIMMED=false
RM_SAM=false

# small helper so a flag given without a value fails with a clear message
need_val() { [[ $# -ge 2 && "${2:0:1}" != "-" ]] || { echo "ERROR: option '$1' needs a value."; echo; usage; exit 1; }; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref)                    need_val "$@"; REF="$2";         shift 2 ;;
        --fastq-dir|--fastqdir)   need_val "$@"; FASTQ_DIR="$2";   shift 2 ;;
        --r1|--R1)                need_val "$@"; R1="$2";          shift 2 ;;
        --r2|--R2)                need_val "$@"; R2="$2";          shift 2 ;;
        --outdir|--output-dir)    need_val "$@"; OUTDIR="$2";      shift 2 ;;
        --fwd-keyword|--fwd)      need_val "$@"; FWD_KEYWORD="$2"; shift 2 ;;
        --stranded)               need_val "$@"; STRANDED="$2";    shift 2 ;;
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
CONDA_ENV_SUBREAD=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration/conda-envs/subread

# Threads = the CPUs SLURM allocated to this single task (fallback 20)
THREADS=${SLURM_CPUS_PER_TASK:-20}

# ---------------------------------------------------------------------------
# Derive the replicate (directory) name from R1 by removing the forward-read
# keyword + extension, then trimming any trailing separator.
#   Parental_rep1_R1_001.fastq.gz  --(fwd=R1_001)-->  Parental_rep1
# ---------------------------------------------------------------------------
noext=$(echo "${R1}" | sed -E 's/\.(fastq|fq)\.gz$//; s/\.(fastq|fq)$//')
Sample=$(echo "${noext}" | sed -E "s/(.*)[._-]?${FWD_KEYWORD}.*/\1/" | sed -E 's/[._-]+$//')
[[ -n "${Sample}" ]] || Sample="${noext}"   # fallback if keyword stripped everything

SampleDir="${OUTDIR}/${Sample}"
TrimDir="${SampleDir}/trimmed_fastqs"   # trimmed fastqs live here (removable)
SamDir="${SampleDir}/sam"               # STAR SAM lives here    (removable)
mkdir -p "${SampleDir}" "${TrimDir}" "${SamDir}"
cd "${SampleDir}"

echo "============================================================"
echo "Replicate : ${Sample}"
echo "Reference : ${REF_LC}"
echo "Started   : $(date)"
echo "Node      : $(hostname)"
echo "Threads   : ${THREADS}   Mem: ${SLURM_MEM_PER_NODE:-NA}"
echo "Fastq dir : ${FASTQ_DIR}"
echo "Output    : ${SampleDir}"
echo "Genome    : ${GENOME_DIR}"
echo "GTF       : ${GTF}"
echo "Flags     : fwd_keyword=${FWD_KEYWORD}  stranded=${STRANDED}  rm_trimmed=${RM_TRIMMED}  rm_sam=${RM_SAM}"
echo "============================================================"

# ---------------------------------------------------------------------------
# Filenames (trimmed fastqs -> trimmed_fastqs/ ; SAM -> sam/ ; finals in dir)
# ---------------------------------------------------------------------------
AdaptTrim1=${TrimDir}/${Sample}_R1_adapter.fastq.gz
AdaptTrim2=${TrimDir}/${Sample}_R2_adapter.fastq.gz
QualTrim1=${TrimDir}/${Sample}_R1_qc.fastq.gz
QualTrim2=${TrimDir}/${Sample}_R2_qc.fastq.gz
StarPrefix=${Sample}_
SamFile=${SamDir}/${Sample}_Aligned.out.sam
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
# STAR writes <prefix>Aligned.out.sam into the cwd; move it into sam/
mv "${StarPrefix}Aligned.out.sam" "${SamFile}"

# ---------------------------------------------------------------------------
# 4. Sort + index (needed for htseq pos-mode and for bigwig)
# ---------------------------------------------------------------------------
echo ">> [4/6] sort + index"
module load samtools/1.16.1-gcc
samtools sort -@ ${THREADS} -o "${SortedBam}" "${SamFile}"
samtools index "${SortedBam}"

# ---------------------------------------------------------------------------
# 5. featureCounts gene counts (gene meta-features; raw integers)
#    Output is converted to an htseq-style 2-column file (+ __ summary rows) so
#    the downstream DESeq2 script reads it unchanged (Figure 02 QC preserved).
# ---------------------------------------------------------------------------
# map --stranded to featureCounts -s
case "${STRANDED}" in
    no|unstranded)  STRAND_S=0 ;;
    yes|forward)    STRAND_S=1 ;;
    reverse)        STRAND_S=2 ;;
    *) echo "ERROR: --stranded must be no|yes|reverse (got '${STRANDED}')"; exit 1 ;;
esac

echo ">> [5/6] featureCounts (-s ${STRAND_S}, meta-feature gene_id)"
module purge
module load anaconda3/2024.06
source $ANACONDA3_SH
conda activate "${CONDA_ENV_SUBREAD}"
unset PYTHONPATH

FCraw=${Sample}_featurecounts.txt
featureCounts -a "${GTF}" -o "${FCraw}" \
    -p --countReadPairs \
    -t exon -g gene_id \
    -s ${STRAND_S} \
    -T ${THREADS} \
    "${SortedBam}"

# --- make an htseq-compatible counts file: gene<TAB>count, then __ summary rows
#     featureCounts columns: Geneid Chr Start End Strand Length <count>  -> keep 1,7
tail -n +3 "${FCraw}" | cut -f1,7 > "${Counts}"
awk -F'\t' '
  $1=="Unassigned_NoFeatures"     {nf=$2}
  $1=="Unassigned_Ambiguity"      {amb=$2}
  $1=="Unassigned_MappingQuality" {mq=$2}
  $1=="Unassigned_Unmapped"       {un=$2}
  $1=="Unassigned_MultiMapping"   {mm=$2}
  END{
    printf "__no_feature\t%d\n",          nf+0
    printf "__ambiguous\t%d\n",           amb+0
    printf "__too_low_aQual\t%d\n",       mq+0
    printf "__not_aligned\t%d\n",         un+0
    printf "__alignment_not_unique\t%d\n", mm+0
  }' "${FCraw}.summary" >> "${Counts}"

conda deactivate

# ---------------------------------------------------------------------------
# 6. bigwig (CPM-normalized, for IGV)
# ---------------------------------------------------------------------------
echo ">> [6/6] bigwig generation"
module purge
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
# remove trimmed_fastqs/ and sam/ if they ended up empty
rmdir "${TrimDir}" 2>/dev/null || true
rmdir "${SamDir}"  2>/dev/null || true

echo "============================================================"
echo "Finished: ${Sample} at $(date)"
echo "Kept    : ${SortedBam}, ${Counts}, ${BigWig}"
echo "============================================================"