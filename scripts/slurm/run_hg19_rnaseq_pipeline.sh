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
# Usage:
#   sbatch run_rnaseq.sh <fastq_dir> <output_dir> <R1.fastq.gz> <R2.fastq.gz> \
#          [--rm_trimmed_fastqs] [--rm_sam]
#
# Stages: adapter trim -> quality trim -> STAR align -> sort/index
#         -> htseq count -> CPM bigwig
#
# Output layout (each sample gets its own clean directory):
#   <output_dir>/<sample>/
#       <sample>_genecounts.tsv     raw counts (DESeq2 input; NOT normalized)
#       <sample>_sorted.bam (+.bai) coordinate-sorted alignments
#       <sample>.cpm.norm.bw        CPM-normalized track for IGV
#       <sample>_Log.final.out      STAR mapping summary (+ other STAR logs)
#       intermediates/              trimmed fastqs + SAM (deletable via flags)
#
# Optional flags:
#   --rm_trimmed_fastqs   delete trimmed fastqs at the end
#   --rm_sam              delete the (large) SAM at the end
#   (default: keep everything; intermediates/ is safe to rm -rf manually)
###############################################################################

# ---------------------------------------------------------------------------
# Positional inputs
# ---------------------------------------------------------------------------
FastqDir=$1
OutDir=$2
Fastq1=$3
Fastq2=$4
shift 4  

# ---------------------------------------------------------------------------
# Optional flags (default: keep all intermediates)
# ---------------------------------------------------------------------------
RM_TRIMMED=false
RM_SAM=false
for arg in "$@"; do
    case "$arg" in
        --rm_trimmed_fastqs) RM_TRIMMED=true ;;
        --rm_sam)            RM_SAM=true ;;
        *) echo "ERROR: unknown flag '$arg'"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# References and tools
# ---------------------------------------------------------------------------
GENOME_DIR=/lustre/isaac24/proj/UTK0009/RNA-Seq/genome/hg19
GTF=/lustre/isaac24/proj/UTK0009/RNA-Seq/reference_genomes/hg19/HomoSapiens_hg19.gtf
ADAPTERS=/lustre/isaac24/proj/UTK0009/RNA-Seq/bbmap/resources/adapters.fa
BBDUK=/lustre/isaac24/proj/UTK0009/RNA-Seq/conda_env/opt/bbmap-38.18/bbduk.sh

# Conda env holding deepTools (bamCoverage).
CONDA_ENV=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration/conda-envs/multiomics-py

# Threads = the CPUs SLURM allocated to this single task (fallback 20)
THREADS=${SLURM_CPUS_PER_TASK:-20}

# ---------------------------------------------------------------------------
# Derive sample name and set up a per-sample output directory
# ---------------------------------------------------------------------------
Sample=$(basename "${Fastq1}" | sed 's/_R1_001\.fastq\.gz//')
SampleDir="${OutDir}/${Sample}"
InterDir="${SampleDir}/intermediates"
mkdir -p "${SampleDir}" "${InterDir}"
cd "${SampleDir}"   

echo "============================================================"
echo "Sample : ${Sample}"
echo "Started: $(date)"
echo "Node   : $(hostname)"
echo "Threads: ${THREADS}   Mem: ${SLURM_MEM_PER_NODE:-NA}"
echo "Output : ${SampleDir}"
echo "Flags  : rm_trimmed=${RM_TRIMMED}  rm_sam=${RM_SAM}"
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
    in="${FastqDir}/${Fastq1}" in2="${FastqDir}/${Fastq2}" \
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