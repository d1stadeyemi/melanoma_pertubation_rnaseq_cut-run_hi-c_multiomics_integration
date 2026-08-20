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
# rnaseq_fastq_to_counts.sh   (v2)   -- ONE sample = ONE paired-end R1/R2 set
#
# Stages: trim (BBDuk, single pass) -> STAR (sorted BAM) -> index
#         -> strandedness detection -> featureCounts -> CPM bigwig -> QC row
#
# Output layout (each fastq PAIR gets its own clean directory):
#   <outdir>/<replicate>/
#       <replicate>_genecounts.tsv        raw counts (htseq-style; DESeq2 input)
#       <replicate>_featurecounts.txt(.summary)   native featureCounts output
#       <replicate>_sorted.bam (+.bai)    coordinate-sorted alignments
#       <replicate>.cpm.norm.bw           CPM-normalized track for IGV
#       <replicate>_library_qc.tsv        ONE-ROW QC RECORD (see below) <== NEW
#       <replicate>_Log.final.out         STAR mapping summary
#       logs/                             BBDuk + tool logs
#       trimmed_fastqs/                   trimmed fastqs (--rm-trimmed-fastqs)
#
###############################################################################
#
#  ==========================================================================
#  PARAMETER TABLE -- every threshold this script uses, and why
#  (read this instead of hunting through the code)
#  ==========================================================================
#
#  --- TRIMMING (BBDuk) -----------------------------------------------------
#  PARAMETER        v1 (OLD)      v2 (NOW)     RATIONALE
#  ktrim            r             r            3' adapter trimming. unchanged.
#  k                23            23           adapter kmer length. unchanged.
#  mink             11            11           shorter kmer at read tip. unchanged.
#  hdist            1             1            1 mismatch allowed in kmer. unchanged.
#  tbo              (absent)      ON           BBDuk author (Bushnell) recommends
#                                              tbo+tpe for "normal paired-end
#                                              fragment libraries". tbo trims
#                                              adapter using pair-overlap, which
#                                              catches adapter that kmer matching
#                                              misses in low-quality reads.
#  tpe              (absent)      ON           trims both mates to equal length so
#                                              one-sided adapter detection does not
#                                              desynchronise the pair.
#  ftm              (absent)      5            drops Illumina's inaccurate extra
#                                              base (151bp -> 150bp). No-op on
#                                              lengths already a multiple of 5.
#  qtrim            r             r            3'-end quality trimming. unchanged.
#  trimq            28            10           *** THE BIG ONE ***
#                                              Q28 is far outside every published
#                                              recommendation. MacManes 2014
#                                              (Front Genet 5:13) found gentle
#                                              trimming at Phred <2 or <5 optimal
#                                              for mRNA-seq; Bushnell calls Q25
#                                              "very stringent" and says Q5-Q10
#                                              suffices for mapping.
#                                              WORSE: on NovaSeq/RTA3 the
#                                              instrument emits only THREE binned
#                                              Q-scores (12, 23, 37). trimq=28
#                                              therefore discards every base that
#                                              is not in the top bin -- it is not
#                                              "strict trimming", it is a
#                                              near-total truncation of Q23 bases.
#                                              Q10 is gentle-to-moderate and
#                                              defensible; use --trimq 5 to follow
#                                              MacManes exactly.
#  minlen           (ABSENT!)     36           *** THE OTHER BIG ONE ***
#                                              v1 set NO minimum length, so BBDuk's
#                                              default minlen=10 applied and 10bp
#                                              fragments went to STAR. Williams
#                                              et al. 2016 (BMC Bioinformatics
#                                              17:103) showed aggressive trimming
#                                              shifted >10% of genes' expression
#                                              estimates, and that MOST of that
#                                              artefact disappears once a minimum
#                                              length filter is applied -- the
#                                              changes were "primarily being
#                                              driven by spurious mapping of short
#                                              reads". 36bp is the field default.
#  passes           2 (adapter    1            BBDuk applies adapter-trim BEFORE
#                   then quality)              quality-trim within a single pass,
#                                              which is the documented recommended
#                                              order. One pass is faster, halves
#                                              the intermediate fastqs, and emits
#                                              KTrimmed/QTrimmed stats separately
#                                              in one log.
#
#  --- ALIGNMENT (STAR) -----------------------------------------------------
#  outFilterScoreMinOverLread   0.2   0.66     *** RESTORED TO STAR DEFAULT ***
#  outFilterMatchNminOverLread  0.2   0.66     v1 relaxed these by 3x, so a read
#                                              only had to align over 20% of its
#                                              length. That inflates the mapping
#                                              rate with partial/spurious
#                                              alignments and multimappers. Dobin
#                                              (STAR author): "The default
#                                              parameters are close to the optimum
#                                              for expression quantification", and
#                                              "if the reads are properly trimmed,
#                                              then the filtering does not have to
#                                              be relaxed." v1 relaxed them
#                                              BECAUSE trimq=28 with no minlen was
#                                              producing very short reads -- fixing
#                                              the trimming removes the reason.
#                                              Both remain settable if a specific
#                                              library truly needs them.
#  outSAMtype        SAM          BAM Sorted   STAR sorts directly; removes the
#                                              multi-100GB intermediate SAM and the
#                                              separate `samtools sort` step.
#  outFilterMultimapNmax  (default 10)  10     unchanged (STAR default).
#  twopassMode       (none)       off          available via --two-pass; doubles
#                                              runtime, improves novel-junction
#                                              sensitivity. Off for DGE.
#
#  --- COUNTING (featureCounts) ---------------------------------------------
#  -p --countReadPairs   ON       ON           counts FRAGMENTS not reads. Since
#                                              subread 2.0.2, -p alone only
#                                              declares the data paired; without
#                                              --countReadPairs counts are ~2x.
#                                              Both flags are required.
#  -t exon -g gene_id    ON       ON           gene-level meta-features. unchanged.
#  -s                    0 (fixed) --stranded  see STRANDEDNESS below.
#
#  --- STRANDEDNESS DETECTION -----------------------------------------------
#  subsample         2%           2%           `samtools view -s 42.02`.
#  call thresholds   >0.8 fwd,    >=0.80 fwd,  v1 had NO "undetermined" branch:
#                    >0.8 rev,    >=0.80 rev,  anything not clearly stranded was
#                    else         |f1-f2|<0.10 called "unstranded", hiding gDNA
#                    "unstranded" -> unstranded contamination and failed libraries.
#                                 else         Thresholds now match nf-core/rnaseq
#                                 "undetermined" (stranded_threshold 0.8,
#                                              unstranded_threshold 0.1).
#  --stranded default  no          no          UNCHANGED ON PURPOSE. Auto-applying
#                                              per-sample detection would count
#                                              samples differently from each other,
#                                              which is worse than counting all of
#                                              them the same "wrong" way. For a
#                                              NEW project: run ONE sample, read
#                                              strand_verdict, then pass that value
#                                              explicitly for the whole batch.
#                                              `--stranded auto` exists but warns.
#  stranded_used     NOT RECORDED  RECORDED    v1 logged the DETECTED strandedness
#                                              but never what -s actually produced
#                                              the counts, so the DESeq2 manifest
#                                              could not state how counts were made.
#
#  --- QC FLAG THRESHOLDS (advisory only; nothing is deleted) ---------------
#  pct_surviving_trim  < 70   -> HIGH_TRIM_LOSS
#  pct_uniq_mapped     < 60   -> LOW_UNIQUE_MAPPING
#  pct_too_short       > 20   -> HIGH_TOO_SHORT      (classic under-trimmed adapter)
#  pct_multimapping    > 30   -> HIGH_MULTIMAPPING
#  pct_assigned        < 50   -> LOW_ASSIGNMENT      (often wrong -s, or gDNA)
#  detected != used           -> STRAND_MISMATCH
#
#  ==========================================================================
#  USAGE
#  ==========================================================================
#  Normally you do NOT call this directly -- use submit_rnaseq_fastq_to_counts.sh
#  to run a whole fastq folder. Call this only to (re)run a single sample.
#
#  MINIMAL:
#    sbatch rnaseq_fastq_to_counts.sh \
#           --ref grch38 \
#           --fastq-dir /lustre/.../fastqs \
#           --r1 Parental_rep1_R1_001.fastq.gz \
#           --r2 Parental_rep1_R2_001.fastq.gz \
#           --outdir /lustre/.../rnaseq_output
#
#  TYPICAL (reverse-stranded kit, clean up intermediates):
#    sbatch rnaseq_fastq_to_counts.sh \
#           --ref grch38 --fastq-dir /lustre/.../fastqs \
#           --r1 S1_R1_001.fastq.gz --r2 S1_R2_001.fastq.gz \
#           --outdir /lustre/.../rnaseq_output \
#           --stranded reverse --rm-trimmed-fastqs
#
#  FOLLOW MacManes EXACTLY (gentlest defensible trimming):
#    ... --trimq 5 --minlen 36
#
#  LEGACY v1 BEHAVIOUR (only to reproduce an old run -- NOT recommended):
#    ... --trimq 28 --minlen 10 \
#        --star-score-min-over-lread 0.2 --star-match-nmin-over-lread 0.2
#
#  DETECT STRANDEDNESS FIRST (recommended for a new project):
#    run ONE sample, then:  cat <outdir>/<rep>/<rep>_library_qc.tsv
#    and pass the strand_verdict as --stranded for the whole batch.
###############################################################################

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

# ---------------------------------------------------------------------------
# Help text
# ---------------------------------------------------------------------------
usage() {
cat <<'EOF'
------------------------------------------------------------------------------
rnaseq_fastq_to_counts.sh (v2)  --  run ONE paired-end sample (R1 + R2)
------------------------------------------------------------------------------
REQUIRED:
  --ref <hg19|grch38>       Reference genome build. Default: grch38
  --fastq-dir <FOLDER>      Folder containing the R1 and R2 fastq.gz files.
  --r1 <FILE>               R1 file NAME (must live in --fastq-dir).
  --r2 <FILE>               R2 file NAME (must live in --fastq-dir).

OUTPUT:
  --outdir <FOLDER>         Where results go. Default: submission folder.
  --fwd-keyword <STR>       Token marking the forward read in the file name; it
                            is stripped to name the output dir. Default: R1_001

LIBRARY:
  --stranded <no|yes|forward|reverse|auto>
                            featureCounts -s 0|1|2. Default: no
                            'auto' uses this sample's own detection -- do NOT use
                            it across a batch (see PARAMETER TABLE at top).

TRIMMING (BBDuk):
  --trimq <INT>             3' quality-trim threshold.        Default: 10
  --minlen <INT>            discard reads shorter than this.  Default: 36
  --qtrim <r|l|rl|f>        which end(s) to quality-trim.     Default: r
  --ktrim <r|l|f>           which end to adapter-trim.        Default: r
  --trim-k <INT>            adapter kmer length.              Default: 23
  --trim-mink <INT>         min kmer at read tip.             Default: 11
  --trim-hdist <INT>        kmer hamming distance.            Default: 1
  --maxns <INT>             discard reads with > this many Ns. Default: -1 (off)
  --ftm <INT>               force-trim length modulo.         Default: 5
  --no-tbo                  disable trim-by-overlap.
  --no-tpe                  disable trim-pairs-evenly.
  --adapters <FILE>         adapter FASTA. Default: bbmap resources/adapters.fa
  --skip-trim               skip trimming entirely (use raw fastqs).

ALIGNMENT (STAR):
  --star-score-min-over-lread <FLOAT>   Default: 0.66  (STAR default)
  --star-match-nmin-over-lread <FLOAT>  Default: 0.66  (STAR default)
  --star-multimap-nmax <INT>            Default: 10    (STAR default)
  --two-pass                            enable --twopassMode Basic (slow).
  --star-extra "<ARGS>"                 extra raw STAR arguments.

COUNTING (featureCounts):
  --fc-feature <STR>        -t value.  Default: exon
  --fc-attribute <STR>      -g value.  Default: gene_id
  --fc-extra "<ARGS>"       extra raw featureCounts arguments.
  --gtf <FILE>              override the annotation GTF chosen by --ref.
  --genome-dir <DIR>        override the STAR index chosen by --ref.

COVERAGE:
  --bigwig-binsize <INT>    bamCoverage --binSize. Default: 10
  --skip-bigwig             do not make the CPM bigwig.

HOUSEKEEPING:
  --rm-trimmed-fastqs       delete trimmed fastqs when finished.
  --rm-sam                  accepted for backwards compatibility (no-op: v2 never
                            writes a SAM).
  --threads <INT>           override $SLURM_CPUS_PER_TASK.
  --dry-run                 resolve + print the full configuration, then STOP.
  -h, --help                show this help.
------------------------------------------------------------------------------
EOF
}

# ---------------------------------------------------------------------------
# Defaults  (mirror the PARAMETER TABLE above)
# ---------------------------------------------------------------------------
REF="grch38"
FASTQ_DIR=""
R1=""
R2=""
OUTDIR=""
FWD_KEYWORD="R1_001"
STRANDED="no"

TRIMQ=10
MINLEN=36
QTRIM="r"
KTRIM="r"
TRIM_K=23
TRIM_MINK=11
TRIM_HDIST=1
MAXNS=-1
FTM=5
USE_TBO=true
USE_TPE=true
ADAPTERS_OVERRIDE=""
SKIP_TRIM=false

STAR_SCORE_MIN=0.66
STAR_MATCH_NMIN=0.66
STAR_MULTIMAP_NMAX=10
TWO_PASS=false
STAR_EXTRA=""

FC_FEATURE="exon"
FC_ATTRIBUTE="gene_id"
FC_EXTRA=""
GTF_OVERRIDE=""
GENOMEDIR_OVERRIDE=""

BIGWIG_BINSIZE=10
SKIP_BIGWIG=false

RM_TRIMMED=false
THREADS_OVERRIDE=""
DRY_RUN=false

# QC flag thresholds (see PARAMETER TABLE)
QC_MIN_SURVIVING=70
QC_MIN_UNIQ=60
QC_MAX_TOOSHORT=20
QC_MAX_MULTI=30
QC_MIN_ASSIGNED=50

# strandedness call thresholds (nf-core/rnaseq defaults)
STRAND_THRESHOLD=0.80
UNSTRAND_THRESHOLD=0.10

need_val() {
    [[ $# -ge 2 && -n "${2:-}" && "${2:0:1}" != "-" ]] || {
        echo "ERROR: option '$1' needs a value."; echo; usage; exit 1; }
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref)                        need_val "$@"; REF="$2";                  shift 2 ;;
        --fastq-dir|--fastqdir)       need_val "$@"; FASTQ_DIR="$2";            shift 2 ;;
        --r1|--R1)                    need_val "$@"; R1="$2";                   shift 2 ;;
        --r2|--R2)                    need_val "$@"; R2="$2";                   shift 2 ;;
        --outdir|--output-dir)        need_val "$@"; OUTDIR="$2";               shift 2 ;;
        --fwd-keyword|--fwd)          need_val "$@"; FWD_KEYWORD="$2";          shift 2 ;;
        --stranded)                   need_val "$@"; STRANDED="$2";             shift 2 ;;

        --trimq)                      need_val "$@"; TRIMQ="$2";                shift 2 ;;
        --minlen|--minlength)         need_val "$@"; MINLEN="$2";               shift 2 ;;
        --qtrim)                      need_val "$@"; QTRIM="$2";                shift 2 ;;
        --ktrim)                      need_val "$@"; KTRIM="$2";                shift 2 ;;
        --trim-k)                     need_val "$@"; TRIM_K="$2";               shift 2 ;;
        --trim-mink)                  need_val "$@"; TRIM_MINK="$2";            shift 2 ;;
        --trim-hdist)                 need_val "$@"; TRIM_HDIST="$2";           shift 2 ;;
        --maxns)                      need_val "$@"; MAXNS="$2";                shift 2 ;;
        --ftm)                        need_val "$@"; FTM="$2";                  shift 2 ;;
        --no-tbo)                     USE_TBO=false;                            shift ;;
        --no-tpe)                     USE_TPE=false;                            shift ;;
        --adapters)                   need_val "$@"; ADAPTERS_OVERRIDE="$2";    shift 2 ;;
        --skip-trim)                  SKIP_TRIM=true;                           shift ;;

        --star-score-min-over-lread)  need_val "$@"; STAR_SCORE_MIN="$2";       shift 2 ;;
        --star-match-nmin-over-lread) need_val "$@"; STAR_MATCH_NMIN="$2";      shift 2 ;;
        --star-multimap-nmax)         need_val "$@"; STAR_MULTIMAP_NMAX="$2";   shift 2 ;;
        --two-pass)                   TWO_PASS=true;                            shift ;;
        --star-extra)                 need_val "$@"; STAR_EXTRA="$2";           shift 2 ;;

        --fc-feature)                 need_val "$@"; FC_FEATURE="$2";           shift 2 ;;
        --fc-attribute)               need_val "$@"; FC_ATTRIBUTE="$2";         shift 2 ;;
        --fc-extra)                   need_val "$@"; FC_EXTRA="$2";             shift 2 ;;
        --gtf)                        need_val "$@"; GTF_OVERRIDE="$2";         shift 2 ;;
        --genome-dir)                 need_val "$@"; GENOMEDIR_OVERRIDE="$2";   shift 2 ;;

        --bigwig-binsize)             need_val "$@"; BIGWIG_BINSIZE="$2";       shift 2 ;;
        --skip-bigwig)                SKIP_BIGWIG=true;                         shift ;;

        --rm-trimmed-fastqs|--rm_trimmed_fastqs) RM_TRIMMED=true;               shift ;;
        --rm-sam|--rm_sam)
            echo "NOTE: --rm-sam is a no-op in v2 (STAR writes sorted BAM directly; no SAM exists)."
            shift ;;
        --threads)                    need_val "$@"; THREADS_OVERRIDE="$2";     shift 2 ;;
        --dry-run)                    DRY_RUN=true;                             shift ;;
        -h|--help)                    usage; exit 0 ;;
        *) echo "ERROR: unknown option '$1'"; echo; usage; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------
missing=()
[[ -n "$FASTQ_DIR" ]] || missing+=("--fastq-dir")
[[ -n "$R1" ]]        || missing+=("--r1")
[[ -n "$R2" ]]        || missing+=("--r2")
if [[ ${#missing[@]} -gt 0 ]]; then
    echo "ERROR: you did not provide the required option(s): ${missing[*]}"
    echo; usage; exit 1
fi

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
        echo; usage; exit 1 ;;
esac
[[ -n "$GTF_OVERRIDE" ]]        && GTF="$GTF_OVERRIDE"
[[ -n "$GENOMEDIR_OVERRIDE" ]]  && GENOME_DIR="$GENOMEDIR_OVERRIDE"

# map --stranded to featureCounts -s
STRANDED_LC=$(echo "$STRANDED" | tr '[:upper:]' '[:lower:]')
case "$STRANDED_LC" in
    no|unstranded)  STRAND_S=0 ;;
    yes|forward)    STRAND_S=1 ;;
    reverse)        STRAND_S=2 ;;
    auto)           STRAND_S=-1 ;;   # resolved after detection
    *) echo "ERROR: --stranded must be no|yes|forward|reverse|auto (got '$STRANDED')"; exit 1 ;;
esac

for pair in "TRIMQ:$TRIMQ" "MINLEN:$MINLEN" "TRIM_K:$TRIM_K" "TRIM_MINK:$TRIM_MINK" \
            "TRIM_HDIST:$TRIM_HDIST" "FTM:$FTM" "STAR_MULTIMAP_NMAX:$STAR_MULTIMAP_NMAX" \
            "BIGWIG_BINSIZE:$BIGWIG_BINSIZE"; do
    nm="${pair%%:*}"; vl="${pair#*:}"
    [[ "$vl" =~ ^[0-9]+$ ]] || { echo "ERROR: --${nm,,} must be a non-negative integer (got '$vl')"; exit 1; }
done
[[ "$MAXNS" =~ ^-?[0-9]+$ ]] || { echo "ERROR: --maxns must be an integer (got '$MAXNS')"; exit 1; }
for pair in "STAR_SCORE_MIN:$STAR_SCORE_MIN" "STAR_MATCH_NMIN:$STAR_MATCH_NMIN"; do
    nm="${pair%%:*}"; vl="${pair#*:}"
    [[ "$vl" =~ ^[0-9]*\.?[0-9]+$ ]] || { echo "ERROR: ${nm} must be a number 0-1 (got '$vl')"; exit 1; }
done

if [[ -z "$OUTDIR" ]]; then
    OUTDIR="$SUBMIT_DIR"
    echo "NOTE: no --outdir given -> writing results into: $OUTDIR"
fi
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

[[ -d "$FASTQ_DIR" ]]     || { echo "ERROR: --fastq-dir folder not found: $FASTQ_DIR"; exit 1; }
[[ -f "$FASTQ_DIR/$R1" ]] || { echo "ERROR: R1 file not found: $FASTQ_DIR/$R1"; exit 1; }
[[ -f "$FASTQ_DIR/$R2" ]] || { echo "ERROR: R2 file not found: $FASTQ_DIR/$R2"; exit 1; }
[[ "$R1" != "$R2" ]]      || { echo "ERROR: --r1 and --r2 are the same file ($R1)."; exit 1; }

# ---------------------------------------------------------------------------
# Tools / env
# ---------------------------------------------------------------------------
ADAPTERS="${ADAPTERS_OVERRIDE:-/lustre/isaac24/proj/UTK0009/RNA-Seq/bbmap/resources/adapters.fa}"
BBDUK=/lustre/isaac24/proj/UTK0009/RNA-Seq/conda_env/opt/bbmap-38.18/bbduk.sh
PROJ_ENVS=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration/conda-envs
CONDA_ENV="${PROJ_ENVS}/multiomics-py"           # deepTools (bamCoverage)
CONDA_ENV_SUBREAD="${PROJ_ENVS}/subread"         # featureCounts

THREADS=${THREADS_OVERRIDE:-${SLURM_CPUS_PER_TASK:-8}}
[[ "$THREADS" =~ ^[0-9]+$ ]] || { echo "ERROR: --threads must be an integer"; exit 1; }

# STAR BAM-sort RAM: ~60% of the job allocation, floor 8G, in bytes.
if [[ -n "${SLURM_MEM_PER_NODE:-}" && "${SLURM_MEM_PER_NODE}" =~ ^[0-9]+$ ]]; then
    BAMSORT_RAM=$(( SLURM_MEM_PER_NODE * 1024 * 1024 * 60 / 100 ))
else
    BAMSORT_RAM=$(( 32 * 1024 * 1024 * 1024 ))
fi
(( BAMSORT_RAM < 8589934592 )) && BAMSORT_RAM=8589934592

# ---------------------------------------------------------------------------
# Derive the replicate (directory) name from R1
# ---------------------------------------------------------------------------
noext=$(echo "${R1}" | sed -E 's/\.(fastq|fq)\.gz$//; s/\.(fastq|fq)$//')
Sample=$(echo "${noext}" | sed -E "s/(.*)[._-]?${FWD_KEYWORD}.*/\1/" | sed -E 's/[._-]+$//')
[[ -n "${Sample}" ]] || Sample="${noext}"

SampleDir="${OUTDIR}/${Sample}"
TrimDir="${SampleDir}/trimmed_fastqs"
LogDir="${SampleDir}/logs"
TmpDir="${SampleDir}/tmp"

BBDUK_ARGS=( -Xmx8g "threads=${THREADS}"
             "ref=${ADAPTERS}" "ktrim=${KTRIM}" "k=${TRIM_K}" "mink=${TRIM_MINK}"
             "hdist=${TRIM_HDIST}" "qtrim=${QTRIM}" "trimq=${TRIMQ}"
             "minlen=${MINLEN}" )
(( FTM > 0 ))     && BBDUK_ARGS+=( "ftm=${FTM}" )
(( MAXNS >= 0 ))  && BBDUK_ARGS+=( "maxns=${MAXNS}" )
$USE_TBO          && BBDUK_ARGS+=( tbo )
$USE_TPE          && BBDUK_ARGS+=( tpe )

echo "============================================================"
echo "rnaseq_fastq_to_counts.sh (v2)"
echo "Replicate  : ${Sample}"
echo "Reference  : ${REF_LC}"
echo "Started    : $(date)"
echo "Node       : $(hostname)"
echo "Threads    : ${THREADS}   Mem: ${SLURM_MEM_PER_NODE:-NA}"
echo "Fastq dir  : ${FASTQ_DIR}"
echo "R1/R2      : ${R1} / ${R2}"
echo "Output     : ${SampleDir}"
echo "Genome idx : ${GENOME_DIR}"
echo "GTF        : ${GTF}"
echo "------------------------------------------------------------"
echo "TRIM   : ${BBDUK_ARGS[*]}"
echo "STAR   : scoreMinOverLread=${STAR_SCORE_MIN} matchNminOverLread=${STAR_MATCH_NMIN}"
echo "         multimapNmax=${STAR_MULTIMAP_NMAX} twoPass=${TWO_PASS}"
echo "COUNT  : -t ${FC_FEATURE} -g ${FC_ATTRIBUTE} -s ${STRAND_S} (--stranded ${STRANDED_LC})"
echo "============================================================"

if $DRY_RUN; then
    echo "--dry-run set: configuration resolved, inputs verified. Stopping."
    exit 0
fi

[[ -d "$GENOME_DIR" ]] || { echo "ERROR: STAR index dir not found: $GENOME_DIR"; exit 1; }
[[ -f "$GTF" ]]        || { echo "ERROR: GTF not found: $GTF"; exit 1; }
if ! $SKIP_TRIM; then
    [[ -x "$BBDUK" ]]      || { echo "ERROR: bbduk.sh not executable: $BBDUK"; exit 1; }
    [[ -f "$ADAPTERS" ]]   || { echo "ERROR: adapter FASTA not found: $ADAPTERS"; exit 1; }
fi

mkdir -p "${SampleDir}" "${TrimDir}" "${LogDir}" "${TmpDir}"
cd "${SampleDir}"

Trim1=${TrimDir}/${Sample}_R1.trimmed.fastq.gz
Trim2=${TrimDir}/${Sample}_R2.trimmed.fastq.gz
StarPrefix=${Sample}_
SortedBam=${Sample}_sorted.bam
Counts=${Sample}_genecounts.tsv
FCraw=${Sample}_featurecounts.txt
BigWig=${Sample}.cpm.norm.bw
QCfile=${Sample}_library_qc.tsv
BBDukLog=${LogDir}/${Sample}_bbduk.log
StarLog=${Sample}_Log.final.out

# ---------------------------------------------------------------------------
# 1. Trimming -- SINGLE BBDuk pass (adapter first, then quality; documented order)
# ---------------------------------------------------------------------------
if $SKIP_TRIM; then
    echo ">> [1/6] trimming SKIPPED (--skip-trim); using raw fastqs"
    Trim1="${FASTQ_DIR}/${R1}"
    Trim2="${FASTQ_DIR}/${R2}"
    RAW_READS="NA"; KEPT_READS="NA"; PCT_SURVIVING="NA"
else
    echo ">> [1/6] BBDuk adapter + quality trimming (single pass)"
    "${BBDUK}" "${BBDUK_ARGS[@]}" \
        in1="${FASTQ_DIR}/${R1}" in2="${FASTQ_DIR}/${R2}" \
        out1="${Trim1}" out2="${Trim2}" \
        stats="${LogDir}/${Sample}_bbduk_stats.txt" \
        2> >(tee "${BBDukLog}" >&2)

    RAW_READS=$(awk '/^Input:/{print $2; exit}'  "${BBDukLog}")
    KEPT_READS=$(awk '/^Result:/{print $2; exit}' "${BBDukLog}")
    RAW_READS=${RAW_READS:-NA}; KEPT_READS=${KEPT_READS:-NA}
    if [[ "$RAW_READS" =~ ^[0-9]+$ && "$KEPT_READS" =~ ^[0-9]+$ && "$RAW_READS" -gt 0 ]]; then
        PCT_SURVIVING=$(awk -v a="$KEPT_READS" -v b="$RAW_READS" 'BEGIN{printf "%.2f", 100*a/b}')
    else
        PCT_SURVIVING="NA"
    fi
    echo "   reads in=${RAW_READS}  kept=${KEPT_READS}  (${PCT_SURVIVING}% surviving)"
fi

# ---------------------------------------------------------------------------
# 2. STAR alignment -> coordinate-sorted BAM (no intermediate SAM)
# ---------------------------------------------------------------------------
echo ">> [2/6] STAR alignment"
module load star/2.7.6a
STAR_ARGS=( --runThreadN "${THREADS}"
            --genomeDir "${GENOME_DIR}"
            --readFilesIn "${Trim1}" "${Trim2}"
            --readFilesCommand zcat
            --outFilterScoreMinOverLread  "${STAR_SCORE_MIN}"
            --outFilterMatchNminOverLread "${STAR_MATCH_NMIN}"
            --outFilterMultimapNmax "${STAR_MULTIMAP_NMAX}"
            --outSAMtype BAM SortedByCoordinate
            --limitBAMsortRAM "${BAMSORT_RAM}"
            --outTmpDir "${TmpDir}/star"
            --outFileNamePrefix "${StarPrefix}" )
$TWO_PASS && STAR_ARGS+=( --twopassMode Basic )
# shellcheck disable=SC2206  # deliberate word-splitting of user-supplied extra args
[[ -n "$STAR_EXTRA" ]] && STAR_ARGS+=( ${STAR_EXTRA} )

rm -rf "${TmpDir}/star"          # STAR refuses to start if outTmpDir exists
STAR "${STAR_ARGS[@]}"
mv "${StarPrefix}Aligned.sortedByCoord.out.bam" "${SortedBam}"

# ---------------------------------------------------------------------------
# 3. Index
# ---------------------------------------------------------------------------
echo ">> [3/6] index BAM"
module load samtools/1.16.1-gcc
samtools index -@ "${THREADS}" "${SortedBam}"

star_stat() {   # star_stat "<label substring>"  -> trimmed value or NA
    awk -F'|' -v k="$1" 'index($0,k){gsub(/^[ \t]+|[ \t]+$/,"",$2); print $2; exit}' "${StarLog}"
}
num_only() { echo "${1:-}" | tr -d '%' | grep -E '^[0-9]+(\.[0-9]+)?$' || echo "NA"; }

STAR_INPUT=$(num_only "$(star_stat 'Number of input reads')")
UNIQ_N=$(num_only     "$(star_stat 'Uniquely mapped reads number')")
PCT_UNIQ=$(num_only   "$(star_stat 'Uniquely mapped reads %')")
PCT_MULTI=$(num_only  "$(star_stat '% of reads mapped to multiple loci')")
PCT_TOOMANY=$(num_only "$(star_stat '% of reads mapped to too many loci')")
PCT_TOOSHORT=$(num_only "$(star_stat '% of reads unmapped: too short')")
MISMATCH_RATE=$(num_only "$(star_stat 'Mismatch rate per base')")
AVG_IN_LEN=$(num_only  "$(star_stat 'Average input read length')")
AVG_MAP_LEN=$(num_only "$(star_stat 'Average mapped length')")
echo "   uniquely mapped ${PCT_UNIQ}%  multi ${PCT_MULTI}%  too-short ${PCT_TOOSHORT}%"

# ---------------------------------------------------------------------------
# 4. Strandedness detection (reports; only changes -s when --stranded auto)
# ---------------------------------------------------------------------------
echo ">> [4/6] strandedness detection (2% subsample, featureCounts -s 0/1/2)"

module purge
module load anaconda3/2024.06
# shellcheck disable=SC1090  # path supplied by the module system at runtime
source "$ANACONDA3_SH"
unset PYTHONPATH                # unset BEFORE activation (project convention)
export PYTHONNOUSERSITE=1
conda activate "${CONDA_ENV_SUBREAD}"

SubBam=${TmpDir}/${Sample}_strandcheck.bam
samtools view -@ "${THREADS}" -s 42.02 -b "${SortedBam}" > "${SubBam}"

declare -a SASSIGN=(0 0 0)
for s in 0 1 2; do
    featureCounts -a "${GTF}" -o "${TmpDir}/sc_s${s}.txt" \
        -p --countReadPairs -t "${FC_FEATURE}" -g "${FC_ATTRIBUTE}" \
        -s "${s}" -T "${THREADS}" "${SubBam}" > /dev/null 2>&1 || true
    v=$(awk -F'\t' '$1=="Assigned"{print $2; exit}' "${TmpDir}/sc_s${s}.txt.summary" 2>/dev/null || true)
    SASSIGN[$s]=${v:-0}
done

read -r FRAC_FWD FRAC_REV STRAND_VERDICT <<< "$(
  awk -v a1="${SASSIGN[1]}" -v a2="${SASSIGN[2]}" \
      -v st="${STRAND_THRESHOLD}" -v ut="${UNSTRAND_THRESHOLD}" 'BEGIN{
    tot = a1 + a2;
    if (tot <= 0) { print "NA NA undetermined"; exit }
    f1 = a1/tot; f2 = a2/tot;
    d  = (f1 > f2) ? f1 - f2 : f2 - f1;
    if      (f1 >= st) v = "forward";
    else if (f2 >= st) v = "reverse";
    else if (d  <  ut) v = "unstranded";
    else               v = "undetermined";
    printf "%.4f %.4f %s", f1, f2, v;
  }')"

echo "   assigned: s0=${SASSIGN[0]} s1=${SASSIGN[1]} s2=${SASSIGN[2]}"
echo "   fwd=${FRAC_FWD} rev=${FRAC_REV} -> ${STRAND_VERDICT}"
if [[ "${STRAND_VERDICT}" == "undetermined" ]]; then
    echo "   WARNING: strandedness is UNDETERMINED - neither clearly stranded nor"
    echo "            clearly unstranded. Common causes: genomic DNA carry-over,"
    echo "            very low library complexity, or a failed library."
fi

# resolve --stranded auto
if [[ "${STRAND_S}" -eq -1 ]]; then
    case "${STRAND_VERDICT}" in
        forward)    STRAND_S=1 ;;
        reverse)    STRAND_S=2 ;;
        unstranded) STRAND_S=0 ;;
        *)          STRAND_S=0
                    echo "   WARNING: --stranded auto could not resolve; falling back to -s 0." ;;
    esac
    echo "   WARNING: --stranded auto resolved THIS sample to -s ${STRAND_S}."
    echo "            Across a batch this can count samples differently from one"
    echo "            another, which is worse than counting all of them the same way."
    echo "            Prefer an explicit --stranded value for the whole project."
fi

case "${STRAND_S}" in
    0) STRANDED_USED="no" ;;
    1) STRANDED_USED="forward" ;;
    2) STRANDED_USED="reverse" ;;
    *) echo "ERROR: unresolved strandedness (-s ${STRAND_S})"; exit 1 ;;
esac

STRAND_MATCH="yes"
if [[ "${STRAND_VERDICT}" == "forward"    && "${STRANDED_USED}" != "forward" ]] ||
   [[ "${STRAND_VERDICT}" == "reverse"    && "${STRANDED_USED}" != "reverse" ]] ||
   [[ "${STRAND_VERDICT}" == "unstranded" && "${STRANDED_USED}" != "no"      ]]; then
    STRAND_MATCH="no"
    echo "   WARNING: detected '${STRAND_VERDICT}' but counting with '${STRANDED_USED}' (-s ${STRAND_S})."
    echo "            This is deliberate when a project mixes preps: all samples must be"
    echo "            counted the SAME way, not each one 'correctly'."
fi
rm -f "${SubBam}" "${TmpDir}"/sc_s*.txt "${TmpDir}"/sc_s*.txt.summary

# ---------------------------------------------------------------------------
# 5. featureCounts
# ---------------------------------------------------------------------------
echo ">> [5/6] featureCounts (-s ${STRAND_S}, -t ${FC_FEATURE}, -g ${FC_ATTRIBUTE})"
FC_ARGS=( -a "${GTF}" -o "${FCraw}" -p --countReadPairs
          -t "${FC_FEATURE}" -g "${FC_ATTRIBUTE}" -s "${STRAND_S}" -T "${THREADS}" )
# shellcheck disable=SC2206  # deliberate word-splitting of user-supplied extra args
[[ -n "$FC_EXTRA" ]] && FC_ARGS+=( ${FC_EXTRA} )
featureCounts "${FC_ARGS[@]}" "${SortedBam}"

# htseq-style 2-column file (+ __ summary rows) for the DESeq2 script
tail -n +3 "${FCraw}" | cut -f1,7 > "${Counts}"
awk -F'\t' '
  $1=="Unassigned_NoFeatures"     {nf=$2}
  $1=="Unassigned_Ambiguity"      {amb=$2}
  $1=="Unassigned_MappingQuality" {mq=$2}
  $1=="Unassigned_Unmapped"       {un=$2}
  $1=="Unassigned_MultiMapping"   {mm=$2}
  END{
    printf "__no_feature\t%d\n",           nf+0
    printf "__ambiguous\t%d\n",            amb+0
    printf "__too_low_aQual\t%d\n",        mq+0
    printf "__not_aligned\t%d\n",          un+0
    printf "__alignment_not_unique\t%d\n", mm+0
  }' "${FCraw}.summary" >> "${Counts}"

ASSIGNED=$(awk -F'\t' '$1=="Assigned"{print $2; exit}' "${FCraw}.summary")
ASSIGNED=${ASSIGNED:-0}
FC_TOTAL=$(awk -F'\t' 'NR>1{s+=$2} END{print s+0}' "${FCraw}.summary")
if [[ "$FC_TOTAL" -gt 0 ]]; then
    PCT_ASSIGNED=$(awk -v a="$ASSIGNED" -v b="$FC_TOTAL" 'BEGIN{printf "%.2f", 100*a/b}')
else
    PCT_ASSIGNED="NA"
fi
# features with a non-zero count (the "counts that survived" per gene)
FEATURES_DETECTED=$(awk -F'\t' '$1 !~ /^__/ && $2+0 > 0 {n++} END{print n+0}' "${Counts}")
FEATURES_TOTAL=$(awk -F'\t' '$1 !~ /^__/ {n++} END{print n+0}' "${Counts}")
echo "   assigned=${ASSIGNED} (${PCT_ASSIGNED}%)  features detected=${FEATURES_DETECTED}/${FEATURES_TOTAL}"

conda deactivate

# ---------------------------------------------------------------------------
# 6. CPM bigwig
# ---------------------------------------------------------------------------
if $SKIP_BIGWIG; then
    echo ">> [6/6] bigwig SKIPPED (--skip-bigwig)"
else
    echo ">> [6/6] bigwig (CPM, binSize ${BIGWIG_BINSIZE})"
    # NOTE: do NOT `module purge` here -- that unloads anaconda3 and can strip
    # conda from PATH on a fresh node. Switch envs within the same conda install.
    conda activate "${CONDA_ENV}"
    bamCoverage --bam "${SortedBam}" --outFileName "${BigWig}" \
        --binSize "${BIGWIG_BINSIZE}" --normalizeUsing CPM \
        --numberOfProcessors "${THREADS}"
    conda deactivate
fi

# ---------------------------------------------------------------------------
# QC record -- ONE row. Read by the DESeq2 script and by the batch submitter
# when it builds the sample sheet.
# ---------------------------------------------------------------------------
qc_flags=""
add_flag() { qc_flags="${qc_flags:+${qc_flags};}$1"; }
gt() { awk -v a="$1" -v b="$2" 'BEGIN{exit !(a+0 >  b+0)}'; }
lt() { awk -v a="$1" -v b="$2" 'BEGIN{exit !(a+0 <  b+0)}'; }

[[ "$PCT_SURVIVING" != "NA" ]] && lt "$PCT_SURVIVING" "$QC_MIN_SURVIVING" && add_flag "HIGH_TRIM_LOSS"
[[ "$PCT_UNIQ"      != "NA" ]] && lt "$PCT_UNIQ"      "$QC_MIN_UNIQ"      && add_flag "LOW_UNIQUE_MAPPING"
[[ "$PCT_TOOSHORT"  != "NA" ]] && gt "$PCT_TOOSHORT"  "$QC_MAX_TOOSHORT"  && add_flag "HIGH_TOO_SHORT"
[[ "$PCT_MULTI"     != "NA" ]] && gt "$PCT_MULTI"     "$QC_MAX_MULTI"     && add_flag "HIGH_MULTIMAPPING"
[[ "$PCT_ASSIGNED"  != "NA" ]] && lt "$PCT_ASSIGNED"  "$QC_MIN_ASSIGNED"  && add_flag "LOW_ASSIGNMENT"
[[ "$STRAND_MATCH"  == "no"  ]] && add_flag "STRAND_MISMATCH"
[[ "$STRAND_VERDICT" == "undetermined" ]] && add_flag "STRAND_UNDETERMINED"
qc_flags="${qc_flags:-PASS}"

TRIM_PARAMS=$($SKIP_TRIM && echo "skipped" || echo "${BBDUK_ARGS[*]}" | sed 's/-Xmx8g //; s/threads=[0-9]* //')
STAR_PARAMS="scoreMinOverLread=${STAR_SCORE_MIN};matchNminOverLread=${STAR_MATCH_NMIN};multimapNmax=${STAR_MULTIMAP_NMAX};twoPass=${TWO_PASS}"
FC_PARAMS="-p --countReadPairs -t ${FC_FEATURE} -g ${FC_ATTRIBUTE} -s ${STRAND_S}"

{
printf "sample\tref\tgtf\traw_reads\tkept_reads\tpct_surviving_trim\tstar_input_reads\t"
printf "uniq_mapped\tpct_uniq_mapped\tpct_multimapping\tpct_too_many_loci\tpct_too_short\t"
printf "mismatch_rate_per_base\tavg_input_read_len\tavg_mapped_len\t"
printf "strand_s0\tstrand_s1\tstrand_s2\tfrac_forward\tfrac_reverse\tstrand_verdict\t"
printf "stranded_used\tstrand_s_flag\tstrand_matches_detected\t"
printf "assigned_fragments\tpct_assigned\tfeatures_detected\tfeatures_total\t"
printf "trim_params\tstar_params\tfc_params\tqc_flag\n"
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t" \
  "${Sample}" "${REF_LC}" "${GTF}" "${RAW_READS}" "${KEPT_READS}" "${PCT_SURVIVING}" "${STAR_INPUT}"
printf "%s\t%s\t%s\t%s\t%s\t" \
  "${UNIQ_N}" "${PCT_UNIQ}" "${PCT_MULTI}" "${PCT_TOOMANY}" "${PCT_TOOSHORT}"
printf "%s\t%s\t%s\t" "${MISMATCH_RATE}" "${AVG_IN_LEN}" "${AVG_MAP_LEN}"
printf "%s\t%s\t%s\t%s\t%s\t%s\t" \
  "${SASSIGN[0]}" "${SASSIGN[1]}" "${SASSIGN[2]}" "${FRAC_FWD}" "${FRAC_REV}" "${STRAND_VERDICT}"
printf "%s\t%s\t%s\t" "${STRANDED_USED}" "${STRAND_S}" "${STRAND_MATCH}"
printf "%s\t%s\t%s\t%s\t" "${ASSIGNED}" "${PCT_ASSIGNED}" "${FEATURES_DETECTED}" "${FEATURES_TOTAL}"
printf "%s\t%s\t%s\t%s\n" "${TRIM_PARAMS}" "${STAR_PARAMS}" "${FC_PARAMS}" "${qc_flags}"
} > "${QCfile}"

# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------
if [[ "${RM_TRIMMED}" == true ]] && ! $SKIP_TRIM; then
    echo ">> removing trimmed fastqs"
    rm -f "${Trim1}" "${Trim2}"
fi
rm -rf "${TmpDir}"
rmdir "${TrimDir}" 2>/dev/null || true

echo "============================================================"
echo "Finished  : ${Sample} at $(date)"
echo "QC flag   : ${qc_flags}"
echo "Kept      : ${SortedBam}, ${Counts}, ${QCfile}$( $SKIP_BIGWIG || echo ", ${BigWig}")"
echo "============================================================"