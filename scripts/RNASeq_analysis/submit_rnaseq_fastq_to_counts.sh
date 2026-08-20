#!/bin/bash
set -euo pipefail

###############################################################################
# submit_rnaseq_fastq_to_counts.sh   (v2)   -- BATCH submitter + SAMPLE SHEET
#
# 1. Finds every R1/R2 pair in a fastq folder and submits ONE SLURM job per pair
#    (rnaseq_fastq_to_counts.sh), all writing into a shared --outdir.
# 2. Submits a final "collector" job that waits for all of them (--dependency)
#    and writes a DESeq2-ready SAMPLE SHEET with the `batch` column left BLANK
#    for you to fill in, plus per-sample QC columns.
#
#  ==========================================================================
#  WHAT YOU GET IN <outdir>
#  ==========================================================================
#   samples_manifest.tsv      written at SUBMIT time: replicate -> sample, R1, R2
#   samplesheet.csv           written at COLLECT time. THIS is what you hand to
#                             rnaseq_deseq2_dge.R via --samplesheet.
#   sample_qc_summary.tsv     every QC column from every sample, one row each
#   SAMPLESHEET_README.txt    what to edit and what each column means
#
#  samplesheet.csv columns
#  -----------------------
#    sample        replicate name = the per-sample output directory name.
#                  rnaseq_deseq2_dge.R resolves count files from this.
#    condition     GUESSED by stripping the --rep-keyword token from `sample`
#                  (Parental_rep1 -> Parental). *** ALWAYS CHECK THIS ***
#                  Sanitised to [A-Za-z0-9._] because DESeq2 factor levels may
#                  contain only letters, numbers, underscore and period; DESeq2
#                  silently renames anything else, after which
#                  contrast=c("condition", treat, control) no longer matches.
#    batch         LEFT EMPTY ON PURPOSE -- you fill it in. Any grouping labels
#                  work: 1/2/3, a/b/c, kitA/kitB, Jan/Feb. One label per group.
#                  Then run the DESeq2 script with --batch batch.
#                  Leave it empty and simply do not pass --batch to ignore it.
#    path          absolute path to the count file (makes the sheet portable).
#    status        OK | FAILED | MISSING_COUNTS  (see collector logic)
#    qc_flag       PASS, or semicolon-joined warnings from the per-sample run
#    ... plus every QC column: raw_reads, pct_surviving_trim, pct_uniq_mapped,
#        pct_too_short, assigned_fragments, pct_assigned, features_detected,
#        strand_verdict, stranded_used, ...
#    Extra columns are harmless: rnaseq_deseq2_dge.R requires only `sample` and
#    `condition` and carries the rest through untouched.
#
#  ==========================================================================
#  WORKFLOW
#  ==========================================================================
#    1) ./submit_rnaseq_fastq_to_counts.sh --ref grch38 \
#           --fastq-dir /lustre/.../fastqs --outdir /lustre/.../rnaseq_output
#    2) wait  (squeue -u $USER)
#    3) open <outdir>/samplesheet.csv, check `condition`, fill in `batch`
#    4) sbatch submit_rnaseq_deseq2_dge.sh --ref grch38 \
#           --counts-dir /lustre/.../rnaseq_output \
#           --samplesheet /lustre/.../rnaseq_output/samplesheet.csv \
#           --control Parental --treatment Bottom10 --batch batch \
#           --outdir results/exploratory/deseq2_B10_vs_Parental --dry-run
#
#  ==========================================================================
#  USAGE EXAMPLES
#  ==========================================================================
#  # Standard Illumina naming (Sample_rep1_R1_001.fastq.gz), grch38
#  ./submit_rnaseq_fastq_to_counts.sh --ref grch38 \
#        --fastq-dir /lustre/.../fastqs \
#        --outdir /lustre/.../rnaseq_output --rm-trimmed-fastqs
#
#  # Reverse-stranded kit (TruSeq Stranded mRNA / NEBNext Ultra II Directional)
#  #   -- v1 COULD NOT DO THIS: --stranded was not accepted by the submitter, so
#  #      every batch run silently counted at -s 0.
#  ./submit_rnaseq_fastq_to_counts.sh --ref grch38 \
#        --fastq-dir /lustre/.../fastqs --outdir /lustre/.../rnaseq_output \
#        --stranded reverse
#
#  # Vendor naming, e.g. Ctrl_Rep1_R1.fastq.gz
#  ./submit_rnaseq_fastq_to_counts.sh --ref hg19 \
#        --fastq-dir /lustre/.../fastqs \
#        --fwd-keyword R1 --rev-keyword R2 --rep-keyword Rep
#
#  # Follow MacManes 2014 trimming exactly
#  ./submit_rnaseq_fastq_to_counts.sh --ref grch38 --fastq-dir DIR --outdir OUT \
#        --trimq 5 --minlen 36
#
#  # See what would be submitted, submit nothing
#  ./submit_rnaseq_fastq_to_counts.sh --ref grch38 --fastq-dir DIR --dry-run
#
#  # Rebuild the sample sheet later (jobs already finished, or you added samples)
#  ./submit_rnaseq_fastq_to_counts.sh --collect-only --outdir /lustre/.../rnaseq_output
#
#  Every pipeline flag is forwarded unchanged -- run
#  `./rnaseq_fastq_to_counts.sh --help` for the full parameter table.
###############################################################################

SUBMIT_DIR="$PWD"

usage() {
cat <<'EOF'
------------------------------------------------------------------------------
RNA-seq BATCH submitter (v2) -- one job per sample + DESeq2 sample sheet
------------------------------------------------------------------------------
REQUIRED (unless --collect-only):
  --ref <hg19|grch38>       Reference build. Default: grch38
  --fastq-dir <FOLDER>      Folder containing the fastq.gz files.

OUTPUT / NAMING:
  --outdir <FOLDER>         Shared output folder. Default: current folder.
  --fwd-keyword <STR>       Token marking the FORWARD read.   Default: R1_001
  --rev-keyword <STR>       Token marking the REVERSE read.   Default: R2_001
  --rep-keyword <STR>       Token marking the REPLICATE number; stripped to
                            guess `condition`.                Default: rep

SAMPLE SHEET:
  --samplesheet-name <STR>  Default: samplesheet.csv
  --no-samplesheet          Do not submit the collector job.
  --collect-only            Do NOT submit alignment jobs; just (re)build the
                            sample sheet from whatever is already in --outdir.
  --batch-column <STR>      Name of the empty batch column. Default: batch
  --condition-from <rep|dir>
                            rep : strip <sep><rep-keyword><digits> from the
                                  replicate name  (Parental_rep1 -> Parental)
                            dir : use the replicate name unchanged (every sample
                                  becomes its own condition; edit by hand)
                            Default: rep

SLURM (for the collector job only; sample jobs use the pipeline's own header):
  --account <STR>           Default: ACF-UTK0009
  --partition <STR>         Default: condo-rmccord
  --qos <STR>               Default: condo
  --collector-time <STR>    Default: 00:30:00

PIPELINE FLAGS -- all forwarded unchanged to rnaseq_fastq_to_counts.sh:
  library     --stranded <no|yes|forward|reverse|auto>
  trimming    --trimq --minlen --qtrim --ktrim --trim-k --trim-mink
              --trim-hdist --maxns --ftm --no-tbo --no-tpe --adapters --skip-trim
  alignment   --star-score-min-over-lread --star-match-nmin-over-lread
              --star-multimap-nmax --two-pass --star-extra
  counting    --fc-feature --fc-attribute --fc-extra --gtf --genome-dir
  coverage    --bigwig-binsize --skip-bigwig
  cleanup     --rm-trimmed-fastqs --rm-sam
  resources   --threads

OTHER:
  --dry-run                 Resolve everything, print the plan, submit nothing.
  -h, --help                Show this help.

Track your jobs with:  squeue -u $USER
------------------------------------------------------------------------------
EOF
}

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
REF="grch38"
FASTQ_DIR=""
OUTDIR=""
FWD_KEYWORD="R1_001"
REV_KEYWORD="R2_001"
REP_KEYWORD="rep"

SHEET_NAME="samplesheet.csv"
MAKE_SHEET=true
COLLECT_ONLY=false
BATCH_COLUMN="batch"
CONDITION_FROM="rep"

SLURM_ACCOUNT="ACF-UTK0009"
SLURM_PARTITION="condo-rmccord"
SLURM_QOS="condo"
COLLECTOR_TIME="00:30:00"

DRY_RUN=false
PASSTHRU=()          # pipeline flags forwarded verbatim

need_val() {
    [[ $# -ge 2 && -n "${2:-}" && "${2:0:1}" != "-" ]] || {
        echo "ERROR: option '$1' needs a value."; echo; usage; exit 1; }
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        # ---- submitter's own options ----
        --ref)                    need_val "$@"; REF="$2";              shift 2 ;;
        --fastq-dir|--fastqdir)   need_val "$@"; FASTQ_DIR="$2";        shift 2 ;;
        --outdir|--output-dir)    need_val "$@"; OUTDIR="$2";           shift 2 ;;
        --fwd-keyword|--fwd)      need_val "$@"; FWD_KEYWORD="$2";      shift 2 ;;
        --rev-keyword|--rev)      need_val "$@"; REV_KEYWORD="$2";      shift 2 ;;
        --rep-keyword|--rep)      need_val "$@"; REP_KEYWORD="$2";      shift 2 ;;
        --samplesheet-name)       need_val "$@"; SHEET_NAME="$2";       shift 2 ;;
        --no-samplesheet)         MAKE_SHEET=false;                     shift ;;
        --collect-only)           COLLECT_ONLY=true;                    shift ;;
        --batch-column)           need_val "$@"; BATCH_COLUMN="$2";     shift 2 ;;
        --condition-from)         need_val "$@"; CONDITION_FROM="$2";   shift 2 ;;
        --account)                need_val "$@"; SLURM_ACCOUNT="$2";    shift 2 ;;
        --partition)              need_val "$@"; SLURM_PARTITION="$2";  shift 2 ;;
        --qos)                    need_val "$@"; SLURM_QOS="$2";        shift 2 ;;
        --collector-time)         need_val "$@"; COLLECTOR_TIME="$2";   shift 2 ;;
        --dry-run)                DRY_RUN=true;                         shift ;;
        -h|--help)                usage; exit 0 ;;

        # ---- pipeline flags: value-taking, forwarded ----
        --stranded|--trimq|--minlen|--minlength|--qtrim|--ktrim|--trim-k|\
        --trim-mink|--trim-hdist|--maxns|--ftm|--adapters|\
        --star-score-min-over-lread|--star-match-nmin-over-lread|\
        --star-multimap-nmax|--star-extra|\
        --fc-feature|--fc-attribute|--fc-extra|--gtf|--genome-dir|\
        --bigwig-binsize|--threads)
            [[ $# -ge 2 ]] || { echo "ERROR: option '$1' needs a value."; exit 1; }
            PASSTHRU+=( "$1" "$2" ); shift 2 ;;

        # ---- pipeline flags: boolean, forwarded ----
        --no-tbo|--no-tpe|--skip-trim|--two-pass|--skip-bigwig|\
        --rm-trimmed-fastqs|--rm_trimmed_fastqs|--rm-sam|--rm_sam)
            PASSTHRU+=( "$1" ); shift ;;

        *) echo "ERROR: unknown option '$1'"; echo; usage; exit 1 ;;
    esac
done

case "$CONDITION_FROM" in
    rep|dir) : ;;
    *) echo "ERROR: --condition-from must be 'rep' or 'dir' (got '$CONDITION_FROM')"; exit 1 ;;
esac

REF_LC=$(echo "$REF" | tr '[:upper:]' '[:lower:]')
if ! $COLLECT_ONLY; then
    [[ -n "$FASTQ_DIR" ]] || { echo "ERROR: --fastq-dir is required."; echo; usage; exit 1; }
    case "$REF_LC" in
        hg19|grch38) : ;;
        *) echo "ERROR: --ref must be exactly 'hg19' or 'grch38' (you gave '$REF')."; exit 1 ;;
    esac
    [[ -d "$FASTQ_DIR" ]] || { echo "ERROR: --fastq-dir folder not found: $FASTQ_DIR"; exit 1; }
    FASTQ_DIR="$(cd "$FASTQ_DIR" && pwd)"
fi

if [[ -z "$OUTDIR" ]]; then
    OUTDIR="$SUBMIT_DIR"
    echo "NOTE: no --outdir given -> results go into: $OUTDIR"
fi
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

SCRIPT_SELF="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(dirname "$SCRIPT_SELF")"
MAIN_SCRIPT="${SCRIPT_DIR}/rnaseq_fastq_to_counts.sh"

# ---------------------------------------------------------------------------
# Name helpers
# ---------------------------------------------------------------------------
rep_name() {
    local n; n=$(echo "$1" | sed -E 's/\.(fastq|fq)\.gz$//; s/\.(fastq|fq)$//')
    n=$(echo "$n" | sed -E "s/(.*)[._-]?${FWD_KEYWORD}.*/\1/" | sed -E 's/[._-]+$//')
    echo "$n"
}
mate_name() { echo "$1" | sed -E "s/(.*)${FWD_KEYWORD}/\1${REV_KEYWORD}/"; }
samp_name() { echo "$1" | sed -E "s/(.*)[._-]${REP_KEYWORD}[0-9]+.*/\1/I" | sed -E 's/[._-]+$//'; }

# DESeq2 factor levels may contain only letters, numbers, underscore and period.
# Anything else is renamed by DESeq2 itself, after which the contrast built from
# the ORIGINAL string no longer matches any level and results() dies -- typically
# an hour into the run, after DESeq() has already been fitted.
sanitize_level() {
    local s="$1" out
    out=$(echo "$s" | sed -E 's/[^A-Za-z0-9._]+/_/g; s/^_+//; s/_+$//')
    [[ -n "$out" ]] || out="level"
    [[ "$out" =~ ^[0-9] ]] && out="X${out}"
    echo "$out"
}

# ===========================================================================
# COLLECTOR -- build the sample sheet from an --outdir full of finished samples
# ===========================================================================
collect_samplesheet() {
    local sheet="${OUTDIR}/${SHEET_NAME}"
    local qcsum="${OUTDIR}/sample_qc_summary.tsv"
    local readme="${OUTDIR}/SAMPLESHEET_README.txt"

    echo "------------------------------------------------------------"
    echo "Building sample sheet from: ${OUTDIR}"

    # Discover finished replicates: any subdir with a *_genecounts.tsv OR a
    # *_library_qc.tsv (so failures are visible instead of silently dropped).
    local -a reps=()
    local d rep
    for d in "${OUTDIR}"/*/; do
        [[ -d "$d" ]] || continue
        rep="$(basename "$d")"
        [[ "$rep" == "logs" || "$rep" == "tmp" ]] && continue
        if [[ -f "${d}${rep}_genecounts.tsv" || -f "${d}${rep}_library_qc.tsv" ]]; then
            reps+=( "$rep" )
        fi
    done

    if [[ ${#reps[@]} -eq 0 ]]; then
        echo "ERROR: no per-sample result directories found under ${OUTDIR}"
        echo "       expected <outdir>/<replicate>/<replicate>_genecounts.tsv"
        return 1
    fi
    IFS=$'\n' read -r -d '' -a reps < <(printf '%s\n' "${reps[@]}" | sort && printf '\0')
    echo "Found ${#reps[@]} replicate director(ies)."

    # ---- QC header: take it from the first QC file that exists ----
    local qc_header="" first_qc=""
    for rep in "${reps[@]}"; do
        if [[ -f "${OUTDIR}/${rep}/${rep}_library_qc.tsv" ]]; then
            first_qc="${OUTDIR}/${rep}/${rep}_library_qc.tsv"
            qc_header=$(head -1 "$first_qc"); break
        fi
    done

    # ---- sample_qc_summary.tsv ----
    if [[ -n "$qc_header" ]]; then
        { echo "$qc_header"
          for rep in "${reps[@]}"; do
              local f="${OUTDIR}/${rep}/${rep}_library_qc.tsv"
              [[ -f "$f" ]] && tail -n +2 "$f"
          done
        } > "$qcsum"
        echo "QC summary -> ${qcsum}"
    else
        echo "WARNING: no *_library_qc.tsv found; QC columns will be empty."
        : > "$qcsum"
    fi

    # ---- which QC columns to carry into the sheet, in this order ----
    local -a want=( raw_reads pct_surviving_trim star_input_reads pct_uniq_mapped
                    pct_multimapping pct_too_short assigned_fragments pct_assigned
                    features_detected features_total strand_verdict stranded_used
                    strand_matches_detected qc_flag )
    local -a have=() idx=()
    if [[ -n "$qc_header" ]]; then
        local w i
        for w in "${want[@]}"; do
            i=$(awk -F'\t' -v k="$w" 'NR==1{for(j=1;j<=NF;j++) if($j==k){print j; exit}}' <<< "$qc_header")
            if [[ -n "$i" ]]; then have+=( "$w" ); idx+=( "$i" ); fi
        done
    fi

    # ---- write the sheet ----
    { printf "sample,condition,%s,path,status" "$BATCH_COLUMN"
      for w in "${have[@]}"; do printf ",%s" "$w"; done
      printf "\n"
    } > "$sheet"

    local n_ok=0 n_bad=0 renamed=0
    local -A cond_seen=()
    for rep in "${reps[@]}"; do
        local counts="${OUTDIR}/${rep}/${rep}_genecounts.tsv"
        local qcf="${OUTDIR}/${rep}/${rep}_library_qc.tsv"
        local status path cond raw_cond

        if [[ -s "$counts" ]]; then
            if [[ -f "$qcf" ]]; then status="OK"; else status="OK_NO_QC"; fi
            path="$counts"; n_ok=$((n_ok+1))
        else
            status="MISSING_COUNTS"; path=""; n_bad=$((n_bad+1))
        fi

        if [[ "$CONDITION_FROM" == "dir" ]]; then raw_cond="$rep"; else raw_cond=$(samp_name "$rep"); fi
        [[ -n "$raw_cond" ]] || raw_cond="$rep"
        cond=$(sanitize_level "$raw_cond")
        if [[ "$cond" != "$raw_cond" ]]; then
            echo "  NOTE: condition '${raw_cond}' -> '${cond}' (DESeq2-safe characters only)"
            renamed=$((renamed+1))
        fi
        cond_seen["$cond"]=$(( ${cond_seen["$cond"]:-0} + 1 ))

        # QC values for this replicate (blank when absent)
        local qcvals=""
        if [[ ${#idx[@]} -gt 0 && -f "$qcf" ]]; then
            qcvals=$(awk -F'\t' -v cols="$(IFS=,; echo "${idx[*]}")" '
                NR==2{ n=split(cols,c,","); out="";
                       for(i=1;i<=n;i++){ v=$(c[i]); gsub(/,/,";",v); out=out "," v }
                       print out; exit }' "$qcf")
        elif [[ ${#idx[@]} -gt 0 ]]; then
            local n; for ((n=0; n<${#have[@]}; n++)); do qcvals="${qcvals},"; done
        fi

        printf "%s,%s,,%s,%s%s\n" "$rep" "$cond" "$path" "$status" "$qcvals" >> "$sheet"
    done

    # ---- README ----
    cat > "$readme" <<EOF
samplesheet.csv -- generated $(date) by submit_rnaseq_fastq_to_counts.sh (v2)

BEFORE YOU RUN THE DESeq2 SCRIPT, DO TWO THINGS:

 1. CHECK the 'condition' column.
    It is a GUESS, made by stripping '${REP_KEYWORD}<number>' from the replicate
    name (--condition-from ${CONDITION_FROM}). If two replicates that belong to
    the same group got different condition values -- or two different groups got
    the same value -- fix it here by hand.

 2. FILL IN the '${BATCH_COLUMN}' column. It is intentionally EMPTY.
    Put one label per group. Any labels work, they are just group names:
        1 / 2 / 3        a / b / c        kitA / kitB        Jan / Feb
    Use it for anything that could shift counts independently of the biology:
    library prep kit or lot, prep date, sequencing run or flow cell, operator.
    Then run the DESeq2 script with:  --batch ${BATCH_COLUMN}

    A batch term only works if it is NOT confounded with condition: each batch
    label must contain BOTH conditions. If every Parental is batch 1 and every
    Bottom10 is batch 2, the two effects cannot be separated and the DESeq2
    script will stop and tell you so (this is a real property of the design, not
    a bug). Leave the column empty and omit --batch to fit ~ condition only.

COLUMN MEANINGS
  sample                 replicate name = its output directory name. The DESeq2
                         script finds the count file from this (or from 'path').
  condition              experimental group. DESeq2 factor levels may contain
                         only letters, numbers, underscore and period, so any
                         other character was replaced with '_'.
  ${BATCH_COLUMN}$(printf '%*s' $(( 23 - ${#BATCH_COLUMN} )) '')YOU FILL THIS IN. Empty = no batch term.
  path                   absolute path to <replicate>_genecounts.tsv
  status                 OK             counts present, QC present
                         OK_NO_QC       counts present, QC file missing
                         MISSING_COUNTS job failed or is still running -- DELETE
                                        the row or re-run that sample
  raw_reads              read pairs entering BBDuk (x2 = reads)
  pct_surviving_trim     % of reads kept by trimming
  star_input_reads       read pairs entering STAR
  pct_uniq_mapped        STAR uniquely mapped %
  pct_multimapping       STAR multi-locus %
  pct_too_short          STAR 'unmapped: too short' % -- high means adapter or
                         over-trimming, not biology
  assigned_fragments     featureCounts Assigned (fragments, not reads)
  pct_assigned           Assigned / all fragments featureCounts saw
  features_detected      genes with a non-zero count  <-- "counts that survived"
  features_total         genes in the annotation that were counted over
  strand_verdict         detected: forward | reverse | unstranded | undetermined
  stranded_used          what -s ACTUALLY produced these counts
  strand_matches_detected  no => counted deliberately against the detection
  qc_flag                PASS, or semicolon-joined warnings

NEXT STEP
  sbatch submit_rnaseq_deseq2_dge.sh \\
      --ref ${REF_LC} \\
      --counts-dir ${OUTDIR} \\
      --samplesheet ${OUTDIR}/${SHEET_NAME} \\
      --control <CONTROL_GROUP> --treatment <TREATMENT_GROUP> \\
      --batch ${BATCH_COLUMN} \\
      --outdir <RESULTS_DIR> \\
      --dry-run          # ALWAYS dry-run first: it validates the design matrix
EOF

    echo "Sample sheet -> ${sheet}"
    echo "README       -> ${readme}"
    echo "  ${n_ok} sample(s) with counts, ${n_bad} missing/failed, ${renamed} condition name(s) sanitised."
    echo "  conditions guessed:"
    local c
    for c in "${!cond_seen[@]}"; do echo "    ${c}  (${cond_seen[$c]} replicate(s))"; done
    if [[ ${#cond_seen[@]} -lt 2 ]]; then
        echo "  WARNING: fewer than 2 distinct conditions were guessed. Edit the"
        echo "           'condition' column by hand, or re-run with --rep-keyword"
        echo "           set to the token your file names actually use."
    fi
    for c in "${!cond_seen[@]}"; do
        if [[ ${cond_seen[$c]} -lt 2 ]]; then
            echo "  WARNING: condition '${c}' has only ${cond_seen[$c]} replicate. DESeq2 needs"
            echo "           >= 2 per group on each side of a contrast."
        fi
    done
    if [[ $n_bad -gt 0 ]]; then
        echo "  WARNING: ${n_bad} row(s) are MISSING_COUNTS. Re-run or delete those rows"
        echo "           before handing the sheet to the DESeq2 script."
    fi
    # flag non-PASS QC
    if [[ -s "$qcsum" ]]; then
        awk -F'\t' 'NR==1{for(j=1;j<=NF;j++) if($j=="qc_flag") q=j; next}
                    q && $q!="PASS" && $q!="" {print "  QC: " $1 " -> " $q}' "$qcsum"
    fi
    echo "------------------------------------------------------------"
}

if $COLLECT_ONLY; then
    collect_samplesheet
    exit 0
fi

# ===========================================================================
# SUBMIT MODE
# ===========================================================================
[[ -f "${MAIN_SCRIPT}" ]] || {
    echo "ERROR: cannot find rnaseq_fastq_to_counts.sh next to this submitter."
    echo "       expected here: ${MAIN_SCRIPT}"; exit 1; }

shopt -s nullglob
fastqs=( "${FASTQ_DIR}"/*"${FWD_KEYWORD}"*.fastq.gz "${FASTQ_DIR}"/*"${FWD_KEYWORD}"*.fq.gz )
shopt -u nullglob
if [[ ${#fastqs[@]} -eq 0 ]]; then
    echo "ERROR: no forward-read fastqs matching *${FWD_KEYWORD}*.fastq.gz (or .fq.gz) in ${FASTQ_DIR}"
    echo "       if your files use a different token, set --fwd-keyword accordingly."
    exit 1
fi

# ---------------------------------------------------------------------------
# PRE-FLIGHT: pair every R1, and REFUSE to run if two inputs collapse to the
# same replicate name.
#
# The replicate name IS the output directory. Two inputs that differ only after
# the forward keyword -- or X_R1_001.fastq.gz alongside X_R1_001.fq.gz -- give
# the same name, so two concurrent jobs would write the same BAM and the same
# _genecounts.tsv. The result is a plausible-looking count file assembled from
# interleaved partial writes, with no error anywhere. Caught here instead.
# ---------------------------------------------------------------------------
declare -a JOB_R1=() JOB_R2=() JOB_REP=()
declare -A REP_SOURCE=()
declare -a DUPES=()
skipped=0

for R1path in "${fastqs[@]}"; do
    R1name=$(basename "${R1path}")
    R2name=$(mate_name "${R1name}")

    if [[ "${R2name}" == "${R1name}" ]]; then
        echo "WARNING: could not derive R2 from ${R1name} (fwd keyword '${FWD_KEYWORD}' not found) -- skipping"
        skipped=$((skipped+1)); continue
    fi
    if [[ ! -f "${FASTQ_DIR}/${R2name}" ]]; then
        echo "WARNING: no R2 mate for ${R1name} (expected ${R2name}) -- skipping"
        skipped=$((skipped+1)); continue
    fi

    Rep=$(rep_name "${R1name}")
    if [[ -z "$Rep" ]]; then
        echo "WARNING: empty replicate name derived from ${R1name} -- skipping"
        skipped=$((skipped+1)); continue
    fi
    if [[ -n "${REP_SOURCE[$Rep]:-}" ]]; then
        DUPES+=( "${Rep}  <-  ${REP_SOURCE[$Rep]}  AND  ${R1name}" )
        continue
    fi
    REP_SOURCE[$Rep]="${R1name}"
    JOB_R1+=( "${R1name}" ); JOB_R2+=( "${R2name}" ); JOB_REP+=( "${Rep}" )
done

if [[ ${#DUPES[@]} -gt 0 ]]; then
    echo
    echo "ERROR: these input files collapse to the SAME replicate name:"
    printf '         %s\n' "${DUPES[@]}"
    echo
    echo "       The replicate name is also the output directory name, so those jobs"
    echo "       would overwrite each other's BAM and count file simultaneously and"
    echo "       produce a corrupt count matrix with no error message."
    echo "       Fix: rename the fastqs, or choose a --fwd-keyword that keeps them"
    echo "       distinct (the keyword and everything after it is stripped)."
    exit 1
fi
[[ ${#JOB_REP[@]} -gt 0 ]] || { echo "ERROR: no usable R1/R2 pairs found."; exit 1; }

# ---------------------------------------------------------------------------
# Plan
# ---------------------------------------------------------------------------
echo "------------------------------------------------------------"
echo "Reference   : ${REF_LC}"
echo "Fastq dir   : ${FASTQ_DIR}"
echo "Output      : ${OUTDIR}"
echo "Pipeline    : ${MAIN_SCRIPT}"
echo "Pairs found : ${#JOB_REP[@]} (skipped ${skipped})"
if [[ ${#PASSTHRU[@]} -gt 0 ]]; then
    echo "Forwarding  : ${PASSTHRU[*]}"
else
    echo "Forwarding  : <pipeline defaults>  (trimq=10 minlen=36 STAR 0.66/0.66 -s 0)"
    echo "              NOTE: -s 0 (unstranded). If your kit is directional, pass"
    echo "              --stranded reverse (or forward). Run one sample first and"
    echo "              read strand_verdict from its _library_qc.tsv if unsure."
fi
echo "------------------------------------------------------------"

MANIFEST="${OUTDIR}/samples_manifest.tsv"
if $DRY_RUN; then
    echo "DRY RUN -- would submit:"
    for i in "${!JOB_REP[@]}"; do
        cond_guess=$( [[ "$CONDITION_FROM" == "dir" ]] && echo "${JOB_REP[$i]}" || samp_name "${JOB_REP[$i]}" )
        printf "  %-28s R1=%-34s condition guess=%s\n" \
               "${JOB_REP[$i]}" "${JOB_R1[$i]}" "$(sanitize_level "${cond_guess:-${JOB_REP[$i]}}")"
    done
    $MAKE_SHEET && echo "  + 1 collector job -> ${OUTDIR}/${SHEET_NAME}"
    echo "Nothing submitted (--dry-run)."
    exit 0
fi

printf "replicate\tsample\tR1\tR2\n" > "${MANIFEST}"

# ---------------------------------------------------------------------------
# Submit
# ---------------------------------------------------------------------------
declare -a JOBIDS=()
submitted=0
for i in "${!JOB_REP[@]}"; do
    Rep="${JOB_REP[$i]}"; R1name="${JOB_R1[$i]}"; R2name="${JOB_R2[$i]}"
    cond_guess=$( [[ "$CONDITION_FROM" == "dir" ]] && echo "$Rep" || samp_name "$Rep" )
    cond_guess=$(sanitize_level "${cond_guess:-$Rep}")

    printf "%s\t%s\t%s\t%s\n" "${Rep}" "${cond_guess}" "${R1name}" "${R2name}" >> "${MANIFEST}"

    jid=$(sbatch --parsable -J "rna_${Rep}" "${MAIN_SCRIPT}" \
            --ref "${REF_LC}" \
            --fastq-dir "${FASTQ_DIR}" \
            --outdir "${OUTDIR}" \
            --fwd-keyword "${FWD_KEYWORD}" \
            --r1 "${R1name}" --r2 "${R2name}" \
            ${PASSTHRU[@]+"${PASSTHRU[@]}"})
    JOBIDS+=( "$jid" )
    echo "Submitted: ${Rep}  (job ${jid}, condition guess=${cond_guess})"
    submitted=$((submitted + 1))
done

# ---------------------------------------------------------------------------
# Collector job -- waits for ALL sample jobs (afterany, so a failure still
# produces a sheet with MISSING_COUNTS rows rather than no sheet at all).
# ---------------------------------------------------------------------------
COLLECT_JOB=""
if $MAKE_SHEET && [[ ${#JOBIDS[@]} -gt 0 ]]; then
    dep=$(IFS=:; echo "${JOBIDS[*]}")
    mkdir -p "${OUTDIR}/logs"
    COLLECT_JOB=$(sbatch --parsable \
        --dependency="afterany:${dep}" \
        --kill-on-invalid-dep=yes \
        -J "rna_samplesheet" \
        -A "${SLURM_ACCOUNT}" -p "${SLURM_PARTITION}" -q "${SLURM_QOS}" \
        -N 1 -n 1 -c 1 --mem=4G -t "${COLLECTOR_TIME}" \
        --output="${OUTDIR}/logs/samplesheet-%j.out" \
        --error="${OUTDIR}/logs/samplesheet-%j.err" \
        --wrap "bash '${SCRIPT_SELF}' --collect-only --outdir '${OUTDIR}' \
                --ref '${REF_LC}' --rep-keyword '${REP_KEYWORD}' \
                --condition-from '${CONDITION_FROM}' \
                --batch-column '${BATCH_COLUMN}' --samplesheet-name '${SHEET_NAME}'")
fi

echo "------------------------------------------------------------"
echo "Submitted ${submitted} sample job(s)."
echo "Reference : ${REF_LC}"
echo "Output    : ${OUTDIR}"
echo "Manifest  : ${MANIFEST}"
if [[ -n "$COLLECT_JOB" ]]; then
    echo "Collector : job ${COLLECT_JOB} (runs after all samples finish)"
    echo "            -> ${OUTDIR}/${SHEET_NAME}"
    echo "            -> ${OUTDIR}/sample_qc_summary.tsv"
    echo "            -> ${OUTDIR}/SAMPLESHEET_README.txt"
    echo
    echo "WHEN IT FINISHES: open ${SHEET_NAME}, check 'condition',"
    echo "                  fill in the empty '${BATCH_COLUMN}' column, then run the"
    echo "                  DESeq2 script with --samplesheet and --batch ${BATCH_COLUMN}."
elif ! $MAKE_SHEET; then
    echo "Collector : skipped (--no-samplesheet). Build it later with:"
    echo "            ${SCRIPT_SELF} --collect-only --outdir '${OUTDIR}'"
fi
echo "Track with:  squeue -u \$USER"
echo "------------------------------------------------------------"