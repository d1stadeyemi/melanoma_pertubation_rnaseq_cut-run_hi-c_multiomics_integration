#!/usr/bin/env bash
# =============================================================================
# kmet_spikein_count.sh
#
# Count SNAP-CUTANA K-MetStat Panel (EpiCypher 19-1002) spike-in barcode reads
# directly from paired-end CUT&RUN FASTQ files, and report per-PTM recovery,
# antibody specificity, and spike-in fraction of the library.
#
# This is a from-scratch reimplementation of the counting logic described in the
# EpiCypher SNAP-CUTANA Spike-in User Guide v1.2 (Section 3.4), written for
# batch use on an HPC cluster instead of one hand-edited loop per reaction.
#
# Differences vs. the vendor's shell script (all deliberate, all documented):
#   1. Batch: reads a manifest, loops over every sample, no manual editing.
#   2. Works directly on .fastq.gz (no decompression to disk).
#   3. Single streaming pass per sample instead of 32 separate greps per FASTQ.
#   4. Searches BOTH the barcode and its reverse complement in BOTH mates.
#      The vendor script greps the forward barcode only. If your library was
#      sequenced such that a barcode lands on the minus strand of a mate, the
#      vendor script silently misses it. Orientation counts are reported
#      separately so you can see whether this matters for your data.
#   5. Reports a PAIR-level count (a fragment is counted once, no matter how
#      many mates/orientations carry the barcode) in addition to the vendor's
#      R1 and R2 counts. Pair counts are the correct unit for normalization;
#      R1/R2 counts are kept so you can paste them straight into EpiCypher's
#      Panel Analysis .xlsx if you want the vendor heatmap.
#
# Author: (McCord Lab multi-omics project)
# =============================================================================

set -o errexit
set -o nounset
set -o pipefail

VERSION="1.1.0"

# ----------------------------------------------------------------------------- 
# Usage
# -----------------------------------------------------------------------------
usage() {
cat <<'EOF'
kmet_spikein_count.sh  --  count SNAP-CUTANA K-MetStat spike-in barcodes from FASTQ

USAGE
  Batch mode (recommended):
    kmet_spikein_count.sh --manifest samples.tsv --outdir OUTDIR [options]

  Single-sample mode:
    kmet_spikein_count.sh --sample NAME --r1 R1.fastq.gz --r2 R2.fastq.gz \
                          --outdir OUTDIR [--target H3K27me3] [--bam FILE] [options]

REQUIRED
  --manifest FILE     Tab-separated, with a header line. Required columns:
                        sample   r1   r2
                      Optional columns (any order, matched by header name):
                        target   the on-target PTM for this reaction. Must be one of:
                                 Unmodified, H3K4me1, H3K4me2, H3K4me3,
                                 H3K9me1,  H3K9me2,  H3K9me3,
                                 H3K27me1, H3K27me2, H3K27me3,
                                 H3K36me1, H3K36me2, H3K36me3,
                                 H4K20me1, H4K20me2, H4K20me3,
                                 IgG      (negative control; no on-target PTM)
                                 NA       (target not in the panel, e.g. H3K27ac)
                        bam      path to the final filtered BAM for this sample.
                                 If given, the script uses samtools to count the
                                 aligned fragments actually used downstream, and
                                 reports spike-in % against that denominator.
                      Lines starting with '#' are ignored.
  --outdir DIR        Output directory (created if absent).

  ...or, instead of --manifest:
  --sample NAME       Sample name (single-sample mode).
  --r1 FILE           R1 FASTQ (.fastq or .fastq.gz).
  --r2 FILE           R2 FASTQ (.fastq or .fastq.gz).

OPTIONS
  --target PTM        On-target PTM in single-sample mode. Allowed values as above.
                      Default: NA
  --bam FILE          Final filtered BAM in single-sample mode. Default: none.
  --mapq INT          MAPQ threshold used when counting aligned fragments from
                      --bam. Set this to whatever your csaw/deepTools filter uses
                      so the denominators match. Default: 20
  --samtools-flags S  Extra flags passed to `samtools view -c`. Default: '-f 66 -F 3852'
                      (-f 66  = paired + first-in-pair, i.e. count fragments once
                       -F 3852 = drop unmapped, mate-unmapped, secondary,
                                 supplementary, QC-fail, duplicate)
  --threads INT       Threads for decompression / samtools. Default: 4
  --keep-tmp          Keep per-sample temporary files.
  --version           Print version and exit.
  -h, --help          Print this help and exit.

OUTPUTS  (all written under --outdir)
  per_sample/<sample>.barcode_counts.tsv
        One row per barcode (32 rows): ptm, replicate (A/B), barcode,
        r1_fwd, r1_rev, r2_fwd, r2_rev, r1_total, r2_total, pair_count
        -> paste r1_total into Column C and r2_total into Column D of
           EpiCypher's Panel Analysis .xlsx if you want their heatmap.

  kmet_barcode_matrix.tsv
        PTM x sample matrix of PAIR counts (16 PTMs, A+B barcodes summed).
        This is the input for kmet_spikein_normfactors.R.

  kmet_specificity.tsv
        PTM x sample matrix of recovery as % of the on-target PTM.
        EpiCypher's pass criterion: every off-target PTM < 20%.

  kmet_summary.tsv
        Per-sample QC: total read pairs, spike-in pairs, % spike-in of total,
        % spike-in of aligned (if --bam), on-target PTM, % on-target of all
        spike-in reads, worst off-target %, and a PASS/FLAG call.

  MANIFEST_kmet.txt
        Provenance: version, command line, date, host, tool versions, barcode
        table checksum.

NOTES
  - EpiCypher requires >= 50 bp paired-end reads to read the 22 bp barcodes.
  - Run this on ADAPTER-TRIMMED reads, not raw and not quality-trimmed.
    See the project README / discussion for why.
  - Expected spike-in fraction: ~0.5-5% of total reads. High-abundance targets
    (H3K27me3) run 0.1-1%; IgG and low-abundance targets run 10-20%.
    Thousands of barcode reads are needed for the panel to be usable.
EOF
}

# -----------------------------------------------------------------------------
# Defaults
# -----------------------------------------------------------------------------
MANIFEST=""
MANIFEST_CSV=""
ROOT=""
TRIMDIR="trimmed"
FACTOR_FILTER=""
CSV_TARGET="H3K27me3"
R1_GLOB='*_1.paired.fastq.gz'
R2_GLOB='*_2.paired.fastq.gz'
OUTDIR=""
S_SAMPLE=""
S_R1=""
S_R2=""
S_TARGET="NA"
S_BAM=""
MAPQ=20
SAMFLAGS="-f 66 -F 3852"
THREADS=4
KEEP_TMP=0

# -----------------------------------------------------------------------------
# Parse args
# -----------------------------------------------------------------------------
if [[ $# -eq 0 ]]; then usage; exit 1; fi
while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest)        MANIFEST="$2"; shift 2 ;;
    --manifest-csv)    MANIFEST_CSV="$2"; shift 2 ;;
    --root)            ROOT="$2"; shift 2 ;;
    --trimdir)         TRIMDIR="$2"; shift 2 ;;
    --factor)          FACTOR_FILTER="$2"; shift 2 ;;
    --csv-target)      CSV_TARGET="$2"; shift 2 ;;
    --r1-glob)         R1_GLOB="$2"; shift 2 ;;
    --r2-glob)         R2_GLOB="$2"; shift 2 ;;
    --outdir)          OUTDIR="$2"; shift 2 ;;
    --sample)          S_SAMPLE="$2"; shift 2 ;;
    --r1)              S_R1="$2"; shift 2 ;;
    --r2)              S_R2="$2"; shift 2 ;;
    --target)          S_TARGET="$2"; shift 2 ;;
    --bam)             S_BAM="$2"; shift 2 ;;
    --mapq)            MAPQ="$2"; shift 2 ;;
    --samtools-flags)  SAMFLAGS="$2"; shift 2 ;;
    --threads)         THREADS="$2"; shift 2 ;;
    --keep-tmp)        KEEP_TMP=1; shift ;;
    --version)         echo "kmet_spikein_count.sh v${VERSION}"; exit 0 ;;
    -h|--help)         usage; exit 0 ;;
    *) echo "ERROR: unknown argument '$1'" >&2; echo "Run with --help." >&2; exit 1 ;;
  esac
done

[[ -z "$OUTDIR" ]] && { echo "ERROR: --outdir is required." >&2; exit 1; }

if [[ -n "$MANIFEST_CSV" ]]; then
  [[ -n "$MANIFEST" ]] && { echo "ERROR: use --manifest-csv OR --manifest, not both." >&2; exit 1; }
  [[ -n "$S_SAMPLE" || -n "$S_R1" ]] && { echo "ERROR: use --manifest-csv OR single-sample flags, not both." >&2; exit 1; }
  [[ -r "$MANIFEST_CSV" ]] || { echo "ERROR: cannot read '$MANIFEST_CSV'." >&2; exit 1; }
  [[ -z "$ROOT" ]] && { echo "ERROR: --manifest-csv also needs --root (the raw_dirs directory)." >&2; exit 1; }
  [[ -d "$ROOT" ]] || { echo "ERROR: --root '$ROOT' is not a directory." >&2; exit 1; }
  [[ -z "$FACTOR_FILTER" ]] && { echo "ERROR: --manifest-csv also needs --factor (e.g. --factor H3K27me3)." >&2; exit 1; }
  command -v python3 >/dev/null 2>&1 || { echo "ERROR: --manifest-csv mode needs python3 on PATH." >&2; exit 1; }
elif [[ -n "$MANIFEST" ]]; then
  [[ -n "$S_SAMPLE" || -n "$S_R1" ]] && { echo "ERROR: use --manifest OR single-sample flags, not both." >&2; exit 1; }
  [[ -r "$MANIFEST" ]] || { echo "ERROR: cannot read manifest '$MANIFEST'." >&2; exit 1; }
else
  [[ -z "$S_SAMPLE" || -z "$S_R1" || -z "$S_R2" ]] && \
    { echo "ERROR: single-sample mode needs --sample, --r1 and --r2." >&2; exit 1; }
fi

# -----------------------------------------------------------------------------
# Allowed PTM names
# -----------------------------------------------------------------------------
PTM_ORDER=(Unmodified H3K4me1 H3K4me2 H3K4me3 H3K9me1 H3K9me2 H3K9me3 \
           H3K27me1 H3K27me2 H3K27me3 H3K36me1 H3K36me2 H3K36me3 \
           H4K20me1 H4K20me2 H4K20me3)

validate_target() {
  local t="$1"
  [[ "$t" == "NA" || "$t" == "IgG" ]] && return 0
  local p
  for p in "${PTM_ORDER[@]}"; do [[ "$t" == "$p" ]] && return 0; done
  echo "ERROR: --target '$t' is not a K-MetStat PTM." >&2
  echo "       Allowed: ${PTM_ORDER[*]} IgG NA" >&2
  exit 1
}

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
mkdir -p "$OUTDIR/per_sample"
TMPROOT="$(mktemp -d "${TMPDIR:-/tmp}/kmet.XXXXXX")"
cleanup() { [[ "$KEEP_TMP" -eq 0 ]] && rm -rf "$TMPROOT"; }
trap cleanup EXIT

# ---- barcode table (EpiCypher SNAP-CUTANA User Guide v1.2, Table 3) ----------
BCFILE="$OUTDIR/kmet_barcodes.tsv"
cat > "$BCFILE" <<'EOF'
ptm	replicate	barcode
Unmodified	A	TTCGCGCGTAACGACGTACCGT
Unmodified	B	CGCGATACGACCGCGTTACGCG
H3K4me1	A	CGACGTTAACGCGTTTCGTACG
H3K4me1	B	CGCGACTATCGCGCGTAACGCG
H3K4me2	A	CCGTACGTCGTGTCGAACGACG
H3K4me2	B	CGATACGCGTTGGTACGCGTAA
H3K4me3	A	TAGTTCGCGACACCGTTCGTCG
H3K4me3	B	TCGACGCGTAAACGGTACGTCG
H3K9me1	A	TTATCGCGTCGCGACGGACGTA
H3K9me1	B	CGATCGTACGATAGCGTACCGA
H3K9me2	A	CGCATATCGCGTCGTACGACCG
H3K9me2	B	ACGTTCGACCGCGGTCGTACGA
H3K9me3	A	ACGATTCGACGATCGTCGACGA
H3K9me3	B	CGATAGTCGCGTCGCACGATCG
H3K27me1	A	CGCCGATTACGTGTCGCGCGTA
H3K27me1	B	ATCGTACCGCGCGTATCGGTCG
H3K27me2	A	CGTTCGAACGTTCGTCGACGAT
H3K27me2	B	TCGCGATTACGATGTCGCGCGA
H3K27me3	A	ACGCGAATCGTCGACGCGTATA
H3K27me3	B	CGCGATATCACTCGACGCGATA
H3K36me1	A	CGCGAAATTCGTATACGCGTCG
H3K36me1	B	CGCGATCGGTATCGGTACGCGC
H3K36me2	A	GTGATATCGCGTTAACGTCGCG
H3K36me2	B	TATCGCGCGAAACGACCGTTCG
H3K36me3	A	CCGCGCGTAATGCGCGACGTTA
H3K36me3	B	CCGCGATACGACTCGTTCGTCG
H4K20me1	A	GTCGCGAACTATCGTCGATTCG
H4K20me1	B	CCGCGCGTATAGTCCGAGCGTA
H4K20me2	A	CGATACGCCGATCGATCGTCGG
H4K20me2	B	CCGCGCGATAAGACGCGTAACG
H4K20me3	A	CGATTCGACGGTCGCGACCGTA
H4K20me3	B	TTTCGACGCGTCGATTCGGCGA
EOF

# ---- expand to forward + revcomp pattern list --------------------------------
PATFILE="$TMPROOT/patterns.tsv"     # pattern <TAB> ptm <TAB> rep <TAB> orient
awk 'BEGIN{FS=OFS="\t"}
function revcomp(s,   i,c,o){
  o=""
  for(i=length(s); i>=1; i--){
    c=substr(s,i,1)
    if(c=="A")o=o"T"; else if(c=="T")o=o"A";
    else if(c=="G")o=o"C"; else if(c=="C")o=o"G"; else o=o"N"
  }
  return o
}
NR==1{next}
{
  print $3, $1, $2, "fwd"
  print revcomp($3), $1, $2, "rev"
}' "$BCFILE" > "$PATFILE"

cut -f1 "$PATFILE" > "$TMPROOT/patterns.txt"

NPAT=$(wc -l < "$TMPROOT/patterns.txt")
if [[ "$NPAT" -ne 64 ]]; then
  echo "ERROR: expected 64 patterns (32 barcodes x 2 orientations), got $NPAT." >&2
  exit 1
fi
# sanity: all patterns unique
if [[ "$(sort -u "$TMPROOT/patterns.txt" | wc -l)" -ne 64 ]]; then
  echo "ERROR: barcode patterns are not unique after revcomp expansion." >&2
  exit 1
fi

# -----------------------------------------------------------------------------
# awk program that assigns prefiltered read pairs to barcodes
# -----------------------------------------------------------------------------
ASSIGN="$TMPROOT/assign.awk"
cat > "$ASSIGN" <<'AWKEOF'
BEGIN{
  FS=OFS="\t"
  # load pattern table
  np=0
  while((getline line < patfile) > 0){
    split(line, f, "\t")
    np++
    pat[np]=f[1]; ppmt[np]=f[2]; prep[np]=f[3]; pori[np]=f[4]
    key=f[2] "\t" f[3]
    seen[key]=1
  }
  close(patfile)
}
{
  s1=$1; s2=$2
  nhit=0
  delete hitkey
  for(i=1;i<=np;i++){
    h1 = index(s1, pat[i]) > 0
    h2 = index(s2, pat[i]) > 0
    if(!h1 && !h2) continue
    key = ppmt[i] "\t" prep[i]
    if(h1){
      if(pori[i]=="fwd") r1f[key]++; else r1r[key]++
      r1any[key]=1
    }
    if(h2){
      if(pori[i]=="fwd") r2f[key]++; else r2r[key]++
      r2any[key]=1
    }
    if(!(key in hitkey)){ hitkey[key]=1; nhit++ }
  }
  if(nhit==0) next
  # R1/R2 presence counted once per mate per barcode (EpiCypher columns C and D)
  for(k in r1any){ r1tot[k]++; delete r1any[k] }
  for(k in r2any){ r2tot[k]++; delete r2any[k] }
  # pair-level: fragment counted once per barcode
  for(k in hitkey) pairs[k]++
  if(nhit>1) ambig++
  total_spike_pairs++
}
END{
  for(k in seen){
    split(k, f, "\t")
    printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
      f[1], f[2],
      (k in r1f ? r1f[k] : 0), (k in r1r ? r1r[k] : 0),
      (k in r2f ? r2f[k] : 0), (k in r2r ? r2r[k] : 0),
      (k in r1tot ? r1tot[k] : 0), (k in r2tot ? r2tot[k] : 0),
      (k in pairs ? pairs[k] : 0) > countfile
  }
  close(countfile)
  printf "%d\t%d\n", (total_spike_pairs+0), (ambig+0) > statfile
  close(statfile)
}
AWKEOF

# -----------------------------------------------------------------------------
# helpers
# -----------------------------------------------------------------------------
opencat() {
  # stream a fastq(.gz) to stdout
  case "$1" in
    *.gz) if command -v pigz >/dev/null 2>&1; then pigz -dc -p "$THREADS" "$1"; else gzip -dc "$1"; fi ;;
    *)    cat "$1" ;;
  esac
}

count_aligned_fragments() {
  # $1 = bam
  local bam="$1"
  command -v samtools >/dev/null 2>&1 || { echo "NA"; return; }
  samtools view -c -@ "$THREADS" -q "$MAPQ" $SAMFLAGS "$bam" 2>/dev/null || echo "NA"
}

# -----------------------------------------------------------------------------
# Core: process one sample
# -----------------------------------------------------------------------------
process_sample() {
  local sample="$1" r1="$2" r2="$3" target="$4" bam="$5"

  [[ -r "$r1" ]] || { echo "ERROR [$sample]: cannot read R1 '$r1'." >&2; exit 1; }
  [[ -r "$r2" ]] || { echo "ERROR [$sample]: cannot read R2 '$r2'." >&2; exit 1; }
  validate_target "$target"

  echo "[$(date +%H:%M:%S)] $sample  (target=$target)" >&2

  local sd="$TMPROOT/$sample"
  mkdir -p "$sd"

  # One streaming pass:
  #   paste the two sequence-line streams
  #   -> pass-through awk that also records the total read-pair count
  #   -> grep -F -f (multi-pattern, drops ~99% of read pairs)
  #   -> awk assigns surviving pairs to barcodes
  paste <(opencat "$r1" | awk 'NR%4==2') \
        <(opencat "$r2" | awk 'NR%4==2') \
  | awk -v cf="$sd/total_pairs.txt" '{ print } END { print NR+0 > cf; close(cf) }' \
  | { grep -F -f "$TMPROOT/patterns.txt" || true; } \
  | awk -v patfile="$PATFILE" \
        -v countfile="$sd/raw_counts.tsv" \
        -v statfile="$sd/spike_stats.tsv" \
        -f "$ASSIGN"

  # grep may emit nothing; make sure the awk END block still ran
  if [[ ! -s "$sd/raw_counts.tsv" ]]; then
    awk 'BEGIN{FS=OFS="\t"} NR>1{print $1,$2,0,0,0,0,0,0,0}' "$BCFILE" > "$sd/raw_counts.tsv"
    printf "0\t0\n" > "$sd/spike_stats.tsv"
  fi

  local total_pairs spike_pairs ambig
  total_pairs=$(cat "$sd/total_pairs.txt")
  spike_pairs=$(cut -f1 "$sd/spike_stats.tsv")
  ambig=$(cut -f2 "$sd/spike_stats.tsv")

  if [[ "$total_pairs" -eq 0 ]]; then
    echo "ERROR [$sample]: R1 has zero reads." >&2; exit 1
  fi

  # ---- per-sample barcode table, in panel order ----
  local out="$OUTDIR/per_sample/${sample}.barcode_counts.tsv"
  {
    printf "sample\tptm\treplicate\tbarcode\tr1_fwd\tr1_rev\tr2_fwd\tr2_rev\tr1_total\tr2_total\tpair_count\n"
    awk -v S="$sample" 'BEGIN{FS=OFS="\t"}
      NR==FNR && FNR>1 { bc[$1"\t"$2]=$3; next }
      NR!=FNR { k=$1"\t"$2; print S,$1,$2,(k in bc?bc[k]:"NA"),$3,$4,$5,$6,$7,$8,$9 }' \
      "$BCFILE" "$sd/raw_counts.tsv" \
    | awk -v order="$(printf '%s ' "${PTM_ORDER[@]}")" 'BEGIN{FS=OFS="\t"; n=split(order,o," "); for(i=1;i<=n;i++) rank[o[i]]=i}
      { print rank[$2], $0 }' \
    | sort -k1,1n -k4,4 \
    | cut -f2-
  } > "$out"

  # ---- collapse A+B to per-PTM pair counts ----
  awk -v S="$sample" 'BEGIN{FS=OFS="\t"} NR>1 { p[$2]+=$11 } END{ for(k in p) print S,k,p[k] }' "$out" \
    > "$sd/ptm_pairs.tsv"

  # ---- QC metrics ----
  local pct_total on_target_reads pct_on_target worst_off worst_off_ptm aligned pct_aligned
  pct_total=$(awk -v a="$spike_pairs" -v b="$total_pairs" 'BEGIN{printf "%.4f", (b>0? 100*a/b : 0)}')

  aligned="NA"; pct_aligned="NA"
  if [[ -n "$bam" && "$bam" != "NA" && -r "$bam" ]]; then
    aligned=$(count_aligned_fragments "$bam")
    if [[ "$aligned" != "NA" && "$aligned" -gt 0 ]]; then
      pct_aligned=$(awk -v a="$spike_pairs" -v b="$aligned" 'BEGIN{printf "%.4f", 100*a/b}')
    fi
  fi

  on_target_reads="NA"; pct_on_target="NA"; worst_off="NA"; worst_off_ptm="NA"
  if [[ "$target" != "NA" && "$target" != "IgG" ]]; then
    on_target_reads=$(awk -v t="$target" 'BEGIN{FS="\t"} $2==t{print $3}' "$sd/ptm_pairs.tsv")
    on_target_reads=${on_target_reads:-0}
    if [[ "$spike_pairs" -gt 0 ]]; then
      pct_on_target=$(awk -v a="$on_target_reads" -v b="$spike_pairs" 'BEGIN{printf "%.2f", 100*a/b}')
    fi
    if [[ "$on_target_reads" -gt 0 ]]; then
      read -r worst_off worst_off_ptm < <(awk -v t="$target" -v ot="$on_target_reads" \
        'BEGIN{FS="\t"; m=-1; mp="NA"} $2!=t { v=100*$3/ot; if(v>m){m=v; mp=$2} } END{printf "%.2f %s\n", (m<0?0:m), mp}' \
        "$sd/ptm_pairs.tsv")
    fi
  fi

  # ---- PASS / FLAG ----
  local flags=""
  awk -v p="$pct_total" 'BEGIN{exit !(p<0.1)}' && flags="${flags}LOW_SPIKEIN_FRACTION;"
  awk -v p="$pct_total" 'BEGIN{exit !(p>20)}'  && flags="${flags}HIGH_SPIKEIN_FRACTION;"
  [[ "$spike_pairs" -lt 1000 ]] && flags="${flags}FEW_BARCODE_READS(<1000);"
  if [[ "$worst_off" != "NA" ]]; then
    awk -v w="$worst_off" 'BEGIN{exit !(w>=20)}' && flags="${flags}OFF_TARGET>=20%;"
  fi
  [[ -z "$flags" ]] && flags="PASS"

  printf "%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\n" \
    "$sample" "$target" "$total_pairs" "$spike_pairs" "$pct_total" \
    "$aligned" "$pct_aligned" "$on_target_reads" "$pct_on_target" \
    "$worst_off" "$worst_off_ptm" "$ambig" "$flags" \
    >> "$TMPROOT/summary_body.tsv"

  cat "$sd/ptm_pairs.tsv" >> "$TMPROOT/all_ptm_pairs.tsv"
}


# -----------------------------------------------------------------------------
# --manifest-csv mode: resolve the project manifest.csv into an internal TSV.
#
# Naming contract (matches generate_manifest.R):
#   * The USAGE NAME is always the manifest's `folder` column. It is the sample
#     identity that carries into the sample sheet and every downstream analysis.
#   * The REAL FILE is whatever is actually on disk under that folder. It is
#     DISCOVERED by glob, never constructed from the folder name. Real example
#     from this project: folder 'B10-UNC-H3K27Me3-R1' contains
#     'trimmed/B10-UNC-H3K27Me3_1.paired.fastq.gz' -- no '-R1'. Assuming
#     folder == filename silently breaks on rows like that.
#   * bam_path is taken verbatim from the manifest (already disk-discovered).
# -----------------------------------------------------------------------------
resolve_csv_manifest() {
  python3 - "$MANIFEST_CSV" "$ROOT" "$TRIMDIR" "$FACTOR_FILTER" "$CSV_TARGET" \
             "$R1_GLOB" "$R2_GLOB" <<'PYEOF'
import csv, glob, os, sys

csv_path, root, trimdir, factors, target, r1g, r2g = sys.argv[1:8]
want = {f.strip().lower() for f in factors.split(",") if f.strip()}

def die(msg):
    sys.stderr.write("ERROR: " + msg + "\n"); sys.exit(1)

with open(csv_path, newline="") as fh:
    rows = list(csv.DictReader(fh))
if not rows:
    die("manifest CSV '%s' has no rows." % csv_path)
for col in ("folder", "Factor", "bam_path"):
    if col not in rows[0]:
        die("manifest CSV is missing required column '%s'. Found: %s"
            % (col, ", ".join(rows[0].keys())))

seen_factors = sorted({(r.get("Factor") or "").strip() for r in rows})
sel = [r for r in rows if (r.get("Factor") or "").strip().lower() in want]
if not sel:
    die("no rows matched --factor '%s'.\n       Factor values present in %s:\n         %s"
        % (factors, csv_path, "\n         ".join(seen_factors)))

out, problems = [], []
for r in sel:
    folder = (r.get("folder") or "").strip()
    fdir = os.path.join(root, folder, trimdir)
    if not os.path.isdir(fdir):
        problems.append("%s: no '%s/' directory under %s" % (folder, trimdir, os.path.join(root, folder)))
        continue
    r1 = sorted(glob.glob(os.path.join(fdir, r1g)))
    r2 = sorted(glob.glob(os.path.join(fdir, r2g)))
    if len(r1) != 1 or len(r2) != 1:
        problems.append("%s: expected exactly 1 R1 and 1 R2 matching '%s'/'%s' in %s, found %d/%d%s"
                        % (folder, r1g, r2g, fdir, len(r1), len(r2),
                           (" -> " + ", ".join(os.path.basename(x) for x in r1 + r2)) if (r1 or r2) else ""))
        continue
    bam = (r.get("bam_path") or "").strip()
    bam = os.path.normpath(bam) if bam else ""
    if bam and not os.path.exists(bam):
        problems.append("%s: bam_path does not exist: %s" % (folder, bam))
        bam = ""
    fq_base = os.path.basename(r1[0]).replace("_1.paired.fastq.gz", "")
    if fq_base != folder:
        sys.stderr.write("[INFO] %s: fastq basename is '%s' (differs from folder; using folder as the name)\n"
                         % (folder, fq_base))
    out.append((folder, r1[0], r2[0], target, bam))

if problems:
    sys.stderr.write("\n[ACTION NEEDED] %d row(s) could not be resolved:\n" % len(problems))
    for p in problems:
        sys.stderr.write("   " + p + "\n")
    sys.stderr.write("\nRefusing to run on a partial sample set -- a silently missing sample\n"
                     "becomes a silently missing normalization factor. Fix or filter, then rerun.\n")
    sys.exit(1)

print("sample\tr1\tr2\ttarget\tbam")
for row in out:
    print("\t".join(row))
sys.stderr.write("[INFO] resolved %d sample(s) from %s (Factor in {%s}, trimdir=%s)\n"
                 % (len(out), csv_path, ", ".join(sorted(want)), trimdir))
PYEOF
}

# -----------------------------------------------------------------------------
# Dispatch
# -----------------------------------------------------------------------------
: > "$TMPROOT/summary_body.tsv"
: > "$TMPROOT/all_ptm_pairs.tsv"
SAMPLE_LIST=()

if [[ -n "$MANIFEST_CSV" ]]; then
  MANIFEST="$TMPROOT/resolved_manifest.tsv"
  resolve_csv_manifest > "$MANIFEST"
  cp "$MANIFEST" "$OUTDIR/resolved_manifest.tsv"
  echo "[INFO] resolved manifest kept at $OUTDIR/resolved_manifest.tsv" >&2
fi

if [[ -n "$MANIFEST" ]]; then
  # locate columns by header name
  IFS=$'\t' read -r -a HDR < <(grep -v '^#' "$MANIFEST" | head -n1)
  ci_sample=-1; ci_r1=-1; ci_r2=-1; ci_target=-1; ci_bam=-1
  for i in "${!HDR[@]}"; do
    case "$(echo "${HDR[$i]}" | tr '[:upper:]' '[:lower:]' | tr -d '[:space:]')" in
      sample) ci_sample=$i ;;
      r1|fastq_1|fastq1) ci_r1=$i ;;
      r2|fastq_2|fastq2) ci_r2=$i ;;
      target|ptm) ci_target=$i ;;
      bam) ci_bam=$i ;;
    esac
  done
  [[ $ci_sample -lt 0 || $ci_r1 -lt 0 || $ci_r2 -lt 0 ]] && \
    { echo "ERROR: manifest must have columns 'sample', 'r1', 'r2' in the header." >&2; exit 1; }

  while IFS=$'\t' read -r -a F; do
    [[ -z "${F[0]:-}" ]] && continue
    [[ "${F[0]}" == \#* ]] && continue
    s="${F[$ci_sample]}"; r1="${F[$ci_r1]}"; r2="${F[$ci_r2]}"
    t="NA"; [[ $ci_target -ge 0 ]] && t="${F[$ci_target]:-NA}"
    b="";   [[ $ci_bam    -ge 0 ]] && b="${F[$ci_bam]:-}"
    [[ -z "$t" ]] && t="NA"
    process_sample "$s" "$r1" "$r2" "$t" "$b"
    SAMPLE_LIST+=("$s")
  done < <(grep -v '^#' "$MANIFEST" | tail -n +2)
else
  process_sample "$S_SAMPLE" "$S_R1" "$S_R2" "$S_TARGET" "$S_BAM"
  SAMPLE_LIST+=("$S_SAMPLE")
fi

# -----------------------------------------------------------------------------
# Aggregate
# -----------------------------------------------------------------------------
{
  printf "sample\ttarget\ttotal_read_pairs\tspikein_pairs\tpct_spikein_of_total\t"
  printf "aligned_fragments\tpct_spikein_of_aligned\ton_target_pairs\tpct_on_target_of_spikein\t"
  printf "worst_offtarget_pct\tworst_offtarget_ptm\tambiguous_pairs\tqc\n"
  cat "$TMPROOT/summary_body.tsv"
} > "$OUTDIR/kmet_summary.tsv"

# PTM x sample pair-count matrix
awk -v order="$(printf '%s ' "${PTM_ORDER[@]}")" -v samples="$(printf '%s ' "${SAMPLE_LIST[@]}")" '
BEGIN{
  FS=OFS="\t"
  np=split(order, P, " "); ns=split(samples, S, " ")
  while((getline line < ARGV[1]) > 0){ split(line,f,"\t"); v[f[1] SUBSEP f[2]]=f[3] }
  printf "ptm"
  for(j=1;j<=ns;j++) if(S[j]!="") printf "\t%s", S[j]
  printf "\n"
  for(i=1;i<=np;i++){
    if(P[i]=="") continue
    printf "%s", P[i]
    for(j=1;j<=ns;j++){ if(S[j]=="") continue; k=S[j] SUBSEP P[i]; printf "\t%d", (k in v ? v[k] : 0) }
    printf "\n"
  }
  exit
}' "$TMPROOT/all_ptm_pairs.tsv" > "$OUTDIR/kmet_barcode_matrix.tsv"

# specificity matrix: each PTM as % of that sample's on-target PTM
awk 'BEGIN{FS=OFS="\t"}
NR==FNR{ if(FNR>1) tgt[$1]=$2; next }
FNR==1{ n=NF; for(i=2;i<=NF;i++) col[i]=$i; print; next }
{ ptm[FNR]=$1; for(i=2;i<=NF;i++) M[FNR,i]=$i; rows=FNR; next }
END{
  for(r=2;r<=rows;r++){}
}' "$OUTDIR/kmet_summary.tsv" "$OUTDIR/kmet_barcode_matrix.tsv" > /dev/null 2>&1 || true

python3 - "$OUTDIR" <<'PY' 2>/dev/null || echo "NOTE: python3 not found; skipping kmet_specificity.tsv" >&2
import sys, os, csv
out = sys.argv[1]
mat = os.path.join(out, "kmet_barcode_matrix.tsv")
summ = os.path.join(out, "kmet_summary.tsv")
tgt = {}
with open(summ) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        tgt[row["sample"]] = row["target"]
with open(mat) as fh:
    rows = [l.rstrip("\n").split("\t") for l in fh]
hdr, body = rows[0], rows[1:]
samples = hdr[1:]
with open(os.path.join(out, "kmet_specificity.tsv"), "w") as fh:
    fh.write("\t".join(hdr) + "\n")
    denom = []
    for s in samples:
        t = tgt.get(s, "NA")
        d = None
        for r in body:
            if r[0] == t:
                d = float(r[1 + samples.index(s)])
        denom.append(d)
    for r in body:
        vals = []
        for i, s in enumerate(samples):
            d = denom[i]
            v = float(r[1 + i])
            vals.append("NA" if not d else f"{100*v/d:.2f}")
        fh.write("\t".join([r[0]] + vals) + "\n")
PY

# -----------------------------------------------------------------------------
# Provenance
# -----------------------------------------------------------------------------
{
  echo "kmet_spikein_count.sh v${VERSION}"
  echo "date:      $(date -Iseconds)"
  echo "host:      $(hostname)"
  echo "user:      ${USER:-unknown}"
  echo "pwd:       $(pwd)"
  echo "command:   $0 $*"
  echo "bash:      ${BASH_VERSION}"
  echo "awk:       $(awk --version 2>&1 | head -n1 || awk -W version 2>&1 | head -n1)"
  echo "grep:      $(grep --version | head -n1)"
  echo "samtools:  $(command -v samtools >/dev/null 2>&1 && samtools --version | head -n1 || echo 'not found')"
  echo "mapq:      ${MAPQ}"
  echo "samflags:  ${SAMFLAGS}"
  echo "barcodes:  ${BCFILE}  md5=$(md5sum "$BCFILE" | cut -d' ' -f1)"
  echo "panel:     SNAP-CUTANA K-MetStat (EpiCypher 19-1002), User Guide v1.2 Table 3"
  echo "samples:   ${SAMPLE_LIST[*]}"
} > "$OUTDIR/MANIFEST_kmet.txt"

echo >&2
echo "Done. Outputs in: $OUTDIR" >&2
echo "  kmet_summary.tsv           <- read this first" >&2
echo "  kmet_barcode_matrix.tsv    <- input for kmet_spikein_normfactors.R" >&2
echo "  kmet_specificity.tsv       <- off-target recovery (<20% = pass)" >&2
echo "  per_sample/*.tsv           <- vendor-xlsx-compatible R1/R2 counts" >&2