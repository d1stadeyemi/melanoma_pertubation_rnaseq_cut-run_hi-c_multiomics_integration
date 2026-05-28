#!/usr/bin/env bash
#SBATCH --job-name=average_bigwigs
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/average_bigwigs_%j.out
#SBATCH --error=logs/average_bigwigs_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration
# =============================================================================
# average_replicate_bigwigs.sh
#   Average REPLICATE bigWigs into one track per mark/condition/treatment group,
#   for building ML feature signal. (A375 parent-vs-B10 constricted-migration.)
# =============================================================================
#
# WHAT IT DOES ----------------------------------------------------------------
#   Reads the sample sheet, finds each sample's CPM-normalized bigWig inside its
#   peakcalling directory, groups samples by (Factor, Condition, Treatment), and
#   runs deepTools `bigwigAverage` over the REPLICATES of each group to produce
#   one averaged bigWig per group. Groups with a single replicate are copied
#   (nothing to average). Works whether a group has 2 or 3 replicates.
#
# WHY NAME OUTPUTS FROM THE SHEET, NOT THE bigWig FILE ------------------------
#   Several inner bigWig filenames disagree with their folder/identity (e.g.
#   B10-UNC-H3K27Me3-R2's track is named ...Unmigrated...; Bottom10-H3K27Ac-R2's
#   is named BottomR5-2...). The sample sheet (built from the proofread manifest)
#   holds the TRUE identity, so output names are built from the sheet's
#   Factor/Condition/Treatment — NOT from the bigWig filename. Treatment is
#   included ONLY when it is a real perturbation (not 'none'). Output pattern:
#       <Condition>_<Factor>.mean.bw            e.g. B10_H3K27me3.mean.bw  (none)
#       <Condition>_<Treatment>_<Factor>.mean.bw e.g. B10_UNC_H3K27me3.mean.bw
#
# WHY THE bigWig FOLDER (broad vs narrow) DOES NOT MATTER ----------------------
#   The .cpm.norm.bw is genome-wide COVERAGE; it is the same track regardless of
#   whether it sits under macs2.broad/ or macs2.narrow/ (those folders differ
#   only in PEAK calls, not coverage). So we just search both subfolders and use
#   whichever holds the sample's .cpm.norm.bw. (For H3K27Ac too — same track.)
#
# INPUT bigWig DISCOVERY ------------------------------------------------------
#   The sheet's bamReads path tells us the sample FOLDER. From there we look for:
#       <folder>/peakcalling/macs2.broad/*.cpm.norm.bw
#       <folder>/peakcalling/macs2.narrow/*.cpm.norm.bw
#       <folder>/*.cpm.norm.bw            (fallback: top-level)
#   and take the first match. If a sample has no .cpm.norm.bw, it is skipped
#   with a warning (so the group average uses whatever replicates do have one).
#
# FILTERING (all optional, COMBINE as AND; all case-insensitive) --------------
#   --condition parental   only samples whose Condition matches (e.g. parental/B10)
#   --mark H3K27me3         only samples whose Factor matches
#   --treatment none        only samples whose Treatment matches; use 'none' for
#                           untreated (empty Treatment also counts as none)
#   With several flags, a sample must satisfy ALL of them to be included, e.g.
#     --condition parental --treatment none   -> only untreated parental groups.
#   With no filter flags, ALL samples are processed (original behaviour).
#
# REQUIRES: deepTools (bigwigAverage) on PATH — it's in your multiomics-py env.
#
# USAGE -----------------------------------------------------------------------
#   bash average_replicate_bigwigs.sh \
#       --sheet data/processed/sample_sheet.csv \
#       --outdir data/processed/mean_bigwigs
#   # filtered examples:
#   #   --condition parental                  only parental groups
#   #   --mark H3K27me3 --treatment unc        only UNC-treated H3K27me3 groups
#   #   --condition B10 --treatment none       only untreated B10 groups
#   # other options:
#   #   --bin 50        bigwigAverage bin size in bp (default 50)
#   #   --dry-run       print what would happen, do nothing
#
# OUTPUT: one <Condition>[_<Treatment>]_<Factor>.mean.bw per group in --outdir,
#   plus a groups.tsv recording which input bigWigs went into each average
#   (provenance for the methods section).
# =============================================================================

module load anaconda3/2024.06
source $ANACONDA3_SH
conda activate multiomics-py

set -euo pipefail

# ---- defaults ----
SHEET=""; OUTDIR="mean_bigwigs"; BIN=50; DRYRUN=0
ONLY_MARK=""; ONLY_COND=""; ONLY_TREAT=""

# ---- parse args ----
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sheet)     SHEET="$2"; shift 2;;
    --outdir)    OUTDIR="$2"; shift 2;;
    --bin)       BIN="$2"; shift 2;;
    --mark)      ONLY_MARK="$2"; shift 2;;
    --condition) ONLY_COND="$2"; shift 2;;
    --treatment) ONLY_TREAT="$2"; shift 2;;
    --dry-run)   DRYRUN=1; shift;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -z "$SHEET" ]] && { echo "ERROR: --sheet is required." >&2; exit 1; }
command -v bigwigAverage >/dev/null 2>&1 || {
  echo "ERROR: bigwigAverage not found. Activate the deepTools env (multiomics-py)." >&2; exit 1; }

# Lower-case the filter values once, for case-insensitive comparison.
ONLY_MARK_LC="${ONLY_MARK,,}"; ONLY_COND_LC="${ONLY_COND,,}"; ONLY_TREAT_LC="${ONLY_TREAT,,}"

# Announce active filters so the run log is self-documenting.
[[ -n "$ONLY_COND"  ]] && echo "[INFO] Filter: condition == $ONLY_COND"
[[ -n "$ONLY_MARK"  ]] && echo "[INFO] Filter: mark == $ONLY_MARK"
[[ -n "$ONLY_TREAT" ]] && echo "[INFO] Filter: treatment == $ONLY_TREAT"

mkdir -p "$OUTDIR"
GROUPS_TSV="$OUTDIR/groups.tsv"
echo -e "group_output\tinput_bigwigs" > "$GROUPS_TSV"

# ---- locate the sheet's columns by header name (robust to column order) ------
header=$(head -1 "$SHEET")
col() { echo "$header" | tr ',' '\n' | sed 's/"//g' | grep -nx "$1" | cut -d: -f1; }
C_FACTOR=$(col Factor); C_COND=$(col Condition); C_TREAT=$(col Treatment)
C_BAM=$(col bamReads);  C_ID=$(col SampleID)
[[ -z "$C_FACTOR$C_COND$C_TREAT$C_BAM" ]] && {
  echo "ERROR: sheet missing required columns." >&2; exit 1; }

# ---- find a sample's .cpm.norm.bw from its sample folder ---------------------
find_bw() {
  local bam="$1"
  local folder; folder=$(dirname "$(dirname "$(dirname "$bam")")")
  local hit=""
  for sub in peakcalling/macs2.broad peakcalling/macs2.narrow .; do
    hit=$(ls "$folder/$sub"/*.cpm.norm.bw 2>/dev/null | head -1 || true)
    [[ -n "$hit" ]] && { echo "$hit"; return; }
  done
  echo ""
}

# ---- normalise a Treatment value: empty -> 'none' (so --treatment none works)
norm_treat() { local t="${1,,}"; [[ -z "$t" ]] && echo "none" || echo "$t"; }

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# ---- read the sheet, apply filters (AND), collect bigWigs per group ----------
tail -n +2 "$SHEET" | sed 's/"//g' | while IFS=',' read -r -a F; do
  factor="${F[$((C_FACTOR-1))]}"
  cond="${F[$((C_COND-1))]}"
  treat="${F[$((C_TREAT-1))]}"
  bam="${F[$((C_BAM-1))]}"
  id="${F[$((C_ID-1))]}"

  # --- apply the optional filters; a sample must pass ALL active ones --------
  treat_norm="$(norm_treat "$treat")"
  [[ -n "$ONLY_MARK_LC"  && "${factor,,}"  != "$ONLY_MARK_LC"  ]] && continue
  [[ -n "$ONLY_COND_LC"  && "${cond,,}"    != "$ONLY_COND_LC"  ]] && continue
  [[ -n "$ONLY_TREAT_LC" && "$treat_norm"  != "$ONLY_TREAT_LC" ]] && continue

  bw=$(find_bw "$bam")
  if [[ -z "$bw" ]]; then
    echo "[WARN] no .cpm.norm.bw for $id (folder of $bam) — skipped." >&2
    continue
  fi

  # group key: include Treatment ONLY when it is a real perturbation (not none).
  if [[ "$treat_norm" == "none" ]]; then
    key="${cond}_${factor}"
  else
    key="${cond}_${treat}_${factor}"
  fi
  echo "$bw" >> "$TMP/$key"
  echo "[INFO] $id -> group $key : $(basename "$bw")"
done

# ---- average each group's replicates ----------------------------------------
shopt -s nullglob
group_files=("$TMP"/*)
if [[ "${#group_files[@]}" -eq 0 ]]; then
  echo "[WARN] No samples passed the filters — nothing to average." >&2
fi
for gf in "${group_files[@]}"; do
  key=$(basename "$gf")
  mapfile -t bws < "$gf"
  out="$OUTDIR/${key}.mean.bw"
  echo -e "${out}\t$(printf '%s;' "${bws[@]}")" >> "$GROUPS_TSV"

  if [[ "${#bws[@]}" -eq 0 ]]; then
    continue
  elif [[ "${#bws[@]}" -eq 1 ]]; then
    echo "[INFO] $key has 1 replicate — copying (nothing to average)."
    [[ "$DRYRUN" -eq 1 ]] || cp "${bws[0]}" "$out"
  else
    echo "[INFO] Averaging ${#bws[@]} replicates -> $out"
    if [[ "$DRYRUN" -eq 0 ]]; then
      bigwigAverage -b "${bws[@]}" -o "$out" --binSize "$BIN" -p max
    fi
  fi
done

echo "[INFO] Done. Averaged tracks in $OUTDIR ; provenance in $GROUPS_TSV"
[[ "$DRYRUN" -eq 1 ]] && echo "[INFO] (dry-run: no files written)"

conda deactivate
exit 0