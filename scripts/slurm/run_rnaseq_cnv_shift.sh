#!/bin/bash
#SBATCH --job-name=rnaseq_cnv
#SBATCH --account=ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# RNA-seq confirmation of Hi-C CNV calls: 2 test contrasts (B10 vs Parental, DMSO & UNC)
# + 2 null contrasts (UNC vs DMSO within Parental and within B10).
#
# Submit from the project root:
#   sbatch scripts/analysis/run_rnaseq_cnv_shift.sh [hic_resolution=250000]
#
# DESeq2 tables are found automatically under DESEQ_DIR by their contrast names.
# If a name differs, set the path explicitly in the OVERRIDES block below.
set -euo pipefail

RES="${1:-250000}"
DESEQ_DIR="/lustre/isaac24/proj/UTK0009/Taiwo/2026JUL-Taiwo-A375UNC-RNA/New_Deseq2_Output"
GTF="data/refs/hg19/gencode.v19.annotation.gtf.gz"
SEGMENTS="results/exploratory/hic_cnv/res_${RES}/cnv_segments.tsv"
OUT="results/exploratory/hic_cnv/rnaseq_check_res_${RES}"

# ---- OVERRIDES (leave empty to auto-detect) ----
TEST_DMSO=""   # Bottom-10-DMSO vs Parental-DMSO
TEST_UNC=""    # Bottom-10-UNC  vs Parental-UNC
NULL_PAR=""    # Parental-UNC   vs Parental-DMSO
NULL_B10=""    # Bottom-10-UNC  vs Bottom-10-DMSO

CONDA_ENV="multiomics-py"
if [[ "${CONDA_DEFAULT_ENV:-}" != "$CONDA_ENV" && "${CONDA_PREFIX:-}" != *"/$CONDA_ENV" ]]; then
    if command -v conda >/dev/null 2>&1; then
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate "$CONDA_ENV"
    else
        echo "ERROR: conda not found. Run 'conda activate $CONDA_ENV' before sbatch." >&2; exit 1
    fi
fi
mkdir -p logs "$OUT"

# find exactly one "*<numerator>_vs_<denominator>_full.csv"; refuse if only the reversed one exists
find_table() {
    local num="$1" den="$2" hits
    mapfile -t hits < <(find "$DESEQ_DIR" -type f -name "*${num}_vs_${den}_full.csv" 2>/dev/null)
    if (( ${#hits[@]} == 1 )); then echo "${hits[0]}"; return 0; fi
    if (( ${#hits[@]} > 1 )); then
        echo "ERROR: several tables match ${num}_vs_${den}:" >&2; printf '  %s\n' "${hits[@]}" >&2
        echo "Set the path in the OVERRIDES block." >&2; return 1
    fi
    if find "$DESEQ_DIR" -type f -name "*${den}_vs_${num}_full.csv" | grep -q .; then
        echo "ERROR: only the reversed contrast ${den}_vs_${num} exists. Test contrasts need B10 in the" >&2
        echo "numerator; for a null contrast the direction does not matter - set it in OVERRIDES." >&2
        return 1
    fi
    echo "ERROR: no table matching *${num}_vs_${den}_full.csv under $DESEQ_DIR" >&2
    echo "Available tables:" >&2; find "$DESEQ_DIR" -type f -name "*_full.csv" | sed 's/^/  /' >&2
    return 1
}

[[ -z "$TEST_DMSO" ]] && TEST_DMSO=$(find_table "Bottom-10-DMSO" "Parental-DMSO")
[[ -z "$TEST_UNC"  ]] && TEST_UNC=$(find_table  "Bottom-10-UNC"  "Parental-UNC")
[[ -z "$NULL_PAR"  ]] && NULL_PAR=$(find_table  "Parental-UNC"   "Parental-DMSO")
[[ -z "$NULL_B10"  ]] && NULL_B10=$(find_table  "Bottom-10-UNC"  "Bottom-10-DMSO")

for f in "$TEST_DMSO" "$TEST_UNC" "$NULL_PAR" "$NULL_B10" "$GTF" "$SEGMENTS"; do
    [[ -s "$f" ]] || { echo "ERROR: missing input: $f" >&2; exit 1; }
done

echo "test  B10vsPar_DMSO : $TEST_DMSO"
echo "test  B10vsPar_UNC  : $TEST_UNC"
echo "null  Par_UNCvsDMSO : $NULL_PAR"
echo "null  B10_UNCvsDMSO : $NULL_B10"
echo "segments            : $SEGMENTS"
echo "output              : $OUT"

python scripts/analysis/rnaseq_cnv_shift.py \
  --test "B10vsPar_DMSO:$TEST_DMSO" \
  --test "B10vsPar_UNC:$TEST_UNC" \
  --null "Par_UNCvsDMSO:$NULL_PAR" \
  --null "B10_UNCvsDMSO:$NULL_B10" \
  --gtf "$GTF" \
  --segments "$SEGMENTS" \
  --outdir "$OUT"

echo "done -> $OUT"
