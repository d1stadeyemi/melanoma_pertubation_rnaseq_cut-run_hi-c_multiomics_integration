#!/bin/bash
#SBATCH --job-name=vp2mcool
#SBATCH --account=ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo-
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# HiC-Pro *.allValidPairs (already deduplicated) -> multi-resolution,
# ICE-balanced .mcool with cooler.
#
# Usage (submit from the project root):
#   sbatch --job-name=mcool_B10_DMSO scripts/preprocessing/validpairs_to_mcool.sh \
#       <sample_name> <path/to/sample.allValidPairs> [min_mapq=0] [base_res=10000]
# ---------------------------------------------------------------------------
set -euo pipefail

SAMPLE="${1:?sample name required}"
VALIDPAIRS="${2:?path to allValidPairs required}"
MIN_MAPQ="${3:-0}"       # 0 = keep HiC-Pro's MAPQ>=10; 30 = stricter
BASE_RES="${4:-10000}"   # same base for every sample keeps bins identical
CONDA_ENV="multiomics-py"

# Same chrom sizes HiC-Pro used (hg19) -> identical chromosome order/naming
CHROMSIZES="/lustre/isaac24/proj/UTK0009/HiCPro/files/hg19/chrom_hg19.sizes"

# Every level must be a multiple of BASE_RES; includes the lab's standard bins
RESOLUTIONS="${BASE_RES}"
for r in 10000 20000 40000 50000 100000 250000 500000 1000000 2500000; do
    if (( r > BASE_RES && r % BASE_RES == 0 )); then RESOLUTIONS+=",${r}"; fi
done

OUTDIR="data/processed/hic/mcool"
WORKDIR="${OUTDIR}/tmp_${SAMPLE}_${SLURM_JOB_ID:-local}"   # on Lustre, not node /tmp
mkdir -p "$OUTDIR" "$WORKDIR" logs

# Use the env already active in the submitting shell; otherwise try to activate it.
if [[ "${CONDA_DEFAULT_ENV:-}" != "$CONDA_ENV" && "${CONDA_PREFIX:-}" != *"/$CONDA_ENV" ]]; then
    if command -v conda >/dev/null 2>&1; then
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate "$CONDA_ENV"
    else
        echo "ERROR: conda not found. Run 'conda activate $CONDA_ENV' before sbatch." >&2; exit 1
    fi
fi

echo "[$(date)] sample=$SAMPLE base=$BASE_RES res=$RESOLUTIONS mapq>=$MIN_MAPQ"
echo "pairs: $VALIDPAIRS"
cooler --version

# --- 1. Sanity checks on the first 5M lines (full pass on 48 GB is wasteful) ---
head -n 5000000 "$VALIDPAIRS" > "$WORKDIR/head.pairs"
NCOL=$(awk -F'\t' '{print NF}' "$WORKDIR/head.pairs" | sort -u | tr '\n' ' ')
echo "column count(s): $NCOL"   # HiC-Pro v3 writes 12 fields + 1 trailing -> 13
MINCOL=$(awk -F'\t' 'NR==1||NF<m{m=NF} END{print m}' "$WORKDIR/head.pairs")
if (( MIN_MAPQ > 0 && MINCOL < 12 )); then
    echo "ERROR: MAPQ filter requested but <12 columns present" >&2; exit 1
fi
cut -f2 "$WORKDIR/head.pairs" | sort -u > "$WORKDIR/pair_chroms.txt"
cut -f1 "$CHROMSIZES" | sort -u > "$WORKDIR/size_chroms.txt"
MISSING=$(comm -23 "$WORKDIR/pair_chroms.txt" "$WORKDIR/size_chroms.txt" || true)
if [[ -n "$MISSING" ]]; then echo "WARNING: chroms not in chrom sizes (dropped):"; echo "$MISSING"; fi
rm -f "$WORKDIR/head.pairs"

# --- 2. Filter (chrom whitelist + optional MAPQ) and bin at base resolution ---
BASE_COOL="$WORKDIR/${SAMPLE}.${BASE_RES}.cool"
echo "[$(date)] cooler cload pairs -> $BASE_COOL"
awk -v q="$MIN_MAPQ" 'BEGIN{FS=OFS="\t"}
     NR==FNR {keep[$1]=1; next}
     ($2 in keep) && ($5 in keep) && (q==0 || ($11+0>=q && $12+0>=q))' \
     "$CHROMSIZES" "$VALIDPAIRS" \
| cooler cload pairs \
    -c1 2 -p1 3 -c2 5 -p2 6 \
    --assembly hg19 \
    "${CHROMSIZES}:${BASE_RES}" \
    - \
    "$BASE_COOL"
cooler info "$BASE_COOL"

# --- 3. Multi-resolution + ICE balancing ("weight") at every level ---
MCOOL="$OUTDIR/${SAMPLE}.hg19.mcool"
echo "[$(date)] cooler zoomify -> $MCOOL"
cooler zoomify \
    --nproc "${SLURM_CPUS_PER_TASK:-8}" \
    --resolutions "$RESOLUTIONS" \
    --balance \
    -o "$MCOOL" \
    "$BASE_COOL"

# --- 4. Verify ---
cooler ls "$MCOOL"
cooler info "${MCOOL}::resolutions/250000"

rm -rf "$WORKDIR"
echo "[$(date)] done: $MCOOL"
