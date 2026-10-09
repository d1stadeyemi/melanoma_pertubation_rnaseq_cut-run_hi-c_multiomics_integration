#!/bin/bash
#SBATCH --job-name=hic_cnv
#SBATCH --account=ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Differential CNV (B10 vs Parental) from Hi-C coverage. Submit from the project root
# with `conda activate multiomics-py` already done in your shell (sbatch inherits it).
set -euo pipefail

M=data/processed/hic/mcool
RES="${1:-250000}"
OUT="results/exploratory/hic_cnv/res_${RES}"
SWITCH_BED="${2:-}"     # optional: BED of compartment-switch bins

CONDA_ENV="multiomics-py"
# Use the env already active in the submitting shell; otherwise try to activate it.
if [[ "${CONDA_DEFAULT_ENV:-}" != "$CONDA_ENV" && "${CONDA_PREFIX:-}" != *"/$CONDA_ENV" ]]; then
    if command -v conda >/dev/null 2>&1; then
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate "$CONDA_ENV"
    else
        echo "ERROR: conda not found. Run 'conda activate $CONDA_ENV' before sbatch." >&2; exit 1
    fi
fi
mkdir -p logs "$OUT"
which python; python -c "import cooler, scipy, matplotlib; print('cooler', cooler.__version__)"

EXTRA=()
[[ -n "$SWITCH_BED" ]] && EXTRA+=(--switch-bed "$SWITCH_BED")

python scripts/analysis/hic_cnv_ratio.py \
  --resolution "$RES" \
  --test "B10vsPar_DMSO:$M/B10_DMSO.hg19.mcool:$M/Parental_DMSO.hg19.mcool" \
  --test "B10vsPar_UNC:$M/B10_UNC.hg19.mcool:$M/Parental_UNC.hg19.mcool" \
  --null "Par_UNCvsDMSO:$M/Parental_UNC.hg19.mcool:$M/Parental_DMSO.hg19.mcool" \
  --null "B10_UNCvsDMSO:$M/B10_UNC.hg19.mcool:$M/B10_DMSO.hg19.mcool" \
  --outdir "$OUT" \
  "${EXTRA[@]}"
