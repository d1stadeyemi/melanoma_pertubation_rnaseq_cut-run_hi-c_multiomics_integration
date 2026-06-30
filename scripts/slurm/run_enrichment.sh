#!/bin/bash
# =============================================================================
# run_enrichment.slurm
#   Submit region_pathway_enrichment.R on ISAAC condo-rmccord.
#   Your compute node HAS outbound internet (verified), so KEGG (live enrichKEGG)
#   and any other online lookups work fine here. Everything else is local after
#   you have run 00_setup_enrichment.sh once.
#
#   Submit:  sbatch run_enrichment.slurm
# =============================================================================
#SBATCH --job-name=enrich_parunc
#SBATCH --account=ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8          
#SBATCH --mem=64G                  
#SBATCH --time=04:00:00
#SBATCH --output=logs/enrich_parunc_%j.out
#SBATCH --error=logs/enrich_parunc_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

set -euo pipefail
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}   # quiet the sbatch OpenMP notice

# ---- clean shell + activate the env -----------------
module purge 2>/dev/null || true
unset PYTHONPATH || true
unset PYTHONHOME || true
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate multiomics-r
export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"

PROJ=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration
LOLA="$PROJ/data/refs/LOLACore/hg19"      
# keep R/Bioconductor caches on lustre (NOT $HOME) — same as setup
export XDG_CACHE_HOME="$PROJ/data/.cache"
export R_USER_CACHE_DIR="$PROJ/data/.cache/R"
mkdir -p "$R_USER_CACHE_DIR"
SCRIPT="$PROJ/scripts/analysis/pathway_enrichment_regions.R"
cd "$PROJ"
mkdir -p logs

# =============================================================================
# EXAMPLE A — CUT&RUN broad H3K27me3 (csaw), lost vs gained, with LOLA
#   significant.csv = foreground ; all_regions.csv = universe
# =============================================================================
CSAW="results/exploratory/diff_analysis/csaw/Migrated_Parunc/spikein"
Rscript "$SCRIPT" \
  --mode region \
  --regions   "$CSAW/migrated_ParentalUNC_vs_ParentalDMSO_spikein_significant.csv" \
  --background "$CSAW/migrated_ParentalUNC_vs_ParentalDMSO_spikein_all_regions.csv" \
  --select "lost=call:lost,gained=call:resistant" \
  --genome hg19 --locusdef nearest_tss \
  --databases GO_BP,Reactome,Hallmark,KEGG,PRC2 \
  --qc_plots --randomization bylength --cores 8 \
  --lola_db "$LOLA" \
  --outdir "$CSAW/enrichment_h3k27me3_Migrated_Parunc"

# =============================================================================
# EXAMPLE B — Hi-C compartments straight from Reversal_classified.txt
#   (uncomment to run; universe = all rows of the table)
# =============================================================================
# REV=results/exploratory/reversal
# Rscript "$SCRIPT" \
#   --mode bin \
#   --regions "$REV/Reversal_classified.txt" \
#   --select "specific=migration_specific_reversal,reverted=reversal_status:reverted,switch_AtoB=mig_switch:A_to_B,switch_BtoA=mig_switch:B_to_A" \
#   --genome hg19 --databases GO_BP,Reactome,Hallmark,KEGG,PRC2 \
#   --collapse_adjacent --focal specific --cores 8 \
#   --lola_db "$LOLA" \
#   --outdir results/exploratory/enrichment_compartment

echo "[slurm] done."