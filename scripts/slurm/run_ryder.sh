#!/bin/bash
#SBATCH --job-name=ryder_run
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G                       
#SBATCH --time=2:00:00
#SBATCH --output=logs/ryder_run_%j.out
#SBATCH --error=logs/ryder_run_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

set -euo pipefail

log()     { echo -e "[$(date '+%F %T')] $*"; }
die()     { echo -e "[$(date '+%F %T')] ERROR: $*" >&2; exit 1; }
require() { [[ -e "$1" ]] || die "missing required file: $1"; log "  ok: $1"; }
run()     { log "RUN: $*"; "$@"; }

log "=============================================================="
log "Ryder run starting on $(hostname)"
module load anaconda3/2024.06
source "$ANACONDA3_SH"
conda activate ryder
log "env: $(which paw.py) | python: $(python --version 2>&1)"

NP="${SLURM_CPUS_PER_TASK:-10}"
EXT=50000

SHEET=data/processed/sample_sheet.csv
AVG=data/processed/ryder/avg_tracks
CTCF=data/processed/ryder/CTCF_constitutive_hg19.bed
CSF=data/processed/greenlist/hg19.chrom.sizes
PEAKS=data/processed/ryder/H3K27me3_peaks_hg19.bed

PAR_NONE=$AVG/Parental-H3K27Me3_avg_rep.rpm.bw
PAR_UNC=$AVG/Parental-UNC-H3K27Me3_avg_rep.rpm.bw
B10_NONE=$AVG/B10-H3K27Me3_avg_rep.rpm.bw
B10_UNC=$AVG/B10-UNC-H3K27Me3_avg_rep.rpm.bw

OUT=data/processed/ryder/run
mkdir -p "$OUT/parental" "$OUT/b10"

# ---- STEP 1: ROI (skip if already built) ----
if [[ -s "$PEAKS" ]]; then
  log "STEP 1: ROI exists ($(wc -l < "$PEAKS") regions) — skipping rebuild"
else
  log "STEP 1: building H3K27me3 ROI (Factor=H3K27me3, Replicate in {1,2})"
  require "$SHEET"
  SHEET="$SHEET" PEAKS_OUT="$PEAKS" python - <<'PY'
import os, sys
import pandas as pd
sheet=os.environ["SHEET"]; out=os.environ["PEAKS_OUT"]
CANON={f"chr{c}" for c in list(range(1,23))+["X","Y"]}
df=pd.read_csv(sheet); df.columns=[c.strip() for c in df.columns]
df["__rep"]=pd.to_numeric(df["Replicate"],errors="coerce")
sub=df[(df["Factor"].astype(str).str.strip()=="H3K27me3") & (df["__rep"].isin([1,2]))]
peaks=[str(p).strip() for p in sub["Peaks"]]
miss=[p for p in peaks if not os.path.exists(p)]
if miss: sys.exit("ERROR: missing peak files:\n  "+"\n  ".join(miss))
allp=pd.concat([pd.read_csv(p,sep="\t",header=None,usecols=[0,1,2],names=["chrom","start","end"]) for p in peaks],ignore_index=True)
allp=allp[allp["chrom"].isin(CANON)].sort_values(["chrom","start","end"])
merged=[]
for chrom,g in allp.groupby("chrom",sort=True):
    cs=ce=None
    for s,e in zip(g["start"].astype(int),g["end"].astype(int)):
        if cs is None: cs,ce=s,e
        elif s<=ce: ce=max(ce,e)
        else: merged.append((chrom,cs,ce)); cs,ce=s,e
    if cs is not None: merged.append((chrom,cs,ce))
pd.DataFrame(merged,columns=["chrom","start","end"]).to_csv(out,sep="\t",header=False,index=False)
print(f"  WROTE {len(merged)} ROI regions -> {out}",file=sys.stderr)
PY
fi
require "$PEAKS"

# ---- STEP 2: pre-flight ----
log "STEP 2: verifying inputs"
require "$CTCF"; require "$CSF"
require "$PAR_NONE"; require "$PAR_UNC"; require "$B10_NONE"; require "$B10_UNC"

# ---- STEP 3: PARENTAL ----
log "STEP 3: PARENTAL contrast"
PAR_NORM="$OUT/parental/paw_parental_Parental_UNC.bw"
if [[ -s "$PAR_NORM" ]]; then
  log "  SKIP paw parental — $PAR_NORM already exists"
else
  run paw.py -r "$CTCF" -c "$PAR_NONE" -t "$PAR_UNC" \
      -lc Parental_none -lt Parental_UNC -csf "$CSF" -ext "$EXT" -p "$NP" \
      -o "$OUT/parental/paw_parental"
fi
require "$PAR_NORM"
run patrol.py -r "$PEAKS" -c "$PAR_NONE" -t "$PAR_NORM" \
    -lc Parental_none -lt Parental_UNC -mode FC -pcut 0.01 \
    -o "$OUT/parental/patrol_parental"
log "STEP 3 done"

# ---- STEP 4: B10 ----
log "STEP 4: B10 contrast"
B10_NORM="$OUT/b10/paw_b10_B10_UNC.bw"
if [[ -s "$B10_NORM" ]]; then
  log "  SKIP paw b10 — $B10_NORM already exists"
else
  run paw.py -r "$CTCF" -c "$B10_NONE" -t "$B10_UNC" \
      -lc B10_none -lt B10_UNC -csf "$CSF" -ext "$EXT" -p "$NP" \
      -o "$OUT/b10/paw_b10"
fi
require "$B10_NORM"
run patrol.py -r "$PEAKS" -c "$B10_NONE" -t "$B10_NORM" \
    -lc B10_none -lt B10_UNC -mode FC -pcut 0.01 \
    -o "$OUT/b10/patrol_b10"
log "STEP 4 done"

log "ALL DONE. Inspect paw QC (*_4_NoiseCutoff.pdf + aggregate-at-reference) before trusting patrol."