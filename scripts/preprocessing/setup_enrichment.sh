#!/usr/bin/env bash
# =============================================================================
# 00_setup_enrichment.sh
#   ONE-TIME setup for region_pathway_enrichment.R on ISAAC.
#   RUN THIS ON A LOGIN NODE (login nodes have internet; compute nodes do not).
#
#   It (1) opens a clean shell, (2) installs every R package the script needs
#   into the multiomics-r conda env, (3) downloads the LOLACore region database
#   to project storage, and (4) verifies everything loads. After this, the
#   analysis itself runs fully OFFLINE on compute nodes (except KEGG — see notes).
#
#   Usage:   bash 00_setup_enrichment.sh
# =============================================================================
set -euo pipefail

# ---- project paths (EDIT if yours differ) -----------------------------------
PROJ=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration
ENV=multiomics-r
REFS="$PROJ/refs"                  # big downloads live here (NOT in $HOME — quota)
mkdir -p "$REFS"

# ---- 1) clean shell, then activate conda ------------------------------------
# Your known pitfalls: module-loaded spack packages get injected onto PYTHONPATH
# and a stray ~/.local / user R library can shadow the conda-pinned versions.
# So: purge modules, drop PYTHONPATH, force R to install into the ENV library.
module purge 2>/dev/null || true
unset PYTHONPATH || true
unset PYTHONHOME || true
# initialise + activate conda (adjust the first line if conda lives elsewhere)
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV"

# force R package installs into the env library (prevents ~/R/... shadowing)
export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"
# keep any package caches on lustre, not in $HOME (quota safety)
export XDG_CACHE_HOME="$PROJ/.cache"
export R_USER_CACHE_DIR="$PROJ/.cache/R"
mkdir -p "$R_USER_CACHE_DIR"

echo "[setup] R: $(which R)"
echo "[setup] .libPaths()[1] should be inside $CONDA_PREFIX:"
Rscript -e 'cat("   ", .libPaths()[1], "\n")'

# ---- 2) install R packages (skips ones already present) ---------------------
# Bioconductor + CRAN. chipenrich.data / TxDb / reactome.db are the big ones.
cat > "$REFS/_install_pkgs.R" <<'REOF'
options(Ncpus = max(1L, parallel::detectCores() - 1L),
        repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc <- c(
  "chipenrich", "chipenrich.data",                 # Broad-Enrich + its loci/genesets
  "clusterProfiler", "ReactomePA", "reactome.db",  # ORA + Reactome (local)
  "org.Hs.eg.db", "enrichplot", "DOSE",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",             # bin -> gene (hg19)
  "GenomicFeatures", "GenomicRanges", "GenomeInfoDb", "S4Vectors",
  "LOLA")
cran <- c("optparse", "ggplot2", "msigdbr", "data.table", "simpleCache")

need_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
need_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]

if (length(need_cran)) { cat("[install] CRAN:", paste(need_cran, collapse=", "), "\n")
  install.packages(need_cran) }
if (length(need_bioc)) { cat("[install] Bioc:", paste(need_bioc, collapse=", "), "\n")
  BiocManager::install(need_bioc, update = FALSE, ask = FALSE) }

cat("[install] done.\n")
# msigdbr note: if it complains it needs the data package 'msigdbdf', install from
# the repo it prints, e.g.:
#   install.packages("msigdbdf", repos = "https://igordot.r-universe.dev")
REOF
Rscript "$REFS/_install_pkgs.R"

# ---- 3) download LOLACore (only needed for the --lola_db companion) ----------
# Cached DB (~200 MB, loads in ~30 s, needs simpleCache — installed above).
# hg19 + hg38 + mm10 are all inside. Skip this whole block if you won't use LOLA.
LOLA_TGZ="$REFS/LOLACoreCaches_170206.tgz"
if [ ! -e "$REFS/LOLACore_marker" ]; then
  echo "[setup] downloading LOLACore cached DB ..."
  ( cd "$REFS" && wget -c "http://cloud.databio.org.s3.amazonaws.com/regiondb/LOLACoreCaches_170206.tgz" )
  ( cd "$REFS" && tar xzf "$LOLA_TGZ" )
  touch "$REFS/LOLACore_marker"
fi
# find the unpacked hg19 folder and expose it at a stable path
LOLA_HG19="$(find "$REFS" -type d -name hg19 2>/dev/null | grep -i 'LOLACore' | head -1 || true)"
if [ -n "${LOLA_HG19:-}" ]; then
  mkdir -p "$REFS/LOLACore"
  ln -sfn "$LOLA_HG19" "$REFS/LOLACore/hg19"
  echo "[setup] LOLACore hg19 ready -> $REFS/LOLACore/hg19"
else
  echo "[setup][WARN] could not locate the unpacked LOLACore/hg19 folder; check $REFS"
fi

# ---- 4) verify everything loads + the data is reachable ---------------------
cat > "$REFS/_verify.R" <<'REOF'
ok <- TRUE
load1 <- function(p) tryCatch({ suppressMessages(library(p, character.only=TRUE)); TRUE },
                              error=function(e){ cat("  [FAIL]", p, ":", conditionMessage(e), "\n"); FALSE })
pkgs <- c("optparse","ggplot2","msigdbr","data.table","simpleCache",
          "chipenrich","clusterProfiler","ReactomePA","org.Hs.eg.db","enrichplot",
          "GenomicFeatures","GenomicRanges","GenomeInfoDb","S4Vectors","LOLA",
          "TxDb.Hsapiens.UCSC.hg19.knownGene")
for (p in pkgs) ok <- load1(p) & ok

# chipenrich data present? (region mode)
cat("[verify] chipenrich genomes:",
    paste(head(chipenrich::supported_genomes()), collapse=", "), "...\n")
cat("[verify] chipenrich genesets incl GOBP/reactome?:",
    all(c("GOBP","reactome") %in% chipenrich::supported_genesets()$geneset), "\n")

# msigdbr Hallmark + Ben-Porath PRC2 reachable? (positive control)
mh <- tryCatch(nrow(msigdbr::msigdbr(species="Homo sapiens", collection="H")),
               error=function(e) tryCatch(nrow(msigdbr::msigdbr(species="Homo sapiens", category="H")),
                                          error=function(e2){cat("  [FAIL] msigdbr:",conditionMessage(e2),"\n"); 0}))
bp <- tryCatch({ x <- msigdbr::msigdbr(species="Homo sapiens", collection="C2", subcollection="CGP")
                 sum(grepl("BENPORATH|SUZ12|EED|PRC2|H3K27ME3", x$gs_name, ignore.case=TRUE)) },
               error=function(e) tryCatch({ x <- msigdbr::msigdbr(species="Homo sapiens", category="C2", subcategory="CGP")
                 sum(grepl("BENPORATH|SUZ12|EED|PRC2|H3K27ME3", x$gs_name, ignore.case=TRUE)) },
               error=function(e2){cat("  [FAIL] msigdbr C2:CGP:",conditionMessage(e2),"\n"); 0}))
cat(sprintf("[verify] msigdbr Hallmark rows=%d ; Polycomb/PRC2 gene-set rows=%d\n", mh, bp))
if (bp == 0) { cat("  [FAIL] no Ben-Porath/PRC2 sets — positive control will be empty.\n"); ok <- FALSE }

# TxDb gene fetch works? (bin mode)
ng <- tryCatch(length(GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)),
               error=function(e){cat("  [FAIL] TxDb genes:",conditionMessage(e),"\n"); 0})
cat("[verify] hg19 genes from TxDb:", ng, "\n"); if (ng < 1000) ok <- FALSE

# LOLACore loads? (only if the folder exists)
ld <- Sys.getenv("LOLA_HG19", "")
if (nzchar(ld) && dir.exists(ld)) {
  rc <- tryCatch({ db <- LOLA::loadRegionDB(ld); length(db$regionGRL) },
                 error=function(e){cat("  [FAIL] LOLA load:",conditionMessage(e),"\n"); 0})
  cat("[verify] LOLACore hg19 tracks loaded:", rc, "\n"); if (rc < 1) ok <- FALSE
} else cat("[verify] (LOLACore folder not set; skip LOLA check)\n")

cat(if (ok) "\n[VERIFY] ALL GREEN — ready to run.\n" else "\n[VERIFY] problems above — fix before running.\n")
REOF
LOLA_HG19="$REFS/LOLACore/hg19" Rscript "$REFS/_verify.R"

echo
echo "[setup] complete."
echo "       LOLA db path for --lola_db :  $REFS/LOLACore/hg19"
echo "       Now submit analysis jobs with run_enrichment.slurm (compute node, offline)."