#!/usr/bin/env Rscript
# =============================================================================
# liftover_ctcf_hg19.R
# -----------------------------------------------------------------------------
# Lift Ryder's hg38 constitutive-CTCF reference BED to hg19 for use as paw.py -r.
# Reuses the liftOver QC core of prepare_greenlist_hg19.R, but reads a BED (not
# Excel) and DROPS the greenlist-specific filters (gene-distance, IgG), which
# would wrongly remove CTCF sites (they sit AT genes / boundaries, not distant).
#
# PIPELINE
#   1. Read hg38 BED -> GRanges (rtracklayer handles 0-based + any header).
#   2. Canonical chromosomes only (chr1-22,X,Y).
#   3. liftOver hg38->hg19; per region drop if unmapped / multi-chromosome /
#      lifted span width outside [size_lo,size_hi]x original / lifted pieces
#      cover < min_cov_frac of the span (gap-spanning).
#   4. Canonical hg19 chroms; (optionally) merge overlaps.
#   5. (optional) drop regions overlapping an hg19 ENCODE blacklist.
#   6. (optional) pad regions narrower than --min_width (Ryder ignores <100 bp).
#   7. Write hg19 BED (0-based, 3 col) + funnel/provenance report.
#
# RUN in your R env (rtracklayer + GenomicRanges; the env that ran the greenlist):
#   conda activate multiomics-r
#   Rscript scripts/preprocessing/liftover_ctcf_hg19.R \
#       --bed ryder/refData/hg38_CTCF_constitutive.bed \
#       --outdir data/processed/ryder \
#       --out data/processed/ryder/CTCF_constitutive_hg19.bed \
#       --chain data/processed/greenlist/hg38ToHg19.over.chain.gz \
#       --blacklist data/processed/greenlist/hg19-blacklist.v2.bed \
#       --min_width 200
# =============================================================================

suppressWarnings(suppressMessages({
  library(optparse)
  library(rtracklayer)
  library(GenomicRanges)
}))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--bed", type="character", help="REQUIRED. hg38 CTCF BED (>=3 cols, 0-based)."),
  make_option("--outdir", type="character", default=".", help="Output dir. [%default]"),
  make_option("--out", type="character", default=NULL, help="Output BED. [<outdir>/CTCF_constitutive_hg19.bed]"),
  make_option("--chain", type="character", default=NULL, help="hg38ToHg19.over.chain(.gz). Auto-downloaded if absent. [auto]"),
  make_option("--size_lo", type="double", default=0.8, help="Min lifted-span/original width ratio. [%default]"),
  make_option("--size_hi", type="double", default=1.25, help="Max lifted-span/original width ratio. [%default]"),
  make_option("--min_cov_frac", type="double", default=0.8, help="Min summed-piece/span coverage. [%default]"),
  make_option("--blacklist", type="character", default=NULL, help="Optional hg19 ENCODE blacklist BED to subtract. [none]"),
  make_option("--min_width", type="integer", default=0L, help="Pad regions narrower than this (bp), centered. Ryder ignores <100 bp; 200 is safe. 0=off. [%default]"),
  make_option("--no_merge", action="store_true", default=FALSE, help="Do NOT merge overlapping lifted regions (default merges).")
)))
if (is.null(opt$bed)) stop("--bed is required.")
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
out_bed <- if (is.null(opt$out)) file.path(opt$outdir,"CTCF_constitutive_hg19.bed") else opt$out
out_report <- sub("\\.bed$",".report.txt",out_bed); if (out_report==out_bed) out_report <- paste0(out_bed,".report.txt")
CANON <- paste0("chr", c(1:22,"X","Y"))

LOG <- character(0); say <- function(...){l<-paste0(...);cat(l,"\n");LOG<<-c(LOG,l)}
funnel <- character(0); step <- function(label,n){line<-sprintf("  %-46s %7d",label,n);cat(line,"\n");funnel<<-c(funnel,line)}
say("=== liftover_ctcf_hg19.R ===")

# 1. read hg38 BED ----------------------------------------------------------
gl <- import(opt$bed)                                  # rtracklayer: BED 0-based -> 1-based GRanges, handles headers
gl <- granges(gl)                                      # drop name/score/strand; keep ranges only
suppressWarnings(try(seqlevelsStyle(gl) <- "UCSC", silent=TRUE))
say(sprintf("[READ] %s: %d regions; median width %d bp.", basename(opt$bed), length(gl), as.integer(median(width(gl)))))
say("[FUNNEL] region survival:"); step("input (hg38)", length(gl))

# 2. canonical chroms (hg38) ------------------------------------------------
gl <- gl[as.character(seqnames(gl)) %in% CANON]
gl <- keepSeqlevels(gl, intersect(CANON, seqlevelsInUse(gl)), pruning.mode="coarse")
step("after canonical-chrom (hg38)", length(gl))
orig_w <- width(gl)

# 3. chain + liftOver -------------------------------------------------------
resolve_chain <- function(chain_arg, outdir){
  url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
  path <- chain_arg
  if (is.null(path)) {
    path <- file.path(outdir,"hg38ToHg19.over.chain.gz")
    if (!file.exists(path)) {
      say(sprintf("[CHAIN] downloading %s", url))
      ok <- tryCatch({download.file(url,path,mode="wb",quiet=TRUE);TRUE}, error=function(e) FALSE)
      if (!ok || !file.exists(path)) stop("Chain download failed; fetch on a login node and pass --chain.")
    }
  }
  if (!file.exists(path)) stop("--chain not found: ", path)
  if (grepl("\\.gz$", path)) { plain<-tempfile(fileext=".chain"); con<-gzfile(path,"rt"); writeLines(readLines(con,warn=FALSE),plain); close(con); path<-plain }
  import.chain(path)
}
chain <- resolve_chain(opt$chain, opt$outdir)
lifted <- liftOver(gl, chain)
np <- elementNROWS(lifted)
step("mapped by liftOver", sum(np>0))

fl <- unlist(lifted, use.names=FALSE); grp <- rep(seq_along(lifted), np); n_in <- length(lifted)
fchr<-as.character(seqnames(fl)); fst<-start(fl); fen<-end(fl)
g_nchr<-tapply(fchr,grp,function(z)length(unique(z))); g_chr<-tapply(fchr,grp,function(z)z[1])
g_min<-tapply(fst,grp,min); g_max<-tapply(fen,grp,max); g_cov<-tapply(fen-fst+1L,grp,sum)
present<-as.integer(names(g_nchr))
mk<-function(v){o<-rep(NA_real_,n_in);o[present]<-as.numeric(v);o}
v_nchr<-mk(g_nchr); v_chr<-{o<-rep(NA_character_,n_in);o[present]<-as.character(g_chr);o}
v_min<-mk(g_min); v_max<-mk(g_max); v_cov<-mk(g_cov)
span_w<-v_max-v_min+1; ratio<-span_w/orig_w; cov_frac<-v_cov/span_w
keep <- !is.na(v_nchr) & v_nchr==1 & ratio>=opt$size_lo & ratio<=opt$size_hi & cov_frac>=opt$min_cov_frac
say(sprintf("[LIFTOVER] dropped: unmapped=%d, multi-chrom=%d, size-out=%d, low-cov=%d.",
            sum(np==0), sum(!is.na(v_nchr)&v_nchr>1),
            sum(!is.na(ratio)&(ratio<opt$size_lo|ratio>opt$size_hi)),
            sum(!is.na(cov_frac)&cov_frac<opt$min_cov_frac & v_nchr==1)))
gl19 <- GRanges(v_chr[keep], IRanges(as.integer(v_min[keep]), as.integer(v_max[keep])))
step("after liftOver QC (hg19)", length(gl19))

# 4. canonical hg19 + (optional) merge --------------------------------------
gl19 <- gl19[as.character(seqnames(gl19)) %in% CANON]
gl19 <- keepSeqlevels(gl19, intersect(CANON, seqlevelsInUse(gl19)), pruning.mode="coarse")
gl19 <- sort(gl19)
if (!opt$no_merge) { gl19 <- reduce(gl19); step("after merge overlaps", length(gl19)) } else say("[MERGE] skipped (--no_merge).")

# 5. (optional) blacklist ---------------------------------------------------
if (!is.null(opt$blacklist)) {
  if (!file.exists(opt$blacklist)) stop("--blacklist not found: ", opt$blacklist)
  bl <- granges(import(opt$blacklist))
  gl19 <- gl19[!overlapsAny(gl19, bl)]
  step("after blacklist filter", length(gl19))
} else say("[BLACKLIST] none supplied — skipped.")

# 6. (optional) pad narrow regions (Ryder ignores <100 bp) ------------------
if (opt$min_width > 0L) {
  narrow <- width(gl19) < opt$min_width
  if (any(narrow)) {
    gl19[narrow] <- resize(gl19[narrow], width=opt$min_width, fix="center")
    start(gl19)[start(gl19) < 1L] <- 1L                # clamp left edge
    gl19 <- reduce(sort(gl19))                         # re-merge any new overlaps from padding
  }
  step(sprintf("after pad-to >=%d bp", opt$min_width), length(gl19))
}

# 7. write ------------------------------------------------------------------
gl19 <- sort(gl19)
if (!length(gl19)) stop("No regions survived; loosen --size_lo/--size_hi or check chain/input.")
bed <- data.frame(chrom=as.character(seqnames(gl19)), start=start(gl19)-1L, end=end(gl19))
write.table(bed, out_bed, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

git_commit <- tryCatch(system("git rev-parse --short HEAD 2>/dev/null", intern=TRUE), error=function(e) character(0))
report <- c(strrep("=",70), "liftover_ctcf_hg19.R — provenance & funnel", strrep("=",70), "",
  sprintf("timestamp : %s", format(Sys.time(),"%Y-%m-%d %H:%M:%S")),
  sprintf("input bed : %s", normalizePath(opt$bed)),
  sprintf("chain     : %s", if (is.null(opt$chain)) "UCSC hg38ToHg19 (auto)" else opt$chain),
  sprintf("blacklist : %s", if (is.null(opt$blacklist)) "none" else opt$blacklist),
  sprintf("params    : size_lo=%.2f size_hi=%.2f min_cov_frac=%.2f min_width=%d merge=%s",
          opt$size_lo, opt$size_hi, opt$min_cov_frac, opt$min_width, !opt$no_merge),
  sprintf("git       : %s", if (length(git_commit)) git_commit else "NA"),
  sprintf("output    : %s", normalizePath(out_bed)),
  "", "region survival funnel:", funnel, "",
  sprintf("FINAL : %d regions, %d bp, median width %d bp",
          length(gl19), sum(width(gl19)), as.integer(median(width(gl19)))),
  "", "--- log ---", LOG)
writeLines(report, out_report)
say(sprintf("\n[DONE] %d hg19 CTCF regions -> %s", length(gl19), out_bed))
say(sprintf("[DONE] report -> %s", out_report))