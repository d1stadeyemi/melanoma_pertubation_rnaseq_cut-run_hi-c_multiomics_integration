#!/usr/bin/env Rscript
# =============================================================================
# prepare_greenlist_hg19.R
# -----------------------------------------------------------------------------
# Turn the published de Mello et al. (2024) human CUT&RUN GREENLIST (hg38, in an
# Excel sheet) into a CLEAN hg19 BED that drops straight into csaw.R --greenlist.
#
# WHY EACH STEP EXISTS ---------------------------------------------------------
#   The greenlist is a set of gene-distant, low-signal "biology-free" background
#   regions whose read counts track TECHNICAL efficiency. csaw normalizes by
#   computing TMM on counts in these regions. We are NOT regenerating the list
#   (that needs hundreds of IgG samples) — we reuse the validated list and only
#   (1) move it to your genome build and (2) re-verify its integrity afterwards.
#
#   FULLY AUTOMATIC: give it the Excel file; it does the rest. The chain file is
#   downloaded if absent. The IgG background check and the gene-distance check are
#   OPTIONAL and only run if the relevant input is available.
#
# THE PIPELINE (each stage reports how many regions survive) -------------------
#   1. Read Excel  -> hg38 GRanges. Excel Start is treated as 0-based (BED): the
#      exact 1000 bp Start/End differences in the published table are the tell of
#      0-based, end-exclusive coordinates (a 1-based 1 kb bin would differ by 999).
#      Override with --input_base 1 if your copy is 1-based.
#   2. Keep canonical chromosomes only (chr1-22, X, Y) — drop alt/random/Un/M.
#   3. liftOver hg38 -> hg19 (rtracklayer + UCSC hg38ToHg19 chain). Per region:
#        - drop regions that fail to map,
#        - drop regions that map across >1 chromosome (rearranged = doubtful),
#        - keep only regions whose lifted SPAN width stayed close to the ORIGINAL
#          width (--size_lo/--size_hi; RELATIVE to each region, NOT an absolute
#          1 kb — the published regions are merged 1 kb bins of 1-4 kb), and whose
#          lifted pieces cover most of that span (--min_cov_frac; no gap-spanning).
#   4. Keep canonical hg19 chromosomes; collapse any overlaps (reduce) so reads
#      are never double-counted by csaw's regionCounts.
#   5. (optional) Drop regions overlapping an hg19 ENCODE blacklist (--blacklist).
#   6. (optional) Re-impose de Mello's gene-distance rule in hg19: drop regions
#      within --gene_dist bp of a gene (--genes BED/GTF, else TxDb.Hsapiens.UCSC
#      .hg19.knownGene if installed, else skipped with a notice — the published
#      list was already gene-distant in hg38, so this is a safety net).
#   7. (optional) IgG background check (--igg BAM[s]): count IgG fragments in each
#      greenlist region and in genome-wide bins; drop any greenlist region whose
#      IgG density exceeds --igg_fold x the genome-median background (i.e. not
#      background-level in YOUR A375 system). A single representative A375 IgG is
#      adequate here because this tests a region property, not per-sample scale.
#   8. Write the hg19 BED (0-based, 3 columns) + a provenance/funnel report.
#
# OUTPUT -----------------------------------------------------------------------
#   <outdir>/greenlist_hg19.bed         (feed to csaw.R --greenlist)
#   <outdir>/greenlist_hg19.report.txt  (provenance + the survival funnel)
#
# EXAMPLE RUNS -----------------------------------------------------------------
#   # Minimal (conversion + structural cleanup only; chain auto-downloaded):
#   Rscript prepare_greenlist_hg19.R --excel greenlist_hg38.xlsx --outdir data
#
#   # Full clean with the A375 IgG background check and a blacklist:
#   Rscript prepare_greenlist_hg19.R --excel greenlist_hg38.xlsx --outdir data \
#       --igg ../../Taiwo/20231024_Cut_Run/IgG-R1/aligned/dedup/IgG-R1.bam \
#       --blacklist data/hg19-blacklist.v2.bed
#
#   # If compute nodes lack internet, fetch the chain on a login node first:
#   #   wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
#   # then: --chain hg38ToHg19.over.chain.gz
# =============================================================================

suppressWarnings(suppressMessages({
  library(optparse)
  library(csaw)
  library(rtracklayer)     # import.chain / liftOver / import (GTF/BED)
  library(GenomicRanges)
}))

# ---- CLI --------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list = list(
  make_option("--excel", type = "character",
              help = "REQUIRED. .xlsx containing the hg38 greenlist (cols: Chromosome, Start, End)."),
  make_option("--sheet_name", type = "character", default = "1",
              help = "Worksheet index or name holding the greenlist. [default %default]"),
  make_option("--input_base", type = "integer", default = 0L,
              help = "Coordinate base of the Excel Start column: 0 (0-based/BED) or 1 (1-based). [default %default]"),
  make_option("--outdir", type = "character", default = ".",
              help = "Output directory. [default %default]"),
  make_option("--out", type = "character", default = NULL,
              help = "Output BED path. [default <outdir>/greenlist_hg19.bed]"),
  make_option("--chain", type = "character", default = NULL,
              help = "Path to hg38ToHg19.over.chain(.gz). If absent, downloaded from UCSC. [default: auto]"),
  make_option("--size_lo", type = "double", default = 0.8,
              help = "Min lifted-span/original width ratio to keep a region. [default %default]"),
  make_option("--size_hi", type = "double", default = 1.25,
              help = "Max lifted-span/original width ratio to keep a region. [default %default]"),
  make_option("--min_cov_frac", type = "double", default = 0.8,
              help = "Min (sum of lifted piece widths)/(lifted span) — guards against gap-spanning. [default %default]"),
  make_option("--blacklist", type = "character", default = NULL,
              help = "Optional hg19 ENCODE blacklist BED; greenlist regions overlapping it are dropped. [default: none]"),
  make_option("--genes", type = "character", default = NULL,
              help = "Optional hg19 gene annotation (BED or GTF/GFF) for the gene-distance filter. If omitted, TxDb.Hsapiens.UCSC.hg19.knownGene is used when installed; else the filter is skipped. [default: none]"),
  make_option("--gene_dist", type = "integer", default = 5000L,
              help = "Min distance (bp) a greenlist region must keep from any gene (de Mello rule). [default %default]"),
  make_option("--igg", type = "character", default = NULL,
              help = "Optional IgG BAM path(s), comma-separated. Background check only: drops greenlist regions that are not background-level in your IgG. [default: none]"),
  make_option("--igg_fold", type = "double", default = 5.0,
              help = "Drop a greenlist region if its IgG density > this multiple of the genome-median background density. [default %default]"),
  make_option("--frag_max", type = "integer", default = 700L,
              help = "Max fragment length (bp) for IgG counting; match your alignment. [default %default]"),
  make_option("--bg_bin", type = "integer", default = 10000L,
              help = "Genome bin width (bp) for the IgG background-density baseline. [default %default]")
)))
if (is.null(opt$excel)) stop("--excel is required.")
if (!opt$input_base %in% c(0L, 1L)) stop("--input_base must be 0 or 1.")
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
out_bed    <- if (is.null(opt$out)) file.path(opt$outdir, "greenlist_hg19.bed") else opt$out
out_report <- sub("\\.bed$", ".report.txt", out_bed); if (out_report == out_bed) out_report <- paste0(out_bed, ".report.txt")

CANON <- paste0("chr", c(1:22, "X", "Y"))   # canonical chromosomes we allow

# ---- logging (echo + collect for the report) --------------------------------
LOG <- character(0)
say <- function(...) { l <- paste0(...); cat(l, "\n"); LOG <<- c(LOG, l) }
funnel <- character(0)
step <- function(label, n) { line <- sprintf("  %-46s %7d", label, n); cat(line, "\n"); funnel <<- c(funnel, line) }

say("=== prepare_greenlist_hg19.R ===")

# =============================================================================
# 1. Read the Excel greenlist -> hg38 GRanges (1-based internally)
# =============================================================================
read_excel_any <- function(path, sheet) {
  sh <- suppressWarnings(if (grepl("^[0-9]+$", sheet)) as.integer(sheet) else sheet)
  if (requireNamespace("readxl", quietly = TRUE)) {
    as.data.frame(readxl::read_excel(path, sheet = sh))
  } else if (requireNamespace("openxlsx", quietly = TRUE)) {
    openxlsx::read.xlsx(path, sheet = sh)
  } else stop("Need 'readxl' or 'openxlsx' to read .xlsx. Install one: install.packages('readxl').")
}

raw <- read_excel_any(opt$excel, opt$sheet_name)
if (ncol(raw) < 3) stop("Expected >=3 columns (Chromosome, Start, End); got ", ncol(raw), ".")
raw <- raw[, 1:3]; colnames(raw) <- c("chrom", "start", "end")
raw$chrom <- as.character(raw$chrom)
raw$start <- as.numeric(raw$start)
raw$end   <- as.numeric(raw$end)
raw <- raw[!is.na(raw$chrom) & !is.na(raw$start) & !is.na(raw$end), , drop = FALSE]
# Excel 0-based BED -> 1-based GRanges start; 1-based input -> start unchanged.
start1 <- if (opt$input_base == 0L) as.integer(raw$start) + 1L else as.integer(raw$start)
gl_hg38 <- GRanges(seqnames = raw$chrom, ranges = IRanges(start = start1, end = as.integer(raw$end)))
suppressWarnings(try(seqlevelsStyle(gl_hg38) <- "UCSC", silent = TRUE))   # ensure 'chr' prefix
say(sprintf("[READ] %s sheet '%s' (input_base=%d): %d regions; median width %d bp.",
            basename(opt$excel), opt$sheet_name, opt$input_base, length(gl_hg38), as.integer(median(width(gl_hg38)))))
say("[FUNNEL] region survival:")
step("input (hg38)", length(gl_hg38))

# =============================================================================
# 2. Canonical chromosomes only (hg38 side)
# =============================================================================
gl_hg38 <- gl_hg38[as.character(seqnames(gl_hg38)) %in% CANON]
gl_hg38 <- keepSeqlevels(gl_hg38, intersect(CANON, seqlevelsInUse(gl_hg38)), pruning.mode = "coarse")
step("after canonical-chrom (hg38)", length(gl_hg38))
orig_w <- width(gl_hg38)                                # remember each region's hg38 width

# =============================================================================
# 3. Resolve the chain file (download + gunzip if needed) and liftOver
# =============================================================================
resolve_chain <- function(chain_arg, outdir) {
  url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
  path <- chain_arg
  if (is.null(path)) {
    path <- file.path(outdir, "hg38ToHg19.over.chain.gz")
    if (!file.exists(path)) {
      say(sprintf("[CHAIN] downloading %s", url))
      ok <- tryCatch({ download.file(url, path, mode = "wb", quiet = TRUE); TRUE }, error = function(e) FALSE)
      if (!ok || !file.exists(path))
        stop("Chain download failed (compute nodes may lack internet). Fetch it on a login node and pass --chain:\n  wget ", url)
    }
  }
  if (!file.exists(path)) stop("--chain not found: ", path)
  if (grepl("\\.gz$", path)) {                          # import.chain needs a plain text file
    plain <- tempfile(fileext = ".chain")
    con <- gzfile(path, "rt"); writeLines(readLines(con, warn = FALSE), plain); close(con)
    path <- plain
  }
  import.chain(path)
}
chain <- resolve_chain(opt$chain, opt$outdir)

lifted <- liftOver(gl_hg38, chain)                      # GRangesList, one element per input region
np <- elementNROWS(lifted)                              # pieces per region (0 = unmapped)
step("mapped by liftOver", sum(np > 0))

# Per-region aggregates via a flat representation (fast for large lists).
fl  <- unlist(lifted, use.names = FALSE)
grp <- rep(seq_along(lifted), np)
n_in    <- length(lifted)
fchr <- as.character(seqnames(fl)); fst <- start(fl); fen <- end(fl)
g_nchr <- tapply(fchr, grp, function(z) length(unique(z)))   # distinct chromosomes hit
g_chr  <- tapply(fchr, grp, function(z) z[1])
g_min  <- tapply(fst,  grp, min)
g_max  <- tapply(fen,  grp, max)
g_cov  <- tapply(fen - fst + 1L, grp, sum)                   # summed piece widths
present <- as.integer(names(g_nchr))                        # group ids that mapped

# Re-expand to full input length (NA where unmapped).
mk <- function(v) { o <- rep(NA_real_, n_in); o[present] <- as.numeric(v); o }
v_nchr <- mk(g_nchr); v_chr <- { o <- rep(NA_character_, n_in); o[present] <- as.character(g_chr); o }
v_min  <- mk(g_min);  v_max <- mk(g_max); v_cov <- mk(g_cov)
span_w  <- v_max - v_min + 1
ratio   <- span_w / orig_w                                  # lifted span vs original width
cov_frac<- v_cov / span_w

keep <- !is.na(v_nchr) & v_nchr == 1 &                      # single chromosome
        ratio >= opt$size_lo & ratio <= opt$size_hi &       # width preserved (relative)
        cov_frac >= opt$min_cov_frac                        # not spanning a big gap
say(sprintf("[LIFTOVER] dropped: unmapped=%d, multi-chrom=%d, size-out-of-bounds=%d, low-coverage=%d.",
            sum(np == 0),
            sum(!is.na(v_nchr) & v_nchr > 1),
            sum(!is.na(ratio) & (ratio < opt$size_lo | ratio > opt$size_hi)),
            sum(!is.na(cov_frac) & cov_frac < opt$min_cov_frac & v_nchr == 1)))
gl_hg19 <- GRanges(seqnames = v_chr[keep],
                   ranges = IRanges(start = as.integer(v_min[keep]), end = as.integer(v_max[keep])))
step("after liftOver QC (hg19)", length(gl_hg19))

# =============================================================================
# 4. Canonical hg19 chromosomes + collapse overlaps (no double counting)
# =============================================================================
gl_hg19 <- gl_hg19[as.character(seqnames(gl_hg19)) %in% CANON]
gl_hg19 <- keepSeqlevels(gl_hg19, intersect(CANON, seqlevelsInUse(gl_hg19)), pruning.mode = "coarse")
gl_hg19 <- reduce(sort(gl_hg19))
step("after canonical + merge overlaps", length(gl_hg19))

# =============================================================================
# 5. (optional) ENCODE blacklist filter
# =============================================================================
read_bed3 <- function(path) {
  d <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", quote = "")
  GRanges(d[[1]], IRanges(as.integer(d[[2]]) + 1L, as.integer(d[[3]])))
}
if (!is.null(opt$blacklist)) {
  if (!file.exists(opt$blacklist)) stop("--blacklist not found: ", opt$blacklist)
  bl <- read_bed3(opt$blacklist)
  gl_hg19 <- gl_hg19[!overlapsAny(gl_hg19, bl)]
  step("after blacklist filter", length(gl_hg19))
} else say("[BLACKLIST] none supplied — skipped.")

# =============================================================================
# 6. (optional) Gene-distance filter (>= gene_dist bp from any gene)
# =============================================================================
load_genes <- function(opt) {
  if (!is.null(opt$genes)) {
    if (!file.exists(opt$genes)) stop("--genes not found: ", opt$genes)
    if (grepl("\\.(gtf|gff|gff3)(\\.gz)?$", opt$genes, ignore.case = TRUE)) {
      g <- import(opt$genes); g <- g[tolower(as.character(g$type)) == "gene"]
      if (!length(g)) g <- import(opt$genes)             # fall back to all features
      return(g)
    }
    return(read_bed3(opt$genes))
  }
  if (requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE)) {
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    return(suppressWarnings(GenomicFeatures::genes(txdb)))
  }
  NULL
}
genes <- load_genes(opt)
if (!is.null(genes)) {
  suppressWarnings(try(seqlevelsStyle(genes) <- "UCSC", silent = TRUE))
  d2g <- distanceToNearest(gl_hg19, genes, ignore.strand = TRUE)
  dist_vec <- rep(Inf, length(gl_hg19))                  # regions with no gene on their chrom = far
  dist_vec[queryHits(d2g)] <- mcols(d2g)$distance
  gl_hg19 <- gl_hg19[dist_vec >= opt$gene_dist]
  step(sprintf("after gene-distance (>=%d bp)", opt$gene_dist), length(gl_hg19))
} else {
  say("[GENES] no --genes and TxDb.Hsapiens.UCSC.hg19.knownGene not installed — gene-distance filter skipped (published list was already gene-distant in hg38).")
}

# =============================================================================
# 7. (optional) IgG background check
# =============================================================================
if (!is.null(opt$igg)) {
  igg_bams <- trimws(strsplit(opt$igg, ",")[[1]])
  miss <- igg_bams[!file.exists(igg_bams)]; if (length(miss)) stop("IgG BAM(s) not found: ", paste(miss, collapse = ", "))
  if (!requireNamespace("csaw", quietly = TRUE)) stop("--igg requires the 'csaw' package for fragment counting.")
  param <- csaw::readParam(minq = 20L, dedup = TRUE, pe = "both", max.frag = opt$frag_max)
  say(sprintf("[IGG] background check using %d BAM(s): %s", length(igg_bams), paste(basename(igg_bams), collapse = ", ")))

  reg_cnt <- csaw::regionCounts(igg_bams, regions = gl_hg19, param = param)   # IgG fragments per greenlist region
  reg_density <- rowSums(SummarizedExperiment::assay(reg_cnt)) / width(gl_hg19)

  bins <- csaw::windowCounts(igg_bams, bin = TRUE, width = opt$bg_bin, param = param)   # genome background
  bin_density <- rowSums(SummarizedExperiment::assay(bins)) / opt$bg_bin
  bg_med <- median(bin_density[bin_density > 0], na.rm = TRUE)                # genome-median background density
  thresh <- opt$igg_fold * bg_med

  drop_hi <- reg_density > thresh
  say(sprintf("[IGG] genome-median background density = %.4g frag/bp; threshold = %.1fx = %.4g; dropping %d region(s) above it.",
              bg_med, opt$igg_fold, thresh, sum(drop_hi)))
  if (sum(drop_hi) > 0.2 * length(gl_hg19))
    say("[IGG] WARNING: >20% of regions flagged — inspect the liftOver/IgG; this is a lot to remove.")
  gl_hg19 <- gl_hg19[!drop_hi]
  step("after IgG background check", length(gl_hg19))
} else say("[IGG] none supplied — background check skipped (list still valid).")

# =============================================================================
# 8. Write the hg19 BED (0-based, 3 columns) + provenance report
# =============================================================================
gl_hg19 <- reduce(sort(gl_hg19))
if (!length(gl_hg19)) stop("No regions survived. Loosen --size_lo/--size_hi or check the input/chain.")
bed <- data.frame(chrom = as.character(seqnames(gl_hg19)),
                  start = start(gl_hg19) - 1L,           # GRanges 1-based -> BED 0-based
                  end   = end(gl_hg19))
write.table(bed, out_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

git_commit <- tryCatch(system("git rev-parse --short HEAD 2>/dev/null", intern = TRUE), error = function(e) character(0))
git_dirty  <- tryCatch(length(system("git status --porcelain 2>/dev/null", intern = TRUE)) > 0, error = function(e) NA)
report <- c(
  strrep("=", 70), "prepare_greenlist_hg19.R — provenance & funnel", strrep("=", 70), "",
  sprintf("timestamp     : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("input excel   : %s (sheet '%s', input_base=%d)", normalizePath(opt$excel), opt$sheet_name, opt$input_base),
  sprintf("chain         : %s", if (is.null(opt$chain)) "UCSC hg38ToHg19 (auto)" else opt$chain),
  sprintf("blacklist     : %s", if (is.null(opt$blacklist)) "none" else opt$blacklist),
  sprintf("genes         : %s", if (!is.null(opt$genes)) opt$genes else if (!is.null(genes)) "TxDb.Hsapiens.UCSC.hg19.knownGene" else "skipped"),
  sprintf("igg           : %s", if (is.null(opt$igg)) "none" else opt$igg),
  sprintf("params        : size_lo=%.2f size_hi=%.2f min_cov_frac=%.2f gene_dist=%d igg_fold=%.1f frag_max=%d bg_bin=%d",
          opt$size_lo, opt$size_hi, opt$min_cov_frac, opt$gene_dist, opt$igg_fold, opt$frag_max, opt$bg_bin),
  sprintf("git           : %s%s", if (length(git_commit)) git_commit else "NA", if (isTRUE(git_dirty)) " (dirty)" else ""),
  sprintf("output BED    : %s", normalizePath(out_bed)),
  "", "region survival funnel:", funnel, "",
  sprintf("FINAL         : %d regions, %d total bp, median width %d bp",
          length(gl_hg19), sum(width(gl_hg19)), as.integer(median(width(gl_hg19)))),
  "", "--- full log ---", LOG
)
writeLines(report, out_report)

say(sprintf("\n[DONE] %d hg19 greenlist regions -> %s", length(gl_hg19), out_bed))
say(sprintf("[DONE] report -> %s", out_report))
say(sprintf("[NEXT] csaw.R ... --norm libsize,background,greenlist --greenlist %s", out_bed))