#!/usr/bin/env Rscript
# =============================================================================
# greenlist_deseq2_diffbind.R
#   de Mello et al. 2024 greenlist normalization + DESeq2 differential binding,
#   ADAPTED for a BROAD histone mark (H3K27me3), as an INDEPENDENT cross-check
#   of the csaw (edgeR + TMM-on-greenlist) differential-binding run.
#
# THE PAPER'S METHOD, EXACTLY -------------------------------------------------
#   1. Quantify read pile-up in the greenlist regions per sample.
#   2. Compute DESeq2 size factors on those greenlist counts
#        (DESeq2::estimateSizeFactorsForMatrix = median-of-ratios). Exact fn.
#   3. INJECT those size factors into a DESeq2 differential-binding fit over the
#        target signal, OVERRIDING DESeq2's default (which would estimate size
#        factors from the target counts themselves and so re-introduce the
#        "most regions unchanged" assumption that fails under a global shift).
#   4. DESeq() (Wald) -> results(contrast) -> per-region log2FC + BH padj.
#
#   The only deviation from the punctate-TF example in the paper is the COUNTING
#   UNIT: a broad mark is not summarised by narrow peaks. Two units, each in its
#   own sub-directory:
#     peaks   : consensus BROAD peaks (union of the reference reps' Peaks column)
#               -> region-level, directly comparable to csaw's in_peak regions.
#     windows : contiguous NON-overlapping bins (width = --window), enrichment-
#               filtered over background -> genome-wide; parallels the csaw window
#               run (DESeq2 needs independent features, so bins are contiguous,
#               NOT the sliding/merged windows edgeR/csaw uses).
#
# CONSISTENCY -----------------------------------------------------------------
#   Greenlist AND target are counted with the SAME readParam. By default the
#   greenlist size factors are recomputed here over EXACTLY the analysis samples
#   (so the median-of-ratios reference matches the dds). Optionally inject the
#   multiBamSummary factors from greenlist_deseq2_sizefactors.py via
#   --size_factors (a consistency check vs the internal ones is printed).
#
# DROP THE OUTLIER ------------------------------------------------------------
#   --drop Parental-H3K27me3-R3   (comma-separated SampleIDs) excludes from BOTH
#   the size-factor estimation and the fit.
#
# EXAMPLE ---------------------------------------------------------------------
#   Rscript greenlist_deseq2_diffbind.R \
#     --sheet data/processed/sample_sheet.csv --mark H3K27me3 \
#     --greenlist data/processed/greenlist/greenlist_hg19.bed \
#     --group1 "Parental|UNC" --group2 "Parental|none" \
#     --drop Parental-H3K27me3-R3 \
#     --unit peaks,windows --peaks sheet \
#     --blacklist data/processed/greenlist/hg19-blacklist.v2.bed \
#     --outdir results/diff_deseq2/H3K27me3_ParentalUNC_vs_Parentalnone
#
# OUTPUTS (per unit, in <outdir>/<unit>/) -------------------------------------
#   <tag>_<unit>_all_regions.csv / _significant.csv
#   <tag>_<unit>_MA.pdf / _volcano.pdf
#   <tag>_<unit>_sizefactors_used.csv
#   <tag>_<unit>_summary.txt / _MANIFEST.txt
# =============================================================================

suppressWarnings(suppressMessages({
  library(optparse)
  library(csaw)
  library(DESeq2)
  library(GenomicRanges)
  library(SummarizedExperiment)
}))

parse_cli <- function() {
  ol <- list(
    make_option("--sheet", type = "character", help = "REQUIRED. Sample sheet CSV."),
    make_option("--mark", type = "character", help = "REQUIRED. Factor, case-insensitive (e.g. H3K27me3)."),
    make_option("--greenlist", type = "character", help = "REQUIRED. Greenlist BED (hg19), 3-column."),
    make_option("--group1", type = "character", default = "Parental|UNC", help = "Numerator 'Condition|Treatment' (treated). [%default]"),
    make_option("--group2", type = "character", default = "Parental|none", help = "Reference 'Condition|Treatment'. +log2FC = higher in group1. [%default]"),
    make_option("--drop", type = "character", default = NULL, help = "Comma-separated SampleIDs to exclude. [none]"),
    make_option("--unit", type = "character", default = "peaks,windows", help = "Counting unit(s), comma list: peaks | windows. Each -> own sub-dir. [%default]"),
    make_option("--peaks", type = "character", default = "sheet", help = "[unit=peaks] 'sheet' (union Peaks of group2 reference rows) OR comma BED/broadPeak paths. [%default]"),
    make_option("--size_factors", type = "character", default = NULL, help = "Optional _size_factors.tsv from greenlist_deseq2_sizefactors.py to inject instead (consistency check printed). [recompute internally]"),
    make_option("--outdir", type = "character", default = "deseq2_diffbind", help = "Base output dir; <outdir>/<unit>/ per unit. [%default]"),
    make_option("--window", type = "integer", default = 2000L, help = "[unit=windows] contiguous bin width (bp). [%default]"),
    make_option("--bg_bin", type = "integer", default = 10000L, help = "Background bin width (bp) for the enrichment filter. [%default]"),
    make_option("--min_lfc_bg", type = "double", default = log2(3), help = "[unit=windows] min log2 enrichment over background to keep a bin. [log2(3)]"),
    make_option("--frag_max", type = "integer", default = 700L, help = "Max fragment length (bp); match Bowtie2 -X. [%default]"),
    make_option("--blacklist", type = "character", default = NULL, help = "Optional ENCODE blacklist BED. [none]"),
    make_option("--fdr", type = "double", default = 0.05, help = "padj cutoff for the significant subset. [%default]"),
    make_option("--lfc", type = "double", default = 1.0, help = "min |log2FC| for the significant subset. [%default]")
  )
  opt <- parse_args(OptionParser(option_list = ol))
  for (r in c("sheet", "mark", "greenlist"))
    if (is.null(opt[[r]])) stop(sprintf("--%s is required.", r))
  opt
}

parse_group <- function(spec) {
  p <- strsplit(spec, "\\|")[[1]]
  if (length(p) != 2) stop("Group spec must be 'Condition|Treatment', got: ", spec)
  list(cond = tolower(trimws(p[1])), treat = tolower(trimws(p[2])))
}
safe_token <- function(s) gsub("[^A-Za-z0-9]+", "", s)

read_bed3 <- function(path) {
  df <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", quote = "")
  GRanges(df[[1]], IRanges(as.integer(df[[2]]) + 1L, as.integer(df[[3]])))
}

log_lines <- new.env(); log_lines$txt <- character(0)
say <- function(...) { l <- paste0(...); cat(l, "\n"); log_lines$txt <- c(log_lines$txt, l) }
git_hash <- function() tryCatch(system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE)[1],
                                error = function(e) "NA")

load_sheet <- function(path) {
  ss <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  req <- c("SampleID", "Factor", "Condition", "Treatment", "Replicate", "Batch", "bamReads")
  miss <- setdiff(req, colnames(ss)); if (length(miss)) stop("Sheet missing: ", paste(miss, collapse = ", "))
  ss$.cond <- tolower(trimws(ss$Condition)); ss$.treat <- tolower(trimws(ss$Treatment)); ss$.fac <- tolower(trimws(ss$Factor))
  ss
}

select_samples <- function(ss, opt) {
  ss <- ss[ss$.fac == tolower(trimws(opt$mark)), , drop = FALSE]
  if (!nrow(ss)) stop("No samples for mark '", opt$mark, "'.")
  g1 <- parse_group(opt$group1); g2 <- parse_group(opt$group2)
  k <- function(g) paste(g$cond, g$treat, sep = "|")
  ss$.key <- paste(ss$.cond, ss$.treat, sep = "|")
  sub <- ss[ss$.key %in% c(k(g1), k(g2)), , drop = FALSE]
  if (!is.null(opt$drop)) {
    drop <- trimws(strsplit(opt$drop, ",")[[1]])
    hit  <- intersect(drop, sub$SampleID)
    if (length(hit)) say(sprintf("[DROP] excluding: %s", paste(hit, collapse = ", ")))
    miss <- setdiff(drop, sub$SampleID)
    if (length(miss)) say(sprintf("[DROP] note: not in this contrast (ignored): %s", paste(miss, collapse = ", ")))
    sub <- sub[!sub$SampleID %in% drop, , drop = FALSE]
  }
  n1 <- sum(sub$.key == k(g1)); n2 <- sum(sub$.key == k(g2))
  if (n1 < 2 || n2 < 2) stop(sprintf("Need >=2 reps per group after --drop (group1=%d, group2=%d).", n1, n2))
  sub$.role <- ifelse(sub$.key == k(g1), "treat", "ref")
  if (any(!file.exists(sub$bamReads))) stop("Missing BAM(s): ", paste(sub$bamReads[!file.exists(sub$bamReads)], collapse = ", "))
  attr(sub, "g1key") <- k(g1); attr(sub, "g2key") <- k(g2)
  sub
}

make_param <- function(opt) {
  discard <- GRanges()
  if (!is.null(opt$blacklist)) { if (!file.exists(opt$blacklist)) stop("Blacklist not found."); discard <- read_bed3(opt$blacklist) }
  readParam(minq = 20L, dedup = TRUE, pe = "both", max.frag = opt$frag_max, discard = discard)
}

resolve_peaks <- function(opt, meta) {
  if (identical(opt$peaks, "sheet")) {
    if (!"Peaks" %in% colnames(meta)) stop("--peaks sheet but no 'Peaks' column.")
    paths <- unique(meta$Peaks[meta$.role == "ref"]); paths <- paths[!is.na(paths) & nzchar(paths)]
    if (!length(paths)) stop("--peaks sheet: no Peaks on the reference rows.")
  } else paths <- trimws(strsplit(opt$peaks, ",")[[1]])
  miss <- paths[!file.exists(paths)]; if (length(miss)) stop("Peak file(s) not found: ", paste(miss, collapse = ", "))
  gr <- reduce(Reduce(c, lapply(paths, read_bed3)))
  say(sprintf("[PEAKS] %d file(s) -> %d consensus regions", length(paths), length(gr)))
  gr
}

greenlist_size_factors <- function(meta, green_gr, param) {
  gc <- regionCounts(meta$bamReads, regions = green_gr, param = param)
  m  <- assay(gc); colnames(m) <- meta$SampleID
  used <- is.finite(rowMeans(log(m)))
  if (sum(used) < 50) say(sprintf("[GREENLIST] WARNING: only %d usable (zero-free) regions -> noisy size factors.", sum(used)))
  sf <- estimateSizeFactorsForMatrix(m)
  say(sprintf("[GREENLIST] %d regions counted, %d usable; size factors: %s",
              nrow(m), sum(used), paste(sprintf("%s=%.3f", meta$SampleID, sf), collapse = ", ")))
  sf
}

read_external_sf <- function(path, meta) {
  tab <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("SampleID", "deseq2_size_factor") %in% colnames(tab)))
    stop("--size_factors needs columns SampleID, deseq2_size_factor.")
  v <- setNames(tab$deseq2_size_factor, tab$SampleID)[meta$SampleID]
  if (any(is.na(v))) stop("--size_factors missing for: ", paste(meta$SampleID[is.na(v)], collapse = ", "))
  v <- v / exp(mean(log(v)))
  setNames(as.numeric(v), meta$SampleID)
}

unit_counts <- function(unit, meta, opt, param) {
  if (unit == "peaks") {
    gr <- resolve_peaks(opt, meta)
    rc <- regionCounts(meta$bamReads, regions = gr, param = param)
    keep <- rowSums(assay(rc)) > 0
    rc <- rc[keep, ]; gr <- rowRanges(rc)
    say(sprintf("[peaks] %d consensus peaks with signal.", length(gr)))
  } else if (unit == "windows") {
    win  <- windowCounts(meta$bamReads, width = opt$window, spacing = opt$window, param = param)
    bins <- windowCounts(meta$bamReads, bin = TRUE, width = opt$bg_bin, param = param)
    fstat <- filterWindowsGlobal(win, bins)
    keep  <- fstat$filter > opt$min_lfc_bg
    say(sprintf("[windows] kept %d / %d contiguous %dbp bins (>%.2f log2 over background).",
                sum(keep), length(keep), opt$window, opt$min_lfc_bg))
    rc <- win[keep, ]; gr <- rowRanges(rc)
  } else stop("Unknown unit: ", unit)
  m <- round(assay(rc)); storage.mode(m) <- "integer"
  colnames(m) <- meta$SampleID
  list(gr = gr, counts = m)
}

deseq2_db <- function(uc, meta, sf, g1key, g2key, opt) {
  coldata <- data.frame(condition = factor(meta$.role, levels = c("ref", "treat")),
                        row.names = meta$SampleID)
  cm <- uc$counts[, rownames(coldata), drop = FALSE]
  dds <- DESeqDataSetFromMatrix(countData = cm, colData = coldata, design = ~ condition)
  sizeFactors(dds) <- sf[colnames(dds)]
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "treat", "ref"))
  gr <- uc$gr
  df <- data.frame(
    chr = as.character(seqnames(gr)), start = start(gr), end = end(gr),
    baseMean = res$baseMean, logFC = res$log2FoldChange, lfcSE = res$lfcSE,
    stat = res$stat, PValue = res$pvalue, FDR = res$padj, stringsAsFactors = FALSE)
  df$FDR[is.na(df$FDR)] <- 1
  df$call <- ifelse(df$FDR < opt$fdr & abs(df$logFC) > opt$lfc,
                    ifelse(df$logFC > 0, "gained", "lost"), "ns")
  df$direction <- ifelse(df$logFC > 0, "up_in_group1", "down_in_group1")
  list(df = df, dds = dds)
}

ma_volcano <- function(df, tag, unit, outdir, opt) {
  sig <- df$call != "ns"
  pdf(file.path(outdir, sprintf("%s_%s_MA.pdf", tag, unit)), width = 6.5, height = 5)
  plot(log2(df$baseMean + 1), df$logFC, pch = 16, cex = 0.3, col = "#88888855",
       xlab = "log2(baseMean+1)", ylab = "log2FC", main = sprintf("%s %s MA (DESeq2/greenlist)", tag, unit))
  points(log2(df$baseMean[sig] + 1), df$logFC[sig], pch = 16, cex = 0.35,
         col = ifelse(df$logFC[sig] > 0, "#C0392B", "#1F6F8B")); abline(h = 0, col = "grey40")
  dev.off()
  pdf(file.path(outdir, sprintf("%s_%s_volcano.pdf", tag, unit)), width = 6.5, height = 5)
  plot(df$logFC, -log10(df$FDR), pch = 16, cex = 0.3, col = "#88888855",
       xlab = "log2FC", ylab = "-log10 padj", main = sprintf("%s %s volcano", tag, unit))
  points(df$logFC[sig], -log10(df$FDR[sig]), pch = 16, cex = 0.35,
         col = ifelse(df$logFC[sig] > 0, "#C0392B", "#1F6F8B"))
  abline(v = c(-opt$lfc, opt$lfc), h = -log10(opt$fdr), col = "grey60", lty = 2)
  n_up <- sum(sig & df$logFC > 0); n_dn <- sum(sig & df$logFC < 0)
  legend("top", horiz = TRUE, bg = "white", cex = 0.95,
         legend = c(sprintf("lost: %d", n_dn), sprintf("gained: %d", n_up)), text.col = c("#1F6F8B", "#C0392B"))
  dev.off()
}

run_unit <- function(unit, meta, opt, param, sf, tag, base_out) {
  outdir <- file.path(base_out, unit); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  say(sprintf("\n========== UNIT: %s -> %s ==========", unit, outdir))
  uc  <- unit_counts(unit, meta, opt, param)
  fit <- deseq2_db(uc, meta, sf, attr(meta, "g1key"), attr(meta, "g2key"), opt)
  df  <- fit$df
  write.csv(df, file.path(outdir, sprintf("%s_%s_all_regions.csv", tag, unit)), row.names = FALSE)
  sg <- df[df$call != "ns", , drop = FALSE]; sg <- sg[order(sg$FDR), , drop = FALSE]
  write.csv(sg, file.path(outdir, sprintf("%s_%s_significant.csv", tag, unit)), row.names = FALSE)
  n_up <- sum(sg$logFC > 0); n_dn <- sum(sg$logFC < 0)
  say(sprintf("[%s] %d regions tested, %d significant (gained=%d, lost=%d).", unit, nrow(df), nrow(sg), n_up, n_dn))
  ma_volcano(df, tag, unit, outdir, opt)

  sfdf <- data.frame(SampleID = meta$SampleID, role = meta$.role,
                     Condition = meta$Condition, Treatment = meta$Treatment,
                     greenlist_size_factor = sf[meta$SampleID], row.names = NULL)
  write.csv(sfdf, file.path(outdir, sprintf("%s_%s_sizefactors_used.csv", tag, unit)), row.names = FALSE)

  writeLines(c(strrep("=", 64), sprintf("DESeq2/greenlist DB | %s | unit=%s", tag, unit), strrep("=", 64), "",
               sprintf("regions tested = %d ; significant (padj<%.2g, |log2FC|>%.2g) = %d (gained=%d, lost=%d)",
                       nrow(df), opt$fdr, opt$lfc, nrow(sg), n_up, n_dn), "",
               "--- log ---", log_lines$txt),
             file.path(outdir, sprintf("%s_%s_summary.txt", tag, unit)))

  writeLines(c("greenlist_deseq2_diffbind.R",
               sprintf("timestamp: %s", as.character(Sys.time())),
               sprintf("git_commit: %s", git_hash()),
               sprintf("unit: %s", unit),
               sprintf("sheet: %s", opt$sheet), sprintf("mark: %s", opt$mark),
               sprintf("greenlist: %s", opt$greenlist),
               sprintf("group1: %s  group2: %s", opt$group1, opt$group2),
               sprintf("dropped: %s", ifelse(is.null(opt$drop), "none", opt$drop)),
               sprintf("samples: %s", paste(meta$SampleID, collapse = ",")),
               sprintf("size_factor_source: %s", ifelse(is.null(opt$size_factors), "internal (csaw counts + estimateSizeFactorsForMatrix)", opt$size_factors)),
               sprintf("n_regions_tested: %d  n_significant: %d", nrow(df), nrow(sg))),
             file.path(outdir, sprintf("%s_%s_MANIFEST.txt", tag, unit)))
  invisible(df)
}

main <- function() {
  opt <- parse_cli()
  units <- unique(trimws(strsplit(opt$unit, ",")[[1]])); units <- units[nzchar(units)]
  bad <- setdiff(units, c("peaks", "windows")); if (length(bad)) stop("Unknown --unit: ", paste(bad, collapse = ", "))
  if ("peaks" %in% units && is.null(opt$peaks)) stop("--unit peaks needs --peaks.")

  ss   <- load_sheet(opt$sheet)
  meta <- select_samples(ss, opt)
  tag  <- sprintf("%s_%s_vs_%s", opt$mark, safe_token(opt$group1), safe_token(opt$group2))
  say(sprintf("=== DESeq2/greenlist DB | mark=%s | %s vs %s | tag=%s ===", opt$mark, opt$group1, opt$group2, tag))
  say(sprintf("[SAMPLES] %d libraries: %s", nrow(meta), paste(meta$SampleID, collapse = ", ")))

  param    <- make_param(opt)
  green_gr <- reduce(read_bed3(opt$greenlist))
  say(sprintf("[GREENLIST] %d normalizer regions loaded.", length(green_gr)))

  sf_internal <- greenlist_size_factors(meta, green_gr, param)
  if (!is.null(opt$size_factors)) {
    sf_ext <- read_external_sf(opt$size_factors, meta)
    pct <- 100 * abs(sf_ext / exp(mean(log(sf_ext))) - sf_internal / exp(mean(log(sf_internal)))) /
      (sf_internal / exp(mean(log(sf_internal))))
    say(sprintf("[GREENLIST] injecting EXTERNAL (multiBamSummary) factors; max |%%diff| vs internal = %.2f%%", max(pct)))
    sf <- sf_ext
  } else sf <- sf_internal
  base_out <- opt$outdir; dir.create(base_out, showWarnings = FALSE, recursive = TRUE)

  for (u in units) run_unit(u, meta, opt, param, sf, tag, base_out)
  say("\n[DONE] DESeq2/greenlist differential binding complete.")
}

main()