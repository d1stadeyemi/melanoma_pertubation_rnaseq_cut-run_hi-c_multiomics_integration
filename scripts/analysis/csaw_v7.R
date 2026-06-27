#!/usr/bin/env Rscript
# =============================================================================
# csaw.R (v7) — Window-based differential binding for BROAD histone marks
# =============================================================================
# PROJECT: A375 melanoma parental-vs-B10 constricted-migration multi-omics.
#
# CHANGES FROM v6 (summary; full notes inline):
#   A. MULTI-WINDOW-SIZE consolidation via --window 1000,2000,5000. Each size is
#      counted, filtered and TESTED independently, then per-window results across
#      sizes are merged into shared regions with csaw::mergeWindowsList(equiweight)
#      + combineTests (one region-level FDR). Single size == v6 behaviour.
#      Normalization factors are computed ONCE (size-independent) and reused, so
#      the global scale stays coherent across sizes.
#   B. PER-STEP DIAGNOSTIC FIGURES: per-size filtering distributions (background
#      bins + median; window filter stat + threshold; window abundance + cutoff),
#      a normalization-factor audit (per-sample factor + effective lib + within-
#      group CV), and a dispersion comparison across normalizations (+ per-norm BCV).
#   C. GROUP-AWARE FILTERING (--filter_mode group_aware) keeps a window if EITHER
#      contrast group's own average clears the threshold (rescues on/off loci).
#      Default stays 'global' (contrast-blind; v6 behaviour).
#
# NORMALIZATION DIRECTION AUDIT (see compute_nf): edgeR effective lib = N*nf.
#   sf: nf = SF/geomean(SF) -> effective lib = N*SF; larger SF lowers normalized
#       counts = global loss (Jin et al. 2020). spikein: effective lib ∝ spikein
#       reads. Both directions are correct in this script.
# =============================================================================

suppressWarnings(suppressMessages({
  library(optparse)
  library(csaw)
  library(edgeR)
  library(GenomicRanges)
  library(S4Vectors)
}))

# =============================================================================
# CLI
# =============================================================================
parse_cli <- function() {
  option_list <- list(
    make_option("--sheet", type = "character",
                help = "REQUIRED. Sample sheet CSV (SampleID, Factor, Condition, Treatment, Replicate, Batch, bamReads; + optional spikein_reads / Peaks / bamControl)."),
    make_option("--mark", type = "character",
                help = "REQUIRED. Factor/antibody, matched case-insensitively (e.g. H3K27me3 | H3K9me3 | H3K4me1)."),
    make_option("--test", type = "character", default = "contrast",
                help = "Test mode. ALLOWED: 'contrast' | 'interaction'. [default %default]"),
    make_option("--group1", type = "character", default = "B10|none",
                help = "[contrast] Numerator 'Condition|Treatment' (treated). [default %default]"),
    make_option("--group2", type = "character", default = "Parental|none",
                help = "[contrast] Reference 'Condition|Treatment'. Positive logFC = higher in group1. [default %default]"),
    make_option("--cond_ref",  type = "character", default = "Parental",
                help = "[interaction] Reference Condition. [default %default]"),
    make_option("--treat_ref", type = "character", default = "none",
                help = "[interaction] Reference Treatment. [default %default]"),
    make_option("--cond_alt",  type = "character", default = "B10",
                help = "[interaction] Non-reference Condition. [default %default]"),
    make_option("--treat_alt", type = "character", default = "UNC",
                help = "[interaction] Non-reference Treatment. [default %default]"),
    make_option("--outdir", type = "character", default = "csaw_out",
                help = "Output directory. [default %default]"),
    make_option("--norm", type = "character", default = "background",
                help = "Comma-separated normalizations. ALLOWED: libsize | background | greenlist | spikein | sf. 2+ -> 'libsize' auto-added. [default %default]"),
    make_option("--greenlist", type = "character", default = NULL,
                help = "Greenlist normalization regions: 3-col BED in hg19. [default: none]"),
    make_option("--spikein", type = "character", default = NULL,
                help = "Spike-in (Orlando/RRPM): 'sheet' OR CSV(SampleID,spikein_reads). [default: none]"),
    make_option("--normfactors", type = "character", default = NULL,
                help = "Precomputed per-sample SF (ChIPseqSpikeInFree). Table with ID col + factor col (SF/factor/ScalingFactor/norm.factors). [default: none]"),
    make_option("--force_norm", action = "store_true", default = FALSE,
                help = "Apply spikein/sf even when batch is non-separable. [default: FALSE]"),
    make_option("--reference_shift", type = "double", default = NULL,
                help = "External global log2FC (treated vs reference) from IF/western; drawn on the comparison figure. [default: none]"),
    make_option("--orthogonal_mark", type = "character", default = NULL,
                help = "Control mark expected NOT to change (e.g. H3K9me3) for the comparison figure. [default: none]"),
    make_option("--window",  type = "character", default = "2000",
                help = "Window width(s) bp. SINGLE (2000) or COMMA LIST (1000,2000,5000) to test each size and merge. [default %default]"),
    make_option("--spacing", type = "integer", default = 500L,
                help = "Spacing between window starts (bp); also the region-merge tolerance. [default %default]"),
    make_option("--bg_bin",  type = "integer", default = 10000L,
                help = "Background bin width (bp) for TMM (background norm) + the enrichment filter. [default %default]"),
    make_option("--frag_max", type = "integer", default = 700L,
                help = "Max fragment length (bp); match Bowtie2 -X. [default %default]"),
    make_option("--filter_mode", type = "character", default = "global",
                help = "Window filter. ALLOWED: 'global' | 'group_aware'. group_aware also keeps a window if EITHER group's own average passes (rescues on/off loci). [default %default]"),
    make_option("--min_lfc_bg", type = "double", default = log2(3),
                help = "Min log2 enrichment over background/IgG to KEEP a window. [default log2(3)=1.585]"),
    make_option("--fdr", type = "double", default = 0.05,
                help = "FDR cutoff for the significant subset. [default %default]"),
    make_option("--lfc", type = "double", default = 1.0,
                help = "Min |log2FC| for the significant subset. [default %default]"),
    make_option("--n_strata", type = "integer", default = 3L,
                help = "Number of baseline-abundance strata (>=2). [default %default]"),
    make_option("--skip_diagnostics", action = "store_true", default = FALSE,
                help = "Skip per-step diagnostic figures (filter / normfactor / dispersion). [default: FALSE]"),
    make_option("--blacklist", type = "character", default = NULL,
                help = "Optional ENCODE blacklist BED. [default: none]"),
    make_option("--igg", type = "character", default = NULL,
                help = "IgG enrichment filter. 'sheet' OR comma-separated IgG BAM path(s). [default: none]"),
    make_option("--peaks", type = "character", default = NULL,
                help = "RTM-immune peak stratifier. 'sheet' OR comma-separated peak file(s). [default: none]"),
    make_option("--baseline_group", type = "character", default = NULL,
                help = "Independent baseline stratifier 'Condition|Treatment' counted in the same windows. [default: none]"),
    make_option("--genes", type = "character", default = NULL,
                help = "Optional hg19 gene annotation for gene mapping (name col + TSS col OR start/end cols + chrom col). [default: none]")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$sheet)) stop("--sheet is required.")
  if (is.null(opt$mark))  stop("--mark is required (e.g. H3K27me3).")
  if (!opt$test %in% c("contrast", "interaction")) stop("--test must be 'contrast' or 'interaction'.")
  if (!opt$filter_mode %in% c("global", "group_aware")) stop("--filter_mode must be 'global' or 'group_aware'.")
  if (opt$n_strata < 2L) stop("--n_strata must be >= 2.")
  opt
}

# ---- pure-R helpers (unit-tested) -------------------------------------------
parse_window_sizes <- function(s) {
  toks <- trimws(strsplit(as.character(s), ",")[[1]]); toks <- toks[nzchar(toks)]
  if (!length(toks)) stop("--window is empty.")
  v <- suppressWarnings(as.integer(toks))
  if (any(is.na(v))) stop("--window has non-integer token(s): ", paste(toks[is.na(v)], collapse = ", "))
  if (any(v <= 0))   stop("--window sizes must be positive: ", paste(v[v <= 0], collapse = ", "))
  sort(unique(v))
}

combine_keep <- function(keep_global, grp_filter_list, thr) {
  keep <- keep_global
  if (length(grp_filter_list)) for (fg in grp_filter_list) keep <- keep | (fg > thr)
  keep
}

mean_within_group_cv <- function(x, group_key) {
  cvs <- tapply(x, group_key, function(z) if (length(z) > 1 && mean(z) != 0) sd(z) / abs(mean(z)) else NA_real_)
  mean(cvs, na.rm = TRUE)
}

parse_group <- function(spec) {
  parts <- strsplit(spec, "\\|")[[1]]
  if (length(parts) != 2) stop("Group spec must be 'Condition|Treatment', got: ", spec)
  list(cond = tolower(trimws(parts[1])), treat = tolower(trimws(parts[2])))
}

safe_token <- function(s) gsub("[^A-Za-z0-9]+", "", s)

log_lines <- new.env(); log_lines$txt <- character(0)
say <- function(...) { line <- paste0(...); cat(line, "\n"); log_lines$txt <- c(log_lines$txt, line) }

read_bed3 <- function(path) {
  df <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", quote = "")
  GRanges(seqnames = df[[1]], ranges = IRanges(start = as.integer(df[[2]]) + 1L, end = as.integer(df[[3]])))
}

# =============================================================================
# Sample sheet + selection (unchanged from v6)
# =============================================================================
load_sheet <- function(path) {
  ss <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  req <- c("SampleID", "Factor", "Condition", "Treatment", "Replicate", "Batch", "bamReads")
  miss <- setdiff(req, colnames(ss))
  if (length(miss)) stop("Sample sheet missing columns: ", paste(miss, collapse = ", "))
  ss$.cond  <- tolower(trimws(ss$Condition))
  ss$.treat <- tolower(trimws(ss$Treatment))
  ss$.fac   <- tolower(trimws(ss$Factor))
  ss$Batch  <- as.character(ss$Batch)
  ss
}

select_samples <- function(ss, opt) {
  mark <- tolower(trimws(opt$mark))
  ss <- ss[ss$.fac == mark, , drop = FALSE]
  if (!nrow(ss)) stop("No samples for mark '", opt$mark, "'.")
  if (opt$test == "contrast") {
    g1 <- parse_group(opt$group1); g2 <- parse_group(opt$group2)
    key <- function(g) paste(g$cond, g$treat, sep = "|")
    ss$.key <- paste(ss$.cond, ss$.treat, sep = "|")
    sub <- ss[ss$.key %in% c(key(g1), key(g2)), , drop = FALSE]
    if (!nrow(sub)) stop("No samples match group1/group2.")
    sub$Group <- factor(sub$.key, levels = c(key(g2), key(g1)))
    if (any(is.na(sub$Group))) stop("Group level mismatch.")
    sub$.role <- ifelse(sub$.key == key(g1), "group1", "group2")
    attr(sub, "g1key") <- key(g1); attr(sub, "g2key") <- key(g2)
  } else {
    cr <- tolower(opt$cond_ref); ca <- tolower(opt$cond_alt)
    tr <- tolower(opt$treat_ref); ta <- tolower(opt$treat_alt)
    sub <- ss[ss$.cond %in% c(cr, ca) & ss$.treat %in% c(tr, ta), , drop = FALSE]
    if (!nrow(sub)) stop("No samples match the interaction cells.")
    sub$Condition <- factor(sub$.cond, levels = c(cr, ca))
    sub$Treatment <- factor(sub$.treat, levels = c(tr, ta))
    if (any(is.na(sub$Condition)) || any(is.na(sub$Treatment))) stop("Interaction factor level mismatch.")
    attr(sub, "peak_ref") <- sub$.cond == ca & sub$.treat == tr
  }
  if (any(!file.exists(sub$bamReads)))
    stop("Missing BAM(s): ", paste(sub$bamReads[!file.exists(sub$bamReads)], collapse = ", "))
  sub
}

ref_rows <- function(meta, opt) {
  if (opt$test == "contrast") meta$.role == "group2" else attr(meta, "peak_ref")
}

# =============================================================================
# Resolve optional inputs (unchanged from v6)
# =============================================================================
resolve_peaks <- function(opt, meta) {
  if (is.null(opt$peaks)) return(NULL)
  if (identical(opt$peaks, "sheet")) {
    if (!"Peaks" %in% colnames(meta)) stop("--peaks sheet requested but no 'Peaks' column.")
    paths <- unique(meta$Peaks[ref_rows(meta, opt)]); paths <- paths[!is.na(paths) & nzchar(paths)]
    if (!length(paths)) stop("--peaks sheet: no Peaks paths on the reference rows.")
  } else paths <- trimws(strsplit(opt$peaks, ",")[[1]])
  miss <- paths[!file.exists(paths)]; if (length(miss)) stop("Peak file(s) not found: ", paste(miss, collapse = ", "))
  gr <- reduce(Reduce(c, lapply(paths, read_bed3)))
  say(sprintf("[PEAKS] %d file(s) -> %d merged peak regions", length(paths), length(gr)))
  gr
}

resolve_igg <- function(opt, meta) {
  if (is.null(opt$igg)) return(NULL)
  if (identical(opt$igg, "sheet")) {
    if (!"bamControl" %in% colnames(meta)) stop("--igg sheet requested but no 'bamControl' column.")
    bams <- unique(meta$bamControl); bams <- bams[!is.na(bams) & nzchar(bams)]
  } else bams <- trimws(strsplit(opt$igg, ",")[[1]])
  miss <- bams[!file.exists(bams)]; if (length(miss)) stop("IgG BAM(s) not found: ", paste(miss, collapse = ", "))
  bams
}

resolve_spikein <- function(opt, meta) {
  if (is.null(opt$spikein)) return(NULL)
  if (identical(opt$spikein, "sheet")) {
    if (!"spikein_reads" %in% colnames(meta)) stop("--spikein sheet requested but no 'spikein_reads' column.")
    v <- setNames(as.numeric(meta$spikein_reads), meta$SampleID)
  } else {
    sp <- read.csv(opt$spikein, stringsAsFactors = FALSE)
    if (!all(c("SampleID", "spikein_reads") %in% colnames(sp))) stop("--spikein CSV needs SampleID,spikein_reads.")
    v <- setNames(as.numeric(sp$spikein_reads), sp$SampleID)[meta$SampleID]
  }
  if (any(is.na(v)) || any(v <= 0))
    stop("spikein_reads NA/non-positive for: ", paste(meta$SampleID[is.na(v) | v <= 0], collapse = ", "))
  v
}

resolve_normfactors <- function(opt, meta) {
  if (is.null(opt$normfactors)) return(NULL)
  if (!file.exists(opt$normfactors)) stop("--normfactors file not found: ", opt$normfactors)
  first <- readLines(opt$normfactors, n = 1L, warn = FALSE)
  sep   <- if (grepl("\t", first)) "\t" else ","
  tab   <- read.table(opt$normfactors, header = TRUE, sep = sep, quote = "\"",
                      stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE)
  fcol <- intersect(c("SF", "factor", "ScalingFactor", "norm.factors", "normfactor"), colnames(tab))
  if (!length(fcol)) stop("--normfactors: no scaling-factor column found. Columns: ", paste(colnames(tab), collapse = ", "))
  icol <- intersect(c("SampleID", "ID", "Sample", "sample"), colnames(tab))
  if (!length(icol)) stop("--normfactors: no ID column found. Columns: ", paste(colnames(tab), collapse = ", "))
  ids <- as.character(tab[[icol[1]]]); sf <- suppressWarnings(as.numeric(tab[[fcol[1]]]))
  bn  <- basename(meta$bamReads)
  out <- setNames(rep(NA_real_, nrow(meta)), meta$SampleID); qc_used <- character(0)
  for (i in seq_len(nrow(meta))) {
    hit <- which(ids == meta$SampleID[i]); if (!length(hit)) hit <- which(ids == bn[i])
    if (!length(hit)) hit <- which(basename(ids) == bn[i])
    if (length(hit)) { out[i] <- sf[hit[1]]; if ("QC" %in% colnames(tab)) qc_used <- c(qc_used, as.character(tab$QC[hit[1]])) }
  }
  if (length(qc_used) && any(!grepl("pass", tolower(qc_used))))
    say("[NORMFACTORS] WARNING: a matched row has QC != 'pass'. Verify before trusting it.")
  if (any(is.na(out))) stop("--normfactors: could not match a factor for: ", paste(meta$SampleID[is.na(out)], collapse = ", "))
  if (any(out <= 0))   stop("--normfactors: non-positive/NA factor for: ", paste(meta$SampleID[out <= 0], collapse = ", "))
  out
}

resolve_greenlist <- function(opt) {
  if (is.null(opt$greenlist)) return(NULL)
  if (!file.exists(opt$greenlist)) stop("--greenlist file not found: ", opt$greenlist)
  gr <- reduce(read_bed3(opt$greenlist))
  say(sprintf("[GREENLIST] %d normalizer regions loaded (assumed hg19).", length(gr)))
  gr
}

resolve_genes <- function(opt) {
  if (is.null(opt$genes)) return(NULL)
  if (!file.exists(opt$genes)) stop("--genes file not found: ", opt$genes)
  first <- readLines(opt$genes, n = 1L, warn = FALSE)
  sep   <- if (grepl("\t", first)) "\t" else ","
  tab   <- read.table(opt$genes, header = TRUE, sep = sep, quote = "\"",
                      stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE)
  name_col <- intersect(c("gene_name", "gene", "symbol", "Symbol", "name", "Gene", "GeneSymbol"), colnames(tab))
  if (!length(name_col)) stop("--genes: no gene-name column. Columns: ", paste(colnames(tab), collapse = ", "))
  chr_col  <- intersect(c("chr", "chrom", "Chr", "Chrom", "seqnames", "chromosome"), colnames(tab))
  if (!length(chr_col)) stop("--genes: no chromosome column. Columns: ", paste(colnames(tab), collapse = ", "))
  tss_col   <- intersect(c("TSS", "tss", "Tss"), colnames(tab))
  start_col <- intersect(c("start", "Start", "txStart", "bin_start", "geneStart"), colnames(tab))
  end_col   <- intersect(c("end", "End", "txEnd", "bin_end", "geneEnd"), colnames(tab))
  chr  <- as.character(tab[[chr_col[1]]]); gene <- as.character(tab[[name_col[1]]])
  if (length(tss_col)) {
    pos  <- suppressWarnings(as.integer(tab[[tss_col[1]]]))
    keep <- !is.na(pos) & !is.na(chr) & nzchar(gene)
    gr   <- GRanges(seqnames = chr[keep], ranges = IRanges(start = pmax(pos[keep], 1L), width = 1L)); gr$gene <- gene[keep]
    say(sprintf("[GENES] %d genes loaded by TSS from %s.", length(gr), basename(opt$genes)))
  } else if (length(start_col) && length(end_col)) {
    s <- suppressWarnings(as.integer(tab[[start_col[1]]])); e <- suppressWarnings(as.integer(tab[[end_col[1]]]))
    keep <- !is.na(s) & !is.na(e) & !is.na(chr) & nzchar(gene)
    gr   <- GRanges(seqnames = chr[keep], ranges = IRanges(start = pmax(pmin(s[keep], e[keep]), 1L), end = pmax(s[keep], e[keep]))); gr$gene <- gene[keep]
    say(sprintf("[GENES] %d genes loaded by body from %s.", length(gr), basename(opt$genes)))
  } else stop("--genes: need a TSS column, or BOTH start and end columns.")
  if (!length(gr)) stop("--genes: no usable gene rows after parsing ", basename(opt$genes))
  gr
}

# =============================================================================
# Batch assessment (affects spikein/sf only) -- unchanged from v6
# =============================================================================
coef_estimable <- function(design, coef) {
  if (!coef %in% colnames(design)) return(FALSE)
  others <- setdiff(colnames(design), coef)
  r_full <- qr(design)$rank
  r_oth  <- if (length(others)) qr(design[, others, drop = FALSE])$rank else 0L
  (r_full - r_oth) == 1L
}

assess_batch <- function(meta, opt) {
  if (opt$test == "contrast") {
    d_bio <- model.matrix(~ Group, meta); tested <- setdiff(colnames(d_bio), "(Intercept)")
    f_batch <- ~ Batch + Group
  } else {
    d_bio <- model.matrix(~ Condition * Treatment, meta); tested <- grep(":", colnames(d_bio), value = TRUE)
    f_batch <- ~ Batch + Condition * Treatment
  }
  nb <- nlevels(droplevels(factor(meta$Batch)))
  if (nb <= 1L) return(list(status = "same", design = d_bio, tested = tested, allow_spikein = TRUE))
  d_full <- model.matrix(f_batch, meta)
  tested_full <- if (opt$test == "contrast")
    setdiff(colnames(d_full), c("(Intercept)", grep("^Batch", colnames(d_full), value = TRUE)))
  else grep(":", colnames(d_full), value = TRUE)
  full_rank <- qr(d_full)$rank == ncol(d_full)
  tested_ok <- length(tested_full) > 0 && all(vapply(tested_full, function(cf) coef_estimable(d_full, cf), logical(1)))
  if (full_rank && tested_ok)      list(status = "crossed",    design = d_full, tested = tested_full, allow_spikein = TRUE)
  else if (tested_ok)              list(status = "partial",    design = d_bio,  tested = tested,      allow_spikein = FALSE)
  else                             list(status = "confounded", design = d_bio,  tested = tested,      allow_spikein = FALSE)
}

# =============================================================================
# Counting + filtering (MULTI-SIZE)
# =============================================================================
make_param <- function(opt) {
  discard <- GRanges()
  if (!is.null(opt$blacklist)) {
    if (!file.exists(opt$blacklist)) stop("Blacklist not found: ", opt$blacklist)
    discard <- read_bed3(opt$blacklist)
  }
  readParam(minq = 20L, dedup = TRUE, pe = "both", max.frag = opt$frag_max, discard = discard)
}

# Count the test mark for EVERY window size, filter each independently, collect
# per-size diagnostics. Background bins (bg_bin) counted ONCE and shared.
# group_key (per-sample 'Condition|Treatment') enables group-aware mode.
count_and_filter_multi <- function(bams, opt, param, igg_bams, sizes, group_key = NULL) {
  bins   <- windowCounts(bams, bin = TRUE, width = opt$bg_bin, param = param)
  bin_ab <- aveLogCPM(asDGEList(bins))
  per_size <- list(); diags <- list()
  for (w in sizes) {
    win <- windowCounts(bams, width = w, spacing = opt$spacing, param = param)
    if (!is.null(igg_bams)) {
      control     <- regionCounts(igg_bams, regions = rowRanges(win), param = param)
      bin.control <- windowCounts(igg_bams, bin = TRUE, width = opt$bg_bin, param = param)
      scinfo <- scaleControlFilter(bins, bin.control)
      fstat  <- filterWindowsControl(win, control, scale.info = scinfo)
      filt_kind <- "IgG-enrichment (filterWindowsControl)"
    } else {
      fstat  <- filterWindowsGlobal(win, bins)
      filt_kind <- "genomic-background (filterWindowsGlobal)"
    }
    keep_global <- fstat$filter > opt$min_lfc_bg
    grp_filt <- list()
    if (identical(opt$filter_mode, "group_aware") && !is.null(group_key) &&
        length(unique(group_key)) >= 2 && is.null(igg_bams)) {
      for (gg in unique(group_key)) {
        cols <- which(group_key == gg)
        grp_filt[[as.character(gg)]] <- filterWindowsGlobal(win[, cols], bins[, cols])$filter
      }
    } else if (identical(opt$filter_mode, "group_aware") && !is.null(igg_bams)) {
      say("[FILTER] group_aware not applied with --igg (control-based filter); using global keep only.")
    }
    keep <- combine_keep(keep_global, grp_filt, opt$min_lfc_bg)
    n_rescue <- sum(keep & !keep_global)
    diags[[as.character(w)]] <- list(
      size = w, kind = filt_kind, bg_bin = opt$bg_bin,
      bin_ab = bin_ab, bin_median = median(bin_ab, na.rm = TRUE),
      win_ab = as.numeric(fstat$abundances), win_filter = as.numeric(fstat$filter),
      back_level = median(as.numeric(fstat$back.abundances), na.rm = TRUE),
      thr = opt$min_lfc_bg, n_total = length(keep), n_keep = sum(keep),
      n_keep_global = sum(keep_global), n_rescue = n_rescue)
    per_size[[as.character(w)]] <- list(size = w, win = win[keep, ], win_gr = rowRanges(win)[keep],
                                        rescued = (keep & !keep_global)[keep])
    say(sprintf("[FILTER w=%d] %s: kept %d / %d windows (>%.2f log2 over background); group-aware rescued +%d.",
                w, filt_kind, sum(keep), length(keep), opt$min_lfc_bg, n_rescue))
  }
  list(per_size = per_size, bins = bins, diags = diags, sizes = sizes)
}

# =============================================================================
# Per-method normalization factors -- unchanged from v6 (audit in header)
# =============================================================================
compute_nf <- function(method, y0, win, bins, green, spike_v, sf_v) {
  cn <- colnames(y0)
  if (method == "libsize") {
    nf <- rep(1, ncol(y0))
  } else if (method == "background") {
    nf <- asDGEList(normFactors(bins,  se.out = win))$samples$norm.factors
  } else if (method == "greenlist") {
    if (is.null(green)) stop("[greenlist] requested but greenlist counts are missing.")
    nf <- asDGEList(normFactors(green, se.out = win))$samples$norm.factors
  } else if (method == "spikein") {
    sr <- spike_v[cn]; nf <- sr / y0$samples$lib.size; nf <- nf / exp(mean(log(nf)))
  } else if (method == "sf") {
    v <- sf_v[cn]; nf <- v / exp(mean(log(v)))
  } else stop("Unknown normalization method: ", method)
  setNames(as.numeric(nf), cn)
}

dir_labels <- function(opt) {
  if (opt$test == "contrast")
    list(up = "gained", down = "lost", legend_up = "gained", legend_down = "lost")
  else
    list(up = sprintf("stronger loss in %s", opt$cond_ref),
         down = sprintf("stronger loss in %s", opt$cond_alt),
         legend_up = sprintf("more loss in %s", opt$cond_ref),
         legend_down = sprintf("more loss in %s", opt$cond_alt))
}

# =============================================================================
# Fit + test + region-merge across one OR several window sizes
# =============================================================================
fit_test_regions_multi <- function(y_list, ab, opt, win_gr_list, peaks_gr,
                                    indep_logcpm_list, sizes, sig_idx, rescued_list = NULL) {
  ranges.list <- vector("list", length(sizes)); tab.list <- vector("list", length(sizes))
  base.list   <- vector("list", length(sizes)); indep.list <- vector("list", length(sizes))
  disp_summ <- list(); y_sig <- NULL; res_sig <- NULL
  for (i in seq_along(sizes)) {
    y   <- estimateDisp(y_list[[i]], ab$design)
    fit <- glmQLFit(y, ab$design, robust = TRUE)
    res <- glmQLFTest(fit, coef = ab$tested)
    wtab <- res$table
    ref_idx <- if (opt$test == "contrast") which(y$samples$group_role == "group2") else which(y$samples$peak_ref)
    base_lcpm <- aveLogCPM(y$counts[, ref_idx, drop = FALSE],
                           lib.size = y$samples$lib.size[ref_idx] * y$samples$norm.factors[ref_idx])
    ranges.list[[i]] <- win_gr_list[[i]]
    tab.list[[i]]    <- data.frame(logFC = wtab$logFC, logCPM = wtab$logCPM, PValue = wtab$PValue)
    base.list[[i]]   <- base_lcpm
    indep.list[[i]]  <- if (!is.null(indep_logcpm_list)) indep_logcpm_list[[i]] else rep(NA_real_, length(win_gr_list[[i]]))
    disp_summ[[as.character(sizes[i])]] <- data.frame(
      size = sizes[i], n_windows = length(win_gr_list[[i]]),
      common_disp = y$common.dispersion, BCV = sqrt(y$common.dispersion),
      prior_df = if (!is.null(fit$df.prior)) stats::median(fit$df.prior) else NA_real_,
      stringsAsFactors = FALSE)
    if (i == sig_idx) { y_sig <- y; res_sig <- res }
  }
  if (length(sizes) == 1L) {
    merged <- mergeWindows(ranges.list[[1]], tol = opt$spacing, max.width = 10000L)
    ids <- merged$ids; reg <- merged$regions; wts <- NULL
  } else {
    ml  <- mergeWindowsList(ranges.list, tol = opt$spacing, max.width = 10000L)
    ids <- ml$ids; reg <- ml$regions; wts <- ml$weights
  }
  tab.all   <- do.call(rbind, tab.list)
  base.all  <- unlist(base.list, use.names = FALSE)
  indep.all <- unlist(indep.list, use.names = FALSE)
  if (length(ids) != nrow(tab.all))
    stop(sprintf("[MERGE] window/cluster length mismatch (ids=%d, rows=%d).", length(ids), nrow(tab.all)))
  comb <- combineTests(ids, tab.all, weights = wts)
  agg <- function(v) as.numeric(tapply(v, ids, mean))
  df <- data.frame(
    chr = as.character(seqnames(reg)), start = start(reg), end = end(reg),
    logFC = agg(tab.all$logFC), logCPM = agg(tab.all$logCPM),
    PValue = comb$PValue, FDR = comb$FDR, direction = comb$direction,
    n_windows = as.integer(tapply(rep(1L, length(ids)), ids, sum)),
    baseline_logCPM = agg(base.all),
    baseline_indep_logCPM = if (all(is.na(indep.all))) rep(NA_real_, length(reg)) else agg(indep.all),
    stringsAsFactors = FALSE)
  cut_strata <- function(x, k) {
    if (all(is.na(x))) return(factor(rep(NA, length(x))))
    br <- quantile(x, probs = seq(0, 1, length.out = k + 1), na.rm = TRUE)
    lab <- if (k == 3) c("low", "mid", "high") else paste0("s", seq_len(k))
    cut(x, breaks = unique(br), labels = lab[seq_len(length(unique(br)) - 1)], include.lowest = TRUE)
  }
  df$baseline_stratum       <- cut_strata(df$baseline_logCPM, opt$n_strata)
  df$baseline_indep_stratum <- cut_strata(df$baseline_indep_logCPM, opt$n_strata)
  df$in_peak <- if (!is.null(peaks_gr)) overlapsAny(reg, peaks_gr) else NA

  # ---- group-aware filter A/B accounting: which kept windows were RESCUED ----
  if (is.null(rescued_list)) rescued.all <- rep(FALSE, length(ids))
  else rescued.all <- unlist(rescued_list, use.names = FALSE)
  if (length(rescued.all) != length(ids))
    stop(sprintf("[A/B] rescued/cluster length mismatch (rescued=%d, ids=%d).", length(rescued.all), length(ids)))
  if (any(is.na(ids))) stop("[A/B] unexpected NA cluster id.")

  df$n_rescued_windows <- as.integer(tapply(rescued.all, ids, sum))
  df$frac_rescued      <- ifelse(df$n_windows > 0, df$n_rescued_windows / df$n_windows, 0)
  df$region_origin     <- ifelse(df$n_rescued_windows == 0L, "global",
                          ifelse(df$n_rescued_windows == df$n_windows, "rescued_only", "mixed"))

  sig_mask  <- df$FDR < opt$fdr & abs(df$logFC) > opt$lfc
  region_row <- match(ids, sort(unique(ids)))   # window -> df row index (robust to non-contiguous ids)
  reg_sig_w <- sig_mask[region_row]              # is each window's region significant?
  reg_dir_w <- sign(df$logFC)[region_row]        # and in which direction?
  ab_report <- list(
    rescued_windows               = sum(rescued.all),
    rescued_windows_in_sig        = sum(rescued.all & reg_sig_w),
    rescued_in_sig_lost           = sum(rescued.all & reg_sig_w & reg_dir_w < 0),
    rescued_in_sig_gained         = sum(rescued.all & reg_sig_w & reg_dir_w > 0),
    regions_with_rescue           = sum(df$n_rescued_windows > 0L),
    regions_rescued_only          = sum(df$region_origin == "rescued_only"),
    sig_regions_rescued_only      = sum(df$region_origin == "rescued_only" & sig_mask),
    sig_regions_rescued_only_lost = sum(df$region_origin == "rescued_only" & sig_mask & df$logFC < 0),
    sig_regions_mixed             = sum(df$region_origin == "mixed" & sig_mask),
    sig_regions_mixed_lost        = sum(df$region_origin == "mixed" & sig_mask & df$logFC < 0))

  list(table = df, res = res_sig, y = y_sig, disp = do.call(rbind, disp_summ), ab = ab_report)
}

# =============================================================================
# Stratified up/down balance + RTM verdict -- unchanged from v6
# =============================================================================
balance_block <- function(df, stratum_col, opt, label) {
  sig <- df$FDR < opt$fdr & abs(df$logFC) > opt$lfc; s <- df[[stratum_col]]
  say(sprintf("[STRAT:%s]", label)); out <- list()
  for (L in levels(factor(s))) {
    m <- !is.na(s) & s == L; ms <- m & sig
    n_up <- sum(ms & df$logFC > 0); n_dn <- sum(ms & df$logFC < 0); n_sig <- n_up + n_dn
    pct_dn <- if (n_sig) 100 * n_dn / n_sig else NA_real_
    say(sprintf("   %-5s n=%-6d sig=%-5d  up=%-5d down=%-5d  %%down(sig)=%s  medianLFC=%+.3f",
                L, sum(m), n_sig, n_up, n_dn, ifelse(is.na(pct_dn), "NA", sprintf("%.0f", pct_dn)),
                median(df$logFC[m], na.rm = TRUE)))
    out[[L]] <- c(n = sum(m), sig = n_sig, up = n_up, down = n_dn, pct_down = pct_dn)
  }
  out
}

balance_peak <- function(df, opt) {
  if (all(is.na(df$in_peak))) { say("[STRAT:peak]   (no --peaks given)"); return(NULL) }
  sig <- df$FDR < opt$fdr & abs(df$logFC) > opt$lfc; say("[STRAT:peak]"); out <- list()
  for (L in c(TRUE, FALSE)) {
    m <- df$in_peak == L; ms <- m & sig
    n_up <- sum(ms & df$logFC > 0); n_dn <- sum(ms & df$logFC < 0); n_sig <- n_up + n_dn
    pct_dn <- if (n_sig) 100 * n_dn / n_sig else NA_real_
    say(sprintf("   %-8s n=%-6d sig=%-5d  up=%-5d down=%-5d  %%down(sig)=%s  medianLFC=%+.3f",
                ifelse(L, "in_peak", "non_peak"), sum(m), n_sig, n_up, n_dn,
                ifelse(is.na(pct_dn), "NA", sprintf("%.0f", pct_dn)), median(df$logFC[m], na.rm = TRUE)))
    out[[ifelse(L, "peak", "nonpeak")]] <- c(pct_down = pct_dn)
  }
  out
}

rtm_verdict <- function(sd_block, peak_block, indep_block, opt) {
  if (opt$test != "contrast") {
    say("[STRAT verdict] (interaction: flattening verdict suppressed — difference of differences)")
    return(invisible(NULL))
  }
  g <- function(b, hi, lo) if (is.null(b) || is.null(b[[hi]]) || is.null(b[[lo]])) NA_real_ else b[[hi]]["pct_down"] - b[[lo]]["pct_down"]
  grad_sd <- g(sd_block, "high", "low")
  grad_pk <- if (!is.null(peak_block)) (peak_block[["peak"]]["pct_down"] - peak_block[["nonpeak"]]["pct_down"]) else NA_real_
  grad_in <- g(indep_block, "high", "low")
  say("[STRAT verdict]")
  say(sprintf("   samedata loss-gradient (high-low %%down) = %s", ifelse(is.na(grad_sd), "NA", sprintf("%+.0f", grad_sd))))
  if (!is.na(grad_pk)) say(sprintf("   peak     loss-gradient (peak-nonpeak)   = %+.0f  -> %s", grad_pk,
                                   ifelse(!is.na(grad_sd) && grad_pk >= 0.5 * grad_sd, "SURVIVES (flattening real)", "COLLAPSES (largely RTM)")))
  if (!is.na(grad_in)) say(sprintf("   indep    loss-gradient (high-low)       = %+.0f  -> %s", grad_in,
                                   ifelse(!is.na(grad_sd) && grad_in >= 0.5 * grad_sd, "SURVIVES", "weak/RTM")))
  if (is.na(grad_pk) && is.na(grad_in)) say("   (no RTM-resistant stratifier supplied; samedata is RTM-prone)")
}

condition_signal <- function(y, cond_key) {
  lc <- cpm(y, log = TRUE); conds <- unique(cond_key)
  m <- sapply(conds, function(cc) rowMeans(lc[, cond_key == cc, drop = FALSE]))
  if (is.null(dim(m))) m <- matrix(m, ncol = length(conds))
  colnames(m) <- conds; m
}

# =============================================================================
# Figure helpers (unchanged from v6)
# =============================================================================
save_fig <- function(fn, outdir, draw, w = 6.5, h = 5) {
  for (dev in c("pdf", "png")) {
    if (dev == "pdf") pdf(file.path(outdir, paste0(fn, ".pdf")), width = w, height = h)
    else png(file.path(outdir, paste0(fn, ".png")), width = w, height = h, units = "in", res = 200)
    draw(); dev.off()
  }
}

mds_plot <- function(y, tag, norm, outdir, color_by = c("treat", "batch")) {
  color_by <- match.arg(color_by)
  mds <- suppressWarnings(plotMDS(cpm(y, log = TRUE), plot = FALSE))
  xx <- as.numeric(mds$x); yy <- as.numeric(mds$y)
  cond <- as.character(y$samples$cond); treat <- as.character(y$samples$treat); batch <- as.character(y$samples$batch)
  ids <- colnames(y)
  conds <- sort(unique(cond)); pch_map <- setNames(c(21,24,22,23,25)[seq_along(conds)], conds); pchs <- pch_map[cond]
  keyvec <- if (color_by == "treat") treat else batch; keys <- sort(unique(keyvec))
  pal <- c("#1F6F8B","#E08A3C","#2A9D8F","#B23A48","#6D5BA6","#9BA7B0")[seq_along(keys)]
  col_map <- setNames(pal, keys); cols <- col_map[keyvec]
  pct <- function(i) { ve <- mds$var.explained; if (!is.null(ve) && length(ve) >= i && is.finite(ve[i])) sprintf(" (%d%%)", round(100*ve[i])) else "" }
  xr <- range(xx); yr <- range(yy); xpad <- max(diff(xr)*0.30, 0.5); ypad <- max(diff(yr)*0.30, 0.5)
  ttl <- sprintf("%s %s MDS  (shape=Condition, colour=%s)", tag, norm, if (color_by=="batch") "Batch" else "Treatment")
  fn <- paste0(tag, "_", norm, "_MDS", if (color_by=="batch") "_batch" else "")
  save_fig(fn, outdir, function() {
    par(mar = c(5,5,4,9), xpd = TRUE)
    plot(xx, yy, pch = pchs, bg = cols, col = "grey20", cex = 2.2, lwd = 0.9,
         xlim = c(xr[1]-xpad, xr[2]+xpad), ylim = c(yr[1]-ypad, yr[2]+ypad),
         xlab = paste0("Leading logFC dim 1", pct(1)), ylab = paste0("Leading logFC dim 2", pct(2)), main = ttl)
    text(xx, yy, labels = ids, pos = 3, offset = 0.7, cex = 0.6, col = "grey30")
    legend("topright", inset = c(-0.30,0), bty = "n", title = "Condition", pch = unname(pch_map[conds]),
           pt.bg = "grey80", col = "grey20", legend = conds, pt.cex = 1.5)
    legend("bottomright", inset = c(-0.30,0), bty = "n", title = if (color_by=="batch") "Batch" else "Treatment",
           pch = 22, pt.bg = unname(col_map[keys]), col = "grey20", legend = keys, pt.cex = 1.5)
  }, w = 8, h = 5.5)
}

std_plots <- function(y, res, df, tag, norm, outdir, opt) {
  labs <- dir_labels(opt)
  mds_plot(y, tag, norm, outdir, "treat"); mds_plot(y, tag, norm, outdir, "batch")
  save_fig(paste0(tag, "_", norm, "_MA"), outdir, function() {
    plot(df$logCPM, df$logFC, pch = 16, cex = 0.3, col = "#88888855", xlab = "logCPM", ylab = "logFC",
         main = paste(tag, norm, "MA")); sig <- df$FDR < opt$fdr & abs(df$logFC) > opt$lfc
    points(df$logCPM[sig], df$logFC[sig], pch = 16, cex = 0.35, col = "#C0392B"); abline(h = 0, col = "grey40")
  })
  save_fig(paste0(tag, "_", norm, "_volcano"), outdir, function() {
    sig <- df$FDR < opt$fdr & abs(df$logFC) > opt$lfc
    plot(df$logFC, -log10(df$FDR), pch = 16, cex = 0.3, col = "#88888855", xlab = "logFC", ylab = "-log10 FDR",
         main = paste(tag, norm, "volcano"))
    points(df$logFC[sig], -log10(df$FDR[sig]), pch = 16, cex = 0.35, col = ifelse(df$logFC[sig] > 0, "#C0392B", "#1F6F8B"))
    abline(v = c(-opt$lfc, opt$lfc), h = -log10(opt$fdr), col = "grey60", lty = 2)
    n_up <- sum(sig & df$logFC > 0); n_dn <- sum(sig & df$logFC < 0)
    legend("top", horiz = TRUE, bg = "white", box.col = "grey80", cex = 0.95,
           legend = c(sprintf("%s: %d", labs$down, n_dn), sprintf("%s: %d", labs$up, n_up)), text.col = c("#1F6F8B","#C0392B"))
  })
  save_fig(paste0(tag, "_", norm, "_repcor"), outdir, function() {
    lc <- cpm(y, log = TRUE)
    suppressWarnings(pairs(lc[sample(nrow(lc), min(2000, nrow(lc))), , drop = FALSE], pch = ".", main = paste(tag, norm, "rep cor")))
  })
}

strata_balance_fig <- function(df, tag, norm, outdir, opt) {
  labs <- dir_labels(opt); sig <- df$FDR < opt$fdr & abs(df$logFC) > opt$lfc
  has_peak <- !all(is.na(df$in_peak)); has_indep <- !all(is.na(df$baseline_indep_logCPM))
  draw_stratum <- function(strat, title) {
    levs <- levels(factor(df[[strat]])); if (!length(levs)) { plot.new(); title(main = title); return() }
    up <- dn <- numeric(length(levs))
    for (i in seq_along(levs)) { m <- !is.na(df[[strat]]) & df[[strat]] == levs[i] & sig; up[i] <- sum(m & df$logFC > 0); dn[i] <- sum(m & df$logFC < 0) }
    tot <- pmax(up + dn, 1)
    barplot(rbind(100*up/tot, 100*dn/tot), names.arg = levs, col = c("#E08A3C","#2A9D8F"),
            ylab = "% of significant", main = title, ylim = c(0,100))
  }
  draw_peak <- function() {
    if (!has_peak) { plot.new(); title(main = "peak (none)"); return() }
    up <- dn <- numeric(2); levs <- c(TRUE, FALSE)
    for (i in seq_along(levs)) { m <- df$in_peak == levs[i] & sig; up[i] <- sum(m & df$logFC > 0); dn[i] <- sum(m & df$logFC < 0) }
    tot <- pmax(up + dn, 1)
    barplot(rbind(100*up/tot, 100*dn/tot), names.arg = c("in_peak","non_peak"),
            col = c("#E08A3C","#2A9D8F"), ylab = "% of significant", main = "peak (RTM-immune)", ylim = c(0,100))
  }
  n_panels <- 1 + as.integer(has_indep) + as.integer(has_peak)
  draw <- function() {
    par(mfrow = c(1, n_panels), mar = c(4,4,3,1), oma = c(0,0,2.4,0))
    draw_stratum("baseline_stratum", "samedata baseline (RTM-prone)")
    if (has_indep) draw_stratum("baseline_indep_stratum", "independent baseline")
    if (has_peak)  draw_peak()
    mtext(sprintf("%s %s  up/down composition by stratum  (orange=%s, teal=%s)", tag, norm, labs$up, labs$down),
          side = 3, outer = TRUE, cex = 0.85)
  }
  save_fig(paste0(tag, "_", norm, "_strata_balance"), outdir, draw, w = max(5, 3.4*n_panels), h = 4)
}

# =============================================================================
# NEW — per-step DIAGNOSTIC figures
# =============================================================================

# (1) FILTERING diagnostic, one figure per window size.
filter_diagnostics_fig <- function(d, tag, outdir, opt) {
  bgkb <- d$bg_bin / 1000
  cutoff_on_abundance <- d$back_level + d$thr
  draw <- function() {
    par(mfrow = c(1, 3), mar = c(4.5, 4.2, 3.6, 1.2), oma = c(0, 0, 2.2, 0))
    # A: background bins + median
    ha <- hist(d$bin_ab, breaks = 60, plot = FALSE)
    plot(ha, col = "#BFC9CA", border = "white",
         xlab = sprintf("%g kb bin avg logCPM", bgkb), main = sprintf("A. background bins (%g kb)", bgkb))
    abline(v = d$bin_median, col = "#C0392B", lwd = 2, lty = 2)
    text(d$bin_median, max(ha$counts) * 0.95, sprintf("median = %.2f", d$bin_median), pos = 4, col = "#C0392B", cex = 0.8)
    # B: per-window increase over background + threshold (csaw-book Fig 4.3 analog)
    hb <- hist(d$win_filter, breaks = 80, plot = FALSE)
    fail <- hb$mids < d$thr
    plot(hb, col = ifelse(fail, "#E8B7B0", "#7FB3A6"), border = "white",
         xlab = "window log2 increase over background", main = "B. window filter statistic")
    abline(v = d$thr, col = "#C0392B", lwd = 2, lty = 2)
    legend("topright", bty = "n", cex = 0.8,
           legend = c(sprintf("threshold = %.2f", d$thr),
                      sprintf("kept   = %d (%.1f%%)", d$n_keep, 100 * d$n_keep / d$n_total),
                      sprintf("dropped= %d (%.1f%%)", d$n_total - d$n_keep, 100 * (d$n_total - d$n_keep) / d$n_total)),
           text.col = c("#C0392B", "#2A6F5E", "#9B3A2E"))
    # C: window abundance with scaled background + cutoff
    hc <- hist(d$win_ab, breaks = 80, plot = FALSE)
    plot(hc, col = "#CBD5DC", border = "white", xlab = "window avg logCPM", main = "C. window abundance vs cutoff")
    abline(v = d$back_level, col = "#6D5BA6", lwd = 2, lty = 3)
    abline(v = cutoff_on_abundance, col = "#C0392B", lwd = 2, lty = 2)
    legend("topright", bty = "n", cex = 0.78,
           legend = c(sprintf("scaled bg = %.2f", d$back_level), sprintf("keep cutoff = %.2f", cutoff_on_abundance)),
           text.col = c("#6D5BA6", "#C0392B"))
    sub <- if (d$n_rescue > 0) sprintf("  (group-aware rescued +%d on/off windows)", d$n_rescue) else ""
    mtext(sprintf("%s  filtering diagnostics  | window = %d bp | %s%s", tag, d$size, d$kind, sub),
          side = 3, outer = TRUE, cex = 0.9, font = 2)
  }
  save_fig(sprintf("%s_filter_w%d", tag, d$size), outdir, draw, w = 12, h = 4.2)
}

# (2) NORM-FACTOR diagnostic: per-sample factor + effective lib + within-group CV.
normfactor_diagnostics_fig <- function(nf_df, methods, group_key, tag, outdir) {
  pal <- c("#1F6F8B","#E08A3C","#2A9D8F","#B23A48","#6D5BA6")[seq_along(methods)]
  nm  <- setNames(pal, methods)
  cvs <- sapply(methods, function(m) mean_within_group_cv(nf_df[[paste0("efflib_", m)]], group_key))
  draw <- function() {
    par(mfrow = c(1, 2), mar = c(9, 4.5, 3.6, 1), oma = c(0, 0, 2.0, 0))
    n <- nrow(nf_df); xx <- seq_len(n)
    yr <- range(unlist(nf_df[paste0("nf_", methods)]), na.rm = TRUE)
    plot(NA, xlim = c(1, n), ylim = yr, xaxt = "n", xlab = "", ylab = "norm.factor", main = "per-sample norm.factor")
    axis(1, at = xx, labels = nf_df$SampleID, las = 2, cex.axis = 0.62)
    abline(h = 1, col = "grey70", lty = 3)
    for (m in methods) { lines(xx, nf_df[[paste0("nf_", m)]], col = nm[m], lwd = 2); points(xx, nf_df[[paste0("nf_", m)]], pch = 19, col = nm[m]) }
    legend("topleft", bty = "n", legend = methods, col = nm, lwd = 2, cex = 0.8)
    bp <- barplot(cvs, col = nm, ylab = "mean within-group CV of effective lib",
                  main = "replicate scaling spread", las = 2, ylim = c(0, max(cvs, na.rm = TRUE) * 1.25 + 1e-6))
    text(bp, cvs, labels = sprintf("%.3f", cvs), pos = 3, cex = 0.8)
    mtext("higher CV -> replicates scattered -> inflated dispersion -> flatter volcano", side = 1, line = 4.2, cex = 0.72, col = "grey30")
    mtext(sprintf("%s  normalization-factor audit", tag), side = 3, outer = TRUE, cex = 0.95, font = 2)
  }
  save_fig(sprintf("%s_normfactor_diagnostics", tag), outdir, draw, w = 11, h = 4.6)
  out <- nf_df
  for (m in methods) out[[paste0("within_grp_efflib_CV_", m)]] <- cvs[m]
  write.csv(out, file.path(outdir, sprintf("%s_normfactor_diagnostics.csv", tag)), row.names = FALSE)
  invisible(cvs)
}

# (3) DISPERSION comparison across normalizations.
dispersion_compare_fig <- function(disp_by_method, tag, outdir) {
  methods <- names(disp_by_method)
  rows <- do.call(rbind, lapply(methods, function(m) data.frame(method = m, disp_by_method[[m]], stringsAsFactors = FALSE)))
  write.csv(rows, file.path(outdir, sprintf("%s_dispersion_compare.csv", tag)), row.names = FALSE)
  rep_bcv <- sapply(methods, function(m) disp_by_method[[m]]$BCV[1])
  rep_pdf <- sapply(methods, function(m) disp_by_method[[m]]$prior_df[1])
  pal <- c("#1F6F8B","#E08A3C","#2A9D8F","#B23A48","#6D5BA6")[seq_along(methods)]
  draw <- function() {
    par(mfrow = c(1, 2), mar = c(7, 4.5, 3.6, 1), oma = c(0, 0, 2.0, 0))
    bp <- barplot(rep_bcv, col = pal, las = 2, ylab = "common BCV (sqrt dispersion)",
                  main = "dispersion by normalization", ylim = c(0, max(rep_bcv, na.rm = TRUE) * 1.25 + 1e-6))
    text(bp, rep_bcv, labels = sprintf("%.3f", rep_bcv), pos = 3, cex = 0.8)
    bp2 <- barplot(rep_pdf, col = pal, las = 2, ylab = "QL prior df (robust median)",
                   main = "QL prior df by normalization", ylim = c(0, max(rep_pdf, na.rm = TRUE, 1) * 1.25))
    text(bp2, rep_pdf, labels = sprintf("%.1f", rep_pdf), pos = 3, cex = 0.8)
    mtext(sprintf("%s  higher BCV under a normalization = it scatters replicates = it flattens the volcano", tag),
          side = 3, outer = TRUE, cex = 0.82, font = 2)
  }
  save_fig(sprintf("%s_dispersion_compare", tag), outdir, draw, w = 10, h = 4.4)
}

bcv_plot <- function(y, tag, norm, outdir) {
  if (is.null(y) || is.null(y$common.dispersion)) return(invisible(NULL))
  save_fig(sprintf("%s_%s_BCV", tag, norm), outdir, function() {
    suppressWarnings(plotBCV(y, main = sprintf("%s %s  BCV", tag, norm)))
  })
}

# =============================================================================
# Cross-method per-region logFC concordance -- unchanged from v6
# =============================================================================
norm_compare_fig <- function(df_by_method, tag, outdir) {
  methods <- names(df_by_method); if (length(methods) < 2) return(invisible(NULL))
  ref <- if ("libsize" %in% methods) "libsize" else methods[1]
  others <- setdiff(methods, ref)
  draw <- function() {
    par(mfrow = c(1, length(others)), mar = c(4,4,3,1))
    for (m in others) {
      x <- df_by_method[[ref]]$logFC; y <- df_by_method[[m]]$logFC
      n <- min(length(x), length(y)); x <- x[seq_len(n)]; y <- y[seq_len(n)]
      rng <- range(c(x, y), finite = TRUE)
      plot(x, y, pch = 16, cex = 0.25, col = "#33333522", xlim = rng, ylim = rng,
           xlab = sprintf("logFC (%s)", ref), ylab = sprintf("logFC (%s)", m), main = sprintf("%s vs %s", m, ref))
      abline(0, 1, col = "grey40"); abline(h = 0, v = 0, col = "grey85", lty = 3)
      rho <- suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
      legend("topleft", bty = "n", cex = 0.85, legend = sprintf("Spearman rho = %.3f", rho))
    }
  }
  save_fig(paste0(tag, "_norm_compare"), outdir, draw, w = max(5, 3.2*length(others)), h = 3.8)
}

# =============================================================================
# Normalization-comparison figures (Egan 2016 Fig 5 analog) -- unchanged from v6
# =============================================================================
.assemble_marks <- function(test_sig, ortho, mark_test, mark_ortho) {
  marks <- list(list(name = mark_test, sig = test_sig, is_test = TRUE))
  if (!is.null(ortho) && !is.null(ortho$signal))
    marks <- c(marks, list(list(name = mark_ortho, sig = ortho$signal, is_test = FALSE)))
  marks
}

comparison_scatter_fig <- function(test_sig, ortho, methods, g1key, g2key, g1lab, g2lab,
                                    mark_test, mark_ortho, reference_shift, tag, outdir) {
  marks <- .assemble_marks(test_sig, ortho, mark_test, mark_ortho)
  nr <- length(marks); nc <- length(methods)
  draw <- function() {
    par(mfrow = c(nr, nc), mar = c(4, 4, 3.5, 1), oma = c(0, 2.2, if (!is.null(reference_shift)) 2.2 else 0.5, 0))
    for (mk in marks) for (m in methods) {
      S <- mk$sig[[m]]
      if (is.null(S) || !all(c(g1key, g2key) %in% colnames(S))) {
        plot.new(); box(col = "grey85"); title(main = m); text(0.5, 0.5, "n/a", col = "grey60"); next
      }
      x <- S[, g2key]; y <- S[, g1key]; rng <- range(c(x, y), finite = TRUE)
      plot(x, y, pch = 16, cex = 0.25, col = "#33333522", xlim = rng, ylim = rng,
           xlab = sprintf("%s (norm logCPM)", g2lab), ylab = sprintf("%s (norm logCPM)", g1lab), main = m)
      abline(0, 1, col = "grey40", lwd = 1.2)
      if (mk$is_test && !is.null(reference_shift)) abline(reference_shift, 1, col = "#C0392B", lwd = 1.6, lty = 2)
      if (m == "libsize") mtext("(no shift by construction)", side = 3, line = 0.2, cex = 0.6, col = "grey45")
    }
    for (i in seq_len(nr)) mtext(marks[[i]]$name, side = 2, outer = TRUE, line = 0.5, at = 1 - (i - 0.5)/nr, cex = 0.85, font = 2)
    if (!is.null(reference_shift))
      mtext(sprintf("dashed red = IF/western global log2FC = %+.2f (treated vs reference)", reference_shift),
            side = 3, outer = TRUE, line = 0.4, cex = 0.7, col = "#C0392B")
  }
  save_fig(paste0(tag, "_norm_comparison_scatter"), outdir, draw, w = max(5, 3.2*nc), h = max(4, 3.2*nr))
}

comparison_box_fig <- function(test_sig, ortho, methods, g1key, g2key, g1lab, g2lab,
                               mark_test, mark_ortho, reference_shift, tag, outdir) {
  marks <- .assemble_marks(test_sig, ortho, mark_test, mark_ortho)
  nr <- length(marks); nc <- length(methods)
  draw <- function() {
    par(mfrow = c(nr, nc), mar = c(6, 4, 3.5, 1), oma = c(0, 2.2, 0.5, 0))
    for (mk in marks) for (m in methods) {
      S <- mk$sig[[m]]
      if (is.null(S) || !all(c(g1key, g2key) %in% colnames(S))) {
        plot.new(); box(col = "grey85"); title(main = m); text(0.5, 0.5, "n/a", col = "grey60"); next
      }
      ref <- S[, g2key]; trt <- S[, g1key]
      boxplot(list(ref, trt), names = c(g2lab, g1lab), col = c("#BFC9CA", "#E08A3C"),
              ylab = "norm logCPM", main = m, las = 2, outline = FALSE)
      if (mk$is_test && !is.null(reference_shift)) {
        expct <- median(ref, na.rm = TRUE) + reference_shift
        segments(1.6, expct, 2.4, expct, col = "#C0392B", lwd = 2, lty = 2)
        text(2, expct, "IF", pos = 3, col = "#C0392B", cex = 0.7)
      }
      if (m == "libsize") mtext("(no shift by construction)", side = 3, line = 0.2, cex = 0.6, col = "grey45")
    }
    for (i in seq_len(nr)) mtext(marks[[i]]$name, side = 2, outer = TRUE, line = 0.5, at = 1 - (i - 0.5)/nr, cex = 0.85, font = 2)
  }
  save_fig(paste0(tag, "_norm_comparison_box"), outdir, draw, w = max(5, 3.0*nc), h = max(4, 3.2*nr))
}

# =============================================================================
# Orthogonal control-mark panel data -- uses ONE representative (signal) size
# =============================================================================
mark_panel_data <- function(mark_name, g1key, g2key, methods, ss, opt, param, green_gr, sig_size) {
  mk  <- tolower(trimws(mark_name))
  sub <- ss[ss$.fac == mk, , drop = FALSE]
  if (!nrow(sub)) { say(sprintf("[ORTHO] no samples for control mark '%s' — control row skipped.", mark_name)); return(NULL) }
  sub$.key <- paste(sub$.cond, sub$.treat, sep = "|")
  sub <- sub[sub$.key %in% c(g1key, g2key), , drop = FALSE]
  n1 <- sum(sub$.key == g1key); n2 <- sum(sub$.key == g2key)
  if (n1 < 1 || n2 < 1) {
    say(sprintf("[ORTHO] control mark '%s' lacks both conditions (%s:%d, %s:%d) — control row skipped.", mark_name, g1key, n1, g2key, n2))
    return(NULL)
  }
  if (any(!file.exists(sub$bamReads))) { say(sprintf("[ORTHO] control mark '%s' BAM(s) missing — control row skipped.", mark_name)); return(NULL) }
  say(sprintf("[ORTHO] control mark '%s': %d samples — building control panels at %d bp.", mark_name, nrow(sub), sig_size))
  bams <- sub$bamReads
  cf   <- count_and_filter_multi(bams, opt, param, igg_bams = NULL, sizes = sig_size, group_key = NULL)
  win  <- cf$per_size[[1]]$win; bins <- cf$bins
  y0   <- asDGEList(win); colnames(y0) <- sub$SampleID
  cond_key <- sub$.key
  green <- if ("greenlist" %in% methods && !is.null(green_gr)) regionCounts(bams, regions = green_gr, param = param) else NULL
  sv <- NULL
  if ("spikein" %in% methods && "spikein_reads" %in% colnames(sub)) {
    cand <- suppressWarnings(as.numeric(sub$spikein_reads)); names(cand) <- sub$SampleID
    if (!any(is.na(cand)) && all(cand > 0)) sv <- cand
    else say("[ORTHO] control-mark spikein_reads missing/invalid — spikein panel skipped for control mark.")
  }
  out <- setNames(vector("list", length(methods)), methods)
  for (m in methods) {
    if (m == "sf")                          { out[[m]] <- NULL; next }
    if (m == "spikein" && is.null(sv))      { out[[m]] <- NULL; next }
    if (m == "greenlist" && is.null(green)) { out[[m]] <- NULL; next }
    nf <- tryCatch(compute_nf(m, y0, win, bins, green, sv, NULL),
                   error = function(e) { say(sprintf("[ORTHO] method '%s' failed for control mark: %s", m, conditionMessage(e))); NULL })
    if (is.null(nf)) { out[[m]] <- NULL; next }
    y <- y0; y$samples$norm.factors <- nf
    out[[m]] <- condition_signal(y, cond_key)
  }
  list(signal = out, g1 = g1key, g2 = g2key)
}

# =============================================================================
# Write per-method region tables; gene mapping -- unchanged from v6
# =============================================================================
write_outputs <- function(df, tag, norm, outdir, opt) {
  labs <- dir_labels(opt)
  df$call <- ifelse(df$FDR < opt$fdr & abs(df$logFC) > opt$lfc, ifelse(df$logFC > 0, labs$up, labs$down), "ns")
  all_path <- file.path(outdir, sprintf("%s_%s_all_regions.csv", tag, norm))
  sig_path <- file.path(outdir, sprintf("%s_%s_significant.csv", tag, norm))
  write.csv(df, all_path, row.names = FALSE)
  sig <- df[df$call != "ns", , drop = FALSE]; sig <- sig[order(sig$FDR), , drop = FALSE]
  write.csv(sig, sig_path, row.names = FALSE)
  n_up <- sum(sig$logFC > 0); n_dn <- sum(sig$logFC < 0)
  say(sprintf("[%s] regions: %d total, %d significant (%s=%d, %s=%d). -> %s",
              norm, nrow(df), nrow(sig), labs$up, n_up, labs$down, n_dn, basename(sig_path)))
  c(total = nrow(df), sig = nrow(sig), up = n_up, down = n_dn)
}

map_genes_to_regions <- function(df, genes_gr, tag, norm, outdir, opt) {
  if (is.null(genes_gr)) return(invisible(NULL))
  labs  <- dir_labels(opt)
  gldir <- file.path(outdir, "genelists"); dir.create(gldir, showWarnings = FALSE, recursive = TRUE)
  reg     <- GRanges(df$chr, IRanges(df$start, df$end))
  ov_all  <- findOverlaps(genes_gr, reg)
  bg_genes <- sort(unique(genes_gr$gene[queryHits(ov_all)]))
  sig    <- df$FDR < opt$fdr & abs(df$logFC) > opt$lfc
  up_idx <- which(sig & df$logFC > 0); dn_idx <- which(sig & df$logFC < 0)
  genes_for <- function(idx) {
    if (!length(idx)) return(character(0))
    sort(unique(genes_gr$gene[queryHits(findOverlaps(genes_gr, reg[idx]))]))
  }
  up_genes <- genes_for(up_idx); dn_genes <- genes_for(dn_idx)
  up_tok <- safe_token(labs$up); dn_tok <- safe_token(labs$down)
  writeLines(up_genes, file.path(gldir, sprintf("%s_%s_%s.txt", tag, norm, up_tok)))
  writeLines(dn_genes, file.path(gldir, sprintf("%s_%s_%s.txt", tag, norm, dn_tok)))
  writeLines(bg_genes, file.path(gldir, "background_universe.txt"))
  if (length(ov_all)) {
    asg <- data.frame(chr = as.character(seqnames(reg))[subjectHits(ov_all)],
                      start = start(reg)[subjectHits(ov_all)], end = end(reg)[subjectHits(ov_all)],
                      gene = genes_gr$gene[queryHits(ov_all)], stringsAsFactors = FALSE)
    write.csv(asg, file.path(gldir, sprintf("%s_%s_gene_assignment.csv", tag, norm)), row.names = FALSE)
  }
  say(sprintf("[GENES:%s] %s=%d, %s=%d, background=%d -> %s",
              norm, labs$up, length(up_genes), labs$down, length(dn_genes), length(bg_genes), file.path(basename(outdir), "genelists")))
  if (!length(bg_genes)) say("[GENES] WARNING: 0 genes overlapped any region — check chromosome naming matches your BAMs.")
  invisible(list(up = up_genes, down = dn_genes, background = bg_genes))
}

# =============================================================================
# Group-aware filter A/B reporter: what did rescuing on/off windows buy?
#   Writes <tag>_<norm>_filter_AB.csv and logs a plain-language summary. Region
#   provenance is also in the region tables via the new columns:
#     n_rescued_windows, frac_rescued, region_origin (global | mixed | rescued_only).
#   'rescued_only' regions literally do not exist under --filter_mode global, so a
#   SIGNIFICANT rescued_only region is a pure additional discovery from group-aware
#   filtering (the headline number for the B10-UNC loss).
# =============================================================================
report_filter_ab <- function(ab, tag, norm, outdir, opt) {
  if (is.null(ab)) return(invisible(NULL))
  labs <- dir_labels(opt)
  abdf <- data.frame(metric = names(ab), value = as.integer(unlist(ab)), stringsAsFactors = FALSE)
  write.csv(abdf, file.path(outdir, sprintf("%s_%s_filter_AB.csv", tag, norm)), row.names = FALSE)
  if (identical(opt$filter_mode, "group_aware") && ab$rescued_windows > 0) {
    say(sprintf("[A/B:%s] group-aware rescued %d windows -> %d fell in significant regions (%s=%d, %s=%d).",
                norm, ab$rescued_windows, ab$rescued_windows_in_sig, labs$down, ab$rescued_in_sig_lost, labs$up, ab$rescued_in_sig_gained))
    say(sprintf("[A/B:%s] regions made ENTIRELY of rescued windows: %d total, %d significant (%d %s) — these do NOT exist under --filter_mode global.",
                norm, ab$regions_rescued_only, ab$sig_regions_rescued_only, ab$sig_regions_rescued_only_lost, labs$down))
    say(sprintf("[A/B:%s] partly-rescued (mixed) significant regions: %d (%d %s) — exist under global but with different boundaries/power.",
                norm, ab$sig_regions_mixed, ab$sig_regions_mixed_lost, labs$down))
  } else if (identical(opt$filter_mode, "group_aware")) {
    say(sprintf("[A/B:%s] group-aware filter rescued 0 windows (nothing on/off to recover at this threshold).", norm))
  } else {
    say(sprintf("[A/B:%s] filter_mode=global (no rescue). Re-run with --filter_mode group_aware to quantify on/off discoveries.", norm))
  }
  invisible(ab)
}

# =============================================================================
# Per-normalization driver (MULTI-SIZE aware)
# =============================================================================
run_norm <- function(method, per_size, bins, green, spike_v, sf_v, peaks_gr, indep_list,
                     ab, tag, outdir, opt, meta, cond_key, genes_gr, sizes, sig_idx,
                     skip_diag) {
  say(sprintf("\n========== NORMALIZATION: %s ==========", method))
  norm_dir <- file.path(outdir, method); dir.create(norm_dir, showWarnings = FALSE, recursive = TRUE)

  # nf computed ONCE from the size-independent inputs, on the signal-size window object.
  win_sig <- per_size[[sig_idx]]$win
  y0_sig  <- asDGEList(win_sig); colnames(y0_sig) <- meta$SampleID
  nf      <- compute_nf(method, y0_sig, win_sig, bins, green, spike_v, sf_v)

  # build a per-size DGEList list, all carrying the SAME nf + sample metadata.
  y_list <- vector("list", length(sizes)); win_gr_list <- vector("list", length(sizes))
  for (i in seq_along(sizes)) {
    yi <- asDGEList(per_size[[i]]$win); colnames(yi) <- meta$SampleID
    yi$samples$norm.factors <- nf
    yi$samples$cond  <- meta$.cond; yi$samples$treat <- meta$.treat; yi$samples$batch <- meta$Batch
    if (opt$test == "contrast") { yi$samples$group_role <- meta$.role; yi$samples$peak_ref <- NA }
    else                        { yi$samples$peak_ref <- attr(meta, "peak_ref"); yi$samples$group_role <- NA }
    y_list[[i]] <- yi; win_gr_list[[i]] <- per_size[[i]]$win_gr
  }
  rescued_list <- lapply(seq_along(sizes), function(i) per_size[[i]]$rescued)

  ftr <- fit_test_regions_multi(y_list, ab, opt, win_gr_list, peaks_gr, indep_list, sizes, sig_idx, rescued_list)
  df  <- ftr$table; yfit <- ftr$y

  counts <- write_outputs(df, tag, method, norm_dir, opt)
  std_plots(yfit, ftr$res, df, tag, method, norm_dir, opt)
  if (!skip_diag) bcv_plot(yfit, tag, method, norm_dir)

  sd_block <- balance_block(df, "baseline_stratum", opt, "samedata(RTM-prone)")
  in_block <- if (!all(is.na(df$baseline_indep_logCPM))) balance_block(df, "baseline_indep_stratum", opt, "indep") else NULL
  pk_block <- balance_peak(df, opt)
  strata_balance_fig(df, tag, method, norm_dir, opt)
  rtm_verdict(sd_block, pk_block, in_block, opt)

  map_genes_to_regions(df, genes_gr, tag, method, norm_dir, opt)
  report_filter_ab(ftr$ab, tag, method, norm_dir, opt)

  sig <- if (opt$test == "contrast") condition_signal(yfit, cond_key) else NULL
  list(df = df, counts = counts, nf = nf, eff_libsize = y0_sig$samples$lib.size * nf,
       signal = sig, disp = ftr$disp, ab = ftr$ab)
}

# =============================================================================
# MAIN
# =============================================================================
main <- function() {
  opt <- parse_cli()
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

  methods_req <- tolower(trimws(strsplit(opt$norm, ",")[[1]])); methods_req <- unique(methods_req[nzchar(methods_req)])
  allowed <- c("libsize", "background", "greenlist", "spikein", "sf")
  bad <- setdiff(methods_req, allowed)
  if (length(bad)) stop("Unknown --norm token(s): ", paste(bad, collapse = ", "), ". Allowed: ", paste(allowed, collapse = ", "))
  if ("greenlist" %in% methods_req && is.null(opt$greenlist))   stop("--norm includes 'greenlist' but --greenlist not given.")
  if ("spikein"   %in% methods_req && is.null(opt$spikein))     stop("--norm includes 'spikein' but --spikein not given.")
  if ("sf"        %in% methods_req && is.null(opt$normfactors)) stop("--norm includes 'sf' but --normfactors not given.")

  sizes   <- parse_window_sizes(opt$window)
  sig_idx <- which.min(abs(sizes - 2000L))            # representative size for QC/signal (nearest 2 kb)

  ss   <- load_sheet(opt$sheet)
  meta <- select_samples(ss, opt)

  tag <- if (opt$test == "contrast")
    sprintf("%s_%s_vs_%s", opt$mark, safe_token(opt$group1), safe_token(opt$group2))
  else sprintf("%s_interaction", opt$mark)
  say(sprintf("=== csaw v7 | mark=%s | test=%s | windows=%s | filter=%s | tag=%s ===",
              opt$mark, opt$test, paste(sizes, collapse = ","), opt$filter_mode, tag))
  say(sprintf("[SAMPLES] %d libraries: %s", nrow(meta), paste(meta$SampleID, collapse = ", ")))

  ab <- assess_batch(meta, opt)
  say(sprintf("[BATCH] status=%s; tested coef=%s; design cols: %s",
              ab$status, paste(ab$tested, collapse = "+"), paste(colnames(ab$design), collapse = ", ")))

  methods <- methods_req
  if (!ab$allow_spikein && !opt$force_norm) {
    drop <- intersect(c("spikein", "sf"), methods)
    if (length(drop)) say(sprintf("[BATCH] dropping %s (batch %s, not separable; override with --force_norm).", paste(drop, collapse = ","), ab$status))
    methods <- setdiff(methods, drop)
  } else if (!ab$allow_spikein && opt$force_norm) {
    say("[BATCH] --force_norm: applying spikein/sf despite non-separable batch.")
  }
  if (length(methods) > 1 && !"libsize" %in% methods) {
    say("[NORM] adding 'libsize' as the assumption-free reference for the comparison figures.")
    methods <- c("libsize", methods)
  }
  methods <- unique(methods)
  say(sprintf("[NORM] methods to run: %s", paste(methods, collapse = ", ")))

  param    <- make_param(opt)
  igg_bams <- resolve_igg(opt, meta)
  peaks_gr <- resolve_peaks(opt, meta)
  green_gr <- if ("greenlist" %in% methods) resolve_greenlist(opt) else NULL
  spike_v  <- if ("spikein"   %in% methods) resolve_spikein(opt, meta) else NULL
  sf_v     <- if ("sf"        %in% methods) resolve_normfactors(opt, meta) else NULL
  genes_gr <- resolve_genes(opt)

  group_key <- paste(meta$.cond, meta$.treat, sep = "|")
  cf <- count_and_filter_multi(meta$bamReads, opt, param, igg_bams, sizes, group_key = group_key)
  per_size <- cf$per_size; bins <- cf$bins
  green <- if ("greenlist" %in% methods) regionCounts(meta$bamReads, regions = green_gr, param = param) else NULL

  # ---- filtering diagnostics (one per size) ----
  if (!opt$skip_diagnostics) {
    diagdir <- file.path(opt$outdir, "diagnostics"); dir.create(diagdir, showWarnings = FALSE, recursive = TRUE)
    for (w in sizes) filter_diagnostics_fig(cf$diags[[as.character(w)]], tag, diagdir, opt)
    say(sprintf("[DIAG] filtering diagnostics written for window size(s): %s", paste(sizes, collapse = ",")))
  }

  # ---- optional independent baseline group (per size) ----
  indep_list <- NULL
  if (!is.null(opt$baseline_group)) {
    bg <- parse_group(opt$baseline_group); bkey <- paste(bg$cond, bg$treat, sep = "|")
    brow <- ss[ss$.fac == tolower(opt$mark) & paste(ss$.cond, ss$.treat, sep = "|") == bkey, , drop = FALSE]
    if (!nrow(brow)) stop("--baseline_group '", opt$baseline_group, "' matched no samples.")
    if (any(!file.exists(brow$bamReads))) stop("baseline_group BAM(s) missing.")
    say(sprintf("[BASELINE] independent stratifier '%s': %d samples.", opt$baseline_group, nrow(brow)))
    indep_list <- vector("list", length(sizes))
    for (i in seq_along(sizes)) {
      bc <- regionCounts(brow$bamReads, regions = per_size[[i]]$win_gr, param = param)
      indep_list[[i]] <- aveLogCPM(assay(bc), lib.size = bc$totals)
    }
  }

  cond_key <- group_key

  df_by_method <- list(); counts_by_method <- list(); nf_by_method <- list()
  eff_by_method <- list(); sig_by_method <- list(); disp_by_method <- list(); ab_by_method <- list()
  for (m in methods) {
    r <- run_norm(m, per_size, bins, green, spike_v, sf_v, peaks_gr, indep_list,
                  ab, tag, opt$outdir, opt, meta, cond_key, genes_gr, sizes, sig_idx, opt$skip_diagnostics)
    df_by_method[[m]] <- r$df; counts_by_method[[m]] <- r$counts
    nf_by_method[[m]] <- r$nf; eff_by_method[[m]] <- r$eff_libsize
    sig_by_method[[m]] <- r$signal; disp_by_method[[m]] <- r$disp; ab_by_method[[m]] <- r$ab
  }

  norm_compare_fig(df_by_method, tag, opt$outdir)

  # ---- consolidated per-sample factor table ----
  nf_df <- data.frame(SampleID = meta$SampleID, Condition = as.character(meta$Condition),
                      Treatment = as.character(meta$Treatment), Batch = meta$Batch,
                      lib.size = asDGEList(per_size[[sig_idx]]$win)$samples$lib.size, stringsAsFactors = FALSE)
  for (m in methods) { nf_df[[paste0("nf_", m)]] <- nf_by_method[[m]]; nf_df[[paste0("efflib_", m)]] <- eff_by_method[[m]] }
  nf_path <- file.path(opt$outdir, sprintf("%s_normfactors_used.csv", tag))
  write.csv(nf_df, nf_path, row.names = FALSE)
  say(sprintf("[NORMFACTORS] per-sample factor + effective library size per method -> %s", basename(nf_path)))

  # ---- normalization-factor + dispersion diagnostics (the flat-volcano audit) ----
  if (!opt$skip_diagnostics && length(methods) >= 1) {
    normfactor_diagnostics_fig(nf_df, methods, group_key, tag, opt$outdir)
    dispersion_compare_fig(disp_by_method, tag, opt$outdir)
    say("[DIAG] normalization-factor + dispersion diagnostics written.")
  }

  # ---- normalization-comparison figure (contrast only) ----
  if (opt$test == "contrast") {
    g1key <- attr(meta, "g1key"); g2key <- attr(meta, "g2key")
    ortho <- if (!is.null(opt$orthogonal_mark))
      mark_panel_data(opt$orthogonal_mark, g1key, g2key, methods, ss, opt, param, green_gr, sizes[sig_idx]) else NULL
    comparison_scatter_fig(sig_by_method, ortho, methods, g1key, g2key, opt$group1, opt$group2,
                           opt$mark, opt$orthogonal_mark, opt$reference_shift, tag, opt$outdir)
    comparison_box_fig(sig_by_method, ortho, methods, g1key, g2key, opt$group1, opt$group2,
                       opt$mark, opt$orthogonal_mark, opt$reference_shift, tag, opt$outdir)
    say("[COMPARE] normalization-comparison scatter + box written.")
    if (is.null(opt$reference_shift))
      say("[COMPARE] NOTE: no --reference_shift; panels show each method's ASSUMED global shift, not a validated magnitude.")
  } else {
    say("[COMPARE] normalization-comparison figure is contrast-only; skipped for interaction.")
  }

  # ---- summary report ----
  summ <- file.path(opt$outdir, sprintf("%s_summary.txt", tag))
  hdr <- c(strrep("=", 72), sprintf("csaw v7 summary | %s", tag), strrep("=", 72), "",
           sprintf("Windows: %s bp | filter_mode: %s | signal size: %d bp", paste(sizes, collapse = ","), opt$filter_mode, sizes[sig_idx]),
           sprintf("Significant region counts per normalization (FDR<%.2g, |logFC|>%.2g):", opt$fdr, opt$lfc))
  rows <- vapply(methods, function(m) { cc <- counts_by_method[[m]]
    sprintf("  %-12s total=%-7d sig=%-6d up=%-6d down=%-6d", m, cc["total"], cc["sig"], cc["up"], cc["down"]) }, character(1))
  ab_block <- character(0)
  if (identical(opt$filter_mode, "group_aware")) {
    ab_block <- c("", "Group-aware filter A/B (rescued = kept here but DROPPED under --filter_mode global):",
      vapply(methods, function(m) { a <- ab_by_method[[m]]
        sprintf("  %-12s rescued_win=%-6d in_sig=%-5d  sig_regions_rescued_only=%-4d  sig_regions_mixed=%-4d",
                m, a$rescued_windows, a$rescued_windows_in_sig, a$sig_regions_rescued_only, a$sig_regions_mixed) }, character(1)))
  }
  writeLines(c(hdr, rows, ab_block, "", "--- full run log ---", log_lines$txt), summ)
  say(sprintf("\n[DONE] summary -> %s", summ))
}

main()