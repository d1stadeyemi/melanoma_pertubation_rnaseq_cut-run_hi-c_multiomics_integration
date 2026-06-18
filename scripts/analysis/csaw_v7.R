#!/usr/bin/env Rscript
# =============================================================================
# csaw.R  (v7) — Window-based differential binding for BROAD histone marks
#                with MULTI-WINDOW (sliding-window ladder) consolidation
# =============================================================================
#
# PROJECT: A375 melanoma parental-vs-B10 constricted-migration multi-omics.
#
# WHAT THIS SCRIPT DOES --------------------------------------------------------
#   Broad marks (H3K27me3, H3K9me3, H3K4me1) form diffuse domains. We TILE the
#   genome into fixed windows, test each for differential coverage with edgeR's
#   NB quasi-likelihood model, merge adjacent windows into regions, and combine
#   their p-values (Simes -> region-level BH FDR): the FDR-correct region workflow
#   of Lun & Smyth (csaw). Punctate H3K27Ac is handled by diffbind.R.
#
# WHAT CHANGED FROM v6 ---------------------------------------------------------
#   MULTI-WINDOW LADDER. --window now accepts a COMMA-SEPARATED LIST of widths.
#   Each width is independently counted, enrichment-filtered, and DB-tested; the
#   per-width results are then CONSOLIDATED into one region set with a single
#   combined FDR via csaw::mergeResultsList(..., equiweight=TRUE). equiweight is
#   mandatory: a smaller width yields more windows, and without equal weighting
#   the smallest width would dominate every region's Simes combination. The payoff
#   is robustness — the loci-level result no longer hinges on one arbitrary width.
#     --window 2000              (default; single width, behaves like v6)
#     --window 1000,2000,5000    (ladder; consolidated to one FDR)
#
#   Two consequences of consolidation, both handled automatically:
#     * Normalization is computed ONCE per method and shared across all widths.
#       Every normalizer (background TMM on bins, greenlist TMM on greenlist
#       regions, spikein, sf) derives from data INDEPENDENT of the test-window
#       width, so the per-sample factors are identical across widths — computing
#       once is both correct and cheaper. Only the enrichment filter and the
#       dispersion/QL fit are genuinely per-width.
#     * Baseline / strata / peak overlap / per-condition signal / MDS now operate
#       on the CONSOLIDATED regions (reads are recounted directly in them). That
#       is strictly better than window-level estimates and matches the units that
#       are actually reported.
#
#   --max_width default raised to 20 kb (broad H3K27me3 domains run large).
#
#   Everything else is v6: true library-size vs background split, --norm selector,
#   greenlist/spikein/sf, test-aware labels, contrast-only strata verdict, the
#   consolidated per-sample factor table, and the Egan-style normalization-
#   comparison figure with --reference_shift (IF anchor) and --orthogonal_mark.
#
# TWO TEST MODES (--test) -----------------------------------------------------
#   contrast    (default) --group1 NUMERATOR --group2 REFERENCE ("Cond|Treat")
#                 Migration:     --group1 "B10|none"     --group2 "Parental|none"
#                 UNC in B10:    --group1 "B10|UNC"      --group2 "B10|none"
#                 UNC in parent: --group1 "Parental|UNC" --group2 "Parental|none"
#   interaction  Factorial Condition*Treatment; tests the interaction coefficient
#                 (does UNC's effect differ between B10 and Parental?). All four
#                 cells; CASE-INSENSITIVE matching.
#
# SAMPLE SHEET ----------------------------------------------------------------
#   Required : SampleID, Factor, Condition, Treatment, Replicate, Batch, bamReads
#   Optional : spikein_reads (--spikein sheet), Peaks (--peaks sheet),
#              bamControl (--igg sheet)
#   Condition Parental|B10 ; Treatment none|UNC (case-insensitive). Condition/
#   Treatment columns are the source of truth for grouping; keep names consistent.
#
# EXAMPLE RUNS ----------------------------------------------------------------
#   # UNC reversal in B10, window ladder, greenlist, IF anchor, control mark:
#   Rscript csaw.R --sheet ss.csv --mark H3K27me3 \
#       --group1 "B10|UNC" --group2 "B10|none" \
#       --window 1000,2000,5000 \
#       --norm libsize,background,greenlist,spikein,sf \
#       --greenlist data/greenlist_hg19.bed \
#       --spikein sheet --normfactors data/spikeinfree_SF.txt \
#       --peaks sheet --baseline_group "Parental|none" \
#       --reference_shift -2.32 --orthogonal_mark H3K9me3 \
#       --outdir results/diff/H3K27me3_B10unc
#
#   # Migration contrast, single default width:
#   Rscript csaw.R --sheet ss.csv --mark H3K27me3 \
#       --group1 "B10|none" --group2 "Parental|none" --peaks sheet \
#       --outdir results/diff/H3K27me3_migration
#
#   # Interaction (no comparison figure; flattening verdict suppressed):
#   Rscript csaw.R --sheet ss.csv --mark H3K27me3 --test interaction \
#       --window 1000,2000,5000 --peaks sheet \
#       --outdir results/diff/H3K27me3_interaction
#
# OUTPUTS (into --outdir) -----------------------------------------------------
#   <tag>_<norm>_all_regions.csv / _significant.csv      (consolidated regions)
#   <tag>_<norm>_{MDS,MDS_batch,MA,volcano,repcor,strata_balance}.{pdf,png}
#   <tag>_<norm>_norm_compare.{pdf,png}
#   <tag>_norm_comparison_scatter.{pdf,png} / _box.{pdf,png}   (contrast only)
#   <tag>_normfactors_used.csv
#   <tag>_summary.txt
# =============================================================================

suppressWarnings(suppressMessages({
  library(optparse)
  library(csaw)            # windowCounts / regionCounts / normFactors / mergeResultsList
  library(edgeR)           # NB-QL GLM, dispersion, DGEList, cpm, aveLogCPM
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
                help = "REQUIRED. Factor/antibody, matched case-insensitively to the Factor column (e.g. H3K27me3 | H3K9me3 | H3K4me1)."),

    make_option("--test", type = "character", default = "contrast",
                help = "Test mode. ALLOWED: 'contrast' | 'interaction'. [default %default]"),

    make_option("--group1", type = "character", default = "B10|none",
                help = "[contrast] Numerator 'Condition|Treatment' (treated). [default %default]"),
    make_option("--group2", type = "character", default = "Parental|none",
                help = "[contrast] Reference 'Condition|Treatment'. Positive logFC = higher in group1. [default %default]"),
    make_option("--cond_ref",  type = "character", default = "Parental", help = "[interaction] Reference Condition (Parental|B10). [default %default]"),
    make_option("--treat_ref", type = "character", default = "none",     help = "[interaction] Reference Treatment (none|UNC). [default %default]"),
    make_option("--cond_alt",  type = "character", default = "B10",      help = "[interaction] Non-reference Condition. [default %default]"),
    make_option("--treat_alt", type = "character", default = "UNC",      help = "[interaction] Non-reference Treatment. [default %default]"),

    make_option("--outdir", type = "character", default = "csaw_out", help = "Output directory. [default %default]"),

    make_option("--norm", type = "character", default = "libsize,background",
                help = "Comma-separated normalizations. ALLOWED: libsize | background | greenlist | spikein | sf. greenlist needs --greenlist; spikein needs --spikein; sf needs --normfactors. [default %default]"),
    make_option("--greenlist", type = "character", default = NULL,
                help = "Greenlist normalization regions: 3-col BED (hg19) of de Mello CUT&RUN greenlist; TMM is computed on counts in these regions. [default: none]"),
    make_option("--spikein", type = "character", default = NULL,
                help = "Spike-in (Orlando/RRPM): norm.factor = spikein_reads/lib.size. 'sheet' OR a CSV with SampleID,spikein_reads. [default: none]"),
    make_option("--normfactors", type = "character", default = NULL,
                help = "Precomputed per-sample factors (e.g. ChIPseqSpikeInFree). Tab/comma table; ID column (SampleID/ID/Sample) + factor column (SF/factor/...). [default: none]"),
    make_option("--force_norm", action = "store_true", default = FALSE,
                help = "Apply spikein/sf even when batch is non-separable. Only with an external magnitude anchor. [default: FALSE]"),

    make_option("--reference_shift", type = "double", default = NULL,
                help = "Externally measured global log2FC of the mark, TREATED (group1) vs REFERENCE (group2), from IF/western. Negative = loss (5-fold = log2(1/5) = -2.32). Drawn as a reference band on the comparison figure. [default: none]"),
    make_option("--orthogonal_mark", type = "character", default = NULL,
                help = "Control mark expected NOT to change (e.g. H3K9me3). If its contrast-condition samples exist in the sheet, the comparison figure gains a control row; else accepted and skipped with a notice. [default: none]"),

    # ---- WINDOW LADDER ----
    make_option("--window", type = "character", default = "2000",
                help = "Window width(s) in bp. ONE value or a COMMA-SEPARATED LADDER (e.g. 1000,2000,5000); each is counted/filtered/tested then consolidated to one combined FDR via mergeResultsList(equiweight=TRUE). [default %default]"),
    make_option("--spacing", type = "integer", default = 500L,
                help = "Spacing between window starts (bp). Capped to the smallest window if larger (else the smallest scan leaves gaps). [default %default]"),
    make_option("--max_width", type = "integer", default = 20000L,
                help = "Max width (bp) of a consolidated region (caps very broad merged domains). [default %default]"),
    make_option("--bg_bin",  type = "integer", default = 10000L,
                help = "Background bin width (bp) for TMM (background norm) + enrichment filter. [default %default]"),
    make_option("--frag_max", type = "integer", default = 700L,
                help = "Max fragment length (bp); match Bowtie2 -X. [default %default]"),

    make_option("--min_lfc_bg", type = "double", default = log2(2), help = "Min log2 enrichment over background/IgG to KEEP a window. [default log2(2)=1.0]"),
    make_option("--fdr", type = "double", default = 0.05, help = "FDR cutoff for the significant subset. [default %default]"),
    make_option("--lfc", type = "double", default = 1.0,  help = "Min |log2FC| for the significant subset. [default %default]"),
    make_option("--n_strata", type = "integer", default = 3L, help = "Number of baseline-abundance strata (>=2). 3 -> low/mid/high. [default %default]"),

    make_option("--blacklist", type = "character", default = NULL, help = "Optional ENCODE blacklist BED (hg19). [default: none]"),
    make_option("--igg", type = "character", default = NULL, help = "IgG enrichment filter. 'sheet' OR comma-separated IgG BAM path(s). [default: none]"),
    make_option("--peaks", type = "character", default = NULL, help = "RTM-immune peak stratifier. 'sheet' OR comma-separated BED/narrowPeak/broadPeak/SEACR path(s). Use a BROAD call for H3K27me3. [default: none]"),
    make_option("--baseline_group", type = "character", default = NULL, help = "Independent baseline stratifier: a 'Condition|Treatment' group counted in the same regions, used ONLY to rank baseline. [default: none]")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$sheet)) stop("--sheet is required.")
  if (is.null(opt$mark))  stop("--mark is required (e.g. H3K27me3).")
  if (!opt$test %in% c("contrast", "interaction")) stop("--test must be 'contrast' or 'interaction'.")
  if (opt$n_strata < 2L) stop("--n_strata must be >= 2.")
  opt
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
# Sample sheet + selection
# =============================================================================
load_sheet <- function(path) {
  ss <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  req <- c("SampleID", "Factor", "Condition", "Treatment", "Replicate", "Batch", "bamReads")
  miss <- setdiff(req, colnames(ss)); if (length(miss)) stop("Sample sheet missing columns: ", paste(miss, collapse = ", "))
  ss$.cond  <- tolower(trimws(ss$Condition)); ss$.treat <- tolower(trimws(ss$Treatment))
  ss$.fac   <- tolower(trimws(ss$Factor));    ss$Batch  <- as.character(ss$Batch)
  ss
}

select_samples <- function(ss, opt) {
  mark <- tolower(trimws(opt$mark)); ss <- ss[ss$.fac == mark, , drop = FALSE]
  if (!nrow(ss)) stop("No samples for mark '", opt$mark, "'.")
  if (opt$test == "contrast") {
    g1 <- parse_group(opt$group1); g2 <- parse_group(opt$group2)
    key <- function(g) paste(g$cond, g$treat, sep = "|")
    ss$.key <- paste(ss$.cond, ss$.treat, sep = "|")
    sub <- ss[ss$.key %in% c(key(g1), key(g2)), , drop = FALSE]
    if (!nrow(sub)) stop("No samples match group1/group2.")
    sub$Group <- factor(sub$.key, levels = c(key(g2), key(g1)))      # reference = group2
    if (any(is.na(sub$Group))) stop("Group level mismatch.")
    sub$.role <- ifelse(sub$.key == key(g1), "group1", "group2")
    attr(sub, "g1key") <- key(g1); attr(sub, "g2key") <- key(g2)
  } else {
    cr <- tolower(opt$cond_ref); ca <- tolower(opt$cond_alt); tr <- tolower(opt$treat_ref); ta <- tolower(opt$treat_alt)
    sub <- ss[ss$.cond %in% c(cr, ca) & ss$.treat %in% c(tr, ta), , drop = FALSE]
    if (!nrow(sub)) stop("No samples match the interaction cells.")
    sub$Condition <- factor(sub$.cond, levels = c(cr, ca)); sub$Treatment <- factor(sub$.treat, levels = c(tr, ta))
    if (any(is.na(sub$Condition)) || any(is.na(sub$Treatment))) stop("Interaction factor level mismatch.")
    attr(sub, "peak_ref") <- sub$.cond == ca & sub$.treat == tr
  }
  if (any(!file.exists(sub$bamReads))) stop("Missing BAM(s): ", paste(sub$bamReads[!file.exists(sub$bamReads)], collapse = ", "))
  sub
}

ref_rows <- function(meta, opt) if (opt$test == "contrast") meta$.role == "group2" else attr(meta, "peak_ref")

# =============================================================================
# Resolve optional inputs
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
  say(sprintf("[PEAKS] %d file(s) -> %d merged peak regions", length(paths), length(gr))); gr
}
resolve_igg <- function(opt, meta) {
  if (is.null(opt$igg)) return(NULL)
  if (identical(opt$igg, "sheet")) {
    if (!"bamControl" %in% colnames(meta)) stop("--igg sheet requested but no 'bamControl' column.")
    bams <- unique(meta$bamControl); bams <- bams[!is.na(bams) & nzchar(bams)]
  } else bams <- trimws(strsplit(opt$igg, ",")[[1]])
  miss <- bams[!file.exists(bams)]; if (length(miss)) stop("IgG BAM(s) not found: ", paste(miss, collapse = ", ")); bams
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
  if (any(is.na(v)) || any(v <= 0)) stop("spikein_reads NA/non-positive for: ", paste(meta$SampleID[is.na(v) | v <= 0], collapse = ", ")); v
}
resolve_normfactors <- function(opt, meta) {
  if (is.null(opt$normfactors)) return(NULL)
  if (!file.exists(opt$normfactors)) stop("--normfactors file not found: ", opt$normfactors)
  first <- readLines(opt$normfactors, n = 1L, warn = FALSE); sep <- if (grepl("\t", first)) "\t" else ","
  tab <- read.table(opt$normfactors, header = TRUE, sep = sep, quote = "\"", stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE)
  fcol <- intersect(c("SF", "factor", "ScalingFactor", "norm.factors", "normfactor"), colnames(tab))
  if (!length(fcol)) stop("--normfactors: no scaling-factor column. Columns: ", paste(colnames(tab), collapse = ", "))
  icol <- intersect(c("SampleID", "ID", "Sample", "sample"), colnames(tab))
  if (!length(icol)) stop("--normfactors: no ID column. Columns: ", paste(colnames(tab), collapse = ", "))
  ids <- as.character(tab[[icol[1]]]); sf <- suppressWarnings(as.numeric(tab[[fcol[1]]])); bn <- basename(meta$bamReads)
  out <- setNames(rep(NA_real_, nrow(meta)), meta$SampleID); qc_used <- character(0)
  for (i in seq_len(nrow(meta))) {
    hit <- which(ids == meta$SampleID[i]); if (!length(hit)) hit <- which(ids == bn[i]); if (!length(hit)) hit <- which(basename(ids) == bn[i])
    if (length(hit)) { out[i] <- sf[hit[1]]; if ("QC" %in% colnames(tab)) qc_used <- c(qc_used, as.character(tab$QC[hit[1]])) }
  }
  if (length(qc_used) && any(!grepl("pass", tolower(qc_used)))) say("[NORMFACTORS] WARNING: a matched row has QC != 'pass'. Verify before trusting it.")
  if (any(is.na(out))) stop("--normfactors: could not match a factor for: ", paste(meta$SampleID[is.na(out)], collapse = ", "))
  if (any(out <= 0))   stop("--normfactors: non-positive factor for: ", paste(meta$SampleID[out <= 0], collapse = ", ")); out
}
resolve_greenlist <- function(opt) {
  if (is.null(opt$greenlist)) return(NULL)
  if (!file.exists(opt$greenlist)) stop("--greenlist file not found: ", opt$greenlist)
  gr <- reduce(read_bed3(opt$greenlist))
  say(sprintf("[GREENLIST] %d normalizer regions loaded (hg19).", length(gr))); gr
}

# =============================================================================
# Batch assessment (gate affects spikein/sf only)
# =============================================================================
coef_estimable <- function(design, coef) {
  if (!coef %in% colnames(design)) return(FALSE)
  others <- setdiff(colnames(design), coef)
  (qr(design)$rank - (if (length(others)) qr(design[, others, drop = FALSE])$rank else 0L)) == 1L
}
assess_batch <- function(meta, opt) {
  if (opt$test == "contrast") {
    d_bio <- model.matrix(~ Group, meta); tested <- setdiff(colnames(d_bio), "(Intercept)"); f_batch <- ~ Batch + Group
  } else {
    d_bio <- model.matrix(~ Condition * Treatment, meta); tested <- grep(":", colnames(d_bio), value = TRUE); f_batch <- ~ Batch + Condition * Treatment
  }
  nb <- nlevels(droplevels(factor(meta$Batch)))
  if (nb <= 1L) return(list(status = "same", design = d_bio, tested = tested, allow_spikein = TRUE))
  d_full <- model.matrix(f_batch, meta)
  tested_full <- if (opt$test == "contrast") setdiff(colnames(d_full), c("(Intercept)", grep("^Batch", colnames(d_full), value = TRUE))) else grep(":", colnames(d_full), value = TRUE)
  full_rank <- qr(d_full)$rank == ncol(d_full)
  tested_ok <- length(tested_full) > 0 && all(vapply(tested_full, function(cf) coef_estimable(d_full, cf), logical(1)))
  if (full_rank && tested_ok) list(status = "crossed", design = d_full, tested = tested_full, allow_spikein = TRUE)
  else if (tested_ok)         list(status = "partial", design = d_bio,  tested = tested,      allow_spikein = FALSE)
  else                        list(status = "confounded", design = d_bio, tested = tested,    allow_spikein = FALSE)
}

make_param <- function(opt) {
  discard <- GRanges()
  if (!is.null(opt$blacklist)) { if (!file.exists(opt$blacklist)) stop("Blacklist not found: ", opt$blacklist); discard <- read_bed3(opt$blacklist) }
  readParam(minq = 20L, dedup = TRUE, pe = "both", max.frag = opt$frag_max, discard = discard)
}

# =============================================================================
# Per-WIDTH count + enrichment filter (shared bins / bin_control passed in)
# =============================================================================
count_filter_one <- function(bams, W, opt, param, bins, igg_bams, bin_control) {
  win <- windowCounts(bams, width = W, spacing = opt$spacing, param = param)
  if (!is.null(igg_bams)) {                               # IgG-aware enrichment filter
    control <- regionCounts(igg_bams, regions = rowRanges(win), param = param)
    scinfo  <- scaleControlFilter(bins, bin_control)
    fstat   <- filterWindowsControl(win, control, scale.info = scinfo)
  } else fstat <- filterWindowsGlobal(win, bins)          # genomic-background enrichment filter
  win[fstat$filter > opt$min_lfc_bg, ]
}

# =============================================================================
# Per-method normalization factors — computed ONCE, shared across all widths.
#   Every input (bins, greenlist counts, spikein reads, sf) is independent of the
#   test-window width, so the factors are width-invariant. Effective library size
#   in edgeR = lib.size * norm.factor.
# =============================================================================
compute_nf <- function(method, y_ref, bins, green, spike_v, sf_v) {
  cn <- colnames(y_ref)
  if (method == "libsize")         nf <- rep(1, length(cn))                                          # depth-only
  else if (method == "background") nf <- asDGEList(normFactors(bins,  se.out = bins))$samples$norm.factors  # TMM on bins
  else if (method == "greenlist")  { if (is.null(green)) stop("[greenlist] counts missing."); nf <- asDGEList(normFactors(green, se.out = green))$samples$norm.factors }
  else if (method == "spikein")    { sr <- spike_v[cn]; nf <- sr / y_ref$samples$lib.size; nf <- nf / exp(mean(log(nf))) }  # Orlando/RRPM
  else if (method == "sf")         { v <- sf_v[cn]; nf <- v / exp(mean(log(v))) }                    # injected SF
  else stop("Unknown normalization method: ", method)
  setNames(as.numeric(nf), cn)
}

# Low..high strata from baseline ranks.
cut_strata <- function(x, k) {
  if (all(is.na(x))) return(factor(rep(NA, length(x))))
  br <- quantile(x, probs = seq(0, 1, length.out = k + 1), na.rm = TRUE)
  lab <- if (k == 3) c("low", "mid", "high") else paste0("s", seq_len(k))
  cut(x, breaks = unique(br), labels = lab[seq_len(length(unique(br)) - 1)], include.lowest = TRUE)
}

# =============================================================================
# Consolidate per-width DB results into one region set + one combined FDR.
#   Single width  -> mergeWindows + combineTests (identical to v6).
#   Multiple widths -> mergeResultsList(equiweight=TRUE): clusters windows across
#   widths into regions and combines per-window p-values with equal weight per
#   width (so the most numerous, smallest width cannot dominate).
# =============================================================================
consolidate <- function(per_size, opt) {
  wins <- lapply(per_size, `[[`, "win"); tabs <- lapply(per_size, `[[`, "tab")
  if (length(wins) == 1) {
    merged <- mergeWindows(rowRanges(wins[[1]]), tol = opt$spacing, max.width = opt$max_width)
    list(regions = merged$regions, combined = combineTests(merged$ids, tabs[[1]]))
  } else {
    cons <- mergeResultsList(wins, tab.list = tabs, equiweight = TRUE,
                             merge.args = list(tol = opt$spacing, max.width = opt$max_width))
    list(regions = cons$regions, combined = cons$combined)
  }
}

# =============================================================================
# Direction labels — TEST-AWARE
#   contrast:    +logFC = gained mark, -logFC = lost mark.
#   interaction: coefficient = (UNC effect in cond_alt) - (UNC effect in cond_ref),
#                a difference of differences -> +/- mean which condition LOSES MORE.
# =============================================================================
dir_labels <- function(opt) {
  if (opt$test == "contrast")
    list(up = "gained", down = "lost", legend_up = "gained", legend_down = "lost")
  else
    list(up = sprintf("stronger loss in %s", opt$cond_ref), down = sprintf("stronger loss in %s", opt$cond_alt),
         legend_up = sprintf("more loss in %s", opt$cond_ref), legend_down = sprintf("more loss in %s", opt$cond_alt))
}

# =============================================================================
# Stratified up/down balance + RTM verdict (verdict meaningful only for contrast)
# =============================================================================
balance_block <- function(df, stratum_col, opt, label) {
  sig <- df$FDR < opt$fdr & abs(df$logFC) > opt$lfc; s <- df[[stratum_col]]
  say(sprintf("[STRAT:%s]", label)); out <- list()
  for (L in levels(factor(s))) {
    m <- !is.na(s) & s == L; ms <- m & sig
    n_up <- sum(ms & df$logFC > 0); n_dn <- sum(ms & df$logFC < 0); n_sig <- n_up + n_dn
    pct_dn <- if (n_sig) 100 * n_dn / n_sig else NA_real_
    say(sprintf("   %-5s n=%-6d sig=%-5d  up=%-5d down=%-5d  %%down(sig)=%s  medianLFC=%+.3f",
                L, sum(m), n_sig, n_up, n_dn, ifelse(is.na(pct_dn), "NA", sprintf("%.0f", pct_dn)), median(df$logFC[m], na.rm = TRUE)))
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
                ifelse(L, "in_peak", "non_peak"), sum(m), n_sig, n_up, n_dn, ifelse(is.na(pct_dn), "NA", sprintf("%.0f", pct_dn)), median(df$logFC[m], na.rm = TRUE)))
    out[[ifelse(L, "peak", "nonpeak")]] <- c(pct_down = pct_dn)
  }
  out
}
rtm_verdict <- function(sd_block, peak_block, indep_block, opt) {
  if (opt$test != "contrast") { say("[STRAT verdict] (interaction: flattening verdict suppressed — the coefficient is a difference of differences, not a treatment effect)"); return(invisible(NULL)) }
  g <- function(b, hi, lo) if (is.null(b) || is.null(b[[hi]]) || is.null(b[[lo]])) NA_real_ else b[[hi]]["pct_down"] - b[[lo]]["pct_down"]
  grad_sd <- g(sd_block, "high", "low")
  grad_pk <- if (!is.null(peak_block)) (peak_block[["peak"]]["pct_down"] - peak_block[["nonpeak"]]["pct_down"]) else NA_real_
  grad_in <- g(indep_block, "high", "low")
  say("[STRAT verdict]")
  say(sprintf("   samedata loss-gradient (high-low %%down) = %s", ifelse(is.na(grad_sd), "NA", sprintf("%+.0f", grad_sd))))
  if (!is.na(grad_pk)) say(sprintf("   peak     loss-gradient (peak-nonpeak)   = %+.0f  -> %s", grad_pk, ifelse(!is.na(grad_sd) && grad_pk >= 0.5 * grad_sd, "SURVIVES (flattening real)", "COLLAPSES (largely RTM)")))
  if (!is.na(grad_in)) say(sprintf("   indep    loss-gradient (high-low)       = %+.0f  -> %s", grad_in, ifelse(!is.na(grad_sd) && grad_in >= 0.5 * grad_sd, "SURVIVES", "weak/RTM")))
  if (is.na(grad_pk) && is.na(grad_in)) say("   (no RTM-resistant stratifier supplied; samedata is RTM-prone — interpret with care)")
}

# Per-condition normalized signal on the consolidated regions (for comparison fig).
condition_signal <- function(y, cond_key) {
  lc <- cpm(y, log = TRUE); conds <- unique(cond_key)
  m <- sapply(conds, function(cc) rowMeans(lc[, cond_key == cc, drop = FALSE]))
  if (is.null(dim(m))) m <- matrix(m, ncol = length(conds))
  colnames(m) <- conds; m
}

# =============================================================================
# Standard per-method figures (operate on the region-level DGEList)
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
  mds <- suppressWarnings(plotMDS(cpm(y, log = TRUE), plot = FALSE)); xx <- as.numeric(mds$x); yy <- as.numeric(mds$y)
  cond <- as.character(y$samples$cond); treat <- as.character(y$samples$treat); batch <- as.character(y$samples$batch); ids <- colnames(y)
  conds <- sort(unique(cond)); pch_map <- setNames(c(21,24,22,23,25)[seq_along(conds)], conds); pchs <- pch_map[cond]
  keyvec <- if (color_by == "treat") treat else batch; keys <- sort(unique(keyvec))
  pal <- c("#1F6F8B","#E08A3C","#2A9D8F","#B23A48","#6D5BA6","#9BA7B0")[seq_along(keys)]; col_map <- setNames(pal, keys); cols <- col_map[keyvec]
  pct <- function(i) { ve <- mds$var.explained; if (!is.null(ve) && length(ve) >= i && is.finite(ve[i])) sprintf(" (%d%%)", round(100*ve[i])) else "" }
  xr <- range(xx); yr <- range(yy); xpad <- max(diff(xr)*0.30, 0.5); ypad <- max(diff(yr)*0.30, 0.5)
  ttl <- sprintf("%s %s MDS  (shape=Condition, colour=%s)", tag, norm, if (color_by=="batch") "Batch" else "Treatment")
  fn <- paste0(tag, "_", norm, "_MDS", if (color_by=="batch") "_batch" else "")
  save_fig(fn, outdir, function() {
    par(mar = c(5,5,4,9), xpd = TRUE)
    plot(xx, yy, pch = pchs, bg = cols, col = "grey20", cex = 2.2, lwd = 0.9, xlim = c(xr[1]-xpad, xr[2]+xpad), ylim = c(yr[1]-ypad, yr[2]+ypad),
         xlab = paste0("Leading logFC dim 1", pct(1)), ylab = paste0("Leading logFC dim 2", pct(2)), main = ttl)
    text(xx, yy, labels = ids, pos = 3, offset = 0.7, cex = 0.6, col = "grey30")
    legend("topright", inset = c(-0.30,0), bty = "n", title = "Condition", pch = unname(pch_map[conds]), pt.bg = "grey80", col = "grey20", legend = conds, pt.cex = 1.5)
    legend("bottomright", inset = c(-0.30,0), bty = "n", title = if (color_by=="batch") "Batch" else "Treatment", pch = 22, pt.bg = unname(col_map[keys]), col = "grey20", legend = keys, pt.cex = 1.5)
  }, w = 8, h = 5.5)
}
std_plots <- function(y, res, df, tag, norm, outdir, opt) {
  labs <- dir_labels(opt)
  mds_plot(y, tag, norm, outdir, "treat"); mds_plot(y, tag, norm, outdir, "batch")
  save_fig(paste0(tag, "_", norm, "_MA"), outdir, function() {
    plot(df$logCPM, df$logFC, pch = 16, cex = 0.3, col = "#88888855", xlab = "logCPM", ylab = "logFC", main = paste(tag, norm, "MA"))
    sig <- df$FDR < opt$fdr & abs(df$logFC) > opt$lfc; points(df$logCPM[sig], df$logFC[sig], pch = 16, cex = 0.35, col = "#C0392B"); abline(h = 0, col = "grey40")
  })
  save_fig(paste0(tag, "_", norm, "_volcano"), outdir, function() {
    sig <- df$FDR < opt$fdr & abs(df$logFC) > opt$lfc
    plot(df$logFC, -log10(df$FDR), pch = 16, cex = 0.3, col = "#88888855", xlab = "logFC", ylab = "-log10 FDR", main = paste(tag, norm, "volcano"))
    points(df$logFC[sig], -log10(df$FDR[sig]), pch = 16, cex = 0.35, col = ifelse(df$logFC[sig] > 0, "#C0392B", "#1F6F8B"))
    abline(v = c(-opt$lfc, opt$lfc), h = -log10(opt$fdr), col = "grey60", lty = 2)
    n_up <- sum(sig & df$logFC > 0); n_dn <- sum(sig & df$logFC < 0)
    legend("top", horiz = TRUE, bg = "white", box.col = "grey80", cex = 0.95, legend = c(sprintf("%s: %d", labs$down, n_dn), sprintf("%s: %d", labs$up, n_up)), text.col = c("#1F6F8B","#C0392B"))
  })
  save_fig(paste0(tag, "_", norm, "_repcor"), outdir, function() {
    lc <- cpm(y, log = TRUE); suppressWarnings(pairs(lc[sample(nrow(lc), min(2000, nrow(lc))), , drop = FALSE], pch = ".", main = paste(tag, norm, "rep cor")))
  })
}
strata_balance_fig <- function(df, tag, norm, outdir, opt) {
  labs <- dir_labels(opt); sig <- df$FDR < opt$fdr & abs(df$logFC) > opt$lfc
  has_peak <- !all(is.na(df$in_peak)); has_indep <- !all(is.na(df$baseline_indep_logCPM))
  draw_stratum <- function(strat, title) {
    levs <- levels(factor(df[[strat]])); if (!length(levs)) { plot.new(); title(main = title); return() }
    up <- dn <- numeric(length(levs))
    for (i in seq_along(levs)) { m <- !is.na(df[[strat]]) & df[[strat]] == levs[i] & sig; up[i] <- sum(m & df$logFC > 0); dn[i] <- sum(m & df$logFC < 0) }
    tot <- pmax(up + dn, 1); barplot(rbind(100*up/tot, 100*dn/tot), names.arg = levs, col = c("#E08A3C","#2A9D8F"), ylab = "% of significant", main = title, ylim = c(0,100))
  }
  draw_peak <- function() {
    if (!has_peak) { plot.new(); title(main = "peak (none)"); return() }
    up <- dn <- numeric(2); levs <- c(TRUE, FALSE)
    for (i in seq_along(levs)) { m <- df$in_peak == levs[i] & sig; up[i] <- sum(m & df$logFC > 0); dn[i] <- sum(m & df$logFC < 0) }
    tot <- pmax(up + dn, 1); barplot(rbind(100*up/tot, 100*dn/tot), names.arg = c("in_peak","non_peak"), col = c("#E08A3C","#2A9D8F"), ylab = "% of significant", main = "peak (RTM-immune)", ylim = c(0,100))
  }
  n_panels <- 1 + as.integer(has_indep) + as.integer(has_peak)
  draw <- function() {
    par(mfrow = c(1, n_panels), mar = c(4,4,3,1), oma = c(0,0,2.4,0))
    draw_stratum("baseline_stratum", "samedata baseline (RTM-prone)")
    if (has_indep) draw_stratum("baseline_indep_stratum", "independent baseline")
    if (has_peak)  draw_peak()
    mtext(sprintf("%s %s  up/down composition by stratum  (orange=%s, teal=%s)", tag, norm, labs$up, labs$down), side = 3, outer = TRUE, cex = 0.85)
  }
  save_fig(paste0(tag, "_", norm, "_strata_balance"), outdir, draw, w = max(5, 3.4*n_panels), h = 4)
}
norm_compare_fig <- function(df_by_method, tag, outdir) {
  methods <- names(df_by_method); if (length(methods) < 2) return(invisible(NULL))
  ref <- if ("libsize" %in% methods) "libsize" else methods[1]; others <- setdiff(methods, ref)
  draw <- function() {
    par(mfrow = c(1, length(others)), mar = c(4,4,3,1))
    for (m in others) {
      x <- df_by_method[[ref]]$logFC; y <- df_by_method[[m]]$logFC; rng <- range(c(x, y), finite = TRUE)
      plot(x, y, pch = 16, cex = 0.25, col = "#33333522", xlim = rng, ylim = rng, xlab = sprintf("logFC (%s)", ref), ylab = sprintf("logFC (%s)", m), main = sprintf("%s vs %s", m, ref))
      abline(0, 1, col = "grey40"); abline(h = 0, v = 0, col = "grey85", lty = 3)
      rho <- suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs")); legend("topleft", bty = "n", cex = 0.85, legend = sprintf("Spearman rho = %.3f", rho))
    }
  }
  save_fig(paste0(tag, "_norm_compare"), outdir, draw, w = max(5, 3.2*length(others)), h = 3.8)
}

# =============================================================================
# Normalization-comparison figures (Egan Fig 5 analog; contrast only)
#   The shift each method shows is its scaling ASSUMPTION, not a validated
#   magnitude. The dashed red IF line / band is the independent ground truth; an
#   orthogonal control-mark row tests whether a method manufactures shifts.
# =============================================================================
.assemble_marks <- function(test_sig, ortho, mark_test, mark_ortho) {
  marks <- list(list(name = mark_test, sig = test_sig, is_test = TRUE))
  if (!is.null(ortho) && !is.null(ortho$signal)) marks <- c(marks, list(list(name = mark_ortho, sig = ortho$signal, is_test = FALSE)))
  marks
}
comparison_scatter_fig <- function(test_sig, ortho, methods, g1key, g2key, g1lab, g2lab, mark_test, mark_ortho, reference_shift, tag, outdir) {
  marks <- .assemble_marks(test_sig, ortho, mark_test, mark_ortho); nr <- length(marks); nc <- length(methods)
  draw <- function() {
    par(mfrow = c(nr, nc), mar = c(4, 4, 3.5, 1), oma = c(0, 2.2, if (!is.null(reference_shift)) 2.2 else 0.5, 0))
    for (mk in marks) for (m in methods) {
      S <- mk$sig[[m]]
      if (is.null(S) || !all(c(g1key, g2key) %in% colnames(S))) { plot.new(); box(col = "grey85"); title(main = m); text(0.5, 0.5, "n/a", col = "grey60"); next }
      x <- S[, g2key]; y <- S[, g1key]; rng <- range(c(x, y), finite = TRUE)
      plot(x, y, pch = 16, cex = 0.25, col = "#33333522", xlim = rng, ylim = rng, xlab = sprintf("%s (norm logCPM)", g2lab), ylab = sprintf("%s (norm logCPM)", g1lab), main = m)
      abline(0, 1, col = "grey40", lwd = 1.2)
      if (mk$is_test && !is.null(reference_shift)) abline(reference_shift, 1, col = "#C0392B", lwd = 1.6, lty = 2)
      if (m == "libsize") mtext("(no shift by construction)", side = 3, line = 0.2, cex = 0.6, col = "grey45")
    }
    for (i in seq_len(nr)) mtext(marks[[i]]$name, side = 2, outer = TRUE, line = 0.5, at = 1 - (i - 0.5)/nr, cex = 0.85, font = 2)
    if (!is.null(reference_shift)) mtext(sprintf("dashed red = IF global log2FC = %+.2f (treated vs reference); the method whose cloud sits on it matches the orthogonal truth", reference_shift), side = 3, outer = TRUE, line = 0.4, cex = 0.7, col = "#C0392B")
  }
  save_fig(paste0(tag, "_norm_comparison_scatter"), outdir, draw, w = max(5, 3.2*nc), h = max(4, 3.2*nr))
}
comparison_box_fig <- function(test_sig, ortho, methods, g1key, g2key, g1lab, g2lab, mark_test, mark_ortho, reference_shift, tag, outdir) {
  marks <- .assemble_marks(test_sig, ortho, mark_test, mark_ortho); nr <- length(marks); nc <- length(methods)
  draw <- function() {
    par(mfrow = c(nr, nc), mar = c(6, 4, 3.5, 1), oma = c(0, 2.2, 0.5, 0))
    for (mk in marks) for (m in methods) {
      S <- mk$sig[[m]]
      if (is.null(S) || !all(c(g1key, g2key) %in% colnames(S))) { plot.new(); box(col = "grey85"); title(main = m); text(0.5, 0.5, "n/a", col = "grey60"); next }
      ref <- S[, g2key]; trt <- S[, g1key]
      boxplot(list(ref, trt), names = c(g2lab, g1lab), col = c("#BFC9CA", "#E08A3C"), ylab = "norm logCPM", main = m, las = 2, outline = FALSE)
      if (mk$is_test && !is.null(reference_shift)) { expct <- median(ref, na.rm = TRUE) + reference_shift; segments(1.6, expct, 2.4, expct, col = "#C0392B", lwd = 2, lty = 2); text(2, expct, "IF", pos = 3, col = "#C0392B", cex = 0.7) }
      if (m == "libsize") mtext("(no shift by construction)", side = 3, line = 0.2, cex = 0.6, col = "grey45")
    }
    for (i in seq_len(nr)) mtext(marks[[i]]$name, side = 2, outer = TRUE, line = 0.5, at = 1 - (i - 0.5)/nr, cex = 0.85, font = 2)
  }
  save_fig(paste0(tag, "_norm_comparison_box"), outdir, draw, w = max(5, 3.0*nc), h = max(4, 3.2*nr))
}

# =============================================================================
# Orthogonal control-mark panel data (single representative width = win_sizes[1]).
#   Data-derived methods (libsize/background/greenlist) are recomputed from the
#   control mark; spikein only if its spikein_reads are in the sheet; sf skipped
#   (test-mark-specific). NULL with a notice if its samples are insufficient.
# =============================================================================
mark_panel_data <- function(mark_name, g1key, g2key, methods, ss, opt, param, green_gr, W) {
  mk <- tolower(trimws(mark_name)); sub <- ss[ss$.fac == mk, , drop = FALSE]
  if (!nrow(sub)) { say(sprintf("[ORTHO] no samples for control mark '%s' — control row skipped.", mark_name)); return(NULL) }
  sub$.key <- paste(sub$.cond, sub$.treat, sep = "|"); sub <- sub[sub$.key %in% c(g1key, g2key), , drop = FALSE]
  n1 <- sum(sub$.key == g1key); n2 <- sum(sub$.key == g2key)
  if (n1 < 1 || n2 < 1) { say(sprintf("[ORTHO] control mark '%s' lacks samples for both conditions (%s:%d, %s:%d) — control row skipped (flag accepted).", mark_name, g1key, n1, g2key, n2)); return(NULL) }
  if (any(!file.exists(sub$bamReads))) { say(sprintf("[ORTHO] control mark '%s' BAM(s) missing — control row skipped.", mark_name)); return(NULL) }
  say(sprintf("[ORTHO] control mark '%s': %d samples — building control panels at width %d.", mark_name, nrow(sub), W))
  bams <- sub$bamReads
  bins_o <- windowCounts(bams, bin = TRUE, width = opt$bg_bin, param = param)
  win    <- count_filter_one(bams, W, opt, param, bins_o, igg_bams = NULL, bin_control = NULL)
  y0 <- asDGEList(win); colnames(y0) <- sub$SampleID; cond_key <- sub$.key
  y_ref_o <- asDGEList(bins_o); colnames(y_ref_o) <- sub$SampleID
  green <- if ("greenlist" %in% methods && !is.null(green_gr)) regionCounts(bams, regions = green_gr, param = param) else NULL
  sv <- NULL
  if ("spikein" %in% methods && "spikein_reads" %in% colnames(sub)) { cand <- suppressWarnings(as.numeric(sub$spikein_reads)); names(cand) <- sub$SampleID; if (!any(is.na(cand)) && all(cand > 0)) sv <- cand else say("[ORTHO] control-mark spikein_reads missing — spikein panel skipped for control mark.") }
  out <- setNames(vector("list", length(methods)), methods)
  for (m in methods) {
    if (m == "sf") { out[[m]] <- NULL; next }
    if (m == "spikein" && is.null(sv)) { out[[m]] <- NULL; next }
    if (m == "greenlist" && is.null(green)) { out[[m]] <- NULL; next }
    nf <- tryCatch(compute_nf(m, y_ref_o, bins_o, green, sv, NULL), error = function(e) NULL)
    if (is.null(nf)) { out[[m]] <- NULL; next }
    y <- y0; y$samples$norm.factors <- nf[colnames(y)]; out[[m]] <- condition_signal(y, cond_key)
  }
  list(signal = out, g1 = g1key, g2 = g2key)
}

# =============================================================================
# Write per-method region tables; return counts for the summary.
# =============================================================================
write_outputs <- function(df, tag, norm, outdir, opt) {
  labs <- dir_labels(opt)
  df$call <- ifelse(df$FDR < opt$fdr & abs(df$logFC) > opt$lfc, ifelse(df$logFC > 0, labs$up, labs$down), "ns")
  all_path <- file.path(outdir, sprintf("%s_%s_all_regions.csv", tag, norm)); sig_path <- file.path(outdir, sprintf("%s_%s_significant.csv", tag, norm))
  write.csv(df, all_path, row.names = FALSE)
  sig <- df[df$call != "ns", , drop = FALSE]; sig <- sig[order(sig$FDR), , drop = FALSE]; write.csv(sig, sig_path, row.names = FALSE)
  n_up <- sum(sig$logFC > 0); n_dn <- sum(sig$logFC < 0)
  say(sprintf("[%s] regions: %d total, %d significant (%s=%d, %s=%d). -> %s", norm, nrow(df), nrow(sig), labs$up, n_up, labs$down, n_dn, basename(sig_path)))
  c(total = nrow(df), sig = nrow(sig), up = n_up, down = n_dn)
}

# =============================================================================
# Analyze one normalization across the window ladder, then recount in regions.
# =============================================================================
analyze_method <- function(method, meta, opt, param, igg_bams, peaks_gr, bins, bin_control, green, indep_bams, nf, win_sizes, design, tested, cond_key) {
  per_size <- list()
  for (W in win_sizes) {                                   # test each width independently
    win <- count_filter_one(meta$bamReads, W, opt, param, bins, igg_bams, bin_control)
    y <- asDGEList(win); colnames(y) <- meta$SampleID; y$samples$norm.factors <- nf[colnames(y)]
    y <- estimateDisp(y, design); fit <- glmQLFit(y, design, robust = TRUE); res <- glmQLFTest(fit, coef = tested)
    per_size[[as.character(W)]] <- list(win = win, tab = res$table)
    say(sprintf("   [window %5d] %d enriched windows tested", W, nrow(res$table)))
  }
  cons <- consolidate(per_size, opt); regions <- cons$regions; comb <- cons$combined
  say(sprintf("   consolidated -> %d regions (combined FDR)", length(regions)))

  rc <- regionCounts(meta$bamReads, regions = regions, param = param)   # recount reads in the final regions
  yreg <- asDGEList(rc); colnames(yreg) <- meta$SampleID; yreg$samples$norm.factors <- nf[colnames(yreg)]
  yreg$samples$cond <- meta$.cond; yreg$samples$treat <- meta$.treat; yreg$samples$batch <- meta$Batch
  if (opt$test == "contrast") { yreg$samples$group_role <- meta$.role; yreg$samples$peak_ref <- NA } else { yreg$samples$peak_ref <- attr(meta, "peak_ref"); yreg$samples$group_role <- NA }
  eff <- yreg$samples$lib.size * yreg$samples$norm.factors
  logcpm_all <- aveLogCPM(yreg$counts, lib.size = eff)
  ref_idx <- if (opt$test == "contrast") which(yreg$samples$group_role == "group2") else which(yreg$samples$peak_ref)
  baseline <- aveLogCPM(yreg$counts[, ref_idx, drop = FALSE], lib.size = eff[ref_idx])
  indep_lc <- NULL
  if (!is.null(indep_bams)) { irc <- regionCounts(indep_bams, regions = regions, param = param); indep_lc <- aveLogCPM(assay(irc), lib.size = irc$totals) }

  rep_lfc <- if ("rep.logFC" %in% colnames(comb)) comb$rep.logFC else comb$logFC   # representative (best-window) logFC
  df <- data.frame(chr = as.character(seqnames(regions)), start = start(regions), end = end(regions),
                   logFC = rep_lfc, logCPM = logcpm_all, PValue = comb$PValue, FDR = comb$FDR, direction = comb$direction,
                   baseline_logCPM = baseline, baseline_indep_logCPM = if (!is.null(indep_lc)) indep_lc else rep(NA_real_, length(regions)),
                   nWindows = if ("num.tests" %in% colnames(comb)) comb$num.tests else NA_integer_, stringsAsFactors = FALSE)
  df$baseline_stratum <- cut_strata(df$baseline_logCPM, opt$n_strata)
  df$baseline_indep_stratum <- cut_strata(df$baseline_indep_logCPM, opt$n_strata)
  df$in_peak <- if (!is.null(peaks_gr)) overlapsAny(regions, peaks_gr) else NA
  signal <- if (opt$test == "contrast") condition_signal(yreg, cond_key) else NULL
  list(df = df, yreg = yreg, signal = signal)
}

run_norm <- function(method, meta, opt, param, igg_bams, peaks_gr, bins, bin_control, green, indep_bams, nf, win_sizes, design, tested, tag, outdir, cond_key) {
  say(sprintf("\n========== NORMALIZATION: %s ==========", method))
  a <- analyze_method(method, meta, opt, param, igg_bams, peaks_gr, bins, bin_control, green, indep_bams, nf, win_sizes, design, tested, cond_key)
  df <- a$df; yreg <- a$yreg
  counts <- write_outputs(df, tag, method, outdir, opt)
  std_plots(yreg, NULL, df, tag, method, outdir, opt)
  sd_block <- balance_block(df, "baseline_stratum", opt, "samedata(RTM-prone)")
  in_block <- if (!all(is.na(df$baseline_indep_logCPM))) balance_block(df, "baseline_indep_stratum", opt, "indep") else NULL
  pk_block <- balance_peak(df, opt)
  strata_balance_fig(df, tag, method, outdir, opt); rtm_verdict(sd_block, pk_block, in_block, opt)
  list(df = df, counts = counts, nf = nf, eff_libsize = yreg$samples$lib.size * nf[colnames(yreg)], signal = a$signal)
}

# =============================================================================
# MAIN
# =============================================================================
main <- function() {
  opt <- parse_cli(); dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

  # ---- window ladder ----
  win_sizes <- as.integer(trimws(strsplit(opt$window, ",")[[1]])); win_sizes <- sort(unique(win_sizes[!is.na(win_sizes) & win_sizes > 0]))
  if (!length(win_sizes)) stop("--window must be one or more positive integers.")
  if (opt$spacing > min(win_sizes)) opt$spacing <- min(win_sizes)        # avoid gaps in the smallest scan

  # ---- normalization list ----
  methods_req <- tolower(trimws(strsplit(opt$norm, ",")[[1]])); methods_req <- unique(methods_req[nzchar(methods_req)])
  allowed <- c("libsize", "background", "greenlist", "spikein", "sf"); bad <- setdiff(methods_req, allowed)
  if (length(bad)) stop("Unknown --norm token(s): ", paste(bad, collapse = ", "), ". Allowed: ", paste(allowed, collapse = ", "))
  if ("greenlist" %in% methods_req && is.null(opt$greenlist))   stop("--norm includes 'greenlist' but --greenlist BED not given.")
  if ("spikein"   %in% methods_req && is.null(opt$spikein))     stop("--norm includes 'spikein' but --spikein not given.")
  if ("sf"        %in% methods_req && is.null(opt$normfactors)) stop("--norm includes 'sf' but --normfactors not given.")

  ss <- load_sheet(opt$sheet); meta <- select_samples(ss, opt)
  tag <- if (opt$test == "contrast") sprintf("%s_%s_vs_%s", opt$mark, safe_token(opt$group1), safe_token(opt$group2)) else sprintf("%s_interaction", opt$mark)
  say(sprintf("=== csaw v7 | mark=%s | test=%s | tag=%s ===", opt$mark, opt$test, tag))
  say(sprintf("[SAMPLES] %d libraries: %s", nrow(meta), paste(meta$SampleID, collapse = ", ")))
  say(sprintf("[WINDOW] ladder=%s; spacing=%d; max region width=%d", paste(win_sizes, collapse = ","), opt$spacing, opt$max_width))

  ab <- assess_batch(meta, opt)
  say(sprintf("[BATCH] status=%s; tested coef=%s; design cols: %s", ab$status, paste(ab$tested, collapse = "+"), paste(colnames(ab$design), collapse = ", ")))

  methods <- methods_req
  if (!ab$allow_spikein && !opt$force_norm) {
    drop <- intersect(c("spikein", "sf"), methods); if (length(drop)) say(sprintf("[BATCH] dropping %s (batch %s, not separable; override with --force_norm).", paste(drop, collapse = ","), ab$status)); methods <- setdiff(methods, drop)
  } else if (!ab$allow_spikein && opt$force_norm) say("[BATCH] --force_norm: applying spikein/sf despite non-separable batch.")
  if (!"libsize" %in% methods) { say("[NORM] forcing 'libsize' as the assumption-free baseline."); methods <- c("libsize", methods) }
  methods <- unique(methods); say(sprintf("[NORM] methods to run: %s", paste(methods, collapse = ", ")))

  param    <- make_param(opt)
  igg_bams <- resolve_igg(opt, meta); peaks_gr <- resolve_peaks(opt, meta)
  green_gr <- if ("greenlist" %in% methods) resolve_greenlist(opt) else NULL
  spike_v  <- if ("spikein"   %in% methods) resolve_spikein(opt, meta) else NULL
  sf_v     <- if ("sf"        %in% methods) resolve_normfactors(opt, meta) else NULL

  # ---- shared, width-independent counts ----
  bins <- windowCounts(meta$bamReads, bin = TRUE, width = opt$bg_bin, param = param); colnames(bins) <- meta$SampleID
  bin_control <- if (!is.null(igg_bams)) windowCounts(igg_bams, bin = TRUE, width = opt$bg_bin, param = param) else NULL
  green <- if ("greenlist" %in% methods) { g <- regionCounts(meta$bamReads, regions = green_gr, param = param); colnames(g) <- meta$SampleID; g } else NULL
  y_ref <- asDGEList(bins)                                 # colnames=SampleID, lib.size=library totals (for compute_nf)

  indep_bams <- NULL
  if (!is.null(opt$baseline_group)) {
    bg <- parse_group(opt$baseline_group); bkey <- paste(bg$cond, bg$treat, sep = "|")
    brow <- ss[ss$.fac == tolower(opt$mark) & paste(ss$.cond, ss$.treat, sep = "|") == bkey, , drop = FALSE]
    if (!nrow(brow)) stop("--baseline_group '", opt$baseline_group, "' matched no samples.")
    if (any(!file.exists(brow$bamReads))) stop("baseline_group BAM(s) missing.")
    say(sprintf("[BASELINE] independent stratifier '%s': %d samples (recounted per method in the consolidated regions).", opt$baseline_group, nrow(brow)))
    indep_bams <- brow$bamReads
  }
  cond_key <- if (opt$test == "contrast") meta$.key else paste(meta$.cond, meta$.treat, sep = "|")

  # ---- run each normalization (factors computed once, applied to every width) ----
  df_by_method <- list(); counts_by_method <- list(); nf_by_method <- list(); eff_by_method <- list(); sig_by_method <- list()
  for (m in methods) {
    nf <- compute_nf(m, y_ref, bins, green, spike_v, sf_v)
    r <- run_norm(m, meta, opt, param, igg_bams, peaks_gr, bins, bin_control, green, indep_bams, nf, win_sizes, ab$design, ab$tested, tag, opt$outdir, cond_key)
    df_by_method[[m]] <- r$df; counts_by_method[[m]] <- r$counts; nf_by_method[[m]] <- r$nf; eff_by_method[[m]] <- r$eff_libsize; sig_by_method[[m]] <- r$signal
  }

  norm_compare_fig(df_by_method, tag, opt$outdir)

  # ---- consolidated per-sample factor table ----
  nf_df <- data.frame(SampleID = meta$SampleID, Condition = as.character(meta$Condition), Treatment = as.character(meta$Treatment), Batch = meta$Batch, lib.size = y_ref$samples$lib.size, stringsAsFactors = FALSE)
  for (m in methods) { nf_df[[paste0("nf_", m)]] <- nf_by_method[[m]]; nf_df[[paste0("efflib_", m)]] <- eff_by_method[[m]] }
  nf_path <- file.path(opt$outdir, sprintf("%s_normfactors_used.csv", tag)); write.csv(nf_df, nf_path, row.names = FALSE)
  say(sprintf("[NORMFACTORS] per-sample norm.factor + effective library size per method -> %s", basename(nf_path)))

  # ---- normalization-comparison figure (contrast only) ----
  if (opt$test == "contrast") {
    g1key <- attr(meta, "g1key"); g2key <- attr(meta, "g2key")
    ortho <- if (!is.null(opt$orthogonal_mark)) mark_panel_data(opt$orthogonal_mark, g1key, g2key, methods, ss, opt, param, green_gr, win_sizes[1]) else NULL
    comparison_scatter_fig(sig_by_method, ortho, methods, g1key, g2key, opt$group1, opt$group2, opt$mark, opt$orthogonal_mark, opt$reference_shift, tag, opt$outdir)
    comparison_box_fig(sig_by_method, ortho, methods, g1key, g2key, opt$group1, opt$group2, opt$mark, opt$orthogonal_mark, opt$reference_shift, tag, opt$outdir)
    say("[COMPARE] normalization-comparison scatter + box written.")
    if (is.null(opt$reference_shift)) say("[COMPARE] NOTE: no --reference_shift; panels show each method's ASSUMED global shift, not a validated magnitude. Anchor with the IF log2FC.")
  } else say("[COMPARE] normalization-comparison figure is contrast-only; skipped for interaction. Rerun --test contrast to compare global-shift recovery.")

  # ---- summary ----
  summ <- file.path(opt$outdir, sprintf("%s_summary.txt", tag))
  hdr <- c(strrep("=", 72), sprintf("csaw v7 summary | %s", tag), strrep("=", 72), "",
           sprintf("window ladder: %s ; consolidated FDR (equiweight across widths)", paste(win_sizes, collapse = ",")),
           sprintf("Significant region counts per normalization (FDR<%.2g, |logFC|>%.2g):", opt$fdr, opt$lfc))
  rows <- vapply(methods, function(m) { cc <- counts_by_method[[m]]; sprintf("  %-12s total=%-7d sig=%-6d up=%-6d down=%-6d", m, cc["total"], cc["sig"], cc["up"], cc["down"]) }, character(1))
  writeLines(c(hdr, rows, "", "--- full run log ---", log_lines$txt), summ)
  say(sprintf("\n[DONE] summary -> %s", summ))
}

main()