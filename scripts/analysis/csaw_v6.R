#!/usr/bin/env Rscript
# =============================================================================
# csaw.R  (v6) — Window-based differential binding for BROAD histone marks
# =============================================================================
#
# PROJECT: A375 melanoma parental-vs-B10 constricted-migration multi-omics.
#
# WHAT THIS SCRIPT DOES --------------------------------------------------------
#   Broad chromatin marks (H3K27me3, H3K9me3, H3K4me1) form diffuse domains, not
#   sharp peaks. We TILE the genome into fixed windows and test each window for
#   differential coverage with edgeR's negative-binomial quasi-likelihood (QL)
#   model; replicates supply the within-group variance that turns a fold-change
#   into an FDR-controlled p-value. Adjacent windows are merged into regions and
#   their p-values combined (Simes -> region-level BH FDR), exactly the
#   FDR-correct region-level workflow of Lun & Smyth (csaw). Punctate H3K27Ac is
#   handled separately by diffbind.R.
#
# WHAT CHANGED FROM v5 (read this) --------------------------------------------
#   1. NORMALIZATION IS NOW EXPLICIT AND CORRECTLY NAMED.
#        v5 called its only data-free normalization "libsize", but it actually
#        ran normFactors(bins) = TMM-on-bins = BACKGROUND normalization. v6 keeps
#        BOTH and names them correctly:
#          libsize    = TRUE depth-only (norm.factors = 1). No composition/
#                       "majority-unchanged" assumption. Cannot show a global
#                       shift by construction (it forces equal totals).
#          background = TMM on genome-wide bins (the old "libsize"). Assumes most
#                       bins are unchanged -> a real caveat for a broad mark under
#                       a global shift, because the bins are full of the mark.
#      Pick methods with --norm (a COMMA-SEPARATED LIST, not an integer):
#          --norm libsize,background          (default; both need no extra data)
#          --norm libsize,background,greenlist  ... etc.
#   2. GREENLIST normalization added (--greenlist BED). Same TMM machinery as
#      background, but restricted to de Mello's pre-validated, biology-free
#      regions -> the "stable majority" TMM wants, cleaner for a broad mark.
#   3. SPIKE-IN (--spikein) and ChIPseqSpikeInFree SF (--normfactors) unchanged.
#   4. TEST-AWARE LABELS. For --test contrast, up/down = "gained"/"lost" mark.
#      For --test interaction, the tested coefficient is a DIFFERENCE OF
#      DIFFERENCES, so up/down are relabeled "stronger loss in <cond_ref>" /
#      "stronger loss in <cond_alt>" and the flattening/RTM verdict (a treatment
#      concept) is SUPPRESSED for interaction.
#   5. CONSOLIDATED PER-SAMPLE FACTOR TABLE. normfactors_used.csv lists, for every
#      sample, the norm.factor AND effective library size under EACH method run.
#   6. NEW NORMALIZATION-COMPARISON FIGURE (Egan et al. 2016, Fig 5 analog).
#      Scatter (reference vs treated signal, one panel per method) + box plots
#      (signal per condition, one panel per method) showing how much GLOBAL shift
#      each normalization implies. CRITICAL: the amount of shift a method shows is
#      a restatement of its scaling ASSUMPTION, not its correctness — so the
#      figure is only interpretable when anchored:
#        --reference_shift   draws the IF/western-measured global log2FC as a
#                            reference band (the independent ground truth).
#        --orthogonal_mark   re-runs the comparison on a control mark expected NOT
#                            to change (e.g. H3K9me3); a shift in the test mark but
#                            not the control mark = mark-specific = real biology.
#      The libsize panel is annotated "(no shift by construction)" so a flat box is
#      never misread as "no biology". Contrast-only (a 2-condition treated-vs-
#      reference axis); skipped with a notice for --test interaction.
#
# TWO TEST MODES (--test) -----------------------------------------------------
#   contrast    (default) Two-group comparison:
#                 --group1 NUMERATOR --group2 REFERENCE  ("Condition|Treatment")
#                 Migration:     --group1 "B10|none"     --group2 "Parental|none"
#                 UNC in B10:    --group1 "B10|UNC"      --group2 "B10|none"
#                 UNC in parent: --group1 "Parental|UNC" --group2 "Parental|none"
#   interaction  Factorial Condition*Treatment; tests the INTERACTION coefficient
#                 (does UNC's effect differ between B10 and Parental?). Needs all
#                 four cells. Matching is CASE-INSENSITIVE.
#
# OUTPUT TAG ------------------------------------------------------------------
#   contrast    -> <mark>_<group1>_vs_<group2>   (e.g. H3K27me3_B10UNC_vs_B10none)
#   interaction -> <mark>_interaction
#
# SAMPLE SHEET ----------------------------------------------------------------
#   Required : SampleID, Factor, Condition, Treatment, Replicate, Batch, bamReads
#   Optional : spikein_reads (--spikein sheet), Peaks (--peaks sheet),
#              bamControl (--igg sheet)
#   Allowed  : Condition Parental|B10 ; Treatment none|UNC (case-insensitive)
#   NOTE     : Condition/Treatment columns are the SOURCE OF TRUTH for grouping.
#              Sample NAMES are only labels; keep them consistent to avoid swaps.
#
# EXAMPLE RUNS ----------------------------------------------------------------
#   # (1) Migration contrast, default normalizations (libsize + background):
#   Rscript csaw.R --sheet ss.csv --mark H3K27me3 \
#       --group1 "B10|none" --group2 "Parental|none" --peaks sheet \
#       --outdir results/diff/H3K27me3_migration
#
#   # (2) UNC reversal in B10 — flagship for the compartment overlay. All five
#   #     normalizations, RTM-immune peak stratifier, an independent baseline
#   #     group, the IF anchor and the (future) orthogonal control mark:
#   Rscript csaw.R --sheet ss.csv --mark H3K27me3 \
#       --group1 "B10|UNC" --group2 "B10|none" \
#       --norm libsize,background,greenlist,spikein,sf \
#       --greenlist data/greenlist_hg19.bed \
#       --spikein sheet --normfactors data/spikeinfree_SF.txt \
#       --peaks sheet --baseline_group "Parental|none" \
#       --reference_shift -2.32 --orthogonal_mark H3K9me3 \
#       --outdir results/diff/H3K27me3_B10unc
#
#   # (3) Interaction (no comparison figure; flattening verdict suppressed):
#   Rscript csaw.R --sheet ss.csv --mark H3K27me3 --test interaction \
#       --peaks sheet --outdir results/diff/H3K27me3_interaction
#
#   # (4) Comparison-figure-focused contrast (just the data-free pair + greenlist,
#   #     anchored to IF):
#   Rscript csaw.R --sheet ss.csv --mark H3K27me3 \
#       --group1 "B10|UNC" --group2 "B10|none" \
#       --norm libsize,background,greenlist --greenlist data/greenlist_hg19.bed \
#       --reference_shift -2.32 --outdir results/diff/H3K27me3_normcompare
#
# OUTPUTS (into --outdir) -----------------------------------------------------
#   <tag>_<norm>_all_regions.csv / _significant.csv
#   <tag>_<norm>_{MDS,MDS_batch,MA,volcano,repcor,strata_balance}.{pdf,png}
#   <tag>_<norm>_norm_compare.{pdf,png}     (per-region logFC vs libsize)
#   <tag>_norm_comparison_scatter.{pdf,png} (NEW; contrast only)
#   <tag>_norm_comparison_box.{pdf,png}     (NEW; contrast only)
#   <tag>_normfactors_used.csv              (NEW; per-sample factor per method)
#   <tag>_summary.txt
# =============================================================================

suppressWarnings(suppressMessages({
  library(optparse)        # command-line option parsing
  library(csaw)            # windowCounts / regionCounts / normFactors / mergeWindows
  library(edgeR)           # NB-QL GLM, dispersion, DGEList, cpm
  library(GenomicRanges)   # GRanges for peaks / greenlist / windows
  library(S4Vectors)       # low-level Bioconductor helpers
}))

# =============================================================================
# CLI — every option carries its allowed values / format in the help string
# =============================================================================
parse_cli <- function() {
  option_list <- list(
    make_option("--sheet", type = "character",
                help = "REQUIRED. Sample sheet CSV (cols: SampleID, Factor, Condition, Treatment, Replicate, Batch, bamReads; + optional spikein_reads / Peaks / bamControl)."),
    make_option("--mark", type = "character",
                help = "REQUIRED. Factor/antibody, matched case-insensitively to the Factor column (e.g. H3K27me3 | H3K9me3 | H3K4me1)."),

    make_option("--test", type = "character", default = "contrast",
                help = "Test mode. ALLOWED: 'contrast' | 'interaction'. [default %default]"),

    # ---- contrast specification ----
    make_option("--group1", type = "character", default = "B10|none",
                help = "[contrast] Numerator 'Condition|Treatment' (treated). One of Parental|none | B10|none | Parental|UNC | B10|UNC. [default %default]"),
    make_option("--group2", type = "character", default = "Parental|none",
                help = "[contrast] Reference 'Condition|Treatment'. Positive logFC = higher in group1. [default %default]"),

    # ---- interaction specification ----
    make_option("--cond_ref",  type = "character", default = "Parental",
                help = "[interaction] Reference Condition. ALLOWED: 'Parental' | 'B10'. [default %default]"),
    make_option("--treat_ref", type = "character", default = "none",
                help = "[interaction] Reference Treatment. ALLOWED: 'none' | 'UNC'. [default %default]"),
    make_option("--cond_alt",  type = "character", default = "B10",
                help = "[interaction] Non-reference Condition. [default %default]"),
    make_option("--treat_alt", type = "character", default = "UNC",
                help = "[interaction] Non-reference Treatment. [default %default]"),

    make_option("--outdir", type = "character", default = "csaw_out",
                help = "Output directory. [default %default]"),

    # ---- NORMALIZATION SELECTOR (string list, NOT an integer) ----
    make_option("--norm", type = "character", default = "libsize,background",
                help = paste0("Comma-separated normalizations to run. ALLOWED tokens: ",
                              "libsize | background | greenlist | spikein | sf. ",
                              "libsize & background need no external data (defaults). ",
                              "greenlist needs --greenlist; spikein needs --spikein; sf needs --normfactors. ",
                              "[default %default]")),

    # ---- normalization data inputs ----
    make_option("--greenlist", type = "character", default = NULL,
                help = paste0("Greenlist normalization regions: a 3-column BED (chrom<TAB>start<TAB>end) of ",
                              "de Mello CUT&RUN greenlist regions IN YOUR GENOME BUILD (hg19; liftOver the ",
                              "published hg38 list first). TMM is computed on read counts in these regions. [default: none]")),
    make_option("--spikein", type = "character", default = NULL,
                help = "Spike-in (Orlando/RRPM) from READ COUNTS: norm.factor = spikein_reads/lib.size. 'sheet' (read spikein_reads column) OR a CSV with columns SampleID,spikein_reads. [default: none]"),
    make_option("--normfactors", type = "character", default = NULL,
                help = "Precomputed per-sample scaling factors (e.g. ChIPseqSpikeInFree). Tab/comma table with an ID column (SampleID/ID/Sample) and a factor column (SF/factor/ScalingFactor/norm.factors). 'QC' column checked; NA/non-positive rejected. [default: none]"),
    make_option("--force_norm", action = "store_true", default = FALSE,
                help = "Apply spikein/sf EVEN when batch is non-separable (confounded/partial). Only with an external magnitude anchor. [default: FALSE]"),

    # ---- NEW: comparison-figure anchors ----
    make_option("--reference_shift", type = "double", default = NULL,
                help = paste0("Externally measured global log2 fold-change of the mark, TREATED (group1) ",
                              "RELATIVE TO REFERENCE (group2), from IF/western/mass-spec. Negative = loss ",
                              "(e.g. a 5-fold reduction = log2(1/5) = -2.32). Drawn as a reference band on the ",
                              "normalization-comparison figure (the independent ground truth). [default: none]")),
    make_option("--orthogonal_mark", type = "character", default = NULL,
                help = paste0("Name of a CONTROL mark expected NOT to change under the perturbation (e.g. H3K9me3), ",
                              "matched to the Factor column. If its samples for BOTH contrast conditions are present ",
                              "in the sheet, the comparison figure gains a control-mark row (a shift in the test mark ",
                              "but not here = mark-specific = real). If the samples are absent, the flag is accepted ",
                              "and the control row is skipped with a notice. [default: none]")),

    # ---- counting / windows ----
    make_option("--window",  type = "integer", default = 2000L,
                help = "Window width (bp). Broad marks: 2000. [default %default]"),
    make_option("--spacing", type = "integer", default = 500L,
                help = "Spacing between window starts (bp); also the region-merge tolerance. [default %default]"),
    make_option("--bg_bin",  type = "integer", default = 10000L,
                help = "Background bin width (bp) for TMM (background norm) + enrichment filter. [default %default]"),
    make_option("--frag_max", type = "integer", default = 700L,
                help = "Max fragment length (bp); match Bowtie2 -X. [default %default]"),

    # ---- filtering / significance ----
    make_option("--min_lfc_bg", type = "double", default = log2(3),
                help = "Min log2 enrichment over background/IgG to KEEP a window. [default log2(3)=1.0]"),
    make_option("--fdr", type = "double", default = 0.05,
                help = "FDR cutoff for the significant subset. [default %default]"),
    make_option("--lfc", type = "double", default = 1.0,
                help = "Min |log2FC| for the significant subset. [default %default]"),
    make_option("--n_strata", type = "integer", default = 3L,
                help = "Number of baseline-abundance strata (>=2). 3 -> low/mid/high. [default %default]"),

    # ---- optional filters / stratifiers ----
    make_option("--blacklist", type = "character", default = NULL,
                help = "Optional ENCODE blacklist BED (e.g. hg19-blacklist.v2.bed). [default: none]"),
    make_option("--igg", type = "character", default = NULL,
                help = "IgG enrichment filter. 'sheet' (read bamControl column) OR comma-separated IgG BAM path(s). [default: none]"),
    make_option("--peaks", type = "character", default = NULL,
                help = "RTM-immune peak stratifier. 'sheet' (Peaks column of reference rows, unioned) OR comma-separated BED/narrowPeak/broadPeak/SEACR path(s). Use a BROAD call for H3K27me3. [default: none]"),
    make_option("--baseline_group", type = "character", default = NULL,
                help = "Independent baseline stratifier: a 'Condition|Treatment' group counted in the same windows, used ONLY to rank baseline (must be independent of the tested coefficient). [default: none]")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$sheet)) stop("--sheet is required.")
  if (is.null(opt$mark))  stop("--mark is required (e.g. H3K27me3).")
  if (!opt$test %in% c("contrast", "interaction"))
    stop("--test must be 'contrast' or 'interaction'.")
  if (opt$n_strata < 2L) stop("--n_strata must be >= 2.")
  opt
}

# Split a 'Condition|Treatment' spec into a lowercased list.
parse_group <- function(spec) {
  parts <- strsplit(spec, "\\|")[[1]]
  if (length(parts) != 2) stop("Group spec must be 'Condition|Treatment', got: ", spec)
  list(cond = tolower(trimws(parts[1])), treat = tolower(trimws(parts[2])))
}

# Filesystem-safe token (strip '|', spaces, etc.) for the output tag.
safe_token <- function(s) gsub("[^A-Za-z0-9]+", "", s)

# =============================================================================
# Logging — collect lines so we can also dump them into the summary report
# =============================================================================
log_lines <- new.env(); log_lines$txt <- character(0)
say <- function(...) { line <- paste0(...); cat(line, "\n"); log_lines$txt <- c(log_lines$txt, line) }

# Robust 3-column reader for BED / narrowPeak / broadPeak / SEACR / greenlist.
# BED is 0-based half-open; we add 1 to start to get 1-based GRanges coordinates.
read_bed3 <- function(path) {
  df <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                   comment.char = "#", quote = "")
  GRanges(seqnames = df[[1]],
          ranges = IRanges(start = as.integer(df[[2]]) + 1L, end = as.integer(df[[3]])))
}

# =============================================================================
# Sample sheet + selection
# =============================================================================
load_sheet <- function(path) {
  ss <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  req <- c("SampleID", "Factor", "Condition", "Treatment", "Replicate", "Batch", "bamReads")
  miss <- setdiff(req, colnames(ss))
  if (length(miss)) stop("Sample sheet missing columns: ", paste(miss, collapse = ", "))
  ss$.cond  <- tolower(trimws(ss$Condition))    # lowercased helpers for case-insensitive matching
  ss$.treat <- tolower(trimws(ss$Treatment))
  ss$.fac   <- tolower(trimws(ss$Factor))
  ss$Batch  <- as.character(ss$Batch)           # keep batch as a string, never numeric
  ss
}

# Pick the rows for the requested mark + test, build the model factors, check BAMs.
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
    # Reference level = group2, so a positive logFC means "higher in group1".
    sub$Group <- factor(sub$.key, levels = c(key(g2), key(g1)))
    if (any(is.na(sub$Group))) stop("Group level mismatch.")
    sub$.role <- ifelse(sub$.key == key(g1), "group1", "group2")
    attr(sub, "g1key") <- key(g1); attr(sub, "g2key") <- key(g2)
  } else {
    cr <- tolower(opt$cond_ref); ca <- tolower(opt$cond_alt)
    tr <- tolower(opt$treat_ref); ta <- tolower(opt$treat_alt)
    sub <- ss[ss$.cond %in% c(cr, ca) & ss$.treat %in% c(tr, ta), , drop = FALSE]
    if (!nrow(sub)) stop("No samples match the interaction cells.")
    sub$Condition <- factor(sub$.cond, levels = c(cr, ca))   # reference first
    sub$Treatment <- factor(sub$.treat, levels = c(tr, ta))
    if (any(is.na(sub$Condition)) || any(is.na(sub$Treatment)))
      stop("Interaction factor level mismatch.")
    # Reference rows for sheet-peaks: the cond_alt|treat_ref cell (e.g. B10|none).
    attr(sub, "peak_ref") <- sub$.cond == ca & sub$.treat == tr
  }
  if (any(!file.exists(sub$bamReads)))
    stop("Missing BAM(s): ", paste(sub$bamReads[!file.exists(sub$bamReads)], collapse = ", "))
  sub
}

# Which rows are the "reference" used for sheet-driven peaks / baseline ranking?
ref_rows <- function(meta, opt) {
  if (opt$test == "contrast") meta$.role == "group2" else attr(meta, "peak_ref")
}

# =============================================================================
# Resolve optional inputs (peaks / igg / spikein / sf / greenlist)
# =============================================================================
resolve_peaks <- function(opt, meta) {
  if (is.null(opt$peaks)) return(NULL)
  if (identical(opt$peaks, "sheet")) {
    if (!"Peaks" %in% colnames(meta)) stop("--peaks sheet requested but no 'Peaks' column.")
    paths <- unique(meta$Peaks[ref_rows(meta, opt)])
    paths <- paths[!is.na(paths) & nzchar(paths)]
    if (!length(paths)) stop("--peaks sheet: no Peaks paths on the reference rows.")
  } else paths <- trimws(strsplit(opt$peaks, ",")[[1]])
  miss <- paths[!file.exists(paths)]; if (length(miss)) stop("Peak file(s) not found: ", paste(miss, collapse = ", "))
  gr <- reduce(Reduce(c, lapply(paths, read_bed3)))   # union, then collapse overlaps
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

# Returns a named numeric vector of spike-in read counts keyed by SampleID, or NULL.
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

# Precomputed SF table (e.g. ChIPseqSpikeInFree). Flexible column/key detection.
resolve_normfactors <- function(opt, meta) {
  if (is.null(opt$normfactors)) return(NULL)
  if (!file.exists(opt$normfactors)) stop("--normfactors file not found: ", opt$normfactors)
  first <- readLines(opt$normfactors, n = 1L, warn = FALSE)
  sep   <- if (grepl("\t", first)) "\t" else ","                 # auto-detect tab vs comma
  tab   <- read.table(opt$normfactors, header = TRUE, sep = sep, quote = "\"",
                      stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE)
  fcol <- intersect(c("SF", "factor", "ScalingFactor", "norm.factors", "normfactor"), colnames(tab))
  if (!length(fcol)) stop("--normfactors: no scaling-factor column found. Columns: ", paste(colnames(tab), collapse = ", "))
  icol <- intersect(c("SampleID", "ID", "Sample", "sample"), colnames(tab))
  if (!length(icol)) stop("--normfactors: no ID column found. Columns: ", paste(colnames(tab), collapse = ", "))
  ids <- as.character(tab[[icol[1]]]); sf <- suppressWarnings(as.numeric(tab[[fcol[1]]]))
  bn  <- basename(meta$bamReads)
  out <- setNames(rep(NA_real_, nrow(meta)), meta$SampleID); qc_used <- character(0)
  for (i in seq_len(nrow(meta))) {                                # match SampleID, else BAM basename
    hit <- which(ids == meta$SampleID[i]); if (!length(hit)) hit <- which(ids == bn[i])
    if (!length(hit)) hit <- which(basename(ids) == bn[i])
    if (length(hit)) { out[i] <- sf[hit[1]]; if ("QC" %in% colnames(tab)) qc_used <- c(qc_used, as.character(tab$QC[hit[1]])) }
  }
  if (length(qc_used) && any(!grepl("pass", tolower(qc_used))))
    say("[NORMFACTORS] WARNING: a matched row has QC != 'pass' (ChIPseqSpikeInFree flags poor enrichment / complete loss). Verify before trusting it.")
  if (any(is.na(out))) stop("--normfactors: could not match a factor for: ", paste(meta$SampleID[is.na(out)], collapse = ", "))
  if (any(out <= 0))   stop("--normfactors: non-positive/NA factor for: ", paste(meta$SampleID[out <= 0], collapse = ", "))
  out
}

# Read the greenlist BED (only if greenlist is among the requested methods).
resolve_greenlist <- function(opt) {
  if (is.null(opt$greenlist)) return(NULL)
  if (!file.exists(opt$greenlist)) stop("--greenlist file not found: ", opt$greenlist)
  gr <- reduce(read_bed3(opt$greenlist))
  say(sprintf("[GREENLIST] %d normalizer regions loaded (assumed hg19; verify low signal in your IgG).", length(gr)))
  gr
}

# =============================================================================
# Batch assessment (per test): same / crossed / partial / confounded
#   The WHOLE design must be full rank (edgeR's requirement), not just the tested
#   coefficient; both are checked by a rank test. The gate it returns affects ONLY
#   spikein/sf (cross-batch carryover is not separable); libsize/background/
#   greenlist are immune to it.
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
# Counting + filtering
# =============================================================================
make_param <- function(opt) {
  discard <- GRanges()                                  # empty == readParam default (NOT NULL)
  if (!is.null(opt$blacklist)) {
    if (!file.exists(opt$blacklist)) stop("Blacklist not found: ", opt$blacklist)
    discard <- read_bed3(opt$blacklist)
  }
  readParam(minq = 20L, dedup = TRUE, pe = "both", max.frag = opt$frag_max, discard = discard)
}

# Count reads into windows (the test units) and into background bins (for the
# enrichment filter and the background-TMM normalization); return filtered windows.
count_and_filter <- function(bams, opt, param, igg_bams) {
  win  <- windowCounts(bams, width = opt$window, spacing = opt$spacing, param = param)
  bins <- windowCounts(bams, bin = TRUE, width = opt$bg_bin, param = param)
  if (!is.null(igg_bams)) {                              # IgG-aware filter (keep windows enriched over IgG)
    control     <- regionCounts(igg_bams, regions = rowRanges(win), param = param)
    bin.control <- windowCounts(igg_bams, bin = TRUE, width = opt$bg_bin, param = param)
    scinfo <- scaleControlFilter(bins, bin.control)
    fstat  <- filterWindowsControl(win, control, scale.info = scinfo)
    say("[FILTER] IgG-enrichment filter (filterWindowsControl).")
  } else {                                               # else filter over genomic background
    fstat <- filterWindowsGlobal(win, bins); say("[FILTER] genomic-background enrichment filter (filterWindowsGlobal).")
  }
  keep <- fstat$filter > opt$min_lfc_bg
  say(sprintf("[FILTER] kept %d / %d windows (>%.2f log2 over background).", sum(keep), length(keep), opt$min_lfc_bg))
  list(win = win[keep, ], bins = bins)
}

# =============================================================================
# Per-method normalization factors
#   Every method returns a per-sample norm.factors vector keyed by SampleID.
#   In edgeR the EFFECTIVE library size = lib.size * norm.factor; CPM/logFC use it.
# =============================================================================
compute_nf <- function(method, y0, win, bins, green, spike_v, sf_v) {
  cn <- colnames(y0)
  if (method == "libsize") {
    nf <- rep(1, ncol(y0))                               # TRUE depth-only: no composition correction
  } else if (method == "background") {
    nf <- asDGEList(normFactors(bins,  se.out = win))$samples$norm.factors  # TMM on genome-wide bins
  } else if (method == "greenlist") {
    if (is.null(green)) stop("[greenlist] requested but greenlist counts are missing.")
    nf <- asDGEList(normFactors(green, se.out = win))$samples$norm.factors  # TMM on greenlist regions
  } else if (method == "spikein") {
    sr <- spike_v[cn]; nf <- sr / y0$samples$lib.size    # Orlando/RRPM: scale by exogenous reads
    nf <- nf / exp(mean(log(nf)))                        # geometric-mean centre to 1
  } else if (method == "sf") {
    v <- sf_v[cn]; nf <- v / exp(mean(log(v)))           # injected ChIPseqSpikeInFree SF, centred
  } else stop("Unknown normalization method: ", method)
  setNames(as.numeric(nf), cn)
}

# =============================================================================
# Fit + test + region-merge, then attach baseline / peak strata
# =============================================================================
fit_test_regions <- function(y, ab, opt, win_gr, peaks_gr, indep_logcpm) {
  y   <- estimateDisp(y, ab$design)                      # NB dispersion, trended + tagwise (EB-squeezed)
  fit <- glmQLFit(y, ab$design, robust = TRUE)           # QL fit; robust = guard outlier WINDOWS
  res <- glmQLFTest(fit, coef = ab$tested)               # QL F-test on the tested coefficient(s)
  wtab <- res$table                                      # per-window logFC, logCPM, PValue

  # Reference-group baseline abundance per window (RTM-prone 'samedata' stratifier).
  ref_idx <- if (opt$test == "contrast") which(y$samples$group_role == "group2") else which(y$samples$peak_ref)
  base_lcpm <- aveLogCPM(y$counts[, ref_idx, drop = FALSE],
                         lib.size = y$samples$lib.size[ref_idx] * y$samples$norm.factors[ref_idx])

  merged <- mergeWindows(win_gr, tol = opt$spacing, max.width = 10000L)   # adjacent windows -> regions
  comb   <- combineTests(merged$ids, wtab)               # Simes per region -> region BH FDR
  reg    <- merged$regions

  agg <- function(v) as.numeric(tapply(v, merged$ids, mean))  # mean over member windows
  df <- data.frame(
    chr = as.character(seqnames(reg)), start = start(reg), end = end(reg),
    logFC = agg(wtab$logFC), logCPM = agg(wtab$logCPM),
    PValue = comb$PValue, FDR = comb$FDR, direction = comb$direction,
    baseline_logCPM = agg(base_lcpm),
    baseline_indep_logCPM = if (!is.null(indep_logcpm)) agg(indep_logcpm) else rep(NA_real_, length(reg)),
    stringsAsFactors = FALSE)

  cut_strata <- function(x, k) {                          # low..high strata from baseline ranks
    if (all(is.na(x))) return(factor(rep(NA, length(x))))
    br <- quantile(x, probs = seq(0, 1, length.out = k + 1), na.rm = TRUE)
    lab <- if (k == 3) c("low", "mid", "high") else paste0("s", seq_len(k))
    cut(x, breaks = unique(br), labels = lab[seq_len(length(unique(br)) - 1)], include.lowest = TRUE)
  }
  df$baseline_stratum       <- cut_strata(df$baseline_logCPM, opt$n_strata)
  df$baseline_indep_stratum <- cut_strata(df$baseline_indep_logCPM, opt$n_strata)
  df$in_peak <- if (!is.null(peaks_gr)) overlapsAny(reg, peaks_gr) else NA
  list(table = df, res = res, fit = fit, y = y)
}

# =============================================================================
# Direction labels — TEST-AWARE so the interaction is never misread
#   contrast:    +logFC = "gained" mark,   -logFC = "lost" mark.
#   interaction: the coefficient is (UNC effect in cond_alt) - (UNC effect in
#                cond_ref), a difference of differences. So +/- mean which
#                condition LOSES MORE, not gained/lost mark.
# =============================================================================
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
# Stratified up/down balance + RTM verdict (verdict only meaningful for contrast)
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

# The flattening signature (loss at high-baseline peaks, spread to low) is a
# TREATMENT concept; only emit the verdict for a contrast test.
rtm_verdict <- function(sd_block, peak_block, indep_block, opt) {
  if (opt$test != "contrast") {
    say("[STRAT verdict] (interaction: flattening verdict suppressed — the coefficient is a difference of differences, not a treatment effect)")
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
  if (is.na(grad_pk) && is.na(grad_in)) say("   (no RTM-resistant stratifier supplied; samedata is RTM-prone — interpret with care)")
}

# =============================================================================
# Per-condition signal (for the normalization-comparison figure)
#   For a given normalized y and a per-sample condition key, returns a
#   windows x conditions matrix of mean normalized logCPM (replicates averaged).
# =============================================================================
condition_signal <- function(y, cond_key) {
  lc <- cpm(y, log = TRUE)                                # per-window normalized logCPM
  conds <- unique(cond_key)
  m <- sapply(conds, function(cc) rowMeans(lc[, cond_key == cc, drop = FALSE]))
  if (is.null(dim(m))) m <- matrix(m, ncol = length(conds))
  colnames(m) <- conds; m
}

# =============================================================================
# Standard per-method figures (MDS / MA / volcano / repcor / strata balance)
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
  mds_plot(y, tag, norm, outdir, "treat"); mds_plot(y, tag, norm, outdir, "batch")  # bio + batch-coloured MDS

  save_fig(paste0(tag, "_", norm, "_MA"), outdir, function() {                       # MA: logFC vs logCPM
    plot(df$logCPM, df$logFC, pch = 16, cex = 0.3, col = "#88888855", xlab = "logCPM", ylab = "logFC",
         main = paste(tag, norm, "MA")); sig <- df$FDR < opt$fdr & abs(df$logFC) > opt$lfc
    points(df$logCPM[sig], df$logFC[sig], pch = 16, cex = 0.35, col = "#C0392B"); abline(h = 0, col = "grey40")
  })

  save_fig(paste0(tag, "_", norm, "_volcano"), outdir, function() {                  # volcano with test-aware counts
    sig <- df$FDR < opt$fdr & abs(df$logFC) > opt$lfc
    plot(df$logFC, -log10(df$FDR), pch = 16, cex = 0.3, col = "#88888855", xlab = "logFC", ylab = "-log10 FDR",
         main = paste(tag, norm, "volcano"))
    points(df$logFC[sig], -log10(df$FDR[sig]), pch = 16, cex = 0.35, col = ifelse(df$logFC[sig] > 0, "#C0392B", "#1F6F8B"))
    abline(v = c(-opt$lfc, opt$lfc), h = -log10(opt$fdr), col = "grey60", lty = 2)
    n_up <- sum(sig & df$logFC > 0); n_dn <- sum(sig & df$logFC < 0)
    legend("top", horiz = TRUE, bg = "white", box.col = "grey80", cex = 0.95,
           legend = c(sprintf("%s: %d", labs$down, n_dn), sprintf("%s: %d", labs$up, n_up)), text.col = c("#1F6F8B","#C0392B"))
  })

  save_fig(paste0(tag, "_", norm, "_repcor"), outdir, function() {                   # replicate scatter matrix
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
# Cross-method per-region logFC concordance (the RELATIVE-structure robustness
#   check). All methods share identical windows -> identical regions in identical
#   order, so logFC vectors align row-for-row. A tight diagonal here means the
#   loci-level redistribution is the same regardless of the global-scale
#   assumption; only the genome-wide offset differs between methods.
# =============================================================================
norm_compare_fig <- function(df_by_method, tag, outdir) {
  methods <- names(df_by_method); if (length(methods) < 2) return(invisible(NULL))
  ref <- if ("libsize" %in% methods) "libsize" else methods[1]
  others <- setdiff(methods, ref)
  draw <- function() {
    par(mfrow = c(1, length(others)), mar = c(4,4,3,1))
    for (m in others) {
      x <- df_by_method[[ref]]$logFC; y <- df_by_method[[m]]$logFC
      rng <- range(c(x, y), finite = TRUE)
      plot(x, y, pch = 16, cex = 0.25, col = "#33333522", xlim = rng, ylim = rng,
           xlab = sprintf("logFC (%s)", ref), ylab = sprintf("logFC (%s)", m),
           main = sprintf("%s vs %s", m, ref))
      abline(0, 1, col = "grey40"); abline(h = 0, v = 0, col = "grey85", lty = 3)
      rho <- suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
      legend("topleft", bty = "n", cex = 0.85, legend = sprintf("Spearman rho = %.3f", rho))
    }
  }
  save_fig(paste0(tag, "_norm_compare"), outdir, draw, w = max(5, 3.2*length(others)), h = 3.8)
}

# =============================================================================
# NEW — Normalization-comparison figures (Egan et al. 2016, Fig 5 analog)
#   These visualize, for each normalization, HOW MUCH GLOBAL SHIFT it implies
#   between the reference (group2) and treated (group1) conditions.
#   READ THIS CAREFULLY: the size of the shift a method shows is a restatement
#   of its scaling ASSUMPTION, not evidence it is correct. The figures are
#   interpretable only when anchored:
#     - dashed red IF line / band  = the externally measured log2FC (ground truth)
#     - orthogonal control-mark row = a shift here too would mean the normalization
#                                     manufactures shifts (artifact); a shift in the
#                                     test mark but not the control = real biology.
#   The libsize panel is annotated "(no shift by construction)".
# =============================================================================
.assemble_marks <- function(test_sig, ortho, mark_test, mark_ortho) {
  marks <- list(list(name = mark_test, sig = test_sig, is_test = TRUE))
  if (!is.null(ortho) && !is.null(ortho$signal))
    marks <- c(marks, list(list(name = mark_ortho, sig = ortho$signal, is_test = FALSE)))
  marks
}

# Scatter: reference signal (x) vs treated signal (y) per window, one panel per
# method, rows = marks. Cloud below the diagonal = implied global loss.
comparison_scatter_fig <- function(test_sig, ortho, methods, g1key, g2key, g1lab, g2lab,
                                    mark_test, mark_ortho, reference_shift, tag, outdir) {
  marks <- .assemble_marks(test_sig, ortho, mark_test, mark_ortho)
  nr <- length(marks); nc <- length(methods)
  draw <- function() {
    par(mfrow = c(nr, nc), mar = c(4, 4, 3.5, 1), oma = c(0, 2.2, if (!is.null(reference_shift)) 2.2 else 0.5, 0))
    for (mk in marks) for (m in methods) {
      S <- mk$sig[[m]]
      if (is.null(S) || !all(c(g1key, g2key) %in% colnames(S))) {       # method not computable for this mark
        plot.new(); box(col = "grey85"); title(main = m); text(0.5, 0.5, "n/a", col = "grey60"); next
      }
      x <- S[, g2key]; y <- S[, g1key]; rng <- range(c(x, y), finite = TRUE)
      plot(x, y, pch = 16, cex = 0.25, col = "#33333522", xlim = rng, ylim = rng,
           xlab = sprintf("%s (norm logCPM)", g2lab), ylab = sprintf("%s (norm logCPM)", g1lab), main = m)
      abline(0, 1, col = "grey40", lwd = 1.2)                            # no-change diagonal
      if (mk$is_test && !is.null(reference_shift))
        abline(reference_shift, 1, col = "#C0392B", lwd = 1.6, lty = 2)  # IF ground truth: y = x + log2FC
      if (m == "libsize") mtext("(no shift by construction)", side = 3, line = 0.2, cex = 0.6, col = "grey45")
    }
    for (i in seq_len(nr)) mtext(marks[[i]]$name, side = 2, outer = TRUE, line = 0.5,
                                 at = 1 - (i - 0.5)/nr, cex = 0.85, font = 2)
    if (!is.null(reference_shift))
      mtext(sprintf("dashed red = IF/western global log2FC = %+.2f (treated vs reference); the method whose cloud sits on this line matches the orthogonal truth",
                    reference_shift), side = 3, outer = TRUE, line = 0.4, cex = 0.7, col = "#C0392B")
  }
  save_fig(paste0(tag, "_norm_comparison_scatter"), outdir, draw, w = max(5, 3.2*nc), h = max(4, 3.2*nr))
}

# Box: signal distribution per condition, one panel per method, rows = marks.
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
      if (mk$is_test && !is.null(reference_shift)) {                     # where treated SHOULD sit per IF
        expct <- median(ref, na.rm = TRUE) + reference_shift
        segments(1.6, expct, 2.4, expct, col = "#C0392B", lwd = 2, lty = 2)
        text(2, expct, "IF", pos = 3, col = "#C0392B", cex = 0.7)
      }
      if (m == "libsize") mtext("(no shift by construction)", side = 3, line = 0.2, cex = 0.6, col = "grey45")
    }
    for (i in seq_len(nr)) mtext(marks[[i]]$name, side = 2, outer = TRUE, line = 0.5,
                                 at = 1 - (i - 0.5)/nr, cex = 0.85, font = 2)
  }
  save_fig(paste0(tag, "_norm_comparison_box"), outdir, draw, w = max(5, 3.0*nc), h = max(4, 3.2*nr))
}

# =============================================================================
# Orthogonal control-mark panel data (for the comparison figure)
#   Counts the control mark's BAMs for the two contrast conditions, applies each
#   DATA-DERIVED normalization (libsize / background / greenlist) recomputed from
#   the control mark itself, plus spikein IF its spikein_reads are in the sheet.
#   (sf is test-mark-specific and is skipped for the control mark.) Returns
#   per-method windows x 2 condition-signal matrices, or NULL (with a notice) if
#   the control mark's samples for both conditions are not present.
# =============================================================================
mark_panel_data <- function(mark_name, g1key, g2key, methods, ss, opt, param, green_gr) {
  mk  <- tolower(trimws(mark_name))
  sub <- ss[ss$.fac == mk, , drop = FALSE]
  if (!nrow(sub)) { say(sprintf("[ORTHO] no samples for control mark '%s' — control row skipped.", mark_name)); return(NULL) }
  sub$.key <- paste(sub$.cond, sub$.treat, sep = "|")
  sub <- sub[sub$.key %in% c(g1key, g2key), , drop = FALSE]
  n1 <- sum(sub$.key == g1key); n2 <- sum(sub$.key == g2key)
  if (n1 < 1 || n2 < 1) {
    say(sprintf("[ORTHO] control mark '%s' lacks samples for both contrast conditions (%s:%d, %s:%d) — control row skipped (flag accepted; nothing else affected).",
                mark_name, g1key, n1, g2key, n2)); return(NULL)
  }
  if (any(!file.exists(sub$bamReads))) { say(sprintf("[ORTHO] control mark '%s' BAM(s) missing — control row skipped.", mark_name)); return(NULL) }
  say(sprintf("[ORTHO] control mark '%s': %d samples (%s:%d, %s:%d) — building control panels.", mark_name, nrow(sub), g1key, n1, g2key, n2))

  bams <- sub$bamReads
  cf   <- count_and_filter(bams, opt, param, igg_bams = NULL)   # global-background filter for the control mark
  win  <- cf$win; bins <- cf$bins
  y0   <- asDGEList(win); colnames(y0) <- sub$SampleID
  cond_key <- sub$.key

  green <- if ("greenlist" %in% methods && !is.null(green_gr)) regionCounts(bams, regions = green_gr, param = param) else NULL

  sv <- NULL                                                    # control-mark spike-in reads (only if present in sheet)
  if ("spikein" %in% methods && "spikein_reads" %in% colnames(sub)) {
    cand <- suppressWarnings(as.numeric(sub$spikein_reads)); names(cand) <- sub$SampleID
    if (!any(is.na(cand)) && all(cand > 0)) sv <- cand
    else say("[ORTHO] control-mark spikein_reads missing/invalid — spikein panel skipped for control mark.")
  }

  out <- setNames(vector("list", length(methods)), methods)
  for (m in methods) {
    if (m == "sf")                       { out[[m]] <- NULL; say("[ORTHO] 'sf' is test-mark-specific — skipped for control mark."); next }
    if (m == "spikein" && is.null(sv))   { out[[m]] <- NULL; next }
    if (m == "greenlist" && is.null(green)) { out[[m]] <- NULL; next }
    nf <- tryCatch(compute_nf(m, y0, win, bins, green, sv, NULL), error = function(e) { say(sprintf("[ORTHO] method '%s' failed for control mark: %s", m, conditionMessage(e))); NULL })
    if (is.null(nf)) { out[[m]] <- NULL; next }
    y <- y0; y$samples$norm.factors <- nf
    out[[m]] <- condition_signal(y, cond_key)                   # windows x 2 (ref, treated)
  }
  list(signal = out, g1 = g1key, g2 = g2key)
}

# =============================================================================
# Write per-method region tables; return counts for the summary.
# =============================================================================
write_outputs <- function(df, tag, norm, outdir, opt) {
  labs <- dir_labels(opt)
  df$call <- ifelse(df$FDR < opt$fdr & abs(df$logFC) > opt$lfc,
                    ifelse(df$logFC > 0, labs$up, labs$down), "ns")    # test-aware significance label
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

# =============================================================================
# Per-normalization driver: set factors -> fit/test/merge -> outputs + figures.
#   Returns the region table, counts, the norm.factors, the effective library
#   sizes, and (contrast only) the per-condition signal for the comparison figure.
# =============================================================================
run_norm <- function(method, y0, win, bins, green, spike_v, sf_v, peaks_gr, indep_logcpm,
                      win_gr, ab, tag, outdir, opt, cond_key) {
  say(sprintf("\n========== NORMALIZATION: %s ==========", method))
  nf <- compute_nf(method, y0, win, bins, green, spike_v, sf_v)   # per-sample norm.factors
  y  <- y0; y$samples$norm.factors <- nf                         # only the factors differ between methods
  ftr <- fit_test_regions(y, ab, opt, win_gr, peaks_gr, indep_logcpm)
  df  <- ftr$table; yfit <- ftr$y

  counts <- write_outputs(df, tag, method, outdir, opt)
  std_plots(yfit, ftr$res, df, tag, method, outdir, opt)

  sd_block <- balance_block(df, "baseline_stratum", opt, "samedata(RTM-prone)")
  in_block <- if (!all(is.na(df$baseline_indep_logCPM))) balance_block(df, "baseline_indep_stratum", opt, "indep") else NULL
  pk_block <- balance_peak(df, opt)
  strata_balance_fig(df, tag, method, outdir, opt)
  rtm_verdict(sd_block, pk_block, in_block, opt)

  sig <- if (opt$test == "contrast") condition_signal(yfit, cond_key) else NULL   # comparison-figure data
  list(df = df, counts = counts, nf = nf, eff_libsize = y0$samples$lib.size * nf, signal = sig)
}

# =============================================================================
# MAIN
# =============================================================================
main <- function() {
  opt <- parse_cli()
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

  # ---- parse + validate the requested normalization list ----
  methods_req <- tolower(trimws(strsplit(opt$norm, ",")[[1]])); methods_req <- unique(methods_req[nzchar(methods_req)])
  allowed <- c("libsize", "background", "greenlist", "spikein", "sf")
  bad <- setdiff(methods_req, allowed)
  if (length(bad)) stop("Unknown --norm token(s): ", paste(bad, collapse = ", "), ". Allowed: ", paste(allowed, collapse = ", "))
  if ("greenlist" %in% methods_req && is.null(opt$greenlist))   stop("--norm includes 'greenlist' but --greenlist BED not given.")
  if ("spikein"   %in% methods_req && is.null(opt$spikein))     stop("--norm includes 'spikein' but --spikein not given.")
  if ("sf"        %in% methods_req && is.null(opt$normfactors)) stop("--norm includes 'sf' but --normfactors not given.")

  ss   <- load_sheet(opt$sheet)
  meta <- select_samples(ss, opt)

  tag <- if (opt$test == "contrast")
    sprintf("%s_%s_vs_%s", opt$mark, safe_token(opt$group1), safe_token(opt$group2))
  else sprintf("%s_interaction", opt$mark)
  say(sprintf("=== csaw v6 | mark=%s | test=%s | tag=%s ===", opt$mark, opt$test, tag))
  say(sprintf("[SAMPLES] %d libraries: %s", nrow(meta), paste(meta$SampleID, collapse = ", ")))

  # ---- batch gate (affects spikein/sf only) ----
  ab <- assess_batch(meta, opt)
  say(sprintf("[BATCH] status=%s; tested coef=%s; design cols: %s",
              ab$status, paste(ab$tested, collapse = "+"), paste(colnames(ab$design), collapse = ", ")))

  methods <- methods_req
  if (!ab$allow_spikein && !opt$force_norm) {
    drop <- intersect(c("spikein", "sf"), methods)
    if (length(drop)) say(sprintf("[BATCH] dropping %s (batch %s, not separable; override with --force_norm).", paste(drop, collapse = ","), ab$status))
    methods <- setdiff(methods, drop)
  } else if (!ab$allow_spikein && opt$force_norm) {
    say("[BATCH] --force_norm: applying spikein/sf despite non-separable batch (justify with an external magnitude anchor).")
  }
  if (!"libsize" %in% methods) { say("[NORM] forcing 'libsize' to always run as the assumption-free baseline."); methods <- c("libsize", methods) }
  methods <- unique(methods)
  say(sprintf("[NORM] methods to run: %s", paste(methods, collapse = ", ")))

  # ---- resolve inputs ----
  param    <- make_param(opt)
  igg_bams <- resolve_igg(opt, meta)
  peaks_gr <- resolve_peaks(opt, meta)
  green_gr <- if ("greenlist" %in% methods) resolve_greenlist(opt) else NULL
  spike_v  <- if ("spikein"   %in% methods) resolve_spikein(opt, meta) else NULL
  sf_v     <- if ("sf"        %in% methods) resolve_normfactors(opt, meta) else NULL

  # ---- count + filter the TEST mark (shared across all normalizations) ----
  cf <- count_and_filter(meta$bamReads, opt, param, igg_bams)
  win <- cf$win; bins <- cf$bins; win_gr <- rowRanges(win)
  green <- if ("greenlist" %in% methods) regionCounts(meta$bamReads, regions = green_gr, param = param) else NULL

  # ---- optional independent baseline group (ranks baseline only) ----
  indep_logcpm <- NULL
  if (!is.null(opt$baseline_group)) {
    bg <- parse_group(opt$baseline_group); bkey <- paste(bg$cond, bg$treat, sep = "|")
    brow <- ss[ss$.fac == tolower(opt$mark) & paste(ss$.cond, ss$.treat, sep = "|") == bkey, , drop = FALSE]
    if (!nrow(brow)) stop("--baseline_group '", opt$baseline_group, "' matched no samples.")
    if (any(!file.exists(brow$bamReads))) stop("baseline_group BAM(s) missing.")
    say(sprintf("[BASELINE] independent stratifier '%s': %d samples.", opt$baseline_group, nrow(brow)))
    bc <- regionCounts(brow$bamReads, regions = win_gr, param = param)
    indep_logcpm <- aveLogCPM(assay(bc), lib.size = bc$totals)
  }

  # ---- base DGEList (only norm.factors differ per method) ----
  y0 <- asDGEList(win); colnames(y0) <- meta$SampleID
  y0$samples$cond  <- meta$.cond
  y0$samples$treat <- meta$.treat
  y0$samples$batch <- meta$Batch
  if (opt$test == "contrast") { y0$samples$group_role <- meta$.role; y0$samples$peak_ref <- NA }
  else                        { y0$samples$peak_ref <- attr(meta, "peak_ref"); y0$samples$group_role <- NA }
  cond_key <- if (opt$test == "contrast") meta$.key else paste(meta$.cond, meta$.treat, sep = "|")

  # ---- run each normalization ----
  df_by_method <- list(); counts_by_method <- list(); nf_by_method <- list(); eff_by_method <- list(); sig_by_method <- list()
  for (m in methods) {
    r <- run_norm(m, y0, win, bins, green, spike_v, sf_v, peaks_gr, indep_logcpm, win_gr, ab, tag, opt$outdir, opt, cond_key)
    df_by_method[[m]] <- r$df; counts_by_method[[m]] <- r$counts
    nf_by_method[[m]] <- r$nf; eff_by_method[[m]] <- r$eff_libsize; sig_by_method[[m]] <- r$signal
  }

  # ---- cross-method per-region logFC concordance (relative-structure robustness) ----
  norm_compare_fig(df_by_method, tag, opt$outdir)

  # ---- consolidated per-sample factor table ----
  nf_df <- data.frame(SampleID = meta$SampleID, Condition = as.character(meta$Condition),
                      Treatment = as.character(meta$Treatment), Batch = meta$Batch,
                      lib.size = y0$samples$lib.size, stringsAsFactors = FALSE)
  for (m in methods) { nf_df[[paste0("nf_", m)]] <- nf_by_method[[m]]; nf_df[[paste0("efflib_", m)]] <- eff_by_method[[m]] }
  nf_path <- file.path(opt$outdir, sprintf("%s_normfactors_used.csv", tag))
  write.csv(nf_df, nf_path, row.names = FALSE)
  say(sprintf("[NORMFACTORS] per-sample norm.factor + effective library size per method -> %s", basename(nf_path)))

  # ---- NEW: normalization-comparison figure (contrast only) ----
  if (opt$test == "contrast") {
    g1key <- attr(meta, "g1key"); g2key <- attr(meta, "g2key")
    ortho <- if (!is.null(opt$orthogonal_mark))
      mark_panel_data(opt$orthogonal_mark, g1key, g2key, methods, ss, opt, param, green_gr) else NULL
    comparison_scatter_fig(sig_by_method, ortho, methods, g1key, g2key, opt$group1, opt$group2,
                           opt$mark, opt$orthogonal_mark, opt$reference_shift, tag, opt$outdir)
    comparison_box_fig(sig_by_method, ortho, methods, g1key, g2key, opt$group1, opt$group2,
                       opt$mark, opt$orthogonal_mark, opt$reference_shift, tag, opt$outdir)
    say("[COMPARE] normalization-comparison scatter + box written.")
    if (is.null(opt$reference_shift))
      say("[COMPARE] NOTE: no --reference_shift given; panels show each method's ASSUMED global shift, not a validated magnitude. Anchor with the IF log2FC for interpretability.")
  } else {
    say("[COMPARE] normalization-comparison figure is contrast-only; skipped for interaction. Rerun --test contrast (e.g. --group1 'B10|UNC' --group2 'B10|none') to compare global-shift recovery.")
  }

  # ---- summary report ----
  summ <- file.path(opt$outdir, sprintf("%s_summary.txt", tag))
  hdr <- c(strrep("=", 72), sprintf("csaw v6 summary | %s", tag), strrep("=", 72), "",
           sprintf("Significant region counts per normalization (FDR<%.2g, |logFC|>%.2g):", opt$fdr, opt$lfc))
  rows <- vapply(methods, function(m) { cc <- counts_by_method[[m]]
    sprintf("  %-12s total=%-7d sig=%-6d up=%-6d down=%-6d", m, cc["total"], cc["sig"], cc["up"], cc["down"]) }, character(1))
  writeLines(c(hdr, rows, "", "--- full run log ---", log_lines$txt), summ)
  say(sprintf("\n[DONE] summary -> %s", summ))
}

main()