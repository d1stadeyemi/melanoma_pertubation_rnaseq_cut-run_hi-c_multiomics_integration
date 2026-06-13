# =============================================================================
# csaw.R — Window-based differential binding for BROAD histone marks
#          (v5: + fixed spike-in sample keying, MDS shape/colour redesign with a
#           batch-coloured companion, contrast-specific output tags, a legend
#           slot for the strata figure, and up/down counts on the volcano)
# =============================================================================
#
# PROJECT: A375 melanoma parental-vs-B10 constricted-migration multi-omics.
#
# WHAT THIS SCRIPT IS FOR ------------------------------------------------------
#   Broad chromatin marks (H3K27me3, H3K9me3, H3K4me1) form diffuse domains, not
#   sharp peaks. We TILE the genome into fixed windows and test each window for
#   differential coverage with the edgeR negative-binomial QL model; replicates
#   supply the within-group variance that turns a fold-change into an FDR-
#   controlled p-value. Adjacent windows are then merged into regions and their
#   p-values combined. (Punctate H3K27Ac is handled separately by diffbind.R.)
#
# TWO TEST MODES (--test) -----------------------------------------------------
#   contrast    (default) Two-group comparison:
#                 --group1 NUMERATOR --group2 REFERENCE  (Condition|Treatment)
#                 Migration:     --group1 "B10|none"     --group2 "Parental|none"
#                 UNC in B10:    --group1 "B10|UNC"      --group2 "B10|none"
#                 UNC in parent: --group1 "Parental|UNC" --group2 "Parental|none"
#   interaction  Factorial Condition*Treatment; tests the INTERACTION coefficient
#                 (does UNC's effect differ between B10 and Parental?). Requires
#                 all four cells (only H3K27me3 has them here).
#   The four "Condition|Treatment" cells: Parental|none, B10|none, Parental|UNC,
#   B10|UNC. Matching is CASE-INSENSITIVE ("Parental"=="parental").
#
# OUTPUT TAG (self-identifying figures/files) ---------------------------------
#   contrast    -> <mark>_<group1>_vs_<group2>   (e.g. H3K27me3_B10UNC_vs_B10none)
#   interaction -> <mark>_interaction
#   So two different contrasts no longer write identically-titled figures.
#
# KEYWORD 'sheet' RESOLUTION (peaks / igg / spikein) --------------------------
#   Each of --peaks, --igg, --spikein accepts either the literal word 'sheet'
#   (pull the information straight from the sample sheet, no retyping) OR an
#   explicit external file/list. 'sheet' requires the matching sheet column:
#     --peaks   sheet -> 'Peaks' column of the REFERENCE rows (group2 for a
#                        contrast; the B10|none cell for an interaction), unioned
#               file  -> comma-separated BED/narrowPeak/broadPeak/SEACR paths
#     --igg     sheet -> 'bamControl' column of the test samples
#               file  -> comma-separated IgG BAM paths
#     --spikein sheet -> 'spikein_reads' column of the test samples
#               file  -> CSV with columns SampleID,spikein_reads
#   --baseline_group is always sheet-driven: you name a 'Condition|Treatment'
#   group and its samples (already in the sheet) are counted in the same windows.
#
# BATCH HANDLING (per test, automatic) ----------------------------------------
#   Batch is a categorical label (e.g. b20250401 = library-prep day / ConA bead
#   lot / sequencing run). For EACH test the script inspects the Batch x Group
#   structure of the samples involved and reports one of:
#     [BATCH same]       one batch  -> no batch term (single level carries no
#                        information); spike-in allowed.
#     [BATCH crossed]    the tested coefficient is estimable AND the full design
#                        is full-rank with Batch in the model -> Batch IS added to
#                        the design; spike-in allowed (most reliable WITHIN a
#                        batch).
#     [BATCH partial]    the tested coefficient is estimable with Batch in the
#                        model, but the FULL design is rank-deficient because a
#                        NON-tested nuisance coefficient (e.g. the Treatment main
#                        effect in an interaction model) is aliased with Batch.
#                        edgeR needs the WHOLE design full-rank, so Batch is
#                        DROPPED: the tested effect is still estimable without it
#                        and Batch carried no usable information for it. The
#                        aliased main effect is NOT interpretable; only the tested
#                        (interaction) effect is reported. Spike-in DISABLED (the
#                        interaction is a spike-in-independent relative claim; a
#                        global-magnitude statement needs an orthogonal assay).
#     [BATCH CONFOUNDED] Batch aliases the TESTED effect -> a loud warning; the
#                        run CONTINUES as if a single batch (no batch term) so you
#                        still get a result, BUT spike-in is DISABLED for this test
#                        (cross-batch carryover cannot be separated from biology).
#   The router requires the WHOLE design to be full rank (edgeR's requirement),
#   not merely that the tested coefficient is estimable; both are checked by an
#   actual rank test, correct for the contrast and the interaction designs.
#
# BASELINE STRATIFICATION (the PRC2 "flattening" signature) -------------------
#   PRC2 inhibition erodes high-enrichment peaks while signal spreads into
#   low/background regions (Romero et al. 2024). A single genome-wide test lets
#   those opposite moves cancel, so we tag every region with its baseline
#   abundance, split into strata, and report the up/down balance per stratum.
#   Net LOSS at 'high' + net GAIN at 'low' = flattening.
#     samedata  (default, RTM-PRONE) ranks by the reference group whose noise also
#               enters logFC -> mechanically anti-correlated, can MANUFACTURE the
#               shape via regression-to-the-mean (RTM).
#     peak      (--peaks, RTM-IMMUNE, PRIMARY) peak membership is called vs IgG, a
#               different model, so it does not inherit the contrast's logFC noise,
#               and it marks where the mark IS in the treated cells.
#     indep     (--baseline_group, threshold-free, SECONDARY) an independent group
#               counted in the same windows, used only to rank baseline.
#   The run prints a [STRAT:samedata]/[STRAT:peak]/[STRAT:indep] block plus an
#   [STRAT verdict]: does the high-vs-low loss gradient SURVIVE the RTM-resistant
#   stratifier? If yes the flattening is real; if it collapses it was largely RTM.
#
# SPIKE-IN (E. coli carryover; RRPM/Orlando) ---------------------------------
#   Effective library size is set proportional to spike-in reads (norm.factors =
#   spikein_reads / lib.size, geometric-mean-centred), REPLACING TMM, and written
#   as a second 'spikein' result set. VALID ONLY across samples with comparable
#   cell number + ConA beads -> WITHIN BATCH (see batch handling above). Spike-in
#   fixes the GLOBAL scale (total loss, up/down balance), not which loci change
#   relative to each other; treat libsize-vs-spikein agreement as a robustness
#   check and anchor any global-magnitude claim with an orthogonal assay
#   (e.g. western blot for total H3K27me3).
#
# SAMPLE SHEET ----------------------------------------------------------------
#   Required : SampleID, Factor, Condition, Treatment, Replicate, Batch, bamReads
#   Optional : spikein_reads (for --spikein sheet), Peaks (for --peaks sheet),
#              bamControl (for --igg sheet)
#   Allowed  : Condition Parental|B10 ; Treatment none|UNC (case-insensitive)
#   Batch    : categorical, e.g. b20250401 (string so it is never read as numeric)
#
# EXAMPLE USAGE ---------------------------------------------------------------
#   # (1) Migration:
#   Rscript csaw.R --sheet ss.csv --mark H3K27me3 \
#       --group1 "B10|none" --group2 "Parental|none" --peaks sheet \
#       --outdir results/diff/H3K27me3_migration
#   # (2) UNC reversal in B10 (flagship, everything from the sheet):
#   Rscript csaw.R --sheet ss.csv --mark H3K27me3 \
#       --group1 "B10|UNC" --group2 "B10|none" \
#       --peaks sheet --baseline_group "Parental|none" --spikein sheet \
#       --outdir results/diff/H3K27me3_B10unc
#   # (3) Interaction:
#   Rscript csaw.R --sheet ss.csv --mark H3K27me3 --test interaction \
#       --peaks sheet --spikein sheet \
#       --outdir results/diff/H3K27me3_interaction
#   # (4) External files instead of the sheet:
#   Rscript csaw.R ... --peaks a.broadPeak,b.broadPeak --igg igg1.bam,igg2.bam \
#       --spikein factors.csv ...
#
# OUTPUTS (into --outdir) -----------------------------------------------------
#   <tag>_<norm>_all_regions.csv   coords, logFC, FDR, direction, baseline_logCPM,
#                                  baseline_stratum, baseline_indep_logCPM,
#                                  baseline_indep_stratum, in_peak
#   <tag>_<norm>_significant.csv   subset passing FDR & |log2FC|
#   <tag>_<norm>_{MDS,MDS_batch,MA,volcano,repcor,strata_balance}.{pdf,png}
#       MDS       : shape = Condition, colour = Treatment
#       MDS_batch : shape = Condition, colour = Batch (confound diagnostic)
#       volcano   : up/down significant counts annotated
#   <tag>_norm_compare.{pdf,png}   libsize-vs-spikein logFC (only with spike-in)
#   <tag>_summary.txt              consolidated report (counts, batch verdict,
#                                  stratification, spike-in factors, RTM verdict)
#   <norm> is 'libsize' always, plus 'spikein' when --spikein is given & allowed.
#   <tag>  is <mark>_<group1>_vs_<group2> (contrast) or <mark>_interaction.
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
                help = "REQUIRED. Sample sheet CSV. Columns: SampleID, Factor, Condition, Treatment, Replicate, Batch, bamReads (+ optional spikein_reads / Peaks / bamControl for the 'sheet' keyword)."),
    make_option("--mark", type = "character",
                help = "REQUIRED. Factor/antibody, matched case-insensitively to the Factor column. Examples: H3K27me3 | H3K9me3 | H3K4me1. (Punctate H3K27Ac -> diffbind.R.)"),

    make_option("--test", type = "character", default = "contrast",
                help = "Test mode. ALLOWED: 'contrast' | 'interaction'. [default %default]"),

    make_option("--group1", type = "character", default = "B10|none",
                help = "[contrast] Numerator 'Condition|Treatment'. One of: 'Parental|none' | 'B10|none' | 'Parental|UNC' | 'B10|UNC'. [default %default]"),
    make_option("--group2", type = "character", default = "Parental|none",
                help = "[contrast] Reference 'Condition|Treatment' (same allowed set). Positive logFC = higher in group1. [default %default]"),

    make_option("--cond_ref", type = "character", default = "Parental",
                help = "[interaction] Reference Condition. ALLOWED: 'Parental' | 'B10'. [default %default]"),
    make_option("--treat_ref", type = "character", default = "none",
                help = "[interaction] Reference Treatment. ALLOWED: 'none' | 'UNC'. [default %default]"),
    make_option("--cond_alt", type = "character", default = "B10",
                help = "[interaction] Non-reference Condition. ALLOWED: 'Parental' | 'B10'. [default %default]"),
    make_option("--treat_alt", type = "character", default = "UNC",
                help = "[interaction] Non-reference Treatment. ALLOWED: 'none' | 'UNC'. [default %default]"),

    make_option("--outdir", type = "character", default = "csaw_out",
                help = "Output directory, e.g. results/diff/H3K27me3_B10unc. [default %default]"),

    make_option("--window", type = "integer", default = 2000L,
                help = "Window width (bp). Broad marks: 2000. [default %default]"),
    make_option("--spacing", type = "integer", default = 500L,
                help = "Spacing between window starts (bp); also the region-merge tolerance. [default %default]"),
    make_option("--bg_bin", type = "integer", default = 10000L,
                help = "Background bin width (bp) for TMM + the enrichment filter. [default %default]"),
    make_option("--frag_max", type = "integer", default = 700L,
                help = "Max fragment length (bp); match Bowtie2 -X. [default %default]"),

    make_option("--min_lfc_bg", type = "double", default = log2(3),
                help = "Min log2 enrichment over background/IgG to KEEP a window. [default log2(3)=1.58]"),
    make_option("--fdr", type = "double", default = 0.05,
                help = "FDR cutoff for the significant subset. [default %default]"),
    make_option("--lfc", type = "double", default = 1.0,
                help = "Min |log2FC| for the significant subset. [default %default]"),
    make_option("--n_strata", type = "integer", default = 3L,
                help = "Number of baseline-abundance strata (>=2). 3 -> low/mid/high. [default %default]"),

    make_option("--blacklist", type = "character", default = NULL,
                help = "Optional ENCODE blacklist BED (e.g. hg19-blacklist.v2.bed). [default: none]"),
    make_option("--igg", type = "character", default = NULL,
                help = "IgG enrichment filter. 'sheet' (read bamControl column) OR comma-separated IgG BAM path(s). [default: none]"),
    make_option("--spikein", type = "character", default = NULL,
                help = "Spike-in (Orlando/RRPM) normalization from READ COUNTS: norm.factor = spikein_reads/lib.size. 'sheet' (read spikein_reads column) OR a CSV with columns SampleID,spikein_reads. Writes a 'spikein' result set. [default: none]"),
    make_option("--normfactors", type = "character", default = NULL,
                help = "PRECOMPUTED scaling factors injected directly as norm.factors (e.g. ChIPseqSpikeInFree). Path to a tab/comma table with an ID column (SampleID or ID, matched to SampleID or bamReads basename) and a factor column (SF/factor/ScalingFactor/norm.factors). A 'QC' column (ChIPseqSpikeInFree) is checked; 'NA'/non-positive SF is rejected. Writes an 'sf' result set. [default: none]"),
    make_option("--force_norm", action = "store_true", default = FALSE,
                help = "Apply --spikein / --normfactors EVEN when the batch verdict is non-separable (confounded/partial). Use only when the global magnitude is anchored externally (e.g. ChIPseqSpikeInFree + western blot). [default: FALSE]"),

    make_option("--peaks", type = "character", default = NULL,
                help = "RTM-immune peak stratifier. 'sheet' (read Peaks column of the reference rows, unioned) OR comma-separated BED/narrowPeak/broadPeak/SEACR path(s). Use a BROAD call for H3K27me3. [default: none]"),
    make_option("--baseline_group", type = "character", default = NULL,
                help = "Independent baseline stratifier: a 'Condition|Treatment' group (e.g. 'Parental|none') counted in the same windows and used ONLY to rank baseline. Must be independent of the tested coefficient. [default: none]")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$sheet)) stop("--sheet is required.")
  if (is.null(opt$mark))  stop("--mark is required (e.g. H3K27me3).")
  if (!opt$test %in% c("contrast", "interaction"))
    stop("--test must be 'contrast' or 'interaction'.")
  if (opt$n_strata < 2L) stop("--n_strata must be >= 2.")
  opt
}

parse_group <- function(spec) {
  parts <- strsplit(spec, "\\|")[[1]]
  if (length(parts) != 2)
    stop("Group spec must be 'Condition|Treatment', got: ", spec)
  list(cond = tolower(trimws(parts[1])), treat = tolower(trimws(parts[2])))
}

# filesystem-safe token (strip '|', spaces, etc.) for the output tag
safe_token <- function(s) gsub("[^A-Za-z0-9]+", "", s)

# =============================================================================
# Small helpers
# =============================================================================
log_lines <- new.env()                 # collect lines for the summary report
log_lines$txt <- character(0)
say <- function(...) {
  line <- paste0(...)
  cat(line, "\n")
  log_lines$txt <- c(log_lines$txt, line)
}

read_bed3 <- function(path) {
  # Robust 3-column reader for BED / narrowPeak / broadPeak / SEACR-bed.
  df <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                   comment.char = "#", quote = "")
  GRanges(seqnames = df[[1]],
          ranges = IRanges(start = as.integer(df[[2]]) + 1L, end = as.integer(df[[3]])))
}

# =============================================================================
# Sample selection
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
    # Group factor: reference = group2
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
    if (any(is.na(sub$Condition)) || any(is.na(sub$Treatment)))
      stop("Interaction factor level mismatch.")
    # reference rows for sheet-peaks: the cond_alt|treat_ref cell (e.g. B10|none)
    attr(sub, "peak_ref") <- sub$.cond == ca & sub$.treat == tr
  }
  if (any(!file.exists(sub$bamReads)))
    stop("Missing BAM(s): ", paste(sub$bamReads[!file.exists(sub$bamReads)], collapse = ", "))
  sub
}

# Which rows are the "reference" for sheet-peaks?
ref_rows <- function(meta, opt) {
  if (opt$test == "contrast") meta$.role == "group2" else attr(meta, "peak_ref")
}

# =============================================================================
# Resolve peaks / igg / spikein from 'sheet' or external files
# =============================================================================
resolve_peaks <- function(opt, meta) {
  if (is.null(opt$peaks)) return(NULL)
  if (identical(opt$peaks, "sheet")) {
    if (!"Peaks" %in% colnames(meta))
      stop("--peaks sheet requested but the sample sheet has no 'Peaks' column.")
    paths <- unique(meta$Peaks[ref_rows(meta, opt)])
    paths <- paths[!is.na(paths) & nzchar(paths)]
    if (!length(paths)) stop("--peaks sheet: no Peaks paths on the reference rows.")
  } else {
    paths <- trimws(strsplit(opt$peaks, ",")[[1]])
  }
  miss <- paths[!file.exists(paths)]
  if (length(miss)) stop("Peak file(s) not found: ", paste(miss, collapse = ", "))
  gr <- Reduce(c, lapply(paths, read_bed3))
  gr <- reduce(gr)
  say(sprintf("[PEAKS] %d peak file(s) -> %d merged peak regions", length(paths), length(gr)))
  gr
}

resolve_igg <- function(opt, meta) {
  if (is.null(opt$igg)) return(NULL)
  if (identical(opt$igg, "sheet")) {
    if (!"bamControl" %in% colnames(meta))
      stop("--igg sheet requested but the sample sheet has no 'bamControl' column.")
    bams <- unique(meta$bamControl)
    bams <- bams[!is.na(bams) & nzchar(bams)]
  } else {
    bams <- trimws(strsplit(opt$igg, ",")[[1]])
  }
  miss <- bams[!file.exists(bams)]
  if (length(miss)) stop("IgG BAM(s) not found: ", paste(miss, collapse = ", "))
  bams
}

resolve_spikein <- function(opt, meta) {
  # returns a named numeric vector spikein_reads keyed by SampleID, or NULL
  if (is.null(opt$spikein)) return(NULL)
  if (identical(opt$spikein, "sheet")) {
    if (!"spikein_reads" %in% colnames(meta))
      stop("--spikein sheet requested but the sample sheet has no 'spikein_reads' column.")
    v <- setNames(as.numeric(meta$spikein_reads), meta$SampleID)
  } else {
    sp <- read.csv(opt$spikein, stringsAsFactors = FALSE)
    if (!all(c("SampleID", "spikein_reads") %in% colnames(sp)))
      stop("--spikein CSV needs columns SampleID,spikein_reads.")
    v <- setNames(as.numeric(sp$spikein_reads), sp$SampleID)
    v <- v[meta$SampleID]
  }
  if (any(is.na(v)) || any(v <= 0))
    stop("spikein_reads has NA / non-positive values for: ",
         paste(meta$SampleID[is.na(v) | v <= 0], collapse = ", "))
  v
}

resolve_normfactors <- function(opt, meta) {
  # Precomputed per-sample scaling factors (e.g. ChIPseqSpikeInFree _SF.txt or any
  # CSV). Flexible: auto-detect tab/comma; accept SF/factor/ScalingFactor/norm.factors
  # as the factor column and SampleID/ID/Sample as the key; match the key against
  # SampleID first, then the bamReads basename (ChIPseqSpikeInFree keys by bam name).
  # Returns a numeric vector keyed by meta$SampleID, or NULL.
  if (is.null(opt$normfactors)) return(NULL)
  if (!file.exists(opt$normfactors)) stop("--normfactors file not found: ", opt$normfactors)
  first <- readLines(opt$normfactors, n = 1L, warn = FALSE)
  sep   <- if (grepl("\t", first)) "\t" else ","
  tab   <- read.table(opt$normfactors, header = TRUE, sep = sep, quote = "\"",
                      stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE)

  fcol <- intersect(c("SF", "factor", "ScalingFactor", "norm.factors", "normfactor"),
                    colnames(tab))
  if (!length(fcol))
    stop("--normfactors: no scaling-factor column (looked for SF / factor / ScalingFactor / ",
         "norm.factors). Columns present: ", paste(colnames(tab), collapse = ", "))
  icol <- intersect(c("SampleID", "ID", "Sample", "sample"), colnames(tab))
  if (!length(icol))
    stop("--normfactors: no ID column (looked for SampleID / ID / Sample). Columns present: ",
         paste(colnames(tab), collapse = ", "))
  ids <- as.character(tab[[icol[1]]])
  sf  <- suppressWarnings(as.numeric(tab[[fcol[1]]]))

  # match each of our samples to a row in the SF table
  bn  <- basename(meta$bamReads)
  out <- setNames(rep(NA_real_, nrow(meta)), meta$SampleID)
  qc_used <- character(0)
  for (i in seq_len(nrow(meta))) {
    hit <- which(ids == meta$SampleID[i])
    if (!length(hit)) hit <- which(ids == bn[i])
    if (!length(hit)) hit <- which(basename(ids) == bn[i])
    if (length(hit)) {
      out[i] <- sf[hit[1]]
      if ("QC" %in% colnames(tab)) qc_used <- c(qc_used, as.character(tab$QC[hit[1]]))
    }
  }
  if (length(qc_used) && any(!grepl("pass", tolower(qc_used))))
    say("[NORMFACTORS] WARNING: a matched row has QC != 'pass' (ChIPseqSpikeInFree flags poor ",
        "enrichment / complete loss). Verify before trusting its factor.")
  if (any(is.na(out)))
    stop("--normfactors: could not match a factor for: ",
         paste(meta$SampleID[is.na(out)], collapse = ", "),
         "\n  ID column held (first few): ", paste(head(ids, 8), collapse = ", "),
         "\n  Matching is tried against SampleID and the bamReads basename.")
  if (any(out <= 0))
    stop("--normfactors: non-positive / NA factor for: ",
         paste(meta$SampleID[out <= 0], collapse = ", "),
         "\n  (A ChIPseqSpikeInFree 'NA' SF means complete loss / poor enrichment; resolve per its ",
         "README before injecting.)")
  out
}

# =============================================================================
# Batch assessment (per test): same / crossed / partial / confounded
#   The WHOLE design must be full rank (edgeR's requirement), not just the tested
#   coefficient; both are checked by a rank test.
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
    d_bio   <- model.matrix(~ Group, meta)
    tested  <- setdiff(colnames(d_bio), "(Intercept)")
    f_batch <- ~ Batch + Group
  } else {
    d_bio   <- model.matrix(~ Condition * Treatment, meta)
    tested  <- grep(":", colnames(d_bio), value = TRUE)
    f_batch <- ~ Batch + Condition * Treatment
  }
  nb <- nlevels(droplevels(factor(meta$Batch)))
  if (nb <= 1L)
    return(list(status = "same", design = d_bio, tested = tested, allow_spikein = TRUE))

  d_full <- model.matrix(f_batch, meta)
  tested_full <- if (opt$test == "contrast")
    setdiff(colnames(d_full), c("(Intercept)", grep("^Batch", colnames(d_full), value = TRUE)))
  else grep(":", colnames(d_full), value = TRUE)

  # edgeR needs the WHOLE design full-rank, not just the tested coefficient.
  full_rank <- qr(d_full)$rank == ncol(d_full)
  tested_ok <- length(tested_full) > 0 &&
               all(vapply(tested_full, function(cf) coef_estimable(d_full, cf), logical(1)))

  if (full_rank && tested_ok) {
    # Batch is genuinely blockable (groups present across batches).
    list(status = "crossed", design = d_full, tested = tested_full, allow_spikein = TRUE)
  } else if (tested_ok) {
    # Tested coef estimable WITH Batch, but a NON-tested nuisance coef (e.g. the
    # Treatment main effect) is aliased with Batch -> full design rank-deficient.
    # Drop Batch: the tested (interaction) effect is still estimable without it,
    # and Batch held no information the fit could use for it. The aliased main
    # effect is NOT interpretable; only the interaction is reported.
    list(status = "partial", design = d_bio, tested = tested, allow_spikein = FALSE)
  } else {
    # The tested coefficient itself is aliased with Batch -> cannot be separated.
    list(status = "confounded", design = d_bio, tested = tested, allow_spikein = FALSE)
  }
}

# =============================================================================
# Counting + filtering
# =============================================================================
make_param <- function(opt) {
  discard <- GRanges()                       # empty (readParam's own default); NOT NULL
  if (!is.null(opt$blacklist)) {
    if (!file.exists(opt$blacklist)) stop("Blacklist not found: ", opt$blacklist)
    discard <- read_bed3(opt$blacklist)
  }
  readParam(minq = 20L, dedup = TRUE, pe = "both",
            max.frag = opt$frag_max, discard = discard)
}

count_and_filter <- function(bams, opt, param, igg_bams) {
  win  <- windowCounts(bams, width = opt$window, spacing = opt$spacing, param = param)
  bins <- windowCounts(bams, bin = TRUE, width = opt$bg_bin, param = param)
  if (!is.null(igg_bams)) {
    control     <- regionCounts(igg_bams, regions = rowRanges(win), param = param)
    bin.control <- windowCounts(igg_bams, bin = TRUE, width = opt$bg_bin, param = param)
    scinfo <- scaleControlFilter(bins, bin.control)
    fstat  <- filterWindowsControl(win, control, scale.info = scinfo)
    say("[FILTER] IgG-enrichment filter (filterWindowsControl).")
  } else {
    fstat <- filterWindowsGlobal(win, bins)
    say("[FILTER] genomic-background enrichment filter (filterWindowsGlobal).")
  }
  keep <- fstat$filter > opt$min_lfc_bg
  say(sprintf("[FILTER] kept %d / %d windows (>%.2f log2 over background).",
              sum(keep), length(keep), opt$min_lfc_bg))
  list(win = win[keep, ], bins = bins)
}

# =============================================================================
# Normalization (TMM = libsize, and spike-in = RRPM)
# =============================================================================
norm_tmm <- function(win, bins) {
  win <- normFactors(bins, se.out = win)            # composition/efficiency TMM
  asDGEList(win)
}
set_spikein <- function(y, spikein_reads) {
  sr <- spikein_reads[colnames(y)]                  # keyed by SampleID (= colnames(y))
  bad <- is.na(sr) | !is.finite(sr) | sr <= 0
  if (any(bad))
    stop("[SPIKEIN] spike-in reads missing/invalid for: ",
         paste(colnames(y)[bad], collapse = ", "),
         "\n  Check that sample-sheet SampleID values match the BAM columns and that ",
         "spikein_reads are positive finite numbers.")
  nf <- sr / y$samples$lib.size
  nf <- nf / exp(mean(log(nf)))                     # geometric-mean-centre to 1
  y$samples$norm.factors <- as.numeric(nf)
  y
}
set_normfactors <- function(y, sf) {
  # Inject precomputed scaling factors (e.g. ChIPseqSpikeInFree SF) as edgeR
  # norm.factors. Geometric-mean-centred for numerical hygiene; DE significance and
  # the up/down balance are invariant to a global rescale, so the centring does not
  # change the result — only the absolute CPM baseline.
  v <- sf[colnames(y)]
  bad <- is.na(v) | !is.finite(v) | v <= 0
  if (any(bad))
    stop("[NORMFACTORS] invalid factor for: ", paste(colnames(y)[bad], collapse = ", "),
         "\n  Check the --normfactors ID column maps to these SampleIDs.")
  v <- v / exp(mean(log(v)))
  y$samples$norm.factors <- as.numeric(v)
  y
}

# =============================================================================
# Fit + test + region-merge, then attach baseline / peak strata
# =============================================================================
fit_test_regions <- function(y, ab, opt, win_gr, peaks_gr, indep_logcpm) {
  y   <- estimateDisp(y, ab$design)
  fit <- glmQLFit(y, ab$design, robust = TRUE)
  res <- glmQLFTest(fit, coef = ab$tested)          # Group coef OR interaction coef
  wtab <- res$table                                 # per-window logFC, logCPM, PValue

  # reference-group baseline abundance per window (RTM-prone 'samedata')
  if (opt$test == "contrast") {
    ref_idx <- which(y$samples$group_role == "group2")
  } else {
    ref_idx <- which(y$samples$peak_ref)
  }
  base_lcpm <- aveLogCPM(y$counts[, ref_idx, drop = FALSE],
                         lib.size = y$samples$lib.size[ref_idx] *
                                    y$samples$norm.factors[ref_idx])

  # merge adjacent windows into regions, combine p-values
  merged <- mergeWindows(win_gr, tol = opt$spacing, max.width = 10000L)
  comb   <- combineTests(merged$ids, wtab)          # PValue, FDR, direction, rep.test
  reg    <- merged$regions

  # aggregate per-window stats to regions (mean over member windows)
  agg <- function(v) as.numeric(tapply(v, merged$ids, mean))
  reg_logFC    <- agg(wtab$logFC)
  reg_baseline <- agg(base_lcpm)
  reg_indep    <- if (!is.null(indep_logcpm)) agg(indep_logcpm) else rep(NA_real_, length(reg))

  df <- data.frame(
    chr = as.character(seqnames(reg)), start = start(reg), end = end(reg),
    logFC = reg_logFC, logCPM = agg(wtab$logCPM),
    PValue = comb$PValue, FDR = comb$FDR,
    direction = comb$direction,
    baseline_logCPM = reg_baseline,
    baseline_indep_logCPM = reg_indep,
    stringsAsFactors = FALSE
  )
  # strata (low..high) from baseline ranks
  cut_strata <- function(x, k) {
    if (all(is.na(x))) return(factor(rep(NA, length(x))))
    br <- quantile(x, probs = seq(0, 1, length.out = k + 1), na.rm = TRUE)
    lab <- if (k == 3) c("low", "mid", "high") else paste0("s", seq_len(k))
    cut(x, breaks = unique(br), labels = lab[seq_len(length(unique(br)) - 1)],
        include.lowest = TRUE)
  }
  df$baseline_stratum       <- cut_strata(df$baseline_logCPM, opt$n_strata)
  df$baseline_indep_stratum <- cut_strata(df$baseline_indep_logCPM, opt$n_strata)
  # in_peak via overlap
  if (!is.null(peaks_gr)) {
    df$in_peak <- overlapsAny(reg, peaks_gr)
  } else {
    df$in_peak <- NA
  }
  list(table = df, res = res, fit = fit, y = y)
}

# =============================================================================
# Stratified up/down balance + RTM verdict
# =============================================================================
balance_block <- function(df, stratum_col, opt, label) {
  sig <- df$FDR < opt$fdr & abs(df$logFC) > opt$lfc
  s <- df[[stratum_col]]
  levs <- levels(factor(s))
  say(sprintf("[STRAT:%s]", label))
  out <- list()
  for (L in levs) {
    m  <- !is.na(s) & s == L
    ms <- m & sig
    n_up   <- sum(ms & df$logFC > 0)
    n_dn   <- sum(ms & df$logFC < 0)
    n_sig  <- n_up + n_dn
    pct_dn <- if (n_sig) 100 * n_dn / n_sig else NA_real_
    say(sprintf("   %-5s n=%-6d sig=%-5d  up=%-5d down=%-5d  %%down(sig)=%s  medianLFC=%+.3f",
                L, sum(m), n_sig, n_up, n_dn,
                ifelse(is.na(pct_dn), "NA", sprintf("%.0f", pct_dn)),
                median(df$logFC[m], na.rm = TRUE)))
    out[[L]] <- c(n = sum(m), sig = n_sig, up = n_up, down = n_dn, pct_down = pct_dn)
  }
  out
}

balance_peak <- function(df, opt) {
  if (all(is.na(df$in_peak))) { say("[STRAT:peak]   (no --peaks given)"); return(NULL) }
  sig <- df$FDR < opt$fdr & abs(df$logFC) > opt$lfc
  say("[STRAT:peak]")
  out <- list()
  for (L in c(TRUE, FALSE)) {
    m  <- df$in_peak == L; ms <- m & sig
    n_up <- sum(ms & df$logFC > 0); n_dn <- sum(ms & df$logFC < 0); n_sig <- n_up + n_dn
    pct_dn <- if (n_sig) 100 * n_dn / n_sig else NA_real_
    say(sprintf("   %-8s n=%-6d sig=%-5d  up=%-5d down=%-5d  %%down(sig)=%s  medianLFC=%+.3f",
                ifelse(L, "in_peak", "non_peak"), sum(m), n_sig, n_up, n_dn,
                ifelse(is.na(pct_dn), "NA", sprintf("%.0f", pct_dn)),
                median(df$logFC[m], na.rm = TRUE)))
    out[[ifelse(L, "peak", "nonpeak")]] <- c(pct_down = pct_dn)
  }
  out
}

rtm_verdict <- function(sd_block, peak_block, indep_block) {
  # loss gradient (samedata) = %down(high) - %down(low)
  g <- function(b, hi, lo) {
    if (is.null(b) || is.null(b[[hi]]) || is.null(b[[lo]])) return(NA_real_)
    b[[hi]]["pct_down"] - b[[lo]]["pct_down"]
  }
  grad_sd <- g(sd_block, "high", "low")
  grad_pk <- if (!is.null(peak_block))
    (peak_block[["peak"]]["pct_down"] - peak_block[["nonpeak"]]["pct_down"]) else NA_real_
  grad_in <- g(indep_block, "high", "low")
  say("[STRAT verdict]")
  say(sprintf("   samedata loss-gradient (high-low %%down) = %s",
              ifelse(is.na(grad_sd), "NA", sprintf("%+.0f", grad_sd))))
  if (!is.na(grad_pk))
    say(sprintf("   peak     loss-gradient (peak-nonpeak)   = %+.0f  -> %s",
                grad_pk, ifelse(!is.na(grad_sd) && grad_pk >= 0.5 * grad_sd,
                                "SURVIVES (flattening real)", "COLLAPSES (largely RTM)")))
  if (!is.na(grad_in))
    say(sprintf("   indep    loss-gradient (high-low)       = %+.0f  -> %s",
                grad_in, ifelse(!is.na(grad_sd) && grad_in >= 0.5 * grad_sd,
                                "SURVIVES", "weak/RTM")))
  if (is.na(grad_pk) && is.na(grad_in))
    say("   (no RTM-resistant stratifier supplied; samedata is RTM-prone — interpret with care)")
}

# =============================================================================
# Figures
# =============================================================================
save_fig <- function(fn, outdir, draw, w = 6.5, h = 5) {
  for (dev in c("pdf", "png")) {
    if (dev == "pdf") pdf(file.path(outdir, paste0(fn, ".pdf")), width = w, height = h)
    else png(file.path(outdir, paste0(fn, ".png")), width = w, height = h, units = "in", res = 200)
    draw(); dev.off()
  }
}

# MDS with biological glyphs: shape = Condition, colour = Treatment (default) or
# Batch (the confound diagnostic). Real SampleID labels; padded so nothing clips.
mds_plot <- function(y, tag, norm, outdir, color_by = c("treat", "batch")) {
  color_by <- match.arg(color_by)
  lc  <- cpm(y, log = TRUE)
  mds <- suppressWarnings(plotMDS(lc, plot = FALSE))
  xx  <- as.numeric(mds$x); yy <- as.numeric(mds$y)
  cond  <- as.character(y$samples$cond)
  treat <- as.character(y$samples$treat)
  batch <- as.character(y$samples$batch)
  ids   <- colnames(y)

  conds   <- sort(unique(cond))
  pch_map <- setNames(c(21, 24, 22, 23, 25)[seq_along(conds)], conds)
  pchs    <- pch_map[cond]

  keyvec  <- if (color_by == "treat") treat else batch
  keys    <- sort(unique(keyvec))
  pal     <- c("#1F6F8B", "#E08A3C", "#2A9D8F", "#B23A48", "#6D5BA6", "#9BA7B0")[seq_along(keys)]
  col_map <- setNames(pal, keys)
  cols    <- col_map[keyvec]

  pct <- function(i) {
    ve <- mds$var.explained
    if (!is.null(ve) && length(ve) >= i && is.finite(ve[i]))
      sprintf(" (%d%%)", round(100 * ve[i])) else ""
  }
  xr <- range(xx); yr <- range(yy)
  xpad <- max(diff(xr) * 0.30, 0.5); ypad <- max(diff(yr) * 0.30, 0.5)
  ttl <- sprintf("%s %s MDS  (shape=Condition, colour=%s)",
                 tag, norm, if (color_by == "batch") "Batch" else "Treatment")
  fn  <- paste0(tag, "_", norm, "_MDS", if (color_by == "batch") "_batch" else "")

  save_fig(fn, outdir, function() {
    par(mar = c(5, 5, 4, 9), xpd = TRUE)
    plot(xx, yy, pch = pchs, bg = cols, col = "grey20", cex = 2.2, lwd = 0.9,
         xlim = c(xr[1] - xpad, xr[2] + xpad),
         ylim = c(yr[1] - ypad, yr[2] + ypad),
         xlab = paste0("Leading logFC dim 1", pct(1)),
         ylab = paste0("Leading logFC dim 2", pct(2)),
         main = ttl)
    text(xx, yy, labels = ids, pos = 3, offset = 0.7, cex = 0.6, col = "grey30")
    legend("topright", inset = c(-0.30, 0), bty = "n", title = "Condition",
           pch = unname(pch_map[conds]), pt.bg = "grey80", col = "grey20",
           legend = conds, pt.cex = 1.5)
    legend("bottomright", inset = c(-0.30, 0), bty = "n",
           title = if (color_by == "batch") "Batch" else "Treatment",
           pch = 22, pt.bg = unname(col_map[keys]), col = "grey20",
           legend = keys, pt.cex = 1.5)
  }, w = 8, h = 5.5)
}

std_plots <- function(y, res, df, tag, norm, outdir, opt) {
  # MDS: shape = Condition, colour = Treatment, plus a Batch-coloured companion.
  mds_plot(y, tag, norm, outdir, "treat")
  mds_plot(y, tag, norm, outdir, "batch")

  save_fig(paste0(tag, "_", norm, "_MA"), outdir, function() {
    plot(df$logCPM, df$logFC, pch = 16, cex = 0.3, col = "#88888855",
         xlab = "logCPM", ylab = "logFC", main = paste(tag, norm, "MA"))
    sig <- df$FDR < opt$fdr & abs(df$logFC) > opt$lfc
    points(df$logCPM[sig], df$logFC[sig], pch = 16, cex = 0.35, col = "#C0392B")
    abline(h = 0, col = "grey40")
  })

  save_fig(paste0(tag, "_", norm, "_volcano"), outdir, function() {
    sig <- df$FDR < opt$fdr & abs(df$logFC) > opt$lfc
    plot(df$logFC, -log10(df$FDR), pch = 16, cex = 0.3, col = "#88888855",
         xlab = "logFC", ylab = "-log10 FDR", main = paste(tag, norm, "volcano"))
    points(df$logFC[sig], -log10(df$FDR[sig]), pch = 16, cex = 0.35,
           col = ifelse(df$logFC[sig] > 0, "#C0392B", "#1F6F8B"))
    abline(v = c(-opt$lfc, opt$lfc), h = -log10(opt$fdr), col = "grey60", lty = 2)
    n_up <- sum(sig & df$logFC > 0); n_dn <- sum(sig & df$logFC < 0)
    legend("top", horiz = TRUE, bg = "white", box.col = "grey80", cex = 0.95,
           legend = c(sprintf("down: %d", n_dn), sprintf("up: %d", n_up)),
           text.col = c("#1F6F8B", "#C0392B"))
  })

  save_fig(paste0(tag, "_", norm, "_repcor"), outdir, function() {
    lc <- cpm(y, log = TRUE)
    suppressWarnings(pairs(lc[sample(nrow(lc), min(2000, nrow(lc))), , drop = FALSE],
                           pch = ".", main = paste(tag, norm, "rep cor")))
  })
}

# Strata up/down balance. Only POPULATED panels are drawn, plus a slim dedicated
# legend cell at the far right (so the legend never overlaps the bars, and for the
# interaction — which has no indep panel — it sits right beside the peak panel).
strata_balance_fig <- function(df, tag, norm, outdir, opt) {
  sig <- df$FDR < opt$fdr & abs(df$logFC) > opt$lfc
  has_peak  <- !all(is.na(df$in_peak))
  has_indep <- !all(is.na(df$baseline_indep_logCPM))

  draw_stratum <- function(strat, title) {
    levs <- levels(factor(df[[strat]]))
    if (!length(levs)) { plot.new(); title(main = title); return() }
    up <- dn <- numeric(length(levs))
    for (i in seq_along(levs)) {
      m <- !is.na(df[[strat]]) & df[[strat]] == levs[i] & sig
      up[i] <- sum(m & df$logFC > 0); dn[i] <- sum(m & df$logFC < 0)
    }
    tot <- pmax(up + dn, 1)
    barplot(rbind(100 * up / tot, 100 * dn / tot), names.arg = levs,
            col = c("#E08A3C", "#2A9D8F"), ylab = "% of significant",
            main = title, ylim = c(0, 100))
  }
  draw_peak <- function() {
    lv <- c("non_peak", "in_peak"); up <- dn <- numeric(2)
    for (i in 1:2) {
      m <- df$in_peak == (lv[i] == "in_peak") & sig
      up[i] <- sum(m & df$logFC > 0); dn[i] <- sum(m & df$logFC < 0)
    }
    tot <- pmax(up + dn, 1)
    barplot(rbind(100 * up / tot, 100 * dn / tot), names.arg = lv,
            col = c("#E08A3C", "#2A9D8F"), ylab = "% of significant",
            main = "peak (RTM-immune)", ylim = c(0, 100))
  }

  panels <- list(list(kind = "stratum", col = "baseline_stratum",
                      title = "samedata (RTM-prone)"))
  if (has_peak)  panels <- c(panels, list(list(kind = "peak")))
  if (has_indep) panels <- c(panels, list(list(kind = "stratum",
                      col = "baseline_indep_stratum", title = "indep group")))
  np <- length(panels)

  save_fig(paste0(tag, "_", norm, "_strata_balance"), outdir, function() {
    layout(matrix(1:(np + 1), nrow = 1), widths = c(rep(1, np), 0.55))
    par(mar = c(5, 4, 4, 1))
    for (p in panels) {
      if (identical(p$kind, "peak")) draw_peak() else draw_stratum(p$col, p$title)
    }
    par(mar = c(5, 0, 4, 0)); plot.new()
    legend("left", fill = c("#E08A3C", "#2A9D8F"),
           legend = c("gained", "lost"), bty = "n", cex = 1.15)
  }, w = 2.7 * (np + 0.6), h = 4)
}

norm_compare_fig <- function(df_lib, df_other, tag, norm, outdir, opt) {
  key <- function(d) paste(d$chr, d$start, d$end)
  m <- merge(data.frame(k = key(df_lib), lf_lib = df_lib$logFC),
             data.frame(k = key(df_other), lf_oth = df_other$logFC), by = "k")
  save_fig(paste0(tag, "_", norm, "_norm_compare"), outdir, function() {
    plot(m$lf_lib, m$lf_oth, pch = 16, cex = 0.3, col = "#1F6F8B55",
         xlab = "logFC (libsize/TMM)", ylab = sprintf("logFC (%s)", norm),
         main = paste(tag, sprintf("libsize vs %s (per region)", norm)))
    abline(0, 1, col = "grey40"); abline(h = 0, v = 0, col = "grey80", lty = 3)
    off <- median(m$lf_oth - m$lf_lib, na.rm = TRUE)
    legend("topleft", bty = "n",
           legend = sprintf("median global shift = %+.2f log2", off))
  })
}

# =============================================================================
# Per-normalization driver
# =============================================================================
write_outputs <- function(df, tag, norm, outdir, opt) {
  write.csv(df, file.path(outdir, sprintf("%s_%s_all_regions.csv", tag, norm)), row.names = FALSE)
  sig <- df[df$FDR < opt$fdr & abs(df$logFC) > opt$lfc, , drop = FALSE]
  write.csv(sig, file.path(outdir, sprintf("%s_%s_significant.csv", tag, norm)), row.names = FALSE)
  n_up <- sum(sig$logFC > 0); n_dn <- sum(sig$logFC < 0)
  say(sprintf("[%s] significant regions: %d  (up=%d, down=%d)  of %d tested",
              norm, nrow(sig), n_up, n_dn, nrow(df)))
  invisible(c(n_sig = nrow(sig), up = n_up, down = n_dn, total = nrow(df)))
}

run_norm <- function(y, ab, opt, win_gr, peaks_gr, indep_logcpm, tag, norm, outdir) {
  say(sprintf("=== normalization: %s ===", norm))
  ft <- fit_test_regions(y, ab, opt, win_gr, peaks_gr, indep_logcpm)
  df <- ft$table
  counts <- write_outputs(df, tag, norm, outdir, opt)
  std_plots(ft$y, ft$res, df, tag, norm, outdir, opt)
  strata_balance_fig(df, tag, norm, outdir, opt)
  sd_block <- balance_block(df, "baseline_stratum", opt, "samedata")
  pk_block <- balance_peak(df, opt)
  in_block <- if (!all(is.na(df$baseline_indep_logCPM)))
    balance_block(df, "baseline_indep_stratum", opt, "indep") else { say("[STRAT:indep]   (no --baseline_group given)"); NULL }
  rtm_verdict(sd_block, pk_block, in_block)
  list(df = df, counts = counts)
}

# =============================================================================
# Main
# =============================================================================
main <- function() {
  opt <- parse_cli()
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

  # self-identifying tag: contrast encodes the actual comparison
  tag <- if (opt$test == "contrast")
    sprintf("%s_%s_vs_%s", opt$mark, safe_token(opt$group1), safe_token(opt$group2))
  else
    paste0(opt$mark, "_interaction")

  ss   <- load_sheet(opt$sheet)
  meta <- select_samples(ss, opt)
  say(sprintf("[INFO] mark=%s test=%s tag=%s  %d samples:", opt$mark, opt$test, tag, nrow(meta)))
  for (i in seq_len(nrow(meta)))
    say(sprintf("   %-22s  %s|%s  batch=%s", meta$SampleID[i], meta$Condition[i],
                meta$Treatment[i], meta$Batch[i]))

  # ---- batch assessment (per test) ----
  ab <- assess_batch(meta, opt)
  if (ab$status == "same") {
    say("[BATCH same] all test samples share one batch -> no batch term; spike-in allowed.")
  } else if (ab$status == "crossed") {
    say("[BATCH crossed] tested effect is estimable with a full-rank Batch design -> Batch ADDED.")
    say("               (note: spike-in is most reliable WITHIN a single batch.)")
  } else if (ab$status == "partial") {
    say("[BATCH partial] *** Batch aliases a NON-tested coefficient (here the Treatment main effect). ***")
    say("               The full ~Batch+Condition*Treatment design is rank-deficient, but the tested")
    say("               INTERACTION is estimable. Batch is DROPPED (it carried no usable information for")
    say("               the interaction); fitting ~Condition*Treatment. The Treatment main effect is")
    say("               CONFOUNDED with batch and is NOT interpretable -- only the interaction is reported.")
    say("               (Within-treatment batch structure, e.g. a lone extra-batch sample, is unmodelled.)")
    say("               Spike-in DISABLED: the interaction is a spike-in-independent relative claim;")
    say("               any global-magnitude statement needs an orthogonal assay (western/dot blot).")
  } else {
    say("[BATCH CONFOUNDED] *** Batch aliases the TESTED effect. ***")
    say("               Continuing AS IF a single batch (no batch term); results are batch-uncorrected.")
    say("               Spike-in is DISABLED for this test (cross-batch carryover cannot be separated).")
  }

  # ---- resolve optional inputs ----
  peaks_gr <- resolve_peaks(opt, meta)
  igg_bams <- resolve_igg(opt, meta)
  spike_v  <- resolve_spikein(opt, meta)
  sf_v     <- resolve_normfactors(opt, meta)

  gate_ok     <- ab$allow_spikein || isTRUE(opt$force_norm)
  use_spikein <- !is.null(spike_v) && gate_ok
  use_sf      <- !is.null(sf_v)    && gate_ok
  if (!is.null(spike_v) && !gate_ok)
    say("[SPIKEIN] requested but DISABLED (batch not separable; pass --force_norm to override).")
  if (!is.null(sf_v) && !gate_ok)
    say("[NORMFACTORS] requested but DISABLED (batch not separable; pass --force_norm to override).")
  if ((use_spikein || use_sf) && !ab$allow_spikein && isTRUE(opt$force_norm)) {
    say("[FORCE_NORM] *** Applying supplied normalization despite a non-separable batch. ***")
    say("            The GLOBAL magnitude is only as good as what you anchored it to (e.g.")
    say("            ChIPseqSpikeInFree + western); it is NOT separable from batch by these data.")
    say("            Loci-level / relative structure is the interpretable part; treat the global")
    say("            up/down balance as externally anchored, not as proven by the sequencing.")
  }

  # ---- count + filter ----
  param <- make_param(opt)
  cf <- count_and_filter(meta$bamReads, opt, param, igg_bams)
  win <- cf$win; bins <- cf$bins
  win_gr <- rowRanges(win)

  # ---- independent baseline group counted in the SAME windows ----
  indep_logcpm <- NULL
  if (!is.null(opt$baseline_group)) {
    bg <- parse_group(opt$baseline_group)
    bgrows <- ss[ss$.fac == tolower(opt$mark) & ss$.cond == bg$cond & ss$.treat == bg$treat, , drop = FALSE]
    if (!nrow(bgrows)) stop("--baseline_group '", opt$baseline_group, "' matched no samples.")
    say(sprintf("[STRAT] independent baseline group '%s': %d samples.",
                opt$baseline_group, nrow(bgrows)))
    bgc <- regionCounts(bgrows$bamReads, regions = win_gr, param = param)
    indep_logcpm <- aveLogCPM(assay(bgc))
  }

  # ---- build DGEList; name columns by SampleID and carry plotting metadata ----
  y0 <- norm_tmm(win, bins)
  colnames(y0)      <- meta$SampleID          # real IDs -> fixes spike-in keying + MDS labels
  y0$samples$cond   <- meta$.cond             # for MDS shape
  y0$samples$treat  <- meta$.treat            # for MDS colour
  y0$samples$batch  <- as.character(meta$Batch)  # for MDS batch companion
  if (opt$test == "contrast") {
    y0$samples$group <- meta$Group
    y0$samples$group_role <- meta$.role
  } else {
    y0$samples$group <- factor(paste(meta$Condition, meta$Treatment, sep = "_"))
    y0$samples$peak_ref <- attr(meta, "peak_ref")
  }

  # ---- run libsize (always) ----
  res_lib <- run_norm(y0, ab, opt, win_gr, peaks_gr, indep_logcpm, tag, "libsize", opt$outdir)

  # ---- run spike-in (if requested & allowed) ----
  res_spk <- NULL
  if (use_spikein) {
    y1 <- set_spikein(y0, spike_v)
    say("[SPIKEIN] norm.factors (spikein_reads / lib.size, centred):")
    for (i in seq_len(ncol(y1)))
      say(sprintf("   %-22s reads=%-12.0f nf=%.3f", colnames(y1)[i],
                  spike_v[colnames(y1)[i]], y1$samples$norm.factors[i]))
    res_spk <- run_norm(y1, ab, opt, win_gr, peaks_gr, indep_logcpm, tag, "spikein", opt$outdir)
    norm_compare_fig(res_lib$df, res_spk$df, tag, "spikein", opt$outdir, opt)
    say("[COMPARE] libsize vs spike-in: loci-level structure should agree; a global")
    say("          up/down shift is the spike-in correction (see *_spikein_norm_compare).")
  }

  # ---- run precomputed scaling factors (ChIPseqSpikeInFree etc.) if requested & allowed ----
  res_sf <- NULL
  if (use_sf) {
    y2 <- set_normfactors(y0, sf_v)
    say("[NORMFACTORS] norm.factors (supplied scaling factors, geo-mean-centred):")
    for (i in seq_len(ncol(y2)))
      say(sprintf("   %-22s SF=%-10.4f nf=%.3f", colnames(y2)[i],
                  sf_v[colnames(y2)[i]], y2$samples$norm.factors[i]))
    res_sf <- run_norm(y2, ab, opt, win_gr, peaks_gr, indep_logcpm, tag, "sf", opt$outdir)
    norm_compare_fig(res_lib$df, res_sf$df, tag, "sf", opt$outdir, opt)
    say("[COMPARE] libsize vs sf: loci-level structure should agree; the global up/down shift")
    say("          is the supplied (e.g. ChIPseqSpikeInFree) correction (see *_sf_norm_compare).")
  }

  # ---- summary report ----
  writeLines(log_lines$txt, file.path(opt$outdir, paste0(tag, "_summary.txt")))
  say(sprintf("[DONE] outputs in %s", normalizePath(opt$outdir)))
}

main()