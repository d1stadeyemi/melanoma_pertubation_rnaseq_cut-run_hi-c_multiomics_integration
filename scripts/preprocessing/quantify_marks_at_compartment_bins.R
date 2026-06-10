# =============================================================================
# quantify_marks_at_compartment_bins.R
#   Per-250kb-bin histone-mark signal, normalized consistently with csaw.R,
#   for integration with the compartment-reversal analysis.
# =============================================================================
#
# PROJECT: A375 melanoma parental-vs-B10 constricted-migration multi-omics.
#
# WHAT THIS SCRIPT IS FOR ------------------------------------------------------
#   Track 1 (compartment <-> H3K27me3 reversal) and Track 3 (predict compartment
#   change from multi-mark changes) both need a CONTINUOUS, per-bin mark signal
#   that lines up 1:1 with the 250kb compartment bins used for PC1. Your existing
#   signal tables are CPM-normalized, which is a one-way street: you cannot re-
#   apply TMM (to match the csaw differential calls) or spike-in factors (the
#   flag you are keeping) to an already-CPM table. So we recount RAW from the
#   BAMs at the exact compartment-bin coordinates and normalize here, leaving a
#   single place to swap in spike-in later.
#
# WHY THIS MATCHES csaw.R / the PRC2 paper -----------------------------------
#   Normalization is the SAME recipe used for the differential calls and in
#   Romero et al. 2024 (Nat Commun): TMM (trimmed mean of M-values) composition
#   factors computed on 10kb genome-wide background bins, applied to the signal.
#   Putting the continuous deltas on the same basis as the csaw FDR calls removes
#   the "why is signal CPM but tests TMM?" question.
#
# KEY DESIGN POINTS -----------------------------------------------------------
#   * Counting is done over YOUR EXACT compartment bins (regionCounts over the
#     coordinates in --compartments), NOT a fresh genome tiling, so the output
#     merges 1:1 with PC1 / Reversal_classified.txt. A generic 250kb tiling would
#     be offset and lose bins to partial overlaps.
#   * Normalization is computed WITHIN EACH MARK separately. You never TMM one
#     antibody against another (H3K27me3 vs H3K9me3 have different dynamics);
#     each mark gets its own background bins and its own factors.
#   * Sex chromosomes (+ chrM) are excluded from the NORMALIZATION background by
#     default (as deepTools --ignoreForNormalization chrX chrY does), but the
#     compartment bins themselves are all retained in the output.
#   * Replicate signal is averaged in LOG space (= geometric mean of CPM), and a
#     contrast delta is the difference of those means (a log2 fold-change):
#       migration  = mean logCPM(B10|none)      - mean logCPM(Parental|none)
#       b10_treat  = mean logCPM(B10|UNC)        - mean logCPM(B10|none)
#       par_treat  = mean logCPM(Parental|UNC)   - mean logCPM(Parental|none)
#   * Group matching is CASE-INSENSITIVE on Condition|Treatment, exactly like
#     csaw.R.
#
# SPIKE-IN (future-proof) -----------------------------------------------------
#   --spikein factors.csv (columns: SampleID,factor). When supplied, those
#   factors REPLACE the TMM factors (scaled to geometric mean 1), per sample.
#   Because we kept raw counts, this is the only change needed to move any/all
#   marks onto an absolute (spike-in) basis later. Until then, log-CPM/TMM
#   answers RELATIVE/redistribution questions only (see csaw.R note).
#
# OUTPUTS (into --outdir) -----------------------------------------------------
#   <mark>_logCPM_250kb.txt          coordinate + one column per SampleID
#   marks_condition_mean_logCPM.txt  coordinate + <mark>_<cond>_<treat> means
#                                    (use the *_parental_none column as the
#                                     baseline H3K27me3 for peak/background
#                                     stratification in Track 1)
#   marks_delta_250kb.txt            coordinate + <mark>_migration / _b10_treat /
#                                    _par_treat  (the centerpiece: Track 1 reads
#                                    the H3K27me3_* columns; Track 3 reads the
#                                    *_b10_treat columns as ML features)
#   normalization_factors.txt        per-sample lib.size + norm.factor (record)
#
# EXAMPLE ---------------------------------------------------------------------
#   Rscript quantify_marks_at_compartment_bins.R \
#       --sheet        data/processed/sample_sheet.csv \
#       --compartments data/processed/A375-CompartmentsAcrossConditions.txt \
#       --coord_col    coordinate \
#       --outdir       results/exploratory/mark_signal_250kb
#   # later, with spike-in:
#   #   ... --spikein data/processed/spikein_factors.csv
# =============================================================================

suppressWarnings(suppressMessages({
  library(optparse)
}))

# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
parse_cli <- function() {
  option_list <- list(
    make_option("--sheet", type = "character",
                help = "Sample sheet CSV (same one csaw.R uses: needs columns "),
    make_option("--compartments", type = "character",
                help = "Compartment table whose --coord_col defines the 250kb bins"),
    make_option("--coord_col", type = "character", default = "coordinate",
                help = "Coordinate column in --compartments [default %default]"),
    make_option("--outdir", type = "character", default = "mark_signal_250kb"),
    make_option("--marks", type = "character", default = NULL,
                help = "Comma-separated Factor values to quantify; default = all distinct Factor values in the sheet"),
    make_option("--frag_max", type = "integer", default = 700L,
                help = "Max fragment length (match Bowtie2 --X / csaw.R) [default %default]"),
    make_option("--bg_bin", type = "integer", default = 10000L,
                help = "Background bin width (bp) for TMM normalization [default %default]"),
    make_option("--prior_count", type = "double", default = 2.0,
                help = "edgeR prior.count for log-CPM (stabilizes low bins) [default %default]"),
    make_option("--blacklist", type = "character", default = NULL,
                help = "Optional ENCODE blacklist BED of regions to discard"),
    make_option("--spikein", type = "character", default = NULL,
                help = "Optional CSV (SampleID,factor); REPLACES TMM when given"),
    make_option("--norm_exclude", type = "character", default = "chrX,chrY,chrM",
                help = "Chroms excluded from the TMM background only [default %default]"),
    # condition/treatment label strings (override if your sheet uses other casing/words)
    make_option("--cond_parental", type = "character", default = "Parental"),
    make_option("--cond_b10",      type = "character", default = "B10"),
    make_option("--treat_none",    type = "character", default = "none"),
    make_option("--treat_unc",     type = "character", default = "UNC")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$sheet))        stop("--sheet is required.")
  if (is.null(opt$compartments)) stop("--compartments is required.")
  opt
}

# -----------------------------------------------------------------------------
# Parse a coordinate string ("chr1:1-250000", "chr1-1-250000", "chr1_1_250000",
# etc.) into a one-row data.frame(chrom,start,end). GRanges is 1-based inclusive;
# we take the parsed start/end as-is (any 1bp edge convention is negligible at
# 250kb). Keeps the ORIGINAL string so the output merges back exactly.
# -----------------------------------------------------------------------------
parse_coord <- function(x) {
  x <- as.character(x)
  m <- regmatches(x, regexec("(chr[0-9A-Za-z]+)[: _-]+([0-9]+)[: _-]+([0-9]+)", x))
  chrom <- vapply(m, function(z) if (length(z) == 4) z[2] else NA_character_, character(1))
  start <- vapply(m, function(z) if (length(z) == 4) as.numeric(z[3]) else NA_real_, numeric(1))
  end   <- vapply(m, function(z) if (length(z) == 4) as.numeric(z[4]) else NA_real_, numeric(1))
  data.frame(coordinate = x, chrom = chrom, start = start, end = end,
             stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------
# Build the GRanges of compartment bins from the compartment table.
# -----------------------------------------------------------------------------
load_bins <- function(opt) {
  library(GenomicRanges)
  comp <- read.delim(opt$compartments, stringsAsFactors = FALSE, check.names = FALSE)
  if (!opt$coord_col %in% colnames(comp))
    stop("--coord_col '", opt$coord_col, "' not in ", opt$compartments,
         ". Columns: ", paste(head(colnames(comp), 40), collapse = ", "))
  cd <- parse_coord(comp[[opt$coord_col]])
  bad <- is.na(cd$chrom) | is.na(cd$start) | is.na(cd$end)
  if (any(bad)) {
    message("[WARN] ", sum(bad), " coordinate(s) could not be parsed and are dropped. ",
            "First unparsed: ", paste(head(cd$coordinate[bad], 3), collapse = " | "))
    cd <- cd[!bad, ]
  }
  message("[INFO] Parsed ", nrow(cd), " compartment bins. Examples:")
  print(utils::head(cd, 3))
  gr <- GRanges(seqnames = cd$chrom,
                ranges   = IRanges(start = pmax(cd$start, 1), end = cd$end))
  gr$coordinate <- cd$coordinate           # carry the exact key
  gr
}

# -----------------------------------------------------------------------------
# Load the sample sheet and check required columns.
# -----------------------------------------------------------------------------
load_sheet <- function(opt) {
  ss <- read.csv(opt$sheet, stringsAsFactors = FALSE)
  need <- c("SampleID", "Factor", "Condition", "Treatment", "Replicate", "bamReads")
  miss <- setdiff(need, colnames(ss))
  if (length(miss))
    stop("Sample sheet missing column(s): ", paste(miss, collapse = ", "),
         "\nHas: ", paste(colnames(ss), collapse = ", "))
  ss$Condition <- trimws(ss$Condition)
  ss$Treatment <- trimws(ss$Treatment)
  ss
}

# -----------------------------------------------------------------------------
# Normalize ONE mark and return per-bin log-CPM (bins x samples) plus a record
# of the factors used.
#   * regionCounts over the compartment bins  -> region counts; colData$totals
#     holds full library sizes (used as lib.size by asDGEList -> proper CPM).
#   * windowCounts(bin=TRUE, bg_bin) on AUTOSOMES -> TMM factors via normFactors,
#     unless spike-in factors are supplied (then those replace TMM).
# -----------------------------------------------------------------------------
quantify_mark <- function(mark, ss, bins, opt, spikein_df, norm_excl) {
  library(csaw); library(edgeR); library(GenomicRanges)

  sm <- ss[tolower(ss$Factor) == tolower(mark), ]
  if (nrow(sm) == 0) { message("[WARN] no samples for mark ", mark, "; skipping."); return(NULL) }
  sm <- sm[order(sm$Condition, sm$Treatment, sm$Replicate), ]
  message("[INFO] === ", mark, " : ", nrow(sm), " samples ===")

  param <- readParam(pe = "both", max.frag = opt$frag_max)
  if (!is.null(opt$blacklist)) {
    library(rtracklayer)
    param <- reform(param, discard = import(opt$blacklist))
  }

  # Signal: counts in the exact compartment bins.
  region.counts <- regionCounts(sm$bamReads, regions = bins, param = param)
  colnames(region.counts) <- sm$SampleID

  y <- asDGEList(region.counts)          # lib.size = full library totals
  y$samples$group <- paste(tolower(sm$Condition), tolower(sm$Treatment), sep = "|")

  # Normalization factors.
  if (!is.null(spikein_df)) {
    facs <- spikein_df$factor[match(sm$SampleID, spikein_df$SampleID)]
    if (any(is.na(facs)))
      stop("Spike-in file missing factors for: ",
           paste(sm$SampleID[is.na(facs)], collapse = ", "))
    y$samples$norm.factors <- facs / exp(mean(log(facs)))   # geo-mean 1
    message("[INFO] ", mark, ": spike-in normalization.")
  } else {
    binned <- windowCounts(sm$bamReads, bin = TRUE, width = opt$bg_bin, param = param)
    # exclude sex/M chroms from the background used to estimate composition
    keep_bg <- !(as.character(seqnames(rowRanges(binned))) %in% norm_excl)
    binned  <- binned[keep_bg, ]
    y$samples$norm.factors <- normFactors(binned, se.out = TRUE)$norm.factors
    message("[INFO] ", mark, ": TMM (", sum(keep_bg), " autosomal ", opt$bg_bin,
            "bp background bins).")
  }

  logcpm <- cpm(y, log = TRUE, prior.count = opt$prior_count)   # bins x samples
  colnames(logcpm) <- sm$SampleID

  list(mark = mark, samples = sm,
       logcpm = logcpm,
       coord  = bins$coordinate,
       norm   = data.frame(SampleID = sm$SampleID, Factor = mark,
                           lib.size = y$samples$lib.size,
                           norm.factor = y$samples$norm.factors,
                           stringsAsFactors = FALSE))
}

# -----------------------------------------------------------------------------
# Condition-mean log-CPM and the three standard contrast deltas for one mark.
# group keys are lowercased "condition|treatment".
# -----------------------------------------------------------------------------
summarize_mark <- function(q, opt) {
  sm <- q$samples
  grp <- paste(tolower(sm$Condition), tolower(sm$Treatment), sep = "|")
  conds <- unique(grp)

  # mean log-CPM per condition (geometric mean of CPM)
  cond_mean <- sapply(conds, function(g) {
    cols <- which(grp == g)
    if (length(cols) == 1) q$logcpm[, cols] else rowMeans(q$logcpm[, cols, drop = FALSE])
  })
  colnames(cond_mean) <- conds

  key <- function(cond, treat) paste(tolower(cond), tolower(treat), sep = "|")
  g <- function(cond, treat) {
    k <- key(cond, treat)
    if (!k %in% colnames(cond_mean)) {
      message("[WARN] ", q$mark, ": missing group '", k,
              "' -> its contrasts will be NA. Present: ",
              paste(colnames(cond_mean), collapse = ", "))
      return(rep(NA_real_, nrow(cond_mean)))
    }
    cond_mean[, k]
  }

  delta <- data.frame(
    migration = g(opt$cond_b10, opt$treat_none) - g(opt$cond_parental, opt$treat_none),
    b10_treat = g(opt$cond_b10, opt$treat_unc)  - g(opt$cond_b10, opt$treat_none),
    par_treat = g(opt$cond_parental, opt$treat_unc) - g(opt$cond_parental, opt$treat_none),
    check.names = FALSE)

  list(cond_mean = cond_mean, delta = delta)
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
main <- function() {
  opt <- parse_cli()
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
  norm_excl <- trimws(strsplit(opt$norm_exclude, ",")[[1]])

  ss   <- load_sheet(opt)
  bins <- load_bins(opt)

  marks <- if (is.null(opt$marks)) sort(unique(ss$Factor)) else
           trimws(strsplit(opt$marks, ",")[[1]])
  message("[INFO] Marks to quantify: ", paste(marks, collapse = ", "))

  spikein_df <- NULL
  if (!is.null(opt$spikein)) {
    spikein_df <- read.csv(opt$spikein, stringsAsFactors = FALSE)
    if (!all(c("SampleID", "factor") %in% colnames(spikein_df)))
      stop("--spikein CSV must have columns: SampleID, factor")
  } else {
    message("[NOTE] No --spikein; TMM/library-size only. ",
            "Cross-sample deltas are RELATIVE (redistribution), not absolute amounts.")
  }

  coord <- bins$coordinate
  delta_all <- data.frame(coordinate = coord, stringsAsFactors = FALSE)
  mean_all  <- data.frame(coordinate = coord, stringsAsFactors = FALSE)
  norm_rec  <- list()

  for (mk in marks) {
    q <- quantify_mark(mk, ss, bins, opt, spikein_df, norm_excl)
    if (is.null(q)) next

    # per-sample log-CPM table for this mark
    lc <- data.frame(coordinate = q$coord, q$logcpm, check.names = FALSE)
    write.table(lc, file.path(opt$outdir, paste0(mk, "_logCPM_250kb.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)

    s <- summarize_mark(q, opt)
    # append condition means with mark prefix
    cm <- s$cond_mean
    colnames(cm) <- paste0(mk, "_", gsub("[|]", "_", colnames(cm)))
    mean_all <- cbind(mean_all, cm)
    # append the three deltas with mark prefix
    dl <- s$delta
    colnames(dl) <- paste0(mk, "_", colnames(dl))
    delta_all <- cbind(delta_all, dl)

    norm_rec[[mk]] <- q$norm
  }

  write.table(mean_all, file.path(opt$outdir, "marks_condition_mean_logCPM.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(delta_all, file.path(opt$outdir, "marks_delta_250kb.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  if (length(norm_rec))
    write.table(do.call(rbind, norm_rec),
                file.path(opt$outdir, "normalization_factors.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)

  message("\n[INFO] Done -> ", normalizePath(opt$outdir))
  message("[INFO] Centerpiece for Track 1/3: ",
          file.path(opt$outdir, "marks_delta_250kb.txt"))
}

main()