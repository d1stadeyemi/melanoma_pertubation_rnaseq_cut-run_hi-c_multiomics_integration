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
#   READ FILTERING now also matches csaw.R exactly: minq=20, dedup=TRUE,
#   pe="both", max.frag=frag_max. Putting BOTH the read filter AND the
#   normalization on the same basis as the csaw FDR calls removes the
#   "why is signal counted/normalized differently from the tests?" question.
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
# READ FILTERING PARITY (the fix) ---------------------------------------------
#   readParam(minq=20, dedup=TRUE, pe="both", max.frag=frag_max) is applied to
#   BOTH the signal counts (regionCounts over compartment bins) and the
#   background counts (windowCounts for TMM), identical to csaw.R::make_param.
#   These filters are safe even if your BAMs were already MAPQ/dup-filtered
#   (re-applying is idempotent). If your BAMs are single-end, override with
#   --pe none (then max.frag is ignored by csaw).
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
#   <mark>_QC.pdf                    MA plot (migration contrast) + replicate
#                                    correlation heat (normalization sanity)
#   QC_replicate_correlation.txt     within-group replicate Pearson r per mark
#   MANIFEST.md                      provenance: args, inputs+md5, git, versions
#
# EXAMPLE ---------------------------------------------------------------------
#   Rscript quantify_marks_at_compartment_bins.R \
#       --sheet        data/processed/sample_sheet.csv \
#       --compartments data/processed/A375-CompartmentsAcrossConditions.txt \
#       --coord_col    coordinate \
#       --expect_bins   10576 \
#       --reversal_check results/exploratory/reversal/Reversal_classified.txt \
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
                help = "Sample sheet CSV (same one csaw.R uses): SampleID,Factor,Condition,Treatment,Replicate,bamReads"),
    make_option("--compartments", type = "character",
                help = "Compartment table whose --coord_col defines the 250kb bins"),
    make_option("--coord_col", type = "character", default = "coordinate",
                help = "Coordinate column in --compartments [default %default]"),
    make_option("--outdir", type = "character", default = "mark_signal_250kb"),
    make_option("--marks", type = "character", default = NULL,
                help = "Comma-separated Factor values to quantify; default = all distinct Factor values in the sheet"),
    make_option("--pe", type = "character", default = "both",
                help = "Paired-end mode for readParam: 'both'/'none'/'first'/'second' (match csaw.R) [default %default]"),
    make_option("--minq", type = "integer", default = 20L,
                help = "Minimum MAPQ (match csaw.R) [default %default]"),
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
    # integrity / provenance
    make_option("--expect_bins", type = "integer", default = NULL,
                help = "If set, ASSERT the final bin count equals this (e.g. 10576). [default: none]"),
    make_option("--reversal_check", type = "character", default = NULL,
                help = "Optional Reversal_classified.txt; assert its coordinates join 1:1 with the output bins."),
    make_option("--no_qc", action = "store_true", default = FALSE,
                help = "Disable per-mark QC plots (MA + replicate correlation)."),
    # condition/treatment label strings (override if your sheet uses other casing/words)
    make_option("--cond_parental", type = "character", default = "Parental"),
    make_option("--cond_b10",      type = "character", default = "B10"),
    make_option("--treat_none",    type = "character", default = "none"),
    make_option("--treat_unc",     type = "character", default = "UNC")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$sheet))        stop("--sheet is required.")
  if (is.null(opt$compartments)) stop("--compartments is required.")
  if (!opt$pe %in% c("both", "none", "first", "second"))
    stop("--pe must be one of both/none/first/second.")
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
  if (any(duplicated(cd$coordinate)))
    stop("[INTEGRITY] --compartments has duplicated coordinate strings (", 
         sum(duplicated(cd$coordinate)), "); bins must be unique.")
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
  miss_bam <- ss$bamReads[!file.exists(ss$bamReads)]
  if (length(miss_bam))
    stop("bamReads not found:\n  ", paste(head(miss_bam, 10), collapse = "\n  "))
  ss
}

# -----------------------------------------------------------------------------
# Normalize ONE mark and return per-bin log-CPM (bins x samples) plus a record
# of the factors used.
#   * regionCounts over the compartment bins  -> region counts; colData$totals
#     holds full library sizes (used as lib.size by asDGEList -> proper CPM).
#   * windowCounts(bin=TRUE, bg_bin) on AUTOSOMES -> TMM factors via normFactors,
#     unless spike-in factors are supplied (then those replace TMM).
#   readParam matches csaw.R exactly: minq, dedup=TRUE, pe, max.frag, discard.
# -----------------------------------------------------------------------------
quantify_mark <- function(mark, ss, bins, opt, spikein_df, norm_excl, discard_gr) {
  library(csaw); library(edgeR); library(GenomicRanges)

  sm <- ss[tolower(ss$Factor) == tolower(mark), ]
  if (nrow(sm) == 0) { message("[WARN] no samples for mark ", mark, "; skipping."); return(NULL) }
  sm <- sm[order(sm$Condition, sm$Treatment, sm$Replicate), ]
  message("[INFO] === ", mark, " : ", nrow(sm), " samples ===")

  # readParam IDENTICAL to csaw.R::make_param (minq=20, dedup=TRUE, pe, max.frag,
  # discard). Applied to BOTH signal and background counting below.
  param <- readParam(minq = as.integer(opt$minq), dedup = TRUE, pe = opt$pe,
                     max.frag = opt$frag_max, discard = discard_gr)

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
# QC for one mark: replicate correlation (normalization sanity) + an MA plot of
# the migration contrast. Written to <mark>_QC.pdf. Returns the within-group
# replicate-correlation rows for the combined QC table. Guarded so QC never
# kills the run.
# -----------------------------------------------------------------------------
qc_mark <- function(q, opt) {
  out_rows <- list()
  tryCatch({
    sm  <- q$samples
    lc  <- q$logcpm
    grp <- paste(tolower(sm$Condition), tolower(sm$Treatment), sep = "|")

    # within-group replicate Pearson r (only groups with >= 2 reps)
    for (gg in unique(grp)) {
      cols <- which(grp == gg)
      if (length(cols) >= 2) {
        cm <- cor(lc[, cols, drop = FALSE], use = "pairwise.complete.obs")
        r  <- mean(cm[upper.tri(cm)])
        out_rows[[length(out_rows) + 1]] <- data.frame(
          Factor = q$mark, group = gg, n_rep = length(cols),
          mean_rep_pearson = round(r, 4), stringsAsFactors = FALSE)
      }
    }

    pdf(file.path(opt$outdir, paste0(q$mark, "_QC.pdf")), width = 9, height = 4.5)
    op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

    # MA plot of the migration contrast (B10|none vs Parental|none) if available
    key <- function(cond, treat) paste(tolower(cond), tolower(treat), sep = "|")
    kb <- key(opt$cond_b10, opt$treat_none); kp <- key(opt$cond_parental, opt$treat_none)
    mean_grp <- function(k) {
      cols <- which(grp == k)
      if (!length(cols)) return(NULL)
      if (length(cols) == 1) lc[, cols] else rowMeans(lc[, cols, drop = FALSE])
    }
    mb <- mean_grp(kb); mp <- mean_grp(kp)
    if (!is.null(mb) && !is.null(mp)) {
      A <- 0.5 * (mb + mp); M <- mb - mp
      smoothScatter(A, M, nbin = 220, xlab = "A  0.5*(B10+Parental) logCPM",
                    ylab = "M  B10 - Parental (migration)",
                    main = paste0(q$mark, "  MA (migration)"))
      abline(h = 0, col = "red", lwd = 1.5)
      abline(h = median(M, na.rm = TRUE), col = "blue", lty = 2)
      legend("topright", bty = "n", cex = 0.8,
             legend = sprintf("median M = %+.3f", median(M, na.rm = TRUE)))
    } else {
      plot.new(); title(paste0(q$mark, ": migration groups absent for MA"))
    }

    # replicate correlation heat (all samples)
    if (ncol(lc) >= 2) {
      cm <- cor(lc, use = "pairwise.complete.obs")
      image(seq_len(ncol(cm)), seq_len(ncol(cm)), cm, axes = FALSE,
            xlab = "", ylab = "", main = paste0(q$mark, "  sample logCPM correlation"),
            zlim = c(min(cm), 1), col = hcl.colors(32, "Blues", rev = TRUE))
      axis(1, at = seq_len(ncol(cm)), labels = colnames(cm), las = 2, cex.axis = 0.55)
      axis(2, at = seq_len(ncol(cm)), labels = colnames(cm), las = 2, cex.axis = 0.55)
    } else {
      plot.new(); title("single sample: no correlation")
    }
    par(op); dev.off()
    message("[QC] ", q$mark, " -> ", paste0(q$mark, "_QC.pdf"))
  }, error = function(e) message("[QC] ", q$mark, " QC skipped: ", conditionMessage(e)))

  if (length(out_rows)) do.call(rbind, out_rows) else NULL
}

# -----------------------------------------------------------------------------
# Provenance manifest: args, input files (+ md5), git commit, package versions.
# -----------------------------------------------------------------------------
write_manifest <- function(opt, marks, norm_rec) {
  path <- file.path(opt$outdir, "MANIFEST.md")
  git_commit <- tryCatch(
    trimws(system("git rev-parse HEAD 2>/dev/null", intern = TRUE))[1],
    error = function(e) NA_character_)
  inputs <- c(opt$sheet, opt$compartments, opt$blacklist, opt$spikein, opt$reversal_check)
  inputs <- inputs[!is.null(inputs) & !is.na(inputs)]
  md5 <- tryCatch(tools::md5sum(inputs), error = function(e) setNames(rep(NA, length(inputs)), inputs))
  pkg_ver <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "NA")

  con <- file(path, open = "wt")
  writeLines(c(
    "# quantify_marks_at_compartment_bins.R - provenance",
    paste0("- timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("- git_commit: ", ifelse(length(git_commit) && !is.na(git_commit), git_commit, "NA (not a git repo)")),
    paste0("- Rscript: ", R.version.string),
    paste0("- packages: csaw ", pkg_ver("csaw"), ", edgeR ", pkg_ver("edgeR"),
           ", GenomicRanges ", pkg_ver("GenomicRanges")),
    "",
    "## parameters",
    paste0("- pe=", opt$pe, "  minq=", opt$minq, "  dedup=TRUE  max.frag=", opt$frag_max),
    paste0("- bg_bin=", opt$bg_bin, "  prior_count=", opt$prior_count,
           "  norm_exclude=", opt$norm_exclude),
    paste0("- normalization=", ifelse(is.null(opt$spikein), "TMM (composition)", "spike-in (supplied factors)")),
    paste0("- marks=", paste(marks, collapse = ", ")),
    "",
    "## inputs (md5)",
    if (length(inputs)) paste0("- ", basename(inputs), ": ", md5, "  (", inputs, ")") else "- none",
    "",
    "## read-filter parity",
    "- readParam(minq, dedup=TRUE, pe, max.frag, discard) is IDENTICAL to csaw.R::make_param,",
    "  applied to BOTH regionCounts (signal) and windowCounts (TMM background)."
  ), con)
  close(con)
  message("[INFO] provenance -> ", path)
}

# -----------------------------------------------------------------------------
# Integrity check: output bins unique, count as expected, and (optionally) a 1:1
# coordinate join with Reversal_classified.txt.
# -----------------------------------------------------------------------------
integrity_check <- function(delta_all, opt) {
  n <- nrow(delta_all)
  if (any(duplicated(delta_all$coordinate)))
    stop("[INTEGRITY] marks_delta_250kb has duplicated coordinates.")
  message("[INTEGRITY] marks_delta_250kb rows (bins): ", n,
          "  (+1 header line = ", n + 1, " in `wc -l`)")
  if (!is.null(opt$expect_bins)) {
    if (n != opt$expect_bins)
      stop("[INTEGRITY] expected ", opt$expect_bins, " bins but produced ", n, ".")
    message("[INTEGRITY] bin count matches --expect_bins (", opt$expect_bins, "). OK")
  }
  if (!is.null(opt$reversal_check)) {
    if (!file.exists(opt$reversal_check))
      stop("[INTEGRITY] --reversal_check file not found: ", opt$reversal_check)
    rv <- read.delim(opt$reversal_check, stringsAsFactors = FALSE, check.names = FALSE)
    cc <- opt$coord_col
    if (!cc %in% colnames(rv)) {
      # fall back to a column literally named 'coordinate'
      cc <- if ("coordinate" %in% colnames(rv)) "coordinate" else
            stop("[INTEGRITY] coord column '", opt$coord_col,
                 "' not in reversal file; columns: ", paste(head(colnames(rv), 20), collapse = ", "))
    }
    rc <- as.character(rv[[cc]]); dc <- as.character(delta_all$coordinate)
    only_marks <- setdiff(dc, rc); only_rev <- setdiff(rc, dc)
    message("[INTEGRITY] reversal join: ", length(intersect(dc, rc)), " shared | ",
            length(only_marks), " only-in-marks | ", length(only_rev), " only-in-reversal")
    if (length(only_marks) || length(only_rev)) {
      message("  first only-in-marks: ", paste(head(only_marks, 3), collapse = " | "))
      message("  first only-in-reversal: ", paste(head(only_rev, 3), collapse = " | "))
      message("[INTEGRITY][WARN] coordinate sets are NOT identical 1:1. Downstream merges ",
              "will drop the non-shared bins. Confirm both came from the same compartment table.")
    } else {
      message("[INTEGRITY] coordinates join 1:1 with the reversal table. OK")
    }
  }
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

  # build the discard GRanges once (blacklist), matching csaw.R's discard handling
  discard_gr <- GenomicRanges::GRanges()
  if (!is.null(opt$blacklist)) {
    if (!file.exists(opt$blacklist)) stop("Blacklist not found: ", opt$blacklist)
    suppressWarnings(suppressMessages(library(rtracklayer)))
    discard_gr <- import(opt$blacklist)
    message("[INFO] blacklist: ", length(discard_gr), " regions discarded (matches csaw.R --blacklist).")
  }

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
  qc_rows   <- list()

  for (mk in marks) {
    q <- quantify_mark(mk, ss, bins, opt, spikein_df, norm_excl, discard_gr)
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
    if (!isTRUE(opt$no_qc)) {
      qr <- qc_mark(q, opt)
      if (!is.null(qr)) qc_rows[[mk]] <- qr
    }
  }

  write.table(mean_all, file.path(opt$outdir, "marks_condition_mean_logCPM.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(delta_all, file.path(opt$outdir, "marks_delta_250kb.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  if (length(norm_rec))
    write.table(do.call(rbind, norm_rec),
                file.path(opt$outdir, "normalization_factors.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  if (length(qc_rows))
    write.table(do.call(rbind, qc_rows),
                file.path(opt$outdir, "QC_replicate_correlation.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)

  # integrity + provenance
  integrity_check(delta_all, opt)
  write_manifest(opt, marks, norm_rec)

  message("\n[INFO] Done -> ", normalizePath(opt$outdir))
  message("[INFO] Centerpiece for Track 1/3: ",
          file.path(opt$outdir, "marks_delta_250kb.txt"))
}

main()