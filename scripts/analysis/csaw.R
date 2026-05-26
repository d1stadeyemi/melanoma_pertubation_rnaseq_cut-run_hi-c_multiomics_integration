# =============================================================================
# csaw.R — Window-based differential binding for BROAD histone marks
# =============================================================================
#
# PROJECT: A375 melanoma parental-vs-B10 constricted-migration multi-omics.
#
# WHAT THIS SCRIPT IS FOR ------------------------------------------------------
#   Broad chromatin marks (H3K27me3, H3K9me3, H3K4me1) form diffuse domains, not
#   sharp peaks. Peak-based differential tools depend on peak calls that are
#   unreliable for such domains, so we instead TILE the genome into fixed windows
#   and test each window for differential coverage between two conditions, using
#   the edgeR negative-binomial model. Replicates supply the within-group
#   variance that turns a fold-change into a real, FDR-controlled p-value.
#   (The punctate mark H3K27Ac is handled separately by diffbind.R.)
#
# WHY A STATISTICAL TEST (not a difference of means) --------------------------
#   A difference/ratio of averaged signal has no measure of reproducibility — it
#   returns the same value whether replicates agree or disagree. The NB model
#   divides the effect by replicate-derived variability, which is what a Q1
#   journal requires to claim a region "changed".
#
# CONTRASTS (explicit, over Condition|Treatment) ------------------------------
#   Migration:      --group1 "B10|none"      --group2 "Parental|none"
#   UNC in parent:  --group1 "Parental|UNC"  --group2 "Parental|none"
#   UNC in B10:     --group1 "B10|UNC"       --group2 "B10|none"
#   group1 = numerator (positive log2FC => higher in group1).
#   Matching is CASE-INSENSITIVE: "Parental"=="parental", "UNC"=="unc".
#
# BATCH -----------------------------------------------------------------------
#   If the sheet's Batch column has >1 level among the contrast samples, the
#   design becomes ~Batch + group (batch as a blocking covariate). If batch is
#   perfectly confounded with the contrast, the design is rank-deficient and
#   edgeR will (correctly) refuse — that situation cannot be corrected, only
#   acknowledged.
#
# SPIKE-IN (future-proof) -----------------------------------------------------
#   --spikein factors.csv  (columns: SampleID,factor). Until supplied, TMM /
#   library-size normalization runs. For broad REPRESSIVE marks, library-size can
#   MASK a true genome-wide change; such results are provisional until spike-in.
#
# EXAMPLE USAGE ---------------------------------------------------------------
#   Rscript csaw.R \
#       --sheet data/processed/sample_sheet.csv \
#       --outdir results/diff_analysis/csaw \
#       --mark H3K27me3 \
#       --group1 "B10|none" --group2 "Parental|none"
#
#   # with spike-in factors and a bigger window for a very broad mark:
#   Rscript csaw.R --sheet ss.csv --mark H3K9me3 \
#       --group1 "B10|none" --group2 "Parental|none" \
#       --window 5000 --spacing 1000 --spikein spikein_factors.csv
#
# OUTPUTS (into --outdir) -----------------------------------------------------
#   <tag>_<norm>_all_regions.csv   all merged regions: coords, logFC, PValue, FDR
#   <tag>_<norm>_significant.csv   subset passing FDR & |log2FC| thresholds
#   <tag>_<norm>_QC.pdf            MDS (sample clustering) + smear (MA) plots
#   where <tag> encodes mark + contrast, <norm> is 'libsize' or 'spikein'.
# =============================================================================

suppressWarnings(suppressMessages({
  library(optparse)
}))

# -----------------------------------------------------------------------------
# Command-line interface
# -----------------------------------------------------------------------------
parse_cli <- function() {
  option_list <- list(
    make_option("--sheet",   type = "character",
                help = "Sample sheet CSV from build_sample_sheet.R"),
    make_option("--mark",    type = "character",
                help = "Broad Factor to analyze: H3K27me3 | H3K9me3 | H3K4me1"),
    make_option("--group1",  type = "character", default = "B10|none",
                help = "Numerator group as Condition|Treatment [default %default]"),
    make_option("--group2",  type = "character", default = "Parental|none",
                help = "Denominator/reference group [default %default]"),
    make_option("--outdir",  type = "character", default = "csaw_out"),
    make_option("--window",  type = "integer",   default = 2000L,
                help = "Window width in bp; widen for broader marks [default %default]"),
    make_option("--spacing", type = "integer",   default = 500L,
                help = "Spacing between window starts in bp [default %default]"),
    make_option("--bg_bin",  type = "integer",   default = 10000L,
                help = "Background bin width in bp for the enrichment filter [default %default]"),
    make_option("--min_lfc_bg", type = "double", default = log2(3),
                help = "Min log2 enrichment over background to keep a window [default log2(3)]"),
    make_option("--frag_max", type = "integer",  default = 700L,
                help = "Max fragment length (match Bowtie2 --X) [default %default]"),
    make_option("--fdr",     type = "double",    default = 0.05,
                help = "FDR threshold for the significant subset [default %default]"),
    make_option("--lfc",     type = "double",    default = 1.0,
                help = "Min |log2FC| for the significant subset [default %default]"),
    make_option("--blacklist", type = "character", default = NULL,
                help = "Optional ENCODE blacklist BED of regions to discard"),
    make_option("--spikein", type = "character", default = NULL,
                help = "Optional CSV (SampleID,factor) of spike-in scale factors")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$sheet)) stop("--sheet is required.")
  if (is.null(opt$mark))  stop("--mark is required (e.g. H3K27me3).")
  opt
}

# -----------------------------------------------------------------------------
# Parse a "Condition|Treatment" group spec into its two lowercased parts.
# Lowercasing here is what makes contrast matching case-insensitive.
# -----------------------------------------------------------------------------
parse_group <- function(spec) {
  parts <- strsplit(spec, "\\|")[[1]]
  if (length(parts) != 2)
    stop("Group spec must be 'Condition|Treatment', got: ", spec)
  list(cond = tolower(trimws(parts[1])),
       treat = tolower(trimws(parts[2])))
}

# -----------------------------------------------------------------------------
# Load the sheet and subset to the two contrast groups for the requested mark.
# All comparisons are done on lowercased values so the sheet's casing and the
# flag's casing never need to agree.
# -----------------------------------------------------------------------------
load_contrast_samples <- function(sheet_path, mark, g1, g2) {
  full <- read.csv(sheet_path, stringsAsFactors = FALSE)

  # Case-insensitive Factor match.
  ss <- full[tolower(full$Factor) == tolower(mark), ]
  if (nrow(ss) == 0) stop("No samples with Factor == ", mark, " in the sheet.")

  # Case-insensitive membership test for each group.
  in_group <- function(rows, g)
    tolower(rows$Condition) == g$cond & tolower(rows$Treatment) == g$treat

  is_g1 <- in_group(ss, g1)
  is_g2 <- in_group(ss, g2)
  ss <- ss[is_g1 | is_g2, ]
  if (nrow(ss) == 0)
    stop("No samples matched groups for ", mark,
         ". Check that --group1/--group2 Condition|Treatment exist in the sheet.")

  # Label each retained row g1 (numerator) or g2 (reference) and sort.
  ss$grp <- ifelse(in_group(ss, g1), "g1", "g2")
  ss <- ss[order(ss$grp, ss$Replicate), ]

  # Guard against the rank-deficiency error: need >=2 replicates per group.
  tab <- table(ss$grp)
  message("[INFO] ", mark, ": ", nrow(ss), " samples (",
          "g1=", sum(ss$grp == "g1"), ", g2=", sum(ss$grp == "g2"), ")")
  if (any(tab < 2))
    stop("A group has <2 replicates (", paste(names(tab), tab, collapse=", "),
         "). The NB model cannot estimate dispersion. Check your --group1/--group2 ",
         "casing and that both conditions exist for this mark.")
  ss
}

# -----------------------------------------------------------------------------
# Count reads into genome-wide windows, then filter to windows enriched over a
# local large-bin background. Returns list(win, bins) of RangedSummarizedExp.
# -----------------------------------------------------------------------------
count_and_filter_windows <- function(bam_files, sample_ids, opt) {
  library(csaw)
  param <- readParam(pe = "both", max.frag = opt$frag_max)
  if (!is.null(opt$blacklist)) {
    library(rtracklayer)
    param <- reform(param, discard = import(opt$blacklist))
    message("[INFO] Blacklist applied: ", opt$blacklist)
  }

  message("[INFO] Counting windows (width=", opt$window,
          ", spacing=", opt$spacing, ")...")
  win <- windowCounts(bam_files, width = opt$window,
                      spacing = opt$spacing, param = param)
  colnames(win) <- sample_ids

  message("[INFO] Counting background bins (", opt$bg_bin, " bp) and filtering...")
  bins <- windowCounts(bam_files, bin = TRUE, width = opt$bg_bin, param = param)
  keep <- filterWindowsGlobal(win, bins)$filter > opt$min_lfc_bg
  win <- win[keep, ]
  message("[INFO] Windows retained after background filter: ",
          sum(keep), " of ", length(keep))
  list(win = win, bins = bins)
}

# -----------------------------------------------------------------------------
# Build the model design. Adds a Batch blocking term only if >1 batch is present
# among the contrast samples; otherwise a simple ~group design.
# group is releveled so g2 is the reference (=> last coef tests g1 vs g2).
# -----------------------------------------------------------------------------
build_design <- function(ss) {
  group <- factor(ss$grp, levels = c("g2", "g1"))   # g2 = reference
  batch <- factor(ss$Batch)
  if (nlevels(batch) > 1) {
    message("[INFO] >1 batch among contrast samples -> design = ~Batch + group.")
    list(design = model.matrix(~ batch + group), group = group)
  } else {
    list(design = model.matrix(~ group), group = group)
  }
}

# -----------------------------------------------------------------------------
# Run one normalization scheme end-to-end: normalize -> estimate dispersion ->
# QL F-test -> merge windows into regions -> combine per-region stats -> save
# the full table, the significant subset, and QC plots.
# -----------------------------------------------------------------------------
run_one_normalization <- function(win, bins, ss, design_info, norm_tag,
                                   opt, tag, spikein_df = NULL) {
  library(csaw); library(edgeR)
  message("[INFO] === Normalization: ", norm_tag, " ===")

  y      <- asDGEList(win)
  design <- design_info$design
  group  <- design_info$group

  # Normalization factors: spike-in (user-supplied) or TMM-on-background.
  if (!is.null(spikein_df)) {
    facs <- spikein_df$factor[match(ss$SampleID, spikein_df$SampleID)]
    if (any(is.na(facs)))
      stop("Spike-in file missing factors for: ",
           paste(ss$SampleID[is.na(facs)], collapse = ", "))
    y$samples$norm.factors <- facs / exp(mean(log(facs)))  # scale to geo-mean 1
  } else {
    y$samples$norm.factors <- normFactors(bins, se.out = TRUE)$norm.factors
  }

  # Negative-binomial fit + quasi-likelihood F-test on the group coefficient.
  y      <- estimateDisp(y, design)
  fit    <- glmQLFit(y, design, robust = TRUE)
  result <- glmQLFTest(fit, coef = ncol(design))   # last coef = g1 vs g2

  # Merge adjacent significant windows into regions; combine their statistics.
  merged <- mergeWindows(rowRanges(win), tol = opt$spacing, max.width = 50000)
  tab    <- combineTests(merged$ids, result$table)
  out    <- cbind(as.data.frame(granges(merged$regions)), as.data.frame(tab))

  base <- file.path(opt$outdir, paste0(tag, "_", norm_tag))
  write.csv(out, paste0(base, "_all_regions.csv"), row.names = FALSE)

  # Significant subset: FDR threshold, plus a |log2FC| threshold if a logFC
  # column is available from combineTests.
  lfc_col <- intersect(c("rep.logFC", "best.logFC"), colnames(out))
  lfc_col <- if (length(lfc_col)) lfc_col[1] else NA
  sig <- out[!is.na(out$FDR) & out$FDR < opt$fdr, ]
  if (!is.na(lfc_col)) sig <- sig[abs(sig[[lfc_col]]) >= opt$lfc, ]
  write.csv(sig, paste0(base, "_significant.csv"), row.names = FALSE)
  message("[INFO] ", norm_tag, ": ", nrow(sig), " significant regions of ", nrow(out),
          " (FDR<", opt$fdr,
          if (!is.na(lfc_col)) paste0(" & |log2FC|>=", opt$lfc) else "", ")")

  # Direction balance among significant regions (gained vs lost in group1).
  # A wildly lopsided split can indicate a global shift that library-size
  # normalization mishandles — worth noting in the log automatically.
  if (!is.na(lfc_col) && nrow(sig) > 0) {
    up   <- sum(sig[[lfc_col]] > 0)
    down <- sum(sig[[lfc_col]] < 0)
    message(sprintf("[INFO] %s: direction balance — up(in g1)=%d (%.0f%%), down=%d (%.0f%%)",
                    norm_tag, up, 100*up/(up+down), down, 100*down/(up+down)))
  }

  # QC plots for the supplement (PDF for the paper + PNG for fast preview).
  save_qc_plots(y, result, group, ss, opt, base, tag, norm_tag)
  invisible(out)
}

# -----------------------------------------------------------------------------
# render_plot(): draw the SAME plot to BOTH a PDF (vector, for the paper) and a
# PNG (raster, for fast on-screen preview). `draw_fn` is a zero-arg function
# that issues the plotting calls; we just open/close the two devices around it.
# -----------------------------------------------------------------------------
render_plot <- function(draw_fn, base_path, w = 7, h = 6) {
  # PDF
  tryCatch({ pdf(paste0(base_path, ".pdf"), width = w, height = h); draw_fn() },
           error = function(e) message("  [QC] PDF skipped: ", conditionMessage(e)),
           finally = if (!is.null(dev.list())) dev.off())
  # PNG (150 dpi is plenty for a quick look)
  tryCatch({ png(paste0(base_path, ".png"), width = w*150, height = h*150, res = 150); draw_fn() },
           error = function(e) message("  [QC] PNG skipped: ", conditionMessage(e)),
           finally = if (!is.null(dev.list())) dev.off())
}

# -----------------------------------------------------------------------------
# QC plots, each written as both .pdf and .png:
#   *_MDS  : sample clustering. Points coloured by CONDITION with a legend, and
#            short rep labels offset below each point — so replicate labels no
#            longer overprint each other (the original bug). Lets you SEE whether
#            replicates of a condition cluster together (good) or scatter (a
#            reproducibility concern, as seen for Parental H3K27me3).
#   *_MA   : smear/MA — log2FC vs average abundance, with +/- lfc guide lines.
# -----------------------------------------------------------------------------
save_qc_plots <- function(y, result, group, ss, opt, base, tag, norm_tag) {
  library(edgeR)

  # Colour by CONDITION (not the g1/g2 grouping) so replicates of the same
  # biological condition share a colour and you can judge their concordance.
  cond     <- factor(ss$Condition)
  cond_col <- as.integer(cond)
  # Self-documenting labels: condition + replicate (e.g. "parental R1"), so the
  # plot is readable without cross-referencing the colour legend.
  full_lab <- paste(ss$Condition, paste0("R", ss$Replicate))

  render_plot(function() {
    op <- par(mar = c(5, 5, 3, 9), xpd = TRUE); on.exit(par(op))
    mds <- plotMDS(cpm(y, log = TRUE), plot = FALSE)   # compute coords only
    # tiny jitter so coincident replicate points (e.g. concordant B10) don't
    # stack their labels on top of each other.
    jx <- jitter(mds$x, amount = diff(range(mds$x)) * 0.01)
    plot(jx, mds$y, col = cond_col, pch = 19, cex = 1.6,
         xlab = "Leading logFC dim 1", ylab = "Leading logFC dim 2",
         main = paste(tag, norm_tag, "MDS"))
    text(jx, mds$y, labels = full_lab, pos = 1, cex = 0.7)
    legend("topright", inset = c(-0.32, 0), legend = levels(cond),
           col = seq_along(levels(cond)), pch = 19, bty = "n", title = "Condition")
  }, paste0(base, "_MDS"))

  render_plot(function() {
    plotSmear(result, main = paste(tag, norm_tag, "smear (MA)"))
    abline(h = c(-opt$lfc, opt$lfc), col = "blue", lty = 2)
  }, paste0(base, "_MA"))
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
main <- function() {
  opt <- parse_cli()
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

  g1 <- parse_group(opt$group1)
  g2 <- parse_group(opt$group2)
  ss <- load_contrast_samples(opt$sheet, opt$mark, g1, g2)

  # Human-readable tag for output filenames, e.g. H3K27me3_B10-none_vs_Parental-none
  tag <- paste0(opt$mark, "_",
                gsub("[|]", "-", opt$group1), "_vs_",
                gsub("[|]", "-", opt$group2))

  counts      <- count_and_filter_windows(ss$bamReads, ss$SampleID, opt)
  design_info <- build_design(ss)

  # Always run library-size/TMM normalization.
  run_one_normalization(counts$win, counts$bins, ss, design_info,
                        "libsize", opt, tag, spikein_df = NULL)

  # Optionally also run spike-in normalization.
  if (!is.null(opt$spikein)) {
    sf <- read.csv(opt$spikein, stringsAsFactors = FALSE)
    if (!all(c("SampleID", "factor") %in% colnames(sf)))
      stop("--spikein CSV must have columns: SampleID, factor")
    run_one_normalization(counts$win, counts$bins, ss, design_info,
                          "spikein", opt, tag, spikein_df = sf)
  } else {
    message("[NOTE] --spikein not supplied; library-size only. ",
            "Global-change results for broad repressive marks are provisional.")
  }

  message("[INFO] csaw complete: ", tag, " -> ", opt$outdir)
}

main()