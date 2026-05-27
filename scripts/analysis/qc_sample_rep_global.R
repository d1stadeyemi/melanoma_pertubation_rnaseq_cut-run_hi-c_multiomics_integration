# =============================================================================
# qc_sample_rep_global.R — Global, all-sample QC (run BEFORE any differential contrast)
# =============================================================================
#
# PROJECT: A375 melanoma parental-vs-B10 constricted-migration multi-omics.
#
# WHAT THIS IS FOR ------------------------------------------------------------
#   Per-contrast MDS plots only show 2 groups of 1 mark at a time, and only
#   after a contrast has run. This script instead looks at ALL samples across
#   ALL marks AT ONCE, summarizes their genome-wide signal into a single matrix,
#   and produces:
#     (1) a sample-sample CORRELATION heatmap, and
#     (2) a PCA,
#   so you can judge global structure in one figure. This is the standard
#   "do my replicates and samples look sane?" supplementary QC (the R/edgeR
#   analogue of deepTools multiBigwigSummary + plotCorrelation/plotPCA).
#
# TWO THINGS IT CATCHES -------------------------------------------------------
#   (a) REPLICATE CONCORDANCE: replicates of a given mark+condition should sit
#       together; ones that scatter (e.g. the Parental H3K27me3 outlier seen in
#       the per-contrast MDS) stand out here in full context.
#   (b) PROVENANCE / MISLABELS: samples should cluster FIRST by mark. If an
#       ambiguous file (e.g. 'BottomR5-2-H3K27Ac', the 'Unmigrated' BAMs) does
#       NOT cluster with its labelled mark, that is a red flag that the file may
#       not be what its folder name claims — exactly the files we flagged in the
#       manifest. Use this plot to sanity-check them before trusting any contrast.
#
# WHY COUNT FROM BAMs (not the bigWigs) --------------------------------------
#   The sample sheet stores BAM paths, and counting from the SAME BAMs the
#   differential scripts use guarantees the QC reflects the actual analysis
#   inputs (including the real, possibly-mismatched filenames resolved via the
#   manifest). We count into large bins (default 10 kb) genome-wide, which is
#   ample resolution for a clustering/correlation overview and is fast.
#   (deepTools multiBigwigSummary on the .cpm.norm.bw files is an equally valid,
#    commonly-cited alternative; this R version is used here so it can read the
#    sheet directly and include the ambiguous files automatically.)
#
# EXAMPLE USAGE ---------------------------------------------------------------
#   Rscript qc_sample_rep_global.R \
#       --sheet data/processed/sample_sheet.csv \
#       --outdir results/exploratory/qc_sample_rep_global
#
#   # only the broad marks, bigger bins, Spearman correlation:
#   Rscript qc_sample_rep_global.R --sheet ss.csv --outdir qc/ \
#       --marks H3K27me3,H3K9me3,H3K4me1 --bin 25000 --cor_method spearman
#
# OUTPUTS (into --outdir, each as .pdf AND .png) ------------------------------
#   global_correlation   sample-sample correlation heatmap (clustered)
#   global_pca           PCA scatter (PC1 vs PC2), coloured by mark, shaped by condition
#   global_counts.csv    the bin x sample log-CPM matrix actually used (for re-use)
# =============================================================================

suppressWarnings(suppressMessages({
  library(optparse)
}))

# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
parse_cli <- function() {
  option_list <- list(
    make_option("--sheet",  type = "character",
                help = "Sample sheet CSV from build_sample_sheet.R"),
    make_option("--outdir", type = "character", default = "qc_sample_rep_global",
                help = "Output directory for QC plots and matrices [default %default]"),
    make_option("--marks",  type = "character", default = NULL,
                help = "Comma-sep Factors to include (default: all in the sheet)"),
    make_option("--bin",    type = "integer",   default = 10000L,
                help = "Genome bin width in bp for counting [default %default]"),
    make_option("--frag_max", type = "integer", default = 700L,
                help = "Max fragment length (match Bowtie2 --X) [default %default]"),
    make_option("--cor_method", type = "character", default = "spearman",
                help = "Correlation method: spearman (robust) | pearson [default %default]"),
    make_option("--min_count", type = "integer", default = 10L,
                help = "Drop bins whose total count across samples is below this [default %default]")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$sheet)) stop("--sheet is required.")
  opt
}

# -----------------------------------------------------------------------------
# Load the sheet; optionally restrict to certain marks. Returns the data.frame.
# -----------------------------------------------------------------------------
load_samples <- function(sheet_path, marks_csv) {
  ss <- read.csv(sheet_path, stringsAsFactors = FALSE)
  if (!is.null(marks_csv)) {
    want <- tolower(trimws(strsplit(marks_csv, ",")[[1]]))
    ss <- ss[tolower(ss$Factor) %in% want, ]
    if (nrow(ss) == 0) stop("No samples match --marks: ", marks_csv)
  }
  # A readable label per sample for plotting: Factor.Condition.Treatment.Rrep
  ss$label <- paste(ss$Factor, ss$Condition, ss$Treatment,
                    paste0("R", ss$Replicate), sep = ".")
  message("[INFO] Samples included: ", nrow(ss))
  print(table(ss$Factor, ss$Condition))
  ss
}

# -----------------------------------------------------------------------------
# Count all samples into genome-wide bins, filter very low bins, return a
# log2-CPM matrix (bins x samples). Uses csaw for paired-end binned counting.
# -----------------------------------------------------------------------------
build_logcpm_matrix <- function(ss, opt) {
  suppressWarnings(suppressMessages({ library(csaw); library(edgeR) }))
  param <- readParam(pe = "both", max.frag = opt$frag_max)

  message("[INFO] Counting all samples into ", opt$bin, " bp bins (paired-end)...")
  binned <- windowCounts(ss$bamReads, bin = TRUE, width = opt$bin, param = param)
  colnames(binned) <- ss$label

  # Drop near-empty bins so correlation/PCA are not dominated by zeros.
  cmat <- assay(binned)
  keep <- rowSums(cmat) >= opt$min_count
  cmat <- cmat[keep, , drop = FALSE]
  message("[INFO] Bins retained: ", nrow(cmat), " of ", length(keep))

  # log2-CPM (library-size normalized) for a fair sample-to-sample comparison.
  dge <- DGEList(counts = cmat)
  dge <- calcNormFactors(dge)            # TMM, robust to composition differences
  logcpm <- cpm(dge, log = TRUE, prior.count = 1)
  colnames(logcpm) <- ss$label
  logcpm
}

# -----------------------------------------------------------------------------
# render to BOTH pdf and png around one drawing function.
# -----------------------------------------------------------------------------
render_plot <- function(draw_fn, base_path, w = 8, h = 7) {
  tryCatch({ pdf(paste0(base_path, ".pdf"), width = w, height = h); draw_fn() },
           error = function(e) message("  [QC] PDF skipped: ", conditionMessage(e)),
           finally = if (!is.null(dev.list())) dev.off())
  tryCatch({ png(paste0(base_path, ".png"), width = w*150, height = h*150, res = 150); draw_fn() },
           error = function(e) message("  [QC] PNG skipped: ", conditionMessage(e)),
           finally = if (!is.null(dev.list())) dev.off())
}

# -----------------------------------------------------------------------------
# Correlation heatmap: pairwise sample correlation of the log-CPM matrix,
# clustered so similar samples sit together. Replicates of the same mark should
# form tight high-correlation blocks.
# -----------------------------------------------------------------------------
plot_correlation <- function(logcpm, opt, base) {
  cmat <- cor(logcpm, method = opt$cor_method)
  render_plot(function() {
    # Base-R heatmap with a clustering dendrogram; values annotated.
    op <- par(mar = c(10, 10, 3, 2))
    on.exit(par(op))
    n <- ncol(cmat)
    image(1:n, 1:n, cmat[, n:1], axes = FALSE, xlab = "", ylab = "",
          main = paste0("Sample correlation (", opt$cor_method, ")"),
          col = colorRampPalette(c("#2166ac", "white", "#b2182b"))(50),
          zlim = c(min(cmat), 1))
    axis(1, at = 1:n, labels = colnames(cmat), las = 2, cex.axis = 0.7)
    axis(2, at = 1:n, labels = rev(colnames(cmat)), las = 2, cex.axis = 0.7)
    for (i in 1:n) for (j in 1:n)
      text(i, n - j + 1, sprintf("%.2f", cmat[i, j]), cex = 0.6)
  }, base, w = 9, h = 9)
  # Also save the clustered heatmap if pheatmap is available (nicer dendrogram).
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    render_plot(function()
      pheatmap::pheatmap(cmat, main = paste0("Sample correlation (", opt$cor_method, ")"),
                         display_numbers = TRUE, number_format = "%.2f"),
      paste0(base, "_clustered"), w = 9, h = 8)
  }
  invisible(cmat)
}

# -----------------------------------------------------------------------------
# PCA on the log-CPM matrix (samples as observations). Colour by MARK (samples
# should separate by mark first), shape by CONDITION. A sample landing away from
# its mark's cluster is a provenance/mislabel red flag.
# -----------------------------------------------------------------------------
plot_pca <- function(logcpm, ss, base) {
  pca <- prcomp(t(logcpm), scale. = TRUE)
  ve  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)   # variance explained
  mark <- factor(ss$Factor)
  cond <- factor(ss$Condition)
  cols <- as.integer(mark)
  pch  <- c(16, 17, 15, 18, 1, 2)[as.integer(cond)]

  render_plot(function() {
    op <- par(mar = c(5, 5, 3, 9), xpd = TRUE); on.exit(par(op))
    plot(pca$x[, 1], pca$x[, 2], col = cols, pch = pch, cex = 1.8,
         xlab = paste0("PC1 (", ve[1], "%)"),
         ylab = paste0("PC2 (", ve[2], "%)"),
         main = "Global PCA — colour=mark, shape=condition")
    text(pca$x[, 1], pca$x[, 2], labels = paste0("R", ss$Replicate),
         pos = 1, cex = 0.7)
    legend("topright", inset = c(-0.25, 0), legend = levels(mark),
           col = seq_along(levels(mark)), pch = 16, bty = "n", title = "Mark")
    legend("bottomright", inset = c(-0.25, 0), legend = levels(cond),
           pch = c(16, 17, 15, 18)[seq_along(levels(cond))], bty = "n", title = "Condition")
  }, base, w = 9, h = 7)
  invisible(pca)
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
main <- function() {
  opt <- parse_cli()
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

  ss     <- load_samples(opt$sheet, opt$marks)
  logcpm <- build_logcpm_matrix(ss, opt)

  # Save the matrix so it can be reused without recounting.
  write.csv(data.frame(logcpm, check.names = FALSE),
            file.path(opt$outdir, "global_counts.csv"))

  plot_correlation(logcpm, opt, file.path(opt$outdir, "global_correlation"))
  plot_pca(logcpm, ss, file.path(opt$outdir, "global_pca"))

  message("[INFO] Global QC complete -> ", opt$outdir)
  message("[INFO] Check: do replicates of each mark cluster together, and do ",
          "the ambiguous files (BottomR5-2, *Unmigrated*) sit with their labelled mark?")
}

main()