#!/usr/bin/env Rscript
# =============================================================================
# qc_samples.R — Per-sample CUT&RUN QC in ONE pass: depth + duplication + FRiP
#                (optionally) figures. Replaces qc_depth_dup.R + qc_frip.R +
#                qc_plots.R with one flag-driven script sharing one input.
# =============================================================================
#
# PROJECT: A375 melanoma parental-vs-B10 constricted-migration multi-omics.
#
# WHY ONE SCRIPT --------------------------------------------------------------
#   depth/duplication and FRiP are all PER-SAMPLE metrics keyed on the same
#   identity columns (SampleID/Factor/Condition/Treatment/Replicate). Computing
#   them together yields ONE non-redundant table (no post-hoc merge), and the
#   plots can be drawn directly from it. (The genome-wide CLUSTERING/PCA QC is
#   intentionally kept in its own script, qc_sample_rep_global.R, because it is
#   much heavier compute and answers a different question — "do samples cluster
#   sanely as a group" vs "is each individual sample good quality".)
#
# METRICS (all computed by default) -------------------------------------------
#   depth : total/mapped/properly-paired reads + mapped% (from `samtools flagstat`)
#   dup   : Picard PERCENT_DUPLICATION (from metrics.*.txt; bundled with depth
#           since both describe library quality and are cheap)
#   frip  : Fraction of Reads in Peaks (reads-in-peaks / total; from csaw counts)
#   --metrics is OPTIONAL. With NO flag, ALL metrics run. If --metrics is given,
#   at least one of {depth,frip} must be named (dup always rides along with
#   depth). Example: --metrics depth   (fast quick-look, skips the slow FRiP).
#
# WHAT EACH METRIC DIAGNOSES --------------------------------------------------
#   Low depth / high duplication  -> weak/over-amplified library (noise source).
#   Low FRiP                      -> weak antibody enrichment (IP didn't work).
#   Good depth+dup+FRiP but replicates still scramble on MDS -> suspect a SWAP /
#     batch effect rather than quality (chase provenance, e.g. BottomR5-2).
#   Outlier flags are computed WITHIN each mark (broad marks naturally have
#   lower FRiP than punctate marks, so only like-with-like comparison is valid).
#
# OPTIONAL FIGURES (--plot) ---------------------------------------------------
#   (a) frip_by_mark      FRiP bars, faceted per mark, coloured by condition,
#                         weak-FRiP bars outlined. The headline QC figure.
#   (b) qc_three_metrics  stacked FRiP / depth / duplication dot panels sharing
#                         a sample x-axis, so a sample weak on MULTIPLE metrics
#                         (e.g. low depth AND low FRiP) is obvious at a glance.
#   Figures need FRiP+depth, so --plot implies all metrics. ggplot2 is used if
#   present, else a base-R fallback. Each figure is written as PDF + PNG.
#
# EXAMPLE USAGE ---------------------------------------------------------------
#   # everything (all metrics + figures) — the usual call:
#   Rscript qc_samples.R --sheet data/processed/sample_sheet.csv \
#       --outdir results/exploratory/qc_samples --plot
#
#   # fast quality-only look, no FRiP, no plots:
#   Rscript qc_samples.R --sheet ss.csv --outdir qc/ --metrics depth
#
#   # one mark while investigating it:
#   Rscript qc_samples.R --sheet ss.csv --outdir qc/ --mark H3K27ac --plot
#
# OUTPUTS (into --outdir) -----------------------------------------------------
#   qc_samples.csv      ONE per-sample table with all requested metrics + flags
#   frip_by_mark.*      (if --plot) faceted FRiP bar figure (pdf+png)
#   qc_three_metrics.*  (if --plot) stacked three-metric dot figure (pdf+png)
# =============================================================================

suppressWarnings(suppressMessages(library(optparse)))

# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
parse_cli <- function() {
  ol <- list(
    make_option("--sheet", type = "character",
                help = "Sample sheet CSV from build_sample_sheet.R"),
    make_option("--outdir", type = "character", default = "qc_samples",
                help = "Output directory [default %default]"),
    make_option("--metrics", type = "character", default = NULL,
                help = paste("Comma-sep subset of {depth,frip} to compute.",
                             "Omit for ALL. dup always accompanies depth.",
                             "[default: all]")),
    make_option("--mark", type = "character", default = NULL,
                help = "Optional: restrict to one Factor (e.g. H3K27ac)"),
    make_option("--plot", action = "store_true", default = FALSE,
                help = "Also draw QC figures (implies all metrics)"),
    make_option("--samtools", type = "character", default = "samtools",
                help = "samtools binary [default %default]"),
    make_option("--frag_max", type = "integer", default = 700L,
                help = "Max fragment length, match Bowtie2 --X [default %default]"),
    # Within-mark outlier thresholds (consistent with the prior separate scripts)
    make_option("--depth_frac", type = "double", default = 0.5,
                help = "Flag depth < this * within-mark median [default %default]"),
    make_option("--dup_excess", type = "double", default = 15.0,
                help = "Flag dup_pct > within-mark median + this pts [default %default]"),
    make_option("--frip_frac", type = "double", default = 0.5,
                help = "Flag FRiP < this * within-mark median [default %default]")
  )
  opt <- parse_args(OptionParser(option_list = ol))
  if (is.null(opt$sheet)) stop("--sheet is required.")
  opt
}

# -----------------------------------------------------------------------------
# Decide which metrics to run from the (optional) --metrics flag and --plot.
# Rules: no flag -> all. With flag -> must name >=1 of {depth,frip}; unknown
# tokens error out. --plot forces all (figures need both depth and FRiP).
# 'dup' is not separately selectable: it always rides along with depth.
# -----------------------------------------------------------------------------
resolve_metrics <- function(opt) {
  all_metrics <- c("depth", "frip")
  if (opt$plot) {
    if (!is.null(opt$metrics))
      message("[INFO] --plot given: computing ALL metrics (figures need both).")
    return(all_metrics)
  }
  if (is.null(opt$metrics)) return(all_metrics)         # default = all
  req <- tolower(trimws(strsplit(opt$metrics, ",")[[1]]))
  req <- req[req != ""]
  bad <- setdiff(req, all_metrics)
  if (length(bad))
    stop("Unknown --metrics token(s): ", paste(bad, collapse = ", "),
         ". Valid: depth, frip (dup is bundled with depth).")
  if (length(req) == 0)
    stop("--metrics was given but named no valid metric. ",
         "Name at least one of: depth, frip.")
  req
}

# -----------------------------------------------------------------------------
# depth + duplication for one sample.
#   - flagstat: capture BOTH stdout+stderr, no shQuote, check exit status, so a
#     samtools failure surfaces as a warning instead of silent NA. Handles the
#     newer 'mapped (100.00% : N/A)' format.
#   - dup: read PERCENT_DUPLICATION from the Picard metrics.*.txt in aligned/
#     (one level above the dedup BAM); convert fraction -> percent.
# -----------------------------------------------------------------------------
flagstat_counts <- function(bam, samtools) {
  if (!file.exists(bam)) { warning("BAM not found: ", bam)
    return(list(total = NA, mapped = NA, mapped_pct = NA, pp = NA)) }
  out <- tryCatch(
    suppressWarnings(system2(samtools, c("flagstat", bam),
                             stdout = TRUE, stderr = TRUE)),
    error = function(e) { warning("samtools failed for ", bam, ": ",
                                  conditionMessage(e)); character(0) })
  status <- attr(out, "status")
  if (length(out) == 0 || (!is.null(status) && status != 0)) {
    warning("flagstat gave no usable output for ", bam,
            if (length(out)) paste0(" | ", paste(out, collapse = " ")) else "")
    return(list(total = NA, mapped = NA, mapped_pct = NA, pp = NA))
  }
  num1 <- function(line) as.numeric(sub("^([0-9]+).*", "\\1", line))
  total_line  <- out[1]
  mapped_line <- grep(" mapped \\(", out, value = TRUE)[1]
  pp_line     <- grep("properly paired", out, value = TRUE)[1]
  pct <- NA
  if (!is.na(mapped_line)) {
    mm <- regmatches(mapped_line, regexpr("\\(([0-9.]+)%", mapped_line))
    if (length(mm)) pct <- as.numeric(sub("\\(", "", sub("%", "", mm)))
  }
  list(total = num1(total_line),
       mapped = if (!is.na(mapped_line)) num1(mapped_line) else NA,
       mapped_pct = pct,
       pp = if (!is.na(pp_line)) num1(pp_line) else NA)
}

picard_dup_pct <- function(bam) {
  aligned_dir <- dirname(dirname(bam))           # .../aligned (above dedup/)
  if (!dir.exists(aligned_dir)) return(NA)
  mfiles <- list.files(aligned_dir, pattern = "^metrics.*\\.txt$", full.names = TRUE)
  if (length(mfiles) == 0) return(NA)
  lines <- tryCatch(readLines(mfiles[1], warn = FALSE), error = function(e) character(0))
  hdr_i <- grep("PERCENT_DUPLICATION", lines)
  if (length(hdr_i) == 0) return(NA)
  hdr <- strsplit(lines[hdr_i[1]], "\t")[[1]]
  col <- which(hdr == "PERCENT_DUPLICATION")
  dat <- strsplit(lines[hdr_i[1] + 1], "\t")[[1]]
  val <- suppressWarnings(as.numeric(dat[col]))
  if (is.na(val)) NA else round(100 * val, 2)
}

# -----------------------------------------------------------------------------
# FRiP for one sample: reads-in-peaks / total-reads, via csaw (paired-end).
# Reads a MACS narrow/broadPeak (BED-like) into a GRanges (BED is 0-based start).
# -----------------------------------------------------------------------------
read_peaks <- function(path) {
  if (is.na(path) || path == "" || !file.exists(path)) return(NULL)
  bed <- tryCatch(read.table(path, sep = "\t", stringsAsFactors = FALSE,
                             comment.char = "#", header = FALSE),
                  error = function(e) NULL)
  if (is.null(bed) || ncol(bed) < 3) return(NULL)
  GenomicRanges::GRanges(bed[[1]], IRanges::IRanges(bed[[2]] + 1, bed[[3]]))
}

frip_one <- function(bam, peaks_path, param) {
  peaks <- read_peaks(peaks_path)
  if (is.null(peaks)) { warning("No usable peaks: ", peaks_path)
    return(list(total = NA, inpk = NA, frip = NA, n = 0)) }
  total <- tryCatch(sum(SummarizedExperiment::assay(
                      csaw::windowCounts(bam, bin = TRUE, width = 10000, param = param))),
                    error = function(e) NA)
  inpk  <- tryCatch(sum(SummarizedExperiment::assay(
                      csaw::regionCounts(bam, regions = peaks, param = param))),
                    error = function(e) NA)
  frip  <- if (!is.na(total) && total > 0 && !is.na(inpk)) round(inpk / total, 4) else NA
  list(total = total, inpk = inpk, frip = frip, n = length(peaks))
}

# -----------------------------------------------------------------------------
# Build the single per-sample table for the requested metrics.
# -----------------------------------------------------------------------------
build_table <- function(ss, metrics, opt) {
  need_frip <- "frip" %in% metrics
  if (need_frip)
    suppressWarnings(suppressMessages({
      library(csaw); library(GenomicRanges); library(SummarizedExperiment) }))
  param <- if (need_frip) csaw::readParam(pe = "both", max.frag = opt$frag_max) else NULL

  rows <- vector("list", nrow(ss))
  for (i in seq_len(nrow(ss))) {
    r <- ss[i, ]
    rec <- data.frame(SampleID = r$SampleID, Factor = r$Factor,
                      Condition = r$Condition, Treatment = r$Treatment,
                      Replicate = r$Replicate, stringsAsFactors = FALSE)
    if ("depth" %in% metrics) {
      fs <- flagstat_counts(r$bamReads, opt$samtools)
      rec$total_reads     <- fs$total
      rec$mapped_reads    <- fs$mapped
      rec$mapped_pct      <- fs$mapped_pct
      rec$properly_paired <- fs$pp
      rec$dup_pct         <- picard_dup_pct(r$bamReads)
    }
    if (need_frip) {
      fr <- frip_one(r$bamReads, r$Peaks, param)
      rec$frip_total_reads <- fr$total
      rec$reads_in_peaks   <- fr$inpk
      rec$FRiP             <- fr$frip
      rec$n_peaks          <- fr$n
    }
    msg <- paste0("  ", r$SampleID)
    if ("depth" %in% metrics) msg <- paste0(msg, " | dup%=", rec$dup_pct)
    if (need_frip)            msg <- paste0(msg, " | FRiP=", rec$FRiP)
    message(msg)
    rows[[i]] <- rec
  }
  do.call(rbind, rows)
}

# -----------------------------------------------------------------------------
# Add WITHIN-MARK outlier flags for whichever metrics are present.
# -----------------------------------------------------------------------------
add_flags <- function(tab, opt) {
  has_depth <- "mapped_reads" %in% colnames(tab)
  has_frip  <- "FRiP" %in% colnames(tab)
  if (has_depth) { tab$depth_flag <- ""; tab$dup_flag <- "" }
  if (has_frip)  tab$frip_flag <- ""
  for (mk in unique(tab$Factor)) {
    idx <- which(tab$Factor == mk)
    if (has_depth) {
      md <- median(tab$mapped_reads[idx], na.rm = TRUE)
      mu <- median(tab$dup_pct[idx], na.rm = TRUE)
      tab$depth_flag[idx[!is.na(tab$mapped_reads[idx]) &
                         tab$mapped_reads[idx] < opt$depth_frac * md]] <- "LOW_DEPTH"
      tab$dup_flag[idx[!is.na(tab$dup_pct[idx]) &
                       tab$dup_pct[idx] > mu + opt$dup_excess]] <- "HIGH_DUP"
    }
    if (has_frip) {
      mf <- median(tab$FRiP[idx], na.rm = TRUE)
      tab$frip_flag[idx[!is.na(tab$FRiP[idx]) &
                        tab$FRiP[idx] < opt$frip_frac * mf]] <- "LOW_FRIP"
    }
  }
  tab
}

# -----------------------------------------------------------------------------
# render to BOTH pdf and png around one drawing function.
# -----------------------------------------------------------------------------
render_plot <- function(draw_fn, base_path, w = 10, h = 7) {
  tryCatch({ pdf(paste0(base_path, ".pdf"), width = w, height = h); draw_fn() },
           error = function(e) message("  [plot] PDF skipped: ", conditionMessage(e)),
           finally = if (!is.null(dev.list())) dev.off())
  tryCatch({ png(paste0(base_path, ".png"), width = w*150, height = h*150, res = 150); draw_fn() },
           error = function(e) message("  [plot] PNG skipped: ", conditionMessage(e)),
           finally = if (!is.null(dev.list())) dev.off())
}

# -----------------------------------------------------------------------------
# Figures (ggplot2 preferred, base-R fallback). Needs FRiP + depth + dup.
# -----------------------------------------------------------------------------
make_figures <- function(tab, opt) {
  # Helper columns for plotting.
  tab$depth_M <- tab$mapped_reads / 1e6
  tab$grp <- ifelse(tab$Treatment %in% c("none", "", NA),
                    tab$Condition, paste(tab$Condition, tab$Treatment, sep = "."))
  tab$xlab <- paste0(tab$grp, " R", tab$Replicate)
  tab <- tab[order(tab$Factor, tab$grp, tab$Replicate), ]
  tab$frip_low <- !is.na(tab$frip_flag) & tab$frip_flag == "LOW_FRIP"

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    suppressWarnings(suppressMessages(library(ggplot2)))
    cond_cols <- c("B10" = "#1b9e77", "parental" = "#d95f02")

    pa <- ggplot(tab, aes(xlab, FRiP, fill = Condition)) +
      geom_col(aes(colour = frip_low), linewidth = 0.6) +
      scale_colour_manual(values = c("FALSE" = NA, "TRUE" = "red"), guide = "none") +
      scale_fill_manual(values = cond_cols, name = "Condition") +
      facet_wrap(~ Factor, scales = "free_x") +
      labs(title = "FRiP by sample, faceted by mark",
           subtitle = "Red outline = FRiP < 50% of within-mark median (weak enrichment)",
           x = NULL, y = "FRiP (fraction of reads in peaks)") +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            panel.grid.major.x = element_blank())
    render_plot(function() print(pa), file.path(opt$outdir, "frip_by_mark"))

    long <- rbind(
      data.frame(tab[, c("xlab","Factor","Condition")], metric = "FRiP",            value = tab$FRiP),
      data.frame(tab[, c("xlab","Factor","Condition")], metric = "Depth (M reads)", value = tab$depth_M),
      data.frame(tab[, c("xlab","Factor","Condition")], metric = "Duplication (%)", value = tab$dup_pct))
    long$metric <- factor(long$metric,
                          levels = c("FRiP", "Depth (M reads)", "Duplication (%)"))
    long$sample <- factor(long$xlab, levels = unique(tab$xlab))
    pb <- ggplot(long, aes(sample, value, colour = Factor, shape = Condition)) +
      geom_point(size = 3) +
      facet_grid(metric ~ ., scales = "free_y", switch = "y") +
      labs(title = "Per-sample QC: FRiP, depth, duplication",
           subtitle = "Aligned panels — a sample weak on multiple metrics is a concern",
           x = NULL, y = NULL, colour = "Mark", shape = "Condition") +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
            strip.placement = "outside", panel.grid.major.x = element_blank())
    render_plot(function() print(pb), file.path(opt$outdir, "qc_three_metrics"),
                w = 12, h = 9)

  } else {
    message("[INFO] ggplot2 not found — base-R fallback figures.")
    marks <- unique(tab$Factor)
    cc <- function(cond) ifelse(cond == "B10", "#1b9e77", "#d95f02")
    render_plot(function() {
      op <- par(mfrow = c(1, length(marks)), mar = c(9, 4, 3, 1)); on.exit(par(op))
      for (mk in marks) { sub <- tab[tab$Factor == mk, ]
        barplot(sub$FRiP, names.arg = sub$xlab, las = 2, col = cc(sub$Condition),
                border = ifelse(sub$frip_low, "red", NA), main = mk,
                ylab = "FRiP", cex.names = 0.7) }
    }, file.path(opt$outdir, "frip_by_mark"), w = 13, h = 6)
    render_plot(function() {
      op <- par(mfrow = c(3, 1), mar = c(0.5, 5, 2, 1), oma = c(10, 0, 2, 0))
      on.exit(par(op)); x <- seq_len(nrow(tab)); cols <- as.integer(factor(tab$Factor))
      plot(x, tab$FRiP, col = cols, pch = 19, xaxt = "n", ylab = "FRiP",
           main = "Per-sample QC (FRiP / depth / dup)")
      plot(x, tab$depth_M, col = cols, pch = 19, xaxt = "n", ylab = "Depth (M)")
      plot(x, tab$dup_pct, col = cols, pch = 19, xaxt = "n", ylab = "Dup (%)")
      axis(1, at = x, labels = tab$xlab, las = 2, cex.axis = 0.7)
    }, file.path(opt$outdir, "qc_three_metrics"), w = 12, h = 9)
  }
}

# -----------------------------------------------------------------------------
# Pretty-print the relevant columns and any flagged outliers.
# -----------------------------------------------------------------------------
report <- function(tab) {
  cols <- c("SampleID", "Factor", "Condition", "Replicate")
  if ("mapped_reads" %in% colnames(tab)) cols <- c(cols, "mapped_reads", "dup_pct", "depth_flag", "dup_flag")
  if ("FRiP" %in% colnames(tab))         cols <- c(cols, "FRiP", "n_peaks", "frip_flag")
  message("\n[INFO] QC table:")
  print(tab[order(tab$Factor, tab$Condition, tab$Replicate), cols], row.names = FALSE)

  flag_cols <- intersect(c("depth_flag","dup_flag","frip_flag"), colnames(tab))
  if (length(flag_cols)) {
    flagged <- tab[Reduce(`|`, lapply(flag_cols, function(c) tab[[c]] != "")), ]
    if (nrow(flagged) > 0) {
      message("\n[ACTION] Within-mark outliers:")
      print(flagged[, c("SampleID", "Factor", flag_cols)], row.names = FALSE)
    } else message("\n[INFO] No within-mark outliers flagged.")
  }
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
main <- function() {
  opt <- parse_cli()
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
  metrics <- resolve_metrics(opt)
  message("[INFO] Metrics: ", paste(metrics, collapse = ", "),
          if ("depth" %in% metrics) " (+dup)" else "")

  ss <- read.csv(opt$sheet, stringsAsFactors = FALSE)
  if (!is.null(opt$mark)) ss <- ss[tolower(ss$Factor) == tolower(opt$mark), ]
  if (nrow(ss) == 0) stop("No samples (check --mark).")
  message("[INFO] Samples: ", nrow(ss))

  tab <- build_table(ss, metrics, opt)
  tab <- add_flags(tab, opt)

  out_csv <- file.path(opt$outdir, "qc_samples.csv")
  write.csv(tab, out_csv, row.names = FALSE)
  report(tab)

  if (opt$plot) {
    message("\n[INFO] Drawing figures...")
    make_figures(tab, opt)
  }
  message("\n[INFO] -> ", out_csv,
          if (opt$plot) paste0(" (+ figures in ", opt$outdir, ")") else "")
}
main()