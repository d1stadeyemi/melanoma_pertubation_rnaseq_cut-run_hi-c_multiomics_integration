# =============================================================================
# csaw.R — Window-based differential binding for BROAD histone marks
#          (including: interaction test, baseline stratification, IgG-control
#           filtering, replicate-correlation QC)
# =============================================================================
#
# PROJECT: A375 melanoma parental-vs-B10 constricted-migration multi-omics.
#
# WHAT THIS SCRIPT IS FOR ------------------------------------------------------
#   Broad chromatin marks (H3K27me3, H3K9me3, H3K4me1) form diffuse domains, not
#   sharp peaks. We TILE the genome into fixed windows and test each window for
#   differential coverage with the edgeR negative-binomial QL model; replicates
#   supply the within-group variance that turns a fold-change into an FDR-
#   controlled p-value. (Punctate H3K27Ac is handled separately by diffbind.R.)
#
# TWO TEST MODES (--test) -----------------------------------------------------
#   contrast    (default) Two-group comparison, exactly as before:
#                 --group1 NUMERATOR --group2 REFERENCE  (Condition|Treatment)
#                 Migration:     --group1 "B10|none"     --group2 "Parental|none"
#                 UNC in B10:    --group1 "B10|UNC"      --group2 "B10|none"
#                 UNC in parent: --group1 "Parental|UNC" --group2 "Parental|none"
#
#   interaction  Factorial Condition*Treatment model; tests the INTERACTION
#                 coefficient = "does UNC's effect on the mark DIFFER between B10
#                 and Parental?" This is the migration-acquired-PRC2-dependence
#                 set (Track 2). Requires all four cells (only H3K27me3 has them
#                 in this project). Positive interaction logFC = UNC raises the
#                 mark more (or lowers it less) in B10 than in Parental.
#                 Reference levels set by --cond_ref / --treat_ref.
#
#   Matching is CASE-INSENSITIVE throughout ("Parental"=="parental").
#
# BASELINE STRATIFICATION (the PRC2 "flattening" signature) -------------------
#   PRC2 inhibition does not lose the mark uniformly — it ERODES high-enrichment
#   peaks while signal spreads into low/background regions (Romero et al. 2024,
#   Fig 1g/h). A single genome-wide test lets those opposite moves cancel. So we
#   tag every region with its baseline (reference-group) abundance, split into
#   low/mid/high strata, and report the up/down direction balance per stratum.
#   You should see net LOSS in the high stratum and net GAIN in the low stratum
#   for the UNC contrasts; that bidirectional pattern IS the flattening.
#
# IgG-CONTROL FILTERING (optional, CUT&RUN-appropriate) -----------------------
#   CUT&RUN background is near-empty, so the genomic-background enrichment filter
#   behaves differently than in ChIP. If you have IgG for the relevant samples,
#   pass --igg bam1,bam2,... and windows are kept by enrichment over IgG
#   (filterWindowsControl) instead of over genomic background. Default (no --igg)
#   keeps the original filterWindowsGlobal behavior unchanged.
#   NOTE: this path could not be unit-tested locally; sanity-check the retained-
#   window count the first time you use it.
#
# SPIKE-IN (future-proof) -----------------------------------------------------
#   --spikein factors.csv (columns: SampleID,factor). When supplied, spike-in
#   factors REPLACE TMM and a second 'spikein' result set is written. Until then,
#   TMM/library-size runs and global-magnitude claims for repressive marks stay
#   PROVISIONAL (a true genome-wide loss can be masked).
#
# REQUIRED SAMPLE-SHEET COLUMNS ----------------------------------------------
#   SampleID   unique id per library (shown in QC plots; spike-in CSV key)
#   Factor     antibody/mark; matched to --mark (case-insensitive)
#   Condition  biological condition, e.g. Parental | B10
#   Treatment  drug arm, e.g. none | UNC
#   Replicate  1,2,...  (>=2 per group/cell, or the NB model cannot fit)
#   Batch      blocking covariate; if >1 level is present it is added to the design
#   bamReads   path to the coordinate-sorted, indexed BAM
#
# EXAMPLE USAGE ---------------------------------------------------------------
#   # Migration contrast (uses default groups B10|none vs Parental|none)
#   Rscript csaw.R --sheet ss.csv --mark H3K27me3 \
#       --outdir results/diff/H3K27me3_migration
#
#   # UNC reversal in B10   (B10|UNC vs B10|none)
#   Rscript csaw.R --sheet ss.csv --mark H3K27me3 \
#       --group1 "B10|UNC" --group2 "B10|none" \
#       --outdir results/diff/H3K27me3_B10unc
#
#   # UNC in parental = specificity control   (Parental|UNC vs Parental|none)
#   Rscript csaw.R --sheet ss.csv --mark H3K27me3 \
#       --group1 "Parental|UNC" --group2 "Parental|none" \
#       --outdir results/diff/H3K27me3_Parunc
#
#   # Interaction: does UNC's effect on H3K27me3 differ in B10 vs Parental?
#   Rscript csaw.R --sheet ss.csv --mark H3K27me3 --test interaction \
#       --outdir results/diff/H3K27me3_interaction
#
#   # Broad mark, wider windows   (H3K9me3 migration)
#   Rscript csaw.R --sheet ss.csv --mark H3K9me3 \
#       --window 5000 --spacing 1000 \
#       --outdir results/diff/H3K9me3_migration
#
#   # With IgG-based filtering and/or spike-in normalization (when available)
#   Rscript csaw.R --sheet ss.csv --mark H3K27me3 \
#       --igg igg_R1.bam,igg_R2.bam --spikein spikein_factors.csv \
#       --outdir results/diff/H3K27me3_migration_igg_spikein
#
# OUTPUTS (into --outdir) -----------------------------------------------------
#   <tag>_<norm>_all_regions.csv   merged regions: coords, logFC, FDR, direction,
#                                  baseline_logCPM, baseline_stratum
#   <tag>_<norm>_significant.csv   subset passing FDR & |log2FC|
#   <tag>_<norm>_{MDS,MA,volcano,repcor}.{pdf,png}
#   <norm> is 'libsize' or 'spikein'; <tag> encodes mark + contrast / interaction.
# =============================================================================

suppressWarnings(suppressMessages({
  library(optparse)
}))

# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
parse_cli <- function() {
  option_list <- list(
    # --- inputs (REQUIRED) ----------------------------------------------------
    make_option("--sheet", type = "character",
                help = "REQUIRED. Sample sheet CSV. Columns used: SampleID, Factor, Condition, Treatment, Replicate, Batch, bamReads."),
    make_option("--mark", type = "character",
                help = "REQUIRED. Factor/antibody to analyze, matched case-insensitively to the sheet's Factor column. e.g. H3K27me3 | H3K9me3 | H3K4me1."),

    # --- what to test ---------------------------------------------------------
    make_option("--test", type = "character", default = "contrast",
                help = "Test mode. ALLOWED: 'contrast' (two-group) | 'interaction' (Condition*Treatment; tests whether the treatment effect differs between conditions). [default %default]"),

    # --- contrast-mode groups (IGNORED when --test interaction) ---------------
    make_option("--group1", type = "character", default = "B10|none",
                help = "[contrast only] Numerator group as 'Condition|Treatment' (case-insensitive). e.g. 'B10|none', 'B10|UNC', 'Parental|UNC'. [default %default]"),
    make_option("--group2", type = "character", default = "Parental|none",
                help = "[contrast only] Reference/denominator group as 'Condition|Treatment'. Positive logFC = higher in group1. [default %default]"),

    # --- interaction-mode cells (IGNORED when --test contrast) ----------------
    make_option("--cond_ref", type = "character", default = "Parental",
                help = "[interaction only] Reference Condition; must exist in the sheet's Condition column. [default %default]"),
    make_option("--treat_ref", type = "character", default = "none",
                help = "[interaction only] Reference Treatment; must exist in the sheet's Treatment column. [default %default]"),
    make_option("--cond_alt", type = "character", default = "B10",
                help = "[interaction only] Non-reference Condition to include. [default %default]"),
    make_option("--treat_alt", type = "character", default = "UNC",
                help = "[interaction only] Non-reference Treatment to include. All four cells {cond_ref,cond_alt} x {treat_ref,treat_alt} must be present with >=2 reps each. [default %default]"),

    # --- output ---------------------------------------------------------------
    make_option("--outdir", type = "character", default = "csaw_out",
                help = "Output directory (created if absent). [default %default]"),

    # --- windowing / counting, all in BASE PAIRS ------------------------------
    make_option("--window", type = "integer", default = 2000L,
                help = "Window width (bp). Typical 1000-2000 for H3K27me3/H3K4me1; 2000-5000 for very broad H3K9me3. [default %default]"),
    make_option("--spacing", type = "integer", default = 500L,
                help = "Spacing between window starts (bp); usually window/4..window/2 (also the merge tolerance for adjacent windows). Smaller = finer & slower. [default %default]"),
    make_option("--bg_bin", type = "integer", default = 10000L,
                help = "Background bin width (bp) for TMM normalization and the genomic-background enrichment filter. 10000 is standard. [default %default]"),
    make_option("--frag_max", type = "integer", default = 700L,
                help = "Max fragment length (bp); should match your Bowtie2 -X/--maxins setting. Typically 600-1000. [default %default]"),

    # --- filtering & significance thresholds ----------------------------------
    make_option("--min_lfc_bg", type = "double", default = log2(3),
                help = "Min log2 enrichment over background (or IgG, if --igg) to KEEP a window. log2(3)=1.58 keeps >=3x; raise e.g. to log2(5)=2.32 to be stricter. [default log2(3)=1.58]"),
    make_option("--fdr", type = "double", default = 0.05,
                help = "FDR cutoff (range 0-1) for the significant subset. [default %default]"),
    make_option("--lfc", type = "double", default = 1.0,
                help = "Min |log2FC| for the significant subset (>=0). 1 = 2-fold; set 0 to threshold on FDR only. [default %default]"),
    make_option("--n_strata", type = "integer", default = 3L,
                help = "Number of baseline-abundance strata for the flattening report (>=2). 3 -> low/mid/high. [default %default]"),

    # --- optional controls / normalization (all default OFF) ------------------
    make_option("--blacklist", type = "character", default = NULL,
                help = "Optional ENCODE blacklist BED; reads overlapping these regions are discarded. [default: none]"),
    make_option("--igg", type = "character", default = NULL,
                help = "Optional comma-separated IgG BAM path(s). If set, windows are kept by enrichment over IgG (filterWindowsControl) instead of genomic background. [default: none]"),
    make_option("--spikein", type = "character", default = NULL,
                help = "Optional CSV with columns 'SampleID,factor'. If set, REPLACES TMM with spike-in normalization (factors centered to geometric mean 1) and writes a second 'spikein' result set. [default: none]")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$sheet)) stop("--sheet is required.")
  if (is.null(opt$mark))  stop("--mark is required (e.g. H3K27me3).")
  if (!opt$test %in% c("contrast", "interaction"))
    stop("--test must be 'contrast' or 'interaction'.")
  opt
}

# -----------------------------------------------------------------------------
# Parse a "Condition|Treatment" group spec into its two lowercased parts.
# -----------------------------------------------------------------------------
parse_group <- function(spec) {
  parts <- strsplit(spec, "\\|")[[1]]
  if (length(parts) != 2)
    stop("Group spec must be 'Condition|Treatment', got: ", spec)
  list(cond = tolower(trimws(parts[1])), treat = tolower(trimws(parts[2])))
}

# -----------------------------------------------------------------------------
# Sample selection.
#   contrast    : the two groups, labeled g1/g2 (as before).
#   interaction : the four cells {cond_ref,cond_alt} x {treat_ref,treat_alt}.
# Both return ss with the rows for this analysis, sorted, plus columns used
# downstream. A 'ref_cell' logical flags the baseline cell for stratification.
# -----------------------------------------------------------------------------
load_samples <- function(opt) {
  full <- read.csv(opt$sheet, stringsAsFactors = FALSE)
  ss <- full[tolower(full$Factor) == tolower(opt$mark), ]
  if (nrow(ss) == 0) stop("No samples with Factor == ", opt$mark, " in the sheet.")
  ss$Condition <- trimws(ss$Condition); ss$Treatment <- trimws(ss$Treatment)

  if (opt$test == "contrast") {
    g1 <- parse_group(opt$group1); g2 <- parse_group(opt$group2)
    in_group <- function(rows, g)
      tolower(rows$Condition) == g$cond & tolower(rows$Treatment) == g$treat
    ss <- ss[in_group(ss, g1) | in_group(ss, g2), ]
    if (nrow(ss) == 0)
      stop("No samples matched --group1/--group2 for ", opt$mark, ".")
    ss$grp      <- ifelse(in_group(ss, g1), "g1", "g2")
    ss$ref_cell <- ss$grp == "g2"                      # baseline = reference group
    ss <- ss[order(ss$grp, ss$Replicate), ]
    tab <- table(ss$grp)
    message("[INFO] ", opt$mark, " [contrast]: ", nrow(ss),
            " samples (g1=", sum(ss$grp == "g1"), ", g2=", sum(ss$grp == "g2"), ")")
    if (any(tab < 2))
      stop("A group has <2 replicates; NB dispersion is not estimable. ",
           "Check --group1/--group2 casing and that both cells exist for this mark.")
  } else {
    conds  <- tolower(c(opt$cond_ref,  opt$cond_alt))
    treats <- tolower(c(opt$treat_ref, opt$treat_alt))
    ss <- ss[tolower(ss$Condition) %in% conds & tolower(ss$Treatment) %in% treats, ]
    ss$ref_cell <- tolower(ss$Condition) == tolower(opt$cond_ref) &
                   tolower(ss$Treatment) == tolower(opt$treat_ref)
    ss <- ss[order(ss$Condition, ss$Treatment, ss$Replicate), ]
    cells <- table(paste(tolower(ss$Condition), tolower(ss$Treatment), sep = "|"))
    message("[INFO] ", opt$mark, " [interaction]: ", nrow(ss), " samples across cells:")
    print(cells)
    if (length(cells) < 4)
      stop("Interaction needs all four cells {", opt$cond_ref, ",", opt$cond_alt,
           "} x {", opt$treat_ref, ",", opt$treat_alt, "}; only ", length(cells),
           " present. (In this project only H3K27me3 has a UNC arm.)")
    if (any(cells < 2))
      stop("A cell has <2 replicates; the interaction is not reliably estimable.")
  }
  ss
}

# -----------------------------------------------------------------------------
# Count reads into windows, then filter to enriched windows. Filter is over IgG
# (filterWindowsControl) if --igg given, else over genomic background
# (filterWindowsGlobal, original behavior). Returns list(win, bins).
# -----------------------------------------------------------------------------
count_and_filter_windows <- function(bam_files, sample_ids, opt) {
  library(csaw)
  param <- readParam(pe = "both", max.frag = opt$frag_max)
  if (!is.null(opt$blacklist)) {
    library(rtracklayer)
    param <- reform(param, discard = import(opt$blacklist))
    message("[INFO] Blacklist applied: ", opt$blacklist)
  }

  message("[INFO] Counting windows (width=", opt$window, ", spacing=", opt$spacing, ")...")
  win <- windowCounts(bam_files, width = opt$window, spacing = opt$spacing, param = param)
  colnames(win) <- sample_ids

  bins <- windowCounts(bam_files, bin = TRUE, width = opt$bg_bin, param = param)

  if (!is.null(opt$igg)) {
    igg_bams <- trimws(strsplit(opt$igg, ",")[[1]])
    message("[INFO] IgG-control filtering with ", length(igg_bams), " IgG BAM(s).")
    win_ctrl <- windowCounts(igg_bams, width = opt$window, spacing = opt$spacing, param = param)
    bin_ctrl <- windowCounts(igg_bams, bin = TRUE, width = opt$bg_bin, param = param)
    scale.info <- scaleControlFilter(bins, bin_ctrl)          # corrects ChIP-vs-IgG depth
    keep <- filterWindowsControl(win, win_ctrl, prior.count = 2,
                                 scale.info = scale.info)$filter > opt$min_lfc_bg
  } else {
    message("[INFO] Genomic-background filtering (", opt$bg_bin, " bp bins).")
    keep <- filterWindowsGlobal(win, bins)$filter > opt$min_lfc_bg
  }
  win <- win[keep, ]
  message("[INFO] Windows retained after filter: ", sum(keep), " of ", length(keep))
  list(win = win, bins = bins)
}

# -----------------------------------------------------------------------------
# Build the model design and identify the coefficient to test.
#   contrast    : ~group (or ~Batch + group); test last coef = g1 vs g2.
#   interaction : ~Condition*Treatment (or ~Batch + ...); test the ":" coef.
# Returns design, coef (index), group factor (for QC), ref_cols (baseline cols).
# -----------------------------------------------------------------------------
build_design <- function(ss, opt) {
  batch <- factor(ss$Batch)
  has_batch <- nlevels(batch) > 1
  ref_cols <- which(ss$ref_cell)

  if (opt$test == "contrast") {
    group <- factor(ss$grp, levels = c("g2", "g1"))         # g2 = reference
    design <- if (has_batch) model.matrix(~ batch + group) else model.matrix(~ group)
    if (has_batch) message("[INFO] >1 batch -> ~Batch + group.")
    coef <- ncol(design)                                    # last coef = g1 vs g2
    qc_group <- factor(ss$Condition)
  } else {
    # Match --cond_ref/--treat_ref to the sheet's actual casing (relevel needs an exact level)
    cref <- unique(ss$Condition)[tolower(unique(ss$Condition)) == tolower(opt$cond_ref)]
    tref <- unique(ss$Treatment)[tolower(unique(ss$Treatment)) == tolower(opt$treat_ref)]
    if (length(cref) != 1)
      stop("Need exactly one Condition level matching '", opt$cond_ref, "'. Found: ",
          paste(unique(ss$Condition), collapse = ", "))
    if (length(tref) != 1)
      stop("Need exactly one Treatment level matching '", opt$treat_ref, "'. Found: ",
          paste(unique(ss$Treatment), collapse = ", "))
    Condition <- relevel(factor(ss$Condition), ref = cref)
    Treatment <- relevel(factor(ss$Treatment), ref = tref)
    design <- if (has_batch) model.matrix(~ batch + Condition * Treatment) else
                             model.matrix(~ Condition * Treatment)
    if (has_batch) message("[INFO] >1 batch -> ~Batch + Condition*Treatment.")
    inter <- grep(":", colnames(design))                    # the interaction column(s)
    if (length(inter) != 1)
      stop("Expected exactly one interaction coefficient, found: ",
           paste(colnames(design)[inter], collapse = ", "))
    coef <- inter
    message("[INFO] Design columns: ", paste(colnames(design), collapse = " | "))
    message("[INFO] Testing interaction coef: ", colnames(design)[coef])
    qc_group <- factor(paste(ss$Condition, ss$Treatment))
  }
  list(design = design, coef = coef, group = qc_group, ref_cols = ref_cols)
}

# -----------------------------------------------------------------------------
# Per-region baseline abundance (mean log-CPM over the reference samples),
# aggregated from windows to merged regions, used to stratify the flattening.
# -----------------------------------------------------------------------------
baseline_by_region <- function(y, ref_cols, ids, n_regions) {
  library(edgeR)
  lc <- cpm(y, log = TRUE)                       # windows x samples
  base_win <- if (length(ref_cols) > 1) rowMeans(lc[, ref_cols, drop = FALSE]) else lc[, ref_cols]
  # ids gives region membership per window; mean per region, aligned to 1..n_regions
  agg <- tapply(base_win, ids, mean)
  out <- rep(NA_real_, n_regions)
  out[as.integer(names(agg))] <- as.numeric(agg)
  out
}

# -----------------------------------------------------------------------------
# Run one normalization scheme end-to-end.
# -----------------------------------------------------------------------------
run_one_normalization <- function(win, bins, ss, design_info, norm_tag,
                                   opt, tag, spikein_df = NULL) {
  library(csaw); library(edgeR)
  message("[INFO] === Normalization: ", norm_tag, " ===")

  y      <- asDGEList(win)
  design <- design_info$design
  coef   <- design_info$coef

  if (!is.null(spikein_df)) {
    facs <- spikein_df$factor[match(ss$SampleID, spikein_df$SampleID)]
    if (any(is.na(facs)))
      stop("Spike-in file missing factors for: ",
           paste(ss$SampleID[is.na(facs)], collapse = ", "))
    y$samples$norm.factors <- facs / exp(mean(log(facs)))
  } else {
    y$samples$norm.factors <- normFactors(bins, se.out = TRUE)$norm.factors
  }

  y      <- estimateDisp(y, design)
  fit    <- glmQLFit(y, design, robust = TRUE)
  result <- glmQLFTest(fit, coef = coef)

  merged <- mergeWindows(rowRanges(win), tol = opt$spacing, max.width = 50000)
  tab    <- combineTests(merged$ids, result$table)
  out    <- cbind(as.data.frame(granges(merged$regions)), as.data.frame(tab))

  # ---- baseline stratification (flattening signature) ----------------------
  out$baseline_logCPM <- baseline_by_region(y, design_info$ref_cols,
                                             merged$ids, nrow(out))
  brks <- quantile(out$baseline_logCPM, probs = seq(0, 1, length.out = opt$n_strata + 1),
                   na.rm = TRUE)
  brks[1] <- -Inf; brks[length(brks)] <- Inf
  lab <- if (opt$n_strata == 3) c("low", "mid", "high") else paste0("q", seq_len(opt$n_strata))
  out$baseline_stratum <- as.character(cut(out$baseline_logCPM, breaks = brks,
                                           labels = lab, include.lowest = TRUE))

  base <- file.path(opt$outdir, paste0(tag, "_", norm_tag))
  write.csv(out, paste0(base, "_all_regions.csv"), row.names = FALSE)

  lfc_col <- intersect(c("rep.logFC", "best.logFC"), colnames(out))
  lfc_col <- if (length(lfc_col)) lfc_col[1] else NA
  sig <- out[!is.na(out$FDR) & out$FDR < opt$fdr, ]
  if (!is.na(lfc_col)) sig <- sig[abs(sig[[lfc_col]]) >= opt$lfc, ]
  write.csv(sig, paste0(base, "_significant.csv"), row.names = FALSE)
  message("[INFO] ", norm_tag, ": ", nrow(sig), " significant regions of ", nrow(out),
          " (FDR<", opt$fdr,
          if (!is.na(lfc_col)) paste0(" & |log2FC|>=", opt$lfc) else "", ")")

  # overall and per-stratum direction balance (the flattening check)
  if (!is.na(lfc_col) && nrow(sig) > 0) {
    up <- sum(sig[[lfc_col]] > 0); dn <- sum(sig[[lfc_col]] < 0)
    message(sprintf("[INFO] %s: direction balance overall — up=%d (%.0f%%), down=%d (%.0f%%)",
                    norm_tag, up, 100*up/(up+dn), dn, 100*dn/(up+dn)))
    for (s in lab) {
      idx <- which(sig$baseline_stratum == s)
      if (length(idx)) {
        u <- sum(sig[[lfc_col]][idx] > 0); d <- sum(sig[[lfc_col]][idx] < 0)
        verdict <- if (d > u) "net LOSS" else if (u > d) "net GAIN" else "balanced"
        message(sprintf("[STRAT] %-4s baseline: sig up=%d, down=%d -> %s",
                        s, u, d, verdict))
      }
    }
    message("[STRAT] (UNC contrasts should show net LOSS at 'high' and net GAIN at 'low' = flattening)")
  }

  save_volcano_plot <- function(out, opt, base, tag, norm_tag) {
    lfc_col <- intersect(c("rep.logFC", "best.logFC"), colnames(out))
    if (!length(lfc_col)) { message("  [QC] Volcano skipped: no logFC column."); return(invisible(NULL)) }
    lfc_col <- lfc_col[1]
    x <- out[[lfc_col]]; fdr <- out$FDR
    ok <- is.finite(x) & is.finite(fdr) & fdr > 0
    x <- x[ok]; fdr <- fdr[ok]; y <- -log10(fdr)
    sig <- fdr < opt$fdr & abs(x) >= opt$lfc
    col <- ifelse(sig & x > 0, "firebrick", ifelse(sig & x < 0, "royalblue", "grey70"))
    render_plot(function() {
      plot(x, y, pch = 19, cex = 0.5, col = col,
           xlab = "log2 fold-change", ylab = "-log10(FDR)",
           main = paste(tag, norm_tag, "volcano"))
      abline(h = -log10(opt$fdr), col = "grey40", lty = 2)
      abline(v = c(-opt$lfc, opt$lfc), col = "grey40", lty = 2)
      legend("topleft", bty = "n", pch = 19, col = c("firebrick", "royalblue", "grey70"),
             legend = c(sprintf("up (n=%d)", sum(sig & x > 0)),
                        sprintf("down (n=%d)", sum(sig & x < 0)), "n.s."))
    }, paste0(base, "_volcano"))
  }

  save_qc_plots(y, result, design_info$group, ss, opt, base, tag, norm_tag)
  save_volcano_plot(out, opt, base, tag, norm_tag)
  invisible(out)
}

# -----------------------------------------------------------------------------
# render_plot(): same plot to BOTH a PDF (paper) and PNG (preview).
# -----------------------------------------------------------------------------
render_plot <- function(draw_fn, base_path, w = 7, h = 6) {
  tryCatch({ pdf(paste0(base_path, ".pdf"), width = w, height = h); draw_fn() },
           error = function(e) message("  [QC] PDF skipped: ", conditionMessage(e)),
           finally = if (!is.null(dev.list())) dev.off())
  tryCatch({ png(paste0(base_path, ".png"), width = w*150, height = h*150, res = 150); draw_fn() },
           error = function(e) message("  [QC] PNG skipped: ", conditionMessage(e)),
           finally = if (!is.null(dev.list())) dev.off())
}

# -----------------------------------------------------------------------------
# QC plots (each as .pdf and .png):
#   *_MDS    : sample clustering, coloured by Condition, rep labels offset.
#   *_MA     : smear/MA with +/- lfc guides.
#   *_repcor : NEW — Spearman correlation heatmap of per-window log-CPM across
#              samples, so replicate concordance is an explicit figure (the
#              reproducibility argument a Q1 reviewer expects).
# -----------------------------------------------------------------------------
save_qc_plots <- function(y, result, group, ss, opt, base, tag, norm_tag) {
  library(edgeR)
  cond     <- factor(ss$Condition)
  cond_col <- as.integer(cond)
  full_lab <- paste(ss$Condition, paste0("R", ss$Replicate))

  render_plot(function() {
    op <- par(mar = c(5, 5, 3, 9), xpd = TRUE); on.exit(par(op))
    mds <- plotMDS(cpm(y, log = TRUE), plot = FALSE)
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

  # replicate-correlation heatmap (Spearman on per-window log-CPM)
  render_plot(function() {
    M <- cor(cpm(y, log = TRUE), method = "spearman")
    op <- par(mar = c(9, 9, 3, 2)); on.exit(par(op))
    n <- ncol(M)
    image(seq_len(n), seq_len(n), t(M[n:1, , drop = FALSE]),
          axes = FALSE, xlab = "", ylab = "", zlim = c(min(M), 1),
          main = paste(tag, norm_tag, "replicate Spearman"))
    axis(1, at = seq_len(n), labels = colnames(M), las = 2, cex.axis = 0.7)
    axis(2, at = seq_len(n), labels = rev(colnames(M)), las = 2, cex.axis = 0.7)
    for (i in seq_len(n)) for (j in seq_len(n))
      text(j, n - i + 1, sprintf("%.2f", M[i, j]), cex = 0.6)
  }, paste0(base, "_repcor"))
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
main <- function() {
  opt <- parse_cli()
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

  ss <- load_samples(opt)

  # human-readable tag
  if (opt$test == "contrast") {
    tag <- paste0(opt$mark, "_",
                  gsub("[|]", "-", opt$group1), "_vs_", gsub("[|]", "-", opt$group2))
  } else {
    tag <- paste0(opt$mark, "_interaction_",
                  opt$cond_alt, "x", opt$treat_alt, "_vs_",
                  opt$cond_ref, "x", opt$treat_ref)
  }

  counts      <- count_and_filter_windows(ss$bamReads, ss$SampleID, opt)
  design_info <- build_design(ss, opt)

  run_one_normalization(counts$win, counts$bins, ss, design_info,
                        "libsize", opt, tag, spikein_df = NULL)

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