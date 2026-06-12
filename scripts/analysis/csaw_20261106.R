# =============================================================================
# csaw.R — Window-based differential binding for BROAD histone marks
#          (v3: + RTM-robust erosion stratification folded in)
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
#   contrast    (default) Two-group comparison:
#                 --group1 NUMERATOR --group2 REFERENCE  (Condition|Treatment)
#                 Migration:     --group1 "B10|none"     --group2 "Parental|none"
#                 UNC in B10:    --group1 "B10|UNC"      --group2 "B10|none"
#                 UNC in parent: --group1 "Parental|UNC" --group2 "Parental|none"
#   interaction  Factorial Condition*Treatment; tests the INTERACTION coefficient
#                 (does UNC's effect differ between B10 and Parental?). Requires
#                 all four cells (only H3K27me3 has them here).
#   Matching is CASE-INSENSITIVE throughout ("Parental"=="parental").
#
# BASELINE STRATIFICATION (the PRC2 "flattening" signature) -------------------
#   PRC2 inhibition erodes high-enrichment peaks while signal spreads into
#   low/background regions (Romero et al. 2024). A single genome-wide test lets
#   those opposite moves cancel, so we tag every region with its baseline
#   abundance, split into low/mid/high strata, and report the up/down balance per
#   stratum. Net LOSS at 'high' + net GAIN at 'low' = flattening.
#
#   *** RTM CAVEAT + ROBUST STRATIFIERS (v3) ***
#   The default 'baseline_stratum' ranks by the SAME reference group whose noise
#   also enters the tested logFC (logFC = treated - reference). Baseline and
#   logFC therefore share that noise and are mechanically anti-correlated: a
#   noisily-high reference estimate lands a region in 'high' AND gives it a
#   negative logFC. That manufactures the high-loses/low-gains shape via
#   regression-to-the-mean (RTM), independent of biology. Two RTM-resistant
#   stratifiers are now reported alongside it:
#     --peaks         (PRIMARY, RTM-IMMUNE)  Peaks are called vs IgG/input, a
#                     different model than the per-window edgeR test, so peak
#                     membership does not inherit the contrast's logFC noise. It
#                     also captures where the mark IS in the treated cells, so it
#                     does NOT lose migration-installed signal. Use a BROAD call
#                     for H3K27me3 (macs2 --broad, epic2, or SEACR).
#     --baseline_group (OPTIONAL, threshold-free) An independent biological group
#                     counted in the SAME windows and used ONLY to rank baseline
#                     (not added to the design). Valid when independent of the
#                     tested coefficient, e.g. 'Parental|none' for a B10|UNC-vs-
#                     B10|none contrast. NOTE: parental abundance under-weights
#                     migration-installed loci, so treat it as a sensitivity
#                     check secondary to --peaks, not the primary.
#   The run prints a [STRAT:samedata] / [STRAT:peak] / [STRAT:indep] block each,
#   plus an [STRAT] RTM verdict: does the high-vs-low loss gradient SURVIVE the
#   RTM-resistant stratifier? If yes, the flattening is real; if it collapses,
#   it was largely RTM. in_peak / baseline_indep_stratum are written to
#   all_regions.csv so Track 1 can use the robust flag directly (no second run).
#   (Promoter-vs-intergenic, the caller-free partition, lives in the optional
#   Python post-processor robust_erosion_stratification.py for fast re-slicing
#   without re-counting; this script covers the two stratifiers that need the
#   csaw regions themselves.)
#
# IgG-CONTROL FILTERING (optional, CUT&RUN-appropriate) -----------------------
#   If you have IgG, pass --igg bam1,bam2,... and windows are kept by enrichment
#   over IgG (filterWindowsControl) instead of genomic background. Sanity-check
#   the retained-window count the first time. (Not unit-tested locally.)
#
# SPIKE-IN (E. coli carryover; RRPM/Orlando) ---------------------------------
#   --spikein factors.csv (columns: SampleID,spikein_reads) where spikein_reads
#   is the ABSOLUTE E. coli read count per sample (= total reads x %E.coli; build
#   with build_spikein_csv.py). The effective library size is set proportional to
#   spike-in reads (Orlando RRPM: endogenous / exogenous), which REPLACES TMM and
#   writes a second 'spikein' result set. Use it for the H3K27me3 UNC contrasts,
#   where global loss is the claim and library-size normalization masks it.
#   VALID ONLY across samples with comparable cell number + ConA beads (carryover
#   scales with both); prefer WITHIN-BATCH contrasts. Spike-in fixes the GLOBAL
#   scale (total loss, up/down balance) but not which loci change relative to
#   each other, so the high-baseline erosion and the interaction already hold
#   under library-size; treat agreement between 'libsize' and 'spikein' as a
#   robustness check, and anchor the global claim with an orthogonal assay
#   (e.g. western blot for total H3K27me3).
#
# REQUIRED SAMPLE-SHEET COLUMNS ----------------------------------------------
#   SampleID, Factor, Condition, Treatment, Replicate, Batch, bamReads
#
# EXAMPLE USAGE ---------------------------------------------------------------
#   # UNC reversal in B10, with RTM-robust stratification in ONE run:
#   Rscript csaw.R --sheet ss.csv --mark H3K27me3 \
#       --group1 "B10|UNC" --group2 "B10|none" \
#       --peaks B10none_H3K27me3_R1.broadPeak,B10none_H3K27me3_R2.broadPeak \
#       --baseline_group "Parental|none" \
#       --outdir results/diff/H3K27me3_B10unc
#
#   # Interaction (peaks are the right RTM-resistant stratifier here):
#   Rscript csaw.R --sheet ss.csv --mark H3K27me3 --test interaction \
#       --peaks B10none_H3K27me3.broadPeak \
#       --outdir results/diff/H3K27me3_interaction
#
# OUTPUTS (into --outdir) -----------------------------------------------------
#   <tag>_<norm>_all_regions.csv   coords, logFC, FDR, direction,
#                                  baseline_logCPM, baseline_stratum,
#                                  baseline_indep_logCPM, baseline_indep_stratum,
#                                  in_peak
#   <tag>_<norm>_significant.csv   subset passing FDR & |log2FC|
#   <tag>_<norm>_{MDS,MA,volcano,repcor}.{pdf,png}
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
                help = "REQUIRED. Sample sheet CSV. Columns: SampleID, Factor, Condition, Treatment, Replicate, Batch, bamReads."),
    make_option("--mark", type = "character",
                help = "REQUIRED. Factor/antibody, matched case-insensitively to Factor. e.g. H3K27me3 | H3K9me3 | H3K4me1."),

    make_option("--test", type = "character", default = "contrast",
                help = "Test mode. ALLOWED: 'contrast' | 'interaction'. [default %default]"),

    make_option("--group1", type = "character", default = "B10|none",
                help = "[contrast only] Numerator 'Condition|Treatment' (case-insensitive). [default %default]"),
    make_option("--group2", type = "character", default = "Parental|none",
                help = "[contrast only] Reference 'Condition|Treatment'. Positive logFC = higher in group1. [default %default]"),

    make_option("--cond_ref", type = "character", default = "Parental",
                help = "[interaction only] Reference Condition. [default %default]"),
    make_option("--treat_ref", type = "character", default = "none",
                help = "[interaction only] Reference Treatment. [default %default]"),
    make_option("--cond_alt", type = "character", default = "B10",
                help = "[interaction only] Non-reference Condition. [default %default]"),
    make_option("--treat_alt", type = "character", default = "UNC",
                help = "[interaction only] Non-reference Treatment. [default %default]"),

    make_option("--outdir", type = "character", default = "csaw_out",
                help = "Output directory. [default %default]"),

    make_option("--window", type = "integer", default = 2000L,
                help = "Window width (bp). [default %default]"),
    make_option("--spacing", type = "integer", default = 500L,
                help = "Spacing between window starts (bp); also merge tolerance. [default %default]"),
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
                help = "Optional ENCODE blacklist BED. [default: none]"),
    make_option("--igg", type = "character", default = NULL,
                help = "Optional comma-separated IgG BAM path(s) for filterWindowsControl. [default: none]"),
    make_option("--spikein", type = "character", default = NULL,
                help = "Spike-in (Orlando/RRPM) normalization. Either the literal 'sheet' to read the sample sheet's spikein_reads column (one-sheet workflow), or a CSV path with columns SampleID,spikein_reads (absolute E. coli reads). Sets effective library size proportional to spike-in reads, replacing TMM, and writes a second 'spikein' result set. [default: none]"),

    # --- RTM-robust stratifiers (v3) -----------------------------------------
    make_option("--peaks", type = "character", default = NULL,
                help = "Optional peak file(s): BED/narrowPeak/broadPeak/SEACR-bed, comma-separated (UNIONed). Regions overlapping a peak are flagged 'in_peak' and the flattening is reported peak-vs-nonpeak. RTM-IMMUNE primary robust stratifier; for H3K27me3 use a BROAD call (called on the reference cell, e.g. B10|none). [default: none]"),
    make_option("--baseline_group", type = "character", default = NULL,
                help = "Optional 'Condition|Treatment' group counted in the SAME windows and used ONLY to rank baseline abundance for an RTM-proof tertile (NOT added to the design). Valid when independent of the tested coefficient, e.g. 'Parental|none' for a B10|UNC-vs-B10|none contrast. Secondary to --peaks (under-weights migration-installed loci). [default: none]")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$sheet)) stop("--sheet is required.")
  if (is.null(opt$mark))  stop("--mark is required (e.g. H3K27me3).")
  if (!opt$test %in% c("contrast", "interaction"))
    stop("--test must be 'contrast' or 'interaction'.")
  opt
}

# -----------------------------------------------------------------------------
parse_group <- function(spec) {
  parts <- strsplit(spec, "\\|")[[1]]
  if (length(parts) != 2)
    stop("Group spec must be 'Condition|Treatment', got: ", spec)
  list(cond = tolower(trimws(parts[1])), treat = tolower(trimws(parts[2])))
}

# -----------------------------------------------------------------------------
# Sample selection (analysis samples).
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
    ss$ref_cell <- ss$grp == "g2"
    ss <- ss[order(ss$grp, ss$Replicate), ]
    tab <- table(ss$grp)
    message("[INFO] ", opt$mark, " [contrast]: ", nrow(ss),
            " samples (g1=", sum(ss$grp == "g1"), ", g2=", sum(ss$grp == "g2"), ")")
    if (any(tab < 2))
      stop("A group has <2 replicates; NB dispersion is not estimable.")
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
      stop("Interaction needs all four cells; only ", length(cells), " present.")
    if (any(cells < 2))
      stop("A cell has <2 replicates; the interaction is not reliably estimable.")
  }
  ss
}

# -----------------------------------------------------------------------------
# Independent baseline group (for --baseline_group): the BAMs of an independent
# 'Condition|Treatment' cell of the SAME mark, used only to rank abundance.
# Returns list(bams, ids) or NULL.
# -----------------------------------------------------------------------------
load_baseline_group <- function(opt) {
  if (is.null(opt$baseline_group)) return(NULL)
  g <- parse_group(opt$baseline_group)
  full <- read.csv(opt$sheet, stringsAsFactors = FALSE)
  bg <- full[tolower(full$Factor) == tolower(opt$mark) &
             tolower(trimws(full$Condition)) == g$cond &
             tolower(trimws(full$Treatment)) == g$treat, ]
  if (nrow(bg) == 0)
    stop("--baseline_group '", opt$baseline_group, "' matched no ", opt$mark, " samples.")
  message("[INFO] Independent baseline group '", opt$baseline_group, "': ",
          nrow(bg), " sample(s) (ranking only, excluded from design).")
  list(bams = bg$bamReads, ids = bg$SampleID)
}

# -----------------------------------------------------------------------------
# Count reads into windows, filter to enriched windows. Optionally also count an
# independent baseline group into the SAME windows (subset by the same filter),
# returned as win_base for RTM-proof ranking. Returns list(win, bins, win_base).
# -----------------------------------------------------------------------------
count_and_filter_windows <- function(bam_files, sample_ids, opt, baseline = NULL) {
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

  # independent baseline group counted into the SAME tiling (BEFORE filtering, so
  # rows align 1:1 with 'win'); subset by the same keep vector below.
  win_base <- NULL
  if (!is.null(baseline)) {
    message("[INFO] Counting independent baseline group into the same windows...")
    win_base <- windowCounts(baseline$bams, width = opt$window, spacing = opt$spacing,
                             param = param)
    colnames(win_base) <- baseline$ids
    if (nrow(win_base) != nrow(win))
      stop("Baseline-group windows (", nrow(win_base), ") do not align with main windows (",
           nrow(win), "). Check identical genome/param.")
  }

  if (!is.null(opt$igg)) {
    igg_bams <- trimws(strsplit(opt$igg, ",")[[1]])
    message("[INFO] IgG-control filtering with ", length(igg_bams), " IgG BAM(s).")
    win_ctrl <- windowCounts(igg_bams, width = opt$window, spacing = opt$spacing, param = param)
    bin_ctrl <- windowCounts(igg_bams, bin = TRUE, width = opt$bg_bin, param = param)
    scale.info <- scaleControlFilter(bins, bin_ctrl)
    keep <- filterWindowsControl(win, win_ctrl, prior.count = 2,
                                 scale.info = scale.info)$filter > opt$min_lfc_bg
  } else {
    message("[INFO] Genomic-background filtering (", opt$bg_bin, " bp bins).")
    keep <- filterWindowsGlobal(win, bins)$filter > opt$min_lfc_bg
  }
  win <- win[keep, ]
  if (!is.null(win_base)) win_base <- win_base[keep, ]   # same filter -> still aligned
  message("[INFO] Windows retained after filter: ", sum(keep), " of ", length(keep))
  list(win = win, bins = bins, win_base = win_base)
}

# -----------------------------------------------------------------------------
build_design <- function(ss, opt) {
  batch <- factor(ss$Batch)
  has_batch <- nlevels(batch) > 1
  ref_cols <- which(ss$ref_cell)

  if (opt$test == "contrast") {
    group <- factor(ss$grp, levels = c("g2", "g1"))
    design <- if (has_batch) model.matrix(~ batch + group) else model.matrix(~ group)
    if (has_batch) message("[INFO] >1 batch -> ~Batch + group.")
    coef <- ncol(design)
    qc_group <- factor(ss$Condition)
  } else {
    cond_f <- factor(ss$Condition); treat_f <- factor(ss$Treatment)
    cref <- levels(cond_f)[tolower(levels(cond_f)) == tolower(opt$cond_ref)]
    tref <- levels(treat_f)[tolower(levels(treat_f)) == tolower(opt$treat_ref)]
    if (length(cref) != 1)
      stop("--cond_ref '", opt$cond_ref, "' not found among: ",
           paste(levels(cond_f), collapse = ", "))
    if (length(tref) != 1)
      stop("--treat_ref '", opt$treat_ref, "' not found among: ",
           paste(levels(treat_f), collapse = ", "))
    Condition <- relevel(cond_f, ref = cref)
    Treatment <- relevel(treat_f, ref = tref)
    design <- if (has_batch) model.matrix(~ batch + Condition * Treatment) else
                             model.matrix(~ Condition * Treatment)
    if (has_batch) message("[INFO] >1 batch -> ~Batch + Condition*Treatment.")
    inter <- grep(":", colnames(design))
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
# Mean log-CPM over a set of sample columns, aggregated windows -> regions.
# Used for both the same-data baseline (ref_cols of the contrast) and the
# independent baseline (all columns of win_base).
# -----------------------------------------------------------------------------
region_mean_logcpm <- function(lc, cols, ids, n_regions) {
  base_win <- if (length(cols) > 1) rowMeans(lc[, cols, drop = FALSE]) else lc[, cols]
  agg <- tapply(base_win, ids, mean)
  out <- rep(NA_real_, n_regions)
  out[as.integer(names(agg))] <- as.numeric(agg)
  out
}

# tertile (or n_strata-tile) a numeric vector into low/mid/high labels
strata_from <- function(values, n_strata, lab) {
  brks <- quantile(values, probs = seq(0, 1, length.out = n_strata + 1), na.rm = TRUE)
  brks[1] <- -Inf; brks[length(brks)] <- Inf
  as.character(cut(values, breaks = brks, labels = lab, include.lowest = TRUE))
}

# per-stratum up/down balance; prints and returns named pct_down per level
strat_report <- function(sig, lfc_col, strat, levels, tagname) {
  res <- setNames(rep(NA_real_, length(levels)), levels)
  for (lv in levels) {
    idx <- which(strat == lv)
    u <- sum(sig[[lfc_col]][idx] > 0); d <- sum(sig[[lfc_col]][idx] < 0)
    tot <- u + d
    res[lv] <- if (tot) d / tot else NA_real_
    verdict <- if (d > u) "net LOSS" else if (u > d) "net GAIN" else "balanced"
    message(sprintf("[STRAT:%s] %-9s: sig up=%d, down=%d -> %s",
                    tagname, lv, u, d, verdict))
  }
  res
}

# -----------------------------------------------------------------------------
run_one_normalization <- function(win, bins, ss, design_info, norm_tag,
                                   opt, tag, spikein_df = NULL,
                                   win_base = NULL, peaks_gr = NULL) {
  library(csaw); library(edgeR)
  message("[INFO] === Normalization: ", norm_tag, " ===")

  y      <- asDGEList(win)
  design <- design_info$design
  coef   <- design_info$coef

  if (!is.null(spikein_df)) {
    # Spike-in (Orlando/RRPM) normalization. RRPM scales endogenous signal by
    # 1/spike-in reads, i.e. the EFFECTIVE library size must be proportional to
    # the spike-in read count. edgeR effective = lib.size * norm.factors, and
    # csaw's lib.size = win$totals (genome-wide TARGET fragments), so the correct
    # norm.factor is (spikein_reads / lib.size), centred to geomean 1. (Passing
    # spikein_reads straight in as norm.factors would double-count depth.)
    spike <- spikein_df$spikein_reads[match(ss$SampleID, spikein_df$SampleID)]
    if (any(is.na(spike)))
      stop("Spike-in file missing spikein_reads for: ",
           paste(ss$SampleID[is.na(spike)], collapse = ", "),
           "\n  (SampleID must match the sheet exactly; reconcile names first.)")
    if (any(spike <= 0))
      stop("Non-positive spikein_reads for: ",
           paste(ss$SampleID[spike <= 0], collapse = ", "))
    nf <- spike / y$samples$lib.size
    y$samples$norm.factors <- nf / exp(mean(log(nf)))
    message("[INFO] Spike-in (RRPM) norm: effective library size set proportional ",
            "to E. coli reads (range ", round(min(spike)), "-", round(max(spike)), ").")
  } else {
    y$samples$norm.factors <- normFactors(bins, se.out = TRUE)$norm.factors
  }

  y      <- estimateDisp(y, design)
  fit    <- glmQLFit(y, design, robust = TRUE)
  result <- glmQLFTest(fit, coef = coef)

  merged <- mergeWindows(rowRanges(win), tol = opt$spacing, max.width = 50000)
  tab    <- combineTests(merged$ids, result$table)
  out    <- cbind(as.data.frame(granges(merged$regions)), as.data.frame(tab))
  lab    <- if (opt$n_strata == 3) c("low", "mid", "high") else
              paste0("q", seq_len(opt$n_strata))

  # ---- (a) same-data baseline (RTM-prone) ----------------------------------
  lc <- cpm(y, log = TRUE)
  out$baseline_logCPM  <- region_mean_logcpm(lc, design_info$ref_cols, merged$ids, nrow(out))
  out$baseline_stratum <- strata_from(out$baseline_logCPM, opt$n_strata, lab)

  # ---- (b) independent baseline (RTM-proof tertile), if --baseline_group ---
  has_indep <- !is.null(win_base)
  if (has_indep) {
    lcb <- cpm(asDGEList(win_base), log = TRUE)
    out$baseline_indep_logCPM  <- region_mean_logcpm(lcb, seq_len(ncol(lcb)),
                                                     merged$ids, nrow(out))
    out$baseline_indep_stratum <- strata_from(out$baseline_indep_logCPM, opt$n_strata, lab)
  }

  # ---- (c) peak overlap (RTM-immune), if --peaks ---------------------------
  has_peaks <- !is.null(peaks_gr)
  if (has_peaks) out$in_peak <- overlapsAny(merged$regions, peaks_gr)

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

  # ---- direction balance: overall + each stratifier + RTM verdict ----------
  if (!is.na(lfc_col) && nrow(sig) > 0) {
    up <- sum(sig[[lfc_col]] > 0); dn <- sum(sig[[lfc_col]] < 0)
    message(sprintf("[INFO] %s: direction balance overall — up=%d (%.0f%%), down=%d (%.0f%%)",
                    norm_tag, up, 100*up/(up+dn), dn, 100*dn/(up+dn)))

    res_same  <- strat_report(sig, lfc_col, sig$baseline_stratum, lab, "samedata")
    res_indep <- if (has_indep)
      strat_report(sig, lfc_col, sig$baseline_indep_stratum, lab, "indep") else NULL
    res_peak  <- if (has_peaks)
      strat_report(sig, lfc_col, ifelse(sig$in_peak, "in_peak", "out_peak"),
                   c("in_peak", "out_peak"), "peak") else NULL

    grad <- function(res, hi, lo) if (is.na(res[hi]) || is.na(res[lo])) NA_real_ else
                                   as.numeric(res[hi] - res[lo])
    g_same <- grad(res_same, "high", "low")
    message(sprintf("[STRAT] same-data high-vs-low %%down gradient = %+.0f%%", 100*g_same))
    if (has_indep) {
      g_ind <- grad(res_indep, "high", "low")
      retained <- if (is.finite(g_same) && g_same > 1e-6) g_ind / g_same else NA_real_
      verdict <- if (!is.na(g_ind) && g_ind > 0.15 && (is.na(retained) || retained > 0.5))
                   "REAL (survives independent baseline)" else
                   "LARGELY RTM (collapses under independent baseline)"
      message(sprintf("[STRAT] independent-baseline gradient = %+.0f%%  -> %s",
                      100*g_ind, verdict))
    }
    if (has_peaks) {
      g_pk <- grad(res_peak, "in_peak", "out_peak")
      message(sprintf("[STRAT] peak-vs-nonpeak %%down gradient = %+.0f%%  (in-peak loss = the literal peak-erosion claim)",
                      100*g_pk))
    }
    message("[STRAT] (flattening is REAL if the high-baseline LOSS survives the RTM-resistant stratifier; ",
            "low-baseline 'gain' + the same-data gradient alone are RTM-prone.)")
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
           main = paste(tag, "\n", norm_tag, "volcano"))
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
render_plot <- function(draw_fn, base_path, w = 7, h = 6) {
  tryCatch({ pdf(paste0(base_path, ".pdf"), width = w, height = h); draw_fn() },
           error = function(e) message("  [QC] PDF skipped: ", conditionMessage(e)),
           finally = if (!is.null(dev.list())) dev.off())
  tryCatch({ png(paste0(base_path, ".png"), width = w*150, height = h*150, res = 150); draw_fn() },
           error = function(e) message("  [QC] PNG skipped: ", conditionMessage(e)),
           finally = if (!is.null(dev.list())) dev.off())
}

# -----------------------------------------------------------------------------
save_qc_plots <- function(y, result, group, ss, opt, base, tag, norm_tag) {
  library(edgeR)
  cond_f   <- factor(ss$Condition)
  treat_f  <- factor(ss$Treatment)
  pch_pool <- c(16, 17, 15, 18, 8, 7)
  col_pool <- c("grey30", "firebrick", "royalblue", "darkgreen", "purple", "orange")
  pt_pch   <- pch_pool[as.integer(cond_f)]
  pt_col   <- col_pool[as.integer(treat_f)]
  rep_lab  <- paste0("R", ss$Replicate)
  full_lab <- paste0(ss$Condition, " ", ss$Treatment, " R", ss$Replicate)

  render_plot(function() {
    op <- par(mar = c(5, 5, 3, 11), xpd = TRUE); on.exit(par(op))
    mds <- plotMDS(cpm(y, log = TRUE), plot = FALSE)
    jx <- jitter(mds$x, amount = diff(range(mds$x)) * 0.01)
    plot(jx, mds$y, pch = pt_pch, col = pt_col, cex = 1.8, lwd = 2,
         xlab = "Leading logFC dim 1", ylab = "Leading logFC dim 2",
         main = paste(tag, "\n", norm_tag, "MDS"))
    text(jx, mds$y, labels = rep_lab, pos = 1, cex = 0.7)
    legend("topright", inset = c(-0.34, 0), legend = levels(cond_f),
           pch = pch_pool[seq_along(levels(cond_f))], col = "black",
           bty = "n", title = "Condition")
    legend("bottomright", inset = c(-0.34, 0), legend = levels(treat_f),
           pch = 16, col = col_pool[seq_along(levels(treat_f))],
           bty = "n", title = "Treatment")
  }, paste0(base, "_MDS"))

  render_plot(function() {
    plotSmear(result, main = paste(tag, norm_tag, "\n", "smear (MA)"))
    abline(h = c(-opt$lfc, opt$lfc), col = "blue", lty = 2)
  }, paste0(base, "_MA"))

  render_plot(function() {
    ord  <- order(ss$Condition, ss$Treatment, ss$Replicate)
    M    <- cor(cpm(y, log = TRUE)[, ord], method = "spearman")
    labs <- full_lab[ord]
    op <- par(mar = c(11, 11, 3, 2)); on.exit(par(op))
    n <- ncol(M)
    image(seq_len(n), seq_len(n), t(M[n:1, , drop = FALSE]),
          axes = FALSE, xlab = "", ylab = "", zlim = c(min(M), 1),
          main = paste(tag, norm_tag, "\n", "replicate Spearman"))
    axis(1, at = seq_len(n), labels = labs, las = 2, cex.axis = 0.7)
    axis(2, at = seq_len(n), labels = rev(labs), las = 2, cex.axis = 0.7)
    for (i in seq_len(n)) for (j in seq_len(n))
      text(j, n - i + 1, sprintf("%.2f", M[i, j]), cex = 0.6)
  }, paste0(base, "_repcor"))
}

# -----------------------------------------------------------------------------
# Import + union peak file(s) into a single GRanges (rtracklayer handles
# BED/narrowPeak/broadPeak; SEACR .bed is plain BED). Returns GRanges or NULL.
# -----------------------------------------------------------------------------
load_peaks <- function(opt) {
  if (is.null(opt$peaks)) return(NULL)
  library(rtracklayer); library(GenomicRanges)
  files <- trimws(strsplit(opt$peaks, ",")[[1]])
  grs <- lapply(files, function(f) {
    g <- tryCatch(import(f), error = function(e)
      stop("Could not import peaks '", f, "': ", conditionMessage(e)))
    granges(g)
  })
  peaks <- reduce(do.call(c, grs))                     # union, merge overlaps
  message("[INFO] Peaks: ", length(peaks), " merged intervals from ",
          length(files), " file(s).")
  peaks
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
main <- function() {
  opt <- parse_cli()
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

  ss       <- load_samples(opt)
  baseline <- load_baseline_group(opt)
  peaks_gr <- load_peaks(opt)

  if (opt$test == "contrast") {
    tag <- paste0(opt$mark, "_",
                  gsub("[|]", "-", opt$group1), "_vs_", gsub("[|]", "-", opt$group2))
  } else {
    tag <- paste0(opt$mark, "_interaction_",
                  opt$cond_alt, "x", opt$treat_alt, "_vs_",
                  opt$cond_ref, "x", opt$treat_ref)
  }

  counts      <- count_and_filter_windows(ss$bamReads, ss$SampleID, opt, baseline = baseline)
  design_info <- build_design(ss, opt)

  run_one_normalization(counts$win, counts$bins, ss, design_info,
                        "libsize", opt, tag, spikein_df = NULL,
                        win_base = counts$win_base, peaks_gr = peaks_gr)

  if (!is.null(opt$spikein)) {
    if (tolower(opt$spikein) == "sheet") {
      if (!"spikein_reads" %in% colnames(ss))
        stop("--spikein sheet: the sample sheet has no 'spikein_reads' column. ",
             "Regenerate it with the updated generate_manifest.R + build_sample_sheet.R.")
      sf <- data.frame(SampleID = ss$SampleID,
                       spikein_reads = suppressWarnings(as.numeric(ss$spikein_reads)),
                       stringsAsFactors = FALSE)
      message("[INFO] Spike-in read from the sample sheet's spikein_reads column.")
    } else {
      sf <- read.csv(opt$spikein, stringsAsFactors = FALSE)
      if (!all(c("SampleID", "spikein_reads") %in% colnames(sf)))
        stop("--spikein CSV must have columns: SampleID, spikein_reads ",
             "(absolute E. coli read counts; build with build_spikein_csv.py).")
    }
    run_one_normalization(counts$win, counts$bins, ss, design_info,
                          "spikein", opt, tag, spikein_df = sf,
                          win_base = counts$win_base, peaks_gr = peaks_gr)
  } else {
    message("[NOTE] --spikein not supplied; library-size only. ",
            "Global-change results for broad repressive marks are provisional.")
  }

  message("[INFO] csaw complete: ", tag, " -> ", opt$outdir)
}

main()