# =============================================================================
# diffbind.R — Peak-based differential binding for PUNCTATE histone marks
# =============================================================================
#
# PROJECT: A375 melanoma parental-vs-B10 constricted-migration multi-omics.
#
# WHAT THIS SCRIPT IS FOR ------------------------------------------------------
#   H3K27Ac forms well-defined peaks, so a peak-based differential method is
#   ideal: DiffBind builds a consensus peakset from replicates, counts reads in
#   those peaks, and tests for differential binding with DESeq2's negative-
#   binomial model. (Broad marks H3K27me3/H3K9me3/H3K4me1 use csaw.R instead,
#   because peak calling is unreliable for diffuse domains.)
#
# WHY A STATISTICAL TEST (not a difference of means) --------------------------
#   Replicates provide within-group variance; DESeq2 divides the effect by that
#   variance to yield an FDR-controlled p-value per peak. A difference of
#   averaged signal cannot do this and is not defensible as a "change" claim.
#
# IgG ROLE --------------------------------------------------------------------
#   IgG is NOT a term in the contrast. With --use_greylist it is used only to
#   build a greylist of antibody-independent artifact regions that are EXCLUDED
#   before testing. The contrast itself is a clean 2-vs-2 (or n-vs-m).
#
# CONTRASTS (explicit, over Condition|Treatment; CASE-INSENSITIVE) ------------
#   Migration:  --group1 "B10|none"  --group2 "Parental|none"
#   group1 = numerator (positive Fold => higher in group1).
#
# BATCH -----------------------------------------------------------------------
#   If >1 batch is present among the contrast samples, the script attempts a
#   ~Batch + Condition design; if that fails (e.g. perfect confounding) it warns
#   and falls back to an unblocked contrast (result is then batch-confounded).
#
# SPIKE-IN (future-proof) -----------------------------------------------------
#   --spikein factors.csv (columns: SampleID,factor). Until supplied, library-
#   size normalization runs; global-change results are provisional.
#
# EXAMPLE USAGE ---------------------------------------------------------------
#   Rscript diffbind.R \
#       --sheet data/processed/sample_sheet.csv \
#       --outdir results/diff_analysis/diffbind \
#       --mark H3K27ac \
#       --group1 "B10|none" --group2 "Parental|none"
#
#   # with IgG greylist and spike-in factors:
#   Rscript diffbind.R --sheet ss.csv --mark H3K27ac \
#       --group1 "B10|none" --group2 "Parental|none" \
#       --use_greylist --spikein spikein_factors.csv
#
# OUTPUTS (into --outdir) -----------------------------------------------------
#   sheet_<tag>.csv                 the per-mark DiffBind sample sheet used
#   <tag>_<norm>_all_sites.csv      all consensus peaks: coords, Fold(log2), FDR
#   <tag>_<norm>_significant.csv    subset passing FDR & |log2FC| thresholds
#   <tag>_<norm>_QC.pdf             heatmap, PCA, MA, volcano
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
    make_option("--mark",    type = "character", default = "H3K27ac",
                help = "Punctate Factor to analyze [default %default]"),
    make_option("--group1",  type = "character", default = "B10|none",
                help = "Numerator group Condition|Treatment [default %default]"),
    make_option("--group2",  type = "character", default = "Parental|none",
                help = "Denominator/reference group [default %default]"),
    make_option("--outdir",  type = "character", default = "diffbind_out"),
    make_option("--minOverlap", type = "integer", default = 2L,
                help = "Min samples a peak must appear in for the consensus [default %default]"),
    make_option("--summits", type = "integer", default = 200L,
                help = "Recenter peaks +/- this many bp around summit [default %default]"),
    make_option("--fdr",     type = "double",  default = 0.05),
    make_option("--lfc",     type = "double",  default = 1.0,
                help = "Min |log2FC| for the significant subset [default %default]"),
    make_option("--use_greylist", action = "store_true", default = FALSE,
                help = "Build & apply an IgG greylist (artifact exclusion only)"),
    make_option("--spikein", type = "character", default = NULL,
                help = "Optional CSV (SampleID,factor) of spike-in scale factors")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$sheet)) stop("--sheet is required.")
  opt
}

# -----------------------------------------------------------------------------
# Parse "Condition|Treatment" into lowercased parts (=> case-insensitive match).
# -----------------------------------------------------------------------------
parse_group <- function(spec) {
  parts <- strsplit(spec, "\\|")[[1]]
  if (length(parts) != 2)
    stop("Group spec must be 'Condition|Treatment', got: ", spec)
  list(cond = tolower(trimws(parts[1])), treat = tolower(trimws(parts[2])))
}

# -----------------------------------------------------------------------------
# Load sheet, subset to the mark and the two contrast groups (case-insensitive),
# label numerator/reference, and guard against <2 replicates per group.
# -----------------------------------------------------------------------------
load_contrast_samples <- function(sheet_path, mark, g1, g2) {
  full <- read.csv(sheet_path, stringsAsFactors = FALSE)
  ss <- full[tolower(full$Factor) == tolower(mark), ]
  if (nrow(ss) == 0) stop("No samples with Factor == ", mark, " in the sheet.")

  in_group <- function(rows, g)
    tolower(rows$Condition) == g$cond & tolower(rows$Treatment) == g$treat

  ss <- ss[in_group(ss, g1) | in_group(ss, g2), ]
  if (nrow(ss) == 0)
    stop("No samples matched the requested groups. Check Condition|Treatment values.")

  ss$grp <- ifelse(in_group(ss, g1), "g1", "g2")
  message("[INFO] ", mark, ": ", nrow(ss), " samples (",
          "g1=", sum(ss$grp == "g1"), ", g2=", sum(ss$grp == "g2"), ")")
  if (any(table(ss$grp) < 2))
    stop("A group has <2 replicates — DESeq2 cannot estimate dispersion. ",
         "Check --group1/--group2 casing and that both groups exist for this mark.")
  ss
}

# -----------------------------------------------------------------------------
# DiffBind needs columns named Condition / Treatment / etc. We hand it a sheet
# with a Condition column set to the contrast labels so its own contrast/design
# machinery works. We write that temp sheet and return its path.
# -----------------------------------------------------------------------------
write_diffbind_sheet <- function(ss, outdir, tag) {
  # DiffBind reads: SampleID, Condition, Treatment, Factor, Replicate,
  # bamReads, ControlID, bamControl, Peaks, PeakCaller. The sheet already has
  # these. We add a 'grp' column for clarity but DiffBind keys on Condition.
  path <- file.path(outdir, paste0("sheet_", tag, ".csv"))
  write.csv(ss, path, row.names = FALSE)
  path
}

# -----------------------------------------------------------------------------
# Build the consensus-peak DiffBind object and count reads (summit-recentered,
# paired-end). Optionally apply an IgG greylist first.
# -----------------------------------------------------------------------------
build_and_count <- function(sheet_path, ss, opt) {
  library(DiffBind)
  obj <- dba(sampleSheet = sheet_path, minOverlap = opt$minOverlap)
  message("[INFO] Consensus peakset built (minOverlap=", opt$minOverlap, ").")

  if (opt$use_greylist) {
    if (all(is.na(ss$bamControl))) {
      warning("No bamControl in sheet — skipping greylist.")
    } else {
      message("[INFO] Building IgG greylist (artifact regions excluded)...")
      obj <- dba.blacklist(obj, blacklist = FALSE, greylist = TRUE)
    }
  }

  message("[INFO] Counting reads in consensus peaks (summits=", opt$summits, ")...")
  dba.count(obj, summits = opt$summits, bParallel = TRUE)
}

# -----------------------------------------------------------------------------
# Run one normalization scheme: normalize -> contrast (blocked on Batch if >1
# batch) -> DESeq2 analyze -> save full + significant tables + QC plots.
# -----------------------------------------------------------------------------
run_one_normalization <- function(obj, ss, norm_tag, opt, tag, g1, g2,
                                   spikein_df = NULL) {
  library(DiffBind)
  message("[INFO] === Normalization: ", norm_tag, " ===")

  # Normalize: spike-in factors (user) or DiffBind library-size default.
  if (!is.null(spikein_df)) {
    ids  <- dba.show(obj)$ID
    facs <- spikein_df$factor[match(ids, spikein_df$SampleID)]
    if (any(is.na(facs)))
      stop("Spike-in file missing factors for: ",
           paste(ids[is.na(facs)], collapse = ", "))
    obj <- dba.normalize(obj, normalize = facs)
  } else {
    obj <- dba.normalize(obj, normalize = DBA_NORM_LIB)
  }

  # Contrast. Block on Batch only if >1 batch among contrast samples.
  multi_batch <- length(unique(ss$Batch)) > 1
  if (multi_batch) {
    message("[INFO] >1 batch among contrast samples -> attempting ~Batch + Condition.")
    cond_ref <- g2$cond                       # denominator condition = reference
    obj <- tryCatch(
      dba.contrast(obj, design = "~Batch + Condition",
                   reorderMeta = list(Condition = cond_ref)),
      error = function(e) {
        message("[WARN] ~Batch+Condition failed (", conditionMessage(e),
                "); falling back to unblocked contrast (result is batch-confounded).")
        dba.contrast(obj,
                     group1 = obj$masks[[opt$group1]],
                     group2 = obj$masks[[opt$group2]],
                     name1 = opt$group1, name2 = opt$group2)
      })
  } else {
    obj <- dba.contrast(obj,
                        group1 = obj$masks[[opt$group1]],
                        group2 = obj$masks[[opt$group2]],
                        name1 = opt$group1, name2 = opt$group2)
  }

  obj <- dba.analyze(obj, method = DBA_DESEQ2)

  res <- as.data.frame(dba.report(obj, th = 1, bCounts = TRUE))  # th=1 => all sites
  base <- file.path(opt$outdir, paste0(tag, "_", norm_tag))
  write.csv(res, paste0(base, "_all_sites.csv"), row.names = FALSE)

  sig <- res[!is.na(res$FDR) & res$FDR < opt$fdr & abs(res$Fold) >= opt$lfc, ]
  write.csv(sig, paste0(base, "_significant.csv"), row.names = FALSE)
  message("[INFO] ", norm_tag, ": ", nrow(sig), " significant sites of ", nrow(res),
          " (FDR<", opt$fdr, " & |log2FC|>=", opt$lfc, ")")

  save_qc_plots(obj, base)
  invisible(obj)
}

# -----------------------------------------------------------------------------
# Standard DiffBind diagnostic plots for the supplement, each written as BOTH
# a .pdf (vector, for the paper) and a .png (fast on-screen preview).
# render_one() opens both devices around a single plotting call.
# -----------------------------------------------------------------------------
render_one <- function(draw_fn, base_path, w = 7, h = 6) {
  tryCatch({ pdf(paste0(base_path, ".pdf"), width = w, height = h); draw_fn() },
           error = function(e) message("  [QC] PDF skipped: ", conditionMessage(e)),
           finally = if (!is.null(dev.list())) dev.off())
  tryCatch({ png(paste0(base_path, ".png"), width = w*150, height = h*150, res = 150); draw_fn() },
           error = function(e) message("  [QC] PNG skipped: ", conditionMessage(e)),
           finally = if (!is.null(dev.list())) dev.off())
}

save_qc_plots <- function(obj, base) {
  library(DiffBind)
  render_one(function() dba.plotHeatmap(obj),               paste0(base, "_heatmap"))
  render_one(function() dba.plotPCA(obj, label = DBA_CONDITION), paste0(base, "_PCA"))
  render_one(function() dba.plotMA(obj),                    paste0(base, "_MA"))
  render_one(function() dba.plotVolcano(obj),               paste0(base, "_volcano"))
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

  tag <- paste0(opt$mark, "_",
                gsub("[|]", "-", opt$group1), "_vs_",
                gsub("[|]", "-", opt$group2))

  sheet_path <- write_diffbind_sheet(ss, opt$outdir, tag)
  obj        <- build_and_count(sheet_path, ss, opt)

  run_one_normalization(obj, ss, "libsize", opt, tag, g1, g2, spikein_df = NULL)

  if (!is.null(opt$spikein)) {
    sf <- read.csv(opt$spikein, stringsAsFactors = FALSE)
    if (!all(c("SampleID", "factor") %in% colnames(sf)))
      stop("--spikein CSV must have columns: SampleID, factor")
    run_one_normalization(obj, ss, "spikein", opt, tag, g1, g2, spikein_df = sf)
  } else {
    message("[NOTE] --spikein not supplied; library-size only (global shifts provisional).")
  }

  message("[INFO] DiffBind complete: ", tag, " -> ", opt$outdir)
}

main()