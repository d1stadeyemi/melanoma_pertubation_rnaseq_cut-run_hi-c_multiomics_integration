#!/usr/bin/env Rscript
# =============================================================================
# run_chipseqspikeinfree.R  (v2)
#   Spike-in-FREE global scaling factors for CUT&RUN/ChIP histone marks,
#   driven from the project sample sheet, output formatted for csaw.R --normfactors.
#
# PROJECT: A375 melanoma parental-vs-B10 constricted-migration multi-omics.
# =============================================================================
#
# WHAT THIS DOES --------------------------------------------------------------
#   ChIPseqSpikeInFree (Jin et al. 2019) gives a per-sample scaling factor that
#   captures GLOBAL change in a histone mark WITHOUT exogenous spike-in. For each
#   sample it builds the cumulative read-distribution curve and reduces it to a
#   slope (beta) between two landmarks: the signal-start point (Xa,Ya) and the
#   last summit (Xb,Yb), beta = (Yb-Ya)/(Xb-Xa). SF = beta_ref / beta_sample.
#   A globally eroded mark concentrates its remaining reads into a few retained
#   cores -> last summit shifts to very high CPMW -> flat slope -> large SF.
#   LARGER SF = GREATER GLOBAL LOSS.
#
# WHAT'S NEW IN v2 (vs v1) ----------------------------------------------------
#   1. ONE POOL, ALWAYS. All samples of a mark are scaled together so every SF
#      sits on a single common ruler. This is required for any cross-condition
#      comparison (e.g. the Parental-vs-B10 interaction). Splitting conditions
#      into separate runs would give them different references = incomparable SFs.
#
#   2. POST-HOC --reference RE-ANCHOR. The per-sample slope beta is POOL-
#      INDEPENDENT (it comes only from that sample's own curve). Only the
#      reference (the numerator) is a pool-level choice, and it is a global
#      constant that CANCELS in any edgeR contrast. So we let the tool run, then
#      recompute beta from the TURNS column and re-anchor SF = beta_ref/beta_i to
#      a reference YOU choose. Default reference = the steepest-slope UNTREATED,
#      QC-passing sample (cleanest baseline), so treated samples read as SF>1.
#      This does NOT change DB results; it only makes the SF scale interpretable
#      and avoids the tool's habit of auto-anchoring to whatever sample is
#      steepest (which in this dataset was a treated sample).
#
#   3. SUBSET FLAGS (--condition, --treatment) with a guardrail: a run that
#      contains only ONE treatment level cannot show a global treatment change,
#      so it is refused unless --force. (Recommended invocation subsets NEITHER.)
#
#   4. OPTIONAL --qc ARTIFACT PASS. A flat slope (large SF) can be real biology
#      OR a handful of junk windows (blacklist/repeat pileups, residual dups).
#      ChIPseqSpikeInFree's own "pass" QC checks curve SHAPE, not alignment
#      integrity, so it cannot catch this. With --qc we read the tool's
#      <mark>_rawCounts.txt, measure how concentrated each sample's reads are in
#      its top 1-kb bins, and (with --blacklist) what fraction sit in blacklisted
#      regions. A sample whose large SF is carried by a few blacklist bins is an
#      alignment artifact, not a global loss.
#
# CRITICAL CAVEAT (do not skip) -----------------------------------------------
#   Re-anchoring fixes the SCALE, not the DATA. It cannot un-inflate an artifact-
#   driven slope. Trust these SFs only WITH (a) clean alignments and (b) an
#   orthogonal protein-level confirmation of the global change (western/dot blot
#   for total H3K27me3; UNC1999 lowers the band while total H3 is unchanged).
#
# INSTALL (once) --------------------------------------------------------------
#   if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
#   BiocManager::install(c("Rsamtools","GenomicRanges","GenomicAlignments","SummarizedExperiment"))
#   install.packages(c("optparse","devtools")); devtools::install_github("stjude/ChIPseqSpikeInFree")
#
# USAGE -----------------------------------------------------------------------
#   # recommended: all H3K27me3 samples, one pool, auto untreated reference
#   Rscript run_chipseqspikeinfree.R \
#       --sheet data/processed/sample_sheet.csv \
#       --mark  H3K27me3 \
#       --outdir results/exploratory/spikefree
#
#   # pin the reference to a specific clean untreated replicate, and run the QC pass
#   Rscript run_chipseqspikeinfree.R \
#       --sheet data/processed/sample_sheet.csv --mark H3K27me3 \
#       --reference Parental-H3K27me3-R3 \
#       --qc --blacklist refs/hg19-blacklist.v2.bed \
#       --outdir results/exploratory/spikefree
#
#   # then feed csaw.R (only after alignments are clean + western in hand):
#   Rscript csaw.R --sheet ss.csv --mark H3K27me3 --test interaction \
#       --normfactors results/exploratory/spikefree/H3K27me3_normfactors.csv \
#       --force_norm --peaks sheet --outdir results/diff/H3K27me3_interaction
#
# OUTPUTS (into --outdir) -----------------------------------------------------
#   <mark>_meta.txt          metaFile handed to the tool (ID, ANTIBODY, GROUP)
#   <mark>_SF.txt            native ChIPseqSpikeInFree result
#   <mark>_distribution.pdf  cumulative-curve diagnostic
#   <mark>_boxplot.pdf       SF boxplot
#   <mark>_rawCounts.txt     per-1kb-bin counts (used by --qc)
#   <mark>_normfactors.csv   SampleID, SF (re-anchored), beta, SF_tool, QC, GROUP,
#                            reference, is_reference   <-- feed to csaw.R
#   <mark>_artifactQC.csv    (with --qc) concentration + blacklist metrics + flag
# =============================================================================

suppressWarnings(suppressMessages({
  library(optparse)
}))

# ---------------------------------------------------------------------------- #
#  CLI                                                                         #
# ---------------------------------------------------------------------------- #
parse_cli <- function() {
  option_list <- list(
    make_option("--sheet", type = "character",
                help = "REQUIRED. Sample sheet CSV (the same one csaw.R uses). Needs: SampleID, Factor, Condition, Treatment, bamReads."),
    make_option("--mark", type = "character", default = NULL,
                help = "Factor/antibody to process, e.g. H3K27me3 (comma-list ok). Default: every distinct Factor in the sheet."),
    make_option("--condition", type = "character", default = NULL,
                help = "Optional Condition subset, e.g. parental or B10 (comma-list ok). Default: all conditions in one pool (recommended)."),
    make_option("--treatment", type = "character", default = NULL,
                help = "Optional Treatment subset, e.g. none,unc (comma-list ok). Default: all treatments. A single-treatment pool is refused (see --force)."),
    make_option("--reference", type = "character", default = NULL,
                help = "SampleID to use as the SF reference anchor. Default: steepest-slope UNTREATED QC-passing sample."),
    make_option("--outdir", type = "character", default = "spikefree",
                help = "Output directory. [default %default]"),
    make_option("--chrom", type = "character", default = "hg19",
                help = "Genome for chrom sizes. hg19/mm9/mm10/hg38 are built in; else a chrom.sizes path. [default %default]"),
    make_option("--cutoff_qc", type = "double", default = 1.0,
                help = "ChIPseqSpikeInFree cutoff_QC (broad/unideal marks: 1). [default %default]"),
    make_option("--maxlastturn", type = "double", default = 0.97,
                help = "ChIPseqSpikeInFree maxLastTurn (broad marks: 0.97). [default %default]"),
    make_option("--qc", action = "store_true", default = FALSE,
                help = "Run the artifact QC pass (top-bin read concentration; with --blacklist, blacklist read fraction)."),
    make_option("--blacklist", type = "character", default = NULL,
                help = "BED of blacklisted/repeat regions (e.g. ENCODE hg19) for the --qc pass. Optional."),
    make_option("--qc_top_frac", type = "double", default = 0.05,
                help = "Flag a sample if its top-50 1-kb bins hold more than this read fraction. [default %default]"),
    make_option("--force", action = "store_true", default = FALSE,
                help = "Override the single-treatment-pool guardrail.")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$sheet)) stop("--sheet is required.")
  opt
}

split_csv <- function(x) if (is.null(x)) NULL else trimws(strsplit(x, ",")[[1]])

# ---------------------------------------------------------------------------- #
#  Sample sheet                                                                #
# ---------------------------------------------------------------------------- #
load_sheet <- function(path) {
  ss <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("SampleID", "Factor", "Condition", "Treatment", "bamReads")
  miss <- setdiff(need, colnames(ss))
  if (length(miss))
    stop("Sample sheet missing column(s): ", paste(miss, collapse = ", "),
         "\n  Has: ", paste(colnames(ss), collapse = ", "))
  ss
}

# Filter sheet to one mark (+ optional condition/treatment) and apply guardrails.
subset_samples <- function(ss, mark, opt) {
  keep <- tolower(trimws(ss$Factor)) == tolower(trimws(mark))
  if (!is.null(opt$condition))
    keep <- keep & tolower(trimws(ss$Condition)) %in% tolower(split_csv(opt$condition))
  if (!is.null(opt$treatment))
    keep <- keep & tolower(trimws(ss$Treatment)) %in% tolower(split_csv(opt$treatment))
  sm <- ss[keep, , drop = FALSE]
  if (!nrow(sm)) return(sm)

  # Guardrail: a pool with only one treatment level cannot reveal a global
  # treatment change. Refuse unless --force.
  ntreat <- length(unique(tolower(trimws(sm$Treatment))))
  if (ntreat < 2 && !isTRUE(opt$force))
    stop("Only one treatment level in the ", mark, " pool (",
         paste(sort(unique(sm$Treatment)), collapse = ", "), ").\n",
         "  Spike-in-free SFs are only meaningful when treated AND untreated are scaled together.\n",
         "  Include both treatments, or pass --force if you really intend this.")

  # Soft warning: thin groups weaken the downstream DB test.
  tab <- table(paste(trimws(sm$Condition), trimws(sm$Treatment), sep = "_"))
  thin <- names(tab)[tab < 2]
  if (length(thin))
    message("[WARN] group(s) with <2 replicates (weak for DB later): ",
            paste(thin, collapse = ", "))
  sm
}

# ---------------------------------------------------------------------------- #
#  Run the tool                                                                #
# ---------------------------------------------------------------------------- #
run_spikefree <- function(mark, sm, opt) {
  bams <- sm$bamReads
  missing <- bams[!file.exists(bams)]
  if (length(missing)) stop("BAM(s) not found for ", mark, ": ", paste(missing, collapse = ", "))
  noidx <- bams[!file.exists(paste0(bams, ".bai")) & !file.exists(sub("\\.bam$", ".bai", bams))]
  if (length(noidx))
    message("[WARN] no .bai index found for: ", paste(basename(noidx), collapse = ", "),
            "\n        ChIPseqSpikeInFree needs indexed BAMs (samtools index <bam>).")

  meta <- data.frame(
    ID       = basename(bams),
    ANTIBODY = mark,
    GROUP    = paste(trimws(sm$Condition), trimws(sm$Treatment), sep = "_"),
    stringsAsFactors = FALSE
  )
  meta_path <- file.path(opt$outdir, paste0(mark, "_meta.txt"))
  write.table(meta, meta_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("[INFO] metaFile -> ", meta_path)
  message("[INFO] groups: ", paste(sort(unique(meta$GROUP)), collapse = ", "))

  if (!requireNamespace("ChIPseqSpikeInFree", quietly = TRUE))
    stop("Package 'ChIPseqSpikeInFree' is not installed.\n",
         "  devtools::install_github(\"stjude/ChIPseqSpikeInFree\")")

  # ChIPseqSpikeInFree calls assay()/summarizeOverlaps() without importing the
  # generics; under Rscript they are not attached -> 'could not find function'.
  # Attach the Bioconductor classes it implicitly relies on.
  .need <- c("SummarizedExperiment", "GenomicAlignments", "Rsamtools")
  .miss <- .need[!vapply(.need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(.miss))
    stop("Missing Bioconductor dependency: ", paste(.miss, collapse = ", "),
         "\n  BiocManager::install(c(",
         paste(sprintf('\"%s\"', .miss), collapse = ", "), "))")
  invisible(lapply(.need, function(p)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE)))))

  prefix <- file.path(opt$outdir, mark)
  message("[INFO] running ChIPseqSpikeInFree (chrom=", opt$chrom,
          ", cutoff_QC=", opt$cutoff_qc, ", maxLastTurn=", opt$maxlastturn, ") ...")
  ChIPseqSpikeInFree::ChIPseqSpikeInFree(
    bamFiles    = bams,
    chromFile   = opt$chrom,
    metaFile    = meta_path,
    prefix      = prefix,
    cutoff_QC   = opt$cutoff_qc,
    maxLastTurn = opt$maxlastturn
  )
  prefix
}

# ---------------------------------------------------------------------------- #
#  Re-anchor SF to a chosen reference                                          #
# ---------------------------------------------------------------------------- #

# beta = (Yb - Ya) / (Xb - Xa) from the TURNS string "Xa,Ya,Xb,Yb".
parse_beta <- function(turns) {
  vapply(turns, function(s) {
    if (is.na(s) || !nzchar(trimws(s))) return(NA_real_)
    v <- suppressWarnings(as.numeric(strsplit(trimws(s), ",")[[1]]))
    if (length(v) < 4 || any(is.na(v[1:4])) || v[3] == v[1]) return(NA_real_)
    (v[4] - v[2]) / (v[3] - v[1])
  }, numeric(1), USE.NAMES = FALSE)
}

# Choose the reference row index (into the mapped data.frame `df`).
choose_reference <- function(df, ref_arg) {
  ok <- !is.na(df$beta) & grepl("pass", tolower(df$QC))
  if (!is.null(ref_arg)) {
    i <- which(tolower(trimws(df$SampleID)) == tolower(trimws(ref_arg)))
    if (!length(i))
      stop("--reference '", ref_arg, "' not found. Samples: ",
           paste(df$SampleID, collapse = ", "))
    if (!ok[i])
      stop("--reference '", ref_arg, "' failed QC or has no valid slope; pick a clean sample.")
    return(i)
  }
  # default: steepest-beta untreated, QC-passing (cleanest, most-abundant baseline)
  unt  <- ok & grepl("none", tolower(df$Treatment))
  pool <- if (any(unt)) which(unt) else which(ok)
  if (!length(pool)) stop("No QC-passing sample can serve as a reference.")
  if (!any(unt))
    message("[WARN] no untreated QC-pass sample; falling back to steepest sample overall as reference.")
  pool[which.max(df$beta[pool])]
}

# Build the csaw-ready normfactors table.
process_sf <- function(prefix, sm, mark, opt) {
  sf_path <- paste0(prefix, "_SF.txt")
  if (!file.exists(sf_path)) {
    message("[WARN] ", sf_path, " not found; skipping normfactors.csv."); return(invisible(NULL))
  }
  sf <- read.table(sf_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   fill = TRUE, quote = "", check.names = FALSE)
  if (!"TURNS" %in% colnames(sf))
    stop("No TURNS column in ", sf_path, " - cannot re-anchor. Check tool version/output.")

  # map _SF.txt rows (ID == bam basename) back to sample-sheet rows
  idx <- match(basename(sf$ID), basename(sm$bamReads))
  if (any(is.na(idx)))
    stop("Could not map _SF.txt IDs back to SampleIDs (IDs: ",
         paste(head(sf$ID, 5), collapse = ", "), ").")

  df <- data.frame(
    SampleID  = sm$SampleID[idx],
    Condition = trimws(sm$Condition[idx]),
    Treatment = trimws(sm$Treatment[idx]),
    GROUP     = paste(trimws(sm$Condition[idx]), trimws(sm$Treatment[idx]), sep = "_"),
    QC        = if ("QC" %in% colnames(sf)) sf$QC else NA_character_,
    SF_tool   = suppressWarnings(as.numeric(sf$SF)),
    beta      = parse_beta(sf$TURNS),
    stringsAsFactors = FALSE
  )

  ref_i  <- choose_reference(df, opt$reference)
  ref_id <- df$SampleID[ref_i]
  ref_b  <- df$beta[ref_i]
  df$SF          <- ref_b / df$beta          # re-anchored: reference == 1, loss > 1
  df$reference   <- ref_id
  df$is_reference <- seq_len(nrow(df)) == ref_i

  out <- df[, c("SampleID", "SF", "beta", "SF_tool", "QC", "GROUP",
                "reference", "is_reference")]
  nf_path <- file.path(opt$outdir, paste0(mark, "_normfactors.csv"))
  write.csv(out, nf_path, row.names = FALSE)

  message("[INFO] reference anchor: ", ref_id, " (beta=", signif(ref_b, 4), ")")
  message("[INFO] csaw-ready factors -> ", nf_path)
  message("[INFO] re-anchored SF (1 = reference; larger = greater global loss):")
  print(out[order(out$GROUP), c("SampleID", "GROUP", "QC", "beta", "SF")], row.names = FALSE)
  if (any(!grepl("pass", tolower(as.character(out$QC)))))
    message("[WARN] some samples did NOT pass QC -- inspect before injecting into csaw.")
  message("[NOTE] Re-anchoring fixes the SCALE, not the data. Confirm the global ",
          "change with a western/dot blot before trusting magnitude.")
  invisible(out)
}

# ---------------------------------------------------------------------------- #
#  Optional artifact QC pass                                                   #
# ---------------------------------------------------------------------------- #
read_table_fast <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE))
    as.data.frame(data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE))
  else
    read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

# Is a sample's large SF carried by a few junk windows? Reads <prefix>_rawCounts.txt.
# Non-blocking: any failure here warns and returns, leaving the SFs intact.
qc_artifact_check <- function(prefix, sm, mark, opt) {
  tryCatch({
    raw_path <- paste0(prefix, "_rawCounts.txt")
    if (!file.exists(raw_path)) { message("[WARN] ", raw_path, " not found; skipping --qc."); return(invisible(NULL)) }
    raw <- read_table_fast(raw_path)
    cn  <- colnames(raw)

    # locate coordinate columns (defensive: fall back to first three columns)
    pick <- function(opts) { h <- cn[tolower(cn) %in% opts]; if (length(h)) h[1] else NA_character_ }
    chr_col <- pick(c("chr","seqnames","chrom","chromosome"))
    sta_col <- pick(c("start"))
    end_col <- pick(c("end","stop"))
    if (any(is.na(c(chr_col, sta_col, end_col)))) {
      message("[WARN] could not name coord columns in rawCounts; assuming first 3 are chr,start,end.")
      chr_col <- cn[1]; sta_col <- cn[2]; end_col <- cn[3]
    }
    samp_cols <- setdiff(cn, c(chr_col, sta_col, end_col))
    # keep only sample columns that map to this mark's BAMs
    want <- basename(sm$bamReads)
    keep <- samp_cols[basename(samp_cols) %in% want | samp_cols %in% want]
    if (!length(keep)) { message("[WARN] no rawCounts sample columns matched this mark; skipping --qc."); return(invisible(NULL)) }

    # optional blacklist overlap
    bl_hit <- NULL
    if (!is.null(opt$blacklist)) {
      if (!requireNamespace("GenomicRanges", quietly = TRUE))
        stop("GenomicRanges needed for --blacklist overlap.")
      suppressWarnings(suppressMessages(library(GenomicRanges)))
      bl <- read.table(opt$blacklist, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      bl_gr <- GenomicRanges::makeGRangesFromDataFrame(
        data.frame(seqnames = bl[[1]], start = bl[[2]] + 1L, end = bl[[3]]))  # BED is 0-based
      bins_gr <- GenomicRanges::makeGRangesFromDataFrame(
        data.frame(seqnames = raw[[chr_col]],
                   start = as.integer(raw[[sta_col]]) + 1L,
                   end   = as.integer(raw[[end_col]])))
      bl_hit <- IRanges::overlapsAny(bins_gr, bl_gr)
    }

    rows <- lapply(keep, function(s) {
      x   <- as.numeric(raw[[s]]); x[is.na(x)] <- 0
      tot <- sum(x); if (tot == 0) return(NULL)
      o   <- order(x, decreasing = TRUE)
      n01 <- max(1L, floor(length(x) * 0.001))               # top 0.1% of bins
      sid <- sm$SampleID[match(basename(s), basename(sm$bamReads))]
      data.frame(
        SampleID       = if (length(sid)) sid else s,
        total_reads    = tot,
        frac_top50     = sum(x[o[seq_len(min(50, length(o)))]]) / tot,
        frac_top0.1pct = sum(x[o[seq_len(n01)]]) / tot,
        frac_blacklist = if (is.null(bl_hit)) NA_real_ else sum(x[bl_hit]) / tot,
        top50_in_bl    = if (is.null(bl_hit)) NA_integer_ else sum(bl_hit[o[seq_len(min(50, length(o)))]]),
        stringsAsFactors = FALSE)
    })
    qc <- do.call(rbind, rows)
    qc$artifact_flag <- (qc$frac_top50 > opt$qc_top_frac) |
                        (!is.na(qc$frac_blacklist) & qc$frac_blacklist > 0.10)

    qc_path <- file.path(opt$outdir, paste0(mark, "_artifactQC.csv"))
    write.csv(qc, qc_path, row.names = FALSE)
    message("[INFO] artifact QC -> ", qc_path)
    print(qc, row.names = FALSE)
    flagged <- qc$SampleID[qc$artifact_flag]
    if (length(flagged))
      message("[WARN] artifact-suspect (slope likely driven by a few junk windows): ",
              paste(flagged, collapse = ", "),
              "\n        Inspect/realign these before using their SFs.")
    invisible(qc)
  }, error = function(e) {
    message("[WARN] --qc pass failed (SFs are unaffected): ", conditionMessage(e))
    invisible(NULL)
  })
}

# ---------------------------------------------------------------------------- #
#  Main                                                                        #
# ---------------------------------------------------------------------------- #
main <- function() {
  opt <- parse_cli()
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
  ss <- load_sheet(opt$sheet)
  marks <- if (is.null(opt$mark)) sort(unique(trimws(ss$Factor))) else split_csv(opt$mark)
  message("[INFO] marks to process: ", paste(marks, collapse = ", "))

  for (mk in marks) {
    message("\n[INFO] ================ ", mk, " ================")
    sm <- subset_samples(ss, mk, opt)
    if (!nrow(sm)) { message("[WARN] no samples for ", mk, "; skipping."); next }
    message("[INFO] ", nrow(sm), " samples in pool.")
    prefix <- run_spikefree(mk, sm, opt)
    process_sf(prefix, sm, mk, opt)
    if (isTRUE(opt$qc)) qc_artifact_check(prefix, sm, mk, opt)
  }
  message("\n[INFO] Done -> ", normalizePath(opt$outdir))
}

main()