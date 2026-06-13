# =============================================================================
# run_chipseqspikeinfree.R
#   Spike-in-FREE global scaling factors for CUT&RUN/ChIP histone marks, driven
#   straight from the project sample sheet, with output formatted to drop into
#   csaw.R via --normfactors.
# =============================================================================
#
# PROJECT: A375 melanoma parental-vs-B10 constricted-migration multi-omics.
#
# WHAT THIS DOES --------------------------------------------------------------
#   ChIPseqSpikeInFree (Jin et al. 2019, Bioinformatics 36:1270) estimates a
#   per-sample scaling factor that captures GLOBAL changes in a histone mark
#   (e.g. the genome-wide H3K27me3 loss under UNC1999) WITHOUT an exogenous
#   spike-in. It reads the cumulative read-distribution curve: when a mark is
#   globally eroded, a larger share of the remaining reads concentrates in the
#   few retained high-enrichment regions, and the slope of that curve scales
#   with total mark level. Larger SF = greater global loss.
#
#   This wrapper:
#     1. reads your sample sheet (the SAME csv csaw.R uses),
#     2. builds the 3-column metaFile (ID, ANTIBODY, GROUP) ChIPseqSpikeInFree
#        needs, one run per mark,
#     3. runs ChIPseqSpikeInFree (hg19 chrom sizes are built into the package),
#     4. writes the native <prefix>_SF.txt PLUS a tidy <mark>_normfactors.csv
#        keyed by your SampleID, ready for:  csaw.R --normfactors <mark>_normfactors.csv
#
# IT IS NOT A PYTHON TOOL -----------------------------------------------------
#   ChIPseqSpikeInFree is an R/Bioconductor package. There is no Python CLI,
#   no '-m' flag. Install once:
#     if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
#     BiocManager::install(c("Rsamtools","GenomicRanges","GenomicAlignments"))
#     install.packages("devtools"); devtools::install_github("stjude/ChIPseqSpikeInFree")
#
# PREREQUISITES (per the tool's README) ---------------------------------------
#   * BAMs uniquely-mapped + duplicate-removed (samtools -q 1 ; Picard MarkDuplicates
#     REMOVE_DUPLICATES=true) and INDEXED (.bai beside each .bam).
#   * Any exogenous spike-in reads already removed from the BAMs (you have none).
#   * Genome hg19 (built in). For other genomes pass --chrom a chrom.sizes path.
#
# CRITICAL CAVEAT (do not skip) -----------------------------------------------
#   The tool will "find" a global change if you ask it to. Trust the SF only with
#   independent biological evidence of the global drop -- a western/dot blot for
#   total H3K27me3 (UNC1999 lowers the H3K27me3 band even though total H3 is
#   unchanged). And note: the SF cannot separate a true biological global loss
#   from a batch-induced one, so under a treatment-confounded-with-batch design it
#   anchors MAGNITUDE (with the western), not attribution.
#
# BROAD-MARK SETTINGS ---------------------------------------------------------
#   For broad/diffuse or unideally-enriched marks (H3K27me3, H3K9me3) the README
#   recommends cutoff_QC=1 and maxLastTurn=0.97 (the defaults below).
#
# USAGE -----------------------------------------------------------------------
#   Rscript run_chipseqspikeinfree.R \
#       --sheet  data/processed/sample_sheet.csv \
#       --mark   H3K27me3 \
#       --outdir results/exploratory/spikefree
#   # all marks in the sheet at once:
#   Rscript run_chipseqspikeinfree.R --sheet data/processed/sample_sheet.csv --outdir results/exploratory/spikefree
#   # then feed csaw.R (anchored by your western):
#   Rscript csaw.R --sheet ss.csv --mark H3K27me3 --test interaction \
#       --normfactors results/exploratory/spikefree/H3K27me3_normfactors.csv \
#       --force_norm --peaks sheet --outdir results/diff/H3K27me3_interaction
#
# OUTPUTS (into --outdir) -----------------------------------------------------
#   <mark>_meta.txt           the metaFile handed to the tool (ID, ANTIBODY, GROUP)
#   <mark>_SF.txt             native ChIPseqSpikeInFree result (ID,GROUP,ANTIBODY,COLOR,QC,SF,TURNS)
#   <mark>_distribution.pdf   cumulative-curve diagnostic
#   <mark>_boxplot.pdf        SF boxplot
#   <mark>_normfactors.csv    SampleID,SF,QC,GROUP  <-- feed to csaw.R --normfactors
# =============================================================================

suppressWarnings(suppressMessages({
  library(optparse)
}))

# -----------------------------------------------------------------------------
parse_cli <- function() {
  option_list <- list(
    make_option("--sheet", type = "character",
                help = "REQUIRED. Sample sheet CSV (same one csaw.R uses). Needs: SampleID, Factor, Condition, Treatment, bamReads."),
    make_option("--mark", type = "character", default = NULL,
                help = "Factor/antibody to process (e.g. H3K27me3). Default: every distinct Factor in the sheet."),
    make_option("--outdir", type = "character", default = "spikefree",
                help = "Output directory. [default %default]"),
    make_option("--chrom", type = "character", default = "hg19",
                help = "Genome for chrom sizes. 'hg19'/'mm9'/'mm10'/'hg38' are built in; else a chrom.sizes path. [default %default]"),
    make_option("--cutoff_qc", type = "double", default = 1.0,
                help = "ChIPseqSpikeInFree cutoff_QC (broad/unideal marks: 1). [default %default]"),
    make_option("--maxlastturn", type = "double", default = 0.97,
                help = "ChIPseqSpikeInFree maxLastTurn (broad marks: 0.97). [default %default]")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$sheet)) stop("--sheet is required.")
  opt
}

# -----------------------------------------------------------------------------
load_sheet <- function(path) {
  ss <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("SampleID", "Factor", "Condition", "Treatment", "bamReads")
  miss <- setdiff(need, colnames(ss))
  if (length(miss))
    stop("Sample sheet missing column(s): ", paste(miss, collapse = ", "),
         "\n  Has: ", paste(colnames(ss), collapse = ", "))
  ss
}

# Run ChIPseqSpikeInFree on one mark and write a csaw-ready normfactors.csv.
run_one_mark <- function(mark, ss, opt) {
  sm <- ss[tolower(trimws(ss$Factor)) == tolower(trimws(mark)), , drop = FALSE]
  if (!nrow(sm)) { message("[WARN] no samples for mark ", mark, "; skipping."); return(invisible(NULL)) }
  message("\n[INFO] === ", mark, " : ", nrow(sm), " samples ===")

  bams <- sm$bamReads
  missing <- bams[!file.exists(bams)]
  if (length(missing)) stop("BAM(s) not found for ", mark, ": ", paste(missing, collapse = ", "))
  noidx <- bams[!file.exists(paste0(bams, ".bai")) & !file.exists(sub("\\.bam$", ".bai", bams))]
  if (length(noidx))
    message("[WARN] no .bai index found for: ", paste(basename(noidx), collapse = ", "),
            "\n        ChIPseqSpikeInFree needs indexed BAMs (samtools index <bam>).")

  # metaFile: ID must EXACTLY match the bamFiles strings passed below.
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

  # ChIPseqSpikeInFree's CountingReads() calls assay()/summarizeOverlaps() WITHOUT
  # importing those generics in its NAMESPACE -- it assumes they are attached. Under
  # Rscript they are not, giving 'could not find function "assay"'. Attach the
  # Bioconductor classes it implicitly relies on so the generics resolve.
  .need <- c("SummarizedExperiment", "GenomicAlignments", "Rsamtools")
  .miss <- .need[!vapply(.need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(.miss))
    stop("Missing Bioconductor dependency: ", paste(.miss, collapse = ", "),
         "\n  Install with:  BiocManager::install(c(",
         paste(sprintf('\"%s\"', .miss), collapse = ", "), "))",
         "\n  (or:  conda install -c bioconda ",
         paste0("bioconductor-", tolower(.miss), collapse = " "), ")")
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

  # Reformat the native <prefix>_SF.txt into a csaw-ready normfactors.csv.
  sf_path <- paste0(prefix, "_SF.txt")
  if (!file.exists(sf_path)) {
    message("[WARN] expected ", sf_path, " not found; skipping normfactors.csv.")
    return(invisible(NULL))
  }
  sf <- read.table(sf_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   fill = TRUE, quote = "", check.names = FALSE)
  # map the SF table's ID (== the bam strings we passed) back to SampleID
  idx <- match(sf$ID, sm$bamReads)
  if (any(is.na(idx))) idx <- match(basename(sf$ID), basename(sm$bamReads))
  if (any(is.na(idx)))
    stop("Could not map _SF.txt IDs back to SampleIDs for ", mark,
         " (IDs: ", paste(head(sf$ID, 5), collapse = ", "), ").")
  out <- data.frame(
    SampleID = sm$SampleID[idx],
    SF       = suppressWarnings(as.numeric(sf$SF)),
    QC       = if ("QC" %in% colnames(sf)) sf$QC else NA,
    GROUP    = if ("GROUP" %in% colnames(sf)) sf$GROUP else meta$GROUP[idx],
    stringsAsFactors = FALSE
  )
  nf_path <- file.path(opt$outdir, paste0(mark, "_normfactors.csv"))
  write.csv(out, nf_path, row.names = FALSE)
  message("[INFO] csaw-ready factors -> ", nf_path)
  message("[INFO] scaling factors (larger = greater global loss):")
  print(out[, c("SampleID", "GROUP", "QC", "SF")], row.names = FALSE)
  if (any(!grepl("pass", tolower(as.character(out$QC)))))
    message("[WARN] some samples did NOT pass QC (poor enrichment / complete loss) -- ",
            "inspect before injecting into csaw.")
  message("[NOTE] Trust these only WITH a western/dot blot confirming the global change.")
  invisible(out)
}

# -----------------------------------------------------------------------------
main <- function() {
  opt <- parse_cli()
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
  ss <- load_sheet(opt$sheet)
  marks <- if (is.null(opt$mark)) sort(unique(trimws(ss$Factor)))
           else trimws(strsplit(opt$mark, ",")[[1]])
  message("[INFO] marks to process: ", paste(marks, collapse = ", "))
  for (mk in marks) run_one_mark(mk, ss, opt)
  message("\n[INFO] Done -> ", normalizePath(opt$outdir))
}

main()