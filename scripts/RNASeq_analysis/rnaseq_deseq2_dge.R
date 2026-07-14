#!/usr/bin/env Rscript
###############################################################################
# rnaseq_deseq2_dge.R
# -----------------------------------------------------------------------------
# Differential gene expression (DGE) for ONE two-group contrast, from htseq-count
# outputs produced by the upstream `rnaseq_fastq_to_counts.sh` pipeline, all the
# way to annotated DEG tables + a full battery of QC / diagnostic figures.
#
# The unit of analysis here is a CONTRAST (group A vs group B), not a sample.
# You therefore tell the script which replicates belong to which group, either
# with a sample sheet (recommended, reproducible) or on the command line.
#
# WHAT IT PRODUCES (everything is written to disk; nothing is silent):
#   tables/   assignment-rate QC, library sizes, detected genes, size factors,
#             filtering summary, PCA coordinates, sample distances, full &
#             significant results, threshold-sensitivity grid, annotation report
#   figures/  htseq assignment bars, count distributions (raw/norm), RLE,
#             size-factor diagnostics, dispersion plot, mean-SD (VST/rlog),
#             PCA + scree, sample-distance & correlation heatmaps, Cook's
#             distance, p-value histogram, independent-filtering curve, MA
#             (raw + shrunken), volcano, top-gene count plots, top-DEG heatmap
#   MANIFEST.txt, sessionInfo.txt, run.log   (provenance)
#
# DESeq2 reference (attached vignette): Love, Huber & Anders 2014; the model
# uses raw counts (NEVER pre-normalized values) and corrects for library size
# internally via median-of-ratios size factors.
#
# ---------------------------------------------------------------------------
# QUICK EXAMPLES
# ---------------------------------------------------------------------------
# 1) Simplest: keyword-scan the counts dir; grch38; Parental as control.
#    Every replicate whose folder name contains "bottom10"/"parental" is picked.
#    Rscript rnaseq_deseq2_dge.R \
#        --ref grch38 \
#        --counts-dir /lustre/.../rnaseq_output \
#        --group1 Bottom10 --group2 Parental \
#        --control Parental \
#        --outdir results/exploratory/deseq2_B10_vs_Parental
#    # optional: assert 3 reps each and set clean labels
#    #   --group1-reps rep1,rep2,rep3 --group2-reps rep1,rep2,rep3
#
# 2) Reproducible: drive everything from a sample sheet, add a batch covariate
#    Rscript rnaseq_deseq2_dge.R \
#        --ref grch38 \
#        --counts-dir /lustre/.../rnaseq_output \
#        --samplesheet samples.csv \
#        --control Parental --batch batch \
#        --outdir results/exploratory/deseq2_B10_vs_Parental
#
#    samples.csv (header required; 'path' optional; extra cols allowed):
#        sample,condition,batch
#        Par_r1,Parental,1
#        Par_r2,Parental,2
#        B10_r1,Bottom10,1
#        B10_r2,Bottom10,2
#
# 3) Just check the setup without running DESeq2 (validate files/design):
#    Rscript rnaseq_deseq2_dge.R --ref hg19 --counts-dir DIR \
#        --samplesheet samples.csv --control Parental --dry-run
#
# Run `Rscript rnaseq_deseq2_dge.R --help` for the full flag list.
###############################################################################

suppressWarnings(suppressMessages({
  ok_core <- requireNamespace("DESeq2", quietly = TRUE)
}))

## ===========================================================================
## 0. TINY ARG PARSER  (zero external dependencies, so it runs in any R)
## ===========================================================================
# Supports:  --key value   |   --key=value   |   --boolean-flag
# Unknown keys are rejected with a clear message. Booleans are listed below.
.bool_flags <- c("help", "h", "no-annotation", "strip-gene-version", "dry-run", "collapse-to-gene")

parse_args <- function(args) {
  opt <- list()
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (!startsWith(a, "--") && !(a %in% c("-h"))) {
      stop("Unexpected argument (not a --flag): '", a, "'. See --help.", call. = FALSE)
    }
    a2 <- sub("^--?", "", a)               # strip leading -- (or -h)
    val <- NULL
    if (grepl("=", a2, fixed = TRUE)) {
      key <- sub("=.*$", "", a2)
      val <- sub("^[^=]*=", "", a2)
    } else {
      key <- a2
    }
    if (key %in% .bool_flags) {             # boolean: no value consumed
      opt[[key]] <- TRUE
      i <- i + 1L
      next
    }
    if (is.null(val)) {                     # value is the next token
      if (i + 1L > length(args))
        stop("Option --", key, " requires a value.", call. = FALSE)
      val <- args[[i + 1L]]
      i <- i + 2L
    } else {
      i <- i + 1L
    }
    opt[[key]] <- val
  }
  opt
}

print_help <- function() {
  cat(readLines_help())
}
readLines_help <- function() {
'
rnaseq_deseq2_dge.R  —  DESeq2 differential expression for one two-group contrast

REQUIRED
  --ref <hg19|grch38>         Reference build. Selects the default annotation GTF
                              (used only to add gene symbols/biotypes to results).

INPUT — counts location
  --counts-dir <DIR>          Output dir of rnaseq_fastq_to_counts.sh. Each sample
                              is looked up as <DIR>/<sample>/<sample><suffix> then
                              <DIR>/<sample><suffix>.
  --counts-suffix <STR>       Count-file suffix. Default: _genecounts.tsv

INPUT — group assignment (choose ONE of A or B)
  A) Keyword scanning (quickest; pairs with rnaseq_fastq_to_counts.sh output):
     --group1 <KEYWORD>       Case-insensitive substring matched against the
     --group2 <KEYWORD>       replicate names discovered in --counts-dir. Every
                              count file whose name contains the keyword joins
                              that group; the keyword is also the group label.
                              (A file matching BOTH keywords is an error — make
                              them more specific. Sample sheet = surgical control.)
     --group1-reps <a,b,...>  OPTIONAL. Its LENGTH asserts how many replicates
     --group2-reps <a,b,...>  must be found (errors on mismatch); the labels
                              become the column/figure labels as <group>_<label>,
                              mapped to the discovered files in sorted order.
                              Omit to use the discovered replicate names as labels.
  B) Sample sheet (most control / multi-factor designs):
     --samplesheet <CSV/TSV>  Columns: sample,condition[,batch,path,...]

CONTRAST / DESIGN
  --control <NAME>            Reference level = denominator of the log2 fold change
                              (e.g. Parental). LFC>0 means UP in the OTHER group.
                              Default: group2 name (inline) or first level (sheet, warns).
  --batch <COL>               Add this sample-sheet column as a covariate:
                              design becomes ~ <batch> + condition.
  --design <FORMULA>          Advanced: override the full design, e.g. "~ batch + condition".
  --contrast-name <STR>       Label used in filenames/plots. Default: <treat>_vs_<control>.

TESTING PARAMETERS
  --alpha <FLOAT>             FDR threshold for calling DEGs. Default: 0.05
  --lfc-threshold <FLOAT>     Test |LFC| > threshold (Wald against a threshold). Default: 0
  --alt-hypothesis <greaterAbs|lessAbs|greater|less>  Used when lfc-threshold>0. Default: greaterAbs
  --shrink <apeglm|ashr|normal|none>   LFC shrinkage for ranking/plots. Default: apeglm
  --fit-type <parametric|local|mean>   Dispersion fit. Default: parametric

FILTERING
  --min-count <INT>           Keep genes with >= this count ... Default: 10
  --min-samples <INT>         ... in at least this many samples. Default: size of smallest group

TRANSFORMS / PLOTS
  --transform <vst|rlog|both> Transform for QC/clustering. Default: vst (rlog is slow)
  --top-n <INT>               #genes for count plots & top-DEG heatmap. Default: 30
  --biotype-filter <STR>      Keep only this gene biotype before testing (e.g. protein_coding).
                              Requires annotation. Default: none.

ANNOTATION
  --gtf <FILE>                Override annotation GTF (else chosen by --ref).
  --id-map <FILE>             CSV/TSV mapping gene_id -> gene_name[,biotype]. Overrides
                              GTF symbols. Use this to add symbols for hg19 UCSC
                              knownGene (which has none): supply a UCSC kgXref table
                              with columns  gene_id,gene_name  (kgID -> geneSymbol).
  --no-annotation             Skip gene-symbol annotation entirely.
  --strip-gene-version        Strip ".NN" version from Ensembl gene IDs (both counts & GTF).
  --collapse-to-gene          Also run a GENE-LEVEL analysis: transcript-model counts
                              are summed per gene symbol. Outputs split into
                              <outdir>/transcript_level/ and <outdir>/gene_level/.
                              (Most useful for hg19 knownGene, whose IDs are
                              transcript models; for grch38/gencode the counts are
                              already gene-level so the two levels are near-identical.)

MISC
  --outdir <DIR>              Output directory. Default: ./deseq2_out
  --threads <INT>             BiocParallel workers. Default: $SLURM_CPUS_PER_TASK or 1
  --seed <INT>                RNG seed (apeglm/ashr). Default: 1
  --rlib <DIR>                Prepend an R library path (e.g. the multiomics-r library,
                              from `Rscript -e ".libPaths()[1]"`), so annotation packages
                              are found even outside the env. NOTE: the R running this
                              script must be ABI-compatible with packages built there;
                              the safe path is to just run inside multiomics-r (the
                              submit_ wrapper activates it by absolute path).
  --dry-run                   Resolve config, validate files & design, then STOP.
  -h, --help                  Show this help.
'
}

## ===========================================================================
## 1. RESOLVE + VALIDATE CONFIG
## ===========================================================================
raw <- commandArgs(trailingOnly = TRUE)
opt <- tryCatch(parse_args(raw), error = function(e) { message("ERROR: ", conditionMessage(e)); quit(status = 2) })

if (isTRUE(opt[["help"]]) || isTRUE(opt[["h"]]) || length(raw) == 0) {
  print_help(); quit(status = 0)
}

# Reject unknown / misspelled flags (protects against silent typos like
# --controll that would otherwise fall back to a wrong default).
.known_keys <- c(
  "ref","gtf","no-annotation","strip-gene-version","id-map",
  "counts-dir","counts-suffix","samplesheet",
  "group1","group2","group1-reps","group2-reps",
  "control","batch","design","contrast-name",
  "alpha","lfc-threshold","alt-hypothesis","shrink","fit-type",
  "min-count","min-samples","transform","top-n","biotype-filter",
  "outdir","threads","seed","dry-run","help","h","rlib","collapse-to-gene")
# Friendly message for the pre-scan flag names that were renamed
for (.old in c("group1-samples","group2-samples")) if (!is.null(opt[[.old]])) {
  .grp <- sub("-samples","",.old)
  message("ERROR: --", .old, " no longer exists. Group membership is now auto-scanned from ",
          "the --", .grp, " keyword (case-insensitive substring of the replicate names). ",
          "Use --", sub("samples","reps",.old), " only to (optionally) assert the replicate ",
          "count and set figure labels.")
  quit(status = 2)
}
.unknown <- setdiff(names(opt), .known_keys)
if (length(.unknown) > 0) {
  message("ERROR: unknown option(s): ", paste0("--", .unknown, collapse = ", "),
          "\n       run with --help to see valid flags.")
  quit(status = 2)
}

# ---- helper accessors with defaults + coercion -----------------------------
getc <- function(k, default = NULL) if (!is.null(opt[[k]])) opt[[k]] else default
getn <- function(k, default) { v <- opt[[k]]; if (is.null(v)) return(default); as.numeric(v) }
geti <- function(k, default) { v <- opt[[k]]; if (is.null(v)) return(default); as.integer(v) }
getb <- function(k) isTRUE(opt[[k]])
split_csv <- function(s) if (is.null(s)) NULL else trimws(strsplit(s, ",", fixed = TRUE)[[1]])

# --ref is always required (it is the flag you specifically asked to keep) -----
ref <- getc("ref")
if (is.null(ref)) { message("ERROR: --ref <hg19|grch38> is required."); quit(status = 2) }
ref <- tolower(ref)
if (!ref %in% c("hg19", "grch38")) {
  message("ERROR: --ref must be 'hg19' or 'grch38' (got '", ref, "')."); quit(status = 2)
}

# Default annotation GTFs (match the upstream preprocessing scripts) ----------
default_gtf <- c(
  hg19   = "/lustre/isaac24/proj/UTK0009/RNA-Seq/reference_genomes/hg19/HomoSapiens_hg19.gtf",
  grch38 = "/lustre/isaac24/proj/UTK0009/RNA-Seq/reference_genomes/grch38/gencode.v49.annotation.gtf"
)
gtf_path      <- getc("gtf", default_gtf[[ref]])
id_map        <- getc("id-map")
rlib          <- getc("rlib")
if (!is.null(rlib)) {
  if (!dir.exists(rlib)) { message("ERROR: --rlib directory not found: ", rlib); quit(status = 2) }
  .libPaths(c(rlib, .libPaths()))   # prepend so packages here are found first
}
do_annotate   <- !getb("no-annotation")
strip_version <- getb("strip-gene-version")

counts_dir    <- getc("counts-dir")
counts_suffix <- getc("counts-suffix", "_genecounts.tsv")
samplesheet   <- getc("samplesheet")
outdir        <- getc("outdir", "deseq2_out")
alpha         <- getn("alpha", 0.05)
lfc_thr       <- getn("lfc-threshold", 0)
alt_hyp       <- getc("alt-hypothesis", "greaterAbs")
shrink        <- getc("shrink", "apeglm")
fit_type      <- getc("fit-type", "parametric")
min_count     <- geti("min-count", 10L)
min_samples   <- if (!is.null(opt[["min-samples"]])) geti("min-samples", NA) else NA_integer_
transform     <- getc("transform", "vst")
top_n         <- geti("top-n", 30L)
biotype_filt  <- getc("biotype-filter")
batch_col     <- getc("batch")
design_override <- getc("design")
control       <- getc("control")
contrast_name <- getc("contrast-name")
threads       <- geti("threads", { s <- Sys.getenv("SLURM_CPUS_PER_TASK"); if (nzchar(s)) as.integer(s) else 1L })
seed          <- geti("seed", 1L)
dry_run       <- getb("dry-run")
collapse_gene <- getb("collapse-to-gene")

if (is.null(counts_dir)) { message("ERROR: --counts-dir is required."); quit(status = 2) }
if (!dir.exists(counts_dir)) { message("ERROR: --counts-dir not found: ", counts_dir); quit(status = 2) }
if (!shrink %in% c("apeglm","ashr","normal","none")) { message("ERROR: --shrink invalid: ", shrink); quit(status = 2) }
if (!transform %in% c("vst","rlog","both")) { message("ERROR: --transform invalid: ", transform); quit(status = 2) }
if (!fit_type %in% c("parametric","local","mean")) { message("ERROR: --fit-type invalid: ", fit_type); quit(status = 2) }
if (!alt_hyp %in% c("greaterAbs","lessAbs","greater","less")) { message("ERROR: --alt-hypothesis invalid: ", alt_hyp); quit(status = 2) }

## ---------------------------------------------------------------------------
## 1b. BUILD THE SAMPLE TABLE (coldata) from sheet OR inline groups
## ---------------------------------------------------------------------------
# Discover every htseq count file under --counts-dir. Handles the standard
# rnaseq_fastq_to_counts.sh layout  <dir>/<rep>/<rep><suffix>  and a flat
# <dir>/<rep><suffix>. The replicate NAME is the basename minus the suffix.
discover_counts <- function(counts_dir, suffix) {
  files <- unique(c(
    Sys.glob(file.path(counts_dir, "*", paste0("*", suffix))),
    Sys.glob(file.path(counts_dir, paste0("*", suffix)))
  ))
  if (length(files) == 0)
    return(data.frame(name = character(0), path = character(0), stringsAsFactors = FALSE))
  nm <- sub(paste0(suffix, "$"), "", basename(files))
  d  <- data.frame(name = nm, path = files, stringsAsFactors = FALSE)
  d[!duplicated(d$name), , drop = FALSE]      # one row per replicate name
}

# Scanning mode: --group1/--group2 are case-insensitive substring keywords matched
# against discovered replicate names. --group*-reps optionally assert the replicate
# COUNT and relabel the columns/figures (as <group>_<label>).
build_scan <- function() {
  g1 <- getc("group1"); g2 <- getc("group2")
  if (is.null(g1) || is.null(g2))
    stop("Provide either --samplesheet OR both --group1 and --group2 (scan keywords).")
  disc <- discover_counts(counts_dir, counts_suffix)
  if (nrow(disc) == 0)
    stop("no count files (*", counts_suffix, ") found under --counts-dir: ", counts_dir)
  pick <- function(kw) disc[grepl(tolower(kw), tolower(disc$name), fixed = TRUE), , drop = FALSE]
  m1 <- pick(g1); m2 <- pick(g2)
  if (nrow(m1) == 0) stop("no count files matched --group1 '", g1, "'. Discovered: ", paste(disc$name, collapse = ", "))
  if (nrow(m2) == 0) stop("no count files matched --group2 '", g2, "'. Discovered: ", paste(disc$name, collapse = ", "))
  ov <- intersect(m1$name, m2$name)
  if (length(ov) > 0)
    stop("these replicates matched BOTH groups (use more specific --group keywords): ", paste(ov, collapse = ", "))
  m1 <- m1[order(m1$name), , drop = FALSE]; m2 <- m2[order(m2$name), , drop = FALSE]
  lab <- function(m, reps, gg, flag) {
    if (!is.null(reps) && length(reps) != nrow(m))
      stop(sprintf("--%s lists %d label(s) but %d count file(s) matched '%s': %s",
                   flag, length(reps), nrow(m), gg, paste(m$name, collapse = ", ")))
    m$condition <- gg
    m$sample    <- if (is.null(reps)) m$name else make.unique(paste(gg, reps, sep = "_"))
    m
  }
  m1 <- lab(m1, split_csv(getc("group1-reps")), g1, "group1-reps")
  m2 <- lab(m2, split_csv(getc("group2-reps")), g2, "group2-reps")
  both <- rbind(m1, m2)
  data.frame(sample = both$sample, condition = both$condition,
             count_file = both$path, source_name = both$name, stringsAsFactors = FALSE)
}

build_coldata <- function() {
  if (!is.null(samplesheet)) {
    if (!file.exists(samplesheet)) stop("--samplesheet not found: ", samplesheet)
    sep <- if (grepl("\\.tsv$", samplesheet)) "\t" else ","
    df <- read.table(samplesheet, header = TRUE, sep = sep,
                     stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("sample","condition") %in% colnames(df)))
      stop("Sample sheet must have at least 'sample' and 'condition' columns.")
    df$sample    <- as.character(df$sample)
    df$condition <- as.character(df$condition)
    return(df)
  }
  build_scan()   # keyword-scanning mode
}

coldata <- tryCatch(build_coldata(), error = function(e) { message("ERROR: ", conditionMessage(e)); quit(status = 2) })

# Restrict to exactly the two conditions of interest, in a defined order -------
conds <- unique(coldata$condition)
if (length(conds) < 2) { message("ERROR: need >=2 conditions; found: ", paste(conds, collapse=", ")); quit(status = 2) }
if (length(conds) > 2 && is.null(design_override) && is.null(control)) {
  message("NOTE: >2 conditions present; --control is required to define the contrast baseline.")
}

# Reference/control level (denominator). LFC>0 => higher in the treatment level.
if (is.null(control)) {
  control <- if (!is.null(getc("group2"))) getc("group2") else sort(conds)[1]
  message("NOTE: --control not given; using '", control, "' as the reference (LFC denominator).")
}
if (!control %in% conds) { message("ERROR: --control '", control, "' is not one of the conditions."); quit(status = 2) }
treat <- setdiff(conds, control)
if (length(treat) != 1) treat <- treat[1]   # first non-control if >2 levels
if (is.null(contrast_name)) contrast_name <- paste0(treat, "_vs_", control)

## ===========================================================================
## 2. OUTPUT SCAFFOLDING + LOGGING
## ===========================================================================
# When --collapse-to-gene is set, the transcript-level analysis lives in a
# transcript_level/ subdir (gene_level/ is added later); otherwise everything
# goes straight into <outdir> exactly as before.
base_dir <- if (isTRUE(collapse_gene)) file.path(outdir, "transcript_level") else outdir
fig_dir <- file.path(base_dir, "figures")
tab_dir <- file.path(base_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

logfile <- file.path(outdir, "run.log")
log_con <- file(logfile, open = "wt")
logmsg <- function(...) {
  line <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(line, "\n"); cat(line, "\n", file = log_con, append = TRUE)
}
save_table <- function(df, name) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, optional = TRUE)
  # Flatten any list columns (e.g. multi-value annotation fields) so write.csv
  # never hits "unimplemented type 'list' in 'EncodeElement'".
  listcols <- names(df)[vapply(df, is.list, logical(1))]
  for (j in listcols)
    df[[j]] <- vapply(df[[j]], function(x) paste(as.character(x), collapse = ";"), character(1))
  if (length(listcols) > 0)
    logmsg("note: flattened list column(s) [", paste(listcols, collapse = ", "), "] in ", name)
  path <- file.path(tab_dir, name)
  write.csv(df, path, row.names = FALSE)
  logmsg("table  -> ", path)
}
# Save a base-graphics plot to BOTH pdf (vector) and png (quick view)
save_base <- function(name, fun, width = 7, height = 6) {
  tryCatch({
    pdf(file.path(fig_dir, paste0(name, ".pdf")), width = width, height = height); fun(); dev.off()
    png(file.path(fig_dir, paste0(name, ".png")), width = width*150, height = height*150, res = 150); fun(); dev.off()
    logmsg("figure -> ", name, ".{pdf,png}")
  }, error = function(e) {
    while (grDevices::dev.cur() > 1L) grDevices::dev.off()   # close any open device
    logmsg("WARNING: figure '", name, "' failed (", conditionMessage(e), "); continuing.")
  })
}
save_gg <- function(name, gg, width = 7, height = 6) {
  tryCatch({
    ggplot2::ggsave(file.path(fig_dir, paste0(name, ".pdf")), gg, width = width, height = height)
    ggplot2::ggsave(file.path(fig_dir, paste0(name, ".png")), gg, width = width, height = height, dpi = 150)
    logmsg("figure -> ", name, ".{pdf,png}")
  }, error = function(e) {
    while (grDevices::dev.cur() > 1L) grDevices::dev.off()
    logmsg("WARNING: figure '", name, "' failed (", conditionMessage(e), "); continuing.")
  })
}
has <- function(pkg) requireNamespace(pkg, quietly = TRUE)

logmsg("rnaseq_deseq2_dge.R starting")
logmsg("ref=", ref, "  contrast=", contrast_name, "  control(ref level)=", control)
logmsg("outdir=", outdir)

## ---------------------------------------------------------------------------
## 2b. RESOLVE COUNT FILES for every sample
## ---------------------------------------------------------------------------
resolve_count_file <- function(sample, path_col = NULL) {
  if (!is.null(path_col) && nzchar(path_col) && file.exists(path_col)) return(path_col)
  cand <- c(
    file.path(counts_dir, sample, paste0(sample, counts_suffix)),
    file.path(counts_dir, paste0(sample, counts_suffix))
  )
  hit <- cand[file.exists(cand)]
  if (length(hit) >= 1) return(hit[1])
  g <- Sys.glob(file.path(counts_dir, sample, paste0("*", counts_suffix)))
  if (length(g) >= 1) return(g[1])
  g2 <- Sys.glob(file.path(counts_dir, paste0(sample, "*", counts_suffix)))
  if (length(g2) >= 1) return(g2[1])
  NA_character_
}
# Scanning mode already populated count_file; only resolve for sample-sheet mode.
if (!"count_file" %in% colnames(coldata)) {
  path_col_vals <- if ("path" %in% colnames(coldata)) coldata$path else rep(NA_character_, nrow(coldata))
  coldata$count_file <- vapply(seq_len(nrow(coldata)),
                               function(i) resolve_count_file(coldata$sample[i], path_col_vals[i]),
                               character(1))
}
missing_files <- coldata$sample[is.na(coldata$count_file)]
if (length(missing_files) > 0) {
  message("ERROR: could not locate count files for: ", paste(missing_files, collapse = ", "))
  message("       looked under: ", counts_dir, " with suffix '", counts_suffix, "'")
  quit(status = 2)
}

# Design formula -------------------------------------------------------------
if (!is.null(design_override)) {
  design_formula <- stats::as.formula(design_override)
} else if (!is.null(batch_col)) {
  if (!batch_col %in% colnames(coldata)) { message("ERROR: --batch column '", batch_col, "' not in sample sheet."); quit(status = 2) }
  design_formula <- stats::as.formula(paste0("~ ", batch_col, " + condition"))
} else {
  design_formula <- ~ condition
}

# Report the resolved configuration (and stop here if --dry-run) --------------
cfg <- data.frame(
  key = c("ref","gtf","annotate","counts_dir","suffix","control(ref)","treatment",
          "contrast_name","design","alpha","lfc_threshold","alt_hypothesis","shrink",
          "fit_type","min_count","transform","top_n","biotype_filter","threads","outdir"),
  value = c(ref, gtf_path, do_annotate, counts_dir, counts_suffix, control, treat,
            contrast_name, paste(deparse(design_formula), collapse=""), alpha, lfc_thr, alt_hyp, shrink,
            fit_type, min_count, transform, top_n, ifelse(is.null(biotype_filt),"<none>",biotype_filt),
            threads, outdir),
  stringsAsFactors = FALSE)
save_table(cfg, "00_run_config.csv")
save_table(coldata, "01_sample_table.csv")
logmsg("resolved ", nrow(coldata), " replicates across conditions: ", paste(conds, collapse = ", "))
for (cc in conds) {
  mm <- coldata[coldata$condition == cc, , drop = FALSE]
  src <- if ("source_name" %in% colnames(mm)) mm$source_name else mm$sample
  logmsg("  group '", cc, "' (", nrow(mm), " reps): ", paste(src, collapse = ", "))
}

if (dry_run) {
  logmsg("--dry-run set: configuration validated, all count files found. Stopping before DESeq2.")
  cat("\nDRY RUN OK. Config written to ", tab_dir, ". Remove --dry-run to execute.\n", sep = "")
  quit(status = 0)
}

if (!ok_core) ok_core <- requireNamespace("DESeq2", quietly = TRUE)  # re-check (in case --rlib added it)
if (!ok_core) {
  message("ERROR: DESeq2 is not installed in this R environment (activate multiomics-r). ",
          "Config and sample validation still succeeded above.")
  quit(status = 3)
}
suppressWarnings(suppressMessages({
  library(DESeq2)
  if (has("ggplot2")) library(ggplot2)
}))
set.seed(seed)
BPPARAM <- if (threads > 1 && has("BiocParallel")) BiocParallel::MulticoreParam(threads) else NULL
if (!is.null(BPPARAM)) logmsg("BiocParallel: ", threads, " workers")

## ===========================================================================
## 3. READ htseq COUNTS  (+ assignment-rate QC from the __ summary rows)
## ===========================================================================
# htseq-count appends 5 summary rows starting with "__". They are NOT genes and
# must be removed before modelling, but they are excellent library QC.
read_htseq <- function(path) {
  df <- read.table(path, sep = "\t", stringsAsFactors = FALSE,
                   col.names = c("gene_id", "count"), quote = "", comment.char = "")
  is_special <- startsWith(df$gene_id, "__")
  list(genes = df[!is_special, , drop = FALSE],
       special = setNames(df$count[is_special], df$gene_id[is_special]))
}

logmsg("reading ", nrow(coldata), " htseq count files ...")
reads <- lapply(coldata$count_file, read_htseq)

# --- assignment QC table ----------------------------------------------------
special_names <- c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
qc <- data.frame(sample = coldata$sample, condition = coldata$condition, stringsAsFactors = FALSE)
qc$assigned <- vapply(reads, function(r) sum(r$genes$count), numeric(1))
for (sn in special_names)
  qc[[sub("^__","",sn)]] <- vapply(reads, function(r) { v <- r$special[[sn]]; if (is.null(v)) 0 else v }, numeric(1))
qc$total <- rowSums(qc[, c("assigned", sub("^__","",special_names))])
qc$assignment_rate <- round(qc$assigned / qc$total, 4)
save_table(qc, "02_htseq_assignment_qc.csv")

# --- build the integer count matrix (aligned on common gene_ids) ------------
gene_lists <- lapply(reads, function(r) r$genes$gene_id)
common_genes <- Reduce(intersect, gene_lists)
if (length(common_genes) == 0) { message("ERROR: no gene IDs shared across samples."); quit(status = 2) }
if (any(vapply(gene_lists, length, integer(1)) != length(common_genes)))
  logmsg("WARNING: samples had differing gene sets; using the ", length(common_genes),
         " gene IDs common to all.")
cts <- vapply(reads, function(r) {
  m <- match(common_genes, r$genes$gene_id); r$genes$count[m]
}, numeric(length(common_genes)))
rownames(cts) <- common_genes
colnames(cts) <- coldata$sample
if (strip_version) rownames(cts) <- sub("\\.\\d+$", "", rownames(cts))
storage.mode(cts) <- "integer"
logmsg("count matrix: ", nrow(cts), " genes x ", ncol(cts), " samples")

# --- library-level QC figures ----------------------------------------------
save_base("02_htseq_assignment_bars", function() {
  mat <- t(as.matrix(qc[, c("assigned", sub("^__","",special_names))]))
  colnames(mat) <- qc$sample
  # meaningful palette: assigned = green; the concerning multimappers red; rest muted
  asg_cols <- c(assigned = "#3B8C4E", no_feature = "#F2A65A", ambiguous = "#F2CC8F",
                too_low_aQual = "#BDBDBD", not_aligned = "#7F9DB9",
                alignment_not_unique = "#C1483B")[rownames(mat)]
  par(mar = c(9,5,3,9), xpd = TRUE)
  barplot(mat, col = asg_cols, border = NA, las = 2, ylab = "reads",
          main = "htseq read assignment per sample")
  legend("topright", inset = c(-0.30,0), legend = rownames(mat),
         fill = asg_cols, border = NA, bty = "n", cex = 0.75)
}, width = 9, height = 6)

save_base("03_library_size", function() {
  par(mar = c(9,5,3,2)); barplot(colSums(cts), las = 2, ylab = "assigned counts",
    main = "Library size (assigned counts)", col = as.integer(factor(coldata$condition)))
}, width = 8, height = 5)

detected <- colSums(cts > 0)
save_table(data.frame(sample = names(detected), genes_detected = detected), "03_genes_detected.csv")
save_base("03_genes_detected", function() {
  par(mar = c(9,5,3,2)); barplot(detected, las = 2, ylab = "# genes with count > 0",
    main = "Detected genes per sample", col = as.integer(factor(coldata$condition)))
}, width = 8, height = 5)

save_base("04_count_distribution_raw", function() {
  par(mar = c(9,5,3,2)); boxplot(log2(cts + 1), las = 2, ylab = "log2(count + 1)",
    main = "Raw count distribution", col = as.integer(factor(coldata$condition)))
}, width = 8, height = 5)

## ===========================================================================
## 4. ANNOTATION  (gene_id -> symbol / biotype / chrom)  [optional]
## ===========================================================================
annotate_genes <- function(gtf, strip_version) {
  if (!file.exists(gtf)) { logmsg("WARNING: GTF not found (", gtf, "); skipping annotation."); return(NULL) }
  ann <- NULL
  if (has("rtracklayer")) {
    ann <- tryCatch({
      gr <- rtracklayer::import(gtf)
      md <- as.data.frame(S4Vectors::mcols(gr))
      # Prefer 'gene' feature rows (GENCODE); if none exist (e.g. UCSC knownGene,
      # which only has exon/CDS), keep every row and collapse by gene_id below.
      if ("type" %in% colnames(md) && any(md$type == "gene")) {
        keep <- md$type == "gene"; gr <- gr[keep]; md <- md[keep, , drop = FALSE]
      }
      gname <- if ("gene_name" %in% colnames(md)) md$gene_name else md$gene_id
      btype <- if ("gene_type" %in% colnames(md)) md$gene_type else if ("gene_biotype" %in% colnames(md)) md$gene_biotype else NA
      data.frame(gene_id = as.character(md$gene_id),
                 gene_name = ifelse(is.na(gname), as.character(md$gene_id), as.character(gname)),
                 biotype = btype,
                 chrom = as.character(GenomeInfoDb::seqnames(gr)), stringsAsFactors = FALSE)
    }, error = function(e) { logmsg("rtracklayer import failed (", conditionMessage(e), "); trying manual parse."); NULL })
  }
  if (is.null(ann)) {   # base-R fallback parser (no Bioconductor needed)
    grab <- function(line, key) {
      m <- regmatches(line, regexpr(paste0(key, ' "[^"]+"'), line))
      if (length(m)) sub(paste0(key, ' "([^"]+)"'), "\\1", m) else NA_character_
    }
    # Stream the file; keep either 'gene' rows only, or (if feature=NULL) all rows.
    stream_parse <- function(feature) {
      con <- file(gtf, "r"); on.exit(close(con))
      ids <- character(0); nms <- character(0); bts <- character(0); chr <- character(0)
      repeat {
        chunk <- readLines(con, n = 100000L)
        if (length(chunk) == 0) break
        chunk <- chunk[!startsWith(chunk, "#")]
        if (length(chunk) == 0) next
        parts <- strsplit(chunk, "\t", fixed = TRUE)
        f1 <- vapply(parts, function(x) if (length(x) >= 1) x[1] else "", character(1))
        f3 <- vapply(parts, function(x) if (length(x) >= 3) x[3] else "", character(1))
        if (!is.null(feature)) { sel <- f3 == feature; chunk <- chunk[sel]; f1 <- f1[sel] }
        if (length(chunk) == 0) next
        ids <- c(ids, vapply(chunk, grab, character(1), key = "gene_id"))
        nms <- c(nms, vapply(chunk, grab, character(1), key = "gene_name"))
        bts <- c(bts, vapply(chunk, function(l) { b <- grab(l,"gene_type"); if (is.na(b)) grab(l,"gene_biotype") else b }, character(1)))
        chr <- c(chr, f1)
      }
      data.frame(gene_id = ids, gene_name = nms, biotype = bts, chrom = chr, stringsAsFactors = FALSE)
    }
    ann <- tryCatch({
      g <- stream_parse("gene")               # GENCODE path
      if (nrow(g) == 0) g <- stream_parse(NULL)  # UCSC knownGene: no 'gene' rows -> use all
      g$gene_name <- ifelse(is.na(g$gene_name), g$gene_id, g$gene_name)
      g
    }, error = function(e) { logmsg("manual GTF parse failed: ", conditionMessage(e)); NULL })
  }
  if (is.null(ann) || nrow(ann) == 0) return(NULL)
  if (strip_version) ann$gene_id <- sub("\\.\\d+$", "", ann$gene_id)
  ann[!duplicated(ann$gene_id), , drop = FALSE]
}

# --- external gene_id -> gene_name[/biotype] map (e.g. UCSC kgXref for hg19) --
read_id_map <- function(path, strip_version) {
  sep <- if (grepl("\\.tsv$", path)) "\t" else ","
  m <- read.table(path, header = TRUE, sep = sep, stringsAsFactors = FALSE, quote = "", check.names = FALSE)
  if (!"gene_id"   %in% colnames(m)) stop("--id-map must have a 'gene_id' column")
  if (!"gene_name" %in% colnames(m)) stop("--id-map must have a 'gene_name' column")
  if (strip_version) m$gene_id <- sub("\\.\\d+$", "", m$gene_id)
  m[!duplicated(m$gene_id), , drop = FALSE]
}

# Overlay an id-map onto GTF-derived annotation (id-map wins on gene_name/biotype)
apply_id_map <- function(ann, m) {
  if (is.null(ann)) {
    return(data.frame(gene_id = m$gene_id, gene_name = m$gene_name,
                      biotype = if ("biotype" %in% colnames(m)) m$biotype else NA,
                      chrom = NA_character_, stringsAsFactors = FALSE))
  }
  idx <- match(ann$gene_id, m$gene_id); hit <- !is.na(idx)
  ann$gene_name[hit] <- m$gene_name[idx[hit]]
  if ("biotype" %in% colnames(m)) ann$biotype[hit] <- m$biotype[idx[hit]]
  extra <- m[!m$gene_id %in% ann$gene_id, , drop = FALSE]
  if (nrow(extra) > 0)
    ann <- rbind(ann, data.frame(gene_id = extra$gene_id, gene_name = extra$gene_name,
                                 biotype = if ("biotype" %in% colnames(extra)) extra$biotype else NA,
                                 chrom = NA_character_, stringsAsFactors = FALSE))
  ann
}

# Optional: auto-map UCSC knownGene (uc########.#) IDs to symbols via Bioconductor
map_ucsc_to_symbol <- function(ann) {
  if (is.null(ann)) return(ann)
  uc <- grepl("^uc[0-9]", ann$gene_id)
  if (mean(uc) < 0.5) return(ann)   # not a knownGene-style ID set
  need <- c("TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db", "AnnotationDbi")
  if (!all(vapply(need, has, logical(1)))) {
    logmsg("hg19 knownGene IDs detected but ", paste(need[!vapply(need, has, logical(1))], collapse=", "),
           " not installed. Symbols will be UCSC accessions. Supply --id-map (UCSC kgXref) for symbols.")
    return(ann)
  }
  tryCatch({
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    ent  <- AnnotationDbi::mapIds(txdb, keys = ann$gene_id, column = "GENEID",
                                  keytype = "TXNAME", multiVals = "first")
    sym  <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = as.character(ent),
                                  column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
    ann$gene_name <- ifelse(is.na(sym), ann$gene_name, sym)
    logmsg("mapped UCSC knownGene -> symbols via TxDb + org.Hs.eg.db (",
           sum(!is.na(sym)), "/", nrow(ann), " symbols)")
    ann
  }, error = function(e) { logmsg("UCSC->symbol mapping failed: ", conditionMessage(e)); ann })
}

ann <- if (do_annotate) annotate_genes(gtf_path, strip_version) else NULL

# Enrich with symbols: --id-map takes precedence; otherwise try UCSC auto-map.
if (do_annotate) {
  if (!is.null(id_map)) {
    m <- tryCatch(read_id_map(id_map, strip_version),
                  error = function(e) { message("ERROR: ", conditionMessage(e)); quit(status = 2) })
    ann <- apply_id_map(ann, m)
    logmsg("applied --id-map: ", id_map, " (", nrow(m), " entries)")
  } else {
    ann <- map_ucsc_to_symbol(ann)
  }
}

# Coerce annotation columns to ATOMIC character. rtracklayer::import and
# AnnotationDbi::mapIds can return list / CharacterList columns; those break both
# write.csv AND ggplot facet/label sorting (e.g. the top-gene count facets and
# the top-DEG heatmap rownames). Fixing it here fixes every downstream use.
if (!is.null(ann)) {
  atomize <- function(x) {
    if (is.list(x)) vapply(x, function(v) if (length(v)) as.character(v)[[1]] else NA_character_, character(1))
    else as.character(x)
  }
  for (.c in intersect(c("gene_id","gene_name","biotype","chrom"), names(ann))) ann[[.c]] <- atomize(ann[[.c]])
  # Never leave gene_name empty/NA: fall back to the gene_id so no figure or
  # table ever shows "NA" (e.g. the few uc IDs with no UCSC symbol).
  if ("gene_name" %in% names(ann))
    ann$gene_name <- ifelse(is.na(ann$gene_name) | !nzchar(ann$gene_name), ann$gene_id, ann$gene_name)
}

if (!is.null(ann)) {
  n_match <- sum(rownames(cts) %in% ann$gene_id)
  matched_ann <- ann[ann$gene_id %in% rownames(cts), , drop = FALSE]
  # "symbols available" = gene_name differs from gene_id for a real fraction of genes
  sym_frac <- if (nrow(matched_ann) > 0)
    mean(!is.na(matched_ann$gene_name) & matched_ann$gene_name != matched_ann$gene_id) else 0
  rep_df <- data.frame(gtf = gtf_path, id_map = ifelse(is.null(id_map), "<none>", id_map),
                       n_genes_counts = nrow(cts), n_annotated = n_match,
                       frac_annotated = round(n_match / nrow(cts), 3),
                       symbols_available = sym_frac > 0.1,
                       frac_with_symbol = round(sym_frac, 3))
  save_table(rep_df, "05_annotation_report.csv")
  if ("biotype" %in% colnames(ann) && any(!is.na(ann$biotype))) {
    bt <- as.data.frame(table(ann$biotype[ann$gene_id %in% rownames(cts)]), stringsAsFactors = FALSE)
    if (ncol(bt) == 2 && nrow(bt) > 0) { colnames(bt) <- c("biotype","n_genes"); save_table(bt, "05_biotype_breakdown.csv") }
  }
  logmsg("annotation: ", n_match, "/", nrow(cts), " gene IDs matched; symbols on ",
         round(100*sym_frac), "% of matched genes")
  if (n_match / nrow(cts) < 0.5)
    logmsg("WARNING: <50% of gene IDs matched the annotation — check --ref / --strip-gene-version / --gtf.")
  if (sym_frac <= 0.1)
    logmsg("NOTE: no usable gene symbols (results will show gene IDs). For hg19 knownGene, ",
           "supply --id-map with a UCSC kgXref (kgID -> geneSymbol) table, or install ",
           "TxDb.Hsapiens.UCSC.hg19.knownGene + org.Hs.eg.db.")
} else if (do_annotate) {
  logmsg("annotation unavailable; results will use gene IDs only.")
}

# Optional biotype filter (before modelling) --------------------------------
if (!is.null(biotype_filt)) {
  if (is.null(ann)) { message("ERROR: --biotype-filter requires annotation."); quit(status = 2) }
  if (!("biotype" %in% colnames(ann)) || !any(!is.na(ann$biotype))) {
    message("ERROR: --biotype-filter requested but the annotation carries no biotype ",
            "(e.g. hg19 UCSC knownGene). Use grch38, or drop --biotype-filter."); quit(status = 2)
  }
  keep_bt <- rownames(cts) %in% ann$gene_id[!is.na(ann$biotype) & ann$biotype == biotype_filt]
  if (sum(keep_bt) == 0) { message("ERROR: --biotype-filter '", biotype_filt, "' matched 0 genes."); quit(status = 2) }
  logmsg("biotype filter '", biotype_filt, "': keeping ", sum(keep_bt), "/", nrow(cts), " genes")
  cts <- cts[keep_bt, , drop = FALSE]
}

## ===========================================================================
## 5-9. ANALYSIS ENGINE  (run once per level: transcript, and gene if collapsing)
##      cts/ann are parameters; save_table/save_base/save_gg write into level_dir.
## ===========================================================================
dge_and_figures <- function(cts, ann, level_dir, level_lab) {
fig_dir <<- file.path(level_dir, "figures")
tab_dir <<- file.path(level_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
logmsg("==== analysis level: ", level_lab, "  ->  ", level_dir, " ====")

## 5. BUILD DESeqDataSet, PRE-FILTER, RUN DESeq
coldata$condition <- factor(coldata$condition)
coldata$condition <- relevel(coldata$condition, ref = control)   # control = reference
if (!is.null(batch_col)) coldata[[batch_col]] <- factor(coldata[[batch_col]])
rownames(coldata) <- coldata$sample

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = design_formula)

# Pre-filtering (vignette: >= min_count in >= smallest group size) ------------
group_sizes <- table(coldata$condition)
smallest    <- min(group_sizes)
if (is.na(min_samples)) min_samples <- smallest
# sensitivity table: how many genes survive at a few thresholds
thr_grid <- expand.grid(min_count = c(1,5,10,20), min_samples = c(1, smallest))
thr_grid$genes_kept <- apply(thr_grid, 1, function(r)
  sum(rowSums(counts(dds) >= r[["min_count"]]) >= r[["min_samples"]]))
save_table(thr_grid, "06_prefilter_sensitivity.csv")

keep <- rowSums(counts(dds) >= min_count) >= min_samples
filt_summary <- data.frame(stage = c("input","kept","removed"),
                           n_genes = c(nrow(dds), sum(keep), sum(!keep)),
                           note = c("all genes in matrix",
                                    paste0(">=", min_count, " counts in >=", min_samples, " samples"),
                                    "low-count"))
save_table(filt_summary, "06_prefilter_summary.csv")
dds <- dds[keep, ]
logmsg("pre-filter: kept ", sum(keep), "/", length(keep), " genes")

logmsg("running DESeq() (fitType=", fit_type, ") ...")
dds <- if (is.null(BPPARAM)) {
  DESeq(dds, fitType = fit_type)
} else {
  DESeq(dds, fitType = fit_type, parallel = TRUE, BPPARAM = BPPARAM)
}

# Size factors + normalization diagnostics -----------------------------------
sf <- sizeFactors(dds)
save_table(data.frame(sample = names(sf), size_factor = as.numeric(sf),
                      library_size = colSums(counts(dds))), "07_size_factors.csv")
save_base("07_size_factor_vs_libsize", function() {
  plot(colSums(counts(dds)), sf, xlab = "library size (filtered)", ylab = "size factor",
       main = "Size factor vs library size", pch = 19, col = as.integer(coldata$condition))
  text(colSums(counts(dds)), sf, labels = colnames(dds), pos = 3, cex = 0.6)
})
norm_cts <- counts(dds, normalized = TRUE)
save_base("07_count_distribution_normalized", function() {
  par(mar = c(9,5,3,2)); boxplot(log2(norm_cts + 1), las = 2, ylab = "log2(norm count + 1)",
    main = "Normalized count distribution", col = as.integer(coldata$condition))
}, width = 8, height = 5)
# RLE (relative log expression) before & after normalization -----------------
rle_plot <- function(mat, main) {
  logm <- log2(mat + 1); med <- apply(logm, 1, median); rle <- logm - med
  par(mar = c(9,5,3,2)); boxplot(rle, las = 2, ylab = "RLE", main = main,
    col = as.integer(coldata$condition), outline = FALSE); abline(h = 0, lty = 2)
}
save_base("08_RLE_raw",        function() rle_plot(counts(dds), "RLE (raw counts)"), width = 8, height = 5)
save_base("08_RLE_normalized", function() rle_plot(norm_cts,    "RLE (normalized)"), width = 8, height = 5)

# Dispersion diagnostic ------------------------------------------------------
save_base("09_dispersion_estimates", function() plotDispEsts(dds, main = "Dispersion estimates"))

## ===========================================================================
## 6. TRANSFORMS + SAMPLE-LEVEL QC (PCA, distances, correlation, Cook's)
## ===========================================================================
do_vst  <- transform %in% c("vst","both")
do_rlog <- transform %in% c("rlog","both")
vsd <- if (do_vst)  tryCatch(vst(dds, blind = FALSE),  error = function(e) varianceStabilizingTransformation(dds, blind = FALSE)) else NULL
rld <- if (do_rlog) rlog(dds, blind = FALSE) else NULL
ntd <- normTransform(dds)                  # log2(n+1), for mean-SD baseline
qc_obj <- if (!is.null(vsd)) vsd else rld

# mean-SD plots (variance stabilization check) -------------------------------
if (has("vsn")) {
  save_base("10_meanSD_log2",  function() print(vsn::meanSdPlot(assay(ntd), plot = FALSE)$gg + ggplot2::ggtitle("mean-SD: log2(n+1)")))
  if (!is.null(vsd)) save_base("10_meanSD_vst", function() print(vsn::meanSdPlot(assay(vsd), plot = FALSE)$gg + ggplot2::ggtitle("mean-SD: VST")))
  if (!is.null(rld)) save_base("10_meanSD_rlog", function() print(vsn::meanSdPlot(assay(rld), plot = FALSE)$gg + ggplot2::ggtitle("mean-SD: rlog")))
}

# PCA ------------------------------------------------------------------------
intgroup <- if (!is.null(batch_col)) c("condition", batch_col) else "condition"
pca <- plotPCA(qc_obj, intgroup = intgroup, returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
save_table(pca, "11_pca_coordinates.csv")
if (has("ggplot2")) {
  aes_shape <- if (!is.null(batch_col)) ggplot2::aes(PC1, PC2, color = condition, shape = .data[[batch_col]]) else ggplot2::aes(PC1, PC2, color = condition)
  g <- ggplot2::ggplot(pca, aes_shape) + ggplot2::geom_point(size = 4) +
       ggplot2::xlab(paste0("PC1: ", pv[1], "% var")) + ggplot2::ylab(paste0("PC2: ", pv[2], "% var")) +
       ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.18)) +
       ggplot2::theme_bw() + ggplot2::ggtitle(paste0("PCA (", ifelse(do_vst,"VST","rlog"), ")"))
  # coord_fixed() dropped on purpose: with PC1 ~96% vs PC2 ~3% it squashes the
  # panel into an unreadable strip and clips labels; the axis titles report the
  # variance split, so the scale honesty is preserved.
  if (has("ggrepel")) {
    g <- g + ggrepel::geom_text_repel(ggplot2::aes(label = name), size = 3,
             show.legend = FALSE, max.overlaps = Inf, min.segment.length = 0)
  } else {
    g <- g + ggplot2::geom_text(ggplot2::aes(label = name), vjust = -1, size = 3, show.legend = FALSE)
  }
  save_gg("11_pca", g)
} else {
  save_base("11_pca", function() { plot(pca$PC1, pca$PC2, col = as.integer(factor(pca$condition)), pch = 19,
    xlab = paste0("PC1: ",pv[1],"%"), ylab = paste0("PC2: ",pv[2],"%"), main = "PCA"); text(pca$PC1, pca$PC2, pca$name, pos = 3, cex = 0.6) })
}
# Scree (variance explained by top PCs) --------------------------------------
rv <- matrixStats::rowVars(assay(qc_obj)); sel <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
pcaFull <- prcomp(t(assay(qc_obj)[sel, ]))
scree <- (pcaFull$sdev^2 / sum(pcaFull$sdev^2))[seq_len(min(10, length(pcaFull$sdev)))]
save_table(data.frame(PC = seq_along(scree), pct_variance = round(100*scree,2)), "11_pca_scree.csv")
save_base("11_pca_scree", function() barplot(100*scree, names.arg = paste0("PC",seq_along(scree)),
  ylab = "% variance", main = "PCA scree (top-500 variable genes)"))

# Sample-to-sample distance heatmap ------------------------------------------
sampleDists <- dist(t(assay(qc_obj)))
distMat <- as.matrix(sampleDists)
rownames(distMat) <- colnames(distMat) <- paste(colnames(qc_obj), coldata$condition, sep = " | ")
if (has("pheatmap")) {
  cols <- if (has("RColorBrewer")) colorRampPalette(rev(RColorBrewer::brewer.pal(9,"Blues")))(255) else NULL
  save_base("12_sample_distance_heatmap", function() pheatmap::pheatmap(distMat,
    clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists,
    color = cols, display_numbers = TRUE, number_format = "%.0f", number_color = "grey15",
    fontsize_number = 9, main = "Sample-to-sample Euclidean distance (VST)"))
} else {
  save_base("12_sample_distance_heatmap", function() heatmap(distMat, main = "Sample distances"))
}
# Sample correlation heatmap --------------------------------------------------
corMat <- cor(assay(qc_obj), method = "spearman")
save_table(as.data.frame(corMat), "12_sample_correlation.csv")
if (has("pheatmap")) save_base("12_sample_correlation_heatmap", function()
  pheatmap::pheatmap(corMat, display_numbers = TRUE, number_format = "%.2f",
    number_color = "grey15", fontsize_number = 9, main = "Sample Spearman correlation (VST)"))

# Cook's distance (outlier diagnostic) ---------------------------------------
save_base("13_cooks_distance", function() {
  par(mar = c(9,5,3,2)); boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2,
    main = "Cook's distance per sample", ylab = "log10(Cook's)", col = as.integer(coldata$condition))
}, width = 8, height = 5)

# Top-variable-gene heatmap (unsupervised) -----------------------------------
topVar <- head(order(matrixStats::rowVars(assay(qc_obj)), decreasing = TRUE), top_n)
if (has("pheatmap")) {
  anndf <- data.frame(condition = coldata$condition); rownames(anndf) <- colnames(qc_obj)
  mat <- assay(qc_obj)[topVar, ]; mat <- mat - rowMeans(mat)
  save_base("14_top_variable_genes_heatmap", function()
    pheatmap::pheatmap(mat, annotation_col = anndf, show_rownames = (top_n <= 60),
      main = paste0("Top ", top_n, " variable genes (VST, centered)")), width = 8, height = 9)
}

## ===========================================================================
## 7. DIFFERENTIAL TESTING + RESULTS
## ===========================================================================
logmsg("extracting results: ", treat, " vs ", control, " (alpha=", alpha, ", lfcThreshold=", lfc_thr, ")")
res <- if (lfc_thr > 0) {
  results(dds, contrast = c("condition", treat, control), alpha = alpha,
          lfcThreshold = lfc_thr, altHypothesis = alt_hyp)
} else {
  results(dds, contrast = c("condition", treat, control), alpha = alpha)
}

# Shrunken LFC (for ranking / MA / volcano) ----------------------------------
coef_name <- paste0("condition_", treat, "_vs_", control)
if (!coef_name %in% resultsNames(dds)) {          # labels with '-', spaces, etc.
  .cc <- grep("^condition_.*_vs_", resultsNames(dds), value = TRUE)
  if (length(.cc) == 1) coef_name <- .cc
}
resShrink <- NULL
shrink_used <- shrink
if (shrink != "none") {
  attempt <- function(type) {
    if (type == "apeglm") {
      if (!coef_name %in% resultsNames(dds)) stop("coef ", coef_name, " absent")
      lfcShrink(dds, coef = coef_name, type = "apeglm")
    } else if (type == "ashr") {
      lfcShrink(dds, contrast = c("condition", treat, control), type = "ashr")
    } else {
      lfcShrink(dds, contrast = c("condition", treat, control), type = "normal")
    }
  }
  resShrink <- tryCatch(attempt(shrink), error = function(e) {
    # 'normal' is built into DESeq2 (no extra package), so it is the safe fallback.
    if (shrink != "normal") {
      logmsg("NOTE: shrinkage '", shrink, "' unavailable (", conditionMessage(e),
             "); falling back to 'normal'. Install the '", shrink,
             "' package for the preferred estimator.")
      shrink_used <<- "normal"
      tryCatch(attempt("normal"), error = function(e2) {
        logmsg("WARNING: 'normal' shrinkage also failed: ", conditionMessage(e2))
        shrink_used <<- "none"; NULL })
    } else {
      logmsg("WARNING: 'normal' shrinkage failed: ", conditionMessage(e))
      shrink_used <<- "none"; NULL
    }
  })
}

# Assemble annotated results table -------------------------------------------
res_df <- as.data.frame(res); res_df$gene_id <- rownames(res_df)
if (!is.null(resShrink)) {
  sh <- as.data.frame(resShrink)[, "log2FoldChange", drop = FALSE]
  colnames(sh) <- "log2FoldChange_shrunk"; sh$gene_id <- rownames(sh)
  res_df <- merge(res_df, sh, by = "gene_id", all.x = TRUE)
}
if (!is.null(ann)) res_df <- merge(res_df, ann, by = "gene_id", all.x = TRUE)
res_df <- res_df[order(res_df$padj), ]
front <- intersect(c("gene_id","gene_name","biotype","chrom","baseMean",
                     "log2FoldChange","log2FoldChange_shrunk","lfcSE","stat","pvalue","padj"),
                   colnames(res_df))
res_df <- res_df[, c(front, setdiff(colnames(res_df), front))]
save_table(res_df, paste0("20_results_", contrast_name, "_full.csv"))

sig <- res_df[!is.na(res_df$padj) & res_df$padj < alpha & abs(res_df$log2FoldChange) >= lfc_thr, ]
save_table(sig, paste0("20_results_", contrast_name, "_significant.csv"))
save_table(as.data.frame(norm_cts) |> (\(d){ d$gene_id <- rownames(d); d })(),
           "21_normalized_counts.csv")
if (!is.null(vsd)) save_table(as.data.frame(assay(vsd)) |> (\(d){ d$gene_id <- rownames(d); d })(), "21_vst_matrix.csv")

# Summary + NA accounting ----------------------------------------------------
summ <- data.frame(
  metric = c("genes_tested","up (padj<alpha & LFC>0)","down (padj<alpha & LFC<0)",
             "NA_padj (low count / outlier)","alpha","lfc_threshold","shrink"),
  value  = c(sum(!is.na(res$pvalue)),
             sum(!is.na(res_df$padj) & res_df$padj < alpha & res_df$log2FoldChange > 0),
             sum(!is.na(res_df$padj) & res_df$padj < alpha & res_df$log2FoldChange < 0),
             sum(is.na(res_df$padj)), alpha, lfc_thr, shrink_used))
save_table(summ, paste0("20_results_", contrast_name, "_summary.csv"))
logmsg("DEGs @ padj<", alpha, ", |LFC|>=", lfc_thr, ": ", nrow(sig),
       " (up=", summ$value[2], ", down=", summ$value[3], ")")

# Threshold-sensitivity grid (robustness of DEG count) -----------------------
grid <- expand.grid(alpha = c(0.01,0.05,0.1), lfc = c(0,0.585,1,2))
grid$n_DEG <- apply(grid, 1, function(r)
  sum(!is.na(res_df$padj) & res_df$padj < r[["alpha"]] & abs(res_df$log2FoldChange) >= r[["lfc"]]))
save_table(grid, paste0("22_threshold_sensitivity_", contrast_name, ".csv"))

## ===========================================================================
## 8. RESULTS FIGURES
## ===========================================================================
# p-value histogram (calibration diagnostic)
save_base("30_pvalue_histogram", function()
  hist(res$pvalue[res$baseMean > metadata(res)$filterThreshold], breaks = 0:20/20,
       col = "grey80", main = "p-value histogram (passing filter)", xlab = "p-value"))

# Independent-filtering curve
if (!is.null(metadata(res)$filterNumRej))
  save_base("31_independent_filtering", function() {
    plot(metadata(res)$filterNumRej, type = "b", ylab = "# rejections",
         xlab = "quantiles of mean-normalized-count filter", main = "Independent filtering")
    lines(metadata(res)$lo.fit, col = "red"); abline(v = metadata(res)$filterTheta, lty = 2)
  })

# MA plots: raw MLE and shrunken
save_base("32_MA_raw", function() plotMA(res, ylim = c(-4,4), main = paste0("MA (MLE): ", contrast_name)))
if (!is.null(resShrink))
  save_base("32_MA_shrunk", function() plotMA(resShrink, ylim = c(-4,4), main = paste0("MA (", shrink_used, "-shrunk): ", contrast_name)))

# Volcano — colour by direction: Up (red) / Down (blue) / NS (grey)
volc <- res_df[!is.na(res_df$pvalue), ]
vthr <- max(lfc_thr, 0.585)                                  # visual >=1.5x guide lines
.sig <- !is.na(volc$padj) & volc$padj < alpha & abs(volc$log2FoldChange) >= max(lfc_thr, 0)
volc$regulation <- factor(ifelse(!.sig, "NS", ifelse(volc$log2FoldChange > 0, "Up", "Down")),
                          levels = c("Up","Down","NS"))
n_up <- sum(volc$regulation == "Up"); n_dn <- sum(volc$regulation == "Down")
lab_col <- if ("gene_name" %in% colnames(volc)) "gene_name" else "gene_id"
vcols   <- c(Up = "firebrick", Down = "royalblue3", NS = "grey75")
if (has("ggplot2")) {
  ord    <- volc[order(volc$padj), ]
  toplab <- rbind(head(ord[ord$regulation == "Up",   ], 10),
                  head(ord[ord$regulation == "Down", ], 10))
  g <- ggplot2::ggplot(volc, ggplot2::aes(log2FoldChange, -log10(pvalue), color = regulation)) +
    ggplot2::geom_point(alpha = 0.6, size = 1) +
    ggplot2::scale_color_manual(values = vcols, drop = FALSE,
      name   = paste0("padj<", alpha, if (lfc_thr > 0) paste0(" & |LFC|>=", lfc_thr) else ""),
      labels = c(paste0("Up (", n_up, ")"), paste0("Down (", n_dn, ")"), "Not sig.")) +
    ggplot2::geom_vline(xintercept = c(-1, 1) * vthr, lty = 2, col = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(alpha), lty = 3, col = "grey60") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = paste0("Volcano: ", contrast_name),
                  subtitle = paste0("Up ", n_up, "  |  Down ", n_dn,
                                    "   (padj<", alpha, ", |LFC|>=", max(lfc_thr, 0), ")"))
  if (has("ggrepel")) g <- g + ggrepel::geom_text_repel(data = toplab,
      ggplot2::aes(label = .data[[lab_col]]), size = 3, max.overlaps = 20, show.legend = FALSE)
  save_gg("33_volcano", g)
} else {
  save_base("33_volcano", function() {
    with(volc, plot(log2FoldChange, -log10(pvalue), pch = 19, cex = 0.4,
      col = vcols[as.character(regulation)], main = paste0("Volcano: ", contrast_name),
      xlab = "log2FoldChange", ylab = "-log10(pvalue)"))
    abline(v = c(-1, 1) * vthr, lty = 2); abline(h = -log10(alpha), lty = 3)
    legend("topright", legend = c(paste0("Up (", n_up, ")"), paste0("Down (", n_dn, ")"), "Not sig."),
           col = vcols, pch = 19, bty = "n")
  })
}

# Top-gene count plots
top_genes <- head(res_df$gene_id[!is.na(res_df$padj)], min(top_n, 12))
if (length(top_genes) > 0 && has("ggplot2")) {
  dlist <- lapply(top_genes, function(g) {
    d <- plotCounts(dds, gene = g, intgroup = "condition", returnData = TRUE)
    d$gene <- if (!is.null(ann)) { nm <- ann$gene_name[ann$gene_id == g]; if (length(nm)) as.character(nm)[1] else g } else g
    d })
  dd <- do.call(rbind, dlist)
  dd$gene <- factor(as.character(dd$gene), levels = unique(as.character(dd$gene)))
  g <- ggplot2::ggplot(dd, ggplot2::aes(condition, count, color = condition)) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.15)) +
    ggplot2::scale_y_log10() + ggplot2::facet_wrap(~ gene, scales = "free_y") +
    ggplot2::theme_bw() + ggplot2::ggtitle("Normalized counts: top genes") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_gg("34_top_gene_counts", g, width = 10, height = 8)
}

# Top-DEG heatmap (VST, z-scored) --------------------------------------------
sig_ids <- head(sig$gene_id, top_n)
if (length(sig_ids) >= 2 && has("pheatmap") && !is.null(qc_obj)) {
  mat <- assay(qc_obj)[sig_ids, , drop = FALSE]
  mat <- t(scale(t(mat)))                          # z-score per gene
  rn  <- if (!is.null(ann)) { m <- ann$gene_name[match(sig_ids, ann$gene_id)]; ifelse(is.na(m), sig_ids, m) } else sig_ids
  rownames(mat) <- make.unique(rn)
  anndf <- data.frame(condition = coldata$condition); rownames(anndf) <- colnames(qc_obj)
  save_base("35_top_DEG_heatmap", function()
    pheatmap::pheatmap(mat, annotation_col = anndf, show_rownames = (length(sig_ids) <= 60),
      main = paste0("Top ", length(sig_ids), " DEGs (VST z-score)")), width = 8, height = 10)
}

## ===========================================================================
## 9. PROVENANCE
## ===========================================================================
git_hash <- tryCatch(system2("git", c("-C", dirname(normalizePath(sub("--file=", "",
              grep("--file=", commandArgs(FALSE), value = TRUE)[1]))), "rev-parse","HEAD"),
              stdout = TRUE, stderr = FALSE), error = function(e) NA)
manifest <- c(
  "rnaseq_deseq2_dge.R — run manifest",
  paste0("timestamp:   ", format(Sys.time())),
  paste0("host:        ", Sys.info()[["nodename"]]),
  paste0("git_commit:  ", if (length(git_hash)) git_hash[1] else NA),
  paste0("R:           ", R.version.string),
  paste0("DESeq2:      ", as.character(packageVersion("DESeq2"))),
  paste0("level:       ", level_lab),
  paste0("ref:         ", ref),
  paste0("gtf:         ", gtf_path),
  paste0("contrast:    ", treat, " vs ", control, "  (", contrast_name, ")"),
  paste0("design:      ", paste(deparse(design_formula), collapse = "")),
  paste0("alpha:       ", alpha, "   lfcThreshold: ", lfc_thr, "   shrink: ", shrink_used),
  paste0("n_samples:   ", ncol(dds), "   n_genes_tested: ", sum(!is.na(res$pvalue))),
  paste0("n_DEG:       ", nrow(sig)),
  paste0("args:        ", paste(raw, collapse = " ")),
  "outputs:", paste0("  tables:  ", tab_dir), paste0("  figures: ", fig_dir))
writeLines(manifest, file.path(level_dir, "MANIFEST.txt"))
writeLines(capture.output(sessionInfo()), file.path(level_dir, "sessionInfo.txt"))
logmsg("DONE (", level_lab, "). DEGs=", nrow(sig), ". See ", level_dir)
}  # ---- end dge_and_figures ----

## ===========================================================================
## DISPATCH: transcript level always; gene level if --collapse-to-gene
## ===========================================================================
dge_and_figures(cts, ann, base_dir, "transcript")

if (isTRUE(collapse_gene)) {
  if (is.null(ann)) {
    logmsg("WARNING: --collapse-to-gene needs annotation (gene symbols); skipping gene level.")
  } else {
    # Sum transcript-model counts per gene symbol. gene_name already falls back
    # to the transcript id when no symbol exists, so every row maps to a key.
    # NOTE: this consolidates isoform-split signal but cannot recover reads the
    # counter already discarded as ambiguous -- true gene-level counting needs a
    # gene-level GTF (gencode) or union/featureCounts upstream.
    key <- ann$gene_name[match(rownames(cts), ann$gene_id)]
    key[is.na(key) | !nzchar(key)] <- rownames(cts)[is.na(key) | !nzchar(key)]
    cts_gene <- rowsum(cts, group = key); storage.mode(cts_gene) <- "integer"
    ann_gene <- data.frame(gene_id = rownames(cts_gene), gene_name = rownames(cts_gene),
                           biotype = NA_character_,
                           chrom = ann$chrom[match(rownames(cts_gene), ann$gene_name)],
                           stringsAsFactors = FALSE)
    dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
    write.csv(data.frame(level = c("transcript", "gene"),
                         n_features = c(nrow(cts), nrow(cts_gene))),
              file.path(outdir, "tables", "gene_collapse_report.csv"), row.names = FALSE)
    logmsg("gene-level collapse: ", nrow(cts), " transcript models -> ", nrow(cts_gene), " genes")
    dge_and_figures(cts_gene, ann_gene, file.path(outdir, "gene_level"), "gene")
  }
}

close(log_con)
cat("\n==> Finished.",
    if (isTRUE(collapse_gene)) paste0(" See ", outdir, "/{transcript_level,gene_level}/")
    else paste0(" See ", outdir), "\n", sep = "")