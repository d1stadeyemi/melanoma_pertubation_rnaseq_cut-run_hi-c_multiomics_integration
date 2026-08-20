#!/usr/bin/env Rscript
###############################################################################
# deseq2_01_setup.R   --  STAGE 1 of 4:  DEFINE AND VALIDATE THE EXPERIMENT
#
# Nothing statistical happens here. This stage answers four questions and stops
# if it cannot answer any of them:
#     WHICH samples?      WHICH groups?      WHICH model?      Is the model FITTABLE?
#
# It reads no count data, so it runs in seconds. Run it alone with --dry-run
# before committing a cluster job: almost every failure in a DESeq2 analysis is
# a design failure, and this stage finds them all before the expensive steps.
#
# WRITES
#   tables/00_run_config.{csv,tsv}            every resolved parameter
#   tables/01_sample_table.{csv,tsv}          the samples that will be analysed
#   tables/01b_covariate_by_condition.{csv,tsv}   batch x condition balance
#   tables/01c_samplesheet_template.csv       FILL-IN TEMPLATE (see below)
#   _state/01_state.rds                       handoff to stage 2
#
#  ==========================================================================
#  PARAMETERS OWNED BY THIS STAGE
#  ==========================================================================
#  --ref <hg19|grch38>        Default: grch38. Chooses the default annotation
#                             GTF (stage 2 uses it; recorded here).
#  --counts-dir <DIR>         Output directory of rnaseq_fastq_to_counts.sh.
#  --counts-suffix <STR>      Default: _genecounts.tsv
#
#  HOW SAMPLES ARE ASSIGNED TO GROUPS -- pick ONE:
#   A) --samplesheet <CSV>    columns sample,condition[,batch,path,...]
#                             This is what submit_rnaseq_fastq_to_counts.sh
#                             auto-generates. PREFERRED.
#   B) --group1 <KEY> --group2 <KEY>
#                             case-insensitive substring match against the
#                             replicate names discovered in --counts-dir.
#
#  --control <NAME>           DENOMINATOR of the log2 fold change.
#                             LFC > 0 means UP in the treatment group.
#  --treatment <NAME>         NUMERATOR. Required when >2 conditions exist.
#  --contrast-name <STR>      Label for filenames. Default <treat>_vs_<control>.
#
#  BATCH / COVARIATES -- three routes, in increasing order of effort:
#  --batch-map "rep1=A,rep2=B"
#                             PATTERN=LABEL pairs matched against sample names.
#                             Covers the common case ("rep1s were prepped in
#                             June, rep2s in August") with no file at all.
#                             Creates the column named by --batch.
#  --covariates <CSV>         a `sample` column plus any extra columns, joined
#                             onto the sample table. Use 01c_samplesheet_
#                             template.csv as the starting point -- it is
#                             written with the sample names already filled in.
#  --samplesheet <CSV>        put `batch` straight in the sheet.
#
#  --batch <COL>              Add COL as a covariate: design becomes
#                             ~ <COL> + condition. Coerced to a FACTOR.
#                             Term order is irrelevant because the contrast is
#                             always extracted explicitly by name.
#  --design <FORMULA>         Full override, e.g. "~ batch + condition".
#                             Must contain a term named `condition`.
#                             INTERACTIONS (* or :) are REFUSED -- see below.
#  --numeric-covars <a,b>     Design terms to keep NUMERIC (e.g. RIN). Anything
#                             not listed is a factor. Note: with a continuous
#                             variable in the design, DESeq2 does NOT perform
#                             automatic outlier detection and replacement.
#  --collapse-technical <COL> Sum counts across TECHNICAL replicates (multiple
#                             sequencing runs of the SAME library). Never use
#                             on biological replicates. Applied in stage 3.
#  --subset-to-contrast       With >2 conditions, drop the other groups before
#                             fitting. See "WHEN TO SUBSET" below.
#
#  --outdir <DIR>             Default ./deseq2_out
#  --dry-run                  Stop after this stage.
#  --force                    Reuse stale upstream state (stages 2-4).
#
#  ==========================================================================
#  WHY INTERACTION DESIGNS ARE REFUSED HERE
#  ==========================================================================
#  `--design "~ condition*timepoint"` passes a naive validation: all.vars()
#  finds `condition` and the rank check succeeds. But then
#  contrast=c("condition", treat, control) returns the condition effect AT THE
#  REFERENCE LEVEL of timepoint only -- not an average effect -- and the
#  vignette states the potential for over-shrinking LFC "is also why DESeq2's
#  shrinkage estimator is not recommended for designs with interaction terms".
#  Silently reporting a reference-level-only effect as "the" condition effect is
#  the kind of error that survives to publication, so this stage stops instead.
#  The vignette's own advice for the common case: rather than an interaction
#  term, combine the factors into a single grouping variable and contrast its
#  levels. Do that in the sample sheet's `condition` column.
#
#  ==========================================================================
#  WHEN TO SUBSET (--subset-to-contrast)
#  ==========================================================================
#  Default is OFF: all groups are fitted together, which shares information
#  across samples and is the better estimator in the usual case. Turn it ON only
#  when the PCA or the sample-distance heatmap from stage 3 shows one group with
#  much higher within-group spread than the others. The vignette describes
#  exactly this: comparing two groups after subsetting the DESeqDataSet to only
#  those groups "will be more sensitive than a model including all samples
#  together", while noting such an extreme range of within-group variability is
#  not common. So: fit everything first, look at the EDA, subset only if the
#  picture demands it -- and say so in the methods.
###############################################################################

.args_raw <- commandArgs(trailingOnly = TRUE)
.self <- (function() {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) normalizePath(sub("^--file=", "", a[1]), mustWork = FALSE) else NA_character_
})()
.common <- file.path(dirname(.self), "deseq2_common.R")
if (!file.exists(.common)) { message("ERROR: deseq2_common.R not found next to ", .self); quit(status = 2) }
source(.common)

help_text <- function() {
'
rnaseq DESeq2 pipeline (v3, 4 stages) -- stage 1: setup

  deseq2_01_setup.R    define + validate the experiment   <-- run this first
  deseq2_02_counts.R   load counts, annotate, library QC
  deseq2_03_model.R    prefilter, DESeq() fit, transforms, sample QC
  deseq2_04_results.R  contrasts, shrinkage, DEG tables, figures

Normally you run all four at once:
  sbatch submit_rnaseq_deseq2_dge.sh --ref grch38 --counts-dir DIR \\
         --samplesheet samples.csv --control Parental --treatment Bottom10 \\
         --batch batch --outdir RESULTS

Each stage documents its own parameters in its header. Summary:

STAGE 1  samples + design
  --ref --counts-dir --counts-suffix --samplesheet --covariates
  --group1 --group2 --group1-reps --group2-reps
  --control --treatment --contrast-name
  --batch --batch-map --design --numeric-covars --collapse-technical
  --subset-to-contrast
STAGE 2  counts + annotation
  --gtf --id-map --no-annotation --strip-gene-version --biotype-filter
  --collapse-to-gene
STAGE 3  model
  --min-count (10) --min-samples (smallest group) --test <wald|lrt> --reduced
  --fit-type <parametric|local|mean|glmGamPoi> --transform <vst|rlog|both>
STAGE 4  results
  --alpha (0.05) --lfc-report (0.585 = 1.5-fold) --lfc-threshold (0)
  --alt-hypothesis <greaterAbs|lessAbs|greater|less> --rank-on <mle|shrunk>
  --shrink <apeglm|ashr|normal|none> --no-independent-filtering --top-n (30)
GLOBAL
  --outdir --threads --seed --rlib --stage <1|2|3|4|all> --dry-run --force -h

THRESHOLDS: --lfc-report is the DEG-CALLING cutoff (applied with --alpha).
--lfc-threshold is the TESTING null handed to results(lfcThreshold=); it
changes what the p-value MEANS. They are separate on purpose. See stage 4.
'
}

opt <- tryCatch(parse_args(.args_raw),
                error = function(e) { message("ERROR: ", conditionMessage(e)); quit(status = 2) })
if (isTRUE(opt[["help"]]) || isTRUE(opt[["h"]]) || length(.args_raw) == 0) {
  cat(help_text()); quit(status = 0)
}
.unknown <- setdiff(names(opt), .known_keys)
if (length(.unknown) > 0) {
  message("ERROR: unknown option(s): ", paste0("--", .unknown, collapse = ", "),
          "\n       run with --help to see valid flags.")
  quit(status = 2)
}
g <- mk_getters(opt)

outdir <- g$chr("outdir", "deseq2_out")
log_open(outdir, "1/setup", truncate = TRUE)   # stage 1 starts a fresh log
logbanner(paste0("rnaseq DESeq2 pipeline ", PIPELINE_VERSION, "  --  STAGE 1: SETUP"))
logmsg("args: ", paste(.args_raw, collapse = " "))

## ---------------------------------------------------------------------------
## 1.1  Resolve and validate every stage-1 parameter
## ---------------------------------------------------------------------------
ref <- tolower(g$chr("ref", "grch38"))
if (!ref %in% c("hg19", "grch38")) die("--ref must be 'hg19' or 'grch38' (got '", ref, "').")

counts_dir    <- g$chr("counts-dir")
counts_suffix <- g$chr("counts-suffix", "_genecounts.tsv")
if (is.null(counts_dir))     die("--counts-dir is required.")
if (!dir.exists(counts_dir)) die("--counts-dir not found: ", counts_dir)
counts_dir <- normalizePath(counts_dir)

samplesheet  <- g$chr("samplesheet")
covariates_f <- g$chr("covariates")
batch_col    <- g$chr("batch")
batch_map    <- g$chr("batch-map")
design_over  <- g$chr("design")
num_covars   <- split_csv(g$chr("numeric-covars"))
tech_col     <- g$chr("collapse-technical")
control      <- g$chr("control")
treatment    <- g$chr("treatment")
contrast_nm  <- g$chr("contrast-name")
subset_ctr   <- g$bool("subset-to-contrast")
dry_run      <- g$bool("dry-run")
threads      <- g$int("threads", { s <- Sys.getenv("SLURM_CPUS_PER_TASK")
                                   if (nzchar(s) && grepl("^[0-9]+$", s)) as.integer(s) else 1L },
                      min = 1, max = 512)
seed         <- g$int("seed", 1L)

# Downstream stages own these, but check the SYNTAX here so a typo fails in
# seconds rather than after the count matrix has been built.
for (.k in c("min-count","min-samples","top-n"))
  if (!is.null(opt[[.k]]) && !grepl("^-?[0-9]+$", trimws(opt[[.k]])))
    die("--", .k, " must be a whole number (got '", opt[[.k]], "').")
for (.k in c("alpha","lfc-report","lfc-threshold"))
  if (!is.null(opt[[.k]]) && is.na(suppressWarnings(as.numeric(opt[[.k]]))))
    die("--", .k, " must be a number (got '", opt[[.k]], "').")

if (!is.null(batch_map) && is.null(batch_col)) {
  batch_col <- "batch"
  logmsg("NOTE: --batch-map given without --batch; using column name 'batch'.")
}

## ---------------------------------------------------------------------------
## 1.2  Build the sample table
## ---------------------------------------------------------------------------
discover_counts <- function(dir, suffix) {
  files <- unique(c(Sys.glob(file.path(dir, "*", paste0("*", suffix))),
                    Sys.glob(file.path(dir, paste0("*", suffix)))))
  if (length(files) == 0)
    return(data.frame(name = character(0), path = character(0), stringsAsFactors = FALSE))
  nm <- sub(paste0(suffix, "$"), "", basename(files))
  d <- data.frame(name = nm, path = files, stringsAsFactors = FALSE)
  d[!duplicated(d$name), , drop = FALSE]
}

build_scan <- function() {
  g1 <- g$chr("group1"); g2 <- g$chr("group2")
  if (is.null(g1) || is.null(g2))
    die("provide either --samplesheet OR both --group1 and --group2 (scan keywords).")
  disc <- discover_counts(counts_dir, counts_suffix)
  if (nrow(disc) == 0)
    die("no count files (*", counts_suffix, ") found under --counts-dir: ", counts_dir)
  pick <- function(kw) disc[grepl(tolower(kw), tolower(disc$name), fixed = TRUE), , drop = FALSE]
  m1 <- pick(g1); m2 <- pick(g2)
  if (nrow(m1) == 0) die("no count files matched --group1 '", g1, "'. Discovered: ",
                         paste(disc$name, collapse = ", "))
  if (nrow(m2) == 0) die("no count files matched --group2 '", g2, "'. Discovered: ",
                         paste(disc$name, collapse = ", "))
  ov <- intersect(m1$name, m2$name)
  if (length(ov) > 0)
    die("these replicates matched BOTH groups (use more specific keywords): ",
        paste(ov, collapse = ", "))
  m1 <- m1[order(m1$name), , drop = FALSE]; m2 <- m2[order(m2$name), , drop = FALSE]
  lab <- function(m, reps, gg, flag) {
    if (!is.null(reps) && length(reps) != nrow(m))
      die(sprintf("--%s lists %d label(s) but %d count file(s) matched '%s': %s",
                  flag, length(reps), nrow(m), gg, paste(m$name, collapse = ", ")))
    m$condition <- gg
    m$sample <- if (is.null(reps)) m$name else make.unique(paste(gg, reps, sep = "_"))
    m
  }
  m1 <- lab(m1, split_csv(g$chr("group1-reps")), g1, "group1-reps")
  m2 <- lab(m2, split_csv(g$chr("group2-reps")), g2, "group2-reps")
  both <- rbind(m1, m2)
  data.frame(sample = both$sample, condition = both$condition,
             count_file = both$path, source_name = both$name, stringsAsFactors = FALSE)
}

if (!is.null(samplesheet)) {
  if (!file.exists(samplesheet)) die("--samplesheet not found: ", samplesheet)
  coldata <- read_delim_auto(samplesheet)
  if (!all(c("sample","condition") %in% colnames(coldata)))
    die("sample sheet must have at least 'sample' and 'condition' columns. Found: ",
        paste(colnames(coldata), collapse = ", "))
  coldata$sample    <- trimws(as.character(coldata$sample))
  coldata$condition <- trimws(as.character(coldata$condition))
  if (any(!nzchar(coldata$sample)))    die("sample sheet has an empty 'sample' value.")
  if (any(!nzchar(coldata$condition))) die("sample sheet has an empty 'condition' value.")
  # An auto-generated sheet may still carry failed rows.
  if ("status" %in% colnames(coldata)) {
    bad <- coldata$status %in% c("MISSING_COUNTS", "FAILED")
    if (any(bad))
      die(sum(bad), " row(s) in the sample sheet have status MISSING_COUNTS/FAILED: ",
          paste(coldata$sample[bad], collapse = ", "),
          "\n       Re-run those samples, or delete those rows from the sheet.")
  }
  logmsg("sample sheet: ", samplesheet, " (", nrow(coldata), " rows)")
} else {
  coldata <- build_scan()
  logmsg("keyword scan of ", counts_dir, " -> ", nrow(coldata), " replicates")
}

dup_s <- unique(coldata$sample[duplicated(coldata$sample)])
if (length(dup_s) > 0)
  die("duplicate sample name(s): ", paste(dup_s, collapse = ", "),
      "\n       Sample names must be unique -- they become count-matrix column names.")

## ---------------------------------------------------------------------------
## 1.3  Attach covariates:  --batch-map, then --covariates
## ---------------------------------------------------------------------------
if (!is.null(batch_map)) {
  pairs <- split_csv(batch_map)
  if (!all(grepl("=", pairs, fixed = TRUE)))
    die("--batch-map must be PATTERN=LABEL pairs, e.g. \"rep1=A,rep2=B\" (got '", batch_map, "').")
  pat <- sub("=.*$", "", pairs); lab <- sub("^[^=]*=", "", pairs)
  if (any(!nzchar(pat)) || any(!nzchar(lab)))
    die("--batch-map has an empty pattern or label: '", batch_map, "'.")
  assigned <- rep(NA_character_, nrow(coldata))
  hits <- integer(length(pat))
  for (i in seq_along(pat)) {
    m <- grepl(pat[i], coldata$sample, fixed = TRUE)
    hits[i] <- sum(m)
    if (any(m & !is.na(assigned)))
      die("--batch-map pattern '", pat[i], "' overlaps an earlier pattern; ",
          "each sample must match exactly one pattern.")
    assigned[m] <- lab[i]
  }
  if (any(hits == 0))
    die("--batch-map pattern(s) matched no sample: ", paste(pat[hits == 0], collapse = ", "),
        "\n       Sample names are: ", paste(coldata$sample, collapse = ", "))
  if (any(is.na(assigned)))
    die("--batch-map left sample(s) unassigned: ", paste(coldata$sample[is.na(assigned)], collapse = ", "),
        "\n       Every sample needs a batch label. Add a pattern, or use --covariates.")
  coldata[[batch_col]] <- assigned
  logmsg("--batch-map -> column '", batch_col, "': ",
         paste(sprintf("%s=%s(n=%d)", pat, lab, hits), collapse = "  "))
}

if (!is.null(covariates_f)) {
  if (!file.exists(covariates_f)) die("--covariates file not found: ", covariates_f)
  cv <- read_delim_auto(covariates_f)
  if (!"sample" %in% colnames(cv)) die("--covariates file must have a 'sample' column.")
  cv$sample <- trimws(as.character(cv$sample))
  if (anyDuplicated(cv$sample))
    die("--covariates has duplicate 'sample' rows: ",
        paste(unique(cv$sample[duplicated(cv$sample)]), collapse = ", "))
  miss <- setdiff(coldata$sample, cv$sample)
  if (length(miss) > 0)
    die("--covariates is missing row(s) for: ", paste(miss, collapse = ", "),
        "\n       Sample names must match exactly. Start from the template this\n",
        "       stage writes: ", file.path(outdir, "tables", "01c_samplesheet_template.csv"))
  add <- setdiff(colnames(cv), c("sample", colnames(coldata)))
  idx <- match(coldata$sample, cv$sample)
  for (cn in add) coldata[[cn]] <- cv[[cn]][idx]
  logmsg("--covariates added column(s): ", paste(add, collapse = ", "))
}

## ---------------------------------------------------------------------------
## 1.4  Resolve the contrast
## ---------------------------------------------------------------------------
raw_conds <- unique(coldata$condition)
clean <- sanitize_level(coldata$condition)
if (!identical(clean, coldata$condition)) {
  chg <- unique(data.frame(from = coldata$condition, to = clean,
                           stringsAsFactors = FALSE))
  chg <- chg[chg$from != chg$to, , drop = FALSE]
  for (i in seq_len(nrow(chg)))
    logmsg("NOTE: condition '", chg$from[i], "' -> '", chg$to[i],
           "' (DESeq2 factor levels allow only letters, numbers, _ and .)")
  coldata$condition <- clean
  if (!is.null(control))   control   <- sanitize_level(control)
  if (!is.null(treatment)) treatment <- sanitize_level(treatment)
}
conds <- unique(coldata$condition)
if (length(conds) < 2) die("need >= 2 conditions; found: ", paste(conds, collapse = ", "))

if (is.null(control)) {
  control <- if (!is.null(g$chr("group2"))) sanitize_level(g$chr("group2")) else sort(conds)[1]
  logmsg("NOTE: --control not given; using '", control, "' as the reference (LFC denominator).")
}
if (!control %in% conds)
  die("--control '", control, "' is not one of the conditions: ", paste(conds, collapse = ", "))

if (is.null(treatment)) {
  others <- setdiff(conds, control)
  if (length(others) != 1)
    die("there are ", length(conds), " conditions (", paste(conds, collapse = ", "), "). ",
        "--treatment is required to say which one is the fold-change numerator.")
  treatment <- others
} else {
  if (!treatment %in% conds)
    die("--treatment '", treatment, "' is not one of the conditions: ", paste(conds, collapse = ", "))
  if (identical(treatment, control)) die("--treatment and --control are the same ('", control, "').")
}
treat <- treatment
if (is.null(contrast_nm)) contrast_nm <- paste0(treat, "_vs_", control)
contrast_nm <- gsub("[^A-Za-z0-9._-]+", "_", contrast_nm)

if (subset_ctr && length(conds) > 2) {
  keep <- coldata$condition %in% c(control, treat)
  logmsg("--subset-to-contrast: dropping ", sum(!keep), " sample(s) from ",
         paste(setdiff(conds, c(control, treat)), collapse = ", "))
  coldata <- coldata[keep, , drop = FALSE]
  conds <- unique(coldata$condition)
} else if (subset_ctr) {
  logmsg("NOTE: --subset-to-contrast has no effect (only 2 conditions present).")
}

n_ctrl <- sum(coldata$condition == control); n_treat <- sum(coldata$condition == treat)
if (n_ctrl < 2 || n_treat < 2)
  die("each side of the contrast needs >= 2 replicates to estimate dispersion (",
      control, "=", n_ctrl, ", ", treat, "=", n_treat, ").")
if (n_ctrl == 2 || n_treat == 2)
  logmsg("NOTE: only 2 replicates in a group (", control, "=", n_ctrl, ", ", treat, "=", n_treat,
         "). Dispersion estimates will be imprecise; treat the DEG list as exploratory.")

## ---------------------------------------------------------------------------
## 1.5  Resolve count file paths
## ---------------------------------------------------------------------------
resolve_count_file <- function(sample, path_col = NULL) {
  if (!is.null(path_col) && !is.na(path_col) && nzchar(path_col) && file.exists(path_col))
    return(normalizePath(path_col))
  cand <- c(file.path(counts_dir, sample, paste0(sample, counts_suffix)),
            file.path(counts_dir, paste0(sample, counts_suffix)))
  hit <- cand[file.exists(cand)]
  if (length(hit) >= 1) return(normalizePath(hit[1]))
  gg <- Sys.glob(file.path(counts_dir, sample, paste0("*", counts_suffix)))
  if (length(gg) >= 1) return(normalizePath(gg[1]))
  g2 <- Sys.glob(file.path(counts_dir, paste0(sample, "*", counts_suffix)))
  if (length(g2) >= 1) return(normalizePath(g2[1]))
  NA_character_
}
if (!"count_file" %in% colnames(coldata)) {
  pc <- if ("path" %in% colnames(coldata)) as.character(coldata$path) else rep(NA_character_, nrow(coldata))
  coldata$count_file <- vapply(seq_len(nrow(coldata)),
                               function(i) resolve_count_file(coldata$sample[i], pc[i]), character(1))
} else {
  coldata$count_file <- normalizePath(coldata$count_file, mustWork = FALSE)
}
missing_files <- coldata$sample[is.na(coldata$count_file) | !file.exists(coldata$count_file)]
if (length(missing_files) > 0) {
  logmsg("looked under: ", counts_dir, " with suffix '", counts_suffix, "'")
  die("could not locate count files for: ", paste(missing_files, collapse = ", "))
}

## ---------------------------------------------------------------------------
## 1.6  Build the design and TYPE every variable explicitly
## ---------------------------------------------------------------------------
if (!is.null(design_over)) {
  # Interaction guard -- see the header block for the full reasoning.
  if (grepl("[*:]", design_over))
    die("--design contains an interaction term (* or :): ", design_over, "\n",
        "       This pipeline extracts contrast=c(\"condition\", treatment, control),\n",
        "       which under an interaction returns the condition effect at the\n",
        "       REFERENCE level of the other factor only -- not an average effect.\n",
        "       DESeq2's shrinkage is also not recommended for interaction designs.\n",
        "       Instead, combine the factors into ONE grouping variable and put it\n",
        "       in the sample sheet's `condition` column (the vignette's own advice),\n",
        "       e.g. condition = Parental_0h, Parental_24h, B10_0h, B10_24h.")
  design_formula <- tryCatch(stats::as.formula(design_over),
                             error = function(e) die("--design is not a valid formula: ", design_over))
} else if (!is.null(batch_col)) {
  if (!batch_col %in% colnames(coldata))
    die("--batch column '", batch_col, "' is not in the sample table.\n",
        "       Available columns: ", paste(colnames(coldata), collapse = ", "), "\n",
        "       Supply it via --samplesheet, --covariates, or --batch-map.")
  if (make.names(batch_col) != batch_col)
    die("--batch column name '", batch_col, "' is not a valid R name; rename it.")
  design_formula <- stats::as.formula(paste0("~ ", batch_col, " + condition"))
} else design_formula <- ~ condition

design_vars <- all.vars(design_formula)
if (!"condition" %in% design_vars)
  die("the design must contain a term named `condition` (got: ",
      paste(deparse(design_formula), collapse = ""), ").")
miss_vars <- setdiff(design_vars, colnames(coldata))
if (length(miss_vars) > 0)
  die("design term(s) not in the sample table: ", paste(miss_vars, collapse = ", "),
      "\n       Available columns: ", paste(colnames(coldata), collapse = ", "))

coldata$condition <- factor(coldata$condition, levels = c(control, setdiff(conds, control)))
covar_vars  <- setdiff(design_vars, "condition")
covar_types <- character(0)
has_continuous <- FALSE
for (v in covar_vars) {
  raw <- coldata[[v]]
  # NA check where we can still name the offending sample.
  blank <- is.na(raw) | !nzchar(trimws(as.character(raw)))
  if (any(blank))
    die("design variable '", v, "' is empty/NA for sample(s): ",
        paste(coldata$sample[blank], collapse = ", "),
        "\n       Fill every row of that column, or drop the term from the design.")
  if (!is.null(num_covars) && v %in% num_covars) {
    n <- suppressWarnings(as.numeric(raw))
    if (any(is.na(n))) die("--numeric-covars lists '", v, "' but it has non-numeric values.")
    coldata[[v]] <- n; covar_types[v] <- "numeric"; has_continuous <- TRUE
  } else {
    coldata[[v]] <- factor(sanitize_level(as.character(raw)))
    covar_types[v] <- paste0("factor(", nlevels(coldata[[v]]), ": ",
                             paste(levels(coldata[[v]]), collapse = "/"), ")")
    if (nlevels(coldata[[v]]) < 2)
      die("design term '", v, "' has only one level -- it cannot be fit. Drop it.")
  }
}
if (!is.null(tech_col)) {
  if (!tech_col %in% colnames(coldata))
    die("--collapse-technical column '", tech_col, "' is not in the sample table.")
  coldata[[tech_col]] <- as.character(coldata[[tech_col]])
}

if (has_continuous)
  logmsg("NOTE: a continuous covariate is in the design (",
         paste(names(covar_types)[covar_types == "numeric"], collapse = ", "), "). ",
         "DESeq2 does NOT perform automatic outlier detection and replacement with ",
         "continuous variables in the design, because its within-group variance ",
         "estimate does not extend to them. Stage 3 still plots Cook's distances, but ",
         "treat that figure as manual inspection only -- nothing was auto-flagged.")

## ---------------------------------------------------------------------------
## 1.7  FULL-RANK / CONFOUNDING CHECK
## ---------------------------------------------------------------------------
# Report typing BEFORE the rank check: when a design is refused, HOW a variable
# was typed is usually the reason (a continuous variable left as a factor gets
# one level per sample and cannot be fit).
if (length(covar_types) > 0)
  logmsg("covariate typing: ", paste(sprintf("%s=%s", names(covar_types), covar_types), collapse = "; "),
         "   (numeric = a LINEAR term, not group effects)")

rank_info <- design_rank_report(coldata, design_formula)
cross_tabs <- list()
for (v in covar_vars)
  if (is.factor(coldata[[v]]))
    cross_tabs[[v]] <- table(coldata[[v]], coldata$condition, dnn = c(v, "condition"))

io <- make_io(outdir)

if (length(cross_tabs) > 0) {
  bal <- do.call(rbind, lapply(names(cross_tabs), function(v) {
    tb <- cross_tabs[[v]]
    data.frame(covariate = v, level = rownames(tb), as.data.frame.matrix(tb),
               balanced = all(tb > 0), check.names = FALSE, stringsAsFactors = FALSE)
  }))
  io$table(bal, "01b_covariate_by_condition.csv")
  for (v in names(cross_tabs)) {
    tb <- cross_tabs[[v]]
    for (rr in rownames(tb))
      logmsg("  balance ", v, "=", rr, ": ",
             paste(sprintf("%s=%d", colnames(tb), tb[rr, ]), collapse = "  "))
    if (any(tb == 0))
      logmsg("WARNING: covariate '", v, "' is NOT balanced across condition -- some cells ",
             "are empty. Condition and '", v, "' cannot be separated in those cells.")
  }
}

if (!rank_info$ok) {
  logmsg("design: ", paste(deparse(design_formula), collapse = ""),
         "  rank=", rank_info$rank, " of ", rank_info$ncol, " columns")
  die("the design is NOT full rank -- a covariate is confounded with condition.\n",
      "       DESeq2 cannot separate the two effects, so this model cannot be fit.\n",
      "       Options:\n",
      "         (a) drop the covariate (run without --batch / simplify --design), or\n",
      "         (b) rebalance so each covariate level contains BOTH conditions, or\n",
      "         (c) if the covariate marks TECHNICAL replicates of the same library,\n",
      "             use --collapse-technical <COL> instead of a batch term.\n",
      "       See the DESeq2 vignette, 'Model matrix not full rank'.")
}

## ---------------------------------------------------------------------------
## 1.8  Sample-sheet TEMPLATE  (removes the two-pass CSV ritual entirely)
## ---------------------------------------------------------------------------
# The sample names are known only after discovery, which is why hand-writing a
# --covariates CSV is error-prone: one wrong name and the join fails. Emitting
# the template with the names already filled in makes the round trip
# open -> type A/B -> pass back, with nothing to mistype.
tmpl <- data.frame(sample = coldata$sample,
                   condition = as.character(coldata$condition),
                   stringsAsFactors = FALSE)
tmpl[[if (is.null(batch_col)) "batch" else batch_col]] <-
  if (!is.null(batch_col) && batch_col %in% colnames(coldata))
    as.character(coldata[[batch_col]]) else ""
tmpl$path <- coldata$count_file
write.csv(tmpl, file.path(io$tab_dir, "01c_samplesheet_template.csv"), row.names = FALSE)
logmsg("template -> tables/01c_samplesheet_template.csv",
       if (all(tmpl[[3]] == "")) "  (batch column EMPTY -- fill it in and pass with --samplesheet)" else "")

## ---------------------------------------------------------------------------
## 1.9  Config table + state handoff
## ---------------------------------------------------------------------------
cfg <- data.frame(
  key = c("pipeline_version","stage","ref","counts_dir","counts_suffix","n_samples",
          "conditions","control(reference)","treatment","contrast_name","design",
          "covariate_types","subset_to_contrast","collapse_technical","threads","seed","outdir"),
  value = c(PIPELINE_VERSION, "1/setup", ref, counts_dir, counts_suffix, nrow(coldata),
            paste(levels(coldata$condition), collapse = "|"), control, treat, contrast_nm,
            paste(deparse(design_formula), collapse = ""),
            if (length(covar_types) == 0) "<none>"
              else paste(sprintf("%s=%s", names(covar_types), covar_types), collapse = "; "),
            subset_ctr, if (is.null(tech_col)) "<none>" else tech_col, threads, seed, outdir),
  stringsAsFactors = FALSE)
io$table(cfg, "00_run_config.csv")
io$table(coldata, "01_sample_table.csv")

logmsg("resolved ", nrow(coldata), " samples across: ", paste(levels(coldata$condition), collapse = ", "))
for (cc in levels(coldata$condition)) {
  mm <- coldata[coldata$condition == cc, , drop = FALSE]
  src <- if ("source_name" %in% colnames(mm)) mm$source_name else mm$sample
  logmsg("  '", cc, "' (", nrow(mm), "): ", paste(src, collapse = ", "))
}
logmsg("design ", paste(deparse(design_formula), collapse = ""),
       " is FULL RANK (", rank_info$rank, "/", rank_info$ncol, " columns)")
logmsg("contrast: ", treat, " vs ", control, "  (", contrast_nm, ")")

save_state(list(coldata = coldata, design_formula = design_formula,
                covar_vars = covar_vars, covar_types = covar_types,
                has_continuous = has_continuous,
                control = control, treat = treat, contrast_nm = contrast_nm,
                conds = levels(coldata$condition), tech_col = tech_col,
                ref = ref, counts_dir = counts_dir, counts_suffix = counts_suffix,
                threads = threads, seed = seed, outdir = outdir,
                args_raw = .args_raw),
           outdir, 1, opt)

if (dry_run) {
  logmsg("--dry-run: setup validated. Stopping before stage 2.")
  cat("\nDRY RUN OK.\n",
      "  samples : ", paste(coldata$sample, collapse = ", "), "\n",
      "  design  : ", paste(deparse(design_formula), collapse = ""), " (full rank)\n",
      "  contrast: ", treat, " vs ", control, "\n",
      "  template: ", file.path(io$tab_dir, "01c_samplesheet_template.csv"), "\n",
      "Remove --dry-run to run stages 2-4.\n", sep = "")
}
logmsg("STAGE 1 COMPLETE")
log_close()