#!/usr/bin/env Rscript
###############################################################################
# deseq2_common.R  --  SHARED LIBRARY for the 4-stage DESeq2 pipeline
#
# This file is NOT a stage. It is sourced by all four stage scripts and holds
# everything they have in common: the argument parser, the flag registry, the
# logger, the table/figure writers, the statistical predicates, and the
# stage-to-stage state handoff.
#
# You should not need to read this file to decide whether a PARAMETER is right.
# Every parameter is documented in the stage that uses it:
#
#   deseq2_01_setup.R    experiment definition: samples, groups, design, batch
#   deseq2_02_counts.R   count loading, annotation, library-composition QC
#   deseq2_03_model.R    dispersion, DESeq() fit, transforms, sample-level QC
#   deseq2_04_results.R  contrasts, shrinkage, thresholds, DEG tables, figures
#
# WHY THE PIPELINE IS SPLIT
# -------------------------
# Four ~400-line files with a documented header each, instead of one 1200-line
# file. Three concrete gains beyond readability:
#   1. RESUMABILITY. DESeq() is the slow step and lives in stage 3. Changing
#      --alpha or --lfc-report only re-runs stage 4 (seconds, not an hour).
#   2. FAILURE ISOLATION. A crash names the stage, so you know whether the
#      problem is your sample sheet, your annotation, the fit, or the contrast.
#   3. INSPECTABLE INTERMEDIATES. Each stage writes an .rds you can load in R.
#
# THE COST, AND HOW IT IS PAID
# ----------------------------
# Splitting means state crosses file boundaries, so a stale intermediate could
# silently be paired with fresh arguments. That is the one real hazard of this
# design, and it is closed by a STAGE SIGNATURE: each stage records a canonical
# string of exactly the flags that could change its output, and every stage
# refuses to run if its upstream signature does not match. Re-running stage 4
# with a new --alpha is allowed precisely because --alpha is not in stage 3's
# key set. Nothing is trusted implicitly.
###############################################################################

PIPELINE_VERSION <- "v3 (4-stage)"

## ===========================================================================
## FLAG REGISTRY
## Single source of truth. Adding a flag here makes it valid in every stage.
## `stage_keys` below decides which flags invalidate which stage.
## ===========================================================================
.bool_flags <- c("help", "h", "no-annotation", "strip-gene-version", "dry-run",
                 "collapse-to-gene", "no-independent-filtering",
                 "subset-to-contrast", "force")

.known_keys <- c(
  # stage 1 -- experiment definition
  "ref", "counts-dir", "counts-suffix", "samplesheet", "covariates",
  "group1", "group2", "group1-reps", "group2-reps",
  "control", "treatment", "batch", "batch-map", "design", "numeric-covars",
  "collapse-technical", "contrast-name", "subset-to-contrast",
  # stage 2 -- counts + annotation
  "gtf", "id-map", "no-annotation", "strip-gene-version", "biotype-filter",
  # stage 3 -- model
  "min-count", "min-samples", "test", "reduced", "fit-type", "transform",
  # stage 4 -- results
  "alpha", "lfc-report", "lfc-threshold", "alt-hypothesis", "rank-on",
  "shrink", "no-independent-filtering", "top-n",
  # global
  "outdir", "threads", "seed", "rlib", "collapse-to-gene", "dry-run",
  "stage", "force", "help", "h")

# Which flags could change a stage's OUTPUT. A stage's signature is built from
# its own keys plus every upstream stage's keys.
stage_keys <- list(
  s1 = c("ref","counts-dir","counts-suffix","samplesheet","covariates",
         "group1","group2","group1-reps","group2-reps","control","treatment",
         "batch","batch-map","design","numeric-covars","collapse-technical",
         "contrast-name","subset-to-contrast"),
  s2 = c("gtf","id-map","no-annotation","strip-gene-version","biotype-filter"),
  s3 = c("min-count","min-samples","test","reduced","fit-type","transform"),
  s4 = c("alpha","lfc-report","lfc-threshold","alt-hypothesis","rank-on",
         "shrink","no-independent-filtering","top-n")
)

## ===========================================================================
## ARGUMENT PARSING
## ===========================================================================
parse_args <- function(args) {
  opt <- list(); i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (!startsWith(a, "--") && !(a %in% c("-h")))
      stop("Unexpected argument (not a --flag): '", a, "'. See --help.", call. = FALSE)
    a2 <- sub("^--?", "", a); val <- NULL
    if (grepl("=", a2, fixed = TRUE)) {
      key <- sub("=.*$", "", a2); val <- sub("^[^=]*=", "", a2)
    } else key <- a2
    if (key %in% .bool_flags) { opt[[key]] <- TRUE; i <- i + 1L; next }
    if (is.null(val)) {
      if (i + 1L > length(args)) stop("Option --", key, " requires a value.", call. = FALSE)
      val <- args[[i + 1L]]; i <- i + 2L
    } else i <- i + 1L
    opt[[key]] <- val
  }
  opt
}

# --- accessors --------------------------------------------------------------
mk_getters <- function(opt) {
  list(
    chr = function(k, default = NULL) if (!is.null(opt[[k]])) as.character(opt[[k]]) else default,
    bool = function(k) isTRUE(opt[[k]]),
    # num/int VALIDATE. v2 let `--min-count ten` become NA and fail later inside
    # rowSums(counts(dds) >= NA) with an unreadable subsetting error.
    num = function(k, default, min = -Inf, max = Inf) {
      v <- opt[[k]]; if (is.null(v)) return(default)
      n <- suppressWarnings(as.numeric(v))
      if (is.na(n)) die("--", k, " must be a number (got '", v, "').")
      if (n < min || n > max) die("--", k, " must be between ", min, " and ", max, " (got ", n, ").")
      n
    },
    int = function(k, default, min = -Inf, max = Inf) {
      v <- opt[[k]]; if (is.null(v)) return(default)
      if (!grepl("^-?[0-9]+$", trimws(v))) die("--", k, " must be a whole number (got '", v, "').")
      n <- as.integer(v)
      if (n < min || n > max) die("--", k, " must be between ", min, " and ", max, " (got ", n, ").")
      n
    }
  )
}
split_csv <- function(s) if (is.null(s)) NULL else trimws(strsplit(s, ",", fixed = TRUE)[[1]])

## ===========================================================================
## LOGGING -- every stage APPENDS to one <outdir>/run.log
## ===========================================================================
.log_con   <- NULL
.log_path  <- NULL
.log_stage <- ""

# NOTE: do NOT reg.finalizer() here. reg.finalizer(environment(), ...) attaches
# to log_open's LOCAL frame; that frame is garbage-collected moments later, the
# finalizer fires, and the log connection closes mid-run. The console keeps
# printing while run.log silently stops -- and because it depends on GC timing
# it does not reproduce reliably. Instead the path is remembered and logmsg
# reopens on demand, so the log cannot be lost.
log_open <- function(outdir, stage, truncate = FALSE) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  .log_stage <<- stage
  .log_path  <<- file.path(outdir, "run.log")
  if (truncate) try(writeLines(character(0), .log_path), silent = TRUE)
  .log_con <<- tryCatch(file(.log_path, open = "at"), error = function(e) NULL)
  invisible(NULL)
}
log_close <- function() {
  if (!is.null(.log_con) && inherits(.log_con, "connection")) {
    try(close(.log_con), silent = TRUE)
  }
  .log_con <<- NULL
}
logmsg <- function(...) {
  line <- paste0("[", format(Sys.time(), "%H:%M:%S"), "][", .log_stage, "] ", ...)
  cat(line, "\n")
  # Reopen if the connection was lost for any reason; never drop log lines.
  if (is.null(.log_con) && !is.null(.log_path))
    .log_con <<- tryCatch(file(.log_path, open = "at"), error = function(e) NULL)
  if (!is.null(.log_con)) {
    ok <- tryCatch({ cat(line, "\n", file = .log_con); flush(.log_con); TRUE },
                   error = function(e) FALSE)
    if (!ok) {
      .log_con <<- NULL
      try({ cat(line, "\n", file = .log_path, append = TRUE) }, silent = TRUE)
    }
  }
}
logbanner <- function(title) {
  logmsg(strrep("=", 66)); logmsg(title); logmsg(strrep("=", 66))
}
die <- function(...) {
  msg <- paste0(...)
  if (!is.null(.log_con)) logmsg("ERROR: ", msg) else message("ERROR: ", msg)
  log_close(); quit(status = 2)
}
has <- function(pkg) requireNamespace(pkg, quietly = TRUE)

## ===========================================================================
## TABLE / FIGURE WRITERS  (one closure per output directory)
## ===========================================================================
make_io <- function(level_dir) {
  fig_dir <- file.path(level_dir, "figures"); tab_dir <- file.path(level_dir, "tables")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

  save_table <- function(df, name) {
    df <- as.data.frame(df, stringsAsFactors = FALSE, optional = TRUE)
    listcols <- names(df)[vapply(df, is.list, logical(1))]
    for (j in listcols)
      df[[j]] <- vapply(df[[j]], function(x) paste(as.character(x), collapse = ";"), character(1))
    stem <- sub("\\.csv$", "", name)
    write.csv(df, file.path(tab_dir, paste0(stem, ".csv")), row.names = FALSE)
    # TSV twin: opened via Excel's import dialog it will not convert SEPT9 /
    # MARCH1 / DEC1 into dates, and downstream tools prefer tabs.
    write.table(df, file.path(tab_dir, paste0(stem, ".tsv")), sep = "\t",
                row.names = FALSE, quote = FALSE, na = "NA")
    logmsg("table  -> ", file.path(basename(level_dir), "tables", paste0(stem, ".{csv,tsv}")))
  }
  save_base <- function(name, fun, width = 7, height = 6) {
    tryCatch({
      pdf(file.path(fig_dir, paste0(name, ".pdf")), width = width, height = height); fun(); dev.off()
      png(file.path(fig_dir, paste0(name, ".png")), width = width*150, height = height*150, res = 150); fun(); dev.off()
      logmsg("figure -> ", name, ".{pdf,png}")
    }, error = function(e) {
      while (grDevices::dev.cur() > 1L) grDevices::dev.off()
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
  list(fig_dir = fig_dir, tab_dir = tab_dir, dir = level_dir,
       table = save_table, base = save_base, gg = save_gg)
}

## ===========================================================================
## STAGE STATE  --  the handoff, and the guard that makes splitting safe
## ===========================================================================
state_dir <- function(outdir) file.path(outdir, "_state")

# Canonical string of the flags that could change this stage's output.
# Sorted, so flag ORDER on the command line never causes a false mismatch.
signature_for <- function(opt, upto) {
  keys <- unlist(stage_keys[seq_len(upto)], use.names = FALSE)
  vals <- vapply(sort(keys), function(k) {
    v <- opt[[k]]
    if (is.null(v)) "" else if (isTRUE(v)) "TRUE" else as.character(v)
  }, character(1))
  paste(paste0(sort(keys), "=", vals), collapse = ";")
}

save_state <- function(obj, outdir, stage_num, opt) {
  d <- state_dir(outdir); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  obj$.meta <- list(stage = stage_num,
                    signature = signature_for(opt, stage_num),
                    version = PIPELINE_VERSION,
                    time = format(Sys.time()),
                    r_version = R.version.string)
  f <- file.path(d, sprintf("%02d_state.rds", stage_num))
  saveRDS(obj, f)
  logmsg("state  -> _state/", basename(f))
  invisible(f)
}

load_state <- function(outdir, stage_num, opt, force = FALSE) {
  f <- file.path(state_dir(outdir), sprintf("%02d_state.rds", stage_num))
  if (!file.exists(f))
    die("stage ", stage_num, " output not found: ", f, "\n",
        "       Run the earlier stage(s) first. The submitter does this for you:\n",
        "         sbatch submit_rnaseq_deseq2_dge.sh <your flags>")
  st <- readRDS(f)
  want <- signature_for(opt, stage_num)
  if (!identical(st$.meta$signature, want)) {
    if (force) {
      logmsg("WARNING: --force set; using stage ", stage_num, " output that was ",
             "produced with DIFFERENT arguments. Results may not match your flags.")
    } else {
      old <- strsplit(st$.meta$signature, ";", fixed = TRUE)[[1]]
      new <- strsplit(want, ";", fixed = TRUE)[[1]]
      diffs <- setdiff(union(old, new), intersect(old, new))
      die("stage ", stage_num, " output on disk was produced with different arguments.\n",
          "       Changed flag(s): ", paste(head(diffs, 12), collapse = "  "), "\n",
          "       Re-run from stage ", stage_num, " (the submitter does this by default),\n",
          "       or pass --force to reuse the existing output anyway.")
    }
  }
  st
}

## ===========================================================================
## STATISTICAL HELPERS
## ===========================================================================

# ONE definition of "clears the fold-change cutoff", derived from
# --alt-hypothesis. Used by the DEG table, the summary, the volcano and the
# biotype breakdown so they can never disagree.
# Vignette: greaterAbs |b|>x (two-tailed); lessAbs |b|<x (p = max of the upper
# and lower tests); greater b>x; less b<-x.
make_lfc_pass <- function(alt, thr) {
  if (is.na(thr) || thr <= 0) return(function(x) !is.na(x))
  switch(alt,
    greaterAbs = function(x) !is.na(x) & abs(x) >= thr,
    lessAbs    = function(x) !is.na(x) & abs(x) <   thr,
    greater    = function(x) !is.na(x) & x     >=  thr,
    less       = function(x) !is.na(x) & x     <= -thr,
    stop("unknown altHypothesis: ", alt))
}

# Design-matrix rank check. Vignette, "Model matrix not full rank": when two
# variables carry the same information the software cannot fit an effect for
# both and a unique solution of coefficients is not possible.
design_rank_report <- function(cd, fml) {
  mm <- tryCatch(stats::model.matrix(fml, data = cd), error = function(e) NULL)
  if (is.null(mm)) return(list(ok = FALSE, reason = "model matrix could not be built",
                               rank = NA_integer_, ncol = NA_integer_))
  r <- qr(mm)$rank
  list(ok = r >= ncol(mm), reason = if (r >= ncol(mm)) "full rank" else "rank deficient",
       rank = r, ncol = ncol(mm), mm = mm)
}

# DESeq2 factor levels may contain only letters, numbers, underscore and period.
# DESeq2 renames anything else itself, after which a contrast built from the
# ORIGINAL string matches no level and results() dies -- an hour into the run.
sanitize_level <- function(s) {
  out <- gsub("[^A-Za-z0-9._]+", "_", s)
  out <- sub("^_+", "", sub("_+$", "", out))
  out[!nzchar(out)] <- "level"
  ifelse(grepl("^[0-9]", out), paste0("X", out), out)
}

read_delim_auto <- function(path) {
  sep <- if (grepl("\\.tsv$|\\.txt$", path, ignore.case = TRUE)) "\t" else ","
  read.table(path, header = TRUE, sep = sep, stringsAsFactors = FALSE,
             quote = "\"", check.names = FALSE, comment.char = "")
}

## ===========================================================================
## PROVENANCE
## ===========================================================================
git_info <- function(script_path) {
  out <- list(commit = "NA", dirty = NA)
  if (is.na(script_path) || !nzchar(Sys.which("git"))) return(out)
  gd <- dirname(script_path)
  h <- suppressWarnings(tryCatch(system2("git", c("-C", gd, "rev-parse", "HEAD"),
                                         stdout = TRUE, stderr = FALSE),
                                 error = function(e) character(0)))
  if (length(h) >= 1 && nzchar(h[1])) {
    out$commit <- h[1]
    st <- suppressWarnings(tryCatch(system2("git", c("-C", gd, "status", "--porcelain"),
                                            stdout = TRUE, stderr = FALSE),
                                    error = function(e) character(0)))
    out$dirty <- length(st) > 0
  }
  out
}

this_script_path <- function() {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) normalizePath(sub("^--file=", "", a[1]), mustWork = FALSE) else NA_character_
}

# Locate deseq2_common.R relative to the calling stage script.
source_common <- function() invisible(NULL)   # placeholder; stages source directly