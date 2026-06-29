#!/usr/bin/env Rscript
# =============================================================================
# build_sample_sheet.R  —  CONSOLIDATED sample-sheet builder
# =============================================================================
# Merges the former two-stage pipeline (generate_manifest.R  ->  build_sample_sheet.R)
# into ONE script that prints ALL the diagnostics of both stages but writes only
# ONE file: the analysis sample sheet consumed by run_csaw.R / run_diffbind.R.
#
# WHY THE TWO STAGES WERE MERGED (and what is preserved)
#   The old manifest stage existed mainly as a HUMAN PROOFREADING CHECKPOINT:
#   it discovered files on disk (handling folder!=filename naming and spike-in
#   alignment failures), then left Batch and IgG BLANK for you to fill before the
#   sheet was built. That checkpoint is still available here — it is just no
#   longer a separate file/invocation. You can:
#     * run --root for a one-shot disk discovery -> sheet (all flags printed), or
#     * run --manifest <proofread.csv> to build from a manifest you curated, or
#     * generate the sheet, edit its Batch/IgG columns by hand, and re-run with
#       --merge_existing <that sheet> to re-resolve the IgG control BAMs.
#   Nothing is ever silently zeroed: name mismatches and spike-in failures are
#   surfaced as [ACTION NEEDED]; rows are always written so you can eyeball them.
#
# INPUT MODES (exactly one of):
#   (A) --root raw_dir [--recursive]
#       Walk every sample folder, DISCOVER the real analysis BAM inside (its inner
#       name may disagree with the folder), the peak file, and the E. coli spike-in
#       bowtie2 log; take biological labels (Condition/Treatment/Factor/Replicate)
#       from the FOLDER name (the corrected source of truth). Batch defaults to
#       'batch1'; IgG is blank unless supplied via --merge_existing.
#   (B) --manifest manifest.csv
#       Build from a proofread manifest (Batch/IgG already filled). Inherits the
#       spike-in counts and any name-mismatch flags recorded in the manifest.
#
# CASE-INSENSITIVE biological labels: 'Parental'=='parental', 'H3K27Me3'=='H3K27me3'.
# BATCH & IgG are OPTIONAL. spikein_reads may be empty for samples whose spike-in
#   alignment FAILED — run_csaw.R only errors if such a sample is in a contrast
#   run with --spikein sheet.
#
# SINGLE OUTPUT FILE — the sample sheet, columns (in this order):
#   SampleID, Condition, Treatment, Factor, Replicate, Batch, IgG,
#   bamReads, ControlID, bamControl, Peaks, PeakCaller,
#   spikein_reads, spikein_pct, spikein_status
#
# -----------------------------------------------------------------------------
# USAGE EXAMPLES
# -----------------------------------------------------------------------------
#   # See every flag:
#   Rscript build_sample_sheet.R --help
#
#   # (A) One-shot disk discovery -> sheet (prints mismatch + spike-in QC):
#   Rscript build_sample_sheet.R \
#       --root data/raw/raw_dirs/ \
#       --out  data/processed/sample_sheet.csv
#
#   # (A) Carry over manually-filled Batch/IgG from a prior sheet/manifest, and
#   #     resolve the IgG control BAMs from a separate IgG tree:
#   Rscript build_sample_sheet.R \
#       --root data/raw/raw_dirs/ \
#       --merge_existing data/processed/sample_sheet.csv \
#       --igg_root data/raw/raw_dirs/ \
#       --require_igg \
#       --out data/processed/sample_sheet.csv
#
#   # (B) Build from a proofread manifest (Batch/IgG already filled in it):
#   Rscript build_sample_sheet.R \
#       --manifest data/raw/manifest.csv \
#       --out data/processed/sample_sheet.csv
#
#   # Narrow-mark project where H3K27Ac uses MACS2 narrow peaks, BordR8 present:
#   Rscript build_sample_sheet.R --root data/raw/raw_dirs/ \
#       --treatments UNC,EPZ,DMSO --narrow_marks H3K27AC --out sheet.csv
#
# AFTER RUNNING
#   1. Read the [ACTION NEEDED] blocks (name mismatches AND spike-in failures).
#   2. If Batch/IgG are blank and you need them: fill them in the OUTPUT sheet,
#      then re-run with --merge_existing <sheet> to resolve the control BAMs.
#   3. Inspect the printed Factor x Condition table before running analyses.
# =============================================================================

suppressWarnings(suppressMessages(library(optparse)))

# =============================================================================
# CLI
# =============================================================================
opt <- parse_args(OptionParser(option_list = list(
  make_option("--root", type = "character", default = NULL,
              help = "Disk-discovery mode: raw_dir containing ALL sample folders."),
  make_option("--manifest", type = "character", default = NULL,
              help = "Manifest mode: proofread manifest.csv (Batch/IgG filled)."),
  make_option("--recursive", action = "store_true", default = FALSE,
              help = "[--root] Recurse into sub-directories when listing sample folders."),
  make_option("--out", type = "character", default = "sample_sheet.csv",
              help = "Output sample sheet CSV path. [default %default]"),
  make_option("--bam", type = "character", default = "dedup",
              help = "Which BAM dir to record: dedup|dup. [default %default]"),
  make_option("--igg_root", type = "character", default = NULL,
              help = "Dir containing IgG folders named in the IgG column (default: --root, or the manifest's dir)."),
  make_option("--treatments", type = "character", default = "UNC,EPZ,DMSO",
              help = "Treatment tokens to detect in names. [default %default]"),
  make_option("--narrow_marks", type = "character", default = "H3K27AC",
              help = "Marks that use MACS2 narrow peaks. [default %default]"),
  make_option("--sample_regex", type = "character",
              default = "^(Bottom10|Bottom|Control|Parental|B10|B8|Bottom8)-(?:UNC-|EPZ-|DMSO-)?([A-Za-z0-9]+)-R([0-9]+)$",
              help = "Regex parsing Condition / Mark / Replicate from a folder (or BAM) name."),
  make_option("--igg_regex", type = "character", default = "^IgG-R([0-9]+)$",
              help = "[--root] Folders matching this are treated as IgG controls (not samples)."),
  make_option("--bam_suffix_strip", type = "character",
              default = "\\.(dedup|dup\\.marked|sorted|filtered|markdup|rmdup)$",
              help = "Regex of trailing pipeline suffixes stripped from a BAM basename before the mismatch check."),
  make_option("--require_igg", action = "store_true", default = FALSE,
              help = "Error (instead of warn) when a named IgG control BAM cannot be found."),
  make_option("--merge_existing", type = "character", default = NULL,
              help = "Optional prior sheet/manifest: carry over Batch & IgG by SampleID/folder.")
)))

if (is.null(opt$root) && is.null(opt$manifest))
  stop("Provide --root (disk discovery) or --manifest (proofread manifest).")
if (!is.null(opt$root) && !is.null(opt$manifest))
  stop("Provide only ONE of --root or --manifest, not both.")
if (!opt$bam %in% c("dedup", "dup")) stop("--bam must be 'dedup' or 'dup'.")

# Resolve input roots to ABSOLUTE paths so the discovered bamReads/Peaks in the
# sheet are absolute and resolve from any working directory the downstream tools
# (csaw.R, run_chipseqspikeinfree.R, ...) are launched from. Without this, a sheet
# built with '--root ./' stores './/<folder>/...' paths that only resolve from the
# directory the sheet was built in.
if (!is.null(opt$root))     opt$root     <- normalizePath(opt$root,     mustWork = FALSE)
if (!is.null(opt$igg_root)) opt$igg_root <- normalizePath(opt$igg_root, mustWork = FALSE)

treatments   <- toupper(trimws(strsplit(opt$treatments, ",")[[1]])); treatments <- treatments[treatments != ""]
narrow_marks <- toupper(trimws(strsplit(opt$narrow_marks, ",")[[1]]))
bam_subdir   <- if (opt$bam == "dup") "dup.marked" else "dedup"

# Expand the optional "(?:UNC-|EPZ-|DMSO-)?" treatment alternation in the regex so
# it tracks --treatments (keeps folder/BAM parsing in sync with the token list).
treat_alt <- paste(treatments, collapse = "|")
sample_re <- opt$sample_regex
if (grepl("\\(\\?:UNC-\\|EPZ-\\|DMSO-\\)\\?", sample_re) && nzchar(treat_alt)) {
  sample_re <- sub("\\(\\?:UNC-\\|EPZ-\\|DMSO-\\)\\?", paste0("(?:(?:", treat_alt, ")-)?"), sample_re)
} else if (grepl("\\(\\?:UNC-\\|EPZ-\\)\\?", sample_re) && nzchar(treat_alt)) {
  sample_re <- sub("\\(\\?:UNC-\\|EPZ-\\)\\?", paste0("(?:(?:", treat_alt, ")-)?"), sample_re)
}

# =============================================================================
# Label helpers (case-insensitive; identical semantics to the original scripts)
# =============================================================================
blank <- function(x) is.na(x) | trimws(as.character(x)) == ""

# Canonical Condition value written to the sheet (lowercase 'parental', 'B10', 'B8').
cond_norm <- function(x) {
  x <- tolower(trimws(x))
  if (x %in% c("control", "parental"))            "parental"
  else if (x %in% c("bottom", "bottom10", "b10")) "B10"
  else if (x %in% c("bottom8", "b8"))             "B8"
  else x
}
# Canonicalize a mark token. Histone marks -> e.g. 'H3K27me3'; non-histone tokens
# (e.g. a cosmetic 'Migrated' Factor) are returned upper-cased UNCHANGED (no
# mangling), so run_csaw.R --mark matches case-insensitively.
canon_mark <- function(m) {
  m <- toupper(trimws(m))
  if (!grepl("^H[0-9]+K[0-9]+", m)) return(m)
  core <- sub("^(H[0-9]+K[0-9]+).*$", "\\1", m)
  rest <- tolower(sub("^H[0-9]+K[0-9]+", "", m))
  paste0(core, rest)
}
# Treatment token found between dashes in a name (UNC/EPZ/DMSO/...), else 'none'.
treat_of <- function(name) {
  hit <- treatments[vapply(treatments, function(tk) grepl(paste0("-", tk, "-"), toupper(name)), logical(1))]
  if (length(hit) >= 1) hit[1] else "none"
}
# Parse Condition/Treatment/Mark/Replicate out of a folder (or BAM) name.
parse_fields <- function(name) {
  if (is.na(name) || !nzchar(name)) return(NULL)
  if (nzchar(opt$bam_suffix_strip)) name <- sub(opt$bam_suffix_strip, "", name)
  if (!grepl(sample_re, name, perl = TRUE)) return(NULL)
  list(cond  = cond_norm(sub(sample_re, "\\1", name, perl = TRUE)),
       treat = treat_of(name),
       mark  = canon_mark(sub(sample_re, "\\2", name, perl = TRUE)),
       rep   = as.integer(sub(sample_re, "\\3", name, perl = TRUE)))
}

# =============================================================================
# Disk-discovery helpers (folded in from generate_manifest.R)
# =============================================================================
# Find the real analysis BAM inside a folder (its inner name may differ from the
# folder name). Returns base (no .bam) + full path, or NA/NA if none.
find_bam <- function(folder_path) {
  d <- file.path(folder_path, "aligned", bam_subdir)
  if (!dir.exists(d)) return(list(base = NA, path = NA))
  bams <- list.files(d, pattern = "\\.bam$", full.names = TRUE)
  bams <- bams[!grepl("\\.bai$", bams)]
  if (length(bams) == 0) return(list(base = NA, path = NA))
  if (length(bams) > 1) warning("Multiple BAMs in ", d, " — picking first: ", basename(bams[1]))
  list(base = sub("\\.bam$", "", basename(bams[1])), path = bams[1])
}
# Find the peak file (MACS2 narrow for narrow marks, broad otherwise), preferring
# a blacklist_filtered copy if present.
find_peak <- function(folder_path, is_narrow) {
  sub  <- if (is_narrow) "macs2.narrow" else "macs2.broad"
  patt <- if (is_narrow) "_peaks\\.narrowPeak$" else "_peaks\\.broadPeak$"
  base <- file.path(folder_path, "peakcalling", sub)
  if (!dir.exists(base)) return(NA)
  cands <- c(list.files(file.path(base, "blacklist_filtered"), patt, full.names = TRUE),
             list.files(base, patt, full.names = TRUE),
             list.files(base, patt, full.names = TRUE, recursive = TRUE))
  cands <- cands[!is.na(cands)]
  if (length(cands) == 0) NA else cands[1]
}
# Locate logs/*.spikein.bowtie2 for a folder.
find_spikein <- function(folder_path) {
  d <- file.path(folder_path, "logs")
  if (!dir.exists(d)) return(NA_character_)
  f <- list.files(d, pattern = "\\.spikein\\.bowtie2$", full.names = TRUE)
  if (length(f) == 0) return(NA_character_)
  f[1]
}
# Parse a bowtie2 spike-in alignment summary. spikein_reads = total_pairs * rate%.
# Failure modes: NOT_FOUND | FAILED_ALIGNMENT | UNPARSED | OK.
parse_spikein_log <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path))
    return(list(status = "NOT_FOUND", reads = NA, pct = NA, total = NA))
  txt <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
  if (length(txt) == 0) return(list(status = "UNPARSED", reads = NA, pct = NA, total = NA))
  failed   <- any(grepl("is not a Bowtie 2 index|does not exist|Exiting now|^\\(ERR\\)", txt, perl = TRUE))
  has_rate <- any(grepl("overall alignment rate", txt))
  if (failed && !has_rate) return(list(status = "FAILED_ALIGNMENT", reads = NA, pct = NA, total = NA))
  tline <- txt[grepl("reads; of these", txt)]
  rline <- txt[grepl("overall alignment rate", txt)]
  if (length(tline) == 0 || length(rline) == 0)
    return(list(status = "UNPARSED", reads = NA, pct = NA, total = NA))
  total <- suppressWarnings(as.numeric(gsub("[^0-9]", "", sub("reads;.*", "", tline[1]))))
  pct   <- suppressWarnings(as.numeric(sub("%.*", "", trimws(rline[1]))))
  if (is.na(total) || is.na(pct)) return(list(status = "UNPARSED", reads = NA, pct = NA, total = NA))
  list(status = "OK", reads = round(total * pct / 100), pct = pct, total = total)
}
# Does the inner BAM name disagree with the folder's biological identity?
field_mismatch <- function(folder, bam_base, f_cond, f_treat, f_mark, f_rep) {
  if (is.na(bam_base) || !nzchar(bam_base)) return("NO_BAM_FOUND")
  bf <- parse_fields(bam_base)
  if (is.null(bf))
    return(paste0("CHECK: bam='", bam_base, "' does not parse into Condition/Treatment/Mark/Replicate (folder='", folder, "')"))
  diffs <- c()
  if (!identical(tolower(bf$cond),  tolower(f_cond)))  diffs <- c(diffs, paste0("Condition(folder=", f_cond, ",bam=", bf$cond, ")"))
  if (!identical(tolower(bf$treat), tolower(f_treat))) diffs <- c(diffs, paste0("Treatment(folder=", f_treat, ",bam=", bf$treat, ")"))
  if (!identical(tolower(bf$mark),  tolower(f_mark)))  diffs <- c(diffs, paste0("Factor(folder=", f_mark, ",bam=", bf$mark, ")"))
  if (!isTRUE(bf$rep == f_rep))                        diffs <- c(diffs, paste0("Replicate(folder=", f_rep, ",bam=", bf$rep, ")"))
  if (length(diffs) == 0) "" else paste0("CHECK: ", paste(diffs, collapse = "; "))
}

# =============================================================================
# IgG control-BAM resolution (folded in from build_sample_sheet.R)
# Given an IgG FOLDER name, find a control BAM under igg_root/<folder>/aligned/*.
# =============================================================================
resolve_control_bam <- function(igg_folder, igg_root, require_igg, sample_folder) {
  if (blank(igg_folder)) return(list(ControlID = NA, bamControl = NA))
  igg_folder <- trimws(as.character(igg_folder))
  for (sd in c("dedup", "dup.marked")) {
    d <- file.path(igg_root, igg_folder, "aligned", sd)
    if (dir.exists(d)) {
      b <- list.files(d, "\\.bam$", full.names = TRUE); b <- b[!grepl("\\.bai$", b)]
      if (length(b) >= 1) return(list(ControlID = igg_folder, bamControl = b[1]))
    }
  }
  msg <- paste0("IgG '", igg_folder, "' for ", sample_folder, " not found under ", igg_root)
  if (isTRUE(require_igg)) stop(msg) else warning(msg, " — bamControl blank.")
  list(ControlID = igg_folder, bamControl = NA)
}

# =============================================================================
# --merge_existing loader: carry manually-filled Batch / IgG by SampleID/folder.
# =============================================================================
prior <- NULL
if (!is.null(opt$merge_existing) && file.exists(opt$merge_existing)) {
  prior <- read.csv(opt$merge_existing, stringsAsFactors = FALSE)
  # accept either a sheet (SampleID) or a manifest (folder) as the key column
  prior$.key <- if ("SampleID" %in% names(prior)) prior$SampleID
                else if ("folder" %in% names(prior)) prior$folder else NA
  message("[INFO] Will carry over Batch/IgG by SampleID/folder from ", opt$merge_existing)
}
get_prior <- function(key, col) {
  if (is.null(prior) || !(col %in% names(prior))) return("")
  v <- prior[[col]][!is.na(prior$.key) & prior$.key == key]
  if (length(v) >= 1 && !is.na(v[1])) as.character(v[1]) else ""
}

# =============================================================================
# STAGE 1 — build the internal record set (the in-memory "manifest")
#   Columns: folder, bam_basename, Condition, Treatment, Factor, Replicate,
#            Batch, IgG, bam_path, peak_path, peak_caller,
#            spikein_reads, spikein_pct, spikein_status, spikein_logfile, name_mismatch
# =============================================================================
igg_root <- if (!is.null(opt$igg_root)) opt$igg_root else if (!is.null(opt$root)) opt$root else dirname(opt$manifest)

igg_folders_seen <- character(0)   # for the [INFO] hint in --root mode

if (!is.null(opt$root)) {
  # ---------------------------------------------------------------------------
  # MODE A: walk the raw_dir and DISCOVER everything from disk.
  # ---------------------------------------------------------------------------
  if (opt$recursive) {
    ap <- list.dirs(opt$root, full.names = TRUE, recursive = TRUE); ap <- ap[basename(ap) != ""]
  } else {
    d  <- list.dirs(opt$root, full.names = FALSE, recursive = FALSE); d <- d[d != ""]
    ap <- file.path(opt$root, d)
  }
  rows <- list()
  for (fp in ap) {
    folder <- basename(fp)
    if (grepl(opt$igg_regex, folder, perl = TRUE)) { igg_folders_seen <- c(igg_folders_seen, folder); next }
    if (!grepl(sample_re, folder, perl = TRUE)) { warning("Folder did not match sample pattern, skipped: ", folder); next }

    cond <- cond_norm(sub(sample_re, "\\1", folder, perl = TRUE))
    mark <- canon_mark(sub(sample_re, "\\2", folder, perl = TRUE))
    rep  <- as.integer(sub(sample_re, "\\3", folder, perl = TRUE))
    trt  <- tolower(treat_of(folder))

    bam   <- find_bam(fp)
    isnar <- toupper(mark) %in% narrow_marks
    peak  <- find_peak(fp, isnar)
    mism  <- field_mismatch(folder, bam$base, cond, trt, mark, rep)

    spk_path <- find_spikein(fp); spk <- parse_spikein_log(spk_path)
    spk_file <- if (is.na(spk_path)) "" else sub("\\.spikein\\.bowtie2$", "", basename(spk_path))

    rows[[length(rows) + 1]] <- data.frame(
      folder = folder, bam_basename = ifelse(is.na(bam$base), "", bam$base),
      Condition = cond, Treatment = trt, Factor = mark, Replicate = rep,
      Batch = get_prior(folder, "Batch"), IgG = get_prior(folder, "IgG"),
      bam_path = ifelse(is.na(bam$path), "", bam$path),
      peak_path = ifelse(is.na(peak), "", peak), peak_caller = ifelse(isnar, "narrow", "broad"),
      spikein_reads = ifelse(is.na(spk$reads), "", spk$reads),
      spikein_pct   = ifelse(is.na(spk$pct),   "", spk$pct),
      spikein_status = spk$status, spikein_logfile = spk_file,
      name_mismatch = mism, stringsAsFactors = FALSE)
  }
  if (!length(rows)) stop("No sample folders matched under --root.")
  man <- do.call(rbind, rows)

} else {
  # ---------------------------------------------------------------------------
  # MODE B: read a proofread manifest; re-canonicalize labels for the sheet.
  # ---------------------------------------------------------------------------
  raw <- read.csv(opt$manifest, stringsAsFactors = FALSE)
  orig_cols <- colnames(raw)
  req <- c("folder", "bam_basename", "Condition", "Treatment", "Factor", "Replicate")
  miss <- setdiff(req, colnames(raw))
  if (length(miss)) stop("Manifest missing columns: ", paste(miss, collapse = ", "))
  if (!"bam_path" %in% colnames(raw) || any(blank(raw$bam_path)))
    stop("Manifest needs a non-empty 'bam_path' for every row. Re-generate the manifest.")
  for (cc in c("peak_path", "peak_caller", "Batch", "IgG", "spikein_reads", "spikein_pct",
               "spikein_status", "spikein_logfile", "name_mismatch"))
    if (!cc %in% colnames(raw)) raw[[cc]] <- if (cc == "spikein_status") "NA" else ""
  if (!"spikein_reads" %in% orig_cols)
    message("[NOTE] manifest had no spikein_reads column; sheet spikein_reads will be blank.")

  man <- data.frame(
    folder = raw$folder, bam_basename = raw$bam_basename,
    Condition = vapply(raw$Condition, cond_norm, character(1)),
    Treatment = tolower(trimws(raw$Treatment)),
    Factor = vapply(raw$Factor, canon_mark, character(1)),
    Replicate = suppressWarnings(as.integer(raw$Replicate)),
    Batch = as.character(raw$Batch),
    IgG = as.character(raw$IgG),
    bam_path = as.character(raw$bam_path),
    peak_path = as.character(raw$peak_path),
    peak_caller = as.character(raw$peak_caller),
    spikein_reads = as.character(raw$spikein_reads),
    spikein_pct = as.character(raw$spikein_pct),
    spikein_status = as.character(raw$spikein_status),
    spikein_logfile = as.character(raw$spikein_logfile),
    name_mismatch = as.character(raw$name_mismatch),
    stringsAsFactors = FALSE)
  # carry Batch/IgG from --merge_existing where the manifest left them blank
  for (i in seq_len(nrow(man))) {
    if (blank(man$Batch[i])) { pv <- get_prior(man$folder[i], "Batch"); if (nzchar(pv)) man$Batch[i] <- pv }
    if (blank(man$IgG[i]))   { pv <- get_prior(man$folder[i], "IgG");   if (nzchar(pv)) man$IgG[i]   <- pv }
  }
}

# Default Batch -> 'batch1' where still blank (after any merge).
man$Batch[blank(man$Batch)] <- "batch1"
# derive peak_caller if absent
need_caller <- blank(man$peak_caller)
if (any(need_caller))
  man$peak_caller[need_caller] <- ifelse(toupper(man$Factor[need_caller]) %in% narrow_marks, "narrow", "broad")

# =============================================================================
# STAGE 1 DIAGNOSTICS  (generate_manifest.R style) — printed, never written
# =============================================================================
message("[INFO] Discovered ", nrow(man), " sample(s) (", 
        if (!is.null(opt$root)) "disk-discovery / --root" else "manifest / --manifest", " mode).")
message("[INFO] Factor/Treatment x Condition:")
print(table(paste(man$Factor, man$Treatment, sep = "/"), man$Condition))

# name mismatches (inner BAM disagrees with folder, or unparseable)
flagged <- man[!blank(man$name_mismatch), c("folder", "bam_basename", "name_mismatch")]
if (nrow(flagged) > 0) {
  message("\n[ACTION NEEDED] ", nrow(flagged),
          " row(s) where the inner BAM name disagrees with the folder (or could not be parsed) — VERIFY:")
  print(flagged, row.names = FALSE)
} else message("\n[OK] No folder/BAM naming mismatches detected.")

# spike-in failures
spk_bad <- man[!blank(man$spikein_status) & !man$spikein_status %in% c("OK", "NA"),
               c("folder", "spikein_status", "spikein_logfile")]
if (nrow(spk_bad) > 0) {
  message("\n[ACTION NEEDED] ", nrow(spk_bad), " sample(s) without a usable spike-in count:")
  message("   FAILED_ALIGNMENT = bowtie2 errored (e.g. missing E. coli index); re-run vs a valid index.")
  message("   NOT_FOUND/UNPARSED = no parseable .spikein.bowtie2 in logs/.")
  print(spk_bad, row.names = FALSE)
} else message("[OK] All samples have a parseable spike-in status.")

# spike-in log-name vs folder mismatch (informational)
spk_mm <- man[!blank(man$spikein_logfile) & man$spikein_logfile != man$folder, c("folder", "spikein_logfile")]
if (nrow(spk_mm) > 0) {
  message("\n[INFO] spike-in log filename differs from folder for ", nrow(spk_mm),
          " sample(s) (recorded by FOLDER; verify identity):")
  print(spk_mm, row.names = FALSE)
}

# spike-in QC: usable window + cross-sample spread (the carryover guard)
ok <- man[man$spikein_status == "OK", ]
if (nrow(ok) > 0) {
  pct <- suppressWarnings(as.numeric(ok$spikein_pct))
  oor <- ok$folder[!is.na(pct) & (pct < 0.2 | pct > 5)]
  if (length(oor) > 0)
    message("\n[FLAG] ", length(oor), " sample(s) outside the 0.2-5% usable spike-in window (EpiCypher); ",
            "scale factors suspect: ", paste(oor, collapse = ", "))
  if (sum(!is.na(pct)) > 1) {
    sp <- range(pct, na.rm = TRUE); ratio <- if (sp[1] > 0) sp[2] / sp[1] else Inf
    message(sprintf("[INFO] spike-in fraction span %.2f%%-%.2f%% (max/min = %.0fx).", sp[1], sp[2], ratio))
    if (ratio > 5)
      message("[FLAG] >5x spread: a single spike-in scalar across all samples is NOT defensible — ",
              "group by batch (matched cells + beads) and normalize WITHIN batch.")
  }
}

# IgG folder hint (--root mode): show what could be filled into the IgG column
if (!is.null(opt$root) && length(unique(igg_folders_seen)) > 0)
  message("\n[INFO] IgG control folder(s) present under --root: ",
          paste(sort(unique(igg_folders_seen)), collapse = ", "),
          "\n       (fill the sheet's IgG column with one of these, then re-run with --merge_existing to resolve bamControl.)")

# =============================================================================
# STAGE 2 — transform records into the analysis SHEET (resolve IgG control BAMs)
# =============================================================================
sheet_rows <- list()
for (i in seq_len(nrow(man))) {
  r <- man[i, ]
  ctrl <- resolve_control_bam(r$IgG, igg_root, opt$require_igg, r$folder)
  sp_reads <- suppressWarnings(as.numeric(r$spikein_reads)); if (is.na(sp_reads)) sp_reads <- NA
  sp_pct   <- suppressWarnings(as.numeric(r$spikein_pct));   if (is.na(sp_pct))   sp_pct   <- NA
  sheet_rows[[i]] <- data.frame(
    SampleID = r$folder, Condition = r$Condition, Treatment = r$Treatment,
    Factor = r$Factor, Replicate = as.integer(r$Replicate), Batch = r$Batch,
    IgG = ifelse(blank(r$IgG), NA, as.character(r$IgG)),
    bamReads = as.character(r$bam_path),
    ControlID = ctrl$ControlID, bamControl = ctrl$bamControl,
    Peaks = ifelse(blank(r$peak_path), NA, as.character(r$peak_path)),
    PeakCaller = r$peak_caller,
    spikein_reads = sp_reads, spikein_pct = sp_pct,
    spikein_status = ifelse(blank(r$spikein_status), NA, as.character(r$spikein_status)),
    stringsAsFactors = FALSE)
}
sheet <- do.call(rbind, sheet_rows); rownames(sheet) <- NULL
sheet$bamReads   <- as.character(sheet$bamReads)
sheet$bamControl <- as.character(sheet$bamControl)
sheet$Peaks      <- as.character(sheet$Peaks)

# =============================================================================
# STAGE 2 DIAGNOSTICS  (build_sample_sheet.R style)
# =============================================================================
# spike-in availability summary (samples lacking a usable count)
bad <- sheet[!is.na(sheet$spikein_status) & !sheet$spikein_status %in% c("OK", "NA"),
             c("SampleID", "Factor", "Treatment", "Condition", "spikein_status")]
if (nrow(bad) > 0) {
  message("\n[NOTE] ", nrow(bad), " sample(s) lack a usable spike-in count; --spikein sheet will ",
          "error only if one of them is in the contrast you run:")
  print(bad, row.names = FALSE)
}

if (length(unique(sheet$Batch)) == 1) message("[INFO] Single batch — no batch term will be used.")

# identity collision: two folders that map to the SAME biological identity
id_key <- with(sheet, paste(Condition, Treatment, Factor, Replicate, sep = "|"))
dup_keys <- id_key[duplicated(id_key)]
if (length(dup_keys) > 0) {
  clash <- sheet[id_key %in% dup_keys, c("SampleID", "Condition", "Treatment", "Factor", "Replicate")]
  clash <- clash[order(paste(clash$Condition, clash$Treatment, clash$Factor, clash$Replicate)), ]
  message("\n[ERROR] Multiple folders map to the SAME biological identity ",
          "(Condition/Treatment/Factor/Replicate). Verify they are not the same sample under two ",
          "folder names, or correct the Replicate numbers:")
  print(clash, row.names = FALSE)
  stop("Identity collision in sample sheet; resolve before running analyses.")
}

# path validation (existing, non-blank paths must point at real files)
chk <- function(v, l) { v <- v[!is.na(v) & v != ""]; m <- v[!file.exists(v)]
  if (length(m)) stop("Missing ", l, " file(s):\n  ", paste(m, collapse = "\n  ")) }
chk(sheet$bamReads, "bamReads"); chk(sheet$bamControl, "bamControl"); chk(sheet$Peaks, "Peaks")

# Absolutize the file paths we just validated, so the sheet is portable across
# working directories regardless of input mode.
abspath <- function(v) { v <- as.character(v); ok <- !is.na(v) & v != "" & file.exists(v)
  v[ok] <- normalizePath(v[ok]); v }
sheet$bamReads   <- abspath(sheet$bamReads)
sheet$bamControl <- abspath(sheet$bamControl)
sheet$Peaks      <- abspath(sheet$Peaks)

# =============================================================================
# SINGLE OUTPUT — the sample sheet
# =============================================================================
sheet <- sheet[, c("SampleID", "Condition", "Treatment", "Factor", "Replicate", "Batch", "IgG",
                   "bamReads", "ControlID", "bamControl", "Peaks", "PeakCaller",
                   "spikein_reads", "spikein_pct", "spikein_status")]
write.csv(sheet, opt$out, row.names = FALSE)
message("\n[INFO] -> ", opt$out, " (", nrow(sheet), " samples). Inspect before analyses.")