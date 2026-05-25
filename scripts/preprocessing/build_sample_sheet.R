#!/usr/bin/env Rscript
# =============================================================================
# build_sample_sheet.R  (STAGE 0)  — turn a manifest (or a raw_dir) into the
# analysis sample sheet that run_diffbind.R and run_csaw.R consume.
#
# TWO INPUT MODES:
#   (A) --manifest manifest.csv   [PREFERRED for your data]
#       Reads the auto-generated, hand-proofread manifest. Uses its real
#       bam_basename/bam_path/peak_path, biological labels, and your manually
#       filled Batch & IgG. This is the safe path because your folder names and
#       inner file names disagree — the manifest records the verified truth.
#   (B) --root raw_dir [--recursive]   [fallback]
#       Parses folder names directly (assumes folder==filename). Fine for clean
#       samples, but will MISS samples whose inner BAM name differs. Prefer (A).
#
# CASE-INSENSITIVE: 'Parental'=='parental', 'H3K27Me3'=='H3K27me3'. You may fill
#   the manifest with any casing; everything is normalized here.
#
# BATCH & IgG are OPTIONAL and may be entirely empty or partially filled:
#   - empty Batch  -> all samples treated as one batch ('batch1'); no batch term.
#   - filled Batch -> downstream adds a batch covariate only when a contrast
#                     actually spans >1 batch.
#   - empty IgG    -> bamControl blank; greylist/background simply not used.
#   - filled IgG   -> resolved to that IgG folder's BAM as bamControl.
#
# OUTPUT (DiffBind-compatible + extras the analysis reads):
#   SampleID, Condition, Treatment, Factor, Replicate, Batch,
#   bamReads, ControlID, bamControl, Peaks, PeakCaller
# =============================================================================

suppressWarnings(suppressMessages(library(optparse)))

opt <- parse_args(OptionParser(option_list=list(
  make_option("--manifest", type="character", default=NULL,
              help="Proofread manifest.csv from generate_manifest.R (preferred)"),
  make_option("--root", type="character", default=NULL,
              help="Fallback: raw_dir to parse directly (if no --manifest)"),
  make_option("--recursive", action="store_true", default=FALSE),
  make_option("--out", type="character", default="sample_sheet.csv"),
  make_option("--bam", type="character", default="dedup", help="dedup|dup (fallback mode only)"),
  make_option("--igg_root", type="character", default=NULL,
              help="Dir containing IgG folders named in the manifest's IgG column (default: --root or manifest dir)"),
  make_option("--treatments", type="character", default="UNC,EPZ"),
  make_option("--narrow_marks", type="character", default="H3K27AC"),
  make_option("--sample_regex", type="character",
              default="^(Bottom10|Bottom|Control|Parental|B10)-(?:UNC-|EPZ-)?([A-Za-z0-9]+)-R([0-9]+)$"),
  make_option("--require_igg", action="store_true", default=FALSE)
)))
if (is.null(opt$manifest) && is.null(opt$root))
  stop("Provide --manifest (preferred) or --root (fallback).")

treatments   <- toupper(trimws(strsplit(opt$treatments, ",")[[1]])); treatments<-treatments[treatments!=""]
narrow_marks <- toupper(trimws(strsplit(opt$narrow_marks, ",")[[1]]))

cond_norm <- function(x){ x<-tolower(trimws(x))
  if (x %in% c("control","parental")) "parental"
  else if (x %in% c("bottom","bottom10","b10")) "B10" else x }
canon_mark <- function(m){ m<-toupper(trimws(m))
  core<-sub("^(H[0-9]+K[0-9]+).*$","\\1",m); rest<-tolower(sub("^H[0-9]+K[0-9]+","",m))
  paste0(core,rest) }
blank <- function(x) is.na(x) | trimws(as.character(x))==""

# =============================================================================
# MODE A: manifest
# =============================================================================
if (!is.null(opt$manifest)) {
  man <- read.csv(opt$manifest, stringsAsFactors=FALSE)
  req <- c("folder","bam_basename","Condition","Treatment","Factor","Replicate")
  miss <- setdiff(req, colnames(man))
  if (length(miss)) stop("Manifest missing columns: ", paste(miss, collapse=", "))

  # bam_path/peak_path may be in the manifest; if not, we cannot guess — error.
  if (!"bam_path" %in% colnames(man) || any(blank(man$bam_path)))
    stop("Manifest needs a non-empty 'bam_path' for every row ",
         "(generate_manifest.R fills this). Re-generate the manifest.")
  if (!"peak_path" %in% colnames(man)) man$peak_path <- ""
  if (!"peak_caller" %in% colnames(man)) man$peak_caller <- NA

  igg_root <- if (!is.null(opt$igg_root)) opt$igg_root else dirname(opt$manifest)

  rows <- list()
  for (i in seq_len(nrow(man))) {
    r <- man[i,]
    caller <- if (!is.na(r$peak_caller) && r$peak_caller!="") r$peak_caller
              else if (toupper(canon_mark(r$Factor)) %in% narrow_marks) "narrow" else "broad"

    # Batch: blank -> batch1
    batch <- if (blank(r$Batch)) "batch1" else trimws(as.character(r$Batch))

    # IgG: blank -> no control. Else resolve to that folder's BAM (look it up
    # in raw_dir via the same aligned/<dedup|dup.marked> structure).
    cid <- NA; bc <- NA
    if ("IgG" %in% colnames(man) && !blank(r$IgG)) {
      igg_folder <- trimws(as.character(r$IgG))
      # Try dedup then dup.marked; accept any single BAM inside.
      for (sd in c("dedup","dup.marked")) {
        d <- file.path(igg_root, igg_folder, "aligned", sd)
        if (dir.exists(d)) {
          b <- list.files(d, "\\.bam$", full.names=TRUE)
          b <- b[!grepl("\\.bai$", b)]
          if (length(b)>=1) { cid<-igg_folder; bc<-b[1]; break }
        }
      }
      if (is.na(bc)) {
        m <- paste0("IgG '", igg_folder, "' for ", r$folder, " not found under ", igg_root)
        if (opt$require_igg) stop(m) else warning(m, " — bamControl blank.")
      }
    }

    rows[[length(rows)+1]] <- data.frame(
      SampleID=r$folder, Condition=cond_norm(r$Condition),
      Treatment=tolower(trimws(r$Treatment)), Factor=canon_mark(r$Factor),
      Replicate=as.integer(r$Replicate), Batch=batch,
      bamReads=r$bam_path, ControlID=cid, bamControl=bc,
      Peaks=ifelse(blank(r$peak_path), NA, r$peak_path), PeakCaller=caller,
      stringsAsFactors=FALSE)
  }
  sheet <- do.call(rbind, rows)
  rownames(sheet) <- NULL
  # Ensure character columns stay character (rbind can coerce to factor)
  sheet$bamReads <- as.character(sheet$bamReads)
  sheet$bamControl <- as.character(sheet$bamControl)
  sheet$Peaks <- as.character(sheet$Peaks)

# =============================================================================
# MODE B: direct raw_dir parse (fallback; assumes folder==filename)
# =============================================================================
} else {
  treat_alt <- paste(treatments, collapse="|")
  sample_re <- opt$sample_regex
  if (grepl("\\(\\?:UNC-\\|EPZ-\\)\\?", sample_re))
    sample_re <- sub("\\(\\?:UNC-\\|EPZ-\\)\\?",
                     paste0("(?:(?:",treat_alt,")-)?"), sample_re)
  bam_sub <- if (opt$bam=="dup") "dup.marked" else "dedup"
  treat_of <- function(f){ h<-treatments[vapply(treatments,
                function(tk) grepl(paste0("-",tk,"-"),toupper(f)),logical(1))]
                if (length(h)>=1) h[1] else "none" }

  if (opt$recursive) { ap<-list.dirs(opt$root,full.names=TRUE,recursive=TRUE)
    ad<-basename(ap); names(ad)<-ap
  } else { d<-list.dirs(opt$root,full.names=FALSE,recursive=FALSE); d<-d[d!=""]
    ad<-d; names(ad)<-file.path(opt$root,d) }
  tp<-names(ad)[grepl(sample_re,ad,perl=TRUE)]; tn<-ad[grepl(sample_re,ad,perl=TRUE)]
  if (length(tp)==0) stop("No folders matched in fallback mode.")

  rows<-list()
  for (i in seq_along(tp)) {
    fp<-tp[i]; folder<-tn[i]
    cond<-sub(sample_re,"\\1",folder,perl=TRUE); mark<-sub(sample_re,"\\2",folder,perl=TRUE)
    rep <-sub(sample_re,"\\3",folder,perl=TRUE)
    bam <- file.path(fp,"aligned",bam_sub,paste0(folder,".bam"))
    isnar<- toupper(canon_mark(mark)) %in% narrow_marks
    sub  <- if (isnar) "macs2.narrow" else "macs2.broad"
    patt <- if (isnar) "narrowPeak" else "broadPeak"
    pk   <- file.path(fp,"peakcalling",sub,paste0(folder,"_peaks.",patt))
    rows[[length(rows)+1]] <- data.frame(
      SampleID=folder, Condition=cond_norm(cond), Treatment=treat_of(folder),
      Factor=canon_mark(mark), Replicate=as.integer(rep), Batch="batch1",
      bamReads=bam, ControlID=NA, bamControl=NA,
      Peaks=pk, PeakCaller=if(isnar)"narrow" else "broad", stringsAsFactors=FALSE)
  }
  sheet<-do.call(rbind,rows)
  rownames(sheet) <- NULL
  # Ensure character columns stay character (rbind can coerce to factor)
  sheet$bamReads <- as.character(sheet$bamReads)
  sheet$bamControl <- as.character(sheet$bamControl)
  sheet$Peaks <- as.character(sheet$Peaks)
}

# ---- validate existing paths ------------------------------------------------
chk <- function(v,l){ v<-v[!is.na(v) & v!=""]; m<-v[!file.exists(v)]
  if (length(m)) stop("Missing ",l," file(s):\n  ",paste(m,collapse="\n  ")) }
chk(sheet$bamReads,"bamReads"); chk(sheet$bamControl,"bamControl")
# Peaks may be NA (csaw doesn't require them); only check the non-NA ones.
chk(sheet$Peaks[!is.na(sheet$Peaks)], "Peaks")

message("[INFO] Factor/Treatment x Condition:")
print(table(paste(sheet$Factor,sheet$Treatment,sep="/"), sheet$Condition))
if (length(unique(sheet$Batch))==1) message("[INFO] Single batch — no batch term will be used.")
write.csv(sheet, opt$out, row.names=FALSE)
message("[INFO] -> ", opt$out, " (", nrow(sheet), " samples). Inspect before analyses.")