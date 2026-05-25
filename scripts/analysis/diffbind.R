# =============================================================================
# diffbind.R  — PEAK-BASED differential binding for PUNCTATE marks
# ASSIGNED: H3K27Ac (well-defined peaks; broad marks go to run_csaw.R instead).
#
# CONTRASTS are explicit, over Condition|Treatment, so all designs work:
#   migration:    --group1 "B10|none"     --group2 "parental|none"
#   UNC in parent:--group1 "parental|UNC"  --group2 "parental|none"
#   UNC in B10:   --group1 "B10|UNC"       --group2 "B10|none"
#
# PIPELINE: subset sheet to (mark x the two groups) -> consensus peaks
#   -> optional IgG greylist (artifact exclusion ONLY, not a contrast term)
#   -> count reads (summit-recentered, paired-end) -> normalize (libsize always;
#   spike-in if --spikein) -> DESeq2 NB test -> log2FC + FDR -> save + QC plots.
#
# BATCH: if the sheet's Batch column has >1 level ACROSS the two groups, the
#   design becomes ~Batch + Condition (batch as covariate). If batch is perfectly
#   confounded with the contrast, DESeq2/DiffBind will error (rank-deficient) —
#   that is correct behaviour; it cannot be corrected, only acknowledged.
#
# SPIKE-IN (future): --spikein CSV (SampleID,factor). Until then libsize only;
#   library-size can mask GLOBAL shifts, so such results are provisional.
# =============================================================================
suppressWarnings(suppressMessages(library(optparse)))
opt <- parse_args(OptionParser(option_list=list(
  make_option("--sheet", type="character"),
  make_option("--mark", type="character", default="H3K27ac"),
  make_option("--treatment", type="character", default=NULL,
              help="Optional: restrict to one Treatment (none|UNC) before contrast"),
  make_option("--group1", type="character", default="B10|none", help="numerator Condition|Treatment"),
  make_option("--group2", type="character", default="parental|none", help="denominator Condition|Treatment"),
  make_option("--outdir", type="character", default="diffbind_out"),
  make_option("--minOverlap", type="integer", default=2L),
  make_option("--summits", type="integer", default=200L),
  make_option("--fdr", type="double", default=0.05),
  make_option("--lfc", type="double", default=1.0),
  make_option("--use_greylist", action="store_true", default=FALSE),
  make_option("--spikein", type="character", default=NULL)
)))
if (is.null(opt$sheet)) stop("--sheet required.")
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
suppressWarnings(suppressMessages(library(DiffBind)))

split_grp <- function(s){ p<-strsplit(s,"\\|")[[1]]; list(cond=p[1], treat=p[2]) }
g1<-split_grp(opt$group1); g2<-split_grp(opt$group2)

full <- read.csv(opt$sheet, stringsAsFactors=FALSE)
ss <- full[toupper(full$Factor)==toupper(opt$mark), ]
if (!is.null(opt$treatment)) ss <- ss[ss$Treatment==opt$treatment, ]
# keep only rows belonging to the two contrast groups
in_g <- function(r,g) r$Condition==g$cond & r$Treatment==g$treat
ss <- ss[in_g(ss,g1) | in_g(ss,g2), ]
if (nrow(ss)==0) stop("No samples for ",opt$mark," in groups ",opt$group1,"/",opt$group2)
ss$grp <- ifelse(in_g(ss,g1), opt$group1, opt$group2)
message("[INFO] ",opt$mark,": ",nrow(ss)," samples"); print(table(ss$grp))
if (any(table(ss$grp)<2)) warning("A group has <2 replicates — weak dispersion.")

tag <- paste0(opt$mark,"_",gsub("[|]","-",opt$group1),"_vs_",gsub("[|]","-",opt$group2))
tmp <- file.path(opt$outdir, paste0("sheet_",tag,".csv")); write.csv(ss,tmp,row.names=FALSE)

obj <- dba(sampleSheet=tmp, minOverlap=opt$minOverlap)
if (opt$use_greylist) {
  if (all(is.na(ss$bamControl))) warning("No bamControl — skipping greylist.")
  else { message("[INFO] IgG greylist..."); obj<-dba.blacklist(obj, blacklist=FALSE, greylist=TRUE) }
}
message("[INFO] Counting (summits=",opt$summits,")...")
obj <- dba.count(obj, summits=opt$summits, bParallel=TRUE)

# batch covariate if >1 batch spans the contrast
multi_batch <- length(unique(ss$Batch))>1
opt$cond_ref <- split_grp(opt$group2)$cond   # denominator condition = reference level
run_norm <- function(o, nt, sf=NULL){
  message("[INFO] === norm: ",nt," ===")
  if (!is.null(sf)){ ids<-dba.show(o)$ID; f<-sf$factor[match(ids,sf$SampleID)]
    if (any(is.na(f))) stop("spikein missing: ",paste(ids[is.na(f)],collapse=","))
    o<-dba.normalize(o, normalize=f)
  } else o<-dba.normalize(o, normalize=DBA_NORM_LIB)
  if (multi_batch){
    message("[INFO] >1 batch across contrast -> attempting to block on Batch.")
    message("[WARN] Cross-batch DiffBind blocking depends on the Batch column being")
    message("       ingested as metadata; verify the design printed below. If batch is")
    message("       perfectly confounded with the contrast this WILL fail (correctly).")
    o<-tryCatch(
      dba.contrast(o, design="~Batch + Condition", reorderMeta=list(Condition=opt$cond_ref)),
      error=function(e){
        message("[WARN] ~Batch+Condition design failed (",conditionMessage(e),
                "); falling back to unblocked contrast. Treat result as batch-confounded.")
        dba.contrast(o, group1=o$masks[[opt$group1]], group2=o$masks[[opt$group2]],
                     name1=opt$group1, name2=opt$group2) })
  } else {
    o<-dba.contrast(o, group1=o$masks[[opt$group1]], group2=o$masks[[opt$group2]],
                    name1=opt$group1, name2=opt$group2)
  }
  o<-dba.analyze(o, method=DBA_DESEQ2)
  res<-as.data.frame(dba.report(o, th=1, bCounts=TRUE))
  base<-file.path(opt$outdir, paste0(tag,"_",nt))
  write.csv(res, paste0(base,"_all_sites.csv"), row.names=FALSE)
  sig<-res[!is.na(res$FDR)&res$FDR<opt$fdr&abs(res$Fold)>=opt$lfc,]
  write.csv(sig, paste0(base,"_significant.csv"), row.names=FALSE)
  message("[INFO] ",nt,": ",nrow(sig)," sig of ",nrow(res))
  pdf(paste0(base,"_QC.pdf"),7,6)
  for (f in list(function()dba.plotHeatmap(o),function()dba.plotPCA(o,label=DBA_CONDITION),
                 function()dba.plotMA(o),function()dba.plotVolcano(o)))
    tryCatch(f(), error=function(e)message("  plot skipped"))
  dev.off(); invisible(o)
}
run_norm(obj,"libsize")
if (!is.null(opt$spikein)){ sf<-read.csv(opt$spikein,stringsAsFactors=FALSE)
  if (!all(c("SampleID","factor")%in%colnames(sf))) stop("spikein needs SampleID,factor")
  run_norm(obj,"spikein",sf)
} else message("[NOTE] no --spikein; libsize only (global shifts provisional).")
message("[INFO] DiffBind done: ",tag," -> ",opt$outdir)