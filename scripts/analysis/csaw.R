# =============================================================================
# csaw.R  — WINDOW-BASED differential binding for BROAD marks
# ASSIGNED: H3K27me3, H3K9me3, H3K4me1 (diffuse domains; peak calling
#   unreliable, so window-based testing is the stronger choice). Punctate
#   H3K27Ac -> run_diffbind.R.
#
# CONTRASTS explicit over Condition|Treatment (same scheme as run_diffbind.R):
#   migration:    --group1 "B10|none"     --group2 "parental|none"
#   UNC in parent:--group1 "parental|UNC"  --group2 "parental|none"
#   UNC in B10:   --group1 "B10|UNC"       --group2 "B10|none"
#
# PIPELINE: subset sheet -> count reads into genome WINDOWS (no peaks needed)
#   -> background bins + enrichment filter (IgG optional as background ref)
#   -> normalize (TMM-on-background default; spike-in if --spikein)
#   -> edgeR NB test -> merge windows into regions -> region FDR + log2FC -> save.
#
# BATCH: if Batch has >1 level across the contrast, design = ~Batch + group.
#   Perfect confounding -> model is rank-deficient and edgeR will fail (correct).
#
# SPIKE-IN (future): --spikein CSV (SampleID,factor). For broad REPRESSIVE marks
#   especially, TMM/libsize can MASK genome-wide change; such results provisional.
# =============================================================================
suppressWarnings(suppressMessages(library(optparse)))
opt <- parse_args(OptionParser(option_list=list(
  make_option("--sheet", type="character"),
  make_option("--mark", type="character", help="H3K27me3|H3K9me3|H3K4me1"),
  make_option("--treatment", type="character", default=NULL),
  make_option("--group1", type="character", default="B10|none"),
  make_option("--group2", type="character", default="parental|none"),
  make_option("--outdir", type="character", default="csaw_out"),
  make_option("--window", type="integer", default=2000L),
  make_option("--spacing", type="integer", default=500L),
  make_option("--bg_bin", type="integer", default=10000L),
  make_option("--min_lfc_bg", type="double", default=log2(3)),
  make_option("--frag_max", type="integer", default=700L),
  make_option("--fdr", type="double", default=0.05),
  make_option("--lfc", type="double", default=1.0),
  make_option("--blacklist", type="character", default=NULL),
  make_option("--spikein", type="character", default=NULL)
)))
if (is.null(opt$sheet)) stop("--sheet required.")
if (is.null(opt$mark))  stop("--mark required.")
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
suppressWarnings(suppressMessages({library(csaw);library(edgeR)
  library(GenomicRanges);library(rtracklayer)}))

split_grp<-function(s){p<-strsplit(s,"\\|")[[1]];list(cond=p[1],treat=p[2])}
g1<-split_grp(opt$group1); g2<-split_grp(opt$group2)
full<-read.csv(opt$sheet, stringsAsFactors=FALSE)
ss<-full[toupper(full$Factor)==toupper(opt$mark),]
if (!is.null(opt$treatment)) ss<-ss[ss$Treatment==opt$treatment,]
in_g<-function(r,g) r$Condition==g$cond & r$Treatment==g$treat
ss<-ss[in_g(ss,g1)|in_g(ss,g2),]
if (nrow(ss)==0) stop("No samples for ",opt$mark," in ",opt$group1,"/",opt$group2)
ss$grp<-ifelse(in_g(ss,g1),"g1","g2")              # g1=numerator
ss<-ss[order(ss$grp, ss$Replicate),]
message("[INFO] ",opt$mark,": ",nrow(ss)," samples"); print(table(ss$grp))
if (any(table(ss$grp)<2)) warning("A group has <2 replicates — weak dispersion.")

tag<-paste0(opt$mark,"_",gsub("[|]","-",opt$group1),"_vs_",gsub("[|]","-",opt$group2))
bams<-ss$bamReads
group<-factor(ss$grp, levels=c("g2","g1"))          # g2 reference
batch<-factor(ss$Batch); multi_batch<-nlevels(batch)>1

param<-readParam(pe="both", max.frag=opt$frag_max)
if (!is.null(opt$blacklist)){ param<-reform(param, discard=import(opt$blacklist))
  message("[INFO] blacklist applied.") }

message("[INFO] window counts (w=",opt$window,", s=",opt$spacing,")...")
win<-windowCounts(bams, width=opt$window, spacing=opt$spacing, param=param)
colnames(win)<-ss$SampleID
message("[INFO] background bins (",opt$bg_bin,") + filter...")
bins<-windowCounts(bams, bin=TRUE, width=opt$bg_bin, param=param)
keep<-filterWindowsGlobal(win,bins)$filter > opt$min_lfc_bg
win<-win[keep,]; message("[INFO] kept ",sum(keep)," / ",length(keep)," windows.")

if (multi_batch) message("[INFO] >1 batch across contrast -> design ~Batch + group.")
design <- if (multi_batch) model.matrix(~batch+group) else model.matrix(~group)

run_norm<-function(w,b,nt,sf=NULL){
  message("[INFO] === norm: ",nt," ===")
  y<-asDGEList(w)
  if (!is.null(sf)){ f<-sf$factor[match(ss$SampleID,sf$SampleID)]
    if (any(is.na(f))) stop("spikein missing: ",paste(ss$SampleID[is.na(f)],collapse=","))
    y$samples$norm.factors<-f/exp(mean(log(f)))
  } else y$samples$norm.factors<-normFactors(b, se.out=TRUE)$norm.factors
  y<-estimateDisp(y, design)
  fit<-glmQLFit(y, design, robust=TRUE)
  res<-glmQLFTest(fit, coef=ncol(design))           # last coef = group g1 vs g2
  merged<-mergeWindows(rowRanges(w), tol=opt$spacing, max.width=50000)
  tab<-combineTests(merged$ids, res$table)
  out<-cbind(as.data.frame(granges(merged$regions)), as.data.frame(tab))
  base<-file.path(opt$outdir, paste0(tag,"_",nt))
  write.csv(out, paste0(base,"_all_regions.csv"), row.names=FALSE)
  lcol<-intersect(c("rep.logFC","best.logFC"), colnames(out))[1]
  sig<-out[!is.na(out$FDR)&out$FDR<opt$fdr,]
  if (!is.na(lcol)) sig<-sig[abs(sig[[lcol]])>=opt$lfc,]
  write.csv(sig, paste0(base,"_significant.csv"), row.names=FALSE)
  message("[INFO] ",nt,": ",nrow(sig)," sig regions of ",nrow(out))
  pdf(paste0(base,"_QC.pdf"),7,6)
  tryCatch(plotMDS(cpm(y,log=TRUE), labels=ss$SampleID, col=as.integer(group),
                   main=paste(tag,nt,"MDS")), error=function(e)message("  MDS skip"))
  tryCatch({plotSmear(res, main=paste(tag,nt)); abline(h=c(-opt$lfc,opt$lfc),col="blue",lty=2)},
           error=function(e)message("  smear skip"))
  dev.off(); invisible(out)
}
run_norm(win,bins,"libsize")
if (!is.null(opt$spikein)){ sf<-read.csv(opt$spikein,stringsAsFactors=FALSE)
  if (!all(c("SampleID","factor")%in%colnames(sf))) stop("spikein needs SampleID,factor")
  run_norm(win,bins,"spikein",sf)
} else message("[NOTE] no --spikein; libsize only (broad-mark global change provisional).")
message("[INFO] csaw done: ",tag," -> ",opt$outdir)