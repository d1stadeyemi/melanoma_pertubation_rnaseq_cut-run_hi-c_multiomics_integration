# =============================================================================
# v8 validation harness: stub Bioconductor/edgeR, source the REAL functions,
# exercise the new/changed code paths to a throwaway device, assert behaviour.
# Real Bioconductor calls (windowCounts/getPESizes/estimateDisp/...) are
# cluster-first-run; here we verify the R logic that surrounds them.
# =============================================================================
PASS <- 0L; FAIL <- 0L
ok  <- function(cond, msg) { if (isTRUE(cond)) { PASS<<-PASS+1L; cat(sprintf("  PASS  %s\n", msg)) }
                             else { FAIL<<-FAIL+1L; cat(sprintf("  FAIL  %s\n", msg)) } }
err <- function(expr, msg) { e <- tryCatch({expr; NULL}, error=function(e) conditionMessage(e))
                             if (is.null(e)) { PASS<<-PASS+1L; cat(sprintf("  PASS  %s\n", msg)) }
                             else { FAIL<<-FAIL+1L; cat(sprintf("  FAIL  %s  [%s]\n", msg, e)) } }

# ---- stubs for symbols used INSIDE the functions we test --------------------
cpm <- function(y, log=FALSE, ...) {
  m <- if (!is.null(y$counts)) y$counts else matrix(0,1,1)
  lib <- y$samples$lib.size * y$samples$norm.factors
  v <- sweep(m, 2, lib, "/") * 1e6
  if (log) log2(v + 1) else v
}
aveLogCPM <- function(counts, lib.size=NULL, ...) {
  if (is.matrix(counts)) rowMeans(log2(counts+1)) else log2(as.numeric(counts)+1)
}
estimateDisp <- function(y, design, ...) {
  n <- nrow(y$counts)
  y$common.dispersion  <- 0.08
  y$tagwise.dispersion <- runif(n, 0.04, 0.13)
  y$trended.dispersion <- rep(0.07, n)
  y$AveLogCPM          <- rowMeans(log2(y$counts+1))
  y
}
glmQLFit  <- function(y, design, robust=TRUE, ...) list(df.prior=rep(20,nrow(y$counts)))
glmQLFTest<- function(fit, coef, ...) {
  n <- length(fit$df.prior)
  list(table=data.frame(logFC=rnorm(n,-0.3,1), logCPM=runif(n,-1,5), PValue=runif(n)))
}
# region merging: ids map each window to a region; make 2 windows per region
.make_ids <- function(nwin) rep(seq_len(ceiling(nwin/2)), each=2)[seq_len(nwin)]
mergeWindows <- function(gr, tol, max.width, ...) {
  nwin <- length(gr); ids <- .make_ids(nwin); nreg <- max(ids)
  list(ids=ids, regions=.fakeGR(nreg))
}
mergeWindowsList <- function(grl, tol, max.width, ...) {
  nwin <- sum(vapply(grl, length, integer(1))); ids <- .make_ids(nwin); nreg <- max(ids)
  list(ids=ids, regions=.fakeGR(nreg), weights=rep(1, nwin))
}
combineTests <- function(ids, tab, weights=NULL, ...) {
  nreg <- length(unique(ids))
  data.frame(PValue=runif(nreg), FDR=runif(nreg),
             direction=sample(c("up","down","mixed"),nreg,replace=TRUE))
}
# minimal fake GRanges with seqnames/start/end/length/overlapsAny
.fakeGR <- function(n) structure(list(.chr=rep("chr1",n), .start=seq_len(n)*1000L,
                                      .end=seq_len(n)*1000L+999L, .n=n), class="fakeGR")
seqnames    <- function(x) if (inherits(x,"fakeGR")) x$.chr else base::attr(x,"sn")
start       <- function(x) if (inherits(x,"fakeGR")) x$.start else NextMethod()
end         <- function(x) if (inherits(x,"fakeGR")) x$.end   else NextMethod()
length.fakeGR <- function(x) x$.n
overlapsAny <- function(a, b, ...) rep(c(TRUE,FALSE), length.out=length(a))
GRanges <- function(...) .fakeGR(1L)
IRanges <- function(...) NULL
# csaw QC stubs
reform <- function(param, ...) param
getPESizes <- function(bam, param, ...) {
  list(sizes=as.integer(abs(rnorm(5000, 250, 120))),
       diagnostics=c(total.reads=1e6, mapped.reads=9.5e5, single=2e4,
                     mate.unmapped=1e4, unoriented=5e3, inter.chr=8e3))
}
windowCounts <- function(bams, bin=TRUE, width=1e4, param=NULL, ...) {
  k <- length(bams); structure(list(totals=round(runif(k,8e5,1.2e6))), class="fakeSE")
}
say <- function(...) invisible(NULL)   # silence logging during tests

# Self-contained: read csaw.R, strip the library() block (Bioconductor is not
# needed for these logic tests), and source the rest. The sys.nframe() guard in
# csaw.R already prevents main() from firing on source().
.src <- readLines("csaw.R")
.keep <- TRUE; .out <- character(0)
for (.ln in .src) {
  if (grepl("^suppressWarnings\\(suppressMessages\\(\\{", .ln)) { .keep <- FALSE; next }
  if (!.keep) { if (grepl("^\\}\\)\\)", .ln)) .keep <- TRUE; next }
  .out <- c(.out, .ln)
}
.tmp <- tempfile(fileext = ".R"); writeLines(.out, .tmp)
source(.tmp)

cat("\n== A. dir_labels normalization-awareness ==\n")
optC <- list(test="contrast", cond_ref="Parental", cond_alt="B10")
ok(dir_labels(optC, "spikein")$up == "resistant", "spikein -> 'resistant'")
ok(dir_labels(optC, "sf")$up      == "resistant", "sf -> 'resistant'")
ok(dir_labels(optC, "background")$up == "gained", "background -> 'gained'")
ok(dir_labels(optC, "libsize")$up    == "gained", "libsize -> 'gained'")
ok(dir_labels(optC, NULL)$up          == "gained", "NULL norm -> 'gained' (safe default)")
ok(dir_labels(optC, "spikein")$down  == "lost",   "down stays 'lost'")
optI <- list(test="interaction", cond_ref="Parental", cond_alt="B10")
ok(grepl("stronger loss", dir_labels(optI, "spikein")$up), "interaction labels unaffected by norm")

cat("\n== B. region_sample_lcpm aggregation alignment (the new plumbing) ==\n")
# Mirror the exact aggregation the script does, assert region-level means align.
set.seed(1)
ids <- c(1,1,2,2,3,3); samp.all <- matrix(rnorm(6*4), nrow=6,
         dimnames=list(NULL, c("DMSO_R1","DMSO_R2","UNC_R1","UNC_R2")))
agg <- function(v) as.numeric(tapply(v, ids, mean))
rsl <- do.call(cbind, lapply(seq_len(ncol(samp.all)),
        function(j) as.numeric(tapply(samp.all[,j], ids, mean))))
colnames(rsl) <- colnames(samp.all)
ok(nrow(rsl)==length(unique(ids)), "region_sample_lcpm has one row per region")
ok(all(colnames(rsl)==colnames(samp.all)), "sample columns preserved in order")
# region 1 = mean of windows 1,2 for each sample
ok(isTRUE(all.equal(rsl[1,], colMeans(samp.all[1:2,]), check.attributes=FALSE)),
   "region-1 means equal mean of its two windows (per sample)")
# alignment with a logFC-style agg on the same ids
lfc <- agg(rowMeans(samp.all))
ok(length(lfc)==nrow(rsl), "row-aligned with agg() over same ids")

cat("\n== C. correlation_scatter_fig runs (condition + replicate) ==\n")
tmp <- tempfile("scatter_"); dir.create(tmp)
ncl <- 200
df <- data.frame(FDR=runif(ncl), logFC=rnorm(ncl,-0.4,1.2))
RSL <- matrix(rnorm(ncl*4, 3, 1.5), nrow=ncl,
              dimnames=list(NULL, c("Parental-DMSO-R1","Parental-DMSO-R2","Parental-UNC-R1","Parental-UNC-R2")))
meta <- data.frame(SampleID=colnames(RSL),
                   .role=c("group2","group2","group1","group1"), stringsAsFactors=FALSE)
optS <- list(test="contrast", fdr=0.05, lfc=1.0, group1="Parental|UNC", group2="Parental|DMSO",
             cond_ref="Parental", cond_alt="B10", reference_shift=NULL)
err(correlation_scatter_fig(df, RSL, meta, optS, "TESTTAG", "spikein", tmp, do_replicate=FALSE),
    "condition-only scatter (spikein, no replicate)")
ok(file.exists(file.path(tmp,"TESTTAG_spikein_condition_scatter.png")), "condition scatter PNG written")
err(correlation_scatter_fig(df, RSL, meta, optS, "TESTTAG", "background", tmp, do_replicate=TRUE),
    "condition + replicate scatter (background)")
ok(file.exists(file.path(tmp,"TESTTAG_background_replicate_scatter.png")), "replicate scatter PNG written")
# with reference_shift drawn
optS2 <- modifyList(optS, list(reference_shift=-1.2))
err(correlation_scatter_fig(df, RSL, meta, optS2, "TESTTAG", "spikein", tmp, do_replicate=FALSE),
    "condition scatter with --reference_shift IF line")
# interaction must early-return (no file)
optSi <- modifyList(optS, list(test="interaction"))
err(correlation_scatter_fig(df, RSL, meta, optSi, "TESTTAG", "spikein", tmp, do_replicate=TRUE),
    "interaction early-returns without error")

cat("\n== D. bcv_plot custom render (cloud + lines) ==\n")
ytest <- list(common.dispersion=0.09, tagwise.dispersion=runif(500,0.05,0.14),
              trended.dispersion=rep(0.08,500), AveLogCPM=runif(500,-1,5))
err(bcv_plot(ytest, "TESTTAG", "spikein", tmp), "bcv_plot with full dispersion fields")
ok(file.exists(file.path(tmp,"TESTTAG_spikein_BCV.png")), "BCV PNG written")
# degenerate: collapsed tagwise (sf-like) and missing trend
ysf <- list(common.dispersion=0.16, tagwise.dispersion=rep(0.16,300), AveLogCPM=runif(300,-1,5))
err(bcv_plot(ysf, "TESTTAG", "sf", tmp), "bcv_plot sf-like (collapsed cloud, no trend)")
ok(is.null(bcv_plot(list(), "TESTTAG", "none", tmp)), "bcv_plot no-op when dispersion absent")

cat("\n== E. fragment_qc draws (stubbed getPESizes/windowCounts/reform) ==\n")
metaF <- data.frame(SampleID=c("Par-DMSO-R1","Par-DMSO-R2","Par-UNC-R1","Par-UNC-R2"),
                    bamReads=paste0("/x/",1:4,".bam"), stringsAsFactors=FALSE)
binsF <- structure(list(totals=c(9e5,8.5e5,7e5,7.2e5)), class="fakeSE")
optF <- list(frag_max=700L, bg_bin=10000L)
err(fragment_qc(metaF, binsF, param=NULL, opt=optF, tag="TESTTAG", outdir=tmp),
    "fragment_qc end-to-end (sizes + read-fate + duplication)")
ok(file.exists(file.path(tmp,"TESTTAG_fragment_sizes.png")), "fragment_sizes PNG written")
ok(file.exists(file.path(tmp,"TESTTAG_read_qc.png")),        "read_qc PNG written")
ok(file.exists(file.path(tmp,"TESTTAG_fragment_qc.csv")),    "fragment_qc CSV written")
qc <- read.csv(file.path(tmp,"TESTTAG_fragment_qc.csv"))
ok(all(c("SampleID","total_reads","used_pairs","dup_rate") %in% names(qc)), "QC CSV has expected columns")
ok(nrow(qc)==4, "QC CSV one row per library")

cat("\n== F. mds_plot stagger + title (stub plotMDS via cpm) ==\n")
plotMDS <- function(x, plot=FALSE, ...) {
  list(x=c(-0.95,-0.88,0.85,0.94), y=c(-0.02,0.01,0.96,-0.91),
       var.explained=c(0.50,0.28,0.12,0.10))
}
ymds <- list(counts=matrix(rpois(40,50),nrow=10),
             samples=data.frame(lib.size=rep(1e6,4), norm.factors=rep(1,4),
               cond=rep("parental",4), treat=c("dmso","dmso","unc","unc"),
               batch=rep("b1",4), stringsAsFactors=FALSE))
colnames(ymds$counts) <- c("Parental-DMSO-Migrated-R1","Parental-DMSO-Migrated-R2",
                           "Parental-UNC-Migrated-R1","Parental-UNC-Migrated-R2")
err(mds_plot(ymds, "TESTTAG", "background", tmp, "treat"), "mds_plot runs with overlapping DMSO points")
ok(file.exists(file.path(tmp,"TESTTAG_background_MDS.png")), "MDS PNG written")

cat("\n== G. filter + normfactor diagnostic figures still render ==\n")
d <- list(size=2000L, kind="genomic-background (filterWindowsGlobal)", bg_bin=10000L,
          bin_ab=rnorm(5000,1,1.5), bin_median=0.39,
          win_ab=rnorm(8000,1,1.5), win_filter=rnorm(8000,1.2,1),
          back_level=-1.84, thr=1.585, n_total=8000L, n_keep=5600L,
          n_keep_global=5400L, n_rescue=200L)
err(filter_diagnostics_fig(d, "TESTTAG", tmp, list()), "filter_diagnostics_fig (two-line title)")
ok(file.exists(file.path(tmp,"TESTTAG_filter_w2000.png")), "filter PNG written")
nf_df <- data.frame(SampleID=metaF$SampleID,
                    nf_libsize=rep(1,4), efflib_libsize=runif(4,8e5,1e6),
                    nf_background=runif(4,0.8,1.2), efflib_background=runif(4,8e5,1e6),
                    nf_spikein=runif(4,0.7,1.3), efflib_spikein=runif(4,8e5,1e6))
gk <- c("p|d","p|d","p|u","p|u")
err(normfactor_diagnostics_fig(nf_df, c("libsize","background","spikein"), gk, "TESTTAG", tmp),
    "normfactor_diagnostics_fig (caption in outer margin)")
ok(file.exists(file.path(tmp,"TESTTAG_normfactor_diagnostics.png")), "normfactor PNG written")

cat(sprintf("\n=========== RESULT: %d PASS, %d FAIL ===========\n", PASS, FAIL))
if (FAIL>0) quit(status=1)