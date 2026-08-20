#!/usr/bin/env Rscript
###############################################################################
# deseq2_03_model.R   --  STAGE 3 of 4:  PREFILTER, FIT, TRANSFORM, SAMPLE QC
#
# The expensive stage. Builds the DESeqDataSet, fits dispersions and
# coefficients, and produces every sample-level diagnostic. It does NOT extract
# a contrast -- that is stage 4, deliberately, so that changing --alpha or
# --lfc-report re-runs seconds of work instead of the whole fit.
#
# WRITES
#   tables/06_prefilter_sensitivity, 06_prefilter_summary, 07_size_factors,
#          11_pca_coordinates(+_batch_removed), 11_pca_scree, 12_sample_correlation
#   figures/07_size_factor_vs_libsize, 07/08 count + RLE distributions,
#          09_dispersion_estimates, 10_meanSD_*, 11_pca(+_batch_removed),
#          11_pca_scree, 12_sample_distance_heatmap, 12_sample_correlation_heatmap,
#          13_cooks_distance, 14_top_variable_features_heatmap
#   _state/03_state.rds
#
#  ==========================================================================
#  PARAMETERS OWNED BY THIS STAGE
#  ==========================================================================
#  --min-count <INT>     Default 10. Keep features with >= this many counts...
#  --min-samples <INT>   Default = size of the SMALLEST condition group.
#                        ...in at least this many samples. This is the
#                        vignette's recommendation verbatim: a count of 10 is
#                        reasonable for bulk RNA-seq, and the minimal number of
#                        samples should be the smallest group size.
#                        Prefiltering is NOT required for correctness -- it
#                        reduces memory and speeds up fitting, and results()
#                        does independent filtering of its own. The full
#                        sensitivity grid is written so you can see exactly how
#                        many features each choice would keep.
#
#  --test <wald|lrt>     Default wald. LRT compares the full design against
#                        --reduced and tests ALL dropped terms at once; the LFC
#                        column still shows the named contrast, so the p-value
#                        and the fold change answer different questions. Use LRT
#                        for "does condition matter at all" across >2 levels.
#  --reduced <FORMULA>   Default = design minus condition.
#
#  --fit-type <parametric|local|mean|glmGamPoi>   Default parametric.
#                        Change only if 09_dispersion_estimates.pdf shows the
#                        red trend line missing the point cloud.
#
#  --transform <vst|rlog|both>   Default vst. Both are used for QC/clustering
#                        ONLY -- never for testing. rlog is much slower and is
#                        the better choice only for small n.
#                        blind = FALSE throughout: the vignette notes blind
#                        dispersion estimation is inappropriate when many genes
#                        differ by design, because it attributes design
#                        differences to noise and over-shrinks the transformed
#                        values. blind=FALSE still does not use group identity
#                        to move samples around -- only the fitted mean-
#                        dispersion trend.
#
#  ==========================================================================
#  TECHNICAL REPLICATE COLLAPSING  (--collapse-technical, set in stage 1)
#  ==========================================================================
#  Applied here via DESeq2::collapseReplicates. The vignette's distinction is
#  exact: a technical replicate means multiple sequencing runs of the SAME
#  library, and biological replicates must never be collapsed with this
#  function. This stage enforces that every other design variable is constant
#  within a collapse group, then re-checks design rank and replicate counts
#  afterwards -- collapsing can turn a fittable design into an unfittable one.
#
#  If a library was split across lanes L001/L002, the upstream submitter made it
#  two "replicates". That is the case this flag exists for; it is NOT a batch
#  effect and should not be modelled as one.
#
#  ==========================================================================
#  BATCH-CORRECTED PCA IS A PICTURE, NOT A MODEL
#  ==========================================================================
#  11_pca_batch_removed uses limma::removeBatchEffect on the transformed matrix,
#  with design = model.matrix(~condition) passed in so treatment variation is
#  retained. The vignette FAQ prescribes exactly this, and is explicit that vst
#  does not itself remove batch variation. The corrected matrix is used for
#  VISUALISATION ONLY. DESeq() always sees raw counts plus the covariate term.
###############################################################################

.args_raw <- commandArgs(trailingOnly = TRUE)
.self <- (function() {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) normalizePath(sub("^--file=", "", a[1]), mustWork = FALSE) else NA_character_
})()
source(file.path(dirname(.self), "deseq2_common.R"))

opt <- parse_args(.args_raw)
.unknown <- setdiff(names(opt), .known_keys)
if (length(.unknown) > 0) { message("ERROR: unknown option(s): ",
  paste0("--", .unknown, collapse = ", ")); quit(status = 2) }
g <- mk_getters(opt)

outdir <- g$chr("outdir", "deseq2_out")
log_open(outdir, "3/model")
logbanner("STAGE 3: PREFILTER, DESeq() FIT, TRANSFORMS, SAMPLE QC")

# Validate this stage's own flags FIRST. A typo should be reported as a typo,
# not surface later as an upstream-signature mismatch.
min_count   <- g$int("min-count", 10L, min = 0)
min_samples <- g$int("min-samples", NA_integer_, min = 1)
test_type   <- tolower(g$chr("test", "wald"))
reduced_str <- g$chr("reduced")
fit_type    <- g$chr("fit-type", "parametric")
transform   <- g$chr("transform", "vst")

if (!test_type %in% c("wald","lrt"))                            die("--test must be wald or lrt.")
if (!fit_type %in% c("parametric","local","mean","glmGamPoi"))  die("--fit-type invalid: ", fit_type)
if (!transform %in% c("vst","rlog","both"))                     die("--transform invalid: ", transform)

force <- g$bool("force")
st1 <- load_state(outdir, 1, opt, force = force)
st2 <- load_state(outdir, 2, opt, force = force)

coldata        <- st1$coldata
design_formula <- st1$design_formula
covar_vars     <- st1$covar_vars
has_continuous <- st1$has_continuous
control        <- st1$control;  treat <- st1$treat
tech_col       <- st1$tech_col
cts            <- st2$cts
ann            <- st2$ann
collapse_gene  <- st2$collapse_gene
threads        <- st1$threads

rlib <- g$chr("rlib")
if (!is.null(rlib)) { if (!dir.exists(rlib)) die("--rlib not found: ", rlib); .libPaths(c(rlib, .libPaths())) }
if (!requireNamespace("DESeq2", quietly = TRUE))
  die("DESeq2 is not installed in this R environment (activate multiomics-r).")
suppressWarnings(suppressMessages({ library(DESeq2); if (has("ggplot2")) library(ggplot2) }))
set.seed(st1$seed)

reduced_formula <- NULL
if (test_type == "lrt") {
  reduced_formula <- if (!is.null(reduced_str))
    tryCatch(stats::as.formula(reduced_str), error = function(e) die("--reduced invalid: ", reduced_str))
  else stats::update(design_formula, ~ . - condition)
  if ("condition" %in% all.vars(reduced_formula))
    die("--reduced must NOT contain `condition` (it is the term being tested).")
  rmiss <- setdiff(all.vars(reduced_formula), colnames(coldata))
  if (length(rmiss)) die("--reduced term(s) not in the sample table: ", paste(rmiss, collapse=", "))
  logmsg("LRT: full ", paste(deparse(design_formula), collapse=""),
         "  vs reduced ", paste(deparse(reduced_formula), collapse=""))
}

BPPARAM <- if (threads > 1 && has("BiocParallel")) BiocParallel::MulticoreParam(threads) else NULL
if (!is.null(BPPARAM)) logmsg("BiocParallel: ", threads, " workers")

## ===========================================================================
## THE ENGINE -- run once per level (transcript, and gene if --collapse-to-gene)
## ===========================================================================
fit_level <- function(cts, ann, coldata, level_dir, level_lab) {
  io <- make_io(level_dir)
  logmsg("---- level: ", level_lab, " -> ", level_dir, " ----")

  rownames(coldata) <- coldata$sample
  stopifnot(identical(colnames(cts), rownames(coldata)))
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = design_formula)

  ## --- technical replicate collapsing --------------------------------------
  if (!is.null(tech_col)) {
    grp <- as.character(colData(dds)[[tech_col]])
    for (v in unique(c("condition", covar_vars))) {
      bad <- tapply(as.character(colData(dds)[[v]]), grp, function(x) length(unique(x)) > 1)
      if (any(bad))
        die("--collapse-technical '", tech_col, "': variable '", v, "' is not constant within ",
            "group(s): ", paste(names(bad)[bad], collapse=", "),
            ".\n       Technical replicates must share every experimental variable. If they ",
            "do not, these are BIOLOGICAL replicates and must not be collapsed.")
    }
    n_before <- ncol(dds)
    dds <- collapseReplicates(dds, groupby = factor(grp), run = as.character(colData(dds)$sample))
    coldata <- as.data.frame(colData(dds)); coldata$sample <- rownames(coldata)
    logmsg("collapsed technical replicates on '", tech_col, "': ", n_before, " -> ", ncol(dds), " libraries")
    io$table(coldata, "01c_after_technical_collapse.csv")
    ri <- design_rank_report(coldata, design_formula)
    if (!ri$ok) die("after collapsing, the design is no longer full rank (", ri$rank, "/", ri$ncol,
                    "). Drop the covariate or rebalance.")
    for (cc in c(control, treat)) {
      nn <- sum(as.character(coldata$condition) == cc)
      if (nn < 2) die("after collapsing, group '", cc, "' has ", nn,
                      " biological replicate(s); at least 2 are required.")
    }
  }
  cond_int <- as.integer(coldata$condition)

  ## --- prefilter -----------------------------------------------------------
  group_sizes <- table(droplevels(coldata$condition))
  smallest <- min(group_sizes)
  ms <- if (is.na(min_samples)) as.integer(smallest) else min_samples
  if (ms > ncol(dds)) die("--min-samples (", ms, ") exceeds the number of samples (", ncol(dds), ").")

  thr_grid <- expand.grid(min_count = c(1,5,10,20),
                          min_samples = sort(unique(c(1L, as.integer(smallest), as.integer(ms)))))
  thr_grid$features_kept <- apply(thr_grid, 1, function(r)
    sum(rowSums(counts(dds) >= r[["min_count"]]) >= r[["min_samples"]]))
  thr_grid$is_active <- thr_grid$min_count == min_count & thr_grid$min_samples == ms
  io$table(thr_grid, "06_prefilter_sensitivity.csv")

  keep <- rowSums(counts(dds) >= min_count) >= ms
  io$table(data.frame(stage = c("input","kept","removed"),
                      n_features = c(nrow(dds), sum(keep), sum(!keep)),
                      note = c("all features", paste0(">=", min_count, " counts in >=", ms, " samples"),
                               "low-count"), stringsAsFactors = FALSE), "06_prefilter_summary.csv")
  dds <- dds[keep, ]
  logmsg("prefilter (>=", min_count, " in >=", ms, "): kept ", sum(keep), "/", length(keep))
  if (nrow(dds) < 100)
    logmsg("WARNING: only ", nrow(dds), " features survived prefiltering. Check --min-count / ",
           "--min-samples and that the annotation build matches the counts.")

  ## --- fit -----------------------------------------------------------------
  logmsg("running DESeq() (test=", toupper(test_type), ", fitType=", fit_type, ") ...")
  a <- list(object = dds, fitType = fit_type)
  if (test_type == "lrt") { a$test <- "LRT"; a$reduced <- reduced_formula }
  if (!is.null(BPPARAM)) { a$parallel <- TRUE; a$BPPARAM <- BPPARAM }
  dds <- do.call(DESeq, a)

  sf <- sizeFactors(dds); if (is.null(sf)) sf <- rep(NA_real_, ncol(dds))
  io$table(data.frame(sample = colnames(dds), size_factor = as.numeric(sf),
                      library_size = colSums(counts(dds)), stringsAsFactors = FALSE),
           "07_size_factors.csv")
  io$base("07_size_factor_vs_libsize", function() {
    plot(colSums(counts(dds)), sf, xlab="library size (filtered)", ylab="size factor",
         main="Size factor vs library size", pch=19, col=cond_int)
    text(colSums(counts(dds)), sf, labels=colnames(dds), pos=3, cex=.6) })

  norm_cts <- counts(dds, normalized = TRUE)
  io$base("07_count_distribution_normalized", function() { par(mar=c(9,5,3,2))
    boxplot(log2(norm_cts+1), las=2, ylab="log2(norm count + 1)",
            main="Normalized count distribution", col=cond_int) }, width=8, height=5)
  rle_plot <- function(mat, main) {
    logm <- log2(mat+1); med <- apply(logm,1,median)
    par(mar=c(9,5,3,2)); boxplot(logm-med, las=2, ylab="RLE", main=main, col=cond_int, outline=FALSE)
    abline(h=0, lty=2) }
  io$base("08_RLE_raw", function() rle_plot(counts(dds), "RLE (raw counts)"), width=8, height=5)
  io$base("08_RLE_normalized", function() rle_plot(norm_cts, "RLE (normalized)"), width=8, height=5)
  io$base("09_dispersion_estimates", function() plotDispEsts(dds, main="Dispersion estimates"))

  ## --- transforms ----------------------------------------------------------
  do_vst  <- transform %in% c("vst","both"); do_rlog <- transform %in% c("rlog","both")
  vsd <- if (do_vst) tryCatch(vst(dds, blind=FALSE),
                              error=function(e) varianceStabilizingTransformation(dds, blind=FALSE)) else NULL
  rld <- if (do_rlog) rlog(dds, blind=FALSE) else NULL
  ntd <- normTransform(dds)
  qc_obj <- if (!is.null(vsd)) vsd else rld
  qc_lab <- if (!is.null(vsd)) "VST" else "rlog"

  if (has("vsn")) {
    io$base("10_meanSD_log2", function()
      print(vsn::meanSdPlot(assay(ntd), plot=FALSE)$gg + ggplot2::ggtitle("mean-SD: log2(n+1)")))
    if (!is.null(vsd)) io$base("10_meanSD_vst", function()
      print(vsn::meanSdPlot(assay(vsd), plot=FALSE)$gg + ggplot2::ggtitle("mean-SD: VST")))
    if (!is.null(rld)) io$base("10_meanSD_rlog", function()
      print(vsn::meanSdPlot(assay(rld), plot=FALSE)$gg + ggplot2::ggtitle("mean-SD: rlog")))
  }

  ## --- PCA -----------------------------------------------------------------
  factor_covars <- covar_vars[vapply(covar_vars, function(v) is.factor(coldata[[v]]), logical(1))]
  shape_var <- if (length(factor_covars)) factor_covars[1] else NULL
  intgroup  <- if (!is.null(shape_var)) c("condition", shape_var) else "condition"
  pca <- plotPCA(qc_obj, intgroup=intgroup, returnData=TRUE)
  pv  <- round(100 * attr(pca, "percentVar"))
  io$table(pca, "11_pca_coordinates.csv")

  draw_pca <- function(dat, ttl, fname) {
    if (has("ggplot2")) {
      aes_use <- if (!is.null(shape_var))
        ggplot2::aes(PC1, PC2, color=condition, shape=.data[[shape_var]])
      else ggplot2::aes(PC1, PC2, color=condition)
      gg <- ggplot2::ggplot(dat, aes_use) + ggplot2::geom_point(size=4) +
        ggplot2::xlab(paste0("PC1: ", pv[1], "% var")) + ggplot2::ylab(paste0("PC2: ", pv[2], "% var")) +
        ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=0.18)) +
        ggplot2::theme_bw() + ggplot2::ggtitle(ttl)
      # coord_fixed() dropped on purpose: with PC1 ~96% vs PC2 ~3% it squashes the
      # panel into an unreadable strip; the axis titles preserve scale honesty.
      gg <- if (has("ggrepel"))
        gg + ggrepel::geom_text_repel(ggplot2::aes(label=name), size=3, show.legend=FALSE,
                                      max.overlaps=Inf, min.segment.length=0)
      else gg + ggplot2::geom_text(ggplot2::aes(label=name), vjust=-1, size=3, show.legend=FALSE)
      io$gg(fname, gg)
    } else io$base(fname, function() {
      plot(dat$PC1, dat$PC2, col=as.integer(factor(dat$condition)), pch=19,
           xlab=paste0("PC1: ",pv[1],"%"), ylab=paste0("PC2: ",pv[2],"%"), main=ttl)
      text(dat$PC1, dat$PC2, dat$name, pos=3, cex=.6) })
  }
  draw_pca(pca, paste0("PCA (", qc_lab, ")"), "11_pca")

  if (!is.null(shape_var) && has("limma")) {
    tryCatch({
      mm_keep <- stats::model.matrix(~ condition, data = as.data.frame(colData(qc_obj)))
      qc_bc <- qc_obj
      assay(qc_bc) <- limma::removeBatchEffect(assay(qc_obj),
                        batch = colData(qc_obj)[[shape_var]], design = mm_keep)
      pca_bc <- plotPCA(qc_bc, intgroup=intgroup, returnData=TRUE)
      io$table(pca_bc, "11_pca_coordinates_batch_removed.csv")
      draw_pca(pca_bc, paste0("PCA (", qc_lab, ", '", shape_var, "' removed)"), "11_pca_batch_removed")
      logmsg("batch-corrected PCA written (limma::removeBatchEffect on '", shape_var,
             "'). VISUALISATION ONLY -- DESeq() uses raw counts + the covariate term.")
    }, error = function(e) logmsg("WARNING: batch-corrected PCA failed (", conditionMessage(e), ")."))
  }

  rv <- matrixStats::rowVars(assay(qc_obj))
  sel <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
  pf <- prcomp(t(assay(qc_obj)[sel, , drop=FALSE]))
  scree <- (pf$sdev^2 / sum(pf$sdev^2))[seq_len(min(10, length(pf$sdev)))]
  io$table(data.frame(PC=seq_along(scree), pct_variance=round(100*scree,2)), "11_pca_scree.csv")
  io$base("11_pca_scree", function() barplot(100*scree, names.arg=paste0("PC",seq_along(scree)),
    ylab="% variance", main="PCA scree (top-500 variable features)"))

  ## --- distances / correlation / outliers ----------------------------------
  sampleDists <- dist(t(assay(qc_obj))); distMat <- as.matrix(sampleDists)
  rownames(distMat) <- colnames(distMat) <- paste(colnames(qc_obj), coldata$condition, sep=" | ")
  if (has("pheatmap")) {
    cols <- if (has("RColorBrewer")) colorRampPalette(rev(RColorBrewer::brewer.pal(9,"Blues")))(255) else NULL
    io$base("12_sample_distance_heatmap", function() pheatmap::pheatmap(distMat,
      clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, color=cols,
      display_numbers=TRUE, number_format="%.0f", number_color="grey15", fontsize_number=9,
      main=paste0("Sample-to-sample Euclidean distance (", qc_lab, ")")))
  } else io$base("12_sample_distance_heatmap", function() heatmap(distMat, main="Sample distances"))

  corMat <- cor(assay(qc_obj), method="spearman")
  io$table(data.frame(sample=rownames(corMat), as.data.frame(corMat), check.names=FALSE),
           "12_sample_correlation.csv")
  if (has("pheatmap")) io$base("12_sample_correlation_heatmap", function()
    pheatmap::pheatmap(corMat, display_numbers=TRUE, number_format="%.2f", number_color="grey15",
      fontsize_number=9, main=paste0("Sample Spearman correlation (", qc_lab, ")")))

  cooks_title <- if (has_continuous)
    "Cook's distance per sample\n(MANUAL INSPECTION ONLY: continuous covariate in design)" else
    "Cook's distance per sample"
  io$base("13_cooks_distance", function() { par(mar=c(9,5,3,2))
    boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2, main=cooks_title,
            ylab="log10(Cook's)", col=cond_int) }, width=8, height=5)
  if (has_continuous)
    logmsg("NOTE: 13_cooks_distance is for MANUAL inspection. With a continuous variable in ",
           "the design DESeq2 does not automatically flag or replace outliers, so nothing ",
           "in this run acted on these values.")

  topVar <- head(order(matrixStats::rowVars(assay(qc_obj)), decreasing=TRUE), g$int("top-n", 30L, min=1))
  if (has("pheatmap") && length(topVar) >= 2) {
    anndf <- data.frame(condition=coldata$condition); rownames(anndf) <- colnames(qc_obj)
    mv <- assay(qc_obj)[topVar, , drop=FALSE]; mv <- mv - rowMeans(mv)
    io$base("14_top_variable_features_heatmap", function() pheatmap::pheatmap(mv,
      annotation_col=anndf, show_rownames=(length(topVar) <= 60),
      main=paste0("Top ", length(topVar), " variable features (", qc_lab, ", centered)")),
      width=8, height=9)
  }

  list(dds = dds, vsd = vsd, rld = rld, qc_obj = qc_obj, qc_lab = qc_lab,
       coldata = coldata, ann = ann, norm_cts = norm_cts, min_samples_used = ms,
       level_dir = level_dir, level_lab = level_lab)
}

## ===========================================================================
## DISPATCH
## ===========================================================================
base_dir <- if (collapse_gene) file.path(outdir, "transcript_level") else outdir
levels_out <- list()
levels_out[["primary"]] <- fit_level(cts, ann, coldata, base_dir,
                                     if (collapse_gene) "transcript" else "primary")

if (collapse_gene) {
  if (is.null(ann)) {
    logmsg("WARNING: --collapse-to-gene needs annotation; skipping the gene level.")
  } else {
    # Sum transcript-model counts per gene symbol. This consolidates isoform-split
    # signal but CANNOT recover reads featureCounts already discarded as ambiguous
    # -- only re-counting against a gene-level annotation does that.
    key <- ann$gene_name[match(rownames(cts), ann$gene_id)]
    key[is.na(key) | !nzchar(key)] <- rownames(cts)[is.na(key) | !nzchar(key)]
    cts_gene <- rowsum(cts, group = key); storage.mode(cts_gene) <- "integer"
    modal <- function(x) { x <- x[!is.na(x) & nzchar(x)]
      if (!length(x)) NA_character_ else names(sort(table(x), decreasing = TRUE))[1] }
    bt <- vapply(split(ann$biotype[match(rownames(cts), ann$gene_id)], key), modal, character(1))
    ch <- vapply(split(ann$chrom[match(rownames(cts), ann$gene_id)],  key), modal, character(1))
    ann_gene <- data.frame(gene_id = rownames(cts_gene), gene_name = rownames(cts_gene),
                           biotype = unname(bt[rownames(cts_gene)]),
                           chrom   = unname(ch[rownames(cts_gene)]), stringsAsFactors = FALSE)
    make_io(outdir)$table(data.frame(level=c("transcript","gene"),
                                     n_features=c(nrow(cts), nrow(cts_gene))),
                          "gene_collapse_report.csv")
    logmsg("gene-level collapse: ", nrow(cts), " transcript models -> ", nrow(cts_gene), " genes")
    levels_out[["gene"]] <- fit_level(cts_gene, ann_gene, coldata,
                                      file.path(outdir, "gene_level"), "gene")
  }
}

save_state(list(levels = levels_out, test_type = test_type, reduced_formula = reduced_formula,
                fit_type = fit_type, transform = transform, min_count = min_count),
           outdir, 3, opt)
logmsg("STAGE 3 COMPLETE  (", length(levels_out), " level(s) fitted)")
log_close()