# Smoke-test fit_test_regions_multi by stubbing the Bioconductor calls it makes,
# so the REAL function body executes on synthetic 2-window-size input.

# ---- stubs (minimal shapes the function expects) ----
estimateDisp <- function(y, design) { y$common.dispersion <- 0.15; y }
glmQLFit     <- function(y, design, robust=TRUE) list(df.prior = rep(8, 5), .y = y)
glmQLFTest   <- function(fit, coef) {
  n <- nrow(fit$.y$counts)
  list(table = data.frame(logFC = fit$.y$.logFC, logCPM = rep(0, n), PValue = fit$.y$.pval))
}
aveLogCPM <- function(counts, lib.size=NULL) {
  if (is.matrix(counts)) rowMeans(counts) else rep(0, length(counts))
}
# region clustering stub: 2 windows per cluster, contiguous ids
mergeWindowsList <- function(ranges.list, tol, max.width) {
  total <- sum(vapply(ranges.list, length, integer(1)))
  ids <- rep(seq_len(total %/% 2 + total %% 2), each = 2)[seq_len(total)]
  K <- length(unique(ids))
  list(ids = ids, regions = data.frame(.k = seq_len(K)), weights = rep(1, total))
}
mergeWindows <- function(ranges, tol, max.width) {
  n <- length(ranges); ids <- rep(seq_len(n %/% 2 + n %% 2), each = 2)[seq_len(n)]
  list(ids = ids, regions = data.frame(.k = seq_len(length(unique(ids)))))
}
combineTests <- function(ids, tab, weights=NULL) {
  K <- length(unique(ids))
  data.frame(PValue = rep(0.01, K), FDR = c(0.001, 0.2, rep(0.3, K-2))[seq_len(K)],
             direction = rep("mixed", K))
}
seqnames <- function(x) rep("chr1", nrow(x))
start    <- function(x) seq_len(nrow(x)) * 1000L
end      <- function(x) seq_len(nrow(x)) * 1000L + 999L
overlapsAny <- function(a, b) rep(FALSE, nrow(a))

# fake GRanges = integer vector with a length() (windows). Use simple vectors.
mk_gr <- function(n) structure(seq_len(n), class = "fakegr")
length.fakegr <- function(x) length(unclass(x))

# ---- extract just fit_test_regions_multi from csaw_v7.R and eval it ----
src <- readLines("csaw_v7.R")
s <- grep("^fit_test_regions_multi <- function", src)
# find matching close: next top-level "^[a-zA-Z_.].* <- function" after s, then back up
e <- grep("^balance_block <- function", src) - 1
eval(parse(text = paste(src[s:e], collapse = "\n")))
cat("fit_test_regions_multi sourced. lines", s, "-", e, "\n")

# ---- synthetic two-size input ----
mk_y <- function(n, logFC, pval) {
  y <- list()
  y$counts <- matrix(rpois(n*4, 50), nrow = n)
  y$samples <- data.frame(lib.size = rep(1e6, 4), norm.factors = rep(1, 4),
                          group_role = c("group2","group2","group1","group1"))
  y$.logFC <- logFC; y$.pval <- pval; y
}
set.seed(1)
n1 <- 6; n2 <- 4
y1 <- mk_y(n1, logFC = c(-3,-3, 1, -2,-2, 0.1), pval = rep(0.001, n1))
y2 <- mk_y(n2, logFC = c(-4, 2, -3, 0.2),        pval = rep(0.001, n2))
opt <- list(test="contrast", fdr=0.05, lfc=1.0, n_strata=3, spacing=500)
ab  <- list(design = matrix(0, 4, 2), tested = "x")
win_gr_list <- list(mk_gr(n1), mk_gr(n2))
rescued_list <- list(c(F,F,F,T,T,T), c(F,T,T,F))  # some rescued in each size

out <- fit_test_regions_multi(y_list = list(y1, y2), ab = ab, opt = opt,
        win_gr_list = win_gr_list, peaks_gr = NULL, indep_logcpm_list = NULL,
        sizes = c(1000L, 5000L), sig_idx = 1L, rescued_list = rescued_list)

cat("\n=== fit_test_regions_multi RAN without error ===\n")
cat("regions:", nrow(out$table), "\n")
cat("region_origin values:", paste(unique(out$table$region_origin), collapse=","), "\n")
cat("A/B report:\n"); str(out$ab)
cat("disp summary rows:", nrow(out$disp), "\n")
stopifnot(all(c("n_rescued_windows","frac_rescued","region_origin","n_windows") %in% names(out$table)))
stopifnot(length(out$ab) == 10)
cat("\nSMOKE TEST PASSED\n")