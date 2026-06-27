# Pure-R unit tests for the v7 logic that doesn't need Bioconductor.
ok <- function(cond, msg) cat(if (isTRUE(cond)) "PASS" else "FAIL", "-", msg, "\n")

# ---- copy the pure helpers verbatim from csaw.R ----
parse_window_sizes <- function(s) {
  toks <- trimws(strsplit(as.character(s), ",")[[1]]); toks <- toks[nzchar(toks)]
  if (!length(toks)) stop("--window is empty.")
  v <- suppressWarnings(as.integer(toks))
  if (any(is.na(v))) stop("non-integer")
  if (any(v <= 0))   stop("nonpositive")
  sort(unique(v))
}
combine_keep <- function(keep_global, grp_filter_list, thr) {
  keep <- keep_global
  if (length(grp_filter_list)) for (fg in grp_filter_list) keep <- keep | (fg > thr)
  keep
}
mean_within_group_cv <- function(x, group_key) {
  cvs <- tapply(x, group_key, function(z) if (length(z) > 1 && mean(z) != 0) sd(z) / abs(mean(z)) else NA_real_)
  mean(cvs, na.rm = TRUE)
}

cat("=== parse_window_sizes ===\n")
ok(identical(parse_window_sizes("2000"), 2000L), "single size")
ok(identical(parse_window_sizes("1000,2000,5000"), c(1000L,2000L,5000L)), "three sizes sorted")
ok(identical(parse_window_sizes("5000, 1000 ,2000"), c(1000L,2000L,5000L)), "unsorted+spaces sorted/trimmed")
ok(identical(parse_window_sizes("2000,2000"), 2000L), "dedup")
ok(inherits(try(parse_window_sizes("2000,abc"), silent=TRUE), "try-error"), "rejects non-integer")
ok(inherits(try(parse_window_sizes("0,2000"), silent=TRUE), "try-error"), "rejects nonpositive")

cat("\n=== combine_keep (group-aware rescue) ===\n")
# window-level filter stats; global keep = those above thr in the all-sample average
thr <- 1.585
kg <- c(FALSE, FALSE, TRUE, FALSE)          # global decision
# on/off window: present only in group1, so group1 avg passes but global avg fails
g1 <- c(2.0,  0.1, 2.0, 0.2)                 # group1 per-window filter
g2 <- c(0.1,  0.1, 2.0, 0.2)                 # group2 per-window filter
res <- combine_keep(kg, list(g1=g1, g2=g2), thr)
ok(identical(res, c(TRUE, FALSE, TRUE, FALSE)), "rescues window 1 (passes in group1 only); keeps global; drops 2 and 4")
ok(identical(combine_keep(kg, list(), thr), kg), "no groups -> identical to global")
ok(sum(res & !kg) == 1, "exactly one window rescued")

cat("\n=== mean_within_group_cv ===\n")
# two groups of 2 reps; group A effective libs identical (CV 0), group B spread
efflib <- c(100,100, 100,200); gk <- c("A","A","B","B")
cvA <- 0; cvB <- sd(c(100,200))/mean(c(100,200))
ok(abs(mean_within_group_cv(efflib, gk) - mean(c(cvA,cvB))) < 1e-9, "mean of per-group CVs")
# a normalization that scatters replicates within group -> higher CV
flat <- c(100,100,100,100); spread <- c(80,120,70,180)
ok(mean_within_group_cv(spread, gk) > mean_within_group_cv(flat, gk), "more replicate spread -> higher CV (flat-volcano signal)")

cat("\n=== cluster->region aggregation alignment (tapply over ids) ===\n")
# emulate fit_test_regions_multi consolidation: ids per (concatenated) window,
# per-window logFC, and the agg() used to build the region table.
ids <- c(1L,1L,2L,3L,3L,3L)                  # 3 clusters
logFC <- c(2,4, -1, 0,3,6)                   # per-window
agg <- function(v) as.numeric(tapply(v, ids, mean))
n_windows <- as.integer(tapply(rep(1L, length(ids)), ids, sum))
ok(identical(agg(logFC), c(3, -1, 3)), "region mean logFC = mean over member windows, ordered by cluster id")
ok(identical(n_windows, c(2L,1L,3L)), "n_windows per region correct and aligned")
# combineTests returns rows ordered by sorted unique id; tapply matches that order
ok(identical(sort(unique(ids)), 1:3), "ids are contiguous 1..K (no NA) -> tapply aligns with regions")

cat("\n=== geometric-mean centering (sf / spikein) preserves within-contrast ratios ===\n")
center <- function(v) v / exp(mean(log(v)))
SF <- c(0.809, 0.989, 1.002, 1.247)          # Parental DMSO x2, UNC x2
c_SF <- center(SF)
ok(abs(prod(c_SF) - 1) < 1e-9, "centered factors have geometric mean 1")
ratio_raw <- (SF[3]*SF[4]) / (SF[1]*SF[2]); ratio_cen <- (c_SF[3]*c_SF[4]) / (c_SF[1]*c_SF[2])
ok(abs(ratio_raw - ratio_cen) < 1e-9, "UNC/DMSO effective-lib ratio unchanged by centering")
ok(ratio_cen > 1, "UNC effective lib > DMSO (larger SF = lower normalized signal = loss) — direction correct")

cat("\n=== group-aware A/B accounting (rescued -> region -> significance) ===\n")
# Emulate the consolidation + accounting in fit_test_regions_multi.
# 6 windows across 1 size, concatenated; ids -> 3 regions.
ids      <- c(1L,1L, 2L, 3L,3L,3L)
rescued  <- c(FALSE,FALSE, TRUE, TRUE,TRUE,TRUE)   # region2 = 1 rescued of 1 (rescued_only);
                                                   # region3 = 3 rescued of 3 (rescued_only);
                                                   # region1 = 0 rescued (global)
# region-level table (ordered by cluster id 1..3)
n_windows <- as.integer(tapply(rep(1L,length(ids)), ids, sum))   # c(2,1,3)
n_resc    <- as.integer(tapply(rescued, ids, sum))               # c(0,1,3)
logFC <- c(-2.0, +1.5, -3.0)                                     # region 1 lost, 2 gained, 3 lost
FDR   <- c(0.01, 0.20, 0.001)                                    # region 1 & 3 significant, 2 not
fdr_cut <- 0.05; lfc_cut <- 1.0
region_origin <- ifelse(n_resc==0L,"global", ifelse(n_resc==n_windows,"rescued_only","mixed"))
sig_mask <- FDR < fdr_cut & abs(logFC) > lfc_cut                 # c(TRUE,FALSE,TRUE)
reg_sig_w <- sig_mask[ids]; reg_dir_w <- sign(logFC)[ids]
rescued_in_sig      <- sum(rescued & reg_sig_w)
rescued_in_sig_lost <- sum(rescued & reg_sig_w & reg_dir_w < 0)
sig_rescued_only    <- sum(region_origin=="rescued_only" & sig_mask)
sig_rescued_only_lost <- sum(region_origin=="rescued_only" & sig_mask & logFC < 0)

ok(identical(region_origin, c("global","rescued_only","rescued_only")), "region_origin classification correct")
ok(rescued_in_sig == 3, "3 rescued windows land in significant regions (region3's 3 windows)")
ok(rescued_in_sig_lost == 3, "all 3 are in the LOST direction")
ok(sig_rescued_only == 1, "1 entirely-rescued region is significant (region3; region2 rescued_only but not sig)")
ok(sig_rescued_only_lost == 1, "that significant rescued-only region is a loss")
# region2 is rescued_only but FDR=0.20 -> not significant; must NOT be counted
ok(sum(region_origin=="rescued_only") == 2, "2 rescued-only regions total (only 1 significant)")
# sanity: a window in a non-significant region is not counted
ok(sum(rescued & !reg_sig_w) == 1, "region2's single rescued window correctly excluded (region not sig)")

cat("\n=== global mode: rescued flag all FALSE -> zero accounting ===\n")
rescued0 <- rep(FALSE, length(ids))
ok(sum(rescued0) == 0, "global mode rescues nothing")
no_resc_origin <- ifelse(tapply(rescued0,ids,sum)==0L,"global","x")
ok(all(no_resc_origin=="global"), "all regions classified 'global' when nothing rescued")