#!/usr/bin/env Rscript
# =============================================================================
# csaw_vs_cpm_agreement.R
#   Validate that the DESCRIPTIVE CPM log2-ratio features used for ML actually
#   track the INFERENTIAL, FDR-controlled differential calls from csaw.
#
# WHY THIS EXISTS (the non-redundant bridge)
#   Two different statistical objects describe the same biology:
#     (1) csaw rep.logFC  — edgeR negative-binomial, replicate-aware, SHRUNK
#                           (empirical-Bayes moderated), comes with FDR.
#                           This is the INFERENTIAL differential call.
#     (2) CPM log2ratio    — descriptive log2(B10/Parental) of per-library
#                           CPM-normalized mean signal. No replicate model,
#                           no p-value. This is the ML FEATURE.
#   If (2) tracks (1), the ML features are not a normalization artifact — they
#   encode the same change the rigorous test detects. This is the figure that
#   connects your feature engineering to your differential analysis, and it is
#   NOT redundant with RF/SHAP (which never touch csaw output).
#
# WHAT IT DOES
#   1. Reads one csaw *_all_regions.csv (ONE mark, ONE contrast).
#   2. Reads the 5kb CPM signal table and picks the matching log2ratio column.
#   3. Joins by GENOMIC OVERLAP (not string match): each csaw region is matched
#      to the 5kb bins it physically overlaps (GenomicRanges::findOverlaps).
#   4. Collapses overlapping 5kb bins to ONE CPM value per csaw region, by
#      OVERLAP-WIDTH-WEIGHTED mean (a bin contributing more bp counts more).
#   5. Correlates csaw rep.logFC vs the aggregated CPM log2ratio:
#        - Spearman (PRIMARY: rank agreement; robust to csaw's shrinkage)
#        - Pearson  (secondary; expect <1 slope because CPM is unshrunk)
#   6. Writes the merged per-region table, a stats summary, and a scatter plot
#      (optionally colored by csaw FDR significance).
#
# WHY 5kb, NOT 250kb
#   csaw regions are kb-scale. Joined to 250kb bins, dozens of distinct regions
#   collapse onto one bin and you'd correlate each region against a 250kb mean
#   dominated by OTHER regions. The 5kb grid keeps the CPM value local to each
#   csaw region. (This is the one analysis where the 5kb extraction is the
#   right input.)
#
# DIRECTIONALITY — READ THIS
#   csaw logFC sign depends on which group was numerator (group1). Your run
#   'h3k9me3 B10 vs Parental' has B10 as numerator -> +logFC = higher in B10.
#   The CPM column 'log2ratio_p-b10_<mark>_mean' is log2(B10/Parental) ->
#   +ratio = higher in B10. THESE AGREE. If your csaw contrast had Parental as
#   numerator, the signs would be flipped and the correlation would come out
#   NEGATIVE — use --flip_cpm to align them. The script prints the sign
#   convention it assumed so you can verify.
#
# CAVEAT (state in the paper)
#   csaw normalizes by TMM; the CPM bigWigs are per-library depth only. Residual
#   scatter reflects BOTH csaw's empirical-Bayes shrinkage AND this
#   normalization difference, not just noise. For a migration-axis mark with no
#   global shift (e.g. H3K9me3 B10-vs-Parental) TMM's assumptions are
#   reasonable, so agreement should be good; large divergence is itself
#   informative (possible global shift -> needs spike-in; cf. Orlando 2014).
#
# USAGE
#   Rscript scripts/analysis/csaw_vs_cpm_agreement.R \
#     --csaw   results/diff_analysis/csaw/H3K9me3_B10-none_vs_Parental-none_libsize_all_regions.csv \
#     --cpm    data/processed/cut_run_signal_5kb.txt \
#     --mark   H3K9me3 \
#     --outdir results/exploratory/csaw_vs_cpm
#
#   # if your csaw contrast used Parental as numerator (signs flipped):
#   Rscript ... --flip_cpm
#
#   # override the CPM column name explicitly if auto-detection misses:
#   Rscript ... --cpm_col log2ratio_p-b10_h3k9me3_mean
#
# OUTPUTS (into --outdir, stem = <mark>_<csaw basename>)
#   <stem>_merged.txt        per-region: coords, csaw logFC/FDR, CPM ratio, n_bins
#   <stem>_summary.txt       Spearman/Pearson rho, p, n, slope, sign convention
#   <stem>_scatter.pdf/png   csaw logFC (x) vs CPM log2ratio (y), with fits
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list = list(
  make_option("--csaw", type = "character",
              help = "csaw *_all_regions.csv for ONE mark, ONE contrast (required)"),
  make_option("--cpm", type = "character",
              help = "5kb CPM signal table (cut_run_signal_5kb.txt) (required)"),
  make_option("--mark", type = "character",
              help = "Mark name, e.g. H3K9me3 — used to auto-pick the CPM column (required)"),
  make_option("--cpm_col", type = "character", default = NULL,
              help = "Explicit CPM log2ratio column name (overrides auto-detect)"),
  make_option("--csaw_lfc_col", type = "character", default = "rep.logFC",
              help = "csaw logFC column to use [default %default]"),
  make_option("--fdr_col", type = "character", default = "FDR",
              help = "csaw FDR column for significance coloring [default %default]"),
  make_option("--fdr_thresh", type = "double", default = 0.05,
              help = "FDR threshold for the 'significant' color [default %default]"),
  make_option("--flip_cpm", action = "store_true", default = FALSE,
              help = "Negate the CPM ratio (use if csaw numerator is the OTHER condition)"),
  make_option("--min_bins", type = "integer", default = 1L,
              help = "Min overlapping 5kb bins required to keep a csaw region [default %default]"),
  make_option("--outdir", type = "character", default = "results/exploratory/csaw_vs_cpm")
)))
for (req in c("csaw", "cpm", "mark"))
  if (is.null(opt[[req]])) stop("--", req, " is required.")
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Load csaw regions and build a GRanges. csaw CSV has seqnames/start/end.
# -----------------------------------------------------------------------------
cat("[INFO] Reading csaw regions:", opt$csaw, "\n")
cs <- read.csv(opt$csaw, stringsAsFactors = FALSE)
need <- c("seqnames", "start", "end", opt$csaw_lfc_col)
miss <- setdiff(need, colnames(cs))
if (length(miss)) stop("csaw file missing columns: ", paste(miss, collapse = ", "))
cs$.row <- seq_len(nrow(cs))   # stable id to map overlaps back
cat(sprintf("[INFO]   %d csaw regions\n", nrow(cs)))

cs_gr <- GRanges(seqnames = cs$seqnames,
                 ranges   = IRanges(start = cs$start, end = cs$end))
mcols(cs_gr)$.row <- cs$.row

# -----------------------------------------------------------------------------
# Load 5kb CPM table; pick the log2ratio column for this mark.
# Expected column form: log2ratio_p-b10_<mark>_mean  (migration axis).
# -----------------------------------------------------------------------------
cat("[INFO] Reading 5kb CPM table:", opt$cpm, "\n")
cpm <- read.csv(opt$cpm, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
req_cpm <- c("chr", "start", "end")
miss <- setdiff(req_cpm, colnames(cpm))
if (length(miss)) stop("CPM table missing columns: ", paste(miss, collapse = ", "))

# Auto-detect the migration log2ratio column for this mark (case-insensitive),
# unless the user gave one explicitly. We specifically want the p-b10 migration
# axis, NOT the drug (unc) axis.
pick_cpm_col <- function(cols, mark, explicit) {
  if (!is.null(explicit)) {
    if (!explicit %in% cols) stop("--cpm_col '", explicit, "' not in CPM table.")
    return(explicit)
  }
  want <- tolower(paste0("log2ratio_p-b10_", mark, "_mean"))
  hit  <- cols[tolower(cols) == want]
  if (length(hit) == 1) return(hit)
  # looser fallback: any migration-axis log2ratio containing the mark
  cand <- cols[grepl("^log2ratio_p-b10_", tolower(cols)) &
               grepl(tolower(mark), tolower(cols)) &
               !grepl("unc", tolower(cols))]
  if (length(cand) >= 1) {
    if (length(cand) > 1)
      warning("Multiple candidate CPM columns; using first: ", cand[1])
    return(cand[1])
  }
  stop("Could not find a migration log2ratio column for mark '", mark,
       "'. Looked for '", want, "'. Pass --cpm_col explicitly.\n",
       "  Available log2ratio columns: ",
       paste(cols[grepl("log2ratio", tolower(cols))], collapse = ", "))
}
cpm_col <- pick_cpm_col(colnames(cpm), opt$mark, opt$cpm_col)
cat("[INFO]   CPM column:", cpm_col, "\n")

cpm$.cpm_ratio <- suppressWarnings(as.numeric(cpm[[cpm_col]]))
if (opt$flip_cpm) {
  cpm$.cpm_ratio <- -cpm$.cpm_ratio
  cat("[INFO]   --flip_cpm applied: CPM ratio negated to match csaw numerator.\n")
}

cpm_gr <- GRanges(seqnames = cpm$chr,
                  ranges   = IRanges(start = cpm$start, end = cpm$end))
mcols(cpm_gr)$ratio <- cpm$.cpm_ratio

# -----------------------------------------------------------------------------
# Overlap join: every (csaw region, 5kb bin) intersection.
# For each csaw region, aggregate the overlapping bins' CPM ratio by
# OVERLAP-WIDTH-WEIGHTED mean (width of the actual intersection, not the bin).
# -----------------------------------------------------------------------------
cat("[INFO] Computing genomic overlaps (csaw regions x 5kb bins)...\n")
ov <- findOverlaps(cs_gr, cpm_gr, ignore.strand = TRUE)
if (length(ov) == 0)
  stop("No overlaps between csaw regions and CPM bins. Check that both use the ",
       "same chromosome naming (e.g. 'chr1' vs '1').")

# Width of each intersection = how much of the bin actually falls in the region.
inter <- pintersect(cs_gr[queryHits(ov)], cpm_gr[subjectHits(ov)],
                    ignore.strand = TRUE)
w     <- width(inter)
qid   <- mcols(cs_gr)$.row[queryHits(ov)]      # csaw region id
ratio <- mcols(cpm_gr)$ratio[subjectHits(ov)]  # bin CPM ratio

# Drop intersections where the bin's CPM ratio is NA (can't weight a missing value).
ok    <- !is.na(ratio) & w > 0
qid   <- qid[ok]; ratio <- ratio[ok]; w <- w[ok]

# Width-weighted mean per csaw region, plus a count of contributing bins.
num   <- tapply(ratio * w, qid, sum)
den   <- tapply(w,         qid, sum)
nbin  <- tapply(rep(1L, length(qid)), qid, sum)
agg   <- data.frame(.row = as.integer(names(num)),
                    cpm_ratio_wmean = as.numeric(num / den),
                    n_bins = as.integer(nbin[as.character(names(num))]))

# -----------------------------------------------------------------------------
# Merge aggregated CPM back onto csaw regions; keep regions with enough bins.
# -----------------------------------------------------------------------------
m <- merge(cs, agg, by = ".row", all.x = FALSE)
m <- m[m$n_bins >= opt$min_bins, ]
m$csaw_logFC <- suppressWarnings(as.numeric(m[[opt$csaw_lfc_col]]))
m <- m[is.finite(m$csaw_logFC) & is.finite(m$cpm_ratio_wmean), ]
if (nrow(m) < 10)
  stop("Too few joined regions (", nrow(m), ") for a meaningful correlation.")
cat(sprintf("[INFO] Joined regions used: %d\n", nrow(m)))

# Significance flag for coloring (if FDR column present).
has_fdr <- opt$fdr_col %in% colnames(m)
if (has_fdr)
  m$sig <- ifelse(suppressWarnings(as.numeric(m[[opt$fdr_col]])) < opt$fdr_thresh,
                  paste0("FDR<", opt$fdr_thresh), "n.s.")

# -----------------------------------------------------------------------------
# Correlations. Spearman is PRIMARY (rank agreement; robust to csaw shrinkage,
# which compresses magnitudes nonlinearly). Pearson is secondary.
# -----------------------------------------------------------------------------
sp <- suppressWarnings(cor.test(m$csaw_logFC, m$cpm_ratio_wmean,
                                method = "spearman"))
pe <- suppressWarnings(cor.test(m$csaw_logFC, m$cpm_ratio_wmean,
                                method = "pearson"))
fit   <- lm(cpm_ratio_wmean ~ csaw_logFC, data = m)
slope <- unname(coef(fit)[2])

cat(sprintf("[RESULT] Spearman rho = %.3f  (p = %.3g)\n",
            unname(sp$estimate), sp$p.value))
cat(sprintf("[RESULT] Pearson  r   = %.3f  (p = %.3g)\n",
            unname(pe$estimate), pe$p.value))
cat(sprintf("[RESULT] OLS slope (CPM ~ csaw) = %.3f (expect <1: CPM is unshrunk, csaw is shrunk)\n", slope))

# -----------------------------------------------------------------------------
# Write merged table and summary.
# -----------------------------------------------------------------------------
stem <- paste0(opt$mark, "_", sub("\\.csv$", "", basename(opt$csaw)))
keep_cols <- c("seqnames", "start", "end", "width",
               opt$csaw_lfc_col, opt$fdr_col, "csaw_logFC",
               "cpm_ratio_wmean", "n_bins")
keep_cols <- keep_cols[keep_cols %in% colnames(m)]
write.table(m[, keep_cols], file.path(opt$outdir, paste0(stem, "_merged.txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

summary_df <- data.frame(
  mark            = opt$mark,
  csaw_file       = basename(opt$csaw),
  cpm_column      = cpm_col,
  cpm_flipped     = opt$flip_cpm,
  n_regions       = nrow(m),
  spearman_rho    = unname(sp$estimate),
  spearman_p      = sp$p.value,
  pearson_r       = unname(pe$estimate),
  pearson_p       = pe$p.value,
  ols_slope       = slope,
  sign_convention = "positive = higher in csaw group1 (numerator) = higher in B10 (default)"
)
write.table(summary_df, file.path(opt$outdir, paste0(stem, "_summary.txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("[INFO] Wrote merged table and summary to", opt$outdir, "\n")

# -----------------------------------------------------------------------------
# Scatter plot: csaw logFC (x, inferential) vs CPM log2ratio (y, descriptive).
# Reference y=x (identity) and the OLS fit are both drawn; the gap between them
# visualizes the shrinkage/normalization difference.
# -----------------------------------------------------------------------------
subtitle <- sprintf(
  "Spearman rho = %.3f (p = %.2g)  |  Pearson r = %.3f  |  n = %d regions",
  unname(sp$estimate), sp$p.value, unname(pe$estimate), nrow(m))

p <- ggplot(m, aes(x = csaw_logFC, y = cpm_ratio_wmean)) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.4) +
  {if (has_fdr)
     geom_point(aes(color = sig), size = 0.9, alpha = 0.5)
   else
     geom_point(size = 0.9, alpha = 0.5, color = "#1A5276")} +
  geom_smooth(method = "lm", se = TRUE, color = "#C0392B",
              linewidth = 0.9, formula = y ~ x) +
  scale_color_manual(values = setNames(c("#C0392B", "grey70"),
                       c(paste0("FDR<", opt$fdr_thresh), "n.s.")),
                     name = "csaw") +
  labs(
    title    = sprintf("csaw differential vs CPM log2-ratio: %s (B10 vs Parental)",
                       opt$mark),
    subtitle = subtitle,
    x = sprintf("csaw %s  (inferential, shrunk; +=higher in B10)", opt$csaw_lfc_col),
    y = "CPM log2(B10/Parental), width-weighted  (descriptive ML feature)"
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey30"),
        legend.position = "right")

ggsave(file.path(opt$outdir, paste0(stem, "_scatter.pdf")), p,
       width = 7, height = 5.5, dpi = 600)
ggsave(file.path(opt$outdir, paste0(stem, "_scatter.png")), p,
       width = 7, height = 5.5, dpi = 600)
cat("[INFO] Wrote scatter plot.\n")

cat(sprintf("\n[INFO] Done -> %s\n", normalizePath(opt$outdir)))
cat("[INTERPRET] High Spearman rho => CPM-derived ML features track the ",
    "FDR-controlled differential calls (features are real change, not ",
    "normalization artifact). The dotted line is y=x; the red line is the ",
    "actual fit — the gap reflects csaw shrinkage + TMM-vs-CPM normalization, ",
    "which is expected, not error.\n", sep = "")