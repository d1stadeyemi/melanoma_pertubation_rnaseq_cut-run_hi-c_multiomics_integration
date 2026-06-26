#!/usr/bin/env Rscript
# =============================================================================
# csaw_overlap_venn.R — Cross-condition overlap of csaw differential regions
# =============================================================================
#
# PROJECT: A375 melanoma parental-vs-B10 constricted-migration multi-omics.
#
# THE QUESTION ----------------------------------------------------------------
#   You ran csaw.R twice for the SAME mark under the SAME normalization, once
#   per cell line, e.g.:
#       A = Parental UNC effect :  --group1 "Parental|UNC" --group2 "Parental|none"
#       B = B10      UNC effect :  --group1 "B10|UNC"      --group2 "B10|none"
#   Each run wrote per-region tables. This script asks: WHICH regions respond to
#   UNC1999 in BOTH cell lines, which are cell-line-specific, and — crucially —
#   do the shared regions move the SAME direction (conserved loss) or OPPOSITE
#   directions (discordant)?
#
#   It is a general 2-set comparator: A and B can be any two csaw region tables
#   (two cell lines, two marks, two time points, ...). "up"/"down" are taken from
#   the sign of logFC (relabel with --up-label/--down-label, e.g. gained/lost).
#
# WHY "CONSENSUS LOCI" (read this) --------------------------------------------
#   The merged windows (regions) from two independent csaw runs are rarely
#   base-identical, so "is this region in both?" has no unique answer at the
#   single-region level. We therefore cluster regions from A and B into CONSENSUS
#   LOCI: build a graph whose edges join an A region and a B region that overlap
#   by at least --min-overlap (reciprocally, by default), then take connected
#   components. Each component is ONE consensus locus, classified as:
#       shared  = supported by >=1 A region AND >=1 B region
#       A-only  = supported by A region(s) only
#       B-only  = supported by B region(s) only
#   Every Venn count is then a single, symmetric number (a count of loci), not an
#   asymmetric "A-side vs B-side" pair. Raw region counts are also reported.
#
# WHAT IT PRODUCES ------------------------------------------------------------
#   THREE VENNS (your request), each computed INDEPENDENTLY on its own subset:
#     venn_all   — all significant regions (direction ignored)
#     venn_up    — only up (logFC>0) significant regions
#     venn_down  — only down (logFC<0) significant regions
#     NOTE: a DISCORDANT shared locus (up in A, down in B) is "shared" in
#     venn_all but appears as A-only in venn_up and B-only in venn_down — so
#     all-shared is NOT simply up-shared + down-shared. That gap IS the
#     discordant signal, quantified explicitly below.
#
#   DIRECTION CONCORDANCE (added) — at shared loci (from the all-clustering):
#     concordant_up / concordant_down / discordant, as a bar + a 2x2 contingency.
#     For an EZH2 inhibitor on H3K27me3 you EXPECT concordant loss to dominate;
#     discordant loci are the interesting minority.
#
#   logFC CONCORDANCE SCATTER (added) — A logFC vs B logFC per consensus locus,
#     coloured by category, with the no-change diagonal, quadrant counts, and a
#     Spearman/Pearson correlation on shared loci. This is the quantitative
#     summary of all three Venns at once: a tight y=x cloud = conserved response;
#     off-diagonal = divergence. (Full scatter needs *_all_regions.csv inputs so
#     a non-significant side still has a logFC; with *_significant.csv only, the
#     non-significant side is NA and those points are dropped — see INPUTS.)
#
#   OVERLAP ENRICHMENT (added) — Jaccard index of consensus loci, plus a Fisher
#     test of association on CO-TESTED loci (loci with a tested region in BOTH
#     runs): is "significant in A" associated with "significant in B" beyond
#     chance? Requires *_all_regions.csv for both (the tested universe); skipped
#     with a notice for *_significant.csv inputs.
#
#   GENE LISTS (your request, extended) — with --genes, one-gene-per-line lists
#     for every category, ready for the pathway-enrichment script:
#       shared_all / shared_concordant_up / shared_concordant_down / shared_discordant
#       <A>_only_all / _up / _down ;  <B>_only_all / _up / _down
#     plus background_universe.txt (genes in ANY input region).
#
#   MASTER TABLE — consensus_loci.tsv: every locus with coordinates, per-side
#     logFC and direction, category, concordance, region counts, and mapped genes.
#
# INPUTS ----------------------------------------------------------------------
#   --a / --b : csaw region CSVs for the two conditions. Either the
#               *_all_regions.csv (RECOMMENDED — enables the full scatter and the
#               Fisher enrichment) or the *_significant.csv (Venns + concordance +
#               genes still work; scatter is shared-only; Fisher skipped). Both
#               must be the SAME mark and SAME normalization for a fair compare.
#               Required columns: chr,start,end,logFC,FDR (csaw output has these).
#   Significance is (re)derived as FDR < --fdr AND |logFC| > --lfc, so it is
#   consistent across the two files regardless of which thresholds each run used.
#
# ENV: conda env multiomics-r  (optparse, GenomicRanges, IRanges, S4Vectors).
#      Base graphics only (no Venn package needed). Run on the cluster.
#
# ------------------------------------------------------------------------------
# USAGE EXAMPLES
# ------------------------------------------------------------------------------
#   # (1) Parental vs B10 UNC effect, background norm, 50% reciprocal overlap,
#   #     with gene mapping (recommended: feed the *_all_regions.csv):
#   Rscript csaw_overlap_venn.R \
#       --a results/diff/H3K27me3_Parentalunc/background/H3K27me3_ParentalUNC_vs_Parentalnone_background_all_regions.csv \
#       --b results/diff/H3K27me3_B10unc/background/H3K27me3_B10UNC_vs_B10none_background_all_regions.csv \
#       --a-name Parental --b-name B10 \
#       --up-label gained --down-label lost \
#       --genes data/processed/GeneListTSShg19_cleaned.txt \
#       --outdir results/diff/overlap_H3K27me3_UNC_Parental_vs_B10
#
#   # (2) Stricter overlap (70% reciprocal) and a different FDR/logFC cutoff:
#   Rscript csaw_overlap_venn.R --a A_all_regions.csv --b B_all_regions.csv \
#       --a-name Parental --b-name B10 --min-overlap 0.7 --fdr 0.01 --lfc 1.0 \
#       --outdir results/diff/overlap_strict
#
#   # (3) "Any overlap counts as shared" (1 bp), significant-only inputs:
#   Rscript csaw_overlap_venn.R --a A_significant.csv --b B_significant.csv \
#       --a-name Parental --b-name B10 --overlap-mode any \
#       --outdir results/diff/overlap_anyhit
# ------------------------------------------------------------------------------
# =============================================================================

suppressWarnings(suppressMessages({
  library(optparse)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
}))

# =============================================================================
# CLI
# =============================================================================
parse_cli <- function() {
  option_list <- list(
    make_option("--a", type = "character",
                help = "REQUIRED. csaw region CSV for condition A (*_all_regions.csv recommended, or *_significant.csv). Needs chr,start,end,logFC,FDR."),
    make_option("--b", type = "character",
                help = "REQUIRED. csaw region CSV for condition B (same mark + same normalization as A)."),
    make_option("--a-name", type = "character", default = "A", dest = "a_name",
                help = "Label for condition A in figures/files (e.g. Parental). [default %default]"),
    make_option("--b-name", type = "character", default = "B", dest = "b_name",
                help = "Label for condition B (e.g. B10). [default %default]"),
    make_option("--outdir", type = "character", default = "csaw_overlap_out",
                help = "Output directory. [default %default]"),

    # ---- significance definition (re-applied to both files for consistency) ----
    make_option("--fdr", type = "double", default = 0.05,
                help = "FDR cutoff defining a significant region. [default %default]"),
    make_option("--lfc", type = "double", default = 1.0,
                help = "Min |log2FC| defining a significant region. [default %default]"),

    # ---- overlap / matching ----
    make_option("--min-overlap", type = "double", default = 0.5, dest = "min_overlap",
                help = paste0("Fraction (0-1] of overlap required to join an A region and a B region into the same ",
                              "consensus locus. With --overlap-mode reciprocal (default) BOTH regions must be covered ",
                              ">= this fraction (like bedtools -f X -r). 0.5 = 50%% reciprocal. [default %default]")),
    make_option("--overlap-mode", type = "character", default = "reciprocal", dest = "overlap_mode",
                help = paste0("How --min-overlap is applied. ALLOWED: 'reciprocal' (both regions covered >= frac; ",
                              "recommended), 'any' (>=1 bp, ignores --min-overlap), 'a' (A covered >= frac), ",
                              "'b' (B covered >= frac). [default %default]")),

    # ---- direction labels (display + filenames only; direction itself = logFC sign) ----
    make_option("--up-label", type = "character", default = "up", dest = "up_label",
                help = "Label for logFC>0 regions (e.g. 'gained' for a contrast). [default %default]"),
    make_option("--down-label", type = "character", default = "down", dest = "down_label",
                help = "Label for logFC<0 regions (e.g. 'lost'). [default %default]"),

    # ---- optional gene mapping ----
    make_option("--genes", type = "character", default = NULL,
                help = paste0("Optional hg19 gene table; maps genes to each locus category and writes one-gene-per-line ",
                              "lists into <outdir>/genelists/ (ready for the pathway-enrichment script). Tab/comma table ",
                              "with a gene-name column (gene_name|gene|symbol|name), a chromosome column (chr|chrom|seqnames; ",
                              "MUST match your region coordinates' naming, e.g. 'chr1'), and EITHER a TSS column (TSS) OR ",
                              "start+end columns (start/end|txStart/txEnd|bin_start/bin_end). [default: none]"))
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$a) || is.null(opt$b)) stop("--a and --b are both required.")
  if (!opt$overlap_mode %in% c("reciprocal", "any", "a", "b"))
    stop("--overlap-mode must be one of: reciprocal | any | a | b.")
  if (opt$min_overlap <= 0 || opt$min_overlap > 1) stop("--min-overlap must be in (0, 1].")
  opt
}

# Filesystem-safe token.
safe_token <- function(s) gsub("[^A-Za-z0-9]+", "", s)

# =============================================================================
# Logging
# =============================================================================
log_lines <- new.env(); log_lines$txt <- character(0)
say <- function(...) { line <- paste0(...); cat(line, "\n"); log_lines$txt <- c(log_lines$txt, line) }

# =============================================================================
# Read a csaw region CSV -> GRanges with logFC / FDR / sig / dir; flag whether a
# non-significant ("universe") fraction is present (i.e. an *_all_regions.csv).
# =============================================================================
read_region_csv <- function(path, fdr, lfc, name) {
  if (!file.exists(path)) stop(name, " file not found: ", path)
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("chr", "start", "end", "logFC", "FDR")
  miss <- setdiff(need, colnames(d))
  if (length(miss))
    stop(name, ": CSV missing column(s): ", paste(miss, collapse = ", "),
         ". Expected a csaw *_all_regions.csv or *_significant.csv.")
  d$chr   <- as.character(d$chr)
  d$start <- as.integer(d$start); d$end <- as.integer(d$end)
  d$logFC <- as.numeric(d$logFC); d$FDR <- as.numeric(d$FDR)
  ok <- is.finite(d$start) & is.finite(d$end) & is.finite(d$logFC) & is.finite(d$FDR)
  if (any(!ok)) { say(sprintf("[READ] %s: dropping %d rows with non-finite coords/logFC/FDR.", name, sum(!ok))); d <- d[ok, , drop = FALSE] }
  sig <- d$FDR < fdr & abs(d$logFC) > lfc
  gr <- GRanges(d$chr, IRanges(d$start, d$end))
  gr$logFC <- d$logFC; gr$FDR <- d$FDR; gr$sig <- sig
  gr$dir <- ifelse(d$logFC > 0, "up", "down")
  has_univ <- any(!sig)                       # any non-significant row -> tested universe present
  say(sprintf("[READ] %-10s %s: %d regions; %d significant (%d up, %d down)%s",
              name, basename(path), length(gr), sum(sig),
              sum(sig & gr$dir == "up"), sum(sig & gr$dir == "down"),
              if (has_univ) "  [universe available]" else "  [significant-only]"))
  list(gr = gr, has_universe = has_univ)
}

# =============================================================================
# Reciprocal-overlap graph clustering of two region sets into consensus loci
# =============================================================================

# Qualifying A<->B overlap edges under the chosen mode; returns an n x 2 matrix
# of (A index, B index).
overlap_edges <- function(grA, grB, frac, mode) {
  if (!length(grA) || !length(grB)) return(matrix(integer(0), ncol = 2))
  hits <- findOverlaps(grA, grB)
  if (!length(hits)) return(matrix(integer(0), ncol = 2))
  qa <- queryHits(hits); sb <- subjectHits(hits)
  ov <- pmin(end(grA)[qa], end(grB)[sb]) - pmax(start(grA)[qa], start(grB)[sb]) + 1
  ov <- pmax(ov, 0)
  fa <- ov / width(grA)[qa]                   # fraction of the A region covered
  fb <- ov / width(grB)[sb]                   # fraction of the B region covered
  keep <- switch(mode,
                 reciprocal = fa >= frac & fb >= frac,
                 any        = ov > 0,
                 a          = fa >= frac,
                 b          = fb >= frac,
                 stop("Unknown overlap mode: ", mode))
  cbind(qa[keep], sb[keep])
}

# Union-find (disjoint-set) over n nodes given an m x 2 edge matrix; returns the
# root id of every node (so split(seq_len(n), roots) gives the components).
union_find <- function(n, E) {
  if (n == 0) return(integer(0))
  parent <- seq_len(n)
  find <- function(i) {
    while (parent[i] != i) { parent[i] <<- parent[parent[i]]; i <- parent[i] }  # path halving (<<- targets parent)
    i
  }
  if (length(E) && nrow(E)) for (k in seq_len(nrow(E))) {
    ra <- find(E[k, 1]); rb <- find(E[k, 2]); if (ra != rb) parent[ra] <- rb    # local assign in this scope
  }
  vapply(seq_len(n), find, numeric(1))                                          # numeric: robust to int/double E
}

# Cluster A and B regions into consensus loci. Returns a list of components, each
# a list(a_idx = indices into grA, b_idx = indices into grB).
cluster_two_sets <- function(grA, grB, frac, mode) {
  nA <- length(grA); nB <- length(grB); N <- nA + nB
  if (N == 0) return(list())
  E <- overlap_edges(grA, grB, frac, mode)
  if (nrow(E)) E[, 2] <- E[, 2] + nA            # offset B nodes into the shared id space
  roots <- union_find(N, E)
  grp <- split(seq_len(N), roots)
  lapply(unname(grp), function(nodes)
    list(a_idx = nodes[nodes <= nA], b_idx = nodes[nodes > nA] - nA))
}

# Venn counts (a_only / shared / b_only) from a clustering.
venn_counts <- function(clusters) {
  if (!length(clusters)) return(c(a_only = 0L, shared = 0L, b_only = 0L))
  cat <- vapply(clusters, function(cl)
    if (length(cl$a_idx) && length(cl$b_idx)) "shared"
    else if (length(cl$a_idx)) "a_only" else "b_only", character(1))
  c(a_only = sum(cat == "a_only"), shared = sum(cat == "shared"), b_only = sum(cat == "b_only"))
}

jaccard <- function(v) if (sum(v) == 0) NA_real_ else unname(v["shared"] / sum(v))

# =============================================================================
# Per-locus helpers
# =============================================================================
# Width-weighted mean logFC of every region in `full_gr` overlapping `locus`
# (a single-range GRanges). NA if none. Used so a non-significant side still has
# a logFC for the scatter, when an *_all_regions.csv was supplied.
agg_logfc <- function(locus, full_gr) {
  if (!length(full_gr)) return(NA_real_)
  h <- findOverlaps(locus, full_gr)
  if (!length(h)) return(NA_real_)
  idx <- subjectHits(h)
  ovw <- pmin(end(locus), end(full_gr)[idx]) - pmax(start(locus), start(full_gr)[idx]) + 1
  ovw <- pmax(ovw, 1)
  sum(full_gr$logFC[idx] * ovw) / sum(ovw)
}

# Dominant direction (width-weighted logFC sign) of a set of significant regions.
dir_of <- function(gr) {
  if (!length(gr)) return(NA_character_)
  w <- width(gr); m <- sum(gr$logFC * w) / sum(w)
  if (m > 0) "up" else "down"
}

# =============================================================================
# Optional gene annotation (flexible TSS-or-body reader, like csaw.R)
# =============================================================================
resolve_genes <- function(opt) {
  if (is.null(opt$genes)) return(NULL)
  if (!file.exists(opt$genes)) stop("--genes file not found: ", opt$genes)
  first <- readLines(opt$genes, n = 1L, warn = FALSE)
  sep   <- if (grepl("\t", first)) "\t" else ","
  tab   <- read.table(opt$genes, header = TRUE, sep = sep, quote = "\"",
                      stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE)
  name_col <- intersect(c("gene_name", "gene", "symbol", "Symbol", "name", "Gene", "GeneSymbol"), colnames(tab))
  if (!length(name_col)) stop("--genes: no gene-name column (gene_name/gene/symbol/name). Columns: ", paste(colnames(tab), collapse = ", "))
  chr_col  <- intersect(c("chr", "chrom", "Chr", "Chrom", "seqnames", "chromosome"), colnames(tab))
  if (!length(chr_col)) stop("--genes: no chromosome column (chr/chrom/seqnames). Columns: ", paste(colnames(tab), collapse = ", "))
  tss_col   <- intersect(c("TSS", "tss", "Tss"), colnames(tab))
  start_col <- intersect(c("start", "Start", "txStart", "bin_start", "geneStart"), colnames(tab))
  end_col   <- intersect(c("end", "End", "txEnd", "bin_end", "geneEnd"), colnames(tab))
  chr  <- as.character(tab[[chr_col[1]]]); gene <- as.character(tab[[name_col[1]]])
  if (length(tss_col)) {
    pos  <- suppressWarnings(as.integer(tab[[tss_col[1]]]))
    keep <- !is.na(pos) & !is.na(chr) & nzchar(gene)
    gr   <- GRanges(chr[keep], IRanges(pmax(pos[keep], 1L), width = 1L)); gr$gene <- gene[keep]
    say(sprintf("[GENES] %d genes loaded by TSS from %s.", length(gr), basename(opt$genes)))
  } else if (length(start_col) && length(end_col)) {
    s <- suppressWarnings(as.integer(tab[[start_col[1]]])); e <- suppressWarnings(as.integer(tab[[end_col[1]]]))
    keep <- !is.na(s) & !is.na(e) & !is.na(chr) & nzchar(gene)
    gr   <- GRanges(chr[keep], IRanges(pmax(pmin(s[keep], e[keep]), 1L), pmax(s[keep], e[keep]))); gr$gene <- gene[keep]
    say(sprintf("[GENES] %d genes loaded by body from %s.", length(gr), basename(opt$genes)))
  } else stop("--genes: need a TSS column, or BOTH start and end. Columns: ", paste(colnames(tab), collapse = ", "))
  if (!length(gr)) stop("--genes: no usable rows in ", basename(opt$genes))
  gr
}

# Genes whose feature overlaps any region in `loci_gr` (a GRanges).
genes_in <- function(loci_gr, genes_gr) {
  if (is.null(genes_gr) || !length(loci_gr)) return(character(0))
  sort(unique(genes_gr$gene[queryHits(findOverlaps(genes_gr, loci_gr))]))
}

# =============================================================================
# Figures (base graphics; pdf + png, like csaw.R)
# =============================================================================
save_fig <- function(fn, outdir, draw, w = 6.5, h = 5) {
  for (dev in c("pdf", "png")) {
    if (dev == "pdf") pdf(file.path(outdir, paste0(fn, ".pdf")), width = w, height = h)
    else png(file.path(outdir, paste0(fn, ".png")), width = w, height = h, units = "in", res = 200)
    draw(); dev.off()
  }
}

COL_A <- "#1F6F8B"; COL_B <- "#E08A3C"; COL_SHARED <- "#6D5BA6"
COL_UP <- "#C0392B"; COL_DOWN <- "#1F6F8B"; COL_DISC <- "#7F8C8D"

# Dependency-free 2-set Venn (fixed overlapping circles; counts in the 3 zones).
draw_venn2 <- function(cnt, labA, labB, title, subtitle = "") {
  na <- cnt["a_only"]; ns <- cnt["shared"]; nb <- cnt["b_only"]
  draw <- function() {
    par(mar = c(3, 1, 4, 1))
    plot.new(); plot.window(xlim = c(0, 10), ylim = c(0, 7), asp = 1)
    circ <- function(cx, cy, r, col) { th <- seq(0, 2*pi, length.out = 240)
      polygon(cx + r*cos(th), cy + r*sin(th), col = col, border = "grey25", lwd = 1.5) }
    r <- 2.5; cy <- 3.4; cxA <- 3.7; cxB <- 6.3
    circ(cxA, cy, r, adjustcolor(COL_A, 0.40))
    circ(cxB, cy, r, adjustcolor(COL_B, 0.40))
    text(cxA - 1.25, cy, na, cex = 1.7, font = 2)
    text(cxB + 1.25, cy, nb, cex = 1.7, font = 2)
    text((cxA + cxB)/2, cy, ns, cex = 1.7, font = 2, col = "white")
    text(cxA - 0.4, cy + r + 0.5, sprintf("%s\n(%d loci)", labA, na + ns), col = COL_A, font = 2, cex = 1.05, adj = 1)
    text(cxB + 0.4, cy + r + 0.5, sprintf("%s\n(%d loci)", labB, nb + ns), col = COL_B, font = 2, cex = 1.05, adj = 0)
    title(main = title, line = 2.2)
    if (nzchar(subtitle)) mtext(subtitle, side = 1, line = 1, cex = 0.85, col = "grey30")
  }
  draw
}

# logFC concordance scatter: A vs B per consensus locus, coloured by category.
scatter_fig <- function(master, opt, outdir) {
  m <- master[is.finite(master$a_logFC) & is.finite(master$b_logFC), , drop = FALSE]
  if (nrow(m) < 3) { say("[SCATTER] <3 loci with logFC on both sides — scatter skipped (feed *_all_regions.csv to populate the non-significant side)."); return(invisible(NULL)) }
  colmap <- c(shared = COL_SHARED, a_only = COL_A, b_only = COL_B)
  sh <- m[m$category == "shared", , drop = FALSE]
  rho <- if (nrow(sh) >= 3) suppressWarnings(cor(sh$a_logFC, sh$b_logFC, method = "spearman", use = "complete.obs")) else NA_real_
  rp  <- if (nrow(sh) >= 3) suppressWarnings(cor(sh$a_logFC, sh$b_logFC, method = "pearson",  use = "complete.obs")) else NA_real_
  # shared-locus quadrant counts (concordance read straight off the plot)
  q1 <- sum(sh$a_logFC > 0 & sh$b_logFC > 0); q3 <- sum(sh$a_logFC < 0 & sh$b_logFC < 0)
  q2 <- sum(sh$a_logFC < 0 & sh$b_logFC > 0); q4 <- sum(sh$a_logFC > 0 & sh$b_logFC < 0)
  rng <- range(c(m$a_logFC, m$b_logFC), finite = TRUE)
  draw <- function() {
    par(mar = c(5, 5, 4, 8), xpd = FALSE)
    plot(m$a_logFC, m$b_logFC, pch = 16, cex = 0.5, col = adjustcolor(colmap[m$category], 0.6),
         xlim = rng, ylim = rng, xlab = sprintf("%s  log2FC", opt$a_name),
         ylab = sprintf("%s  log2FC", opt$b_name),
         main = sprintf("UNC effect concordance: %s vs %s", opt$a_name, opt$b_name))
    abline(0, 1, col = "grey40", lwd = 1.3); abline(h = 0, v = 0, col = "grey80", lty = 3)
    usr <- par("usr")
    text(usr[2]*0.92, usr[4]*0.92, sprintf("both %s: %d", opt$up_label,   q1), col = COL_UP,  cex = 0.85, adj = 1)
    text(usr[1]*0.92, usr[3]*0.92, sprintf("both %s: %d", opt$down_label, q3), col = COL_DOWN, cex = 0.85, adj = 0)
    text(usr[1]*0.92, usr[4]*0.92, sprintf("discord: %d", q2), col = COL_DISC, cex = 0.85, adj = 0)
    text(usr[2]*0.92, usr[3]*0.92, sprintf("discord: %d", q4), col = COL_DISC, cex = 0.85, adj = 1)
    lab <- c(sprintf("shared (rho=%s, r=%s)", ifelse(is.na(rho), "NA", sprintf("%.2f", rho)), ifelse(is.na(rp), "NA", sprintf("%.2f", rp))),
             sprintf("%s only", opt$a_name), sprintf("%s only", opt$b_name))
    par(xpd = TRUE)
    legend(usr[2] + 0.02*diff(usr[1:2]), usr[4], bty = "n", pch = 16, cex = 0.8,
           col = c(COL_SHARED, COL_A, COL_B), legend = lab)
  }
  save_fig("logfc_concordance_scatter", outdir, draw, w = 7.5, h = 5.5)
  list(rho = rho, pearson = rp, q_bothup = q1, q_bothdown = q3, q_disc = q2 + q4)
}

# Direction-concordance bar at shared loci.
concordance_bar <- function(master, opt, outdir) {
  sh <- master[master$category == "shared" & !is.na(master$concordance), , drop = FALSE]
  if (!nrow(sh)) { say("[CONC] no shared loci — concordance bar skipped."); return(invisible(NULL)) }
  n_cu <- sum(sh$concordance == "concordant_up")
  n_cd <- sum(sh$concordance == "concordant_down")
  n_di <- sum(sh$concordance == "discordant")
  vals <- c(n_cu, n_cd, n_di)
  labs <- c(sprintf("concordant\n%s", opt$up_label), sprintf("concordant\n%s", opt$down_label), "discordant")
  draw <- function() {
    par(mar = c(4, 5, 4, 1))
    bp <- barplot(vals, names.arg = labs, col = c(COL_UP, COL_DOWN, COL_DISC),
                  ylab = "shared consensus loci", main = sprintf("Direction concordance at shared loci\n(%s vs %s)", opt$a_name, opt$b_name),
                  ylim = c(0, max(vals) * 1.15 + 1))
    text(bp, vals, labels = vals, pos = 3, font = 2)
  }
  save_fig("direction_concordance_bar", outdir, draw, w = 5.5, h = 5)
  c(concordant_up = n_cu, concordant_down = n_cd, discordant = n_di)
}

# =============================================================================
# Overlap enrichment (Fisher) on CO-TESTED consensus loci (needs both universes)
# =============================================================================
overlap_enrichment <- function(A, B, opt, outdir) {
  if (!(A$has_universe && B$has_universe)) {
    say("[ENRICH] Fisher enrichment needs *_all_regions.csv for BOTH inputs (a tested universe) — skipped. Jaccard on significant loci is still reported.")
    return(NULL)
  }
  uni <- cluster_two_sets(A$gr, B$gr, opt$min_overlap, opt$overlap_mode)
  cot <- Filter(function(cl) length(cl$a_idx) && length(cl$b_idx), uni)   # tested in BOTH runs
  if (length(cot) < 4) { say("[ENRICH] too few co-tested loci for Fisher — skipped."); return(NULL) }
  sigA <- vapply(cot, function(cl) any(A$gr$sig[cl$a_idx]), logical(1))
  sigB <- vapply(cot, function(cl) any(B$gr$sig[cl$b_idx]), logical(1))
  tab  <- table(factor(sigA, c(TRUE, FALSE)), factor(sigB, c(TRUE, FALSE)))
  ft   <- suppressWarnings(fisher.test(tab))
  say(sprintf("[ENRICH] co-tested loci=%d | sig both=%d, %s-only=%d, %s-only=%d, neither=%d",
              length(cot), tab["TRUE","TRUE"], opt$a_name, tab["TRUE","FALSE"], opt$b_name, tab["FALSE","TRUE"], tab["FALSE","FALSE"]))
  say(sprintf("[ENRICH] Fisher OR=%.2f, p=%.3g (is significance in %s associated with significance in %s?)",
              unname(ft$estimate), ft$p.value, opt$a_name, opt$b_name))
  enr <- data.frame(
    metric = c("co_tested_loci", "sig_both", paste0("sig_", safe_token(opt$a_name), "_only"),
               paste0("sig_", safe_token(opt$b_name), "_only"), "sig_neither", "odds_ratio", "fisher_p"),
    value  = c(length(cot), tab["TRUE","TRUE"], tab["TRUE","FALSE"], tab["FALSE","TRUE"], tab["FALSE","FALSE"],
               round(unname(ft$estimate), 4), signif(ft$p.value, 4)))
  write.table(enr, file.path(outdir, "overlap_enrichment.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  list(odds_ratio = unname(ft$estimate), p = ft$p.value, table = tab)
}

# =============================================================================
# Build the master consensus-loci table from the ALL-significant clustering
# =============================================================================
build_master <- function(grA_sig, grB_sig, grA_full, grB_full, clusters) {
  if (!length(clusters)) return(data.frame())
  rows <- lapply(clusters, function(cl) {
    rs  <- c(grA_sig[cl$a_idx], grB_sig[cl$b_idx])
    chr <- as.character(seqnames(rs))[1]
    s <- min(start(rs)); e <- max(end(rs))
    locus <- GRanges(chr, IRanges(s, e))
    a_sig <- length(cl$a_idx) > 0; b_sig <- length(cl$b_idx) > 0
    category <- if (a_sig && b_sig) "shared" else if (a_sig) "a_only" else "b_only"
    a_dir <- if (a_sig) dir_of(grA_sig[cl$a_idx]) else NA_character_
    b_dir <- if (b_sig) dir_of(grB_sig[cl$b_idx]) else NA_character_
    conc  <- if (category == "shared") { if (a_dir == b_dir) paste0("concordant_", a_dir) else "discordant" } else NA_character_
    data.frame(chr = chr, start = s, end = e, width = e - s + 1, category = category,
               a_sig = a_sig, b_sig = b_sig,
               a_logFC = agg_logfc(locus, grA_full), b_logFC = agg_logfc(locus, grB_full),
               a_dir = a_dir, b_dir = b_dir, concordance = conc,
               n_a_regions = length(cl$a_idx), n_b_regions = length(cl$b_idx),
               stringsAsFactors = FALSE)
  })
  m <- do.call(rbind, rows)
  m[order(m$chr, m$start), , drop = FALSE]
}

# =============================================================================
# Gene lists per category, ready for the pathway-enrichment script
# =============================================================================
write_gene_lists <- function(master, genes_gr, opt, outdir) {
  if (is.null(genes_gr)) return(invisible(NULL))
  gldir <- file.path(outdir, "genelists"); dir.create(gldir, showWarnings = FALSE, recursive = TRUE)
  gr_of <- function(df) if (nrow(df)) GRanges(df$chr, IRanges(df$start, df$end)) else GRanges()
  aT <- safe_token(opt$a_name); bT <- safe_token(opt$b_name)
  up <- opt$up_label; dn <- opt$down_label

  # background = genes in ANY consensus locus (every locus came from a tested region)
  bg <- genes_in(gr_of(master), genes_gr)
  writeLines(bg, file.path(gldir, "background_universe.txt"))

  emit <- function(df, fname) writeLines(genes_in(gr_of(df), genes_gr), file.path(gldir, fname))

  sh <- master[master$category == "shared", , drop = FALSE]
  emit(sh,                                          sprintf("shared_all.txt"))
  emit(sh[sh$concordance == "concordant_up",   ],  sprintf("shared_concordant_%s.txt", safe_token(up)))
  emit(sh[sh$concordance == "concordant_down", ],  sprintf("shared_concordant_%s.txt", safe_token(dn)))
  emit(sh[sh$concordance == "discordant",      ],  sprintf("shared_discordant.txt"))

  ao <- master[master$category == "a_only", , drop = FALSE]
  emit(ao,                          sprintf("%s_only_all.txt", aT))
  emit(ao[ao$a_dir == "up",   ],    sprintf("%s_only_%s.txt", aT, safe_token(up)))
  emit(ao[ao$a_dir == "down", ],    sprintf("%s_only_%s.txt", aT, safe_token(dn)))

  bo <- master[master$category == "b_only", , drop = FALSE]
  emit(bo,                          sprintf("%s_only_all.txt", bT))
  emit(bo[bo$b_dir == "up",   ],    sprintf("%s_only_%s.txt", bT, safe_token(up)))
  emit(bo[bo$b_dir == "down", ],    sprintf("%s_only_%s.txt", bT, safe_token(dn)))

  # provenance (pathway script ignores files matching 'gene_assignment')
  loci_gr <- gr_of(master)
  if (length(loci_gr) && length(genes_gr)) {
    h <- findOverlaps(genes_gr, loci_gr)
    if (length(h)) {
      asg <- data.frame(chr = as.character(seqnames(loci_gr))[subjectHits(h)],
                        start = start(loci_gr)[subjectHits(h)], end = end(loci_gr)[subjectHits(h)],
                        category = master$category[subjectHits(h)], gene = genes_gr$gene[queryHits(h)],
                        stringsAsFactors = FALSE)
      write.csv(asg, file.path(gldir, "gene_assignment.csv"), row.names = FALSE)
    }
  }
  say(sprintf("[GENES] background=%d genes; per-category lists -> %s", length(bg), gldir))
  if (!length(bg)) say("[GENES] WARNING: 0 genes mapped — check chromosome naming (e.g. 'chr1' vs '1') matches your region coordinates.")
  # also attach a genes column to the master (semicolon-joined per locus)
  if (length(loci_gr) && length(genes_gr)) {
    h <- findOverlaps(genes_gr, loci_gr)
    gl <- tapply(genes_gr$gene[queryHits(h)], subjectHits(h), function(x) paste(sort(unique(x)), collapse = ";"))
    master$genes <- ""; master$genes[as.integer(names(gl))] <- gl
  }
  master
}

# =============================================================================
# MAIN
# =============================================================================
main <- function() {
  opt <- parse_cli()
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
  say(sprintf("=== csaw_overlap_venn | %s vs %s | overlap=%s %.0f%% | FDR<%.3g |logFC|>%.2g ===",
              opt$a_name, opt$b_name, opt$overlap_mode, 100 * opt$min_overlap, opt$fdr, opt$lfc))

  A <- read_region_csv(opt$a, opt$fdr, opt$lfc, opt$a_name)
  B <- read_region_csv(opt$b, opt$fdr, opt$lfc, opt$b_name)
  genes_gr <- resolve_genes(opt)

  grA_sig <- A$gr[A$gr$sig]; grB_sig <- B$gr[B$gr$sig]
  grA_up  <- grA_sig[grA_sig$dir == "up"];   grB_up  <- grB_sig[grB_sig$dir == "up"]
  grA_dn  <- grA_sig[grA_sig$dir == "down"]; grB_dn  <- grB_sig[grB_sig$dir == "down"]

  # ---- three independent clusterings -> three Venns ----
  cl_all <- cluster_two_sets(grA_sig, grB_sig, opt$min_overlap, opt$overlap_mode)
  cl_up  <- cluster_two_sets(grA_up,  grB_up,  opt$min_overlap, opt$overlap_mode)
  cl_dn  <- cluster_two_sets(grA_dn,  grB_dn,  opt$min_overlap, opt$overlap_mode)
  v_all <- venn_counts(cl_all); v_up <- venn_counts(cl_up); v_dn <- venn_counts(cl_dn)

  say(sprintf("[VENN all ] %s-only=%d  shared=%d  %s-only=%d  (Jaccard=%.3f)",
              opt$a_name, v_all["a_only"], v_all["shared"], opt$b_name, v_all["b_only"], jaccard(v_all)))
  say(sprintf("[VENN %-4s] %s-only=%d  shared=%d  %s-only=%d  (Jaccard=%.3f)",
              opt$up_label, opt$a_name, v_up["a_only"], v_up["shared"], opt$b_name, v_up["b_only"], jaccard(v_up)))
  say(sprintf("[VENN %-4s] %s-only=%d  shared=%d  %s-only=%d  (Jaccard=%.3f)",
              opt$down_label, opt$a_name, v_dn["a_only"], v_dn["shared"], opt$b_name, v_dn["b_only"], jaccard(v_dn)))

  save_fig("venn_all",  opt$outdir, draw_venn2(v_all, opt$a_name, opt$b_name,
           "All significant regions", sprintf("consensus loci | %.0f%% %s overlap | Jaccard=%.3f", 100*opt$min_overlap, opt$overlap_mode, jaccard(v_all))), w = 6, h = 5)
  save_fig("venn_up",   opt$outdir, draw_venn2(v_up, opt$a_name, opt$b_name,
           sprintf("%s regions (logFC > 0)", tools::toTitleCase(opt$up_label)), sprintf("Jaccard=%.3f", jaccard(v_up))), w = 6, h = 5)
  save_fig("venn_down", opt$outdir, draw_venn2(v_dn, opt$a_name, opt$b_name,
           sprintf("%s regions (logFC < 0)", tools::toTitleCase(opt$down_label)), sprintf("Jaccard=%.3f", jaccard(v_dn))), w = 6, h = 5)

  # ---- master table (all-clustering) + downstream analyses ----
  master <- build_master(grA_sig, grB_sig, A$gr, B$gr, cl_all)
  if (!is.null(genes_gr) && nrow(master)) master <- write_gene_lists(master, genes_gr, opt, opt$outdir)
  write.table(master, file.path(opt$outdir, "consensus_loci.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  conc <- concordance_bar(master, opt, opt$outdir)
  scat <- scatter_fig(master, opt, opt$outdir)
  enr  <- overlap_enrichment(A, B, opt, opt$outdir)

  # ---- summary report ----
  summ <- file.path(opt$outdir, "overlap_summary.txt")
  L <- c(strrep("=", 70), sprintf("csaw overlap summary | %s vs %s", opt$a_name, opt$b_name), strrep("=", 70), "",
         sprintf("Overlap rule: %s, min fraction %.2f", opt$overlap_mode, opt$min_overlap),
         sprintf("Significance:  FDR < %.3g and |log2FC| > %.2g", opt$fdr, opt$lfc), "",
         "Consensus-locus Venn counts (a_only / shared / b_only | Jaccard):",
         sprintf("  all : %d / %d / %d | %.3f", v_all["a_only"], v_all["shared"], v_all["b_only"], jaccard(v_all)),
         sprintf("  %-4s: %d / %d / %d | %.3f", opt$up_label,   v_up["a_only"], v_up["shared"], v_up["b_only"], jaccard(v_up)),
         sprintf("  %-4s: %d / %d / %d | %.3f", opt$down_label, v_dn["a_only"], v_dn["shared"], v_dn["b_only"], jaccard(v_dn)), "")
  if (!is.null(conc)) L <- c(L, "Direction concordance at shared loci (all-clustering):",
                             sprintf("  concordant %s   : %d", opt$up_label, conc["concordant_up"]),
                             sprintf("  concordant %s   : %d", opt$down_label, conc["concordant_down"]),
                             sprintf("  discordant       : %d", conc["discordant"]), "")
  if (!is.null(scat)) L <- c(L, sprintf("Shared-locus logFC correlation: Spearman=%.3f, Pearson=%.3f",
                                        scat$rho, scat$pearson), "")
  if (!is.null(enr))  L <- c(L, sprintf("Overlap enrichment (Fisher on co-tested loci): OR=%.2f, p=%.3g",
                                        enr$odds_ratio, enr$p), "")
  L <- c(L, "Files:",
         "  venn_all/up/down.{pdf,png}        - the three Venns",
         "  logfc_concordance_scatter.{pdf,png} - A vs B logFC per locus",
         "  direction_concordance_bar.{pdf,png} - concordant/discordant shared loci",
         "  consensus_loci.tsv                - master table (coords, logFC, category, concordance, genes)",
         if (!is.null(enr)) "  overlap_enrichment.tsv            - Fisher 2x2" else NULL,
         if (!is.null(genes_gr)) "  genelists/                        - per-category gene lists + background_universe.txt" else NULL,
         "", "--- run log ---", log_lines$txt)
  writeLines(L, summ)
  say(sprintf("[DONE] summary -> %s", summ))
}

main()