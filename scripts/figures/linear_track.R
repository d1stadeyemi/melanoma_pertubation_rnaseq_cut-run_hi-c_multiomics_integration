#!/usr/bin/env Rscript
# =============================================================================
# Linear stacked genome-browser tracks for EXAMPLE COMPARTMENT-SWITCH LOCI
# A375 melanoma — Parental vs Bottom10
#
# WHY THIS PLOT (vs circos) ---------------------------------------------------
#   The scientific question is VERTICAL co-variation: at a switch locus, do the
#   parental marks, the PC1 compartment flip, and expression move together?
#   A linear stacked-track view answers that directly — you read straight down
#   a genomic position and see every assay at once. Circos bends tracks onto
#   different radii and destroys that read. This is the figure that shows the
#   MECHANISM at a locus (e.g. high parental H3K4me1 -> B->A flip -> expr up).
#
#   Marks are kept at NATIVE 5kb here (fine mark structure is the whole point).
#   PC1 is 250kb (its native resolution). They share the x-axis (bp), so the
#   resolution difference is honest and visible, not silently mixed.
#
# TRACKS (top -> bottom), per locus -------------------------------------------
#   1 Parental H3K27ac      (5kb)   active/enhancer acetylation
#   2 Parental H3K4me1       (5kb)   enhancer priming  (drives B->A in our SHAP)
#   3 Parental H3K27me3      (5kb)   facultative repression (the EZH2i target)
#   4 Parental H3K9me3       (5kb)   constitutive repression (drives A->B)
#   5 Parental PC1          (250kb)  baseline compartment (red A / blue B)
#   6 Bottom10 PC1          (250kb)  post-migration compartment
#   7 Expression log2FC     (gene)   B10/parental, up=red down=blue
#   Switch interval shaded across all tracks.
#
# INPUTS ----------------------------------------------------------------------
#   --cutrun  cut_run_signal_table.txt  (5kb; parental_<mark>_mean columns)
#   --merged  merged_hic_rnaseq_cutrun_parent_bottom10.txt (250kb PC1 + expr)
#   --loci    optional CSV: chr,start,end,name,switch  (regions to plot).
#             If omitted, the script AUTO-PICKS the top example switch loci
#             (largest |PC1 change|) — a few B->A and a few A->B.
#
# USAGE -----------------------------------------------------------------------
#   Rscript scripts/linear_tracks.R \
#       --cutrun data/processed/cut_run_signal_table.txt \
#       --merged data/processed/merged_hic_rnaseq_cutrun_parent_bottom10.txt \
#       --outdir results/final/figures/linear_tracks
#
#   # plot specific loci instead of auto-picking:
#   Rscript scripts/linear_tracks.R ... --loci my_regions.csv
#
# OUTPUT ----------------------------------------------------------------------
#   <outdir>/linear_tracks_<name>.pdf / .png   one figure per locus
#   <outdir>/auto_selected_loci.csv            the loci chosen (if auto)
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(stringr)
  library(GenomicRanges)
  library(karyoploteR)
})

# ── CLI ──────────────────────────────────────────────────────────────────────
opt <- parse_args(OptionParser(option_list = list(
  make_option("--cutrun", type = "character"),
  make_option("--merged", type = "character"),
  make_option("--outdir", type = "character",
              default = "results/final/figures/linear_tracks"),
  make_option("--loci",   type = "character", default = NULL,
              help = "Optional CSV chr,start,end,name,switch. Auto-pick if absent."),
  make_option("--genome", type = "character", default = "hg19"),
  make_option("--flank",  type = "double", default = 1.5,
              help = "Flank multiplier around each locus for context [default 1.5x width]"),
  make_option("--n_each", type = "integer", default = 3,
              help = "If auto-picking, how many B->A and A->B loci each [default 3].")
)))
if (is.null(opt$cutrun) || is.null(opt$merged))
  stop("--cutrun and --merged are required.")
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ── Load merged (PC1 + expression, 250kb) ────────────────────────────────────
cat("[INFO] Loading merged table...\n")
m <- read_tsv(opt$merged, show_col_types = FALSE) %>%
  filter(!is.na(coordinate)) %>%
  mutate(
    chr   = str_extract(coordinate, "^chr[0-9XY]+"),
    start = as.numeric(str_extract(coordinate, "(?<=:)[0-9]+")),
    end   = as.numeric(str_extract(coordinate, "(?<=-)[0-9]+$")),
    Parental_orig = as.numeric(Parental_orig),
    Bottom10_orig = as.numeric(Bottom10_orig),
    log2FC        = suppressWarnings(as.numeric(log2FC))
  ) %>%
  filter(!is.na(chr), !is.na(start), !is.na(end))

standard_chr <- paste0("chr", c(1:22, "X", "Y"))
m <- m %>% filter(chr %in% standard_chr)

# Per-bin PC1 (one value per 250kb bin; genes may repeat a bin)
pc1 <- m %>%
  group_by(chr, start, end) %>%
  summarise(PC1_par = mean(Parental_orig, na.rm = TRUE),
            PC1_b10 = mean(Bottom10_orig, na.rm = TRUE),
            .groups = "drop") %>%
  filter(!is.na(PC1_par), !is.na(PC1_b10))

# Gene-level expression points (keep one row per gene if a gene col exists)
gene_col <- intersect(c("gene", "gene_name", "symbol", "Gene"), colnames(m))
expr <- m %>%
  filter(!is.na(log2FC)) %>%
  select(chr, start, end, log2FC)

# ── Classify switches (same logic as your other scripts) ─────────────────────
classify_switch <- function(both, sw, shift) {
  mapply(function(b, s, h) {
    b <- trimws(as.character(b)); s <- trimws(as.character(s)); h <- trimws(as.character(h))
    if (!b %in% c("-", "", "NA", "nan")) {
      if (grepl("toA|BtoA", b)) return("B_to_A")
      if (grepl("toB|AtoB", b)) return("A_to_B")
    }
    if (!s %in% c("-", "", "NA", "nan")) {
      if (s == "BtoA") return("B_to_A"); if (s == "AtoB") return("A_to_B")
    }
    if (!h %in% c("-", "", "NA", "nan")) {
      if (h == "TowardA") return("B_to_A"); if (h == "TowardB") return("A_to_B")
    }
    "Stable"
  }, both, sw, shift)
}
m$switch_cat <- classify_switch(m$`P->B10_Both_original`,
                                m$`P->B10_Switch_Orig`,
                                m$`P->B10_Shift_Orig`)

# ── Load CUT&RUN (5kb marks, parental baseline) ──────────────────────────────
cat("[INFO] Loading CUT&RUN 5kb signal...\n")
cr <- read_tsv(opt$cutrun, show_col_types = FALSE)
# parse coords from chr/start/end if present, else from region_id
if (all(c("chr", "start", "end") %in% colnames(cr))) {
  cr <- cr %>% mutate(chr = as.character(chr),
                      start = as.numeric(start), end = as.numeric(end))
} else {
  cr <- cr %>% mutate(
    chr   = str_extract(region_id, "^chr[0-9XY]+"),
    start = as.numeric(str_extract(region_id, "(?<=_)[0-9]+(?=_)")),
    end   = as.numeric(str_extract(region_id, "[0-9]+$"))
  )
}
mark_cols <- c(par_h3k27ac  = "parental_h3k27ac_mean",
               par_h3k4me1  = "parental_h3k4me1_mean",
               par_h3k27me3 = "parental_h3k27me3_mean",
               par_h3k9me3  = "parental_h3k9me3_mean")
missing <- setdiff(unname(mark_cols), colnames(cr))
if (length(missing)) stop("Missing CUT&RUN columns: ", paste(missing, collapse=", "))
for (c0 in unname(mark_cols)) cr[[c0]] <- suppressWarnings(as.numeric(cr[[c0]]))

# ── Pick loci ────────────────────────────────────────────────────────────────
if (!is.null(opt$loci)) {
  loci <- read_csv(opt$loci, show_col_types = FALSE)
  if (!"name" %in% colnames(loci))
    loci$name <- sprintf("%s_%d_%d", loci$chr, loci$start, loci$end)
  cat(sprintf("[INFO] Using %d user-supplied loci.\n", nrow(loci)))
} else {
  cat("[INFO] Auto-picking example switch loci by |PC1 change|...\n")
  cand <- m %>%
    filter(switch_cat %in% c("B_to_A", "A_to_B")) %>%
    group_by(chr, start, end, switch_cat) %>%
    summarise(PC1_par = mean(Parental_orig, na.rm = TRUE),
              PC1_b10 = mean(Bottom10_orig, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(dPC1 = PC1_b10 - PC1_par, adp = abs(dPC1))
  pick <- bind_rows(
    cand %>% filter(switch_cat == "B_to_A") %>% arrange(desc(adp)) %>% head(opt$n_each),
    cand %>% filter(switch_cat == "A_to_B") %>% arrange(desc(adp)) %>% head(opt$n_each)
  ) %>%
    mutate(name = sprintf("%s_%s_%d", switch_cat, chr, start),
           switch = switch_cat)
  loci <- pick %>% select(chr, start, end, name, switch)
  write_csv(loci, file.path(opt$outdir, "auto_selected_loci.csv"))
  cat(sprintf("[INFO] Auto-picked %d loci -> auto_selected_loci.csv\n", nrow(loci)))
}

# ── Helper: build a GRanges of one mark within a window ──────────────────────
mark_gr <- function(chrom, from, to, col) {
  sel <- cr %>% filter(chr == chrom, end >= from, start <= to)
  if (nrow(sel) == 0) return(NULL)
  GRanges(seqnames = sel$chr,
          ranges   = IRanges(start = sel$start, end = sel$end),
          score    = sel[[col]])
}
pc1_gr <- function(chrom, from, to, which) {
  sel <- pc1 %>% filter(chr == chrom, end >= from, start <= to)
  if (nrow(sel) == 0) return(NULL)
  GRanges(seqnames = sel$chr, ranges = IRanges(sel$start, sel$end),
          score = sel[[which]])
}
expr_gr <- function(chrom, from, to) {
  sel <- expr %>% filter(chr == chrom, end >= from, start <= to)
  if (nrow(sel) == 0) return(NULL)
  GRanges(seqnames = sel$chr, ranges = IRanges(sel$start, sel$end),
          score = sel$log2FC)
}

# colours
COL <- list(h3k27ac="#27AE60", h3k4me1="#8E44AD", h3k27me3="#E67E22",
            h3k9me3="#2C3E50", A="#C0392B", B="#1A5276",
            up="#C0392B", down="#1A5276", shade="#F7E59A")

# signed two-colour PC1 plotter (positive=A red, negative=B blue)
plot_signed <- function(kp, gr, r0, r1, poscol, negcol, label) {
  if (is.null(gr)) return(invisible())
  rng <- max(abs(range(gr$score, na.rm = TRUE)), 1e-6)
  pos <- gr; pos$score[pos$score < 0] <- 0
  neg <- gr; neg$score[neg$score > 0] <- 0
  kpBars(kp, pos, y1 = pos$score/rng, y0 = 0, r0=r0, r1=r1,
         ymin=-1, ymax=1, col=poscol, border=NA)
  kpBars(kp, neg, y1 = neg$score/rng, y0 = 0, r0=r0, r1=r1,
         ymin=-1, ymax=1, col=negcol, border=NA)
  kpAbline(kp, h=0, r0=r0, r1=r1, ymin=-1, ymax=1, col="black", lwd=0.6)
  kpAddLabels(kp, labels=label, r0=r0, r1=r1, cex=0.7, label.margin=0.02)
}
# single-colour density area for a mark
plot_mark <- function(kp, gr, r0, r1, col, label) {
  if (is.null(gr)) return(invisible())
  mx <- max(gr$score, na.rm = TRUE); if (!is.finite(mx) || mx<=0) mx <- 1
  kpArea(kp, gr, y = gr$score/mx, r0=r0, r1=r1, ymin=0, ymax=1,
         col=col, border=col, base.y=0)
  kpAddLabels(kp, labels=label, r0=r0, r1=r1, cex=0.7, label.margin=0.02)
}

# ── Draw one locus ───────────────────────────────────────────────────────────
draw_locus <- function(row) {
  chrom <- row$chr; w <- row$end - row$start
  from <- max(1, round(row$start - opt$flank * w))
  to   <- round(row$end + opt$flank * w)
  region <- toGRanges(data.frame(chrom, from, to))

  pp <- getDefaultPlotParams(plot.type = 1)
  pp$leftmargin <- 0.18; pp$data1height <- 700

  draw <- function() {
    kp <- plotKaryotype(genome = opt$genome, zoom = region, plot.type = 1,
                        plot.params = pp,
                        main = sprintf("%s   (%s:%s-%s)  %s",
                                       row$name, chrom,
                                       format(row$start, big.mark=","),
                                       format(row$end, big.mark=","),
                                       row$switch))
    kpAddBaseNumbers(kp, tick.dist = max(50000, round(w/2)), cex = 0.5)
    # shade the switch interval across the whole panel
    kpRect(kp, chr=chrom, x0=row$start, x1=row$end, y0=0, y1=1,
           col=COL$shade, border=NA, r0=0, r1=1)

    # track bands (top->bottom)
    plot_mark(kp, mark_gr(chrom,from,to,mark_cols["par_h3k27ac"]),  0.88,1.00, COL$h3k27ac,  "Par H3K27ac")
    plot_mark(kp, mark_gr(chrom,from,to,mark_cols["par_h3k4me1"]),  0.74,0.86, COL$h3k4me1,  "Par H3K4me1")
    plot_mark(kp, mark_gr(chrom,from,to,mark_cols["par_h3k27me3"]), 0.60,0.72, COL$h3k27me3, "Par H3K27me3")
    plot_mark(kp, mark_gr(chrom,from,to,mark_cols["par_h3k9me3"]),  0.46,0.58, COL$h3k9me3,  "Par H3K9me3")
    plot_signed(kp, pc1_gr(chrom,from,to,"PC1_par"), 0.30,0.44, COL$A, COL$B, "Par PC1")
    plot_signed(kp, pc1_gr(chrom,from,to,"PC1_b10"), 0.14,0.28, COL$A, COL$B, "B10 PC1")
    plot_signed(kp, expr_gr(chrom,from,to),          0.00,0.12, COL$up, COL$down, "log2FC")
  }

  base <- file.path(opt$outdir, paste0("linear_tracks_", row$name))
  pdf(paste0(base, ".pdf"), width = 9, height = 7); draw(); dev.off()
  png(paste0(base, ".png"), width = 9, height = 7, units="in", res = 300); draw(); dev.off()
  cat(sprintf("  Saved -> %s.[pdf|png]\n", base))
}

for (i in seq_len(nrow(loci))) draw_locus(loci[i, ])
cat(sprintf("[INFO] Done. %d locus figures -> %s\n", nrow(loci), opt$outdir))