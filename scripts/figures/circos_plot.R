#!/usr/bin/env Rscript
# =============================================================================
# Circos WHOLE-GENOME OVERVIEW — A375 melanoma Parental vs Bottom10
# (a "wow" overview figure; the mechanistic detail lives in linear_tracks.R)
#
# WHAT CHANGED vs the earlier version -----------------------------------------
#   1. Marks come from the 5kb CUT&RUN table and are drawn as DENSE LINE/AREA
#      rings (not ~570k rectangles, which would choke the renderer and bloat
#      the PDF). Dense rings = the look you want, without the cost.
#   2. PC1 is shown as PAIRED rings (Parental then Bottom10, adjacent) so the
#      compartment flip is read WITHIN one figure — no eyeballing two plots.
#   3. The meaningless self-links (region1==region2) are REPLACED by a signed
#      Δ-PC1 ring (Bottom10 - Parental): red = gained-A (B->A), blue = gained-B
#      (A->B). This actually encodes switching, instead of decorative arcs.
#   4. PDF and PNG are produced by ONE draw_circos() function (no copy-paste
#      drift between devices).
#
# RESOLUTION NOTE -------------------------------------------------------------
#   Marks are 5kb (dense, for aesthetics); PC1 is 250kb (native). On a circos
#   overview the viewer does not co-read a mark value against a PC1 value at a
#   single radius, so mixing resolutions is fine HERE. For quantitative
#   locus-level co-interpretation use linear_tracks.R instead.
#   For very dense marks the script optionally bins them to --mark_display_bp
#   purely for drawing speed (default 10kb); set =5000 to keep native 5kb.
#
# INPUTS / USAGE --------------------------------------------------------------
#   Rscript scripts/circos_overview.R \
#       --cutrun data/processed/cut_run_signal_table.txt \
#       --merged data/processed/merged_hic_rnaseq_cutrun_parent_bottom10.txt \
#       --outdir results/final/figures/circos
#
# OUTPUT: <outdir>/circos_overview_A375.pdf / .png
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(circlize)
  library(dplyr)
  library(readr)
  library(stringr)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--cutrun", type="character"),
  make_option("--merged", type="character"),
  make_option("--outdir", type="character",
              default="results/final/figures/circos"),
  make_option("--mark_display_bp", type="integer", default=10000L,
              help="Bin width (bp) for DISPLAY of 5kb marks; 5000 keeps native [default 10000]")
)))
if (is.null(opt$cutrun) || is.null(opt$merged))
  stop("--cutrun and --merged are required.")
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

standard_chr <- paste0("chr", c(1:22, "X", "Y"))

# ── Hi-C PC1 (250kb) ──────────────────────────────────────────────────────────
cat("[INFO] Loading merged (PC1) ...\n")
m <- read_tsv(opt$merged, show_col_types = FALSE) %>%
  filter(!is.na(coordinate)) %>%
  mutate(chr=str_extract(coordinate,"^chr[0-9XY]+"),
         start=as.numeric(str_extract(coordinate,"(?<=:)[0-9]+")),
         end=as.numeric(str_extract(coordinate,"(?<=-)[0-9]+$")),
         Parental_orig=as.numeric(Parental_orig),
         Bottom10_orig=as.numeric(Bottom10_orig)) %>%
  filter(chr %in% standard_chr, !is.na(start), !is.na(end))

pc1 <- m %>% group_by(chr,start,end) %>%
  summarise(PC1_par=mean(Parental_orig,na.rm=TRUE),
            PC1_b10=mean(Bottom10_orig,na.rm=TRUE), .groups="drop") %>%
  filter(!is.na(PC1_par), !is.na(PC1_b10)) %>%
  mutate(dPC1 = PC1_b10 - PC1_par)

# ── CUT&RUN marks (5kb) -> display bins ──────────────────────────────────────
cat("[INFO] Loading CUT&RUN 5kb signal ...\n")
cr <- read_tsv(opt$cutrun, show_col_types = FALSE)
if (!all(c("chr","start","end") %in% colnames(cr))) {
  cr <- cr %>% mutate(chr=str_extract(region_id,"^chr[0-9XY]+"),
                      start=as.numeric(str_extract(region_id,"(?<=_)[0-9]+(?=_)")),
                      end=as.numeric(str_extract(region_id,"[0-9]+$")))
}
cr <- cr %>% filter(chr %in% standard_chr)
mark_cols <- c(h3k27ac="parental_h3k27ac_mean", h3k4me1="parental_h3k4me1_mean",
               h3k27me3="parental_h3k27me3_mean", h3k9me3="parental_h3k9me3_mean")
for (c0 in unname(mark_cols)) cr[[c0]] <- suppressWarnings(as.numeric(cr[[c0]]))

# Aggregate 5kb -> display bins (mean) for speed; keeps the dense look.
binw <- opt$mark_display_bp
if (binw <= 0) {
  cat("[INFO] mark_display_bp<=0 -> using marks as-is (no re-binning).\n")
  marks_disp <- cr %>% select(chr, start, end, all_of(unname(mark_cols)))
} else {
  cat(sprintf("[INFO] Aggregating marks to %d bp for display...\n", binw))
  marks_disp <- cr %>%
    mutate(bin = floor(start / binw) * binw) %>%
    group_by(chr, bin) %>%
    summarise(across(all_of(unname(mark_cols)), ~mean(.x, na.rm = TRUE)),
              .groups = "drop") %>%
    mutate(start = bin + 1, end = bin + binw) %>%
    select(chr, start, end, everything(), -bin)
}

# ── hg19 sizes ────────────────────────────────────────────────────────────────
chr_sizes <- data.frame(
  chr=paste0("chr",c(1:22,"X","Y")),
  length=c(249250621,243199373,198022430,191154276,180915260,171115067,
           159138663,146364022,141213431,135534747,135006516,133851895,
           115169878,107349540,102531392,90354753,81195210,78077248,
           59128983,63025520,48129895,51304566,155270560,59373566),
  stringsAsFactors=FALSE) %>%
  filter(chr %in% unique(pc1$chr))

# ── colours / helpers ─────────────────────────────────────────────────────────
COL <- list(h3k27ac="#27AE60", h3k4me1="#8E44AD", h3k27me3="#E67E22",
            h3k9me3="#2C3E50", A="#C0392B", B="#1A5276")
sel <- function(df, cols) df %>% select(all_of(c("chr","start","end",cols))) %>% as.data.frame()

# signed PC1 / dPC1 ring (pos=A red, neg=B blue)
panel_signed <- function(region, value, poscol, negcol, ...) {
  v <- value[[1]]
  v[!is.finite(v)] <- 0                      # guard NA/Inf bins
  circos.genomicRect(region, value, ytop = pmax(v, 0), ybottom = pmin(v, 0),
                     col = ifelse(v > 0, poscol, negcol), border = NA, ...)
  circos.lines(CELL_META$cell.xlim, c(0, 0), col = "black", lwd = 0.4)
}
# single-colour mark line+area ring
panel_mark <- function(region, value, col, ...) {
  x <- (region[[1]] + region[[2]]) / 2     # region = (start, end), 2 cols
  y <- value[[1]]
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(invisible())
  circos.lines(x[ok], y[ok], area = TRUE, col = paste0(col, "55"),
               border = col, lwd = 0.4)
}

# ── ONE draw function (used by both devices) ─────────────────────────────────
draw_circos <- function() {
  circos.par(start.degree=90, gap.degree=1.5,
             track.margin=c(0.004,0.004), cell.padding=c(0,0,0,0),
             points.overflow.warning=FALSE)
  chr_init <- data.frame(chr=chr_sizes$chr, start=1, end=chr_sizes$length)
  circos.genomicInitialize(chr_init, plotType=c("axis","labels"),
                           axis.labels.cex=0.45, labels.cex=0.65)

  # mark rings (outer -> inner): the dense, decorative layers
  circos.genomicTrack(sel(marks_disp, mark_cols["h3k27ac"]),  track.height=0.075,
    bg.border=NA, panel.fun=function(r,v,...) panel_mark(r,v,COL$h3k27ac,...))
  circos.genomicTrack(sel(marks_disp, mark_cols["h3k4me1"]),  track.height=0.075,
    bg.border=NA, panel.fun=function(r,v,...) panel_mark(r,v,COL$h3k4me1,...))
  circos.genomicTrack(sel(marks_disp, mark_cols["h3k27me3"]), track.height=0.075,
    bg.border=NA, panel.fun=function(r,v,...) panel_mark(r,v,COL$h3k27me3,...))
  circos.genomicTrack(sel(marks_disp, mark_cols["h3k9me3"]),  track.height=0.075,
    bg.border=NA, panel.fun=function(r,v,...) panel_mark(r,v,COL$h3k9me3,...))

  # paired PC1 rings
  circos.genomicTrack(sel(pc1,"PC1_par"), track.height=0.085, bg.border=NA,
    bg.col="#FBFBFB",
    panel.fun=function(r,v,...) panel_signed(r,v,COL$A,COL$B,...))
  circos.genomicTrack(sel(pc1,"PC1_b10"), track.height=0.085, bg.border=NA,
    bg.col="#FBFBFB",
    panel.fun=function(r,v,...) panel_signed(r,v,COL$A,COL$B,...))

  # Δ-PC1 ring (replaces self-links): red = B->A gain, blue = A->B gain
  circos.genomicTrack(sel(pc1,"dPC1"), track.height=0.085, bg.border=NA,
    panel.fun=function(r,v,...) panel_signed(r,v,COL$A,COL$B,...))

  legend("center", bty="n", cex=0.62, x.intersp=0.5,
    legend=c("A compartment (PC1>0)","B compartment (PC1<0)",
             "ΔPC1>0 : B→A gain","ΔPC1<0 : A→B gain"),
    fill=c(COL$A,COL$B,COL$A,COL$B), border=NA)
  legend("topleft", inset=c(0.0,0.0), bty="n", cex=0.55,
    legend=c("Rings outer→inner:","H3K27ac · H3K4me1 · H3K27me3 · H3K9me3 (parental, 5kb)",
             "Parental PC1 · Bottom10 PC1 · ΔPC1 (250kb)"))
  title(main="A375 melanoma: chromatin, 3D compartments and switching\nParental vs Bottom10 (constricted migration)",
        cex.main=1.05, font.main=2, line=-1)
  circos.clear()
}

# ── render both ───────────────────────────────────────────────────────────────
pdf_path <- file.path(opt$outdir, "circos_overview_A375.pdf")
pdf(pdf_path, width=12, height=12); draw_circos(); dev.off()
cat(sprintf("[INFO] PDF -> %s\n", pdf_path))

png_path <- file.path(opt$outdir, "circos_overview_A375.png")
png(png_path, width=12, height=12, units="in", res=600); draw_circos(); dev.off()
cat(sprintf("[INFO] PNG -> %s\n", png_path))