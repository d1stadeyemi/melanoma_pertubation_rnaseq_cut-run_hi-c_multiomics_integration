#!/usr/bin/env Rscript
# =============================================================================
# Circos Plot: Hi-C PC1 + RNA-seq + Compartment Switch Links
# A375 Melanoma — Parental vs Bottom10
#
# Tracks from outside to inside:
#   Ring 1 (outermost): Chromosome ideogram with labels
#   Ring 2: Gene expression — Parental baseMean (log10 scaled)
#   Ring 3: Gene expression — Bottom10 log2FC (colored Up/Down)
#   Ring 4: Hi-C PC1 — Parental_orig (red = A, blue = B)
#   Ring 5: Hi-C PC1 — Bottom10_orig (red = A, blue = B)
#   Inner links: Genes that switched compartment (B->A = red, A->B = blue)
#
# Usage:
#   Rscript scripts/circos_plot.R \
#     outputs/tables/merged_hic_rnaseq_cutrun_parent_bottom10.txt \
#     outputs/figures/circos/
# =============================================================================

suppressPackageStartupMessages({
  library(circlize)
  library(dplyr)
  library(readr)
  library(stringr)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
args        <- commandArgs(trailingOnly = TRUE)
merged_file <- if (length(args) >= 1) args[1] else
  "outputs/tables/merged_hic_rnaseq_cutrun_parent_bottom10.txt"
outdir      <- if (length(args) >= 2) args[2] else
  "outputs/figures/circos"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ── Load data ──────────────────────────────────────────────────────────────────
cat("[INFO] Loading merged table...\n")
df <- read_tsv(merged_file, show_col_types = FALSE)
cat(sprintf("[INFO] %d rows x %d columns\n", nrow(df), ncol(df)))

# ── Parse coordinate column into chr, start, end ──────────────────────────────
# coordinate format: "chr19:58750001-59000000"
# This is the Hi-C bin coordinate — 250kb resolution
cat("[INFO] Parsing coordinates...\n")
df <- df %>%
  filter(!is.na(coordinate)) %>%
  mutate(
    # Extract chromosome
    chr   = str_extract(coordinate, "^chr[0-9XY]+"),
    # Extract start position
    start = as.numeric(str_extract(coordinate, "(?<=:)[0-9]+")),
    # Extract end position
    end   = as.numeric(str_extract(coordinate, "(?<=-)[0-9]+$"))
  ) %>%
  filter(!is.na(chr) & !is.na(start) & !is.na(end))

cat(sprintf("[INFO] %d genes with valid coordinates\n", nrow(df)))

# ── Keep only standard chromosomes ────────────────────────────────────────────
standard_chr <- paste0("chr", c(1:22, "X", "Y"))
df <- df %>% filter(chr %in% standard_chr)

# ── Force numeric columns ─────────────────────────────────────────────────────
df <- df %>%
  mutate(
    Parental_orig  = as.numeric(Parental_orig),
    Bottom10_orig  = as.numeric(Bottom10_orig),
    log2FC         = as.numeric(log2FC),
    baseMean       = as.numeric(baseMean),
    padj           = as.numeric(padj)
  )

cat(sprintf("[INFO] Genes on standard chromosomes: %d\n", nrow(df)))

# ── Compartment switch classification ─────────────────────────────────────────
# Classify each gene into switch category for inner links
classify_switch <- function(both, sw, shift) {
  mapply(function(b, s, h) {
    b <- trimws(as.character(b))
    s <- trimws(as.character(s))
    h <- trimws(as.character(h))
    if (!b %in% c("-", "", "NA", "nan")) {
      if (grepl("toA|BtoA", b)) return("B_to_A")
      if (grepl("toB|AtoB", b)) return("A_to_B")
    }
    if (!s %in% c("-", "", "NA", "nan")) {
      if (s == "BtoA") return("B_to_A")
      if (s == "AtoB") return("A_to_B")
    }
    if (!h %in% c("-", "", "NA", "nan")) {
      if (h == "TowardA") return("B_to_A")
      if (h == "TowardB") return("A_to_B")
    }
    return("Stable")
  }, both, sw, shift)
}

df$switch_cat <- classify_switch(
  df$`P->B10_Both_original`,
  df$`P->B10_Switch_Orig`,
  df$`P->B10_Shift_Orig`
)

# ── Per-bin PC1 averages (one value per 250kb Hi-C bin) ───────────────────────
# Multiple genes can share the same Hi-C bin coordinate
# Average PC1 per bin for smooth circos tracks
bin_parental <- df %>%
  group_by(chr, start, end) %>%
  summarise(
    PC1_parental  = mean(Parental_orig, na.rm = TRUE),
    PC1_bottom10  = mean(Bottom10_orig, na.rm = TRUE),
    mean_log2FC   = mean(log2FC,        na.rm = TRUE),
    mean_baseMean = mean(log10(baseMean + 1), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(PC1_parental))

cat(sprintf("[INFO] Hi-C bins: %d\n", nrow(bin_parental)))

# ── Switched genes for inner links ────────────────────────────────────────────
# Only use the strongest switches (Both_original category)
# for visual clarity — too many links = unreadable plot
switched <- df %>%
  filter(switch_cat %in% c("B_to_A", "A_to_B")) %>%
  filter(!is.na(chr) & !is.na(start) & !is.na(end)) %>%
  # Deduplicate by bin — one link per Hi-C bin
  group_by(chr, start, end, switch_cat) %>%
  slice(1) %>%
  ungroup()

cat(sprintf("[INFO] Switched bins for links: %d\n", nrow(switched)))

# ── Human chromosome sizes (hg19) ─────────────────────────────────────────────
# These are the reference chromosome lengths for hg19
chr_sizes <- data.frame(
  chr    = paste0("chr", c(1:22, "X", "Y")),
  length = c(
    249250621, 243199373, 198022430, 191154276, 180915260,
    171115067, 159138663, 146364022, 141213431, 135534747,
    135006516, 133851895, 115169878, 107349540, 102531392,
    90354753,  81195210,  78077248,  59128983,  63025520,
    48129895,  51304566,  155270560, 59373566
  ),
  stringsAsFactors = FALSE
)

# Keep only chromosomes present in your data
chr_in_data <- unique(bin_parental$chr)
chr_sizes   <- chr_sizes %>% filter(chr %in% chr_in_data)

# ── Color functions ────────────────────────────────────────────────────────────
# PC1: positive = A compartment (red), negative = B compartment (blue)
pc1_color <- function(val) {
  ifelse(is.na(val), "#CCCCCC",
         ifelse(val > 0, "#C0392B", "#1A5276"))
}

# log2FC: up = red, down = blue, NS = grey
lfc_color <- function(val) {
  ifelse(is.na(val), "#CCCCCC",
         ifelse(val > 1,  "#C0392B",
                ifelse(val < -1, "#1A5276", "#AAAAAA")))
}

# ── OPEN PDF ──────────────────────────────────────────────────────────────────
pdf_path <- file.path(outdir, "circos_A375_Parental_vs_Bottom10.pdf")
pdf(pdf_path, width = 12, height = 12)

# ── Initialize circos ─────────────────────────────────────────────────────────
# circos.par sets global parameters before any tracks are drawn
circos.par(
  "start.degree"      = 90,       # start at top (12 o'clock)
  "gap.degree"        = 1.5,      # gap between chromosomes in degrees
  "track.margin"      = c(0.005, 0.005),
  "cell.padding"      = c(0, 0, 0, 0),
  "points.overflow.warning" = FALSE
)

# Initialize the genomic coordinate system using chromosome sizes
# This sets up one sector per chromosome
# circlize requires exactly: chr, start, end column names
# Build a proper genomic data frame with start=1 for each chromosome
chr_init <- data.frame(
  chr   = chr_sizes$chr,
  start = 1,
  end   = chr_sizes$length,
  stringsAsFactors = FALSE
)

circos.genomicInitialize(
  data            = chr_init,
  plotType        = c("axis", "labels"),
  axis.labels.cex = 0.5,
  labels.cex      = 0.7
)

# ── TRACK 1: Parental baseMean (gene expression, log10) ──────────────────────
# Bar chart — height = log10(baseMean + 1)
# Only genes with expression data
expr_data <- bin_parental %>%
  filter(!is.na(mean_baseMean)) %>%
  select(chr, start, end, mean_baseMean) %>%
  as.data.frame()

circos.genomicTrack(
  expr_data,
  track.height = 0.08,
  bg.border    = NA,
  bg.col       = "#F8F8F8",
  panel.fun    = function(region, value, ...) {
    circos.genomicRect(
      region, value,
      ytop        = value[[1]],
      ybottom     = 0,
      col         = "#27AE60",
      border      = NA,
      ...
    )
  }
)
# Add track label
text(-1.25, 0, "Parental\nExpression", cex = 0.55, adj = c(1, 0.5))

# ── TRACK 2: Bottom10 log2FC ──────────────────────────────────────────────────
# Bar chart colored by direction (up = red, down = blue)
lfc_data <- bin_parental %>%
  filter(!is.na(mean_log2FC)) %>%
  select(chr, start, end, mean_log2FC) %>%
  as.data.frame()

circos.genomicTrack(
  lfc_data,
  track.height = 0.08,
  bg.border    = NA,
  bg.col       = "#F8F8F8",
  panel.fun    = function(region, value, ...) {
    # Color each bar by direction
    col_vec <- lfc_color(value[[1]])
    circos.genomicRect(
      region, value,
      ytop    = pmax(value[[1]], 0),
      ybottom = pmin(value[[1]], 0),
      col     = col_vec,
      border  = NA,
      ...
    )
    # Baseline at 0
    circos.lines(
      CELL_META$cell.xlim,
      c(0, 0),
      col = "black", lwd = 0.5
    )
  }
)
text(-1.25, 0, "Bottom10\nlog2FC", cex = 0.55, adj = c(1, 0.5))

# ── TRACK 3: Parental PC1 ─────────────────────────────────────────────────────
# Bar chart: positive (A compartment) = red, negative (B) = blue
pc1_par_data <- bin_parental %>%
  filter(!is.na(PC1_parental)) %>%
  select(chr, start, end, PC1_parental) %>%
  as.data.frame()

circos.genomicTrack(
  pc1_par_data,
  track.height = 0.09,
  bg.border    = NA,
  bg.col       = "#F8F8F8",
  panel.fun    = function(region, value, ...) {
    col_vec <- pc1_color(value[[1]])
    circos.genomicRect(
      region, value,
      ytop    = pmax(value[[1]], 0),
      ybottom = pmin(value[[1]], 0),
      col     = col_vec,
      border  = NA,
      ...
    )
    circos.lines(
      CELL_META$cell.xlim,
      c(0, 0),
      col = "black", lwd = 0.5
    )
  }
)
text(-1.25, 0, "Parental\nPC1", cex = 0.55, adj = c(1, 0.5))

# ── TRACK 4: Bottom10 PC1 ────────────────────────────────────────────────────
pc1_bot_data <- bin_parental %>%
  filter(!is.na(PC1_bottom10)) %>%
  select(chr, start, end, PC1_bottom10) %>%
  as.data.frame()

circos.genomicTrack(
  pc1_bot_data,
  track.height = 0.09,
  bg.border    = NA,
  bg.col       = "#F8F8F8",
  panel.fun    = function(region, value, ...) {
    col_vec <- pc1_color(value[[1]])
    circos.genomicRect(
      region, value,
      ytop    = pmax(value[[1]], 0),
      ybottom = pmin(value[[1]], 0),
      col     = col_vec,
      border  = NA,
      ...
    )
    circos.lines(
      CELL_META$cell.xlim,
      c(0, 0),
      col = "black", lwd = 0.5
    )
  }
)
text(-1.25, 0, "Bottom10\nPC1", cex = 0.55, adj = c(1, 0.5))

# ── INNER LINKS: Compartment switch connections ───────────────────────────────
# Each switched bin is linked from its Parental position to its Bottom10 position
# Since both are at the same genomic location, we draw links between
# the Parental PC1 track and the Bottom10 PC1 track at the switch site.
# Visual effect: colored arcs inside the innermost ring highlighting switch loci.
#
# For genuine inter-chromosomal links you would need inter-chromosomal Hi-C data.
# Here we draw intra-position arcs connecting the same bin across conditions,
# which highlights the switched regions visually — same approach used in
# papers showing compartment changes with circos.

if (nrow(switched) > 0) {
  # Subsample for visual clarity — max 200 links
  set.seed(42)
  n_links <- min(200, nrow(switched))
  sw_plot <- switched[sample(nrow(switched), n_links), ]

  # B->A switches = red links, A->B = blue links
  link_colors <- ifelse(sw_plot$switch_cat == "B_to_A",
                        "#C0392B44",   # red, 44 = ~27% transparency in hex
                        "#1A527644")   # blue

  for (i in seq_len(nrow(sw_plot))) {
    tryCatch({
      circos.genomicLink(
        region1 = sw_plot[i, c("chr", "start", "end")],
        region2 = sw_plot[i, c("chr", "start", "end")],
        col     = link_colors[i],
        lwd     = 0.8,
        border  = NA
      )
    }, error = function(e) NULL)
  }
}

# ── LEGEND ────────────────────────────────────────────────────────────────────
# Draw legend in the center of the circos
legend(
  "center",
  legend = c(
    "A compartment (PC1 > 0)",
    "B compartment (PC1 < 0)",
    "Upregulated (log2FC > 1)",
    "Downregulated (log2FC < -1)",
    "B\u2192A switch",
    "A\u2192B switch"
  ),
  fill = c("#C0392B", "#1A5276", "#C0392B", "#1A5276",
           "#C0392B", "#1A5276"),
  border = NA,
  bty    = "n",
  cex    = 0.65,
  x.intersp = 0.5
)

# Track labels on the left
legend(
  "topleft",
  inset   = c(0.02, 0.02),
  legend  = c("A: Parental expression (log10 baseMean)",
              "B: Bottom10 log2FC",
              "C: Parental Hi-C PC1",
              "D: Bottom10 Hi-C PC1",
              "Links: compartment switches"),
  bty = "n",
  cex = 0.6
)

# Title
title(
  main = "A375 Melanoma: 3D Genome and Transcriptional\nRemodeling Upon Constricted Migration",
  cex.main = 1.1,
  font.main = 2,
  line = -1
)

# ── CLOSE ─────────────────────────────────────────────────────────────────────
circos.clear()
dev.off()

cat(sprintf("[INFO] Circos plot saved -> %s\n", pdf_path))

# Also save as PNG at 600 dpi for submission
png_path <- file.path(outdir, "circos_A375_Parental_vs_Bottom10.png")
png(png_path, width = 12, height = 12, units = "in", res = 600)

circos.par(
  "start.degree"      = 90,
  "gap.degree"        = 1.5,
  "track.margin"      = c(0.005, 0.005),
  "cell.padding"      = c(0, 0, 0, 0),
  "points.overflow.warning" = FALSE
)

# Re-run all tracks for PNG — identical code block
chr_init <- data.frame(
  chr   = chr_sizes$chr,
  start = 1,
  end   = chr_sizes$length,
  stringsAsFactors = FALSE
)

circos.genomicInitialize(
  data            = chr_init,
  plotType        = c("axis", "labels"),
  axis.labels.cex = 0.5,
  labels.cex      = 0.7
)

circos.genomicTrack(expr_data, track.height = 0.08, bg.border = NA,
  bg.col = "#F8F8F8",
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, ytop = value[[1]], ybottom = 0,
                       col = "#27AE60", border = NA, ...)
  })

circos.genomicTrack(lfc_data, track.height = 0.08, bg.border = NA,
  bg.col = "#F8F8F8",
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value,
                       ytop = pmax(value[[1]], 0),
                       ybottom = pmin(value[[1]], 0),
                       col = lfc_color(value[[1]]), border = NA, ...)
    circos.lines(CELL_META$cell.xlim, c(0, 0), col = "black", lwd = 0.5)
  })

circos.genomicTrack(pc1_par_data, track.height = 0.09, bg.border = NA,
  bg.col = "#F8F8F8",
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value,
                       ytop = pmax(value[[1]], 0),
                       ybottom = pmin(value[[1]], 0),
                       col = pc1_color(value[[1]]), border = NA, ...)
    circos.lines(CELL_META$cell.xlim, c(0, 0), col = "black", lwd = 0.5)
  })

circos.genomicTrack(pc1_bot_data, track.height = 0.09, bg.border = NA,
  bg.col = "#F8F8F8",
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value,
                       ytop = pmax(value[[1]], 0),
                       ybottom = pmin(value[[1]], 0),
                       col = pc1_color(value[[1]]), border = NA, ...)
    circos.lines(CELL_META$cell.xlim, c(0, 0), col = "black", lwd = 0.5)
  })

if (nrow(switched) > 0) {
  for (i in seq_len(nrow(sw_plot))) {
    tryCatch({
      circos.genomicLink(
        region1 = sw_plot[i, c("chr", "start", "end")],
        region2 = sw_plot[i, c("chr", "start", "end")],
        col     = link_colors[i], lwd = 0.8, border = NA
      )
    }, error = function(e) NULL)
  }
}

legend("center",
  legend = c("A compartment", "B compartment",
             "Upregulated", "Downregulated",
             "B\u2192A switch", "A\u2192B switch"),
  fill   = c("#C0392B", "#1A5276", "#C0392B",
             "#1A5276", "#C0392B", "#1A5276"),
  border = NA, bty = "n", cex = 0.65)

title(main = "A375 Melanoma: 3D Genome and Transcriptional\nRemodeling Upon Constricted Migration",
      cex.main = 1.1, font.main = 2, line = -1)

circos.clear()
dev.off()

cat(sprintf("[INFO] PNG saved -> %s\n", png_path))