#!/usr/bin/env Rscript
# =============================================================================
# Q1: PLS-DA — mixOmics (Nature/Science standard for supervised multiomics)
# Reference: Rohart et al., PLOS Computational Biology 2017
#
# Usage:
#   Rscript scripts/Q1b_PLSDA.R \
#     outputs/figures/Q1_RF_SHAP/Q1_switch_classified.txt \
#     outputs/figures/Q1_PLSDA
# =============================================================================

suppressPackageStartupMessages({
  library(mixOmics)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
args       <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else
  "outputs/figures/Q1_RF_SHAP/Q1_switch_classified.txt"
outdir     <- if (length(args) >= 2) args[2] else
  "outputs/figures/Q1_PLSDA"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ── Publication theme ─────────────────────────────────────────────────────────
theme_pub <- theme_classic(base_size = 11) +
  theme(
    text             = element_text(family = "sans"),
    axis.text        = element_text(size = 10, color = "black"),
    axis.title       = element_text(size = 12),
    plot.title       = element_text(size = 13, face = "bold", hjust = 0),
    legend.text      = element_text(size = 9),
    legend.title     = element_text(size = 9, face = "bold"),
    panel.grid       = element_blank(),
    axis.ticks       = element_line(color = "black", linewidth = 0.5),
    axis.line        = element_line(color = "black", linewidth = 0.5)
  )

switch_colors <- c(
  "B_to_A"   = "#C0392B",
  "A_to_B"   = "#1A5276",
  "Stable_A" = "#85C1E9",
  "Stable_B" = "#D5D8DC"
)

# ── Load data ──────────────────────────────────────────────────────────────────
cat("[INFO] Loading classified gene table...\n")
df <- read_tsv(input_file, show_col_types = FALSE)
cat(sprintf("[INFO] %d genes loaded\n", nrow(df)))

# Detect mark columns (not gene_name, not switch_cat, not PC1 columns)
mark_cols <- names(df)[names(df) %in%
               c("H3K27me3", "H3K27Ac", "H3K4me1", "H3K9me3")]

if (length(mark_cols) < 2) {
  stop("[ERROR] Need at least 2 mark columns. Check input file.")
}
cat(sprintf("[INFO] Mark columns: %s\n", paste(mark_cols, collapse = ", ")))

# Keep four main switch categories
df <- df %>%
  filter(switch_cat %in% c("B_to_A", "A_to_B", "Stable_A", "Stable_B")) %>%
  filter(complete.cases(across(all_of(mark_cols))))

cat(sprintf("[INFO] Genes for PLS-DA: %d\n", nrow(df)))
print(table(df$switch_cat))

# ── Build X and Y ─────────────────────────────────────────────────────────────
X <- as.matrix(df[, mark_cols])
Y <- factor(df$switch_cat,
            levels = c("B_to_A", "A_to_B", "Stable_A", "Stable_B"))

# ── PLS-DA ────────────────────────────────────────────────────────────────────
set.seed(42)
plsda_res <- plsda(X, Y, ncomp = 3)

# ── Cross-validation (5-fold, 50 repeats — Nature standard) ───────────────────
cat("[INFO] Running 5-fold CV, 50 repeats...\n")
set.seed(42)
perf_res <- perf(
  plsda_res,
  validation  = "Mfold",
  folds       = 5,
  nrepeat     = 50,
  progressBar = FALSE
)

# Save performance plot
pdf(file.path(outdir, "Fig1_PLSDA_performance.pdf"), width = 6, height = 5)
plot(perf_res,
     col             = color.mixo(1:3),
     sd              = TRUE,
     legend.position = "horizontal")
title(
  main     = "PLS-DA Classification Error Rate\n(5-fold CV, 50 repeats)",
  cex.main = 1.1,
  font.main = 2
)
dev.off()
cat("[INFO] Saved Fig1_PLSDA_performance.pdf\n")

# ── VIP scores ────────────────────────────────────────────────────────────────
# VIP > 1 = significant contribution (Chong & Jun 2005 — standard threshold)
vip_mat <- vip(plsda_res)
vip_df  <- as.data.frame(vip_mat) %>%
  rownames_to_column("mark") %>%
  rename(VIP_comp1 = `comp 1`) %>%
  arrange(desc(VIP_comp1)) %>%
  mutate(significant = VIP_comp1 > 1)

write_tsv(vip_df, file.path(outdir, "Q1_PLSDA_VIP_scores.txt"))
cat("[INFO] VIP scores:\n")
print(vip_df)

# VIP bar plot
p_vip <- ggplot(vip_df,
                aes(x = reorder(mark, VIP_comp1),
                    y = VIP_comp1,
                    fill = significant)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "#C0392B", linewidth = 0.8) +
  scale_fill_manual(
    values = c("TRUE" = "#C0392B", "FALSE" = "#85C1E9"),
    labels = c("TRUE" = "VIP > 1", "FALSE" = "VIP \u2264 1")
  ) +
  coord_flip() +
  labs(
    x     = "",
    y     = "VIP Score (Component 1)",
    fill  = "",
    title = "Variable Importance in Projection\nParental Marks \u2192 Compartment Switch Prediction"
  ) +
  annotate("text", x = 0.6, y = 1.05,
           label = "VIP = 1.0 (significance threshold)",
           color = "#C0392B", size = 3, hjust = 0) +
  theme_pub +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "Fig1_PLSDA_VIP.pdf"),
       p_vip, width = 5.5, height = 4.5, dpi = 600)
cat("[INFO] Saved Fig1_PLSDA_VIP.pdf\n")

# ── Score plot (gene/sample separation in PLS-DA space) ──────────────────────
pdf(file.path(outdir, "Fig1_PLSDA_scores.pdf"), width = 6, height = 5.5)
plotIndiv(
  plsda_res,
  comp          = c(1, 2),
  group         = Y,
  ind.names     = FALSE,
  ellipse       = TRUE,        # 95% confidence ellipses — standard
  legend        = TRUE,
  title         = "PLS-DA: Parental Marks vs Compartment Switch",
  col.per.group = unname(switch_colors[levels(Y)]),
  pch           = 19,
  alpha         = 0.6,
  point.lwd     = 0.3
)
dev.off()
cat("[INFO] Saved Fig1_PLSDA_scores.pdf\n")

# ── Loading plot ───────────────────────────────────────────────────────────────
pdf(file.path(outdir, "Fig1_PLSDA_loadings.pdf"), width = 5.5, height = 4.5)
plotLoadings(
  plsda_res,
  comp    = 1,
  method  = "mean",
  contrib = "max",
  title   = "PLS-DA Loadings: Component 1"
)
dev.off()
cat("[INFO] Saved Fig1_PLSDA_loadings.pdf\n")

cat(sprintf("\n[INFO] Q1 PLS-DA complete -> %s\n", normalizePath(outdir)))