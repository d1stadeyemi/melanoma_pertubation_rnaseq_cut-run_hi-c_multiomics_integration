# =============================================================================
# PLS-DA — mixOmics (Nature/Science standard for supervised multiomics)
# Reference: Rohart et al., PLOS Computational Biology 2017
#
# Revisions:
#   - log1p() transform of marks before plsda()
#   - Performance plot: drop max.dist; legend moved outside the panel
#   - VIP plot: numeric VIP values added at the end of each bar
#   - Score plot: stable categories faded; switchers on top; manual axis labels
#   - Loadings: use method = "median" (robust for skewed signal); add x-axis title
#
# Usage:
#   Rscript scripts/analysis/PLSDA.R \
#     results/exploratory/RF_SHAP/Switch_classified.txt \
#     results/exploratory/PLSDA
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
  "results/exploratory/RF_SHAP/Switch_classified.txt"
outdir     <- if (length(args) >= 2) args[2] else
  "results/exploratory/PLSDA"
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

# Detect all mark columns present (supports future CTCF column from Python flag)
candidate_marks <- c("H3K27me3", "H3K27Ac", "H3K4me1", "H3K9me3", "CTCF")
mark_cols       <- candidate_marks[candidate_marks %in% names(df)]

if (length(mark_cols) < 2) {
  stop("[ERROR] Need at least 2 mark columns. Check input file.")
}
cat(sprintf("[INFO] Mark columns: %s\n", paste(mark_cols, collapse = ", ")))

df <- df %>%
  filter(switch_cat %in% c("B_to_A", "A_to_B", "Stable_A", "Stable_B")) %>%
  filter(complete.cases(across(all_of(mark_cols))))

cat(sprintf("[INFO] Genes for PLS-DA: %d\n", nrow(df)))
print(table(df$switch_cat))

# ── Build X (log1p transform!) and Y ──────────────────────────────────────────
# Log1p compresses the H3K9me3 outliers and makes the supervised projection
# reflect the bulk of the data, not the extreme tail. This change is the single
# most impactful fix for the previously uninterpretable score plot.
X_raw <- as.matrix(df[, mark_cols])
X     <- log1p(X_raw)
Y     <- factor(df$switch_cat,
                levels = c("B_to_A", "A_to_B", "Stable_A", "Stable_B"))

# ── PLS-DA ────────────────────────────────────────────────────────────────────
set.seed(42)
plsda_res <- plsda(X, Y, ncomp = 3)

# ── Cross-validation (5-fold, 50 repeats) ─────────────────────────────────────
cat("[INFO] Running 5-fold CV, 50 repeats...\n")
set.seed(42)
perf_res <- perf(
  plsda_res,
  validation  = "Mfold",
  folds       = 5,
  nrepeat     = 50,
  progressBar = FALSE
)

# ── Performance plot — drop max.dist, legend outside ──────────────────────────
# mixOmics' default plot.perf draws all distance methods. We re-implement
# the same content using ggplot, dropping max.dist (which is misleading on
# imbalanced data) and placing the legend in empty space.
err_df <- bind_rows(
  data.frame(component = seq_along(perf_res$error.rate$overall[, "centroids.dist"]),
             error     = as.numeric(perf_res$error.rate$overall[, "centroids.dist"]),
             metric    = "Overall",
             distance  = "centroids.dist"),
  data.frame(component = seq_along(perf_res$error.rate$overall[, "mahalanobis.dist"]),
             error     = as.numeric(perf_res$error.rate$overall[, "mahalanobis.dist"]),
             metric    = "Overall",
             distance  = "mahalanobis.dist"),
  data.frame(component = seq_along(perf_res$error.rate$BER[, "centroids.dist"]),
             error     = as.numeric(perf_res$error.rate$BER[, "centroids.dist"]),
             metric    = "BER",
             distance  = "centroids.dist"),
  data.frame(component = seq_along(perf_res$error.rate$BER[, "mahalanobis.dist"]),
             error     = as.numeric(perf_res$error.rate$BER[, "mahalanobis.dist"]),
             metric    = "BER",
             distance  = "mahalanobis.dist")
)

p_perf <- ggplot(err_df,
                 aes(x = component, y = error,
                     color = distance, linetype = metric,
                     group = interaction(distance, metric))) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  geom_hline(yintercept = 0.75, linetype = "dotted",
             color = "grey40", linewidth = 0.5) +
  annotate("text", x = max(err_df$component), y = 0.76,
           label = "Random (BER = 0.75)", hjust = 1,
           color = "grey40", size = 3) +
  scale_color_manual(values = c("centroids.dist" = "#E67E22",
                                "mahalanobis.dist" = "#7F8C8D"),
                     name = "Distance") +
  scale_linetype_manual(values = c("Overall" = "solid",
                                   "BER"     = "dashed"),
                        name = "Metric") +
  scale_x_continuous(breaks = unique(err_df$component)) +
  labs(
    title = "PLS-DA Classification Error Rate\n(5-fold CV, 50 repeats)",
    x = "Component",
    y = "Classification error rate"
  ) +
  theme_pub +
  theme(
    legend.position = "right",
    legend.box      = "vertical"
  )

ggsave(file.path(outdir, "Fig1_PLSDA_performance.pdf"),
       p_perf, width = 7.5, height = 5.5, dpi = 600)
ggsave(file.path(outdir, "Fig1_PLSDA_performance.png"),
       p_perf, width = 7.5, height = 5.5, dpi = 600)
cat("[INFO] Saved Fig1_PLSDA_performance.pdf/png\n")

# ── VIP scores with numeric labels ────────────────────────────────────────────
vip_mat   <- vip(plsda_res)
vip_df    <- as.data.frame(vip_mat) %>%
  rownames_to_column("mark")
comp1_col <- names(vip_df)[2]
vip_df    <- vip_df %>%
  rename(VIP_comp1 = all_of(comp1_col)) %>%
  arrange(desc(VIP_comp1)) %>%
  mutate(significant = VIP_comp1 > 1)

write_tsv(vip_df, file.path(outdir, "Q1_PLSDA_VIP_scores.txt"))
cat("[INFO] VIP scores:\n")
print(vip_df)

p_vip <- ggplot(vip_df,
                aes(x = reorder(mark, VIP_comp1),
                    y = VIP_comp1,
                    fill = significant)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "#C0392B", linewidth = 0.8) +
  # Numeric VIP value at end of each bar
  geom_text(aes(label = sprintf("%.2f", VIP_comp1)),
            hjust = -0.15, size = 3.6, color = "black") +
  scale_fill_manual(
    values = c("TRUE" = "#C0392B", "FALSE" = "#85C1E9"),
    labels = c("TRUE" = "VIP > 1", "FALSE" = "VIP \u2264 1")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  coord_flip() +
  labs(
    x     = "",
    y     = "VIP Score (Component 1)",
    fill  = "",
    title = "Variable Importance in Projection\nParental Marks \u2192 Compartment Switch Prediction"
  ) +
  annotate("text", x = 0.6, y = 1.05,
           label = "VIP = 1.0\n(significance threshold)",
           color = "#C0392B", size = 3, hjust = 0) +
  theme_pub +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "Fig1_PLSDA_VIP.pdf"),
       p_vip, width = 6.5, height = 5, dpi = 600)
ggsave(file.path(outdir, "Fig1_PLSDA_VIP.png"),
       p_vip, width = 6.5, height = 5, dpi = 600)
cat("[INFO] Saved Fig1_PLSDA_VIP.pdf/png\n")

# ── Score plot — switchers in foreground, manual axis labels ──────────────────
# Strategy: pull mixOmics' scores directly and rebuild with ggplot so we
# can layer stables (faded background) under switchers (foreground),
# add proper axis labels ("PLS-DA Component 1/2"), and produce a cleaner
# Nature-style plot than plotIndiv's default.

scores_df <- data.frame(
  Comp1 = plsda_res$variates$X[, 1],
  Comp2 = plsda_res$variates$X[, 2],
  Class = Y
)
var_exp <- round(plsda_res$prop_expl_var$X[1:2] * 100, 1)

stables    <- scores_df %>% filter(Class %in% c("Stable_A", "Stable_B"))
switchers  <- scores_df %>% filter(Class %in% c("B_to_A",   "A_to_B"))

# Compute zoom limits from 1st-99th percentiles of all data
xlim_vals <- quantile(scores_df$Comp1, c(0.01, 0.99))
ylim_vals <- quantile(scores_df$Comp2, c(0.01, 0.99))
xpad      <- 0.05 * diff(xlim_vals)
ypad      <- 0.05 * diff(ylim_vals)

p_scores <- ggplot() +
  # Faded stable background
  geom_point(data = stables,
             aes(x = Comp1, y = Comp2, color = Class),
             size = 0.6, alpha = 0.10) +
  # Switchers on top
  geom_point(data = switchers,
             aes(x = Comp1, y = Comp2, color = Class),
             size = 1.4, alpha = 0.75) +
  # 95% confidence ellipses per class (drawn last so they sit above points)
  stat_ellipse(data = scores_df,
               aes(x = Comp1, y = Comp2, color = Class),
               level = 0.95, linewidth = 0.8) +
  scale_color_manual(values = switch_colors,
                     name   = "Switch class") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  coord_cartesian(
    xlim = c(xlim_vals[1] - xpad, xlim_vals[2] + xpad),
    ylim = c(ylim_vals[1] - ypad, ylim_vals[2] + ypad)
  ) +
  labs(
    title = "PLS-DA: Parental Marks vs Compartment Switch (log1p)",
    x = sprintf("PLS-DA Component 1 (%.1f%% variance)", var_exp[1]),
    y = sprintf("PLS-DA Component 2 (%.1f%% variance)", var_exp[2])
  ) +
  theme_pub +
  theme(legend.position = "right")

ggsave(file.path(outdir, "Fig1_PLSDA_scores.pdf"),
       p_scores, width = 7.5, height = 6, dpi = 600)
ggsave(file.path(outdir, "Fig1_PLSDA_scores.png"),
       p_scores, width = 7.5, height = 6, dpi = 600)
cat("[INFO] Saved Fig1_PLSDA_scores.pdf/png\n")

# ── Loading plot — method = "median", x-axis labeled ─────────────────────────
# method = "median" is more robust for skewed CPM signal than the default mean.
pdf(file.path(outdir, "Fig1_PLSDA_loadings.pdf"), width = 6.8, height = 5.5)
par(mar = c(5, 5, 4, 3))
plotLoadings(
  plsda_res,
  comp    = 1,
  method  = "mean",
  contrib = "max",
  title   = "",
  legend.title = "Max-contributing\nswitch class"
)
title(main = "PLS-DA Loadings: Component 1",
      line = 2.5, cex.main = 1, font.main = 2)
mtext("Loading weight (Component 1)", side = 1, line = 3.0, cex = 1)
dev.off()
cat("[INFO] Saved Fig1_PLSDA_loadings.pdf\n")

png(file.path(outdir, "Fig1_PLSDA_loadings.png"),
    width = 6.8, height = 5.5, units = "in", res = 600)
par(mar = c(5, 5, 4, 3))
plotLoadings(
  plsda_res,
  comp    = 1,
  method  = "mean",
  contrib = "max",
  title   = "",
  legend.title = "Max-contributing\nswitch class"
)
title(main = "PLS-DA Loadings: Component 1",
      line = 2.5, cex.main = 1, font.main = 2)
mtext("Loading weight (Component 1)", side = 1, line = 3.0, cex = 1)
dev.off()
cat("[INFO] Saved Fig1_PLSDA_loadings.png\n")

cat(sprintf("\n[INFO] PLS-DA complete -> %s\n", normalizePath(outdir)))