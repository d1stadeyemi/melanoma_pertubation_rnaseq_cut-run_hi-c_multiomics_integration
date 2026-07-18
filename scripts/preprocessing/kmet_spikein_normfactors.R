#!/usr/bin/env Rscript
# =============================================================================
# kmet_spikein_normfactors.R
#
# Turn SNAP-CUTANA K-MetStat barcode counts into:
#   (a) csaw/edgeR-ready norm.factors  (effective library size proportional to
#       spike-in recovery), and
#   (b) deepTools bamCoverage --scaleFactor values.
#
# Input:  the outputs of kmet_spikein_count.sh
# Output: normfactors.csv  (drop-in for csaw.R --normfactors)
#
# THE MATH (read this before you trust the numbers)
# -------------------------------------------------
# Let for sample i:
#   S_i = spike-in barcode read pairs (a fixed amount of substrate was added)
#   L_i = library size that csaw will use (total filtered/mapped fragments)
#   c_ij = csaw window counts
#
# A fixed number of spike-in nucleosomes went into every reaction. So S_i is a
# readout of "reads recovered per unit of substrate present" -- i.e. capture
# efficiency x sequencing depth. If the true absolute signal at window j is A_ij,
# then c_ij ~ A_ij * (S_i / N_spike), because the same antibody/MNase/library/
# sequencing chain acted on both spike-in and sample.
#
# We want c_ij / effLib_i ~ A_ij.  Therefore:
#
#       effLib_i  =  L_i * nf_i   MUST BE PROPORTIONAL TO  S_i
#       =>  nf_i  =  k * S_i / L_i
#
# edgeR convention: norm.factors are scaled so their geometric mean is 1, which
# fixes k. That is what this script emits.
#
# Note this is the exact reciprocal of EpiCypher's
#       "Scale Factor = % human reads / % spike-in reads"
# and note that formula simplifies: (H_i/T_i) / (S_i/T_i) = H_i/S_i. The total
# read count T_i cancels. So "%human/%spike" and "human reads / spike reads" are
# the SAME NUMBER -- there is nothing to choose between them.
# =============================================================================

# base R only -- no external package dependencies

# ----------------------------- tiny arg parser -------------------------------
USAGE <- "
kmet_spikein_normfactors.R -- K-MetStat spike-in factors for csaw/edgeR and deepTools

USAGE
  Rscript kmet_spikein_normfactors.R --matrix kmet_barcode_matrix.tsv \\
                                     --summary kmet_summary.tsv \\
                                     --groups manifest_kmet.csv \\
                                     --group-by Condition,Treatment,Factor \\
                                     --plot kmet_qc \\
                                     [--libsizes csaw_totals.tsv] [options]

REQUIRED
  --matrix FILE       kmet_barcode_matrix.tsv from kmet_spikein_count.sh
  --summary FILE      kmet_summary.tsv from kmet_spikein_count.sh

OPTIONS
  --libsizes FILE     2-column TSV/CSV (sample, libsize) giving the library size csaw will
                      actually use, i.e. totals(data). STRONGLY RECOMMENDED. If omitted, the
                      aligned_fragments column of kmet_summary.tsv is used instead.
  --denominator VAL   Which spike-in reads define S_i. One of:
                        on_target   reads on the antibody's own PTM barcode  [default]
                        total       all 16 PTM barcodes summed (EpiCypher's convention)
                        ptm:<NAME>  one specific barcode, e.g. ptm:H3K27me3
  --reference NAME    Anchor deepTools scale factors to 1 for this sample. Does not affect
                      csaw norm.factors. Default: geometric-mean centring.
  --min-spikein INT   Refuse to emit factors when S_i is below this. Default: 1000
  --out FILE          Output CSV. Default: normfactors.csv
  --groups FILE       Optional manifest CSV/TSV giving the REAL replicate groups. Without
                      this the CV check groups by --summary's `target` column, which is a
                      single value when every sample is the same mark -- so the CV it
                      reports is meaningless. Needs a `folder` (or `sample`) column plus
                      the columns named in --group-by. If it also has `kmet_ab_ratio`,
                      an A/B balance panel is added to the plot.
  --group-by COLS     Comma-separated columns in --groups defining a replicate group.
                      Default: Condition,Treatment
  --plot FILE         Optional diagnostic plot. Writes BOTH <file>.pdf and <file>.png
                      (a .pdf/.png extension on FILE is stripped first).
  -h, --help          Print this help and exit.

IMPORTANT -- matching denominators
  edgeR's effective library size is lib.size * norm.factors. csaw sets lib.size from
  totals(data), which counts the fragments that survived YOUR read filter (mapq, dedup,
  proper pairs, blacklist). If the L_i you feed this script is not that same number, the
  norm.factors are wrong by the ratio between them. Either pass --libsizes with totals(data),
  or make sure --mapq / --samtools-flags in kmet_spikein_count.sh reproduce your csaw filter.
"

parse_args_simple <- function(args, defaults) {
  opt <- defaults
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% c("-h", "--help")) { cat(USAGE); quit(status = 0) }
    if (!grepl("^--", a)) { message("ERROR: unexpected argument '", a, "'. Use --help."); quit(status = 1) }
    key <- sub("^--", "", a)
    if (grepl("=", key, fixed = TRUE)) {
      kv <- strsplit(key, "=", fixed = TRUE)[[1]]
      key <- kv[1]; val <- paste(kv[-1], collapse = "=")
    } else {
      if (i == length(args)) { message("ERROR: --", key, " needs a value."); quit(status = 1) }
      i <- i + 1; val <- args[[i]]
    }
    if (!(key %in% names(defaults))) {
      message("ERROR: unknown option --", key, ". Allowed: ",
              paste0("--", names(defaults), collapse = " "), ". Use --help.")
      quit(status = 1)
    }
    opt[[key]] <- val
    i <- i + 1
  }
  opt
}

opt <- parse_args_simple(commandArgs(trailingOnly = TRUE), list(
  matrix = NULL, summary = NULL, libsizes = NULL,
  denominator = "on_target", reference = NULL,
  `min-spikein` = "1000", out = "normfactors.csv", plot = NULL,
  groups = NULL, `group-by` = "Condition,Treatment"
))
opt$`min-spikein` <- suppressWarnings(as.integer(opt$`min-spikein`))
if (is.na(opt$`min-spikein`)) { message("ERROR: --min-spikein must be an integer."); quit(status = 1) }

stopifnot_msg <- function(cond, msg) if (!isTRUE(cond)) { message("ERROR: ", msg); quit(status = 1) }

stopifnot_msg(!is.null(opt$matrix),  "--matrix is required.")
stopifnot_msg(!is.null(opt$summary), "--summary is required.")
stopifnot_msg(file.exists(opt$matrix),  paste("cannot read", opt$matrix))
stopifnot_msg(file.exists(opt$summary), paste("cannot read", opt$summary))

mat  <- read.delim(opt$matrix,  check.names = FALSE, stringsAsFactors = FALSE)
summ <- read.delim(opt$summary, check.names = FALSE, stringsAsFactors = FALSE)

rownames(mat) <- mat$ptm
mat$ptm <- NULL
samples <- colnames(mat)
stopifnot_msg(all(samples %in% summ$sample), "samples in --matrix are missing from --summary.")
summ <- summ[match(samples, summ$sample), ]

# ---- S_i : spike-in reads that define capture efficiency ---------------------
den <- opt$denominator
if (den == "total") {
  S <- colSums(mat)
  den_label <- "all 16 K-MetStat barcodes"
} else if (grepl("^ptm:", den)) {
  p <- sub("^ptm:", "", den)
  stopifnot_msg(p %in% rownames(mat), paste0("PTM '", p, "' is not in the matrix."))
  S <- as.numeric(mat[p, ]); names(S) <- samples
  den_label <- paste0(p, " barcode")
} else if (den == "on_target") {
  tg <- summ$target
  bad <- tg %in% c("NA", "IgG", "") | is.na(tg)
  stopifnot_msg(!any(bad), paste0(
    "--denominator on_target needs a K-MetStat PTM in the 'target' column for every sample.\n",
    "       Offending samples: ", paste(samples[bad], collapse = ", "), "\n",
    "       If your antibody is not in the panel (e.g. H3K27ac), K-MetStat cannot normalize it.\n",
    "       Use --denominator total only if you know what that assumes, or drop those samples."))
  S <- vapply(seq_along(samples), function(i) as.numeric(mat[tg[i], samples[i]]), numeric(1))
  names(S) <- samples
  den_label <- "each sample's own on-target PTM barcode"
} else {
  message("ERROR: --denominator must be 'on_target', 'total', or 'ptm:<NAME>'."); quit(status = 1)
}

# ---- L_i : library size csaw will use ---------------------------------------
if (!is.null(opt$libsizes)) {
  sep <- if (grepl("\\.csv$", opt$libsizes)) "," else "\t"
  lb <- read.delim(opt$libsizes, sep = sep, stringsAsFactors = FALSE)
  stopifnot_msg(ncol(lb) >= 2, "--libsizes needs at least 2 columns (sample, libsize).")
  colnames(lb)[1:2] <- c("sample", "libsize")
  stopifnot_msg(all(samples %in% lb$sample), "--libsizes is missing some samples.")
  L <- as.numeric(lb$libsize[match(samples, lb$sample)])
  L_src <- paste("user-supplied:", opt$libsizes)
} else {
  L <- suppressWarnings(as.numeric(summ$aligned_fragments))
  L_src <- "aligned_fragments column of kmet_summary.tsv (from --bam)"
  stopifnot_msg(!any(is.na(L)), paste0(
    "no library sizes available.\n",
    "       Either rerun kmet_spikein_count.sh with a 'bam' column in the manifest,\n",
    "       or pass --libsizes with totals(data) from csaw."))
}
names(L) <- samples

# ---- guardrails --------------------------------------------------------------
low <- S < opt$`min-spikein`
if (any(low)) {
  message("ERROR: these samples have fewer than ", opt$`min-spikein`, " spike-in reads (S_i):")
  message("       ", paste(sprintf("%s=%d", samples[low], S[low]), collapse = ", "))
  message("       EpiCypher's own guidance is that thousands of barcode reads are needed for the")
  message("       panel to sample reliably. Normalizing on counts this small propagates Poisson")
  message("       noise straight into every logFC. Fix the experiment or use --min-spikein to override.")
  quit(status = 1)
}

frac <- S / L
cv <- function(x) if (length(x) < 2) NA_real_ else stats::sd(x) / mean(x)

# ---- real replicate groups ---------------------------------------------------
# The `target` column is NOT a replicate group when every sample is the same mark.
grp <- setNames(rep(NA_character_, length(samples)), samples)
abr <- setNames(rep(NA_real_, length(samples)), samples)
if (!is.null(opt$groups)) {
  stopifnot_msg(file.exists(opt$groups), paste("cannot read", opt$groups))
  sepg <- if (grepl("\\.tsv$|\\.txt$", opt$groups)) "\t" else ","
  gt <- read.delim(opt$groups, sep = sepg, stringsAsFactors = FALSE, check.names = FALSE)
  key <- if ("folder" %in% names(gt)) "folder" else if ("sample" %in% names(gt)) "sample" else NA
  stopifnot_msg(!is.na(key), "--groups needs a `folder` or `sample` column.")
  gcols <- trimws(strsplit(opt$`group-by`, ",")[[1]])
  miss <- setdiff(gcols, names(gt))
  stopifnot_msg(length(miss) == 0,
    paste0("--groups is missing column(s): ", paste(miss, collapse = ", "),
           "\n       Available: ", paste(names(gt), collapse = ", ")))
  idx <- match(samples, gt[[key]])
  stopifnot_msg(!any(is.na(idx)),
    paste0("these samples are absent from --groups: ",
           paste(samples[is.na(idx)], collapse = ", ")))
  grp <- setNames(apply(gt[idx, gcols, drop = FALSE], 1,
                        function(r) paste(r, collapse = "_")), samples)
  if ("kmet_ab_ratio" %in% names(gt))
    abr <- setNames(suppressWarnings(as.numeric(gt$kmet_ab_ratio[idx])), samples)
} else {
  grp <- setNames(as.character(summ$target), samples)
}

# ---- the factors -------------------------------------------------------------
nf_raw <- S / L
nf <- nf_raw / exp(mean(log(nf_raw)))          # edgeR convention: geometric mean 1

sf_raw <- 1 / S                                 # deepTools multiplies coverage
if (!is.null(opt$reference)) {
  stopifnot_msg(opt$reference %in% samples, paste("--reference", opt$reference, "not found."))
  sf <- sf_raw / sf_raw[opt$reference]
} else {
  sf <- sf_raw / exp(mean(log(sf_raw)))
}

out <- data.frame(
  sample            = samples,
  target            = summ$target,
  spikein_reads_S   = as.integer(S),
  libsize_L         = as.numeric(L),
  spikein_fraction  = signif(frac, 4),
  norm_factor       = signif(nf, 6),            # <- for csaw/edgeR: y$samples$norm.factors
  effective_libsize = round(L * nf),
  deeptools_scaleFactor = signif(sf, 6),
  stringsAsFactors  = FALSE
)

write.csv(out, opt$out, row.names = FALSE, quote = FALSE)

# ---- report ------------------------------------------------------------------
cat("\n=== K-MetStat spike-in normalization ===\n")
cat("S_i defined by : ", den_label, "\n", sep = "")
cat("L_i source     : ", L_src, "\n", sep = "")
cat("\n"); print(out, row.names = FALSE); cat("\n")

cat("Within-group consistency check (the thing that decides whether this is publishable):\n")
if (is.null(opt$groups))
  cat("  !! --groups not given: grouping by `target`. If every sample is the same mark this\n",
      "     is ONE group and the CV below is meaningless. Pass --groups data/raw/manifest.csv\n",
      "     --group-by Condition,Treatment for the real answer.\n", sep = "")
cv_tab <- data.frame()
for (g in unique(grp)) {
  i <- which(grp == g)
  if (length(i) < 2) { 
    cat(sprintf("  %-28s n=1  (no within-group estimate)\n", g)); next 
  }
  cvg <- 100 * cv(frac[i])
  cv_tab <- rbind(cv_tab, data.frame(group = g, n = length(i), cv = cvg))
  cat(sprintf("  %-28s n=%d  CV = %6.1f%%  (%.3g - %.3g)%s\n",
              g, length(i), cvg, min(frac[i]), max(frac[i]),
              if (cvg > 30) "   <-- FLAG" else ""))
}
cat("\n  A CV above ~25-30% within a replicate group means the spike-in is reporting technical\n")
cat("  noise, not capture efficiency, and the resulting factors are not trustworthy for global\n")
cat("  magnitude claims. Check this before you use these factors in a figure.\n\n")

cat("How to use:\n")
cat("  csaw/edgeR : y$samples$norm.factors <- nf[match(colnames(y), samples)]\n")
cat("               (make sure y$samples$lib.size equals libsize_L above)\n")
cat("  deepTools  : bamCoverage --scaleFactor <deeptools_scaleFactor> (do NOT also pass --normalizeUsing)\n\n")

if (!is.null(opt$plot)) {
  base <- sub("\\.(pdf|png)$", "", opt$plot, ignore.case = TRUE)

  ord  <- order(grp, samples)
  gcol <- as.integer(factor(grp[ord]))
  pal  <- grDevices::hcl.colors(max(gcol), "Dark 3")

  has_ab <- any(!is.na(abr))
  npan   <- 3 + as.integer(has_ab)

  draw <- function() {
    op <- par(mfrow = c(2, 2), mar = c(11.5, 4.6, 3, 1.2), mgp = c(3.2, 0.7, 0),
              cex.axis = 0.62, cex.lab = 0.9, cex.main = 0.95)
    on.exit(par(op), add = TRUE)

    # 1. spike-in fraction, log scale (range spans >25x)
    yl <- range(frac) * c(0.85, 6)   # headroom so the legend clears the bars
    barplot(frac[ord], las = 2, log = "y", names.arg = samples[ord], ylim = yl,
            col = pal[gcol], border = NA, ylab = "S_i / L_i  (log scale)",
            main = "K-MetStat spike-in fraction")
    abline(h = exp(mean(log(frac))), lty = 3, col = "grey30")
    legend("top", legend = levels(factor(grp[ord])), fill = pal, bty = "n",
           cex = 0.55, ncol = 4, xpd = NA)

    # 2. norm.factors, log scale, reference at 1
    barplot(nf[ord], las = 2, log = "y", names.arg = samples[ord],
            ylim = range(nf) * c(0.85, 1.6),
            col = pal[gcol], border = NA, ylab = "norm.factor (log scale)",
            main = "csaw/edgeR norm.factors")
    abline(h = 1, lty = 2, col = "red")
    mtext("dashed red = 1.0 (no rescaling)", side = 3, line = -1.1, cex = 0.5,
          col = "red")

    # 3. within-group CV -- the go/no-go
    if (nrow(cv_tab) > 0) {
      cvv <- cv_tab$cv; names(cvv) <- paste0(cv_tab$group, " (n=", cv_tab$n, ")")
      par(mar = c(11.5, 4.6, 3, 1.2))
      barplot(cvv, las = 2, col = ifelse(cvv > 30, "firebrick", "seagreen"),
              border = NA, ylab = "CV of S_i/L_i  (%)",
              main = "Within-group consistency", ylim = c(0, max(35, max(cvv) * 1.15)))
      abline(h = 30, lty = 2, col = "red")
      mtext("above red line = not usable for global magnitude", side = 3, line = -1.1,
            cex = 0.5, col = "firebrick")
    } else {
      plot.new(); title(main = "Within-group consistency")
      text(0.5, 0.5, "no group has n>=2\n(pass --groups)", cex = 0.8)
    }

    # 4. A/B barcode balance -- EpiCypher's within-reaction technical replicate
    if (has_ab) {
      par(mar = c(11.5, 4.6, 3, 1.2))
      barplot(abr[ord], las = 2, names.arg = samples[ord],
              col = ifelse(!is.na(abr[ord]) & (abr[ord] < 0.77 | abr[ord] > 1.3),
                           "firebrick", "grey60"),
              border = NA, ylab = "on-target barcode A / B",
              main = "A/B balance (expect ~1.0)")
      abline(h = 1, lty = 2, col = "red")
      abline(h = c(0.77, 1.3), lty = 3, col = "grey40")
      mtext("red = imbalanced; often marks a different panel lot / run", side = 3,
            line = 0.2, cex = 0.5, col = "firebrick")
    }
  }

  pdf(paste0(base, ".pdf"), width = 12, height = 10); draw(); invisible(dev.off())
  ok_png <- TRUE
  tryCatch({
    if (capabilities("cairo")) {
      png(paste0(base, ".png"), width = 12, height = 10, units = "in", res = 200,
          type = "cairo")
    } else {
      png(paste0(base, ".png"), width = 12 * 200, height = 10 * 200, res = 200)
    }
    draw(); invisible(dev.off())
  }, error = function(e) {
    ok_png <<- FALSE
    message("NOTE: PNG device unavailable (", conditionMessage(e), "); wrote PDF only.")
  })

  cat("Diagnostic plots: ", base, ".pdf", if (ok_png) paste0(" and ", base, ".png") else "",
      "\n\n", sep = "")
}