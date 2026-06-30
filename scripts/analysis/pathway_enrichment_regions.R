#!/usr/bin/env Rscript
# =============================================================================
# region_pathway_enrichment.R
#   Pathway / functional enrichment for GENOMIC REGIONS, with the statistically
#   appropriate method chosen per input type.
# =============================================================================
#
# PROJECT: A375 melanoma parental-vs-B10 constricted-migration multi-omics.
#
# WHY TWO MODES (this is the whole point) -------------------------------------
#   You feed this script genomic regions, not pre-made gene lists, and the two
#   kinds of region you analyse have DIFFERENT statistical hazards, so they are
#   tested with different primary methods:
#
#   --mode region  : CUT&RUN H3K27me3 DIFFERENTIALLY-BOUND regions (broad).
#       Hazard 1: gene-locus-LENGTH bias — long loci collect peaks by chance
#                 (Welch et al. 2014, NAR 42:e105; Fisher's exact test inflated).
#       Hazard 2: H3K27me3 paints BROAD multi-gene domains; point-based tests
#                 are wrong by construction (Cavalcante et al. 2014, Bioinf
#                 30:i393 — Broad-Enrich models the PROPORTION of each gene
#                 locus covered by a peak, with a binomial smoothing spline for
#                 length).
#       PRIMARY METHOD: Broad-Enrich  ->  chipenrich::broadenrich()
#
#   --mode bin     : Hi-C COMPARTMENT-SWITCH bins (fixed 250 kb).
#       Fixed width => locus-length bias largely gone. The real hazards are:
#       (a) the correct BACKGROUND is "bins that COULD have switched" (bins with
#           a confident A/B call in BOTH conditions), NOT the genome
#           (Wijesooriya et al. 2022, PLoS Comput Biol 18:e1009935).
#       (b) adjacent bins are spatially autocorrelated (whole arms move together)
#           -> optional --collapse_adjacent merges contiguous switched bins into
#           one "domain" so a single event is not counted many times.
#       PRIMARY METHOD: clusterProfiler ORA with the bin-restricted universe.
#
# SHARED ACROSS MODES ---------------------------------------------------------
#   * Database blocks: GO_BP, Reactome, Hallmark (+ optional KEGG).
#   * POSITIVE CONTROL: a PRC2/Polycomb block built from MSigDB C2:CGP
#     (Ben-Porath SUZ12/EED/PRC2-target and ES_WITH_H3K27ME3 signatures).
#     If your differential-H3K27me3 / compartment sets do NOT light this up,
#     treat it as a QC red flag before interpreting anything else.
#   * options(enrichment_force_universe = TRUE) so clusterProfiler does NOT
#     silently drop unannotated genes from YOUR background
#     (the "background problem"; see Bioinformatics Advances 2024, vbae159).
#   * Effect size (fold enrichment / odds ratio) reported next to FDR.
#
# A NOTE ON NON-GENIC H3K27me3 ------------------------------------------------
#   H3K27me3 also occurs in intergenic space (broad local enrichments over
#   silent genes AND flanking intergenic regions; Pauler et al. 2009, Genome
#   Res 19:221). The 'nearest_tss' locus definition tiles the whole genome, so
#   intergenic regions are still assigned to a gene locus — but check the
#   qc_plots (distance-to-TSS, gene coverage). For signal that is genuinely
#   region-level, complement this gene-based analysis with a region-overlap test
#   (e.g. LOLA / ChIP-Atlas vs public EZH2/SUZ12 tracks) — that is OUT OF SCOPE
#   here and run separately.
#
# REQUIREMENTS (install once into your R env, e.g. multiomics-r) ---------------
#   BiocManager::install(c(
#     "chipenrich", "chipenrich.data",          # Broad-Enrich (region mode)
#     "clusterProfiler", "ReactomePA", "org.Hs.eg.db",
#     "enrichplot", "DOSE",
#     "TxDb.Hsapiens.UCSC.hg19.knownGene",       # bin->gene (hg19)
#     "TxDb.Hsapiens.UCSC.hg38.knownGene",       # bin->gene (hg38, if needed)
#     "GenomicFeatures", "GenomicRanges", "GenomeInfoDb",
#     "LOLA"))                                    # region-overlap companion
#   install.packages(c("optparse", "ggplot2", "msigdbr", "data.table"))
#   # LOLA also needs a region database on disk (the genome must match --genome):
#   #   * LOLACore (hg19/hg38): download + unpack from the LOLA regionDB page, OR
#   #   * a CUSTOM db: a parent folder with collection subfolders, each holding a
#   #     'regions/' dir of BED files (public EZH2/SUZ12/H3K27me3 tracks from
#   #     ChIP-Atlas/ENCODE) + an 'index.txt'. Files need only chr/start/end.
#
# INPUTS ----------------------------------------------------------------------
#   --regions     A DIRECTORY of BED-like files (chr<TAB>start<TAB>end[, score])
#                 OR a comma-separated list of such files. One file = one set;
#                 the basename (no extension) is the set name. Header optional.
#                 (region mode: differential H3K27me3 regions; bin mode: bins
#                  that switched compartment.)
#   --background  REQUIRED for bin mode, RECOMMENDED for region mode. A BED-like
#                 file of the TESTABLE universe:
#                   bin mode:    every bin with a confident A/B call in BOTH
#                                conditions (the bins that could have switched).
#                   region mode: every region/window the differential test could
#                                have called (used for the clusterProfiler
#                                companion universe; Broad-Enrich itself uses a
#                                genome-wide gene-locus null).
#
# OPTIONAL REGION-OVERLAP COMPANION (LOLA) ------------------------------------
#   --lola_db     Path to a LOLA region database (e.g. LOLACore/hg19, or a custom
#                 folder of public EZH2/SUZ12/H3K27me3 tracks). When supplied,
#                 the SAME region/bin sets are tested for direct overlap with
#                 those reference tracks via Fisher's exact test (Sheffield &
#                 Bock 2016, Bioinformatics 32:587). This is GENE-FREE: it answers
#                 "are these regions canonically Polycomb / EZH2-bound?" and so
#                 captures intergenic H3K27me3 that gene-based tests miss. It is
#                 the REGION-LEVEL twin of the Ben-Porath gene-level control, runs
#                 in BOTH modes, and uses --background (or --lola_universe) as the
#                 LOLA universe. The LOLA db genome MUST match --genome.
#
# EXAMPLE USAGE ---------------------------------------------------------------
#   # 1) CUT&RUN broad H3K27me3 (csaw output), split lost vs gained:
#   #    significant.csv = foreground; all_regions.csv = universe (LOLA/companion).
#   Rscript region_pathway_enrichment.R \
#       --mode region \
#       --regions results/.../csaw/Migrated_Parunc/spikein/migrated_..._significant.csv \
#       --background results/.../csaw/Migrated_Parunc/spikein/migrated_..._all_regions.csv \
#       --select "lost=call:lost,gained=call:gained" \
#       --genome hg19 --locusdef nearest_tss \
#       --databases GO_BP,Reactome,Hallmark,PRC2 \
#       --qc_plots --randomization bylength \
#       --outdir results/exploratory/enrichment_h3k27me3
#
#   # (or one combined set: drop --select; the whole significant.csv = one set)
#
#   # 2) Hi-C compartments straight from Reversal_classified.txt (sets = columns):
#   Rscript region_pathway_enrichment.R \
#       --mode bin \
#       --regions results/exploratory/reversal/Reversal_classified.txt \
#       --select "specific=migration_specific_reversal,reverted=reversal_status:reverted,switch_AtoB=mig_switch:A_to_B,switch_BtoA=mig_switch:B_to_A" \
#       --genome hg19 --databases GO_BP,Reactome,Hallmark,PRC2 \
#       --collapse_adjacent --focal specific \
#       --outdir results/exploratory/enrichment_compartment
#       # universe = ALL rows of Reversal_classified.txt (every callable bin)
#
#   # 3) Add KEGG (optional; live API in bin mode, built-in set in region mode):
#   Rscript region_pathway_enrichment.R ... --databases GO_BP,Reactome,Hallmark,PRC2,KEGG ...
#
#   # 4) Add the LOLA region-overlap companion (runs together, same inputs):
#   Rscript region_pathway_enrichment.R ... --lola_db /lustre/.../LOLACore/hg19 ...
#
# COLUMN / FORMAT NOTES -------------------------------------------------------
#   * Files may be tab OR comma, quoted OR unquoted, header OR not.
#   * chr/start/end are taken BY NAME when a header exists (override with
#     --chr_col/--start_col/--end_col), else by POSITION (cols 1-3). This is why
#     Reversal_classified.txt works even though col 1 is 'coordinate' (chr/start/
#     end are named cols 2-4), and csaw's quoted "chr","start","end" works too.
#
# OUTPUTS (into --outdir) -----------------------------------------------------
#   region mode:
#     broadenrich_<set>_<gs>_results.tab     raw Broad-Enrich tables (+ qc PDFs)
#     broadenrich_all_results.csv            every set/db stacked (+ FDR, effect)
#     broadenrich_<db>_top_<set>.{pdf,png}   top-term plots
#   bin mode:
#     ora_<db>_compare.csv                   compareCluster across sets (+ fold)
#     ora_<db>_compare_dotplot.{pdf,png}     comparative dot plot
#     focal_<db>_<focal>_{cnet,emap}.{pdf,png}
#     ora_all_compare.csv
#   LOLA (if --lola_db):
#     lola_region_overlap_results.csv        full Fisher overlap table (all tracks)
#     lola_polycomb_overlap.csv              EZH2/SUZ12/EED/H3K27me3 hits (control)
#     lola_top_enrichments.{pdf,png}         top reference tracks per set
#   both:
#     positive_control_PRC2.csv              Ben-Porath FDR per set (QC)
#     enrichment_summary.txt                 settings, mapping stats, caveats
# =============================================================================

suppressWarnings(suppressMessages({
  library(optparse)
}))

# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
parse_cli <- function() {
  ol <- list(
    make_option("--mode", type = "character", default = NULL,
      help = paste0("REQUIRED. Which data type you are analysing. ",
        "ALLOWED: region | cutrun | histone | mark  (CUT&RUN histone-mark differential ",
        "regions -> Broad-Enrich)  ||  bin | bins | compartment | hic  (Hi-C compartment ",
        "bins -> clusterProfiler ORA). The two are handled differently on purpose.")),
    make_option("--regions", type = "character", default = NULL,
      help = paste0("REQUIRED. EITHER (a) a directory / comma-separated list of region files ",
        "(each file = one set; basename = set name), OR (b) a SINGLE table from which sets are ",
        "carved with --select (e.g. csaw *_significant.csv, or Reversal_classified.txt). ",
        "Accepts tab/CSV, quoted or unquoted, with or without a header.")),
    make_option("--background", type = "character", default = NULL,
      help = paste0("Testable universe (all tested regions / all callable bins). ",
        "REQUIRED for --mode bin UNLESS --select is used on a table (then all rows of that ",
        "table are the universe). For csaw, this is the *_all_regions.csv. Recommended for region mode.")),

    # --- input parsing flexibility (applies to --regions and --background) ---
    make_option("--chr_col", type = "character", default = "chr",
      help = "Chromosome column NAME if the file has a header. Falls back to column 1 if absent. [default %default]"),
    make_option("--start_col", type = "character", default = "start",
      help = "Start column NAME (header files). Falls back to column 2 if absent. [default %default]"),
    make_option("--end_col", type = "character", default = "end",
      help = "End column NAME (header files). Falls back to column 3 if absent. [default %default]"),
    make_option("--sep", type = "character", default = "auto",
      help = "Field separator. ALLOWED: auto | comma | tab | ws. [default %default]"),
    make_option("--select", type = "character", default = NULL,
      help = paste0("Carve multiple sets out of ONE --regions table. Comma list of ",
        "name=column  (rows where column is TRUE/True/1) OR  name=column:value  (rows where ",
        "column == value). Example (compartments): ",
        "'specific=migration_specific_reversal,reverted=reversal_status:reverted'. ",
        "Example (csaw): 'lost=call:lost,gained=call:gained'.")),

    make_option("--genome", type = "character", default = "hg19",
      help = "ALLOWED: hg19 | hg38. [default %default]"),
    make_option("--databases", type = "character", default = "GO_BP,Reactome,Hallmark,PRC2",
      help = "Comma list. ALLOWED: GO_BP | GO_MF | GO_CC | Reactome | Hallmark | KEGG | PRC2. [default %default]"),

    # --- region mode (Broad-Enrich) ---
    make_option("--locusdef", type = "character", default = "nearest_tss",
      help = "[region mode] chipenrich locus definition. Common: nearest_tss | nearest_gene | 5kb | 10kb | exon | intron. 'nearest_tss' tiles the genome so every region is assigned. [default %default]"),
    make_option("--qc_plots", action = "store_true", default = FALSE,
      help = "[region mode] write Broad-Enrich QC plots (dist-to-TSS, gene coverage). Inspect these for genic/intergenic balance."),
    make_option("--randomization", type = "character", default = "NULL",
      help = "[region mode] Type-I-error check. ALLOWED: NULL | complete | bylength. 'bylength' preserves the locus-length relationship. [default %default]"),
    make_option("--cores", type = "integer", default = 1L,
      help = "[region mode] n_cores for Broad-Enrich. [default %default]"),

    # --- bin mode (clusterProfiler ORA) ---
    make_option("--collapse_adjacent", action = "store_true", default = FALSE,
      help = "[bin mode] merge contiguous switched bins into one domain before gene assignment (mitigates adjacent-bin autocorrelation)."),
    make_option("--bin_assign", type = "character", default = "tss",
      help = "[bin mode] how a gene is assigned to a bin. ALLOWED: tss (gene's TSS falls in the bin; one gene -> one bin) | body (any overlap of gene body). [default %default]"),
    make_option("--focal", type = "character", default = NULL,
      help = "[bin mode] set name used for the single-set emap/cnet figures. [default: auto — a set matching 'migration_specific', else the first]."),

    # --- LOLA region-overlap companion (optional; both modes) ---
    make_option("--lola_db", type = "character", default = NULL,
      help = "[optional] Path to a LOLA region database folder (e.g. LOLACore/hg19, or a custom EZH2/SUZ12/H3K27me3 collection). Genome MUST match --genome. When set, a gene-free region-overlap test runs on the SAME sets."),
    make_option("--lola_universe", type = "character", default = NULL,
      help = "[optional] BED-like file to use as the LOLA universe. [default: --background]. The query sets must be subsets of this universe."),

    # --- shared cutoffs ---
    make_option("--pvalue", type = "double", default = 0.05,
      help = "Raw p-value cutoff. [default %default]"),
    make_option("--qvalue", type = "double", default = 0.05,
      help = "BH q-value (FDR) cutoff. [default %default]"),
    make_option("--min_gs", type = "integer", default = 10L,
      help = "Minimum gene-set size. [default %default]"),
    make_option("--max_gs", type = "integer", default = 500L,
      help = "Maximum gene-set size. [default %default]"),
    make_option("--simplify_cutoff", type = "double", default = 0.7,
      help = "[bin mode] GO redundancy reduction via clusterProfiler::simplify (0 disables; GO_* only). [default %default]"),
    make_option("--top_n", type = "integer", default = 15L,
      help = "Top terms shown per set in figures. [default %default]"),
    make_option("--prc2_regex", type = "character",
      default = "BENPORATH|SUZ12|EED|PRC2|H3K27ME3|POLYCOMB",
      help = "Case-insensitive regex selecting Polycomb/PRC2 signatures from MSigDB C2:CGP by set name. [default '%default']"),

    make_option("--outdir", type = "character", default = "enrichment_out",
      help = "Output directory. [default %default]")
  )
  opt <- parse_args(OptionParser(option_list = ol))

  # normalise friendly mode aliases -> canonical {region, bin}
  if (is.null(opt$mode)) stop("--mode is required (region/cutrun/histone OR bin/compartment/hic).")
  m <- tolower(opt$mode)
  opt$mode <- if (m %in% c("region", "cutrun", "cut&run", "histone", "mark", "marks")) "region"
              else if (m %in% c("bin", "bins", "compartment", "compartments", "hic", "hi-c")) "bin"
              else stop("--mode '", opt$mode, "' not recognised. Use region/cutrun/histone or bin/compartment/hic.")

  if (is.null(opt$regions)) stop("--regions is required.")
  if (!opt$genome %in% c("hg19", "hg38")) stop("--genome must be hg19 or hg38.")
  if (!opt$sep %in% c("auto", "comma", "tab", "ws")) stop("--sep must be auto, comma, tab or ws.")
  if (opt$mode == "bin" && is.null(opt$background) && is.null(opt$select))
    stop("--mode bin needs the callable-bin universe: pass --background, OR use --select on a ",
         "table (e.g. Reversal_classified.txt) so all its rows become the universe.")
  if (!opt$bin_assign %in% c("tss", "body")) stop("--bin_assign must be tss or body.")
  if (!opt$randomization %in% c("NULL", "complete", "bylength"))
    stop("--randomization must be NULL, complete or bylength.")

  allowed <- c("GO_BP", "GO_MF", "GO_CC", "Reactome", "Hallmark", "KEGG", "PRC2")
  dbs <- trimws(strsplit(opt$databases, ",")[[1]])
  bad <- setdiff(dbs, allowed)
  if (length(bad)) stop("Unknown --databases: ", paste(bad, collapse = ", "),
                        ". Allowed: ", paste(allowed, collapse = ", "))
  opt$.dbs <- dbs
  opt
}

# -----------------------------------------------------------------------------
# small logger
# -----------------------------------------------------------------------------
LOG <- new.env(); LOG$txt <- character(0)
say <- function(...) { l <- paste0(...); cat(l, "\n"); LOG$txt <- c(LOG$txt, l) }

# -----------------------------------------------------------------------------
# robust table / region readers
#   Handles: tab OR comma OR whitespace; quoted OR unquoted fields; header OR not;
#   chr/start/end picked BY NAME (header) else by POSITION (cols 1-3); 'chr' prefix
#   enforced. Works on csaw quoted CSVs ("chr","start","end",...) and on the
#   compartment Reversal_classified.txt (coordinate, chr, start, end, ...).
# -----------------------------------------------------------------------------
.detect_sep <- function(first, sep_opt) {
  if (sep_opt == "comma") return(",")
  if (sep_opt == "tab")   return("\t")
  if (sep_opt == "ws")    return("")
  if (grepl("\t", first)) return("\t")
  if (grepl(",",  first)) return(",")
  ""                                   # whitespace
}

read_full_table <- function(path, sep_opt = "auto") {
  if (!file.exists(path)) stop("File not found: ", path)
  first <- readLines(path, n = 1, warn = FALSE)
  sep   <- .detect_sep(first, sep_opt)
  toks  <- strsplit(gsub("[\"']", "", first), if (sep == "") "[[:space:]]+" else sep)[[1]]
  # header present unless the 2nd/3rd tokens are already numeric coordinates
  chk <- toks[2:min(3, length(toks))]
  has_header <- length(toks) < 2 || any(is.na(suppressWarnings(as.numeric(chk))))
  read.table(path, header = has_header, sep = sep, quote = "\"'",
             stringsAsFactors = FALSE, comment.char = "", check.names = FALSE,
             fill = TRUE, na.strings = c("NA", "", "-"))
}

# pull chr/start/end (+ optional extra columns) out of a full table
pick_regions <- function(df, chr_col, start_col, end_col, keep = NULL) {
  cn <- colnames(df)
  if (all(c(chr_col, start_col, end_col) %in% cn)) {
    out <- data.frame(chr = df[[chr_col]], start = df[[start_col]],
                      end = df[[end_col]], stringsAsFactors = FALSE)
  } else {
    if (ncol(df) < 3) stop("Need >=3 columns (chr,start,end) and named cols not found.")
    out <- data.frame(chr = df[[1]], start = df[[2]], end = df[[3]], stringsAsFactors = FALSE)
  }
  out$start <- suppressWarnings(as.integer(out$start))
  out$end   <- suppressWarnings(as.integer(out$end))
  out$chr   <- trimws(as.character(out$chr))
  out$chr   <- ifelse(grepl("^chr", out$chr), out$chr, paste0("chr", out$chr))
  if (!is.null(keep)) for (k in keep) if (k %in% cn) out[[k]] <- df[[k]]
  out[!is.na(out$start) & !is.na(out$end), , drop = FALSE]
}

read_regions <- function(path, opt) {
  unique(pick_regions(read_full_table(path, opt$sep), opt$chr_col, opt$start_col, opt$end_col))
}

# parse --select "name=col[:value], ..." into a list of selectors
parse_select <- function(spec) {
  lapply(trimws(strsplit(spec, ",")[[1]]), function(it) {
    kv <- strsplit(it, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) stop("Bad --select item (need name=col[:value]): ", it)
    nm <- trimws(kv[1]); rhs <- trimws(kv[2])
    if (grepl(":", rhs, fixed = TRUE)) {
      cv <- strsplit(rhs, ":", fixed = TRUE)[[1]]
      list(name = nm, col = trimws(cv[1]), val = trimws(paste(cv[-1], collapse = ":")))
    } else list(name = nm, col = rhs, val = NULL)
  })
}

# rows of a full table matching one selector
.selector_rows <- function(df, sel) {
  if (!sel$col %in% colnames(df))
    stop("--select column '", sel$col, "' not in table. Columns: ",
         paste(colnames(df), collapse = ", "))
  v <- df[[sel$col]]
  if (is.null(sel$val)) {
    truthy <- c("TRUE", "True", "true", "T", "1", "yes", "Yes")
    df[trimws(as.character(v)) %in% truthy | (is.logical(v) & isTRUE(v %in% TRUE)), , drop = FALSE]
  } else {
    df[trimws(as.character(v)) == sel$val, , drop = FALSE]
  }
}

# Build named list of region sets from --regions (dir / list / single table+select)
read_region_sets <- function(opt) {
  spec <- opt$regions
  if (!is.null(opt$select)) {                         # one table -> many sets by column
    full <- read_full_table(spec, opt$sep)
    sels <- parse_select(opt$select)
    sets <- lapply(sels, function(s) {
      sub <- .selector_rows(full, s)
      unique(pick_regions(sub, opt$chr_col, opt$start_col, opt$end_col))
    })
    names(sets) <- vapply(sels, function(s) s$name, character(1))
    attr(sets, "full_table") <- full                  # stash for universe-from-table
    return(sets)
  }
  if (dir.exists(spec)) {                              # directory of per-set files
    files <- list.files(spec, pattern = "\\.(bed|txt|tsv|csv|broadPeak|narrowPeak)$",
                        full.names = TRUE, ignore.case = TRUE)
    files <- files[!grepl("background|universe|callable|tested|all_regions",
                          basename(files), ignore.case = TRUE)]
  } else {
    files <- trimws(strsplit(spec, ",")[[1]])          # comma list of per-set files
  }
  miss <- files[!file.exists(files)]
  if (length(miss)) stop("Region file(s) not found: ", paste(miss, collapse = ", "))
  if (!length(files)) stop("No region files found in: ", spec)
  sets <- lapply(files, read_regions, opt = opt)
  names(sets) <- tools::file_path_sans_ext(basename(files))
  sets
}

# optionally merge contiguous bins into domains (bin mode)
collapse_bins <- function(df) {
  if (!nrow(df)) return(df)
  if (!requireNamespace("GenomicRanges", quietly = TRUE))
    stop("GenomicRanges needed for --collapse_adjacent.")
  gr <- GenomicRanges::makeGRangesFromDataFrame(df[, c("chr", "start", "end")])
  gr <- GenomicRanges::reduce(gr)              # merge touching/overlapping ranges
  as.data.frame(gr)[, c("seqnames", "start", "end")] |>
    setNames(c("chr", "start", "end"))
}

# -----------------------------------------------------------------------------
# MSigDB tables (robust across msigdbr >=10 'collection' and <10 'category')
# -----------------------------------------------------------------------------
get_msig <- function(collection, subcollection = NULL, name_regex = NULL) {
  df <- tryCatch({
    a <- list(species = "Homo sapiens", collection = collection)
    if (!is.null(subcollection)) a$subcollection <- subcollection
    do.call(msigdbr::msigdbr, a)
  }, error = function(e) {
    a <- list(species = "Homo sapiens", category = collection)
    if (!is.null(subcollection)) a$subcategory <- subcollection
    do.call(msigdbr::msigdbr, a)
  })
  gene_col <- if ("ncbi_gene" %in% colnames(df)) "ncbi_gene" else "entrez_gene"
  out <- data.frame(term = df$gs_name, gene = as.character(df[[gene_col]]),
                    stringsAsFactors = FALSE)
  if (!is.null(name_regex))
    out <- out[grepl(name_regex, out$term, ignore.case = TRUE), , drop = FALSE]
  unique(out[!is.na(out$gene) & nzchar(out$gene), ])
}

# write a chipenrich custom geneset file:
#   tab-delimited, WITH header, col1 = geneset id, col2 = Entrez id (one pair/row)
#   (format per chipenrich vignette / read_geneset docs)
write_chipenrich_geneset <- function(t2g, path) {
  write.table(data.frame(gs_id = t2g$term, gene_id = t2g$gene),
              path, sep = "\t", quote = FALSE, row.names = FALSE)
  path
}

# fold enrichment from clusterProfiler GeneRatio / BgRatio
add_fold <- function(df) {
  rr <- function(x) vapply(strsplit(as.character(x), "/"),
                           function(z) as.numeric(z[1]) / as.numeric(z[2]), numeric(1))
  if (nrow(df) && all(c("GeneRatio", "BgRatio") %in% colnames(df)))
    df$FoldEnrichment <- round(rr(df$GeneRatio) / rr(df$BgRatio), 3)
  df
}

save_gg <- function(p, fn, outdir, w = 9, h = 7) {
  ggplot2::ggsave(file.path(outdir, paste0(fn, ".pdf")), p, width = w, height = h, limitsize = FALSE)
  ggplot2::ggsave(file.path(outdir, paste0(fn, ".png")), p, width = w, height = h, dpi = 200, limitsize = FALSE)
}

# =============================================================================
# bin -> gene assignment (bin mode) using a UCSC TxDb
# =============================================================================
bins_to_entrez <- function(df, genome, how) {
  txn <- if (genome == "hg19") "TxDb.Hsapiens.UCSC.hg19.knownGene"
         else                  "TxDb.Hsapiens.UCSC.hg38.knownGene"
  for (p in c(txn, "GenomicFeatures", "GenomicRanges", "GenomeInfoDb", "S4Vectors"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Package needed for bin mode: ", p)
  suppressMessages(library(txn, character.only = TRUE))
  txdb <- get(txn)

  genes <- suppressWarnings(GenomicFeatures::genes(txdb))   # GRanges, names = Entrez
  feat  <- if (how == "tss") GenomicRanges::resize(genes, width = 1, fix = "start") else genes

  bins <- GenomicRanges::makeGRangesFromDataFrame(df[, c("chr", "start", "end")])
  suppressWarnings(GenomeInfoDb::seqlevelsStyle(bins) <- "UCSC")

  ov <- GenomicRanges::findOverlaps(feat, bins, ignore.strand = TRUE)
  unique(names(feat)[S4Vectors::queryHits(ov)])
}

# =============================================================================
# REGION MODE  ->  Broad-Enrich  (primary, length-corrected, broad-aware)
# =============================================================================
run_region_mode <- function(opt, sets, builtin_map, custom_files) {
  if (!requireNamespace("chipenrich", quietly = TRUE))
    stop("chipenrich (and chipenrich.data) required for --mode region.")
  rand <- if (opt$randomization == "NULL") NULL else opt$randomization

  builtin_gs <- unname(unlist(builtin_map))          # e.g. c("GOBP","reactome")
  all_rows <- list()

  for (sn in names(sets)) {
    peaks <- sets[[sn]][, c("chr", "start", "end")]
    say(sprintf("=== [region] set '%s' : %d regions ===", sn, nrow(peaks)))

    # ---- pass 1: built-in databases (GO_BP / Reactome / KEGG) ----
    if (length(builtin_gs)) {
      r <- tryCatch(
        chipenrich::broadenrich(
          peaks = peaks, genome = opt$genome, genesets = builtin_gs,
          locusdef = opt$locusdef, qc_plots = opt$qc_plots, randomization = rand,
          out_name = paste0("broadenrich_", sn, "_builtin"),
          out_path = opt$outdir, n_cores = opt$cores),
        error = function(e) { say(sprintf("  [ERROR] builtin: %s", conditionMessage(e))); NULL })
      if (!is.null(r)) all_rows[[paste0(sn, "_builtin")]] <- cbind(Set = sn, r$results)
    }

    # ---- pass 2: custom databases (Hallmark, PRC2 positive control) ----
    for (cf in names(custom_files)) {
      r <- tryCatch(
        chipenrich::broadenrich(
          peaks = peaks, genome = opt$genome, genesets = custom_files[[cf]],
          locusdef = opt$locusdef, qc_plots = FALSE, randomization = rand,
          out_name = paste0("broadenrich_", sn, "_", cf),
          out_path = opt$outdir, n_cores = opt$cores),
        error = function(e) { say(sprintf("  [ERROR] %s: %s", cf, conditionMessage(e))); NULL })
      if (!is.null(r)) {
        rr <- r$results
        rr$Geneset.Type <- cf
        all_rows[[paste0(sn, "_", cf)]] <- cbind(Set = sn, rr)
      }
    }
  }

  if (!length(all_rows)) { say("[region] no results produced."); return(invisible()) }
  comb <- do.call(rbind, all_rows)
  write.csv(comb, file.path(opt$outdir, "broadenrich_all_results.csv"), row.names = FALSE)

  # effect-size column (Broad-Enrich reports Odds.Ratio / Effect / Status depending on version)
  eff_col <- intersect(c("Odds.Ratio", "Effect"), colnames(comb))
  eff_col <- if (length(eff_col)) eff_col[1] else NA
  fdr_col <- if ("FDR" %in% colnames(comb)) "FDR" else "P.value"

  # top-term plot per database block
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    desc_col <- intersect(c("Description", "Geneset.ID"), colnames(comb))[1]
    for (gt in unique(comb$Geneset.Type)) {
      sub <- comb[comb$Geneset.Type == gt & comb[[fdr_col]] < opt$qvalue, , drop = FALSE]
      if (!nrow(sub)) next
      sub <- do.call(rbind, by(sub, sub$Set, function(d)
        head(d[order(d[[fdr_col]]), ], opt$top_n)))
      sub$neglog10FDR <- -log10(pmax(sub[[fdr_col]], .Machine$double.xmin))
      aes_x <- if (!is.na(eff_col)) eff_col else "neglog10FDR"
      p <- ggplot2::ggplot(sub, ggplot2::aes(x = .data[[aes_x]],
                                             y = stats::reorder(.data[[desc_col]], neglog10FDR),
                                             size = neglog10FDR, colour = Set)) +
        ggplot2::geom_point(alpha = 0.85) +
        ggplot2::labs(title = paste0("Broad-Enrich — ", gt),
                      x = aes_x, y = NULL, size = "-log10 FDR") +
        ggplot2::theme_bw(base_size = 10)
      save_gg(p, paste0("broadenrich_", gt, "_top"), opt$outdir,
              w = 11, h = max(5, 0.32 * nrow(sub)))
    }
  }
  comb
}

# =============================================================================
# BIN MODE  ->  clusterProfiler ORA  (bin-restricted universe, force_universe)
# =============================================================================
run_bin_mode <- function(opt, sets, bg_df, msig) {
  for (p in c("clusterProfiler", "org.Hs.eg.db"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Package needed for bin mode: ", p)
  options(enrichment_force_universe = TRUE)   # keep YOUR background intact

  # build the restricted universe (genes assignable to ANY callable bin)
  universe <- bins_to_entrez(bg_df, opt$genome, opt$bin_assign)
  if (!length(universe)) stop("Background bins mapped to 0 genes — check --background / --genome.")
  say(sprintf("[bin] universe = %d genes assignable to %d callable bins.",
              length(universe), nrow(bg_df)))

  # foreground genes per set
  clusters <- list()
  for (sn in names(sets)) {
    df <- sets[[sn]]
    if (opt$collapse_adjacent) {
      n0 <- nrow(df); df <- collapse_bins(df)
      say(sprintf("[bin] '%s': collapsed %d bins -> %d domains.", sn, n0, nrow(df)))
    }
    g <- intersect(bins_to_entrez(df, opt$genome, opt$bin_assign), universe)
    if (length(g)) { clusters[[sn]] <- g
      say(sprintf("[bin] '%s': %d genes (within universe).", sn, length(g)))
    } else say(sprintf("[bin][WARN] '%s' mapped to 0 in-universe genes; dropped.", sn))
  }
  if (!length(clusters)) stop("No bin sets survived gene assignment.")

  focal <- opt$focal
  if (is.null(focal)) {
    hit <- grep("migration_specific", names(clusters), value = TRUE)
    focal <- if (length(hit)) hit[1] else names(clusters)[1]
  }
  if (!focal %in% names(clusters))
    stop("--focal '", focal, "' not among: ", paste(names(clusters), collapse = ", "))
  say(sprintf("[bin] focal set for emap/cnet = %s", focal))

  run_one_db <- function(db, compare = TRUE, genes = NULL) {
    base <- list(pvalueCutoff = opt$pvalue, qvalueCutoff = opt$qvalue,
                 minGSSize = opt$min_gs, maxGSSize = opt$max_gs)
    fun_args <- function() {
      if (db %in% c("GO_BP", "GO_MF", "GO_CC"))
        list(fun = "enrichGO", OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
             ont = sub("GO_", "", db), readable = TRUE)
      else if (db == "Reactome") list(fun = "enrichPathway", organism = "human", readable = TRUE)
      else if (db == "KEGG")     list(fun = "enrichKEGG", organism = "hsa")
      else                       list(fun = "enricher", TERM2GENE = msig[[db]])  # Hallmark / PRC2
    }
    if (compare) {
      cc <- do.call(clusterProfiler::compareCluster,
                    c(list(geneClusters = clusters, universe = universe), base, fun_args()))
      if (!is.null(cc) && startsWith(db, "GO") && opt$simplify_cutoff > 0 &&
          nrow(as.data.frame(cc)) > 0)
        cc <- tryCatch(clusterProfiler::simplify(cc, cutoff = opt$simplify_cutoff,
                       by = "p.adjust", select_fun = min), error = function(e) cc)
      if (!is.null(cc) && db %in% c("Hallmark", "PRC2", "KEGG"))
        cc <- tryCatch(clusterProfiler::setReadable(cc, org.Hs.eg.db, "ENTREZID"),
                       error = function(e) cc)
      cc
    } else {
      fa <- fun_args(); fn <- fa$fun; fa$fun <- NULL
      r <- do.call(getExportedValue(if (db == "Reactome") "ReactomePA" else "clusterProfiler", fn),
                   c(list(gene = genes, universe = universe), base, fa))
      if (!is.null(r) && startsWith(db, "GO") && opt$simplify_cutoff > 0 &&
          nrow(as.data.frame(r)) > 0)
        r <- tryCatch(clusterProfiler::simplify(r, cutoff = opt$simplify_cutoff,
                      by = "p.adjust", select_fun = min), error = function(e) r)
      if (!is.null(r) && db %in% c("Hallmark", "PRC2", "KEGG"))
        r <- tryCatch(clusterProfiler::setReadable(r, org.Hs.eg.db, "ENTREZID"),
                      error = function(e) r)
      r
    }
  }

  all_tab <- list()
  for (db in opt$.dbs) {
    say(sprintf("=== [bin] %s ===", db))
    cc <- tryCatch(run_one_db(db, compare = TRUE),
                   error = function(e) { say(sprintf("  [ERROR] %s: %s", db, conditionMessage(e))); NULL })
    if (is.null(cc) || nrow(as.data.frame(cc)) == 0) { say(sprintf("  [%s] no enriched terms.", db)); next }
    dff <- add_fold(as.data.frame(cc))
    write.csv(dff, file.path(opt$outdir, sprintf("ora_%s_compare.csv", db)), row.names = FALSE)
    all_tab[[db]] <- cbind(Database = db, dff)

    if (requireNamespace("ggplot2", quietly = TRUE)) {
      nC <- length(unique(dff$Cluster))
      p <- tryCatch(clusterProfiler::dotplot(cc, showCategory = opt$top_n, label_format = 55) +
                      ggplot2::ggtitle(paste0(db, " — enrichment across sets")),
                    error = function(e) NULL)
      if (!is.null(p)) save_gg(p, sprintf("ora_%s_compare_dotplot", db), opt$outdir,
                               w = max(9, 2.4 * nC), h = 8)
    }

    # focal single-set emap + cnet
    rf <- tryCatch(run_one_db(db, compare = FALSE, genes = clusters[[focal]]),
                   error = function(e) NULL)
    if (!is.null(rf) && nrow(as.data.frame(rf)) > 0 && requireNamespace("enrichplot", quietly = TRUE)) {
      tryCatch({
        rs <- enrichplot::pairwise_termsim(rf)
        save_gg(enrichplot::emapplot(rs, showCategory = opt$top_n),
                sprintf("focal_%s_%s_emap", db, focal), opt$outdir, 9, 8)
      }, error = function(e) say(sprintf("   [fig] emap %s skipped: %s", db, conditionMessage(e))))
      tryCatch(save_gg(enrichplot::cnetplot(rf, showCategory = min(opt$top_n, 8)),
                       sprintf("focal_%s_%s_cnet", db, focal), opt$outdir, 10, 8),
               error = function(e) say(sprintf("   [fig] cnet %s skipped: %s", db, conditionMessage(e))))
    }
  }
  if (length(all_tab)) {
    comb <- do.call(rbind, lapply(all_tab, function(d) {
      keep <- c("Database", "Cluster", "ID", "Description", "GeneRatio", "BgRatio",
                "FoldEnrichment", "pvalue", "p.adjust", "qvalue", "Count", "geneID")
      d[, intersect(keep, colnames(d)), drop = FALSE]
    }))
    write.csv(comb, file.path(opt$outdir, "ora_all_compare.csv"), row.names = FALSE)
    return(comb)
  }
  invisible()
}

# =============================================================================
# positive control summary (Ben-Porath PRC2) — works for both modes
# =============================================================================
write_pc <- function(comb, mode, regex, outdir) {
  if (is.null(comb) || !nrow(comb)) return(invisible())
  idcol  <- intersect(c("ID", "Geneset.ID", "Description"), colnames(comb))[1]
  fdrcol <- intersect(c("FDR", "p.adjust", "qvalue", "P.value", "pvalue"), colnames(comb))[1]
  setcol <- intersect(c("Set", "Cluster"), colnames(comb))[1]
  if (is.na(idcol) || is.na(fdrcol)) return(invisible())
  pc <- comb[grepl(regex, comb[[idcol]], ignore.case = TRUE), , drop = FALSE]
  if (!nrow(pc)) { say("[POS-CTRL] no Polycomb/PRC2 term matched — inspect inputs."); return(invisible()) }
  keep <- intersect(c(setcol, idcol, "Description", fdrcol, "Status", "Odds.Ratio",
                      "FoldEnrichment", "Count"), colnames(pc))
  pc <- pc[order(pc[[fdrcol]]), keep, drop = FALSE]
  write.csv(pc, file.path(outdir, "positive_control_PRC2.csv"), row.names = FALSE)
  say(sprintf("[POS-CTRL] %d Polycomb/PRC2 hits written (best %s = %.2g). ",
              nrow(pc), fdrcol, min(pc[[fdrcol]], na.rm = TRUE)),
      "Expect these significant; if not, treat downstream pathways with caution.")
}

# =============================================================================
# LOLA region-overlap companion (gene-free; both modes; optional)
#   Tests each region/bin set for overlap with every reference track in a LOLA
#   region database (Fisher's exact test), using --background as the universe.
#   Captures intergenic H3K27me3 that gene-based methods cannot; the EZH2/SUZ12/
#   EED/H3K27me3 hits are the REGION-LEVEL positive control.
# =============================================================================
df_to_granges <- function(df) {
  GenomicRanges::makeGRangesFromDataFrame(
    df[, c("chr", "start", "end")],
    seqnames.field = "chr", start.field = "start", end.field = "end")
}

run_lola <- function(opt, sets, universe_df) {
  for (p in c("LOLA", "GenomicRanges", "GenomeInfoDb"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Package needed for --lola_db: ", p)
  if (is.null(universe_df))
    stop("LOLA needs a universe: pass --background (testable regions / callable bins) ",
         "or --lola_universe. A union-of-query universe inflates significance.")
  say("=== [LOLA] region-overlap companion ===")
  say(sprintf("[LOLA] loading region DB: %s  (must be %s)", opt$lola_db, opt$genome))

  regionDB <- LOLA::loadRegionDB(opt$lola_db)

  # query sets -> GRangesList (named by set); universe -> single GRanges
  grl <- GenomicRanges::GRangesList(lapply(sets, df_to_granges))
  names(grl) <- names(sets)
  uni <- df_to_granges(universe_df)
  suppressWarnings(GenomeInfoDb::seqlevelsStyle(uni) <- "UCSC")
  suppressWarnings(GenomeInfoDb::seqlevelsStyle(grl) <- "UCSC")  # harmonise chr style

  # restrict query sets to the universe (LOLA-recommended; avoids many-to-many artifacts)
  grl <- tryCatch(LOLA::redefineUserSets(grl, uni), error = function(e) grl)
  tryCatch(LOLA::checkUniverseAppropriateness(grl, uni), error = function(e) NULL)

  res <- LOLA::runLOLA(grl, uni, regionDB, cores = opt$cores)
  rd  <- as.data.frame(res)
  if ("userSet" %in% colnames(rd))
    rd$userSetName <- names(grl)[as.integer(rd$userSet)]
  # order by FDR (qValue) if present, else by ranked p
  ord <- intersect(c("qValue", "pValueLog", "maxRnk"), colnames(rd))[1]
  if (!is.na(ord)) rd <- rd[order(rd[[ord]], decreasing = (ord == "pValueLog")), ]
  write.csv(rd, file.path(opt$outdir, "lola_region_overlap_results.csv"), row.names = FALSE)
  say(sprintf("[LOLA] tested %d set(s) vs %d reference tracks.",
              length(grl), length(regionDB$regionGRL)))

  # region-level positive control: EZH2 / SUZ12 / EED / H3K27me3 tracks
  ann_cols <- intersect(c("antibody", "description", "filename", "cellType",
                          "treatment", "tissue"), colnames(rd))
  if (length(ann_cols)) {
    hay <- apply(rd[, ann_cols, drop = FALSE], 1, paste, collapse = " ")
    pc  <- rd[grepl("EZH2|SUZ12|EED|H3K27me3|PRC2|polycomb", hay, ignore.case = TRUE), , drop = FALSE]
    if (nrow(pc)) {
      write.csv(pc, file.path(opt$outdir, "lola_polycomb_overlap.csv"), row.names = FALSE)
      best <- if ("qValue" %in% colnames(pc)) sprintf("min qValue=%.2g", min(pc$qValue, na.rm = TRUE))
              else sprintf("max pValueLog=%.1f", max(pc$pValueLog, na.rm = TRUE))
      say(sprintf("[LOLA][POS-CTRL] %d EZH2/SUZ12/EED/H3K27me3 overlaps (%s). ", nrow(pc), best),
          "Expect strong overlap; this is the gene-free Polycomb check.")
    } else say("[LOLA][POS-CTRL] no EZH2/SUZ12/H3K27me3 track in this DB matched — ",
               "check the DB has Polycomb tracks (LOLACore does; a minimal custom DB may not).")
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    g <- tryCatch(LOLA::plotTopLOLAEnrichments(res), error = function(e) NULL)
    if (!is.null(g)) tryCatch(save_gg(g, "lola_top_enrichments", opt$outdir, 11, 8),
                              error = function(e) NULL)
  }
  invisible(rd)
}

# =============================================================================
# MAIN
# =============================================================================
main <- function() {
  opt <- parse_cli()
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
  suppressWarnings(suppressMessages(library(msigdbr)))
  if (!requireNamespace("ggplot2", quietly = TRUE)) say("[WARN] ggplot2 missing — no figures.")

  say(sprintf("[INFO] mode=%s genome=%s databases=%s", opt$mode, opt$genome,
              paste(opt$.dbs, collapse = ",")))

  sets <- read_region_sets(opt)
  say(sprintf("[INFO] %d set(s): %s", length(sets),
              paste(sprintf("%s(%d)", names(sets), vapply(sets, nrow, integer(1))), collapse = ", ")))
  if (any(vapply(sets, nrow, integer(1)) == 0))
    say("[WARN] some set(s) are empty — check --select column names/values.")

  # universe: --background if given; else (when --select used) all rows of that table
  bg_df <- NULL
  if (!is.null(opt$background)) {
    bg_df <- read_regions(opt$background, opt)
  } else if (!is.null(attr(sets, "full_table"))) {
    bg_df <- unique(pick_regions(attr(sets, "full_table"),
                                 opt$chr_col, opt$start_col, opt$end_col))
    say("[INFO] universe = all rows of the --select table (no --background given).")
  }
  if (!is.null(bg_df)) say(sprintf("[INFO] background = %d regions/bins.", nrow(bg_df)))

  # ---- precompute MSigDB blocks needed in either mode ----
  msig <- list()
  if ("Hallmark" %in% opt$.dbs) msig$Hallmark <- get_msig("H")
  if ("PRC2" %in% opt$.dbs) {
    msig$PRC2 <- get_msig("C2", "CGP", opt$prc2_regex)
    say(sprintf("[INFO] PRC2 positive control: %d signatures matched '%s'.",
                length(unique(msig$PRC2$term)), opt$prc2_regex))
    if (!length(unique(msig$PRC2$term)))
      say("[WARN] PRC2 regex matched 0 sets — check msigdbr install / --prc2_regex.")
  }
  if ("KEGG" %in% opt$.dbs)
    say("[NOTE] KEGG enabled. Figures from KEGG maps need a free KEGG copyright form; ",
        "MSigDB KEGG sets are frozen-legacy. Reactome is the open, current alternative.")

  comb <- NULL
  if (opt$mode == "region") {
    # built-in chipenrich names for the databases it ships
    bmap <- list()
    if ("GO_BP"    %in% opt$.dbs) bmap$GO_BP    <- "GOBP"
    if ("GO_MF"    %in% opt$.dbs) bmap$GO_MF    <- "GOMF"
    if ("GO_CC"    %in% opt$.dbs) bmap$GO_CC    <- "GOCC"
    if ("Reactome" %in% opt$.dbs) bmap$Reactome <- "reactome"
    if ("KEGG"     %in% opt$.dbs) bmap$KEGG     <- "kegg_pathway"  # confirm via supported_genesets()
    # Hallmark + PRC2 go in as CUSTOM geneset files
    cfiles <- list()
    if ("Hallmark" %in% opt$.dbs)
      cfiles$Hallmark <- write_chipenrich_geneset(msig$Hallmark, file.path(opt$outdir, "_gs_hallmark.txt"))
    if ("PRC2" %in% opt$.dbs)
      cfiles$PRC2 <- write_chipenrich_geneset(msig$PRC2, file.path(opt$outdir, "_gs_prc2.txt"))

    say("[NOTE] region mode: Broad-Enrich uses a GENOME-WIDE gene-locus null and ",
        "corrects for locus length; --background is used only for context, not the test.")
    comb <- run_region_mode(opt, sets, bmap, cfiles)

  } else { # bin mode
    if (is.null(bg_df)) stop("bin mode needs a universe (--background or --select on a table).")
    say(sprintf("[INFO] universe = %d callable bins.", nrow(bg_df)))
    say("[CAVEAT] compartment bins are spatially autocorrelated; ORA assumes independent ",
        "genes. Use --collapse_adjacent and read results as domain-level, not bin-level.")
    comb <- run_bin_mode(opt, sets, bg_df, msig)
  }

  if ("PRC2" %in% opt$.dbs) write_pc(comb, opt$mode, opt$prc2_regex, opt$outdir)

  # ---- optional LOLA region-overlap companion (same inputs, runs together) ----
  if (!is.null(opt$lola_db)) {
    lola_uni <- if (!is.null(opt$lola_universe)) read_regions(opt$lola_universe, opt) else bg_df
    tryCatch(run_lola(opt, sets, lola_uni),
             error = function(e) say(sprintf("[LOLA][ERROR] %s", conditionMessage(e))))
  }

  writeLines(LOG$txt, file.path(opt$outdir, "enrichment_summary.txt"))
  say(sprintf("[DONE] outputs in %s", normalizePath(opt$outdir)))
}

main()