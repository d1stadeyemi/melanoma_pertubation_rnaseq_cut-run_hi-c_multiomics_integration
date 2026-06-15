# =============================================================================
# pathway_enrichment.R — Functional / pathway enrichment for gene lists
#   (Q1-standard over-representation analysis with clusterProfiler)
# =============================================================================
#
# PROJECT: A375 melanoma parental-vs-B10 constricted-migration multi-omics.
#
# WHAT THIS DOES --------------------------------------------------------------
#   Takes the gene lists produced by compartment_reversal_analysis.py
#   (switchers, reverted, migration_specific_reversal, ...) and asks which
#   pathways / signatures they are enriched for, using over-representation
#   analysis (ORA, hypergeometric test) against a PROPER background. It runs all
#   lists together so you get a COMPARATIVE picture (what is specific to the
#   migration-specific set vs general switching), and it splits the analysis into
#   database blocks that answer different questions:
#     GO_BP / Reactome / Hallmark  -> "what processes" (e.g. cell migration, EMT)
#     PRC2                         -> "is this the Polycomb mechanism" (BENPORATH
#                                     PRC2/EED/SUZ12/H3K27ME3 target signatures)
#
# THE THING THAT MAKES OR BREAKS THIS: the BACKGROUND universe -----------------
#   ORA compares your list against a reference set. The correct reference is NOT
#   the whole genome — it is the genes that COULD have been in your list, i.e.
#   the genes assignable to the bins you actually tested (background_universe.txt
#   from the reversal script). Passing the genome inflates significance and a Q1
#   reviewer will reject the panel. --background is therefore effectively
#   required; if you omit it the script runs but prints a loud warning and the
#   result is NOT defensible.
#
# NOTE ON REGION-BASED ENRICHMENT (GREAT) -------------------------------------
#   Your units are genomic bins, so the most rigorous *primary* method is the
#   region-based test in GREAT (binomial over genomic space; handles gene length
#   and adjacent-bin non-independence). This script is the GENE-LIST ORA path you
#   asked for and is the standard human-readable companion to GREAT. If two genes
#   fall in the same bin, prefer collapsing to one-gene-per-bin upstream so
#   clustered duplicates do not inflate counts.
#
# WHAT IS REPORTED ------------------------------------------------------------
#   * FOLD ENRICHMENT (GeneRatio / BgRatio) alongside FDR — effect size first.
#   * BH-adjusted p (p.adjust / qvalue), gene-set size bounds, the member genes.
#   * GO redundancy reduced via clusterProfiler::simplify (semantic similarity).
#
# REQUIREMENTS (install once into the R env, e.g. multiomics-r) ----------------
#   BiocManager::install(c("clusterProfiler","ReactomePA","org.Hs.eg.db",
#                          "enrichplot","DOSE"))
#   install.packages(c("optparse","ggplot2","msigdbr"))
#   (msigdbr supplies MSigDB Hallmark + the C2:CGP PRC2/Polycomb sets; no manual
#    .gmt download. Human only here, via org.Hs.eg.db.)
#
# INPUTS ----------------------------------------------------------------------
#   --genelists  a DIRECTORY of .txt files (one gene symbol/ID per line; the file
#                basename is the list name) OR a comma-separated list of files.
#   --background a .txt of the background universe (one gene per line).
#
# EXAMPLE USAGE ---------------------------------------------------------------
#   # All lists in the reversal genelists dir, all four database blocks:
#   Rscript pathway_enrichment.R \
#       --genelists results/exploratory/reversal/genelists \
#       --background results/exploratory/reversal/genelists/background_universe.txt \
#       --databases GO_BP,Reactome,Hallmark,PRC2 \
#       --focal migration_specific_reversal \
#       --outdir results/exploratory/enrichment
#
#   # Just the PRC2/Polycomb mechanistic block, two explicit lists:
#   Rscript pathway_enrichment.R \
#       --genelists results/exploratory/reversal/genelists/migration_specific_reversal.txt,results/exploratory/reversal/genelists/reverted_under_UNC.txt \
#       --background results/exploratory/reversal/genelists/background_universe.txt \
#       --databases PRC2 --outdir results/exploratory/enrichment_prc2
#
#   # Inputs are Ensembl IDs instead of symbols:
#   Rscript pathway_enrichment.R ... --id_type ENSEMBL ...
#
# OUTPUTS (into --outdir) -----------------------------------------------------
#   enrich_<DB>_compare.csv            all lists, one DB (incl. FoldEnrichment)
#   enrich_<DB>_compare_dotplot.{pdf,png}   comparative dot plot across lists
#   enrich_<DB>_emap.{pdf,png}              enrichment-map network (redundancy)
#   focal_<DB>_<focal>.csv                  the focal list, one DB
#   focal_<DB>_<focal>_cnet.{pdf,png}       gene-concept network (driver genes)
#   enrichment_all_compare.csv         every DB stacked
#   enrichment_summary.txt             ID-mapping stats, counts, settings
# =============================================================================

suppressWarnings(suppressMessages({
  library(optparse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(enrichplot)
  library(ggplot2)
  library(msigdbr)
}))

# =============================================================================
# CLI
# =============================================================================
parse_cli <- function() {
  option_list <- list(
    make_option("--genelists", type = "character",
                help = "REQUIRED. Either a DIRECTORY of .txt gene lists (one gene per line; basename = list name) OR a comma-separated list of .txt files."),
    make_option("--background", type = "character", default = NULL,
                help = "STRONGLY RECOMMENDED. Background-universe .txt (genes assignable to the tested bins). Omitting it is not Q1-defensible; the script warns and falls back to all annotated genes."),

    make_option("--databases", type = "character", default = "GO_BP,Reactome,Hallmark,PRC2",
                help = "Comma-separated DB blocks to run. ALLOWED: GO_BP | GO_MF | GO_CC | Reactome | Hallmark | PRC2. [default %default]"),

    make_option("--id_type", type = "character", default = "SYMBOL",
                help = "Identifier type of the input genes. ALLOWED: SYMBOL | ENSEMBL | ENTREZID. [default %default]"),
    make_option("--organism", type = "character", default = "human",
                help = "Organism. ALLOWED: human (org.Hs.eg.db). Only human is wired up. [default %default]"),

    make_option("--focal", type = "character", default = NULL,
                help = "List name (file basename, no extension) used for the single-list emap/cnet figures. [default: auto — the list matching 'migration_specific', else the first]."),

    make_option("--pvalue", type = "double", default = 0.05,
                help = "Raw p-value cutoff passed to the enrich functions. Typical: 0.05 | 0.01. [default %default]"),
    make_option("--qvalue", type = "double", default = 0.05,
                help = "BH-adjusted q-value cutoff. Typical: 0.05 | 0.1 | 0.2. [default %default]"),
    make_option("--min_gs", type = "integer", default = 10L,
                help = "Minimum gene-set size to test. [default %default]"),
    make_option("--max_gs", type = "integer", default = 500L,
                help = "Maximum gene-set size to test (caps huge unspecific terms). [default %default]"),
    make_option("--simplify_cutoff", type = "double", default = 0.7,
                help = "GO redundancy reduction: semantic-similarity cutoff for clusterProfiler::simplify (0 disables; applies to GO_* only). [default %default]"),
    make_option("--top_n", type = "integer", default = 12L,
                help = "Number of top terms per list shown in the figures. [default %default]"),

    make_option("--prc2_regex", type = "character",
                default = "PRC2|POLYCOMB|SUZ12|EED|EZH2|H3K27ME3",
                help = "Regex (case-insensitive) selecting Polycomb/PRC2 signatures from MSigDB C2:CGP by set name. [default '%default']"),

    make_option("--outdir", type = "character", default = "enrichment_out",
                help = "Output directory. [default %default]")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$genelists)) stop("--genelists is required.")
  if (!opt$id_type %in% c("SYMBOL", "ENSEMBL", "ENTREZID"))
    stop("--id_type must be SYMBOL, ENSEMBL or ENTREZID.")
  if (tolower(opt$organism) != "human")
    stop("Only --organism human is wired up (org.Hs.eg.db).")
  allowed <- c("GO_BP", "GO_MF", "GO_CC", "Reactome", "Hallmark", "PRC2")
  dbs <- trimws(strsplit(opt$databases, ",")[[1]])
  bad <- setdiff(dbs, allowed)
  if (length(bad)) stop("Unknown --databases: ", paste(bad, collapse = ", "),
                        ". Allowed: ", paste(allowed, collapse = ", "))
  opt$.dbs <- dbs
  opt
}

# =============================================================================
# Small helpers (summary log)
# =============================================================================
log_lines <- new.env(); log_lines$txt <- character(0)
say <- function(...) { line <- paste0(...); cat(line, "\n"); log_lines$txt <- c(log_lines$txt, line) }

read_one_list <- function(path) {
  v <- readLines(path, warn = FALSE)
  v <- trimws(v)
  v <- v[nzchar(v) & !startsWith(v, "#")]
  unique(v)
}

read_genelists <- function(spec) {
  if (dir.exists(spec)) {
    files <- list.files(spec, pattern = "\\.txt$", full.names = TRUE)
    files <- files[!grepl("background|gene_assignment", basename(files), ignore.case = TRUE)]
  } else {
    files <- trimws(strsplit(spec, ",")[[1]])
  }
  miss <- files[!file.exists(files)]
  if (length(miss)) stop("Gene-list file(s) not found: ", paste(miss, collapse = ", "))
  if (!length(files)) stop("No gene-list files found in: ", spec)
  out <- lapply(files, read_one_list)
  names(out) <- tools::file_path_sans_ext(basename(files))
  out
}

# =============================================================================
# ID mapping to ENTREZ (clusterProfiler enrich functions key on ENTREZ here)
# =============================================================================
map_ids <- function(genes, from, label) {
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  if (from == "ENTREZID") {
    say(sprintf("[MAP] %-28s %d genes (already ENTREZ)", label, length(genes)))
    return(genes)
  }
  m <- tryCatch(suppressWarnings(suppressMessages(
    bitr(genes, fromType = from, toType = "ENTREZID", OrgDb = org.Hs.eg.db))),
    error = function(e) data.frame())
  if (!nrow(m)) {
    say(sprintf("[MAP] %-28s 0 of %d ids valid as %s — list skipped.", label, length(genes), from))
    return(character(0))
  }
  mapped <- unique(m$ENTREZID)
  say(sprintf("[MAP] %-28s %d input -> %d ENTREZ (%d unmapped)",
              label, length(genes), length(mapped),
              length(genes) - length(unique(m[[from]]))))
  mapped
}

# =============================================================================
# MSigDB term-to-gene tables (robust across msigdbr versions)
# =============================================================================
get_msigdbr_t2g <- function(collection, subcollection = NULL, name_regex = NULL) {
  df <- tryCatch({
    args <- list(species = "Homo sapiens", collection = collection)   # msigdbr >= 10
    if (!is.null(subcollection)) args$subcollection <- subcollection
    do.call(msigdbr::msigdbr, args)
  }, error = function(e) {
    args <- list(species = "Homo sapiens", category = collection)     # msigdbr < 10
    if (!is.null(subcollection)) args$subcategory <- subcollection
    do.call(msigdbr::msigdbr, args)
  })
  gene_col <- if ("ncbi_gene" %in% colnames(df)) "ncbi_gene" else "entrez_gene"
  t2g <- data.frame(term = df$gs_name, gene = as.character(df[[gene_col]]),
                    stringsAsFactors = FALSE)
  if (!is.null(name_regex))
    t2g <- t2g[grepl(name_regex, t2g$term, ignore.case = TRUE), , drop = FALSE]
  unique(t2g[!is.na(t2g$gene) & nzchar(t2g$gene), ])
}

# =============================================================================
# Fold enrichment (effect size) from GeneRatio / BgRatio
# =============================================================================
add_fold_enrichment <- function(df) {
  rr <- function(x) vapply(strsplit(as.character(x), "/"),
                           function(z) as.numeric(z[1]) / as.numeric(z[2]), numeric(1))
  if (nrow(df)) df$FoldEnrichment <- round(rr(df$GeneRatio) / rr(df$BgRatio), 3)
  df
}

# =============================================================================
# Enrichment runners (compareCluster for all lists; enrich_one for the focal)
# =============================================================================
run_compare <- function(clusters, universe, db, opt, msig) {
  base <- list(geneClusters = clusters, universe = universe,
               pvalueCutoff = opt$pvalue, qvalueCutoff = opt$qvalue,
               minGSSize = opt$min_gs, maxGSSize = opt$max_gs)
  cc <- if (db %in% c("GO_BP", "GO_MF", "GO_CC")) {
    do.call(compareCluster, c(base, list(fun = "enrichGO", OrgDb = org.Hs.eg.db,
            keyType = "ENTREZID", ont = sub("GO_", "", db), readable = TRUE)))
  } else if (db == "Reactome") {
    do.call(compareCluster, c(base, list(fun = "enrichPathway",
            organism = "human", readable = TRUE)))
  } else if (db %in% c("Hallmark", "PRC2")) {
    do.call(compareCluster, c(base, list(fun = "enricher", TERM2GENE = msig[[db]])))
  }
  if (!is.null(cc) && startsWith(db, "GO") && opt$simplify_cutoff > 0 &&
      nrow(as.data.frame(cc)) > 0)
    cc <- tryCatch(clusterProfiler::simplify(cc, cutoff = opt$simplify_cutoff,
                   by = "p.adjust", select_fun = min), error = function(e) cc)
  if (!is.null(cc) && db %in% c("Hallmark", "PRC2"))
    cc <- tryCatch(setReadable(cc, org.Hs.eg.db, "ENTREZID"), error = function(e) cc)
  cc
}

enrich_one <- function(genes, universe, db, opt, msig) {
  common <- list(gene = genes, universe = universe, pvalueCutoff = opt$pvalue,
                 qvalueCutoff = opt$qvalue, minGSSize = opt$min_gs, maxGSSize = opt$max_gs)
  r <- if (db %in% c("GO_BP", "GO_MF", "GO_CC")) {
    do.call(enrichGO, c(common, list(OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
            ont = sub("GO_", "", db), readable = TRUE)))
  } else if (db == "Reactome") {
    do.call(enrichPathway, c(common, list(organism = "human", readable = TRUE)))
  } else if (db %in% c("Hallmark", "PRC2")) {
    do.call(enricher, c(common, list(TERM2GENE = msig[[db]])))
  }
  if (!is.null(r) && startsWith(db, "GO") && opt$simplify_cutoff > 0 &&
      nrow(as.data.frame(r)) > 0)
    r <- tryCatch(clusterProfiler::simplify(r, cutoff = opt$simplify_cutoff,
                  by = "p.adjust", select_fun = min), error = function(e) r)
  if (!is.null(r) && db %in% c("Hallmark", "PRC2"))
    r <- tryCatch(setReadable(r, org.Hs.eg.db, "ENTREZID"), error = function(e) r)
  r
}

# =============================================================================
# Figures
# =============================================================================
save_gg <- function(p, fn, outdir, w = 9, h = 7) {
  ggplot2::ggsave(file.path(outdir, paste0(fn, ".pdf")), p, width = w, height = h, limitsize = FALSE)
  ggplot2::ggsave(file.path(outdir, paste0(fn, ".png")), p, width = w, height = h, dpi = 200, limitsize = FALSE)
}

compare_dotplot <- function(cc, db, opt, outdir) {
  n_clusters <- length(unique(as.data.frame(cc)$Cluster))
  w <- max(9, 2.4 * n_clusters)
  p <- tryCatch(
    dotplot(cc, showCategory = opt$top_n, label_format = 55) +
      ggtitle(paste0(db, " — enrichment across gene lists")) +
      theme(axis.text.y = element_text(size = 8)),
    error = function(e) NULL)
  if (!is.null(p)) save_gg(p, paste0("enrich_", db, "_compare_dotplot"), outdir, w = w, h = 8)
}

emap_fig <- function(obj, fn, opt, outdir) {
  tryCatch({
    obj <- enrichplot::pairwise_termsim(obj)
    p <- enrichplot::emapplot(obj, showCategory = opt$top_n)
    save_gg(p, fn, outdir, w = 9, h = 8)
  }, error = function(e) say(sprintf("   [fig] emap %s skipped: %s", fn, conditionMessage(e))))
}

cnet_fig <- function(r, fn, opt, outdir) {
  tryCatch({
    p <- enrichplot::cnetplot(r, showCategory = min(opt$top_n, 8),
                              cex_label_gene = 0.6)
    save_gg(p, fn, outdir, w = 10, h = 8)
  }, error = function(e) say(sprintf("   [fig] cnet %s skipped: %s", fn, conditionMessage(e))))
}

# =============================================================================
# Main
# =============================================================================
main <- function() {
  opt <- parse_cli()
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

  lists_sym <- read_genelists(opt$genelists)
  say(sprintf("[INFO] %d gene lists: %s", length(lists_sym),
              paste(names(lists_sym), collapse = ", ")))

  # background universe
  if (is.null(opt$background)) {
    say("[WARN] *** no --background given. ORA will use all annotated genes as the")
    say("        universe, which INFLATES significance and is NOT Q1-defensible.")
    say("        Pass background_universe.txt (genes assignable to the tested bins). ***")
    universe <- NULL
  } else {
    bg_sym <- read_one_list(opt$background)
    universe <- map_ids(bg_sym, opt$id_type, "background_universe")
    if (!length(universe)) stop("Background mapped to 0 ENTREZ ids — check --id_type.")
  }

  # map each list; drop empties
  clusters <- list()
  for (nm in names(lists_sym)) {
    e <- map_ids(lists_sym[[nm]], opt$id_type, nm)
    if (length(e)) clusters[[nm]] <- e else say(sprintf("[WARN] '%s' mapped to 0 genes; dropped.", nm))
  }
  if (!length(clusters)) stop("No gene lists survived ID mapping.")

  # focal list for single-list figures
  focal <- opt$focal
  if (is.null(focal)) {
    hit <- grep("migration_specific", names(clusters), value = TRUE)
    focal <- if (length(hit)) hit[1] else names(clusters)[1]
  }
  if (!focal %in% names(clusters))
    stop("--focal '", focal, "' is not one of: ", paste(names(clusters), collapse = ", "))
  say(sprintf("[INFO] focal list for emap/cnet = %s", focal))

  # precompute MSigDB tables once
  msig <- list()
  if ("Hallmark" %in% opt$.dbs) msig$Hallmark <- get_msigdbr_t2g("H")
  if ("PRC2" %in% opt$.dbs) {
    msig$PRC2 <- get_msigdbr_t2g("C2", "CGP", opt$prc2_regex)
    say(sprintf("[INFO] PRC2 block: %d Polycomb/PRC2 signatures matched '%s'.",
                length(unique(msig$PRC2$term)), opt$prc2_regex))
  }

  all_tables <- list()
  for (db in opt$.dbs) {
    say(sprintf("=== %s ===", db))
    cc <- tryCatch(run_compare(clusters, universe, db, opt, msig),
                   error = function(e) { say(sprintf("[ERROR] %s compareCluster: %s", db, conditionMessage(e))); NULL })
    if (is.null(cc) || nrow(as.data.frame(cc)) == 0) {
      say(sprintf("[%s] no enriched terms at the chosen cutoffs.", db)); next
    }
    df <- add_fold_enrichment(as.data.frame(cc))
    write.csv(df, file.path(opt$outdir, sprintf("enrich_%s_compare.csv", db)), row.names = FALSE)
    all_tables[[db]] <- cbind(Database = db, df)
    n_terms <- tapply(df$ID, df$Cluster, function(z) length(unique(z)))
    for (nm in names(n_terms)) {
      cnt <- n_terms[[nm]]; if (is.na(cnt)) cnt <- 0L
      say(sprintf("   %-28s %d enriched terms (FDR<%.2g)", nm, cnt, opt$qvalue))
    }

    compare_dotplot(cc, db, opt, opt$outdir)
    emap_fig(cc, sprintf("enrich_%s_emap", db), opt, opt$outdir)

    # focal single-list figures + table
    r <- tryCatch(enrich_one(clusters[[focal]], universe, db, opt, msig),
                  error = function(e) NULL)
    if (!is.null(r) && nrow(as.data.frame(r)) > 0) {
      rf <- add_fold_enrichment(as.data.frame(r))
      write.csv(rf, file.path(opt$outdir, sprintf("focal_%s_%s.csv", db, focal)), row.names = FALSE)
      cnet_fig(r, sprintf("focal_%s_%s_cnet", db, focal), opt, opt$outdir)
      emap_fig(r, sprintf("focal_%s_%s_emap", db, focal), opt, opt$outdir)
    }
  }

  if (length(all_tables)) {
    combined <- do.call(rbind, lapply(all_tables, function(d) {
      keep <- c("Database", "Cluster", "ID", "Description", "GeneRatio", "BgRatio",
                "FoldEnrichment", "pvalue", "p.adjust", "qvalue", "Count", "geneID")
      d[, intersect(keep, colnames(d)), drop = FALSE]
    }))
    write.csv(combined, file.path(opt$outdir, "enrichment_all_compare.csv"), row.names = FALSE)
  }

  writeLines(log_lines$txt, file.path(opt$outdir, "enrichment_summary.txt"))
  say(sprintf("[DONE] outputs in %s", normalizePath(opt$outdir)))
}

main()