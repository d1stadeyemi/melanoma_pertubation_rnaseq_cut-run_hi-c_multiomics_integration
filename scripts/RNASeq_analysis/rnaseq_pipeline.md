# RNA-seq → Differential Expression Pipeline

A375 melanoma, Parental vs Bottom-10 (constricted-migration selected), ± UNC1999.
Two stages, four scripts. Everything takes a `--ref hg19|grch38` flag.

```
 fastq (R1/R2)
     │
     │  Stage 1 — preprocessing  (env: multiomics-py)
     ▼
 rnaseq_fastq_to_counts.sh          ← per-sample worker (trim → STAR → htseq → bigwig)
 submit_rnaseq_fastq_to_counts.sh   ← fans the worker out over a whole fastq directory
     │
     ▼
 per-sample counts:  <outdir>/<sample>/<sample>_genecounts.tsv
     │
     │  Stage 2 — differential expression  (env: multiomics-r)
     ▼
 rnaseq_deseq2_dge.R                ← counts → annotated DEGs + full diagnostics
 submit_rnaseq_deseq2_dge.sh        ← SLURM + conda wrapper for the R script
     │
     ▼
 <outdir>/tables/  <outdir>/figures/  MANIFEST.txt  sessionInfo.txt
```

Put all four scripts in `scripts/`. Each submitter finds its worker by filename in
the SAME directory (`SCRIPT_DIR`), so keep the pairs together.

---

## Scripts

| File | Role | You run | Key inputs | Key outputs |
|------|------|---------|------------|-------------|
| `rnaseq_fastq_to_counts.sh` | One paired-end pair: adapter/quality trim → STAR → sort/index → htseq-count → CPM bigwig | `sbatch` (usually via the batch script) | `--fastq-dir --r1 --r2 --ref [--outdir --fwd-keyword]` | `<replicate>/` dir with counts, sorted BAM, `.cpm.norm.bw`, plus `trimmed_fastqs/` and `sam/` subdirs (removable via flags) |
| `submit_rnaseq_fastq_to_counts.sh` | Finds every R1/R2 pair in a fastq dir, submits one job per pair, writes `samples_manifest.tsv` | `./…` on the login node | `--fastq-dir --ref [--outdir --fwd/--rev/--rep-keyword]` | one SLURM job per pair + `samples_manifest.tsv` |
| `rnaseq_deseq2_dge.R` | One two-group contrast → DEG tables + QC/diagnostic figures | `Rscript` or via the wrapper | `--counts-dir --ref` + `--group1/--group2` keywords (or `--samplesheet`) | `tables/`, `figures/`, `MANIFEST.txt` |
| `submit_rnaseq_deseq2_dge.sh` | Thin SLURM + `multiomics-r` wrapper; forwards all flags to the R script | `sbatch` | same flags as the R script | (runs the R script on a compute node) |

---

## Reference flag (`--ref`) — required in both stages

- **grch38**: gencode.v49 GTF. Counts are keyed on versioned Ensembl gene IDs
  (`ENSG…`); the DESeq2 script self-annotates gene symbol + biotype + chrom from the
  GTF. Nothing extra needed.
- **hg19**: UCSC knownGene GTF (`HomoSapiens_hg19.gtf`). Counts are keyed on `uc…`
  accessions and the GTF carries **no gene symbols**. The DESeq2 script auto-maps
  `uc… → symbol` using `TxDb.Hsapiens.UCSC.hg19.knownGene` + `org.Hs.eg.db`
  (installed in `multiomics-r`). Do **not** pass a modern UCSC `kgXref` id-map — its
  IDs are `ENST…` and will not match `uc…` counts.
- Always sanity-check `tables/05_annotation_report.csv` → `frac_annotated`,
  `symbols_available`, `frac_with_symbol`.

---

## File-naming flexibility (any sequencing vendor)

Fastq names are parsed by editable keywords, so a new vendor's naming needs only
flag changes, not code edits:

- `--fwd-keyword` / `--rev-keyword` (default `R1_001` / `R2_001`) — mark forward/
  reverse reads; used to pair files and to name each replicate's output directory.
- `--rep-keyword` (default `rep`, case-insensitive) — marks the replicate number;
  used only to derive the **biological-sample** name written to
  `samples_manifest.tsv`.

Example for `Ctrl_Rep1_R1.fastq.gz`:
`--fwd-keyword R1 --rev-keyword R2 --rep-keyword Rep`.

`samples_manifest.tsv` (columns: replicate, sample, R1, R2) is the bridge to
Stage 2 — it is exactly the mapping you need to write the DESeq2 `--samplesheet`
(`sample`=replicate, `condition`=sample), or to pick `--group*-samples`.

Output layout is **flat**: one directory per fastq pair, named `<replicate>`,
directly under `--outdir` (e.g. `OUTDIR/Parental_rep1/…`). Trimmed fastqs and the
SAM live in `trimmed_fastqs/` and `sam/` subdirs inside each replicate directory.

---

## Environments

| Stage | conda env | Why |
|-------|-----------|-----|
| preprocessing | `multiomics-py` | `bamCoverage` (deepTools); STAR/samtools/htseq come from `module load` |
| DESeq2 DGE | `multiomics-r` | R 4.4.3 + DESeq2 (+ annotation packages for hg19) |

Both wrappers apply the project env hygiene: `module purge`, `unset PYTHONPATH`,
`PYTHONNOUSERSITE=1`, caches redirected off `$HOME`.

---

## Internal references to keep in sync

Each submitter hardcodes its worker's filename (resolved next to itself):

| Submitter | Line | Must equal |
|-----------|------|------------|
| `submit_rnaseq_fastq_to_counts.sh` | `MAIN_SCRIPT=…` | the preprocessing worker's filename |
| `submit_rnaseq_deseq2_dge.sh` | `R_SCRIPT=…` | `rnaseq_deseq2_dge.R` |

If you rename a worker, update the one line in its submitter. Keep each submitter in
the same folder as its worker.

---

## End-to-end example (grch38)

```bash
# 1. Preprocess every sample in a fastq directory (submits N parallel jobs)
./submit_rnaseq_fastq_to_counts.sh \
    --ref grch38 \
    --fastq-dir /lustre/.../fastqs \
    --outdir    /lustre/.../rnaseq_output \
    --rm-trimmed-fastqs --rm-sam

# 2. When those finish, run one differential-expression contrast.
#    --group1/--group2 are case-insensitive KEYWORDS: every replicate folder in
#    --counts-dir whose name contains the keyword joins that group.
sbatch submit_rnaseq_deseq2_dge.sh \
    --ref grch38 \
    --counts-dir /lustre/.../rnaseq_output \
    --group1 Bottom10 --group2 Parental \
    --control Parental \
    --outdir results/exploratory/deseq2_B10_vs_Parental
    # optional: assert 3 reps each and set clean figure labels
    #   --group1-reps rep1,rep2,rep3 --group2-reps rep1,rep2,rep3
```

hg19 is identical except `--ref hg19` (symbols via auto-map; no `--id-map`).
The keywords are project-agnostic — any student swaps in their own group names
(e.g. `--group1 shCTRL --group2 shGENE`). For surgical control (dropping a bad
replicate, ± drug splits, multi-factor designs) use `--samplesheet samples.csv`
(columns `sample,condition[,batch,path]`, optionally `--batch batch`); the
`samples_manifest.tsv` from Stage 1 is a ready-made starting point.

Add `--dry-run` to the DESeq2 step first to validate files/design without running.

---

## Output layout (per DESeq2 contrast)

```
<outdir>/
  tables/    00_run_config, 01_sample_table, 02_htseq_assignment_qc,
             03_library_size / genes_detected, 05_annotation_report,
             06_prefilter_*, 07_size_factors, 11_pca_*, 12_sample_*,
             20_results_<contrast>_{full,significant,summary}, 22_threshold_sensitivity, …
  figures/   assignment bars, count distributions, RLE, dispersion, mean-SD,
             PCA (+scree), sample-distance & correlation heatmaps, Cook's,
             p-value histogram, independent-filtering, MA (raw+shrunk),
             volcano, top-gene counts, top-DEG heatmap   (each as .pdf + .png)
  MANIFEST.txt   sessionInfo.txt   run.log
```

Suggested placement: point `--outdir` under `results/exploratory/` while iterating,
promote final runs to `results/final/`.