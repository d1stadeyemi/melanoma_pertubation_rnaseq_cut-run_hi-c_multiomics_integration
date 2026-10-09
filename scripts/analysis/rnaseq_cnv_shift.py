#!/usr/bin/env python3
"""
rnaseq_cnv_shift.py — do Hi-C-called CNV segments show the expected expression shift?

Independent confirmation of Hi-C CNV calls using RNA-seq. If B10 lost one copy of
a region, the genes in that region should, ON AVERAGE, be expressed lower in B10
(dosage effect; usually somewhat smaller than the DNA change because of partial
dosage compensation). Individual genes vary; the block shift is the signal.

Input
-----
* DESeq2 results table for the SAME contrast direction as the Hi-C ratio
  (B10 vs Parental -> positive log2FC = higher in B10). Prefer the UNSHRUNK
  log2FoldChange; shrunken LFCs pull low-count genes toward 0 and understate shifts.
* Gene annotation GTF (hg19, same one used for counting) for gene coordinates.
* cnv_segments.tsv from hic_cnv_ratio.py.

Tests, per CNV segment (gain/loss)
----------------------------------
* median log2FC of expressed genes in the segment vs genome background
* Wilcoxon rank-sum p (genes in segment vs genes in CNV-neutral regions)
* block-permutation p: random runs of the same number of consecutive genes taken
  from CNV-neutral regions (keeps local co-regulation of neighbouring genes)
* direction agreement with Hi-C, and correlation of Hi-C vs RNA shift across segments

Example
-------
python scripts/analysis/rnaseq_cnv_shift.py \
  --deseq results/rnaseq/deseq2_B10_vs_Parental.csv \
  --gtf /path/to/hg19.gtf \
  --segments results/exploratory/hic_cnv/res_250000/cnv_segments.tsv \
  --outdir results/exploratory/hic_cnv/rnaseq_check
"""
from __future__ import annotations

import argparse
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

SERIES_GAIN = "#eb6834"
SERIES_LOSS = "#2a78d6"
POINT = "#8a8984"
INK = "#0b0b0b"
INK_2 = "#52514e"
SURFACE = "#fcfcfb"


def read_gtf_genes(path: str) -> pd.DataFrame:
    """Gene-level coordinates from a GTF (feature == gene; falls back to exon spans)."""
    opener = gzip.open if path.endswith(".gz") else open
    rows, exon_rows = [], []
    rx_id = re.compile(r'gene_id "([^"]+)"')
    rx_nm = re.compile(r'gene_name "([^"]+)"')
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "exon"):
                continue
            gid = rx_id.search(f[8])
            gnm = rx_nm.search(f[8])
            rec = (f[0], int(f[3]), int(f[4]), gid.group(1) if gid else "", gnm.group(1) if gnm else "")
            (rows if f[2] == "gene" else exon_rows).append(rec)
    if not rows:  # UCSC-style GTFs have no gene lines
        ex = pd.DataFrame(exon_rows, columns=["chrom", "start", "end", "gene_id", "gene_name"])
        g = ex.groupby(["gene_id", "gene_name", "chrom"], as_index=False).agg(start=("start", "min"),
                                                                             end=("end", "max"))
    else:
        g = pd.DataFrame(rows, columns=["chrom", "start", "end", "gene_id", "gene_name"])
    g["chrom"] = np.where(g["chrom"].str.startswith("chr"), g["chrom"], "chr" + g["chrom"])
    g["gene_id_nov"] = g["gene_id"].str.replace(r"\.\d+$", "", regex=True)
    g["mid"] = (g["start"] + g["end"]) // 2
    return g


def load_deseq(path: str, gene_col: str | None, lfc_col: str, bm_col: str) -> pd.DataFrame:
    sep = "\t" if path.endswith((".tsv", ".txt")) else ","
    d = pd.read_csv(path, sep=sep)
    if gene_col is None:
        first = d.columns[0]
        gene_col = first if (first.startswith("Unnamed") or first in ("", "gene", "gene_id", "X")) else None
        if gene_col is None:
            for c in ("gene_id", "gene", "ensembl_gene_id", "gene_name", "symbol"):
                if c in d.columns:
                    gene_col = c
                    break
    if gene_col is None or gene_col not in d.columns:
        sys.exit(f"ERROR: cannot find gene column; use --gene-col. Columns: {list(d.columns)}")
    for c in (lfc_col, bm_col):
        if c not in d.columns:
            sys.exit(f"ERROR: column '{c}' not in DESeq2 table. Columns: {list(d.columns)}")
    out = d[[gene_col, lfc_col, bm_col]].copy()
    out.columns = ["gene", "lfc", "baseMean"]
    out["gene_nov"] = out["gene"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    return out.dropna(subset=["lfc"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--deseq", required=True, help="DESeq2 results (csv/tsv), B10 vs Parental")
    ap.add_argument("--gtf", required=True, help="hg19 GTF used for counting (.gtf or .gtf.gz)")
    ap.add_argument("--segments", required=True, help="cnv_segments.tsv from hic_cnv_ratio.py")
    ap.add_argument("--gene-col", default=None)
    ap.add_argument("--lfc-col", default="log2FoldChange")
    ap.add_argument("--basemean-col", default="baseMean")
    ap.add_argument("--min-basemean", type=float, default=20, help="expressed-gene filter (default 20)")
    ap.add_argument("--min-genes", type=int, default=15, help="min expressed genes to test a segment")
    ap.add_argument("--n-perm", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--window", type=int, default=41, help="genes in rolling median for the plot")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    genes = read_gtf_genes(args.gtf)
    de = load_deseq(args.deseq, args.gene_col, args.lfc_col, args.basemean_col)

    # match by Ensembl ID (version-stripped) first, then by gene name
    m = de.merge(genes.drop_duplicates("gene_id_nov"), left_on="gene_nov", right_on="gene_id_nov", how="inner")
    if len(m) < 0.5 * len(de):
        m2 = de.merge(genes.drop_duplicates("gene_name"), left_on="gene", right_on="gene_name", how="inner")
        if len(m2) > len(m):
            m = m2
    print(f"DESeq2 genes: {len(de):,}; with coordinates: {len(m):,}")
    m = m[m["baseMean"] >= args.min_basemean].copy()
    m = m[~m["chrom"].isin(["chrY", "chrM", "chrMT"])]
    print(f"expressed (baseMean >= {args.min_basemean}): {len(m):,}")
    m = m.sort_values(["chrom", "mid"]).reset_index(drop=True)
    # centre so the genome-wide median shift is 0 (same convention as the Hi-C ratio)
    m["lfc_c"] = m["lfc"] - m["lfc"].median()

    seg = pd.read_csv(args.segments, sep="\t")
    m["segment_call"] = "neutral"
    m["segment_id"] = -1
    for i, r in seg.iterrows():
        hit = (m["chrom"] == r["chrom"]) & (m["mid"] >= r["start"]) & (m["mid"] < r["end"])
        m.loc[hit, "segment_id"] = i
        if r["call"] != "neutral":
            m.loc[hit, "segment_call"] = r["call"]

    neutral = m[m["segment_call"] == "neutral"]
    neutral_runs = [g["lfc_c"].to_numpy() for _, g in neutral.groupby("chrom")]

    def block_perm(n: int) -> np.ndarray:
        """Medians of random runs of n consecutive CNV-neutral genes."""
        ok = [x for x in neutral_runs if x.size >= n]
        if not ok:
            return np.array([])
        w = np.array([x.size - n + 1 for x in ok], dtype=float)
        w /= w.sum()
        out = np.empty(args.n_perm)
        picks = rng.choice(len(ok), size=args.n_perm, p=w)
        for k, j in enumerate(picks):
            s = rng.integers(ok[j].size - n + 1)
            out[k] = np.median(ok[j][s:s + n])
        return out

    rows = []
    for i, r in seg[seg["call"] != "neutral"].iterrows():
        g = m[m["segment_id"] == i]
        res = {"chrom": r["chrom"], "start": r["start"], "end": r["end"], "call": r["call"],
               "hic_log2": r["mean_all"], "n_genes": len(g)}
        if len(g) >= args.min_genes:
            med = float(np.median(g["lfc_c"]))
            null = block_perm(len(g))
            want_lower = r["call"] == "loss"
            if null.size:
                p_block = (1 + np.sum(null <= med if want_lower else null >= med)) / (1 + null.size)
            else:
                p_block = np.nan
            p_w = mannwhitneyu(g["lfc_c"], neutral["lfc_c"],
                               alternative="less" if want_lower else "greater").pvalue
            res.update({"rna_median_lfc": round(med, 3),
                        "rna_frac_genes_same_dir": round(float(np.mean(np.sign(g["lfc_c"]) ==
                                                                       (-1 if want_lower else 1))), 3),
                        "direction_agrees": (med < 0) == want_lower,
                        "p_wilcoxon_one_sided": p_w, "p_block_perm_one_sided": p_block,
                        "rna_to_dna_ratio": round(med / r["mean_all"], 2) if r["mean_all"] else np.nan})
        rows.append(res)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(args.outdir, "rnaseq_cnv_segment_shift.tsv"), sep="\t", index=False)
    m.to_csv(os.path.join(args.outdir, "rnaseq_genes_with_cnv.tsv"), sep="\t", index=False)

    tested = out.dropna(subset=["rna_median_lfc"]) if "rna_median_lfc" in out else out.iloc[0:0]
    with open(os.path.join(args.outdir, "rnaseq_cnv_summary.txt"), "w") as fh:
        fh.write(f"expressed genes: {len(m):,}; in CNV-changed segments: "
                 f"{int((m['segment_call'] != 'neutral').sum()):,}\n")
        fh.write(f"segments tested (>= {args.min_genes} genes): {len(tested)} of {len(out)}\n")
        if len(tested):
            fh.write(f"direction agrees with Hi-C: {int(tested['direction_agrees'].sum())} / {len(tested)}\n")
            if len(tested) >= 4:
                rho = spearmanr(tested["hic_log2"], tested["rna_median_lfc"]).correlation
                fh.write(f"Spearman rho (Hi-C segment shift vs RNA median shift): {rho:.2f}\n")
            fh.write("\n" + tested.to_string(index=False) + "\n")
        fh.write("\nReading: a CNV is supported when the RNA shift has the same sign, is block-wise\n"
                 "significant, and rna_to_dna_ratio is roughly 0.4-1 (partial dosage compensation).\n"
                 "Caveat: this RNA-seq contrast is confounded with library kit; kit effects are\n"
                 "gene-specific, not positional, so block shifts remain informative.\n")
    print(open(os.path.join(args.outdir, "rnaseq_cnv_summary.txt")).read())

    # ---- plot: rolling median of centred log2FC along each chromosome with CNV calls ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    chroms = [c for c in pd.unique(m["chrom"]) if (m["chrom"] == c).sum() > args.window]
    order = sorted(chroms, key=lambda c: (len(c), c))
    offsets, mids, cum = {}, [], 0
    for c in order:
        L = int(m.loc[m["chrom"] == c, "end"].max())
        offsets[c] = cum
        mids.append(cum + L / 2)
        cum += L
    fig, ax = plt.subplots(figsize=(15, 3.6), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    for c in order:
        g = m[m["chrom"] == c]
        x = offsets[c] + g["mid"].to_numpy()
        ax.scatter(x, np.clip(g["lfc_c"], -2, 2), s=1, color=POINT, alpha=0.25, linewidths=0, rasterized=True)
        roll = g["lfc_c"].rolling(args.window, center=True, min_periods=args.window // 2).median()
        ax.plot(x, roll, color=INK, lw=1)
        ax.axvline(offsets[c], color="#d9d8d3", lw=0.6, zorder=0)
    for r in seg[seg["call"] != "neutral"].itertuples():
        if r.chrom in offsets:
            ax.axvspan(offsets[r.chrom] + r.start, offsets[r.chrom] + r.end, alpha=0.18, lw=0,
                       color=SERIES_GAIN if r.call == "gain" else SERIES_LOSS)
    ax.axhline(0, color=INK_2, lw=0.8)
    ax.set_ylim(-1.5, 1.5)
    ax.set_xlim(0, cum)
    ax.set_xticks(mids)
    ax.set_xticklabels([c.replace("chr", "") for c in order], fontsize=8, color=INK_2)
    ax.set_ylabel("log2FC (centred)", color=INK_2, fontsize=9)
    ax.set_title("RNA-seq B10 vs Parental along the genome — line: rolling median; "
                 "shading: Hi-C CNV calls (orange gain, blue loss)", loc="left", fontsize=10, color=INK)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(args.outdir, f"rnaseq_cnv_genome.{ext}"), dpi=200, facecolor=SURFACE)


if __name__ == "__main__":
    main()
