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

Design (mirrors hic_cnv_ratio.py)
---------------------------------
--test LABEL:results.csv   B10 vs Parental contrasts (DMSO, UNC). A CNV-driven shift must
                           replicate in every test contrast, in the direction of the Hi-C call.
--null LABEL:results.csv   UNC vs DMSO within the same cells. Copy number is identical, so
                           CNV segments should show NO block shift (two-sided test).
(--deseq FILE is still accepted as a single test contrast.)

Example
-------
python scripts/analysis/rnaseq_cnv_shift.py \
  --test B10vsPar_DMSO:tables/20_results_Bottom-10-DMSO_vs_Parental-DMSO_full.csv \
  --test B10vsPar_UNC:tables/..._Bottom-10-UNC_vs_Parental-UNC_full.csv \
  --null Par_UNCvsDMSO:tables/..._Parental-UNC_vs_Parental-DMSO_full.csv \
  --null B10_UNCvsDMSO:tables/..._Bottom-10-UNC_vs_Bottom-10-DMSO_full.csv \
  --gtf data/refs/hg19/gencode.v19.annotation.gtf.gz \
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

# same colours as hic_cnv_ratio.py: test 1 blue, test 2 orange, nulls grey; violet loss, yellow gain
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7"]
NULL_COLOR = "#8a8984"
LOSS_TINT, LOSS_ALPHA = "#4a3aa7", 0.13
GAIN_TINT, GAIN_ALPHA = "#eda100", 0.22
# hg19 chromosome lengths, so the x-axis matches the Hi-C figure exactly
HG19 = {"chr1": 249250621, "chr2": 243199373, "chr3": 198022430, "chr4": 191154276, "chr5": 180915260,
        "chr6": 171115067, "chr7": 159138663, "chr8": 146364022, "chr9": 141213431, "chr10": 135534747,
        "chr11": 135006516, "chr12": 133851895, "chr13": 115169878, "chr14": 107349540, "chr15": 102531392,
        "chr16": 90354753, "chr17": 81195210, "chr18": 78077248, "chr19": 59128983, "chr20": 63025520,
        "chr21": 48129895, "chr22": 51304566, "chrX": 155270560}
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
    for c in (lfc_col, bm_col):
        if c not in d.columns:
            sys.exit(f"ERROR: column '{c}' not in DESeq2 table. Columns: {list(d.columns)}")
    if gene_col is not None and gene_col not in d.columns:
        sys.exit(f"ERROR: --gene-col '{gene_col}' not found. Columns: {list(d.columns)}")
    # every text column is a candidate gene identifier (Ensembl ID, symbol, ...)
    id_cols = [gene_col] if gene_col else [c for c in d.columns
               if d[c].dtype == object or pd.api.types.is_string_dtype(d[c])]
    if not id_cols:
        sys.exit(f"ERROR: no text column to use as gene ID. Columns: {list(d.columns)}")
    # DESeq2 annotation columns (chrom, biotype, ...) are not gene IDs and would collide with
    # the GTF's hg19 coordinates - drop them; prefix the ID columns to avoid name clashes
    not_ids = {"chrom", "chr", "seqnames", "seqname", "chromosome", "strand", "biotype",
               "gene_biotype", "gene_type", "description"}
    id_cols = [c for c in id_cols if c.lower() not in not_ids]
    if not id_cols:
        sys.exit(f"ERROR: no gene-ID column found. Columns: {list(d.columns)}")
    out = d[id_cols + [lfc_col, bm_col]].copy()
    out = out.rename(columns={lfc_col: "lfc", bm_col: "baseMean", **{c: f"de_{c}" for c in id_cols}})
    out.attrs["id_cols"] = [f"de_{c}" for c in id_cols]
    return out.dropna(subset=["lfc"])


def attach_coords(de: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    """Try every (DESeq2 id column) x (GTF key) combination; keep the best-matching one."""
    strip = lambda x: x.astype(str).str.replace(r"\.\d+$", "", regex=True)
    best, best_n, best_desc = None, -1, ""
    id_cols = de.attrs["id_cols"]
    for dc in id_cols:
        key_de = strip(de[dc])
        for gk in ("gene_id_nov", "gene_name", "gene_id"):
            g = genes[genes[gk].astype(str) != ""].drop_duplicates(gk)
            key_g = strip(g[gk]) if gk != "gene_id_nov" else g[gk]
            n = key_de.isin(set(key_g)).sum()
            if n > best_n:
                tmp = de.assign(_k=key_de.values)
                best = tmp.merge(g.assign(_k=key_g.values), on="_k", how="inner")
                best_n, best_desc = n, f"DESeq2 column '{dc}' <-> GTF {gk}"
    if best_n <= 0:
        print("ERROR: no gene identifiers matched between the DESeq2 table and the GTF.", file=sys.stderr)
        for dc in id_cols:
            print(f"  DESeq2 '{dc}' examples: {list(de[dc].astype(str).head(5))}", file=sys.stderr)
        print(f"  GTF gene_id examples:   {list(genes['gene_id'].head(5))}", file=sys.stderr)
        print(f"  GTF gene_name examples: {list(genes['gene_name'].head(5))}", file=sys.stderr)
        sys.exit("Use a GTF whose IDs match the DESeq2 table (Ensembl IDs <-> GENCODE/Ensembl GTF; "
                 "symbols <-> any GTF with gene_name).")
    print(f"matched by {best_desc}: {best_n:,} genes")
    return best


def parse_contrast(spec: str) -> tuple[str, str]:
    if ":" not in spec:
        sys.exit(f"ERROR: expected LABEL:path, got '{spec}'")
    lab, path = spec.split(":", 1)
    if not os.path.exists(path):
        sys.exit(f"ERROR: file not found: {path}")
    return lab, path


def chrom_key(c: str):
    b = c.replace("chr", "")
    return (0, int(b), "") if b.isdigit() else (1, 0, b)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--test", action="append", default=[], help="LABEL:deseq.csv (B10 vs Parental); repeat")
    ap.add_argument("--null", action="append", default=[], help="LABEL:deseq.csv (UNC vs DMSO, same cells); repeat")
    ap.add_argument("--deseq", default=None, help="single test contrast (backward compatible)")
    ap.add_argument("--gtf", required=True, help="hg19 GTF with Ensembl IDs, e.g. GENCODE v19 (.gtf/.gtf.gz)")
    ap.add_argument("--segments", required=True, help="cnv_segments.tsv from hic_cnv_ratio.py")
    ap.add_argument("--gene-col", default=None)
    ap.add_argument("--lfc-col", default="log2FoldChange", help="use UNSHRUNK log2FC (default)")
    ap.add_argument("--basemean-col", default="baseMean")
    ap.add_argument("--min-basemean", type=float, default=20, help="expressed-gene filter (default 20)")
    ap.add_argument("--min-genes", type=int, default=15, help="min expressed genes to test a segment")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--null-frac", type=float, default=0.25,
                    help="a null shift is 'material' if |shift| > this x mean |test shift| (default 0.25)")
    ap.add_argument("--n-perm", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--window", type=int, default=41, help="genes in rolling median for the plot")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    tests = [parse_contrast(x) for x in args.test]
    if args.deseq:
        tests.insert(0, ("B10vsPar", args.deseq))
    nulls = [parse_contrast(x) for x in args.null]
    if not tests:
        sys.exit("ERROR: give at least one --test (or --deseq)")
    contrasts = [(l, p, "test") for l, p in tests] + [(l, p, "null") for l, p in nulls]

    genes = read_gtf_genes(args.gtf)
    seg = pd.read_csv(args.segments, sep="\t")
    altered = seg[seg["call"] != "neutral"]

    per_gene, rows = {}, []
    for lab, path, kind in contrasts:
        print(f"\n== {lab} ({kind}): {path}")
        de = load_deseq(path, args.gene_col, args.lfc_col, args.basemean_col)
        m = attach_coords(de, genes)
        m = m[(m["baseMean"] >= args.min_basemean) & ~m["chrom"].isin(["chrY", "chrM", "chrMT"])].copy()
        print(f"expressed genes with coordinates: {len(m):,}")
        if len(m) < 100:
            sys.exit("ERROR: too few expressed genes with coordinates - check IDs / --min-basemean")
        m = m.sort_values(["chrom", "mid"]).reset_index(drop=True)
        m["lfc_c"] = m["lfc"] - m["lfc"].median()
        m["segment_call"], m["segment_id"] = "neutral", -1
        for i, r in seg.iterrows():
            hit = (m["chrom"] == r["chrom"]) & (m["mid"] >= r["start"]) & (m["mid"] < r["end"])
            m.loc[hit, "segment_id"] = i
            if r["call"] != "neutral":
                m.loc[hit, "segment_call"] = r["call"]
        per_gene[lab] = m
        neutral = m[m["segment_call"] == "neutral"]
        runs = [g["lfc_c"].to_numpy() for _, g in neutral.groupby("chrom")]

        def block_null(n):
            ok = [x for x in runs if x.size >= n]
            if not ok:
                return np.array([])
            w = np.array([x.size - n + 1 for x in ok], float)
            picks = rng.choice(len(ok), size=args.n_perm, p=w / w.sum())
            out = np.empty(args.n_perm)
            for k, j in enumerate(picks):
                s0 = rng.integers(ok[j].size - n + 1)
                out[k] = np.median(ok[j][s0:s0 + n])
            return out

        for i, r in altered.iterrows():
            g = m[m["segment_id"] == i]
            res = {"contrast": lab, "kind": kind, "chrom": r["chrom"], "start": r["start"], "end": r["end"],
                   "call": r["call"], "hic_log2": r["mean_all"], "n_genes": len(g)}
            if len(g) >= args.min_genes:
                med = float(np.median(g["lfc_c"]))
                nul = block_null(len(g))
                if kind == "test":   # one-sided, in the direction of the Hi-C call
                    lower = r["call"] == "loss"
                    p_blk = (1 + np.sum(nul <= med if lower else nul >= med)) / (1 + nul.size) if nul.size else np.nan
                    p_w = mannwhitneyu(g["lfc_c"], neutral["lfc_c"], alternative="less" if lower else "greater").pvalue
                else:                # null contrasts: two-sided, any shift is unexpected
                    p_blk = (1 + np.sum(np.abs(nul) >= abs(med))) / (1 + nul.size) if nul.size else np.nan
                    p_w = mannwhitneyu(g["lfc_c"], neutral["lfc_c"], alternative="two-sided").pvalue
                res.update({"rna_median_lfc": round(med, 3), "p_wilcoxon": p_w, "p_block_perm": p_blk,
                            "rna_to_dna_ratio": round(med / r["mean_all"], 2) if kind == "test" and r["mean_all"] else np.nan})
            rows.append(res)

    long = pd.DataFrame(rows)
    long.to_csv(os.path.join(args.outdir, "rnaseq_cnv_segment_shift_long.tsv"), sep="\t", index=False)
    for lab, m in per_gene.items():
        m.to_csv(os.path.join(args.outdir, f"rnaseq_genes_with_cnv_{lab}.tsv"), sep="\t", index=False)

    # ---- per-segment verdict ----
    verdict = []
    for (c, s0, e0), g in long.groupby(["chrom", "start", "end"], sort=False):
        t = g[g["kind"] == "test"].dropna(subset=["rna_median_lfc"])
        n = g[g["kind"] == "null"].dropna(subset=["rna_median_lfc"])
        call = g["call"].iloc[0]
        row = {"chrom": c, "start": s0, "end": e0, "call": call, "hic_log2": g["hic_log2"].iloc[0],
               "n_genes": int(g["n_genes"].max())}
        for _, x in g.iterrows():
            row[f"lfc_{x['contrast']}"] = x.get("rna_median_lfc", np.nan)
            row[f"pblk_{x['contrast']}"] = x.get("p_block_perm", np.nan)
        if len(t) < len(tests):
            row["verdict"] = "not testable (too few genes)"
        else:
            sign_ok = all((v < 0) == (call == "loss") for v in t["rna_median_lfc"])
            sig_all = (t["p_block_perm"] < args.alpha).all()
            sig_any = (t["p_block_perm"] < args.alpha).any()
            # a null shift only counts if it is significant AND material relative to the test shift
            test_size = float(np.mean(np.abs(t["rna_median_lfc"])))
            null_bad = [(x["p_block_perm"] < args.alpha) and (abs(x["rna_median_lfc"]) > args.null_frac * test_size)
                        for _, x in n.iterrows()]
            null_flat = not any(null_bad)
            row["null_flat"] = null_flat
            tail = "" if null_flat else "; BUT a null contrast shows a material shift"
            if not sign_ok:
                row["verdict"] = "NOT supported (direction disagrees)"
            elif sig_all and len(tests) == 1:
                row["verdict"] = "supported (single test contrast - add the UNC contrast to replicate)" + tail
            elif sig_all:
                row["verdict"] = ("CONFIRMED (replicated in all test contrasts)" if null_flat
                                  else "replicated in test contrasts" + tail)
            elif sig_any:
                row["verdict"] = "supported (right direction, significant in some test contrasts)" + tail
            else:
                row["verdict"] = "weak (right direction, not significant)" + tail
        verdict.append(row)
    vd = pd.DataFrame(verdict)
    vd.to_csv(os.path.join(args.outdir, "rnaseq_cnv_segment_verdict.tsv"), sep="\t", index=False)

    with open(os.path.join(args.outdir, "rnaseq_cnv_summary.txt"), "w") as fh:
        fh.write("test contrasts: " + ", ".join(l for l, _ in tests) + "\n")
        fh.write("null contrasts: " + (", ".join(l for l, _ in nulls) or "none") + "\n")
        fh.write(f"test one-sided in the Hi-C direction; null two-sided; alpha = {args.alpha}\n\n")
        fh.write(vd.round(4).to_string(index=False) + "\n\n")
        for lab, _ in tests:
            t = long[(long["contrast"] == lab)].dropna(subset=["rna_median_lfc"]) if "rna_median_lfc" in long else long.iloc[0:0]
            if len(t) >= 4:
                rho = spearmanr(t["hic_log2"], t["rna_median_lfc"]).correlation
                fh.write(f"Spearman rho, Hi-C vs RNA segment shift ({lab}): {rho:.2f}\n")
        fh.write("\nReading: CNV-driven expression = same sign as Hi-C in every test contrast, block-wise\n"
                 "significant, RNA/DNA ratio ~0.4-1, and NO shift in the null (UNC vs DMSO) contrasts.\n"
                 "RNA/DNA >> 1 suggests copy loss plus additional regulation (e.g. X-inactivation changes).\n"
                 "Caveat: Parental vs B10 libraries differ in kit; kit effects are gene-specific, not positional.\n")
    print("\n" + open(os.path.join(args.outdir, "rnaseq_cnv_summary.txt")).read())

    # ---- plot: one panel per contrast ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ref = per_gene[tests[0][0]]
    order = sorted([c for c in pd.unique(ref["chrom"]) if (ref["chrom"] == c).sum() > args.window], key=chrom_key)
    offsets, mids, cum = {}, [], 0
    for c in order:
        L = HG19.get(c, int(ref.loc[ref["chrom"] == c, "end"].max()))
        offsets[c] = cum
        mids.append(cum + L / 2)
        cum += L
    fig, axes = plt.subplots(len(contrasts), 1, figsize=(15, 2.6 * len(contrasts)), sharex=True, facecolor=SURFACE)
    axes = np.atleast_1d(axes)
    n_test = 0
    for ax, (lab, _, kind) in zip(axes, contrasts):
        m = per_gene[lab]
        pt_color = SERIES[n_test % len(SERIES)] if kind == "test" else NULL_COLOR
        n_test += kind == "test"
        ax.set_facecolor(SURFACE)
        for r in altered.itertuples():
            if r.chrom in offsets:
                ax.axvspan(offsets[r.chrom] + r.start, offsets[r.chrom] + r.end, lw=0, zorder=0,
                           color=GAIN_TINT if r.call == "gain" else LOSS_TINT,
                           alpha=GAIN_ALPHA if r.call == "gain" else LOSS_ALPHA)
        for c in order:
            g = m[m["chrom"] == c]
            if g.empty:
                continue
            x = offsets[c] + g["mid"].to_numpy()
            ax.scatter(x, np.clip(g["lfc_c"], -2, 2), s=2, color=pt_color, alpha=0.35, linewidths=0, rasterized=True)
            roll = g["lfc_c"].rolling(args.window, center=True, min_periods=args.window // 2).median()
            ax.plot(x, roll, color=INK, lw=1)
            ax.axvline(offsets[c], color="#d9d8d3", lw=0.6, zorder=0)
        ax.axhline(0, color=INK_2, lw=0.8)
        ax.set_ylim(-1.5, 1.5)
        ax.set_ylabel("log2FC (centred)", color=INK_2, fontsize=9)
        what = "test: B10 vs Parental" if kind == "test" else "null: same cells, UNC vs DMSO - should be flat"
        ax.set_title(f"{lab}  ({what})", loc="left", fontsize=10, color=INK)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(colors=INK_2, labelsize=8)
    axes[-1].set_xlim(0, cum)
    axes[-1].set_xticks(mids)
    axes[-1].set_xticklabels([c.replace("chr", "") for c in order], fontsize=8, color=INK_2)
    fig.suptitle("RNA-seq log2FC (genes) - black: rolling median; shading: Hi-C CNV calls (violet loss, yellow gain)",
                 x=0.01, ha="left", fontsize=12, color=INK)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(args.outdir, f"rnaseq_cnv_genome.{ext}"), dpi=200, facecolor=SURFACE)


if __name__ == "__main__":
    main()
