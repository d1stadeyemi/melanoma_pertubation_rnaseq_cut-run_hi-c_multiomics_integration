"""
cut_run_pc1_integration.py
================================================================================
The mechanistic spine: does H3K27me3 revert on the SAME loci where the
compartment reverts during constricted migration, and does EZH2i (UNC) erosion
concentrate at the migration-acquired compartment switches?

This first build leads with the headline overlap you asked for:

    csaw INTERACTION-DOWN regions      x      compartment migration-specific reversal
    (loci where UNC strips H3K27me3            (the 437-bin set from
     MORE in B10 than in parental;             compartment_reversal_analysis.py:
     the migration-gated PRC2 set)              reverted in B10 AND parental Stable)

INPUTS
    --reversal     Reversal_classified.txt  (compartment_reversal_analysis.py output:
                   coordinate, _chrom, mig_switch, reversal_status, punc_switch,
                   migration_specific_reversal, dPC1_migration, dPC1_b10_treat, ...)
    --marks_delta  marks_delta_250kb.txt    (quantify_marks_at_compartment_bins.R:
                   coordinate, H3K27me3_migration/_b10_treat/_par_treat, ...)
    --csaw_sig     <tag>_<norm>_significant.csv  (csaw.R interaction result; regions
                   with chr,start,end,logFC,FDR,direction). Overlapped onto the
                   250 kb compartment-bin grid.
    --gene_list    GeneListTSShg19_cleaned.txt   (coordinate, gene_name) for export.

WHAT IT COMPUTES  (effect sizes first; bin-level p-values are PSEUDOREP_, not for inference)
    1. The 2x2 headline: compartment reversal (Y/N) x H3K27me3 erosion (Y/N).
       "Erosion" is defined two ways and BOTH are reported:
         - csaw-down : the bin overlaps an FDR-significant interaction region with
                       logFC<0 (B10 loses more H3K27me3 than parental under UNC).
         - delta     : H3K27me3_b10_treat < 0 (UNC lowered H3K27me3 in B10), the
                       continuous sign from the mark-delta matrix.
       Reported: Haldane odds ratio, P(erosion|reversal) vs P(erosion|not), and a
       CHROMOSOME block-bootstrap 95% CI for each (resample whole chromosomes; the
       same pseudoreplication discipline as the compartment analysis).
    2. Continuous concentration: Cliff's delta of dH3K27me3_b10_treat at reversal
       bins vs the rest (sign = lt-gt; positive = reversal bins are MORE eroded),
       block-bootstrapped.
    3. dH3K27me3 vs dPC1 coupling, both on the migration axis and the UNC-reversal
       axis (Spearman + block-bootstrap).
    4. Gene lists: the RESPONSIVE set (compartment reverts AND H3K27me3 erodes) plus
       the four cross-tab cells and the background universe.

USAGE
    python track1_integrate.py \
        --reversal    results/exploratory/reversal/Reversal_classified.txt \
        --marks_delta results/exploratory/mark_signal_250kb/marks_delta_250kb.txt \
        --csaw_sig    results/exploratory/diff_analysis/csaw/H3K27me3_interaction/H3K27me3_interaction_libsize_significant.csv \
        --csaw_all    results/exploratory/diff_analysis/csaw/H3K27me3_interaction/H3K27me3_interaction_libsize_all_regions.csv \
        --gene_list   data/processed/GeneListTSShg19_cleaned.txt \
        --outdir      results/exploratory/track1
================================================================================
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import spearmanr, fisher_exact


COL = {
    "revert":  "#2A9D8F",   # reverted / responsive (emerald-teal)
    "stay":    "#9BA7B0",   # cool grey
    "erode":   "#B23A48",   # H3K27me3 loss (warm red)
    "both":    "#6A4C93",   # responsive overlap (distinct purple)
    "A":       "#B23A48",
    "B":       "#1F6F8B",
    "accent":  "#E08A3C",
}


# ============================================================================ #
# coordinate parsing
# ============================================================================ #
def parse_coord(series):
    """Parse 'chr1:1-250000' (or '_'/'-' separated) into chrom/start/end arrays."""
    chrom, start, end = [], [], []
    pat = re.compile(r"(chr[0-9XYMxym]+)[:_\- ]+([0-9]+)[:_\- ]+([0-9]+)")
    for s in series.astype(str):
        m = pat.search(s)
        if m:
            chrom.append(m.group(1)); start.append(int(m.group(2))); end.append(int(m.group(3)))
        else:
            chrom.append(None); start.append(np.nan); end.append(np.nan)
    return (np.array(chrom, dtype=object),
            np.array(start, dtype=float), np.array(end, dtype=float))


def chrom_of(coord):
    m = re.search(r"(chr[0-9XYMxym]+)", str(coord))
    return m.group(1) if m else str(coord).split(":")[0]


# ============================================================================ #
# loaders
# ============================================================================ #
def load_reversal(path, coord_col):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if coord_col not in df.columns:
        sys.exit(f"[ERROR] '{coord_col}' not in {path}. Columns: {list(df.columns)[:30]}")
    df[coord_col] = df[coord_col].astype(str).str.strip()
    ch, st, en = parse_coord(df[coord_col])
    df["_chr"], df["_start"], df["_end"] = ch, st, en
    if "_chrom" not in df.columns:
        df["_chrom"] = df[coord_col].map(chrom_of)
    return df


def load_marks_delta(path, coord_col):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if coord_col not in df.columns:
        sys.exit(f"[ERROR] '{coord_col}' not in {path}. Columns: {list(df.columns)[:30]}")
    df[coord_col] = df[coord_col].astype(str).str.strip()
    return df


def load_csaw_sig(path):
    """csaw significant.csv -> regions with chr/start/end/logFC/FDR (robust to col names)."""
    df = pd.read_csv(path)
    cols = {c.lower(): c for c in df.columns}
    def pick(*names):
        for n in names:
            if n in cols:
                return cols[n]
        return None
    c_chr = pick("chr", "seqnames", "chrom")
    c_st  = pick("start"); c_en = pick("end")
    c_lfc = pick("logfc", "log2fc"); c_fdr = pick("fdr", "padj", "adj.p.val")
    if not all([c_chr, c_st, c_en, c_lfc]):
        sys.exit(f"[ERROR] csaw_sig needs chr/start/end/logFC. Columns: {list(df.columns)}")
    out = pd.DataFrame({
        "chr":   df[c_chr].astype(str).str.strip(),
        "start": pd.to_numeric(df[c_st], errors="coerce"),
        "end":   pd.to_numeric(df[c_en], errors="coerce"),
        "logFC": pd.to_numeric(df[c_lfc], errors="coerce"),
    })
    out["FDR"] = pd.to_numeric(df[c_fdr], errors="coerce") if c_fdr else np.nan
    out = out.dropna(subset=["start", "end", "logFC"])
    return out


# ============================================================================ #
# overlap csaw regions onto the 250 kb compartment-bin grid
# ============================================================================ #
def overlap_flags(bins, regions):
    """For each bin, flag whether it overlaps any significant region, and whether
    any overlapping region is down (logFC<0) or up (logFC>0). Bins and regions are
    1-based inclusive; overlap = (bin_start <= reg_end) & (bin_end >= reg_start)."""
    n = len(bins)
    sig  = np.zeros(n, dtype=bool)
    down = np.zeros(n, dtype=bool)
    up   = np.zeros(n, dtype=bool)
    bstart = bins["_start"].to_numpy(float)
    bend   = bins["_end"].to_numpy(float)
    bchr   = bins["_chr"].to_numpy(object)
    for chrom, reg in regions.groupby("chr"):
        bmask = (bchr == chrom)
        if not bmask.any():
            continue
        bi = np.where(bmask)[0]
        rs = reg["start"].to_numpy(float); re_ = reg["end"].to_numpy(float)
        rl = reg["logFC"].to_numpy(float)
        for k in bi:
            ov = (bstart[k] <= re_) & (bend[k] >= rs)
            if ov.any():
                sig[k] = True
                lf = rl[ov]
                if (lf < 0).any(): down[k] = True
                if (lf > 0).any(): up[k]   = True
    return sig, down, up


# ============================================================================ #
# effect-size helpers
# ============================================================================ #
def haldane_or(a, b, c, d):
    a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    return (a * d) / (b * c)


def cliffs_delta(a, b):
    """Cliff's delta with project sign convention (lt - gt): positive => values in
    `a` tend to be LESS than values in `b`. Vectorised O((n+m) log m)."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    na, nb = len(a), len(b)
    if na == 0 or nb == 0:
        return np.nan
    bs = np.sort(b)
    # for each a: #(b > a) and #(b < a)
    gt = nb - np.searchsorted(bs, a, side="right")   # b strictly greater than a
    lt = np.searchsorted(bs, a, side="left")          # b strictly less than a
    # a<b pairs = sum(gt) ; a>b pairs = sum(lt)
    return (gt.sum() - lt.sum()) / (na * nb)          # (a<b) - (a>b), normalised


def spearman_safe(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 10:
        return np.nan
    if np.all(x[m] == x[m][0]) or np.all(y[m] == y[m][0]):
        return np.nan
    return spearmanr(x[m], y[m]).correlation


# ============================================================================ #
# chromosome block-bootstrap
# ============================================================================ #
def block_bootstrap(stat_fn, chrom, n_boot=2000, seed=0):
    """stat_fn(idx)->float on the index subset; resample whole chromosomes."""
    uniq = pd.unique(chrom)
    blocks = {c: np.where(chrom == c)[0] for c in uniq}
    nb = len(uniq)
    obs = stat_fn(np.arange(len(chrom)))
    rng = np.random.default_rng(seed)
    boot = np.full(n_boot, np.nan)
    for it in range(n_boot):
        pick = rng.choice(uniq, size=nb, replace=True)
        idx = np.concatenate([blocks[c] for c in pick])
        boot[it] = stat_fn(idx)
    lo, hi = np.nanpercentile(boot, [2.5, 97.5])
    return obs, lo, hi


# ============================================================================ #
# build the integrated table
# ============================================================================ #
def build_table(args):
    rev = load_reversal(args.reversal, args.coord_col)
    mk  = load_marks_delta(args.marks_delta, args.coord_col)

    # 1:1 inner join on coordinate (the reason exact bins matter)
    n_rev, n_mk = len(rev), len(mk)
    df = rev.merge(mk, on=args.coord_col, how="inner", suffixes=("", "_mk"))
    print(f"[INFO] reversal bins {n_rev:,} | marks_delta bins {n_mk:,} | joined {len(df):,}")
    if len(df) == 0:
        sys.exit("[ERROR] zero bins after coordinate join — check that both files share the "
                 "same coordinate strings (same compartment table).")

    # csaw interaction overlap onto the grid
    reg = load_csaw_sig(args.csaw_sig)
    print(f"[INFO] csaw significant regions: {len(reg):,} "
          f"(down={int((reg['logFC']<0).sum())}, up={int((reg['logFC']>0).sum())})")
    sig, down, up = overlap_flags(df, reg)
    df["h3k27me3_int_sig"]  = sig
    df["h3k27me3_int_down"] = down
    df["h3k27me3_int_up"]   = up
    print(f"[INFO] bins overlapping a sig region: {int(sig.sum()):,}  "
          f"(down={int(down.sum())}, up={int(up.sum())})")

    # OPTIONAL: csaw TESTED universe (from all_regions) -> restrict the csaw-down
    # cross-tab denominator so untested bins don't swell the 'not down' cell.
    if args.csaw_all:
        allreg = load_csaw_sig(args.csaw_all)
        tested, _, _ = overlap_flags(df, allreg)
        df["h3k27me3_tested"] = tested
        print(f"[INFO] csaw all_regions: {len(allreg):,} tested regions | "
              f"{int(tested.sum()):,}/{len(df):,} bins in the tested universe "
              "(csaw-down cross-tab restricted to these).")
    else:
        df["h3k27me3_tested"] = True
        print("[INFO] --csaw_all not given: csaw-down cross-tab uses ALL bins "
              "(untested bins count as 'not down' -> OR is mildly anti-conservative).")

    # compartment reversal flag (headline 437-set)
    if "migration_specific_reversal" in df.columns:
        df["comp_reversal"] = df["migration_specific_reversal"].astype(bool)
    else:
        df["comp_reversal"] = (df.get("reversal_status", "") == "reverted")
        print("[WARN] 'migration_specific_reversal' not found; using reversal_status=='reverted'.")

    # H3K27me3 erosion flags
    db10 = pd.to_numeric(df.get("H3K27me3_b10_treat", np.nan), errors="coerce")
    df["H3K27me3_b10_treat_num"] = db10
    df["erode_csaw"]  = df["h3k27me3_int_down"]        # FDR-significant B10-preferential loss
    df["erode_delta"] = db10 < 0                        # UNC lowered H3K27me3 in B10 (sign)
    return df, reg


# ============================================================================ #
# 2x2 cross-tab with block-bootstrap
# ============================================================================ #
def crosstab_block(df, erosion_col, label, args):
    rev = df["comp_reversal"].to_numpy(bool)
    ero = df[erosion_col].to_numpy(bool)
    chrom = df["_chrom"].to_numpy(object)

    def cells(idx):
        r = rev[idx]; e = ero[idx]
        a = int(np.sum(r & e)); b = int(np.sum(r & ~e))
        c = int(np.sum(~r & e)); d = int(np.sum(~r & ~e))
        return a, b, c, d

    a, b, c, d = cells(np.arange(len(df)))
    obs_or = haldane_or(a, b, c, d)
    p_rev = a / (a + b) if (a + b) else np.nan          # P(erosion | reversal)
    p_non = c / (c + d) if (c + d) else np.nan          # P(erosion | not)

    or_ci  = block_bootstrap(lambda i: haldane_or(*cells(i)), chrom, args.n_boot, args.boot_seed)
    prev_ci = block_bootstrap(
        lambda i: (lambda x: x[0]/(x[0]+x[1]) if (x[0]+x[1]) else np.nan)(cells(i)),
        chrom, args.n_boot, args.boot_seed)
    pnon_ci = block_bootstrap(
        lambda i: (lambda x: x[2]/(x[2]+x[3]) if (x[2]+x[3]) else np.nan)(cells(i)),
        chrom, args.n_boot, args.boot_seed)
    try:
        _, fisher_p = fisher_exact([[a, b], [c, d]])
    except Exception:
        fisher_p = np.nan

    print(f"\n[CROSSTAB:{label}]  erosion = {erosion_col}")
    print(f"   reversal&erode={a}  reversal&stable={b}  other&erode={c}  other&stable={d}")
    print(f"   P(erosion|reversal)={p_rev:.3f} [{prev_ci[1]:.3f},{prev_ci[2]:.3f}]  "
          f"P(erosion|other)={p_non:.3f} [{pnon_ci[1]:.3f},{pnon_ci[2]:.3f}]")
    print(f"   odds ratio={obs_or:.2f} [{or_ci[1]:.2f},{or_ci[2]:.2f}]   "
          f"PSEUDOREP_fisher_p={fisher_p:.2e}")
    return dict(label=label, erosion_col=erosion_col, a=a, b=b, c=c, d=d,
                OR=or_ci, p_rev=prev_ci, p_non=pnon_ci, fisher_p=fisher_p)


# ============================================================================ #
# continuous concentration + coupling
# ============================================================================ #
def concentration_block(df, args):
    rev = df["comp_reversal"].to_numpy(bool)
    db10 = df["H3K27me3_b10_treat_num"].to_numpy(float)
    chrom = df["_chrom"].to_numpy(object)

    def cd(idx):
        return cliffs_delta(db10[idx][rev[idx]], db10[idx][~rev[idx]])

    obs, lo, hi = block_bootstrap(cd, chrom, args.n_boot, args.boot_seed)
    med_rev = np.nanmedian(db10[rev]); med_non = np.nanmedian(db10[~rev])
    print(f"\n[CONCENTRATION] dH3K27me3_b10_treat at reversal vs non-reversal bins")
    print(f"   median(reversal)={med_rev:+.3f}  median(non-reversal)={med_non:+.3f}")
    print(f"   Cliff's delta (lt-gt; + = reversal more eroded)={obs:+.3f} [{lo:+.3f},{hi:+.3f}]")
    return dict(cliffs=(obs, lo, hi), med_rev=med_rev, med_non=med_non)


def coupling_block(df, args):
    chrom = df["_chrom"].to_numpy(object)
    out = {}
    pairs = [
        ("migration", "H3K27me3_migration", "dPC1_migration"),
        ("reversal",  "H3K27me3_b10_treat", "dPC1_b10_treat"),
    ]
    print("\n[COUPLING] dH3K27me3 vs dPC1 (Spearman, block-bootstrap)")
    for name, mcol, pcol in pairs:
        if mcol not in df.columns or pcol not in df.columns:
            print(f"   [skip] {name}: need {mcol} & {pcol}")
            out[name] = None
            continue
        x = pd.to_numeric(df[mcol], errors="coerce").to_numpy(float)
        y = pd.to_numeric(df[pcol], errors="coerce").to_numpy(float)
        obs, lo, hi = block_bootstrap(lambda i: spearman_safe(x[i], y[i]),
                                      chrom, args.n_boot, args.boot_seed)
        print(f"   {name:9s} rho(dH3K27me3, dPC1)={obs:+.3f} [{lo:+.3f},{hi:+.3f}]")
        out[name] = (obs, lo, hi)
    return out


# ============================================================================ #
# figures
# ============================================================================ #
def fig_crosstab(ct, outdir):
    a, b, c, d = ct["a"], ct["b"], ct["c"], ct["d"]
    M = np.array([[a, b], [c, d]])
    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    im = ax.imshow(M, cmap="Purples", aspect="auto")
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f"{M[i, j]:,}", ha="center", va="center",
                    fontsize=14, fontweight="bold",
                    color="white" if M[i, j] > M.max() * 0.6 else "black")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["H3K27me3\neroded", "H3K27me3\nstable"])
    ax.set_yticks([0, 1]); ax.set_yticklabels(["compartment\nreversal", "non-reversal\nbins"])
    o, olo, ohi = ct["OR"]
    ax.set_title(f"Compartment reversal x H3K27me3 erosion ({ct['label']})\n"
                 f"OR={o:.1f} [{olo:.1f}, {ohi:.1f}]   (chrom block-bootstrap)",
                 fontweight="bold", fontsize=10)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="bin count")
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"Fig_crosstab_{ct['label']}.{ext}"),
                    bbox_inches="tight", dpi=600)
    plt.close(fig)


def fig_coupling(df, coup, outdir):
    pairs = [("migration", "H3K27me3_migration", "dPC1_migration",
              "Migration axis"),
             ("reversal", "H3K27me3_b10_treat", "dPC1_b10_treat",
              "UNC reversal axis (B10)")]
    have = [p for p in pairs if p[1] in df.columns and p[2] in df.columns]
    if not have:
        return
    fig, axes = plt.subplots(1, len(have), figsize=(5.8 * len(have), 5.2), squeeze=False)
    for ax, (name, mcol, pcol, title) in zip(axes[0], have):
        x = pd.to_numeric(df[mcol], errors="coerce")
        y = pd.to_numeric(df[pcol], errors="coerce")
        rev = df["comp_reversal"].to_numpy(bool)
        ax.scatter(x[~rev], y[~rev], s=7, alpha=0.20, color=COL["stay"],
                   edgecolors="none", label="other bins")
        ax.scatter(x[rev], y[rev], s=14, alpha=0.75, color=COL["both"],
                   edgecolors="none", label="compartment reversal")
        ax.axhline(0, color="#888", lw=0.7); ax.axvline(0, color="#888", lw=0.7)
        rr = coup.get(name)
        sub = f"rho={rr[0]:+.2f} [{rr[1]:+.2f}, {rr[2]:+.2f}]" if rr else ""
        ax.set_xlabel(f"\u0394H3K27me3 ({name})")
        ax.set_ylabel(f"\u0394PC1 ({name})")
        ax.set_title(f"{title}\n{sub}", fontweight="bold", fontsize=10)
        ax.legend(frameon=False, fontsize=8, loc="best")
    fig.suptitle("H3K27me3 - compartment coupling (block-bootstrap rho)",
                 fontweight="bold", fontsize=11, y=1.03)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"Fig_coupling.{ext}"),
                    bbox_inches="tight", dpi=600)
    plt.close(fig)


def fig_concentration(df, conc, outdir):
    rev = df["comp_reversal"].to_numpy(bool)
    db10 = df["H3K27me3_b10_treat_num"].to_numpy(float)
    a = db10[rev & np.isfinite(db10)]; b = db10[~rev & np.isfinite(db10)]
    if len(a) < 5 or len(b) < 5:
        return
    fig, ax = plt.subplots(figsize=(6.0, 4.8))
    parts = ax.violinplot([b, a], showmedians=True, widths=0.8)
    for pc, col in zip(parts["bodies"], [COL["stay"], COL["both"]]):
        pc.set_facecolor(col); pc.set_alpha(0.6)
    ax.axhline(0, color="#888", lw=0.8, ls="--")
    ax.set_xticks([1, 2]); ax.set_xticklabels(["non-reversal\nbins", "compartment\nreversal"])
    ax.set_ylabel("\u0394H3K27me3  (B10+UNC \u2212 B10)")
    cobs, clo, chi = conc["cliffs"]
    ax.set_title("UNC-driven H3K27me3 change at reversal bins\n"
                 f"Cliff's \u03b4={cobs:+.2f} [{clo:+.2f}, {chi:+.2f}]  "
                 "(+ = reversal more eroded)", fontweight="bold", fontsize=10)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"Fig_concentration.{ext}"),
                    bbox_inches="tight", dpi=600)
    plt.close(fig)


def fig_forest(ct_csaw, ct_delta, conc, coup, outdir):
    rows = []
    for ct in (ct_csaw, ct_delta):
        o, lo, hi = ct["OR"]
        rows.append((f"OR reversal x erosion ({ct['label']})", o, lo, hi, "or"))
    cobs, clo, chi = conc["cliffs"]
    rows.append(("Cliff's \u03b4 (erosion@reversal)", cobs, clo, chi, "d"))
    for name in ("migration", "reversal"):
        rr = coup.get(name)
        if rr:
            rows.append((f"\u03c1 \u0394H3K27me3~\u0394PC1 ({name})", rr[0], rr[1], rr[2], "rho"))
    fig, ax = plt.subplots(figsize=(8.2, 0.7 * len(rows) + 1.5))
    ypos = list(range(len(rows)))[::-1]
    for y, (lab, o, lo, hi, kind) in zip(ypos, rows):
        col = COL["both"] if kind == "or" else (COL["erode"] if kind == "d" else COL["B"])
        ax.errorbar(o, y, xerr=[[o - lo], [hi - o]], fmt="o", color=col,
                    ecolor="#333", capsize=4, ms=7, lw=1.4)
        ax.text(o, y + 0.18, f"{o:+.2f} [{lo:+.2f}, {hi:+.2f}]",
                ha="center", va="bottom", fontsize=8.5)
    ax.axvline(0, color="#888", lw=0.8, ls="--")
    ax.axvline(1, color="#bbb", lw=0.8, ls=":")
    ax.set_yticks(ypos); ax.set_yticklabels([r[0] for r in rows], fontsize=9)
    ax.set_xlabel("effect size (95% chromosome block-bootstrap CI)")
    ax.set_title("Track 1 headline effect sizes\n(OR vs 1 dotted; \u03b4/\u03c1 vs 0 dashed)",
                 fontweight="bold", fontsize=11)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"Fig_track1_forest.{ext}"),
                    bbox_inches="tight", dpi=600)
    plt.close(fig)


# ============================================================================ #
# gene export
# ============================================================================ #
def export_genes(df, args, outdir):
    gdir = os.path.join(outdir, "genelists")
    os.makedirs(gdir, exist_ok=True)
    cc, gn = args.coord_col, args.gene_name_col
    genes = pd.read_csv(args.gene_list, sep="\t", low_memory=False)
    for c in (cc, gn):
        if c not in genes.columns:
            sys.exit(f"[ERROR] '{c}' not in gene_list. Columns: {list(genes.columns)}")
    genes[cc] = genes[cc].astype(str).str.strip()
    genes[gn] = genes[gn].astype(str).str.strip()

    keep = [cc, "comp_reversal", "erode_csaw", "erode_delta",
            "h3k27me3_int_down", "reversal_status"]
    keep = [k for k in keep if k in df.columns]
    g = genes.merge(df[keep], on=cc, how="inner")

    def write(mask, name):
        syms = sorted(set(s for s in g.loc[mask, gn] if s and s.lower() != "nan"))
        with open(os.path.join(gdir, f"{name}.txt"), "w") as fh:
            fh.write("\n".join(syms) + ("\n" if syms else ""))
        print(f"   gene set [{name:42s}] {len(syms):5d} -> {name}.txt")
        return len(syms)

    write(pd.Series(True, index=g.index), "background_universe")
    write(g["comp_reversal"] == True, "compartment_reversal")
    write(g["erode_csaw"] == True, "h3k27me3_eroded_csaw")
    # the headline RESPONSIVE sets (compartment reverts AND H3K27me3 erodes)
    write((g["comp_reversal"] == True) & (g["erode_csaw"] == True),
          "responsive_reversal_and_csawdown")
    if "erode_delta" in g.columns:
        write((g["comp_reversal"] == True) & (g["erode_delta"] == True),
              "responsive_reversal_and_deltaloss")
    # full per-gene assignment table
    g.to_csv(os.path.join(gdir, "track1_gene_assignments.txt"), sep="\t", index=False)
    print(f"   per-gene assignments -> track1_gene_assignments.txt")


# ============================================================================ #
# stats file
# ============================================================================ #
def write_stats(df, ct_csaw, ct_delta, conc, coup, outdir):
    path = os.path.join(outdir, "Track1_stats.txt")
    with open(path, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"joined_bins\t{len(df)}\n")
        fh.write(f"n_chromosomes\t{df['_chrom'].nunique()}\n")
        fh.write(f"compartment_reversal_bins\t{int(df['comp_reversal'].sum())}\n")
        fh.write(f"h3k27me3_int_down_bins\t{int(df['h3k27me3_int_down'].sum())}\n")
        if "h3k27me3_tested" in df.columns:
            fh.write(f"csaw_tested_bins\t{int(df['h3k27me3_tested'].sum())}\n")
        fh.write("# === HEADLINE effect sizes (chromosome block-bootstrap 95% CI) ===\n")
        for ct in (ct_csaw, ct_delta):
            for k in ("a", "b", "c", "d"):
                fh.write(f"crosstab_{ct['label']}_{k}\t{ct[k]}\n")
            for nm in ("OR", "p_rev", "p_non"):
                o, lo, hi = ct[nm]
                fh.write(f"crosstab_{ct['label']}_{nm}\t{o:.5f}\n")
                fh.write(f"crosstab_{ct['label']}_{nm}_CI95_lo\t{lo:.5f}\n")
                fh.write(f"crosstab_{ct['label']}_{nm}_CI95_hi\t{hi:.5f}\n")
        c0, clo, chi = conc["cliffs"]
        fh.write(f"concentration_cliffs_delta\t{c0:.5f}\n")
        fh.write(f"concentration_cliffs_delta_CI95_lo\t{clo:.5f}\n")
        fh.write(f"concentration_cliffs_delta_CI95_hi\t{chi:.5f}\n")
        for name in ("migration", "reversal"):
            rr = coup.get(name)
            if rr:
                o, lo, hi = rr
                fh.write(f"coupling_{name}_rho\t{o:.5f}\n")
                fh.write(f"coupling_{name}_rho_CI95_lo\t{lo:.5f}\n")
                fh.write(f"coupling_{name}_rho_CI95_hi\t{hi:.5f}\n")
        fh.write("# === PSEUDOREPLICATED bin-level p-values - NOT FOR INFERENCE ===\n")
        fh.write(f"PSEUDOREP_crosstab_csaw_fisher_p\t{ct_csaw['fisher_p']:.3e}\n")
        fh.write(f"PSEUDOREP_crosstab_delta_fisher_p\t{ct_delta['fisher_p']:.3e}\n")
    print(f"[INFO] stats -> {path}")


# ============================================================================ #
# CLI + main
# ============================================================================ #
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--reversal", required=True)
    p.add_argument("--marks_delta", required=True)
    p.add_argument("--csaw_sig", required=True)
    p.add_argument("--csaw_all", default=None,
                   help="csaw <tag>_<norm>_all_regions.csv. If given, the csaw-down "
                        "cross-tab is restricted to bins csaw actually TESTED (overlap "
                        "any region here), removing the anti-conservative inflation from "
                        "untested bins defaulting to 'not down'. The delta cross-tab and "
                        "Cliff's delta are unaffected.")
    p.add_argument("--gene_list", default=None)
    p.add_argument("--coord_col", default="coordinate")
    p.add_argument("--gene_name_col", default="gene_name")
    p.add_argument("--outdir", default="results/exploratory/track1")
    p.add_argument("--n_boot", type=int, default=2000)
    p.add_argument("--boot_seed", type=int, default=0)
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df, reg = build_table(args)
    df.to_csv(os.path.join(args.outdir, "Track1_classified.txt"), sep="\t", index=False)

    print("\n[INFO] Headline cross-tabs (chromosome block-bootstrap)...")
    csaw_df = df[df["h3k27me3_tested"]].copy() if args.csaw_all else df
    if args.csaw_all:
        print(f"[INFO] csaw-down cross-tab restricted to {len(csaw_df):,} tested bins.")
    ct_csaw  = crosstab_block(csaw_df, "erode_csaw",  "csawdown", args)
    ct_delta = crosstab_block(df, "erode_delta", "deltaloss", args)
    conc = concentration_block(df, args)
    coup = coupling_block(df, args)

    print("\n[INFO] Figures...")
    fig_crosstab(ct_csaw, args.outdir)
    fig_crosstab(ct_delta, args.outdir)
    fig_coupling(df, coup, args.outdir)
    fig_concentration(df, conc, args.outdir)
    fig_forest(ct_csaw, ct_delta, conc, coup, args.outdir)

    write_stats(df, ct_csaw, ct_delta, conc, coup, args.outdir)

    if args.gene_list:
        print("\n[INFO] Gene lists...")
        export_genes(df, args, args.outdir)
    else:
        print("[INFO] --gene_list not given; skipping gene export.")

    print(f"\n[INFO] Done -> {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()