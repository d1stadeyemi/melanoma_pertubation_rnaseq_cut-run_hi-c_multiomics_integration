"""
=============================================================================
compartment_reversal_analysis.py
  Does UNC (EZH2i) treatment REVERSE migration-acquired compartment switches?
=============================================================================

ALL-DMSO BASELINE (single-experiment consistency)
  Everything is measured WITHIN the drug experiment:
    migration  = Bottom10-DMSO - Parental-DMSO
    B10 treat  = Bottom10-UNC  - Bottom10-DMSO
    par  treat = Parental-UNC  - Parental-DMSO
  Migration and the treatment deltas SHARE the DMSO baseline (opposite signs),
  which injects a regression-to-mean / shared-variable contribution into the
  reversal & control correlations. The artifact diagnostic + attractor-residual
  subset REMOVE this; the partial/residual numbers are the best effect, NOT
  the raw correlation. orig columns are not used anywhere.

INFERENCE NOTE (read before quoting any p-value)
  250 kb bins are NOT independent: compartments are megabase-scale, so adjacent
  bins are strongly autocorrelated and ~2,000 switcher bins really represent only
  a few hundred independent domains / ~23 chromosomes. The per-bin Fisher /
  Spearman p-values are therefore PSEUDOREPLICATED and anti-conservative. They are
  written to the stats file ONLY for completeness, prefixed PSEUDOREP_ and flagged
  not-for-inference. The headline uncertainty is the CHROMOSOME BLOCK-BOOTSTRAP
  95% CI (resamples whole chromosomes, preserving within-chromosome correlation),
  reported for: reversal rate (B10), anti-migration rate (parental control),
  excess reversal, odds ratio, and the attractor-removed residual Spearman rho.
  Lead with effect sizes + these CIs, not the bin-level p.
=============================================================================
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
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.patches import Patch
from scipy.stats import fisher_exact, spearmanr, rankdata


# ---- shared palette -------------------------------------------------------
COL = {"A": "#C0392B", "B": "#1A5276", "switch": "#E67E22",
       "revert": "#27AE60", "stay": "#95A5A6", "further": "#8E44AD",
       "AtoB": "#1A5276", "BtoA": "#C0392B", "Stable": "#D5D8DC"}
ARROW = "\u2192"


# =============================================================================
# Token parsing
#   Each directional change is encoded as a Both / Switch / Shift column triple.
#   parse_switch_token collapses that triple into {A_to_B, B_to_A, Stable}.
#   "Toward"/"shift" calls (sub-threshold movement that does not cross the A/B
#   line) fold into the same direction, so a token is a DIRECTION-of-change call,
#   not a strict sign flip.
# =============================================================================
def parse_switch_token(both, switch, shift):
    for v in (both, switch, shift):
        s = str(v).strip()
        if s in ("-", "", "nan", "None", "NaN"):
            continue
        if "BtoA" in s or "BothtoA" in s or "TowardA" in s or s == "toA":
            return "B_to_A"
        if "AtoB" in s or "BothtoB" in s or "TowardB" in s or s == "toB":
            return "A_to_B"
    return "Stable"


def classify_axis(df, both_c, switch_c, shift_c):
    """Apply parse_switch_token over a Both/Switch/Shift column triple."""
    return df.apply(lambda r: parse_switch_token(
        r.get(both_c, "-"), r.get(switch_c, "-"), r.get(shift_c, "-")), axis=1)


def chrom_of(coord):
    """Extract a chromosome label from a coordinate string for block-bootstrap
    blocking. Handles 'chr1:...', 'chr1-...', 'chr1_...', etc."""
    s = str(coord)
    m = re.search(r"(chr[0-9XYMxym]+)", s)
    if m:
        return m.group(1)
    for sep in (":", "-", "_", "\t", " "):
        if sep in s:
            return s.split(sep)[0]
    return s


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pc1", required=True)
    p.add_argument("--outdir", default="results/exploratory/reversal")
    p.add_argument("--gene_list", default=None)
    p.add_argument("--coord_col", default="coordinate")
    p.add_argument("--gene_name_col", default="gene_name")
    p.add_argument("--chrom_col", default=None,
                   help="Optional explicit chromosome column for block-bootstrap "
                        "blocking. If absent, chromosome is parsed from --coord_col.")
    # ALL-DMSO baseline: migration + treatment both anchored on the DMSO arm
    p.add_argument("--pc1_parental_dmso", default="Parental_DMSO")
    p.add_argument("--pc1_b10_dmso", default="Bottom10-DMSO")
    p.add_argument("--pc1_parental_unc", default="Parental_UNC")
    p.add_argument("--pc1_b10_unc", default="Bottom10-UNC")
    p.add_argument("--mig_cols", nargs=3,
                   default=["P->B10_Both_DMSO", "P->B10_Switch_DMSO", "P->B10_Shift_DMSO"],
                   metavar=("BOTH", "SWITCH", "SHIFT"))
    p.add_argument("--punc_cols", nargs=3,
                   default=["P_DMSO->UNC_Both", "P_DMSO->UNC_switch", "P_DMSO->UNC_shift"],
                   metavar=("BOTH", "SWITCH", "SHIFT"))
    p.add_argument("--b10unc_cols", nargs=3,
                   default=["B10_DMSO->UNC_Both", "B10_DMSO->UNC_switch", "B10_DMSO->UNC_shift"],
                   metavar=("BOTH", "SWITCH", "SHIFT"))
    p.add_argument("--include_stable", action="store_true")
    p.add_argument("--include_unresponsive", action="store_true")
    p.add_argument("--unresp_thresh", type=float, default=None,
                   help="|residual ΔPC1| threshold for the migration-specific subset "
                        "(attractor regressed out). Default: 66th pct of |residual|.")
    # block-bootstrap controls
    p.add_argument("--n_boot", type=int, default=2000,
                   help="Number of chromosome block-bootstrap resamples (default 2000).")
    p.add_argument("--boot_seed", type=int, default=0,
                   help="RNG seed for the block-bootstrap (reproducibility).")
    return p.parse_args()


# =============================================================================
# Build the per-bin reversal table
#   Two INDEPENDENT directional comparisons, each vs its OWN baseline:
#     mig_switch    : Parental-DMSO -> Bottom10-DMSO   (migration axis)
#     b10unc_switch : Bottom10-DMSO -> Bottom10-UNC    (treatment, FROM the B10 state)
#     punc_switch   : Parental-DMSO -> Parental-UNC    (treatment in parental = control)
#   reversal_status combines mig_switch + b10unc_switch for migration switchers:
#     unresponsive : UNC did not move the migrated bin (kept its new compartment)
#     reverted     : UNC moved it OPPOSITE to migration (back toward parental)
#     reinforced   : UNC moved it the SAME way (further); rare - a migrated bin is
#                    already in its new compartment, so this is mostly sub-threshold
# =============================================================================
def build_reversal_table(df, args):
    df["mig_switch"]    = classify_axis(df, *args.mig_cols)
    df["punc_switch"]   = classify_axis(df, *args.punc_cols)
    df["b10unc_switch"] = classify_axis(df, *args.b10unc_cols)

    par = pd.to_numeric(df[args.pc1_parental_dmso], errors="coerce")
    b10 = pd.to_numeric(df[args.pc1_b10_dmso], errors="coerce")
    df["baseline_class"] = np.where(
        (par > 0) & (b10 > 0), "Stable_A",
        np.where((par < 0) & (b10 < 0), "Stable_B",
                 np.where(df["mig_switch"] == "A_to_B", "A_to_B",
                          np.where(df["mig_switch"] == "B_to_A", "B_to_A", "Other"))))

    def reversal_status(r):
        m, t = r["mig_switch"], r["b10unc_switch"]
        if m == "Stable":
            return "not_switcher"
        if t == "Stable":
            return "unresponsive"
        if (m == "A_to_B" and t == "B_to_A") or (m == "B_to_A" and t == "A_to_B"):
            return "reverted"
        if m == t:
            return "reinforced"
        return "other"
    df["reversal_status"] = df.apply(reversal_status, axis=1)
    df["parental_unc_moved"] = df["punc_switch"] != "Stable"
    df["categorical_specific_reversal"] = (
        (df["reversal_status"] == "reverted") & 
        (df["punc_switch"] == "Stable")
    ).fillna(False)

    pu  = pd.to_numeric(df[args.pc1_parental_unc], errors="coerce")
    pdm = pd.to_numeric(df[args.pc1_parental_dmso], errors="coerce")
    bu  = pd.to_numeric(df[args.pc1_b10_unc], errors="coerce")
    bdm = pd.to_numeric(df[args.pc1_b10_dmso], errors="coerce")
    df["dPC1_migration"] = b10 - par      # Bottom10-DMSO - Parental-DMSO
    df["dPC1_b10_treat"] = bu  - bdm       # Bottom10-UNC  - Bottom10-DMSO
    df["dPC1_par_treat"] = pu  - pdm       # Parental-UNC  - Parental-DMSO

    # chromosome label for block-bootstrap blocking
    if args.chrom_col and args.chrom_col in df.columns:
        df["_chrom"] = df[args.chrom_col].astype(str).str.strip()
    else:
        df["_chrom"] = df[args.coord_col].astype(str).map(chrom_of)
    return df


def ctrl_reverted_mask(sw):
    """Parental control: UNC moved the bin OPPOSITE to its migration direction.
    NOT a literal reversal (the parental bin never migrated) - a directional-
    coincidence statistic whose near-zero rate is the specificity proof."""
    m, t = sw["mig_switch"], sw["punc_switch"]
    return (((m == "A_to_B") & (t == "B_to_A")) |
            ((m == "B_to_A") & (t == "A_to_B")))


# =============================================================================
# FIGURE 1 - token-based braided alluvial (switchers only)
#   Consistent with the heatmap/bars: states come from the SAME tokens, not raw
#   PC1 sign. Each band is one bin's compartment PATH across the three stages; a
#   band that leaves its parental compartment and RETURNS (A->B->A or B->A->B) is
#   a reversal (green). Band widths == heatmap counts.
# =============================================================================
def _alluvial_ribbon(ax, x0, y0a, y0b, x1, y1a, y1b, color, alpha):
    """One cubic-bezier flow band from a vertical interval at x0 to one at x1."""
    xm = (x0 + x1) / 2
    verts = [(x0, y0a), (xm, y0a), (xm, y1a), (x1, y1a),
             (x1, y1b), (xm, y1b), (xm, y0b), (x0, y0b), (x0, y0a)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    ax.add_patch(patches.PathPatch(Path(verts, codes), facecolor=color,
                 edgecolor="none", alpha=alpha, zorder=2))


def fig_alluvial(df, outdir, args):
    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])].copy()
    if len(sw) == 0:
        print("  [WARN] no migration switchers; skipping alluvial."); return

    # derive A/B state at each stage from the tokens
    #   mig_switch    -> (Parental, Bottom10) compartments
    #   b10unc_switch -> Bottom10+UNC compartment (Stable = unchanged)
    def _states(r):
        s1, s2 = ("A", "B") if r["mig_switch"] == "A_to_B" else ("B", "A")
        t = r["b10unc_switch"]
        s3 = "B" if t == "A_to_B" else "A" if t == "B_to_A" else s2
        return pd.Series([s1, s2, s3])
    sw[["s1", "s2", "s3"]] = sw.apply(_states, axis=1)
    N = len(sw)

    # unique paths (flows); s3==s1 means it returned to parental = reversal
    flows = (sw.groupby(["s1", "s2", "s3"]).size().reset_index(name="n")
               .sort_values(["s1", "s2", "s3"]).reset_index(drop=True))
    flows["type"] = np.where(flows["s3"] == flows["s1"], "reverted", "stayed")

    state_cols = ["s1", "s2", "s3"]
    stage_names = ["Parental", "Bottom10", "Bottom10+UNC"]
    gap = 0.06
    order = list(flows.index)   # fixed flow order keeps each band contiguous

    # vertical layout per stage: A block at bottom, B block on top
    yint, node_y = {}, {}
    for si, sc in enumerate(state_cols):
        a_h = flows.loc[flows[sc] == "A", "n"].sum() / N
        b_h = flows.loc[flows[sc] == "B", "n"].sum() / N
        node_y[(si, "A")] = (0.0, a_h)
        node_y[(si, "B")] = (a_h + gap, a_h + gap + b_h)
        ya = 0.0
        for fi in order:
            if flows.loc[fi, sc] == "A":
                h = flows.loc[fi, "n"] / N; yint[(fi, si)] = (ya, ya + h); ya += h
        yb = a_h + gap
        for fi in order:
            if flows.loc[fi, sc] == "B":
                h = flows.loc[fi, "n"] / N; yint[(fi, si)] = (yb, yb + h); yb += h

    fig, ax = plt.subplots(figsize=(9, 6))
    for si in range(3):
        for st in ("A", "B"):
            y0, y1 = node_y[(si, st)]
            ax.add_patch(patches.Rectangle((si - 0.03, y0), 0.06, y1 - y0,
                         facecolor=COL[st], edgecolor="white", lw=0.5, zorder=3))
            ax.text(si, (y0 + y1) / 2, st, ha="center", va="center",
                    color="white", fontweight="bold", zorder=4, fontsize=10)
    for fi in order:
        col = COL["revert"] if flows.loc[fi, "type"] == "reverted" else COL["stay"]
        al  = 0.80 if flows.loc[fi, "type"] == "reverted" else 0.35
        for si in (0, 1):
            y0a, y0b = yint[(fi, si)]; y1a, y1b = yint[(fi, si + 1)]
            _alluvial_ribbon(ax, si + 0.03, y0a, y0b, si + 1 - 0.03, y1a, y1b, col, al)

    ax.set_xlim(-0.35, 2.35); ax.set_ylim(-0.04, 1 + gap + 0.06)
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(stage_names, fontsize=11)
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.legend(handles=[
        Patch(facecolor=COL["revert"], alpha=0.8, label="reverted (returned to parental compartment)"),
        Patch(facecolor=COL["stay"],   alpha=0.5, label="stayed in migrated compartment")],
        frameon=False, fontsize=8.5, loc="upper center",
        bbox_to_anchor=(0.5, -0.06), ncol=2)
    ax.set_title("Compartment-state flow of migration switchers (label-based)",
                 fontweight="bold", fontsize=11, pad=10)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"Fig_alluvial_states.{ext}"),
                    bbox_inches="tight", dpi=600)
    plt.close(fig)
    print("  Saved -> Fig_alluvial_states.[pdf|png]")


# =============================================================================
# FIGURE 2 - reversal outcome bars (relabeled parental control)
#   LEFT  B10+UNC      : literal reversal test (reverted / unresponsive / reinforced)
#   RIGHT Parental+UNC : SAME directional rule, but the bin never migrated, so the
#                        categories are relabeled by DIRECTION relative to migration
#                        (anti-migration / no move / pro-migration). Shared colours,
#                        per-panel legends. Migration-specific effect = lots of green
#                        on the LEFT, almost none on the RIGHT.
# =============================================================================
def fig_reversal_bars(df, outdir):
    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])]
    slots = ["reverted", "unresponsive", "reinforced", "other"]
    colmap = {"reverted": COL["revert"], "unresponsive": COL["stay"],
              "reinforced": COL["further"], "other": "#BDC3C7"}
    labels_b10 = {"reverted": "reverted", "unresponsive": "unresponsive",
                  "reinforced": "reinforced", "other": "other"}
    labels_par = {"reverted": "anti-migration move", "unresponsive": "no move",
                  "reinforced": "pro-migration move", "other": "other"}

    def ctrl_status(r):
        m, t = r["mig_switch"], r["punc_switch"]
        if t == "Stable":
            return "unresponsive"
        if (m == "A_to_B" and t == "B_to_A") or (m == "B_to_A" and t == "A_to_B"):
            return "reverted"
        if m == t:
            return "reinforced"
        return "other"

    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2), sharey=True)
    panels = [
        (axes[0], "B10 + UNC  (reversal test)", sw["reversal_status"], labels_b10),
        (axes[1], "Parental + UNC  (specificity control)",
         sw.apply(ctrl_status, axis=1), labels_par),
    ]
    for ax, title, status, labels in panels:
        for direction in ["A_to_B", "B_to_A"]:
            mask = sw["mig_switch"] == direction
            sub = status[mask]
            total = len(sub) if len(sub) else 1
            bottom = 0
            xi = 0 if direction == "A_to_B" else 1
            for c in slots:
                frac = (sub == c).sum() / total * 100
                ax.bar(xi, frac, bottom=bottom, color=colmap[c], width=0.6,
                       edgecolor="white", linewidth=0.6)
                if frac > 6:
                    ax.text(xi, bottom + frac / 2, f"{frac:.0f}%", ha="center",
                            va="center", fontsize=8, color="white", fontweight="bold")
                bottom += frac
        ax.set_xticks([0, 1]); ax.set_xticklabels(["A\u2192B", "B\u2192A"])
        ax.set_title(title, fontweight="bold", fontsize=10)
        ax.set_ylim(0, 100)
        ax.legend(handles=[Patch(facecolor=colmap[c], label=labels[c]) for c in slots],
                  frameon=False, fontsize=8, loc="upper center",
                  bbox_to_anchor=(0.5, -0.08), ncol=2)
    axes[0].set_ylabel("% of migration switchers")
    fig.suptitle("Migration-acquired switches under UNC\n"
                 "left: reversal in B10  |  right: same UNC in PARENTAL cells "
                 "(directions are relative to each bin's migration)",
                 fontweight="bold", fontsize=11, y=1.04)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"Fig_reversal_bars.{ext}"),
                    bbox_inches="tight", dpi=600)
    plt.close(fig)
    print("  Saved -> Fig_reversal_bars.[pdf|png]")


# =============================================================================
# FIGURE 3 - migration x treatment transition heatmap (descriptive only)
#   No p-value in the title: the bin-level Fisher test is pseudoreplicated.
#   Inference lives in the block-bootstrap figure/stats.
# =============================================================================
def fig_transition_heatmap(df, outdir):
    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])]
    mig = ["A_to_B", "B_to_A"]
    trt = ["A_to_B", "B_to_A", "Stable"]
    M = np.zeros((len(mig), len(trt)), dtype=int)
    for i, m in enumerate(mig):
        for j, t in enumerate(trt):
            M[i, j] = ((sw["mig_switch"] == m) & (sw["b10unc_switch"] == t)).sum()

    rev = M[0, 1] + M[1, 0]      # reverted: A->B then B->A, or B->A then A->B
    same = M[0, 0] + M[1, 1]     # reinforced: same direction
    try:
        _, fisher_p = fisher_exact([[M[0, 1], M[0, 0]], [M[1, 1], M[1, 0]]])
    except Exception:
        fisher_p = np.nan

    fig, ax = plt.subplots(figsize=(6, 4.6))
    im = ax.imshow(M, cmap="Oranges", aspect="auto")
    for i in range(len(mig)):
        for j in range(len(trt)):
            ax.text(j, i, f"{M[i, j]:,}", ha="center", va="center",
                    fontsize=11, fontweight="bold",
                    color="white" if M[i, j] > M.max() * 0.6 else "black")
    ax.set_xticks(range(len(trt))); ax.set_xticklabels(
        ["B10+UNC\nA\u2192B", "B10+UNC\nB\u2192A", "B10+UNC\nStable"])
    ax.set_yticks(range(len(mig))); ax.set_yticklabels(
        ["migration\nA\u2192B", "migration\nB\u2192A"])
    ax.set_title("Migration switch vs UNC-treatment response\n"
                 f"(reverted={rev:,}, same-direction={same:,})",
                 fontweight="bold", fontsize=10, pad=10)
    plt.colorbar(im, ax=ax, label="bin count", fraction=0.046, pad=0.04)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"Fig_transition_heatmap.{ext}"),
                    bbox_inches="tight", dpi=600)
    plt.close(fig)
    print("  Saved -> Fig_transition_heatmap.[pdf|png]")
    return rev, same, fisher_p


# =============================================================================
# FIGURE 4 - continuous ΔPC1 scatter (RAW; read diagnostic for artifact-removed)
# =============================================================================
def fig_dpc1_scatter(df, outdir):
    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])].copy()
    needed = ["dPC1_migration", "dPC1_b10_treat", "dPC1_par_treat"]
    sw = sw.dropna(subset=needed)
    if len(sw) < 10:
        print("  [WARN] too few bins for ΔPC1 scatter; skipping.")
        return np.nan, np.nan, np.nan, np.nan

    rho_b10, p_b10 = spearmanr(sw["dPC1_migration"], sw["dPC1_b10_treat"])
    rho_par, p_par = spearmanr(sw["dPC1_migration"], sw["dPC1_par_treat"])

    allx = sw["dPC1_migration"].values
    ally = np.concatenate([sw["dPC1_b10_treat"].values, sw["dPC1_par_treat"].values])
    xlim = np.nanpercentile(np.abs(allx), 99); ylim = np.nanpercentile(np.abs(ally), 99)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.4), sharex=True, sharey=True)
    panels = [
        (axes[0], "dPC1_b10_treat", rho_b10, p_b10,
         "Reversal test:  B10 + UNC",
         "ΔPC1 treatment  (Bottom10+UNC \u2212 Bottom10)"),
        (axes[1], "dPC1_par_treat", rho_par, p_par,
         "Specificity control:  Parental + UNC",
         "ΔPC1 treatment  (Parental+UNC \u2212 Parental)"),
    ]
    for ax, ycol, rho, p, title, ylab in panels:
        for d, c in [("A_to_B", COL["AtoB"]), ("B_to_A", COL["BtoA"])]:
            m = sw["mig_switch"] == d
            ax.scatter(sw.loc[m, "dPC1_migration"], sw.loc[m, ycol],
                       s=12, alpha=0.5, color=c, edgecolors="none",
                       label=f"{d.replace('_', ARROW)} migration")
        ax.axhline(0, color="#888", lw=0.7); ax.axvline(0, color="#888", lw=0.7)
        ax.set_xlim(-xlim, xlim); ax.set_ylim(-ylim, ylim)
        ax.set_xlabel("ΔPC1 migration  (Bottom10 \u2212 Parental)")
        ax.set_ylabel(ylab)
        ax.set_title(f"{title}\nSpearman \u03c1={rho:+.2f}  (bin-level p not shown; "
                     "see bootstrap)", fontweight="bold", fontsize=10, pad=8)
    axes[0].legend(frameon=False, fontsize=9, loc="best")
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"Fig_dPC1_scatter.{ext}"),
                    bbox_inches="tight", dpi=600)
    plt.close(fig)
    print(f"  Saved -> Fig_dPC1_scatter.[pdf|png]  (B10 \u03c1={rho_b10:+.2f}, "
          f"parental-control \u03c1={rho_par:+.2f})")
    return rho_b10, p_b10, rho_par, p_par


def _partial_spearman(x, y, z):
    """Spearman partial correlation of x,y controlling for z (rank-Pearson form)."""
    x = np.asarray(x, float); y = np.asarray(y, float); z = np.asarray(z, float)
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[ok], y[ok], z[ok]
    if len(x) < 10:
        return np.nan
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    def _pear(a, b):
        a = a - a.mean(); b = b - b.mean()
        d = np.sqrt((a @ a) * (b @ b))
        return (a @ b) / d if d > 0 else np.nan
    rxy, rxz, ryz = _pear(rx, ry), _pear(rx, rz), _pear(ry, rz)
    denom = np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
    return (rxy - rxz * ryz) / denom if denom > 0 else np.nan


# =============================================================================
# FIGURE 5 - artifact diagnostic (RTM / partial correlation / convergence)
# =============================================================================
def fig_artifact_diagnostic(df, outdir, args):
    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])].copy()
    par_b = pd.to_numeric(sw[args.pc1_parental_dmso], errors="coerce")
    b10_b = pd.to_numeric(sw[args.pc1_b10_dmso], errors="coerce")
    paru  = pd.to_numeric(sw[args.pc1_parental_unc], errors="coerce")
    b10u  = pd.to_numeric(sw[args.pc1_b10_unc], errors="coerce")
    dMig  = sw["dPC1_migration"]; dB10 = sw["dPC1_b10_treat"]; dPar = sw["dPC1_par_treat"]

    def sp(a, b):
        m = np.isfinite(a) & np.isfinite(b)
        return spearmanr(a[m], b[m])[0] if m.sum() > 10 else np.nan

    rho_rtm_b10 = sp(dB10, b10_b)
    rho_rtm_par = sp(dPar, par_b)
    raw_rev  = sp(dMig, dB10)
    raw_ctrl = sp(dMig, dPar)
    part_rev_b10  = _partial_spearman(dMig, dB10, b10_b)
    part_rev_par  = _partial_spearman(dMig, dB10, par_b)
    part_ctrl_par = _partial_spearman(dMig, dPar, par_b)
    d_untreated = np.nanmedian(np.abs(b10_b - par_b))
    d_treated   = np.nanmedian(np.abs(b10u - paru))

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.9))

    for d, c in [("A_to_B", COL["AtoB"]), ("B_to_A", COL["BtoA"])]:
        m = sw["mig_switch"] == d
        axes[0].scatter(b10_b[m], dB10[m], s=10, alpha=0.45, color=c, edgecolors="none")
    axes[0].axhline(0, color="#888", lw=0.7)
    axes[0].set_xlabel("Bottom10 baseline PC1 (DMSO)")
    axes[0].set_ylabel("\u0394PC1 B10 treatment (B10+UNC \u2212 B10)")
    axes[0].set_title(f"RTM / attractor (B10)\n\u03c1={rho_rtm_b10:+.2f}  "
                      "(neg = pull toward common level)", fontsize=9.5, fontweight="bold")

    for d, c in [("A_to_B", COL["AtoB"]), ("B_to_A", COL["BtoA"])]:
        m = sw["mig_switch"] == d
        axes[1].scatter(par_b[m], dPar[m], s=10, alpha=0.45, color=c, edgecolors="none")
    axes[1].axhline(0, color="#888", lw=0.7)
    axes[1].set_xlabel("Parental baseline PC1 (DMSO)")
    axes[1].set_ylabel("\u0394PC1 parental treatment (Par+UNC \u2212 Par)")
    axes[1].set_title(f"RTM / attractor (Parental)\n\u03c1={rho_rtm_par:+.2f}  "
                      "(neg = pull toward common level)", fontsize=9.5, fontweight="bold")

    labels = ["Reversal\nraw", "Reversal\n| B10 base", "Reversal\n| Par base",
              "Control\nraw", "Control\n| Par base"]
    vals = [raw_rev, part_rev_b10, part_rev_par, raw_ctrl, part_ctrl_par]
    cols = [COL["BtoA"] if v < 0 else COL["AtoB"] for v in vals]
    axes[2].bar(range(len(vals)), vals, color=cols, edgecolor="white", width=0.7)
    axes[2].axhline(0, color="#333", lw=1.0)
    for i, v in enumerate(vals):
        axes[2].text(i, v + (0.03 if v >= 0 else -0.03), f"{v:+.2f}",
                     ha="center", va="bottom" if v >= 0 else "top", fontsize=8.5)
    axes[2].set_xticks(range(len(labels)))
    axes[2].set_xticklabels(labels, fontsize=8)
    axes[2].set_ylabel("Spearman \u03c1")
    axes[2].set_ylim(-1, 1)
    axes[2].set_title("Raw vs PARTIAL correlation", fontsize=9.5, fontweight="bold")

    conv_txt = ("CONVERGES to common state" if d_treated < d_untreated
                else "does NOT converge")
    fig.suptitle("Artifact diagnostic for the reversal/control correlations\n"
                 f"Treated-state distance {d_treated:.4f} vs untreated "
                 f"{d_untreated:.4f}  \u2192 {conv_txt}",
                 fontweight="bold", fontsize=11, y=1.06)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"Fig_artifact_diagnostic.{ext}"),
                    bbox_inches="tight", dpi=600)
    plt.close(fig)
    print("  Saved -> Fig_artifact_diagnostic.[pdf|png]")
    print(f"  [DIAG] RTM \u03c1: B10={rho_rtm_b10:+.2f}, parental={rho_rtm_par:+.2f}")
    print(f"  [DIAG] reversal \u03c1 raw={raw_rev:+.2f} -> partial|B10base={part_rev_b10:+.2f}, "
          f"partial|Parbase={part_rev_par:+.2f}")
    print(f"  [DIAG] control  \u03c1 raw={raw_ctrl:+.2f} -> partial|Parbase={part_ctrl_par:+.2f}")
    print(f"  [DIAG] convergence: treated {d_treated:.4f} vs untreated "
          f"{d_untreated:.4f} -> {conv_txt}")
    return dict(rho_rtm_b10=rho_rtm_b10, rho_rtm_par=rho_rtm_par,
                raw_rev=raw_rev, part_rev_b10=part_rev_b10, part_rev_par=part_rev_par,
                raw_ctrl=raw_ctrl, part_ctrl_par=part_ctrl_par,
                dist_treated=d_treated, dist_untreated=d_untreated)


# =============================================================================
# Attractor-removed migration-specific reversal subset
# =============================================================================
def add_specific_reversal_flag(df, args):
    """Regress (ΔPC1 B10-treat) on the B10 DMSO baseline; the residual is the
    treatment move NOT explained by starting level (attractor removed). Flag bins
    whose residual OPPOSES migration AND whose parental copy is unresponsive."""
    sw_mask = df["mig_switch"].isin(["A_to_B", "B_to_A"])
    dmig = df["dPC1_migration"]; db10 = df["dPC1_b10_treat"]; dpar = df["dPC1_par_treat"]
    b10_base = pd.to_numeric(df[args.pc1_b10_dmso], errors="coerce")

    m = sw_mask & np.isfinite(db10) & np.isfinite(b10_base)
    if m.sum() > 10:
        slope, intercept = np.polyfit(b10_base[m].values, db10[m].values, 1)
        pred = slope * b10_base + intercept
    else:
        pred = pd.Series(0.0, index=df.index)
    resid_b10 = db10 - pred
    df["resid_b10_treat"] = resid_b10

    if args.unresp_thresh is not None:
        tau = float(args.unresp_thresh)
    else:
        tau = float(np.nanpercentile(np.abs(resid_b10[sw_mask].dropna()), 66))
    tau_par = float(np.nanpercentile(np.abs(dpar[sw_mask].dropna()), 50))

    resid_reverts   = (np.sign(resid_b10) == -np.sign(dmig)) & (np.abs(resid_b10) > tau)
    parental_unresp = np.abs(dpar) < tau_par
    df["parental_unresponsive"] = (sw_mask & parental_unresp).fillna(False)
    df["specific_reversal"] = (sw_mask & resid_reverts & parental_unresp).fillna(False)
    return tau


def fig_specific_reversal(df, outdir, tau, args):
    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])].copy()
    spec = sw[sw["specific_reversal"]]
    n_spec = len(spec); n_sw = len(sw)
    med_par_move = float(np.nanmedian(np.abs(spec["dPC1_par_treat"]))) if n_spec else np.nan
    med_b10_move = float(np.nanmedian(np.abs(spec["dPC1_b10_treat"]))) if n_spec else np.nan
    rr = sw.dropna(subset=["dPC1_migration", "resid_b10_treat"])
    rho_resid = spearmanr(rr["dPC1_migration"], rr["resid_b10_treat"])[0] if len(rr) > 10 else np.nan

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0))
    other = sw[~sw["specific_reversal"]]
    axes[0].scatter(other["dPC1_migration"], other["resid_b10_treat"],
                    s=9, alpha=0.25, color="#BBBBBB", edgecolors="none",
                    label="other switchers")
    for d, c in [("A_to_B", COL["AtoB"]), ("B_to_A", COL["BtoA"])]:
        mm = spec["mig_switch"] == d
        axes[0].scatter(spec.loc[mm, "dPC1_migration"], spec.loc[mm, "resid_b10_treat"],
                        s=16, alpha=0.7, color=c, edgecolors="none",
                        label=f"specific reversal ({d.replace('_', ARROW)})")
    axes[0].axhline(0, color="#888", lw=0.7); axes[0].axvline(0, color="#888", lw=0.7)
    axes[0].set_xlabel("ΔPC1 migration  (Bottom10 \u2212 Parental)")
    axes[0].set_ylabel("RESIDUAL ΔPC1 B10-treat\n(attractor regressed out)")
    axes[0].set_title(f"Attractor-removed reversal  (residual \u03c1={rho_resid:+.2f})\n"
                      f"specific subset n={n_spec}/{n_sw} ({n_spec/max(n_sw,1):.1%}); "
                      f"parental |ΔPC1| med={med_par_move:.4f}",
                      fontsize=9.5, fontweight="bold")
    axes[0].legend(frameon=False, fontsize=7.5, loc="best")

    base = tau if tau > 0 else 0.005
    taus = np.linspace(0.4 * base, 2.0 * base, 12)
    dmig = sw["dPC1_migration"]; resid = sw["resid_b10_treat"]; dpar = sw["dPC1_par_treat"]
    tau_par = float(np.nanpercentile(np.abs(dpar.dropna()), 50))
    pu = np.abs(dpar) < tau_par
    sizes = []
    for t in taus:
        rr2 = (np.sign(resid) == -np.sign(dmig)) & (np.abs(resid) > t)
        sizes.append(int((rr2 & pu).sum()))
    axes[1].plot(taus, sizes, "-o", color=COL["revert"], lw=1.6, ms=4)
    axes[1].axvline(tau, color="#333", ls="--", lw=1.0, label=f"chosen tau={tau:.4f}")
    axes[1].set_xlabel("residual threshold tau  (|residual ΔPC1|)")
    axes[1].set_ylabel("No. of migration-specific reversal bins")
    axes[1].set_title("Sensitivity of subset size to tau\n(smooth trend = not cherry-picked)",
                      fontsize=9.5, fontweight="bold")
    axes[1].legend(frameon=False, fontsize=8.5)

    fig.suptitle("Migration-specific reversal (attractor regressed out)\n"
                 "residual-reverting \u2227 parental-unresponsive bins",
                 fontweight="bold", fontsize=11, y=1.04)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"Fig_specific_reversal.{ext}"),
                    bbox_inches="tight", dpi=600)
    plt.close(fig)
    print(f"  Saved -> Fig_specific_reversal.[pdf|png]")
    print(f"  [SPECIFIC] migration-specific reversal: {n_spec}/{n_sw} "
          f"({n_spec/max(n_sw,1):.1%}); tau={tau:.4f}; residual \u03c1={rho_resid:+.2f}; "
          f"parental |ΔPC1| median={med_par_move:.4f}")
    return dict(n_specific=n_spec, n_switchers=n_sw, tau=tau, rho_resid=rho_resid,
                med_par_move=med_par_move, med_b10_move=med_b10_move)


# =============================================================================
# CHROMOSOME BLOCK-BOOTSTRAP  (the best inference for this script)
# -----------------------------------------------------------------------------
# 250 kb bins are autocorrelated, so per-bin Fisher/Spearman p-values are
# pseudoreplicated. We resample WHOLE CHROMOSOMES with replacement (preserving
# each chromosome's internal correlation), recompute every headline statistic on
# each resample, and report percentile 95% CIs. The number of independent units
# is the number of chromosomes (~23) - the best N.
#
# Bootstrapped statistics (over migration switchers):
#   rev_rate_b10 : P(reverted | switcher)            in B10+UNC
#   ctrl_rate    : P(anti-migration | switcher)      in Parental+UNC (control)
#   excess       : rev_rate_b10 - ctrl_rate
#   odds_ratio   : Haldane-corrected OR (reverting in B10 vs control)
#   residual_rho : Spearman(dPC1_migration, residual ΔPC1 B10-treat); the attractor
#                  regression is REFIT inside each resample so the CI includes the
#                  attractor-fit uncertainty.
# =============================================================================
def _haldane_or(a, b, c, d):
    """Odds ratio with Haldane-Anscombe +0.5 correction (handles empty cells)."""
    a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    return (a * d) / (b * c)


def _spearman_safe(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 10:
        return np.nan
    xx, yy = x[m], y[m]
    if np.all(xx == xx[0]) or np.all(yy == yy[0]):
        return np.nan
    return spearmanr(xx, yy).correlation


def _resid_rho(dmig, db10, b10base):
    """Refit attractor (db10 ~ b10base) on THIS resample, then residual rho."""
    m = np.isfinite(db10) & np.isfinite(b10base) & np.isfinite(dmig)
    if m.sum() < 10:
        return np.nan
    x, y, mg = b10base[m], db10[m], dmig[m]
    if np.all(x == x[0]):
        resid = y - y.mean()
    else:
        slope, intercept = np.polyfit(x, y, 1)
        resid = y - (slope * x + intercept)
    return _spearman_safe(mg, resid)


def block_bootstrap(df, args):
    """Observed point estimates + chromosome block-bootstrap 95% CIs, plus the raw
    bootstrap arrays (for the CI figure)."""
    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])].copy()
    n_sw = len(sw)
    if n_sw < 10:
        print("  [WARN] too few switchers for block-bootstrap; skipping.")
        return None

    chrom   = sw["_chrom"].to_numpy()
    is_revB = (sw["reversal_status"] == "reverted").to_numpy()
    is_revC = ctrl_reverted_mask(sw).to_numpy()
    dmig    = pd.to_numeric(sw["dPC1_migration"], errors="coerce").to_numpy(float)
    db10    = pd.to_numeric(sw["dPC1_b10_treat"], errors="coerce").to_numpy(float)
    b10base = pd.to_numeric(sw[args.pc1_b10_dmso], errors="coerce").to_numpy(float)

    uniq = pd.unique(chrom)
    blocks = {ch: np.where(chrom == ch)[0] for ch in uniq}     # resampling blocks
    n_blocks = len(uniq)

    def stats_on(idx):
        rb = is_revB[idx]; rc = is_revC[idx]
        a = int(rb.sum()); b = int((~rb).sum())
        c = int(rc.sum()); d = int((~rc).sum())
        return (rb.mean(), rc.mean(), rb.mean() - rc.mean(),
                _haldane_or(a, b, c, d),
                _resid_rho(dmig[idx], db10[idx], b10base[idx]))

    obs = stats_on(np.arange(n_sw))

    rng = np.random.default_rng(args.boot_seed)
    boot = np.full((args.n_boot, 5), np.nan)
    for it in range(args.n_boot):
        pick = rng.choice(uniq, size=n_blocks, replace=True)   # resample chromosomes
        idx = np.concatenate([blocks[ch] for ch in pick])
        boot[it] = stats_on(idx)

    names = ["rev_rate_b10", "ctrl_rate", "excess", "odds_ratio", "residual_rho"]
    ci = {}
    for k, nm in enumerate(names):
        lo, hi = np.nanpercentile(boot[:, k], [2.5, 97.5])
        ci[nm] = (obs[k], lo, hi)
    ci["_n_blocks"] = n_blocks
    ci["_n_boot"] = args.n_boot
    return ci


def fig_block_bootstrap(ci, outdir):
    """Forest plot of headline effects with chromosome block-bootstrap 95% CIs."""
    if ci is None:
        return
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2),
                             gridspec_kw={"width_ratios": [2, 1.1]})

    propnames = [("rev_rate_b10", "B10 reversal rate"),
                 ("ctrl_rate", "parental control rate"),
                 ("excess", "excess (B10 \u2212 parental)")]
    ypos = list(range(len(propnames)))[::-1]
    for y, (key, lab) in zip(ypos, propnames):
        obs, lo, hi = ci[key]
        axes[0].errorbar(obs * 100, y,
                         xerr=[[(obs - lo) * 100], [(hi - obs) * 100]],
                         fmt="o", color=COL["revert"], ecolor="#333",
                         capsize=4, ms=7, lw=1.4)
        axes[0].text(obs * 100, y + 0.18,
                     f"{obs*100:.1f}%  [{lo*100:.1f}, {hi*100:.1f}]",
                     ha="center", va="bottom", fontsize=8.5)
    axes[0].axvline(0, color="#888", lw=0.8, ls="--")
    axes[0].set_yticks(ypos); axes[0].set_yticklabels([l for _, l in propnames])
    axes[0].set_xlabel("percent of migration switchers (95% CI)")
    axes[0].set_title("Reversal vs specificity control", fontweight="bold", fontsize=10)
    axes[0].set_ylim(-0.6, len(propnames) - 0.2)

    obs, lo, hi = ci["residual_rho"]
    axes[1].errorbar(obs, 0, xerr=[[obs - lo], [hi - obs]], fmt="o",
                     color=COL["BtoA"], ecolor="#333", capsize=4, ms=7, lw=1.4)
    axes[1].text(obs, 0.12, f"{obs:+.2f}  [{lo:+.2f}, {hi:+.2f}]",
                 ha="center", va="bottom", fontsize=8.5)
    axes[1].axvline(0, color="#888", lw=0.8, ls="--")
    axes[1].set_yticks([0]); axes[1].set_yticklabels(["residual \u03c1\n(attractor removed)"])
    axes[1].set_xlim(min(-0.05, lo - 0.05), max(0.05, hi + 0.05))
    axes[1].set_ylim(-0.6, 0.8)
    axes[1].set_xlabel("Spearman \u03c1 (95% CI)")
    o_obs, o_lo, o_hi = ci["odds_ratio"]
    axes[1].set_title(f"Attractor-removed effect\nOR={o_obs:.1f} [{o_lo:.1f}, {o_hi:.1f}]",
                      fontweight="bold", fontsize=10)

    fig.suptitle("Chromosome block-bootstrap 95% CIs  "
                 f"({ci['_n_blocks']} chromosomes, {ci['_n_boot']} resamples)\n"
                 "Blocks = chromosomes, NOT individual bins",
                 fontweight="bold", fontsize=11, y=1.06)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"Fig_block_bootstrap_CI.{ext}"),
                    bbox_inches="tight", dpi=600)
    plt.close(fig)
    print("  Saved -> Fig_block_bootstrap_CI.[pdf|png]")


# =============================================================================
# Gene-list export (enrichment-ready)
# =============================================================================
def export_gene_lists(df, args, outdir):
    gdir = os.path.join(outdir, "genelists")
    os.makedirs(gdir, exist_ok=True)
    cc, gn = args.coord_col, args.gene_name_col

    genes = pd.read_csv(args.gene_list, sep="\t", low_memory=False)
    for c in (cc, gn):
        if c not in genes.columns:
            sys.exit(f"[ERROR] '{c}' not in gene_list. Columns: {list(genes.columns)}")
    genes[cc] = genes[cc].astype(str).str.strip()
    genes[gn] = genes[gn].astype(str).str.strip()
    df[cc] = df[cc].astype(str).str.strip()

    keep = [cc, "mig_switch", "reversal_status", "baseline_class"]
    if "categorical_specific_reversal" in df.columns:
        keep.append("categorical_specific_reversal")
    if "specific_reversal" in df.columns:
        keep.append("specific_reversal")
    g = genes.merge(df[keep], on=cc, how="inner")

    def write_set(symbols, name):
        syms = sorted(set(s for s in symbols if s and s.lower() != "nan"))
        path = os.path.join(gdir, f"{name}.txt")
        with open(path, "w") as fh:
            fh.write("\n".join(syms) + "\n")
        print(f"  gene set [{name}]: {len(syms)} genes -> {path}")

    write_set(g[gn], "background_universe")
    write_set(g.loc[g["mig_switch"] == "A_to_B", gn], "switcher_A_to_B")
    write_set(g.loc[g["mig_switch"] == "B_to_A", gn], "switcher_B_to_A")
    write_set(g.loc[g["reversal_status"] == "reverted", gn], "reverted_under_UNC")
    
    # Categorical-specific reversal (reverted in B10, stable in parental)
    cat_spec_genes = g.loc[g["categorical_specific_reversal"] == True, gn]
    write_set(cat_spec_genes, "categorical_specific_reversal")
    n_cat_spec_bins = (df["categorical_specific_reversal"]).sum()
    n_cat_spec_genes = len(set(s for s in cat_spec_genes if s and s.lower() != "nan"))
    total_sw = len(df[df["mig_switch"].isin(["A_to_B", "B_to_A"])])
    pct_cat_spec = (n_cat_spec_bins / total_sw * 100) if total_sw > 0 else 0
    print(f"  [INFO] Categorical-specific: {n_cat_spec_bins} bins ({pct_cat_spec:.1f}% of switchers), {n_cat_spec_genes} genes")
    
    if "specific_reversal" in g.columns:
        write_set(g.loc[g["specific_reversal"] == True, gn],
                  "migration_specific_reversal")
    if args.include_unresponsive:
        write_set(g.loc[g["reversal_status"] == "unresponsive", gn],
                  "switcher_unresponsive_to_UNC")
    if args.include_stable:
        write_set(g.loc[g["baseline_class"] == "Stable_A", gn], "stable_A")
        write_set(g.loc[g["baseline_class"] == "Stable_B", gn], "stable_B")

    audit_cols = [gn, cc, "mig_switch", "reversal_status", "baseline_class"]
    if "categorical_specific_reversal" in g.columns:
        audit_cols.append("categorical_specific_reversal")
    if "specific_reversal" in g.columns:
        audit_cols.append("specific_reversal")
    g[audit_cols].to_csv(
        os.path.join(gdir, "gene_assignments.txt"), sep="\t", index=False)
    print(f"  Per-gene assignments -> {os.path.join(gdir, 'gene_assignments.txt')}")


# =============================================================================
# Main
# =============================================================================
def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print("[INFO] Loading PC1 / compartment table...")
    df = pd.read_csv(args.pc1, sep="\t", low_memory=False)
    df[args.coord_col] = df[args.coord_col].astype(str).str.strip()
    print(f"  {len(df):,} bins x {len(df.columns)} columns")

    need = ([args.pc1_parental_dmso, args.pc1_b10_dmso, args.pc1_parental_unc,
             args.pc1_b10_unc] + args.mig_cols + args.punc_cols + args.b10unc_cols)
    missing = [c for c in need if c not in df.columns]
    if missing:
        sys.exit("[ERROR] Missing required columns:\n  " + "\n  ".join(missing)
                 + "\n\nAvailable:\n  " + "\n  ".join(map(str, df.columns)))

    df = build_reversal_table(df, args)
    tau_spec = add_specific_reversal_flag(df, args)

    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])]
    n_sw = len(sw)
    rev_rate = (sw["reversal_status"] == "reverted").mean() if n_sw else np.nan
    unresp_rate = (sw["reversal_status"] == "unresponsive").mean() if n_sw else np.nan
    cat_spec_rate = (sw["categorical_specific_reversal"]).mean() if n_sw else np.nan
    print(f"[INFO] Migration switchers: {n_sw:,}  (across {df['_chrom'].nunique()} chromosomes)")
    print(f"[INFO] Reverted under UNC-in-B10: {rev_rate:.1%}")
    print(f"[INFO] Unresponsive to UNC-in-B10: {unresp_rate:.1%}")
    print(f"[INFO] Categorical-specific reversal: {(sw['categorical_specific_reversal']).sum():,} bins ({cat_spec_rate:.1%})")

    df.to_csv(os.path.join(args.outdir, "Reversal_classified.txt"),
              sep="\t", index=False)

    # ---- figures ----------------------------------------------------------
    print("[INFO] Figures...")
    fig_alluvial(df, args.outdir, args)
    fig_reversal_bars(df, args.outdir)
    rev, same, fisher_p = fig_transition_heatmap(df, args.outdir)
    rho_b10, sp_p_b10, rho_par, sp_p_par = fig_dpc1_scatter(df, args.outdir)
    diag = fig_artifact_diagnostic(df, args.outdir, args)
    spec_summary = fig_specific_reversal(df, args.outdir, tau_spec, args)

    # ---- categorical specificity (point estimates) ------------------------
    sw_rev_b10 = (sw["reversal_status"] == "reverted")
    sw_rev_par = ctrl_reverted_mask(sw)
    rev_rate_b10 = sw_rev_b10.mean() if n_sw else np.nan
    ctrl_rate = sw_rev_par.mean() if n_sw else np.nan
    a = int(sw_rev_b10.sum()); b = int((~sw_rev_b10).sum())
    c = int(sw_rev_par.sum()); d_ = int((~sw_rev_par).sum())
    try:
        excess_OR, excess_p = fisher_exact([[a, b], [c, d_]])
    except Exception:
        excess_OR, excess_p = np.nan, np.nan
    excess_pts = rev_rate_b10 - ctrl_rate

    # ---- BEST inference: chromosome block-bootstrap ---------------------
    print("[INFO] Chromosome block-bootstrap (best CIs; blocks = chromosomes)...")
    ci = block_bootstrap(df, args)
    fig_block_bootstrap(ci, args.outdir)
    if ci is not None:
        for nm, lab in [("rev_rate_b10", "B10 reversal rate"),
                        ("ctrl_rate", "parental control rate"),
                        ("excess", "excess reversal"),
                        ("odds_ratio", "odds ratio"),
                        ("residual_rho", "residual rho (attractor-removed)")]:
            o, lo, hi = ci[nm]
            if nm in ("rev_rate_b10", "ctrl_rate", "excess"):
                print(f"  [BOOT] {lab:34s} {o*100:6.1f}%  95% CI [{lo*100:.1f}, {hi*100:.1f}]")
            elif nm == "odds_ratio":
                print(f"  [BOOT] {lab:34s} {o:6.2f}   95% CI [{lo:.2f}, {hi:.2f}]")
            else:
                print(f"  [BOOT] {lab:34s} {o:+6.2f}   95% CI [{lo:+.2f}, {hi:+.2f}]")

    # ---- stats file -------------------------------------------------------
    with open(os.path.join(args.outdir, "Reversal_stats.txt"), "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"migration_switchers\t{n_sw}\n")
        fh.write(f"n_chromosomes\t{df['_chrom'].nunique()}\n")

        fh.write("# === HEADLINE: effect sizes + chromosome block-bootstrap 95% CI ===\n")
        if ci is not None:
            for nm in ("rev_rate_b10", "ctrl_rate", "excess",
                       "odds_ratio", "residual_rho"):
                o, lo, hi = ci[nm]
                fh.write(f"{nm}\t{o:.5f}\n")
                fh.write(f"{nm}_CI95_lo\t{lo:.5f}\n")
                fh.write(f"{nm}_CI95_hi\t{hi:.5f}\n")
            fh.write(f"bootstrap_n_chromosome_blocks\t{ci['_n_blocks']}\n")
            fh.write(f"bootstrap_n_resamples\t{ci['_n_boot']}\n")

        fh.write("# --- descriptive transition counts ---\n")
        fh.write(f"transition_reverted_bins\t{rev}\n")
        fh.write(f"transition_same_direction_bins\t{same}\n")
        fh.write(f"unresponsive_rate_B10UNC\t{unresp_rate:.4f}\n")

        fh.write("# --- artifact diagnostic (shared-variable / regression-to-mean) ---\n")
        fh.write(f"RTM_rho_B10treat_vs_B10baseline\t{diag['rho_rtm_b10']:.4f}\n")
        fh.write(f"RTM_rho_partreat_vs_parbaseline\t{diag['rho_rtm_par']:.4f}\n")
        fh.write(f"reversal_rho_raw\t{diag['raw_rev']:.4f}\n")
        fh.write(f"reversal_rho_partial_ctrl_B10baseline\t{diag['part_rev_b10']:.4f}\n")
        fh.write(f"reversal_rho_partial_ctrl_parbaseline\t{diag['part_rev_par']:.4f}\n")
        fh.write(f"control_rho_raw\t{diag['raw_ctrl']:.4f}\n")
        fh.write(f"control_rho_partial_ctrl_parbaseline\t{diag['part_ctrl_par']:.4f}\n")
        fh.write(f"convergence_dist_treated\t{diag['dist_treated']:.5f}\n")
        fh.write(f"convergence_dist_untreated\t{diag['dist_untreated']:.5f}\n")

        fh.write("# --- migration-specific reversal subset (attractor-removed) ---\n")
        fh.write(f"specific_reversal_tau\t{spec_summary['tau']:.5f}\n")
        fh.write(f"specific_reversal_n_bins\t{spec_summary['n_specific']}\n")
        fh.write(f"specific_reversal_residual_spearman_rho\t{spec_summary['rho_resid']:.4f}\n")
        fh.write(f"specific_reversal_frac_of_switchers\t"
                 f"{spec_summary['n_specific']/max(spec_summary['n_switchers'],1):.4f}\n")
        fh.write(f"specific_reversal_parental_move_median\t{spec_summary['med_par_move']:.5f}\n")
        fh.write(f"specific_reversal_B10_move_median\t{spec_summary['med_b10_move']:.5f}\n")

        fh.write("# === PSEUDOREPLICATED bin-level p-values - NOT FOR INFERENCE ===\n")
        fh.write("# (250kb bins are autocorrelated; these treat bins as independent\n")
        fh.write("#  and are anti-conservative. Use the block-bootstrap CIs above.)\n")
        fh.write(f"PSEUDOREP_excess_reversal_fisher_OR\t{excess_OR:.4f}\n")
        fh.write(f"PSEUDOREP_excess_reversal_fisher_p\t{excess_p:.3e}\n")
        fh.write(f"PSEUDOREP_transition_fisher_p\t{fisher_p:.3e}\n")
        fh.write(f"PSEUDOREP_dPC1_spearman_rho_B10treat\t{rho_b10:.4f}\n")
        fh.write(f"PSEUDOREP_dPC1_spearman_p_B10treat\t{sp_p_b10:.3e}\n")
        fh.write(f"PSEUDOREP_dPC1_spearman_rho_parentalControl\t{rho_par:.4f}\n")
        fh.write(f"PSEUDOREP_dPC1_spearman_p_parentalControl\t{sp_p_par:.3e}\n")
    print(f"[INFO] Stats -> {os.path.join(args.outdir, 'Reversal_stats.txt')}")

    if args.gene_list:
        print("[INFO] Exporting enrichment-ready gene lists...")
        export_gene_lists(df, args, args.outdir)
    else:
        print("[INFO] --gene_list not given; skipping gene-list export.")

    print(f"\n[INFO] Done -> {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()