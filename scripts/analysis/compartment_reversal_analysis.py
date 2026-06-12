"""
compartment_reversal_analysis.py
================================================================================
Does EZH2 inhibition (UNC1999, "UNC") REVERSE the compartment switches that A375
melanoma cells ACQUIRE during constricted migration (Parental -> Bottom10 / B10)?

DESIGN  (every delta is measured inside the single drug experiment, DMSO baseline)
    migration         : Bottom10-DMSO - Parental-DMSO   (what migration changed)
    B10 treatment     : Bottom10-UNC  - Bottom10-DMSO    (does UNC undo it in B10?)
    parental treatment: Parental-UNC  - Parental-DMSO    (does UNC act in non-migrated cells?)
Each comparison is collapsed to one categorical token in {A_to_B, B_to_A, Stable}
from its Both/Switch/Shift column triple (see parse_switch_token).

HEADLINE  --  migration-specific reversal
    A migration switcher bin is "migration-specifically reverted" iff
        (1) UNC moved it back toward parental in Bottom10   -> reverted
        (2) UNC left the parental copy unchanged            -> parental Stable
    i.e. REVERTED in B10 AND parental STABLE under UNC.
    Equivalently: take all reverted bins, then REMOVE the ones whose parental copy
    ALSO switched under UNC (at those loci UNC acts regardless of migration). So
        migration-specific reversal = (reverted bins) - (reverted bins that ALSO
                                       switched in parental)
    This is NOT "(reverted) - (all parental-switched)": most parental-switched bins
    never reverted, so only the OVERLAP is removed. The Venn figure shows this
    directly. Column `migration_specific_reversal`, gene list migration_specific_reversal.txt.

INFERENCE
    250 kb bins are autocorrelated (compartments are megabase-scale), so per-bin
    Fisher/Spearman p-values are pseudoreplicated; they are written only as
    PSEUDOREP_* and are NOT for inference. Headline uncertainty = CHROMOSOME
    block-bootstrap 95% CIs (resample whole chromosomes) for: the B10 reversal rate,
    the parental-switch rate, the migration-specific reversal rate, the specificity
    odds ratio, and the attractor-removed residual rho.

    Specificity odds ratio: its control is the parental ANTI-MIGRATION switch (the
    same reverting move, but in parental), because the OR asks "is the REVERSAL
    event itself rarer in parental?" - the matched control for a reverting move is
    that same reverting move. The SET above instead removes ANY parental switch,
    because the SET asks "at which loci is UNC migration-gated?". Two questions,
    two correct controls.

ROBUSTNESS  (continuous, optional, NOT the headline)
    A fully continuous attractor-removed reversal set is also produced: regress the
    B10 treatment delta on the B10 baseline, keep bins whose residual reverts AND
    whose parental copy is unresponsive. Only a regression-to-the-mean control.
    Column `attractor_removed_reversal`, gene list attractor_removed_reversal_robustness.txt.

USAGE
    python compartment_reversal_analysis.py \
        --pc1 A375-CompartmentsAcrossConditions.txt \
        --gene_list GeneListTSShg19_cleaned.txt \
        --outdir results/exploratory/reversal
    The Venn figure uses matplotlib-venn if installed (recommended, area-proportional):
        pip install matplotlib-venn
    otherwise a built-in two-circle fallback is drawn. Column-name overrides are
    available for every PC1/token column (see parse_args).
================================================================================
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")                 # headless backend (cluster / no display)
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import Circle
import matplotlib.patches as patches
from matplotlib.patches import Patch
from scipy.stats import fisher_exact, spearmanr, rankdata


# ---- shared publication palette ------------------------------------------------
# Refined, muted, print-friendly tones (ColorBrewer-Dark2-inspired).
COL = {
    "A":      "#B23A48",   # compartment A  (warm red)
    "B":      "#1F6F8B",   # compartment B  (teal-blue)
    "revert": "#2A9D8F",   # reverted        (emerald-teal)
    "stay":   "#9BA7B0",   # unresponsive / Stable (cool grey)
    "further":"#6D5BA6",   # reinforced      (muted indigo)
    "switch": "#E08A3C",   # parental switched (warm amber)
    "overlap":"#6A4C93",   # Venn overlap     (distinct purple, not a teal/amber blend)
    "AtoB":   "#1F6F8B",   # scatter: A->B migration
    "BtoA":   "#B23A48",   # scatter: B->A migration
}
ARROW = "\u2192"


# ================================================================================
# Token parsing
#   Each directional change is stored as a Both / Switch / Shift column triple.
#   parse_switch_token collapses that triple into a single {A_to_B, B_to_A, Stable}
#   call. Sub-threshold "Toward"/"Shift" calls fold into their direction, so a token
#   is a DIRECTION-of-change call, not a strict zero-crossing.
# ================================================================================
def parse_switch_token(both, switch, shift):
    """Return 'A_to_B', 'B_to_A' or 'Stable' from one Both/Switch/Shift triple."""
    for v in (both, switch, shift):              # first informative field wins
        s = str(v).strip()
        if s in ("-", "", "nan", "None", "NaN"):
            continue
        if "BtoA" in s or "BothtoA" in s or "TowardA" in s or s == "toA":
            return "B_to_A"
        if "AtoB" in s or "BothtoB" in s or "TowardB" in s or s == "toB":
            return "A_to_B"
    return "Stable"


def classify_axis(df, both_c, switch_c, shift_c):
    """Vectorised parse_switch_token over one Both/Switch/Shift column triple."""
    return df.apply(lambda r: parse_switch_token(
        r.get(both_c, "-"), r.get(switch_c, "-"), r.get(shift_c, "-")), axis=1)


def chrom_of(coord):
    """Extract a chromosome label from a coordinate string, for block-bootstrap
    blocking. Handles 'chr1:..', 'chr1-..', 'chr1_..' and bare 'chr1'."""
    s = str(coord)
    m = re.search(r"(chr[0-9XYMxym]+)", s)
    if m:
        return m.group(1)
    for sep in (":", "-", "_", "\t", " "):
        if sep in s:
            return s.split(sep)[0]
    return s


# ================================================================================
# Command-line interface
# ================================================================================
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pc1", required=True,
                   help="Compartment/PC1 table (TSV) with one row per 250kb bin.")
    p.add_argument("--outdir", default="results/exploratory/reversal")
    p.add_argument("--gene_list", default=None,
                   help="TSV with --coord_col and --gene_name_col; enables gene-list export.")
    p.add_argument("--coord_col", default="coordinate")
    p.add_argument("--gene_name_col", default="gene_name")
    p.add_argument("--chrom_col", default=None,
                   help="Explicit chromosome column for block-bootstrap blocking; "
                        "if absent, chromosome is parsed from --coord_col.")
    # PC1 value columns (DMSO baseline design)
    p.add_argument("--pc1_parental_dmso", default="Parental_DMSO")
    p.add_argument("--pc1_b10_dmso",      default="Bottom10-DMSO")
    p.add_argument("--pc1_parental_unc",  default="Parental_UNC")
    p.add_argument("--pc1_b10_unc",       default="Bottom10-UNC")
    # token column triples for the three directional comparisons
    p.add_argument("--mig_cols", nargs=3,
                   default=["P->B10_Both_DMSO", "P->B10_Switch_DMSO", "P->B10_Shift_DMSO"],
                   metavar=("BOTH", "SWITCH", "SHIFT"))
    p.add_argument("--punc_cols", nargs=3,
                   default=["P_DMSO->UNC_Both", "P_DMSO->UNC_switch", "P_DMSO->UNC_shift"],
                   metavar=("BOTH", "SWITCH", "SHIFT"))
    p.add_argument("--b10unc_cols", nargs=3,
                   default=["B10_DMSO->UNC_Both", "B10_DMSO->UNC_switch", "B10_DMSO->UNC_shift"],
                   metavar=("BOTH", "SWITCH", "SHIFT"))
    # optional extra gene lists
    p.add_argument("--include_stable", action="store_true",
                   help="Also export stable_A / stable_B gene lists.")
    p.add_argument("--include_unresponsive", action="store_true",
                   help="Also export the switcher-unresponsive gene list.")
    p.add_argument("--include_robustness", action="store_true",
                   help="Also export the continuous attractor-removed robustness gene list.")
    # continuous robustness threshold
    p.add_argument("--unresp_thresh", type=float, default=None,
                   help="|residual dPC1| threshold for the continuous robustness subset "
                        "(default: 66th pct of |residual|).")
    # block-bootstrap controls
    p.add_argument("--n_boot", type=int, default=2000)
    p.add_argument("--boot_seed", type=int, default=0)
    return p.parse_args()


# ================================================================================
# Build the per-bin reversal table
#   Three INDEPENDENT directional comparisons, each vs its OWN baseline:
#     mig_switch    : Parental-DMSO -> Bottom10-DMSO   (migration axis)
#     b10unc_switch : Bottom10-DMSO -> Bottom10-UNC    (treatment, from the B10 state)
#     punc_switch   : Parental-DMSO -> Parental-UNC    (treatment in parental = control)
#
#   reversal_status (migration switchers only):
#     unresponsive : UNC left the migrated bin in its new compartment (b10unc Stable)
#     reverted     : UNC moved it OPPOSITE to migration (back toward parental)
#     reinforced   : UNC moved it the SAME way (further into the migrated state)  (rare)
#
#   parental_switched              : punc_switch != Stable  (UNC did anything in parental)
#   migration_specific_reversal    : reverted AND parental Stable    (the HEADLINE set)
# ================================================================================
def build_reversal_table(df, args):
    # --- the three categorical tokens ---
    df["mig_switch"]    = classify_axis(df, *args.mig_cols)
    df["punc_switch"]   = classify_axis(df, *args.punc_cols)
    df["b10unc_switch"] = classify_axis(df, *args.b10unc_cols)

    # --- baseline compartment class (for stable_A / stable_B gene lists) ---
    par = pd.to_numeric(df[args.pc1_parental_dmso], errors="coerce")
    b10 = pd.to_numeric(df[args.pc1_b10_dmso], errors="coerce")
    df["baseline_class"] = np.where(
        (par > 0) & (b10 > 0), "Stable_A",
        np.where((par < 0) & (b10 < 0), "Stable_B",
                 np.where(df["mig_switch"] == "A_to_B", "A_to_B",
                          np.where(df["mig_switch"] == "B_to_A", "B_to_A", "Other"))))

    # --- reversal status from migration + B10 treatment tokens ---
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
        return "other"                      # unreachable for switchers (kept for safety)
    df["reversal_status"] = df.apply(reversal_status, axis=1)

    # --- parental switched, and the HEADLINE migration-specific reversal set ---
    df["parental_switched"] = df["punc_switch"] != "Stable"
    df["migration_specific_reversal"] = (
        (df["reversal_status"] == "reverted") &        # reverted in B10
        (df["punc_switch"] == "Stable")                # parental Stable under UNC
    ).fillna(False)

    # --- continuous deltas (for the robustness/attractor figures) ---
    pu  = pd.to_numeric(df[args.pc1_parental_unc], errors="coerce")
    pdm = pd.to_numeric(df[args.pc1_parental_dmso], errors="coerce")
    bu  = pd.to_numeric(df[args.pc1_b10_unc], errors="coerce")
    bdm = pd.to_numeric(df[args.pc1_b10_dmso], errors="coerce")
    df["dPC1_migration"] = b10 - par      # Bottom10-DMSO - Parental-DMSO
    df["dPC1_b10_treat"] = bu  - bdm      # Bottom10-UNC  - Bottom10-DMSO
    df["dPC1_par_treat"] = pu  - pdm      # Parental-UNC  - Parental-DMSO

    # --- chromosome label for block-bootstrap blocking ---
    if args.chrom_col and args.chrom_col in df.columns:
        df["_chrom"] = df[args.chrom_col].astype(str).str.strip()
    else:
        df["_chrom"] = df[args.coord_col].astype(str).map(chrom_of)
    return df


def anti_migration_mask(sw):
    """Parental UNC switched the bin OPPOSITE to its migration direction (the same
    'reverting' move, but in parental). Used ONLY as the matched control for the
    specificity odds ratio. A near-zero rate here is the proof that the reversal
    EVENT is specific to Bottom10."""
    m, t = sw["mig_switch"], sw["punc_switch"]
    return (((m == "A_to_B") & (t == "B_to_A")) |
            ((m == "B_to_A") & (t == "A_to_B")))


# ================================================================================
# FIGURE 1 - braided alluvial of migration switchers (categorical)
#   Each band follows one bin's compartment PATH across the 3 stages. A band that
#   leaves its parental compartment and RETURNS (A->B->A or B->A->B) is a reversal
#   (green); band widths equal the heatmap counts.
# ================================================================================
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

    def _states(r):
        s1, s2 = ("A", "B") if r["mig_switch"] == "A_to_B" else ("B", "A")
        t = r["b10unc_switch"]
        s3 = "B" if t == "A_to_B" else "A" if t == "B_to_A" else s2
        return pd.Series([s1, s2, s3])
    sw[["s1", "s2", "s3"]] = sw.apply(_states, axis=1)
    N = len(sw)

    flows = (sw.groupby(["s1", "s2", "s3"]).size().reset_index(name="n")
               .sort_values(["s1", "s2", "s3"]).reset_index(drop=True))
    flows["type"] = np.where(flows["s3"] == flows["s1"], "reverted", "stayed")

    state_cols = ["s1", "s2", "s3"]
    stage_names = ["Parental", "Bottom10", "Bottom10+UNC"]
    gap = 0.06
    order = list(flows.index)

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
        al  = 0.85 if flows.loc[fi, "type"] == "reverted" else 0.35
        for si in (0, 1):
            y0a, y0b = yint[(fi, si)]; y1a, y1b = yint[(fi, si + 1)]
            _alluvial_ribbon(ax, si + 0.03, y0a, y0b, si + 1 - 0.03, y1a, y1b, col, al)

    ax.set_xlim(-0.35, 2.35); ax.set_ylim(-0.04, 1 + gap + 0.06)
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(stage_names, fontsize=11)
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.legend(handles=[
        Patch(facecolor=COL["revert"], alpha=0.85, label="reverted (returned to parental compartment)"),
        Patch(facecolor=COL["stay"],   alpha=0.5,  label="stayed in migrated compartment")],
        frameon=False, fontsize=8.5, loc="upper center",
        bbox_to_anchor=(0.5, -0.06), ncol=2)
    ax.set_title("Compartment-state flow of migration switchers (categorical)",
                 fontweight="bold", fontsize=11, pad=10)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"Fig_alluvial_states.{ext}"),
                    bbox_inches="tight", dpi=600)
    plt.close(fig)
    print("  Saved -> Fig_alluvial_states.[pdf|png]")


# ================================================================================
# FIGURE 2 - reversal vs parental-response VENN  (replaces the old stacked bars)
#   Left  circle : bins REVERTED by UNC in Bottom10
#   Right circle : bins SWITCHED by UNC in parental (any direction)
#   Overlap      : reverted AND parental-switched  (removed from the headline set)
#   Left-only    : reverted with parental STABLE  =  migration-specific reversal
#   The left-only crescent IS the subtraction (reverted - overlap), so the figure
#   makes it impossible to mis-read it as (reverted - all-parental-switched).
# ================================================================================
def _draw_venn_manual(ax, left_only, right_only, both, n_left, n_right,
                      label_left, label_right):
    """Dependency-free two-circle Venn (used only if matplotlib-venn is absent)."""
    r, cx_l, cx_r, cy = 1.0, -0.55, 0.55, 0.0
    ax.add_patch(Circle((cx_l, cy), r, facecolor=COL["revert"], alpha=0.50,
                        edgecolor="#444", lw=1.3, zorder=1))
    ax.add_patch(Circle((cx_r, cy), r, facecolor=COL["switch"], alpha=0.50,
                        edgecolor="#444", lw=1.3, zorder=1))
    ax.text(cx_l - 0.42, cy, f"{left_only}", ha="center", va="center",
            fontsize=16, fontweight="bold", color="white", zorder=3)
    ax.text(cx_r + 0.42, cy, f"{right_only}", ha="center", va="center",
            fontsize=15, fontweight="bold", zorder=3)
    ax.text(0.0, cy, f"{both}", ha="center", va="center",
            fontsize=13, fontweight="bold", zorder=3)
    ax.text(cx_l - 0.55, cy + r + 0.22, f"{label_left}\n(n={n_left})", ha="center",
            va="bottom", fontsize=9.5, fontweight="bold")
    ax.text(cx_r + 0.55, cy + r + 0.22, f"{label_right}\n(n={n_right})", ha="center",
            va="bottom", fontsize=9.5, fontweight="bold")
    ax.set_xlim(-2.2, 2.2); ax.set_ylim(-1.7, 2.1); ax.set_aspect("equal"); ax.axis("off")
    return (cx_l - 0.42, cy)          # data coords of the left-only label (for the arrow)


def fig_reversal_venn(df, outdir):
    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])]
    reverted = (sw["reversal_status"] == "reverted")
    par_sw   = (sw["punc_switch"] != "Stable")
    n_rev  = int(reverted.sum())
    n_par  = int(par_sw.sum())
    n_both = int((reverted & par_sw).sum())
    left_only  = n_rev - n_both                       # migration-specific reversal
    right_only = n_par - n_both
    n_sw = int(len(sw))
    n_neither = n_sw - (left_only + right_only + n_both)
    lab_left, lab_right = "Reverted by UNC in B10", "Switched by UNC in parental"

    fig, ax = plt.subplots(figsize=(7.8, 6.4))
    try:
        from matplotlib_venn import venn2, venn2_circles
        v = venn2(subsets=(left_only, right_only, n_both),
                  set_labels=(f"{lab_left}\n(n={n_rev})", f"{lab_right}\n(n={n_par})"),
                  set_colors=(COL["revert"], COL["switch"]), alpha=0.55, ax=ax)
        venn2_circles(subsets=(left_only, right_only, n_both), lw=1.3, color="#444", ax=ax)
        if v.get_patch_by_id("10"):
            v.get_patch_by_id("10").set_alpha(0.78)          # highlight the headline region
        if v.get_patch_by_id("11"):
            v.get_patch_by_id("11").set_color(COL["overlap"])   # overlap gets its own colour
            v.get_patch_by_id("11").set_alpha(0.90)
        for idd in ("10", "01", "11"):
            lab = v.get_label_by_id(idd)
            if lab:
                lab.set_fontsize(14); lab.set_fontweight("bold")
        if v.get_label_by_id("10"):
            v.get_label_by_id("10").set_color("white")        # 437 on teal
        if v.get_label_by_id("11"):
            v.get_label_by_id("11").set_color("white")        # overlap count on purple
        for t in v.set_labels:
            if t:
                t.set_fontsize(10.5); t.set_fontweight("bold")
    except Exception:
        _draw_venn_manual(ax, left_only, right_only, n_both,
                          n_rev, n_par, lab_left, lab_right)

    ax.set_title("Reversal vs parental UNC response\n"
                 f"{n_sw:,} migration switchers  ({n_neither:,} switched in neither)",
                 fontweight="bold", fontsize=11)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"Fig_reversal_venn.{ext}"),
                    bbox_inches="tight", dpi=600)
    plt.close(fig)
    print(f"  Saved -> Fig_reversal_venn.[pdf|png]  (reverted {n_rev}, "
          f"parental-switched {n_par}, overlap {n_both}, migration-specific {left_only})")
    return dict(reverted=n_rev, par_switched=n_par, overlap=n_both, specific=left_only)


# ================================================================================
# FIGURE 3a - migration x B10-treatment transition heatmap (2 rows x 3 cols)
#   Descriptive counts; no p-value in the title (per-bin Fisher is pseudoreplicated).
# ================================================================================
def fig_transition_heatmap(df, outdir):
    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])]
    mig = ["A_to_B", "B_to_A"]
    trt = ["A_to_B", "B_to_A", "Stable"]
    M = np.zeros((len(mig), len(trt)), dtype=int)
    for i, m in enumerate(mig):
        for j, t in enumerate(trt):
            M[i, j] = int(((sw["mig_switch"] == m) & (sw["b10unc_switch"] == t)).sum())

    rev  = M[0, 1] + M[1, 0]      # reverted (off-diagonal)
    same = M[0, 0] + M[1, 1]      # reinforced (same direction)
    try:
        _, fisher_p = fisher_exact([[M[0, 1], M[0, 0]], [M[1, 1], M[1, 0]]])
    except Exception:
        fisher_p = np.nan

    fig, ax = plt.subplots(figsize=(6, 4.6))
    im = ax.imshow(M, cmap="Oranges", aspect="auto")          # matches parental heatmap
    for i in range(len(mig)):
        for j in range(len(trt)):
            ax.text(j, i, f"{M[i, j]:,}", ha="center", va="center",
                    fontsize=11, fontweight="bold",
                    color="white" if M[i, j] > M.max() * 0.6 else "black")
    ax.set_xticks(range(len(trt))); ax.set_xticklabels(
        ["B10+UNC\nA\u2192B", "B10+UNC\nB\u2192A", "B10+UNC\nStable"])
    ax.set_yticks(range(len(mig))); ax.set_yticklabels(
        ["migration\nA\u2192B", "migration\nB\u2192A"])
    ax.set_title("Migration switch vs UNC response in B10\n"
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


# ================================================================================
# FIGURE 3b - parental UNC-response heatmap (1 row x 3 cols)
#   Direct counts of how UNC switched the SAME switcher bins in parental cells.
#   Stable = the parental copy did not move = the basis of the migration-specific set.
# ================================================================================
def fig_parental_heatmap(df, outdir):
    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])]
    cols = ["A_to_B", "B_to_A", "Stable"]
    counts = np.array([[int((sw["punc_switch"] == t).sum()) for t in cols]])
    n_sw = int(len(sw))
    n_switched = int(counts[0, 0] + counts[0, 1])
    n_stable = int(counts[0, 2])

    fig, ax = plt.subplots(figsize=(6.4, 2.5))
    im = ax.imshow(counts, cmap="Oranges", aspect="auto", vmin=0, vmax=max(1, counts.max()))
    for j, t in enumerate(cols):
        c = counts[0, j]
        pct = c / n_sw * 100 if n_sw else 0
        ax.text(j, 0, f"{c:,}\n({pct:.1f}%)", ha="center", va="center",
                fontsize=10, fontweight="bold",
                color="white" if c > counts.max() * 0.6 else "black")
    ax.set_xticks(range(3)); ax.set_xticklabels(
        ["Parental+UNC\nA\u2192B", "Parental+UNC\nB\u2192A", "Parental+UNC\nStable"])
    ax.set_yticks([0]); ax.set_yticklabels(["all migration\nswitchers"])
    ax.set_title("UNC response in PARENTAL cells (same switcher bins)\n"
                 f"switched={n_switched:,}   |   Stable={n_stable:,}  of {n_sw:,}",
                 fontweight="bold", fontsize=9.5, pad=8)
    plt.colorbar(im, ax=ax, label="bin count", fraction=0.046, pad=0.04)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"Fig_parental_response_heatmap.{ext}"),
                    bbox_inches="tight", dpi=600)
    plt.close(fig)
    print("  Saved -> Fig_parental_response_heatmap.[pdf|png]")
    return n_switched, n_stable, n_sw


# ================================================================================
# CONTINUOUS ROBUSTNESS (attractor-removed)  --  NOT the headline
#   Provided only to show the categorical result is not a regression-to-the-mean
#   artifact. Figures 4-6 and the `attractor_removed_reversal` column all live here.
# ================================================================================
def fig_dpc1_scatter(df, outdir):
    """RAW continuous dPC1 scatter (read the artifact diagnostic for the de-RTM'd view)."""
    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])].copy()
    sw = sw.dropna(subset=["dPC1_migration", "dPC1_b10_treat", "dPC1_par_treat"])
    if len(sw) < 10:
        print("  [WARN] too few bins for dPC1 scatter; skipping.")
        return np.nan, np.nan, np.nan, np.nan

    rho_b10, p_b10 = spearmanr(sw["dPC1_migration"], sw["dPC1_b10_treat"])
    rho_par, p_par = spearmanr(sw["dPC1_migration"], sw["dPC1_par_treat"])

    allx = sw["dPC1_migration"].values
    ally = np.concatenate([sw["dPC1_b10_treat"].values, sw["dPC1_par_treat"].values])
    xlim = np.nanpercentile(np.abs(allx), 99); ylim = np.nanpercentile(np.abs(ally), 99)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.4), sharex=True, sharey=True)
    panels = [
        (axes[0], "dPC1_b10_treat", rho_b10, "Reversal test:  B10 + UNC",
         "ΔPC1 treatment  (Bottom10+UNC \u2212 Bottom10)"),
        (axes[1], "dPC1_par_treat", rho_par, "Specificity control:  Parental + UNC",
         "ΔPC1 treatment  (Parental+UNC \u2212 Parental)"),
    ]
    for ax, ycol, rho, title, ylab in panels:
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
          f"parental \u03c1={rho_par:+.2f})")
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


def fig_artifact_diagnostic(df, outdir, args):
    """RTM / partial-correlation / convergence diagnostic (continuous robustness)."""
    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])].copy()
    par_b = pd.to_numeric(sw[args.pc1_parental_dmso], errors="coerce")
    b10_b = pd.to_numeric(sw[args.pc1_b10_dmso], errors="coerce")
    paru  = pd.to_numeric(sw[args.pc1_parental_unc], errors="coerce")
    b10u  = pd.to_numeric(sw[args.pc1_b10_unc], errors="coerce")
    dMig  = sw["dPC1_migration"]; dB10 = sw["dPC1_b10_treat"]; dPar = sw["dPC1_par_treat"]

    def sp(a, b):
        m = np.isfinite(a) & np.isfinite(b)
        return spearmanr(a[m], b[m])[0] if m.sum() > 10 else np.nan

    rho_rtm_b10 = sp(dB10, b10_b); rho_rtm_par = sp(dPar, par_b)
    raw_rev  = sp(dMig, dB10);     raw_ctrl = sp(dMig, dPar)
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
    axes[0].set_title(f"RTM / attractor (B10)\n\u03c1={rho_rtm_b10:+.2f}",
                      fontsize=9.5, fontweight="bold")
    for d, c in [("A_to_B", COL["AtoB"]), ("B_to_A", COL["BtoA"])]:
        m = sw["mig_switch"] == d
        axes[1].scatter(par_b[m], dPar[m], s=10, alpha=0.45, color=c, edgecolors="none")
    axes[1].axhline(0, color="#888", lw=0.7)
    axes[1].set_xlabel("Parental baseline PC1 (DMSO)")
    axes[1].set_ylabel("\u0394PC1 parental treatment (Par+UNC \u2212 Par)")
    axes[1].set_title(f"RTM / attractor (Parental)\n\u03c1={rho_rtm_par:+.2f}",
                      fontsize=9.5, fontweight="bold")

    labels = ["Reversal\nraw", "Reversal\n| B10 base", "Reversal\n| Par base",
              "Control\nraw", "Control\n| Par base"]
    vals = [raw_rev, part_rev_b10, part_rev_par, raw_ctrl, part_ctrl_par]
    cols = [COL["BtoA"] if (v is not None and v < 0) else COL["AtoB"] for v in vals]
    axes[2].bar(range(len(vals)), vals, color=cols, edgecolor="white", width=0.7)
    axes[2].axhline(0, color="#333", lw=1.0)
    for i, v in enumerate(vals):
        if v is None or not np.isfinite(v):
            continue
        axes[2].text(i, v + (0.03 if v >= 0 else -0.03), f"{v:+.2f}",
                     ha="center", va="bottom" if v >= 0 else "top", fontsize=8.5)
    axes[2].set_xticks(range(len(labels))); axes[2].set_xticklabels(labels, fontsize=8)
    axes[2].set_ylabel("Spearman \u03c1"); axes[2].set_ylim(-1, 1)
    axes[2].set_title("Raw vs PARTIAL correlation", fontsize=9.5, fontweight="bold")

    conv_txt = ("CONVERGES to common state" if d_treated < d_untreated else "does NOT converge")
    fig.suptitle("Artifact diagnostic (continuous robustness)\n"
                 f"Treated-state distance {d_treated:.4f} vs untreated {d_untreated:.4f}  "
                 f"\u2192 {conv_txt}", fontweight="bold", fontsize=11, y=1.06)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"Fig_artifact_diagnostic.{ext}"),
                    bbox_inches="tight", dpi=600)
    plt.close(fig)
    print("  Saved -> Fig_artifact_diagnostic.[pdf|png]")
    return dict(rho_rtm_b10=rho_rtm_b10, rho_rtm_par=rho_rtm_par,
                raw_rev=raw_rev, part_rev_b10=part_rev_b10, part_rev_par=part_rev_par,
                raw_ctrl=raw_ctrl, part_ctrl_par=part_ctrl_par,
                dist_treated=d_treated, dist_untreated=d_untreated)


def add_attractor_removed_flag(df, args):
    """Continuous robustness set. Regress the B10 treatment delta on the B10 DMSO
    baseline; the residual is the treatment move NOT explained by starting level
    (attractor removed). Flag bins whose residual OPPOSES migration AND whose
    parental copy is continuously unresponsive (|dPC1_par| below its median)."""
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

    tau = (float(args.unresp_thresh) if args.unresp_thresh is not None
           else float(np.nanpercentile(np.abs(resid_b10[sw_mask].dropna()), 66)))
    tau_par = float(np.nanpercentile(np.abs(dpar[sw_mask].dropna()), 50))

    resid_reverts   = (np.sign(resid_b10) == -np.sign(dmig)) & (np.abs(resid_b10) > tau)
    parental_unresp = np.abs(dpar) < tau_par
    df["attractor_removed_reversal"] = (sw_mask & resid_reverts & parental_unresp).fillna(False)
    return tau


def fig_attractor_removed(df, outdir, tau, args):
    """Continuous robustness figure: residual reversal scatter + tau sensitivity."""
    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])].copy()
    spec = sw[sw["attractor_removed_reversal"]]
    n_spec = len(spec); n_sw = len(sw)
    rr = sw.dropna(subset=["dPC1_migration", "resid_b10_treat"])
    rho_resid = spearmanr(rr["dPC1_migration"], rr["resid_b10_treat"])[0] if len(rr) > 10 else np.nan

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0))
    other = sw[~sw["attractor_removed_reversal"]]
    axes[0].scatter(other["dPC1_migration"], other["resid_b10_treat"],
                    s=9, alpha=0.25, color="#BBBBBB", edgecolors="none", label="other switchers")
    for d, c in [("A_to_B", COL["AtoB"]), ("B_to_A", COL["BtoA"])]:
        mm = spec["mig_switch"] == d
        axes[0].scatter(spec.loc[mm, "dPC1_migration"], spec.loc[mm, "resid_b10_treat"],
                        s=16, alpha=0.75, color=c, edgecolors="none",
                        label=f"attractor-removed reversal ({d.replace('_', ARROW)})")
    axes[0].axhline(0, color="#888", lw=0.7); axes[0].axvline(0, color="#888", lw=0.7)
    axes[0].set_xlabel("ΔPC1 migration  (Bottom10 \u2212 Parental)")
    axes[0].set_ylabel("RESIDUAL ΔPC1 B10-treat\n(attractor regressed out)")
    axes[0].set_title(f"Attractor-removed reversal (residual \u03c1={rho_resid:+.2f})\n"
                      f"n={n_spec}/{n_sw} ({n_spec/max(n_sw,1):.1%})",
                      fontsize=9.5, fontweight="bold")
    axes[0].legend(frameon=False, fontsize=7.5, loc="best")

    base = tau if tau > 0 else 0.005
    taus = np.linspace(0.4 * base, 2.0 * base, 12)
    dmig = sw["dPC1_migration"]; resid = sw["resid_b10_treat"]; dpar = sw["dPC1_par_treat"]
    tau_par = float(np.nanpercentile(np.abs(dpar.dropna()), 50))
    pu = np.abs(dpar) < tau_par
    sizes = [int((((np.sign(resid) == -np.sign(dmig)) & (np.abs(resid) > t)) & pu).sum())
             for t in taus]
    axes[1].plot(taus, sizes, "-o", color=COL["revert"], lw=1.6, ms=4)
    axes[1].axvline(tau, color="#333", ls="--", lw=1.0, label=f"chosen tau={tau:.4f}")
    axes[1].set_xlabel("residual threshold tau"); axes[1].set_ylabel("No. of bins")
    axes[1].set_title("Subset-size sensitivity to tau\n(smooth = not cherry-picked)",
                      fontsize=9.5, fontweight="bold")
    axes[1].legend(frameon=False, fontsize=8.5)

    fig.suptitle("Continuous robustness: attractor-removed reversal\n"
                 "(regression-to-the-mean control; NOT the headline set)",
                 fontweight="bold", fontsize=11, y=1.04)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"Fig_attractor_removed_robustness.{ext}"),
                    bbox_inches="tight", dpi=600)
    plt.close(fig)
    print(f"  Saved -> Fig_attractor_removed_robustness.[pdf|png]  "
          f"(n={n_spec}/{n_sw}, residual \u03c1={rho_resid:+.2f})")
    return dict(n_specific=n_spec, n_switchers=n_sw, tau=tau, rho_resid=rho_resid)


# ================================================================================
# CHROMOSOME BLOCK-BOOTSTRAP  (the headline inference)
# --------------------------------------------------------------------------------
# 250 kb bins are autocorrelated, so we resample WHOLE CHROMOSOMES with replacement
# (preserving within-chromosome correlation) and recompute every statistic on each
# resample -> percentile 95% CIs. Effective N = number of chromosomes (~22).
#
# Statistics (over migration switchers):
#   reversal_rate   : P(reverted | switcher)
#   par_switch_rate : P(parental switched | switcher)
#   specific_rate   : P(reverted AND parental Stable | switcher)   (the headline)
#   odds_ratio      : Haldane OR, reverting-in-B10 vs anti-migration-switch-in-parental
#                     (matched reverting event = the specificity test)
#   residual_rho    : Spearman(dPC1_migration, attractor-removed residual); the attractor
#                     regression is REFIT inside each resample (CI includes fit uncertainty)
# ================================================================================
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
    """Refit the attractor (db10 ~ b10base) on THIS resample, then residual rho."""
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
    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])].copy()
    n_sw = len(sw)
    if n_sw < 10:
        print("  [WARN] too few switchers for block-bootstrap; skipping.")
        return None

    chrom        = sw["_chrom"].to_numpy()
    is_reverted  = (sw["reversal_status"] == "reverted").to_numpy()
    is_par_sw    = (sw["punc_switch"] != "Stable").to_numpy()
    is_specific  = (is_reverted & ~is_par_sw)                   # reverted AND parental Stable
    is_anti_mig  = anti_migration_mask(sw).to_numpy()           # matched OR control
    dmig    = pd.to_numeric(sw["dPC1_migration"], errors="coerce").to_numpy(float)
    db10    = pd.to_numeric(sw["dPC1_b10_treat"], errors="coerce").to_numpy(float)
    b10base = pd.to_numeric(sw[args.pc1_b10_dmso], errors="coerce").to_numpy(float)

    uniq   = pd.unique(chrom)
    blocks = {ch: np.where(chrom == ch)[0] for ch in uniq}
    n_blocks = len(uniq)

    def stats_on(idx):
        rv = is_reverted[idx]; ps = is_par_sw[idx]; sp_ = is_specific[idx]; am = is_anti_mig[idx]
        a = int(rv.sum()); b = int((~rv).sum())          # reverting in B10
        c = int(am.sum()); d = int((~am).sum())          # anti-migration switch in parental
        return (rv.mean(), ps.mean(), sp_.mean(),
                _haldane_or(a, b, c, d),
                _resid_rho(dmig[idx], db10[idx], b10base[idx]))

    obs = stats_on(np.arange(n_sw))
    rng = np.random.default_rng(args.boot_seed)
    boot = np.full((args.n_boot, 5), np.nan)
    for it in range(args.n_boot):
        pick = rng.choice(uniq, size=n_blocks, replace=True)
        idx = np.concatenate([blocks[ch] for ch in pick])
        boot[it] = stats_on(idx)

    names = ["reversal_rate", "par_switch_rate", "specific_rate", "odds_ratio", "residual_rho"]
    ci = {}
    for k, nm in enumerate(names):
        lo, hi = np.nanpercentile(boot[:, k], [2.5, 97.5])
        ci[nm] = (obs[k], lo, hi)
    ci["_n_blocks"] = n_blocks
    ci["_n_boot"] = args.n_boot
    return ci


def fig_block_bootstrap(ci, outdir):
    """Forest plot: the three rates (left) + specificity OR & residual rho (right)."""
    if ci is None:
        return
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2),
                             gridspec_kw={"width_ratios": [2, 1.1]})

    rows = [("reversal_rate",   "Reverted in B10"),
            ("par_switch_rate", "Switched in parental"),
            ("specific_rate",   "Migration-specific reversal")]
    ypos = list(range(len(rows)))[::-1]
    for y, (key, lab) in zip(ypos, rows):
        obs, lo, hi = ci[key]
        colr = COL["revert"] if key == "specific_rate" else "#5B6B73"
        axes[0].errorbar(obs * 100, y,
                         xerr=[[(obs - lo) * 100], [(hi - obs) * 100]],
                         fmt="o", color=colr, ecolor="#333", capsize=4, ms=7, lw=1.4)
        axes[0].text(obs * 100, y + 0.18,
                     f"{obs*100:.1f}%  [{lo*100:.1f}, {hi*100:.1f}]",
                     ha="center", va="bottom", fontsize=8.5)
    axes[0].axvline(0, color="#888", lw=0.8, ls="--")
    axes[0].set_yticks(ypos); axes[0].set_yticklabels([l for _, l in rows])
    axes[0].set_xlabel("percent of migration switchers (95% CI)")
    axes[0].set_title("Reversal, parental switching, and the specific set",
                      fontweight="bold", fontsize=10)
    axes[0].set_ylim(-0.6, len(rows) - 0.2)

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
    o, olo, ohi = ci["odds_ratio"]
    axes[1].set_title(f"Specificity odds ratio\n(reverting in B10 vs same switch in parental)\n"
                      f"OR={o:.1f} [{olo:.1f}, {ohi:.1f}]",
                      fontweight="bold", fontsize=9)

    fig.suptitle("Chromosome block-bootstrap 95% CIs  "
                 f"({ci['_n_blocks']} chromosomes, {ci['_n_boot']} resamples)\n"
                 "Blocks = chromosomes, NOT individual bins",
                 fontweight="bold", fontsize=11, y=1.07)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"Fig_block_bootstrap_CI.{ext}"),
                    bbox_inches="tight", dpi=600)
    plt.close(fig)
    print("  Saved -> Fig_block_bootstrap_CI.[pdf|png]")


# ================================================================================
# Gene-list export (enrichment-ready)
#   background_universe          : every gene with a bin in the table (enrichment universe)
#   switcher_A_to_B / _B_to_A    : migration switchers, by direction
#   reverted_under_UNC           : ALL reverted bins
#   migration_specific_reversal  : the HEADLINE set (reverted AND parental Stable)
#   attractor_removed_reversal_robustness : continuous robustness (optional flag)
# ================================================================================
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

    keep = [cc, "mig_switch", "reversal_status", "punc_switch", "baseline_class",
            "migration_specific_reversal"]
    if "attractor_removed_reversal" in df.columns:
        keep.append("attractor_removed_reversal")
    g = genes.merge(df[keep], on=cc, how="inner")

    def write_set(symbols, name):
        syms = sorted(set(s for s in symbols if s and s.lower() != "nan"))
        path = os.path.join(gdir, f"{name}.txt")
        with open(path, "w") as fh:
            fh.write("\n".join(syms) + ("\n" if syms else ""))
        print(f"  gene set [{name:38s}]: {len(syms):5d} genes -> {path}")
        return len(syms)

    write_set(g[gn], "background_universe")
    write_set(g.loc[g["mig_switch"] == "A_to_B", gn], "switcher_A_to_B")
    write_set(g.loc[g["mig_switch"] == "B_to_A", gn], "switcher_B_to_A")
    write_set(g.loc[g["reversal_status"] == "reverted", gn], "reverted_under_UNC")
    write_set(g.loc[g["migration_specific_reversal"] == True, gn],
              "migration_specific_reversal")

    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])]
    n_rev_bins  = int((sw["reversal_status"] == "reverted").sum())
    n_spec_bins = int(sw["migration_specific_reversal"].sum())
    pct = (n_spec_bins / n_rev_bins * 100) if n_rev_bins else 0
    print(f"  [INFO] migration-specific reversal: {n_spec_bins} bins / {n_rev_bins} reverted "
          f"({pct:.0f}% of reversals have a Stable parental copy)")

    if args.include_robustness and "attractor_removed_reversal" in g.columns:
        write_set(g.loc[g["attractor_removed_reversal"] == True, gn],
                  "attractor_removed_reversal_robustness")
    if args.include_unresponsive:
        write_set(g.loc[g["reversal_status"] == "unresponsive", gn],
                  "switcher_unresponsive_to_UNC")
    if args.include_stable:
        write_set(g.loc[g["baseline_class"] == "Stable_A", gn], "stable_A")
        write_set(g.loc[g["baseline_class"] == "Stable_B", gn], "stable_B")

    audit = [c for c in [gn, cc, "mig_switch", "reversal_status", "punc_switch",
                         "baseline_class", "migration_specific_reversal",
                         "attractor_removed_reversal"] if c in g.columns]
    g[audit].to_csv(os.path.join(gdir, "gene_assignments.txt"), sep="\t", index=False)
    print(f"  Per-gene assignments -> {os.path.join(gdir, 'gene_assignments.txt')}")


# ================================================================================
# Main
# ================================================================================
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
    tau = add_attractor_removed_flag(df, args)            # continuous robustness flag

    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])]
    n_sw = len(sw)
    n_rev    = int((sw["reversal_status"] == "reverted").sum())
    n_unres  = int((sw["reversal_status"] == "unresponsive").sum())
    n_parsw  = int((sw["punc_switch"] != "Stable").sum())
    n_overlap = int(((sw["reversal_status"] == "reverted") &
                     (sw["punc_switch"] != "Stable")).sum())
    n_spec   = int(sw["migration_specific_reversal"].sum())
    pct_stable = (n_spec / n_rev * 100) if n_rev else 0

    print(f"[INFO] Migration switchers             : {n_sw:,}  "
          f"(across {df['_chrom'].nunique()} chromosomes)")
    print(f"[INFO] Reverted by UNC in B10          : {n_rev:,}  ({n_rev/max(n_sw,1):.1%})")
    print(f"[INFO] Unresponsive in B10             : {n_unres:,}  ({n_unres/max(n_sw,1):.1%})")
    print(f"[INFO] Switched by UNC in parental     : {n_parsw:,}  ({n_parsw/max(n_sw,1):.1%})")
    print(f"[INFO]   ...reverted AND parental-switched (overlap to remove): {n_overlap:,}")
    print(f"[INFO] Migration-specific reversal     : {n_spec:,}  "
          f"({n_spec/max(n_sw,1):.1%})   = {n_rev} reverted \u2212 {n_overlap} overlap")
    print(f"[INFO]   (reverted in B10 with parental Stable; "
          f"{n_spec}/{n_rev} = {pct_stable:.0f}% of reversals)")

    df.to_csv(os.path.join(args.outdir, "Reversal_classified.txt"), sep="\t", index=False)

    # ---- figures ----------------------------------------------------------
    print("[INFO] Figures...")
    fig_alluvial(df, args.outdir, args)
    fig_reversal_venn(df, args.outdir)
    rev, same, fisher_p = fig_transition_heatmap(df, args.outdir)
    fig_parental_heatmap(df, args.outdir)
    rho_b10, sp_p_b10, rho_par, sp_p_par = fig_dpc1_scatter(df, args.outdir)   # robustness
    diag = fig_artifact_diagnostic(df, args.outdir, args)                      # robustness
    rob  = fig_attractor_removed(df, args.outdir, tau, args)                   # robustness

    # ---- headline inference: chromosome block-bootstrap -------------------
    print("[INFO] Chromosome block-bootstrap (blocks = chromosomes)...")
    ci = block_bootstrap(df, args)
    fig_block_bootstrap(ci, args.outdir)
    if ci is not None:
        for nm, lab in [("reversal_rate", "Reverted in B10"),
                        ("par_switch_rate", "Switched in parental"),
                        ("specific_rate", "Migration-specific reversal"),
                        ("odds_ratio", "Specificity odds ratio"),
                        ("residual_rho", "Residual rho (robustness)")]:
            o, lo, hi = ci[nm]
            if nm in ("reversal_rate", "par_switch_rate", "specific_rate"):
                print(f"  [BOOT] {lab:30s} {o*100:6.1f}%  95% CI [{lo*100:.1f}, {hi*100:.1f}]")
            elif nm == "odds_ratio":
                print(f"  [BOOT] {lab:30s} {o:6.2f}   95% CI [{lo:.2f}, {hi:.2f}]")
            else:
                print(f"  [BOOT] {lab:30s} {o:+6.2f}   95% CI [{lo:+.2f}, {hi:+.2f}]")

    # ---- per-bin pseudoreplicated p-values (NOT for inference) ------------
    a = n_rev; b = n_sw - n_rev
    is_anti = anti_migration_mask(sw)
    c = int(is_anti.sum()); d_ = n_sw - c
    try:
        ps_or, ps_p = fisher_exact([[a, b], [c, d_]])
    except Exception:
        ps_or, ps_p = np.nan, np.nan

    # ---- stats file -------------------------------------------------------
    with open(os.path.join(args.outdir, "Reversal_stats.txt"), "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"migration_switchers\t{n_sw}\n")
        fh.write(f"n_chromosomes\t{df['_chrom'].nunique()}\n")

        fh.write("# === HEADLINE: migration-specific reversal + chromosome block-bootstrap 95% CI ===\n")
        fh.write(f"reverted_bins\t{n_rev}\n")
        fh.write(f"parental_switched_bins\t{n_parsw}\n")
        fh.write(f"reverted_and_parental_switched_bins\t{n_overlap}\n")
        fh.write(f"migration_specific_reversal_bins\t{n_spec}\n")
        fh.write(f"frac_reversals_with_stable_parental\t{pct_stable/100:.4f}\n")
        if ci is not None:
            for nm in ("reversal_rate", "par_switch_rate", "specific_rate",
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

        fh.write("# --- continuous robustness (attractor-removed) ---\n")
        fh.write(f"robustness_attractor_removed_tau\t{rob['tau']:.5f}\n")
        fh.write(f"robustness_attractor_removed_n_bins\t{rob['n_specific']}\n")
        fh.write(f"robustness_residual_spearman_rho\t{rob['rho_resid']:.4f}\n")
        fh.write(f"robustness_RTM_rho_B10\t{diag['rho_rtm_b10']:.4f}\n")
        fh.write(f"robustness_RTM_rho_parental\t{diag['rho_rtm_par']:.4f}\n")
        fh.write(f"robustness_reversal_rho_raw\t{diag['raw_rev']:.4f}\n")
        fh.write(f"robustness_reversal_rho_partial_B10base\t{diag['part_rev_b10']:.4f}\n")

        fh.write("# === PSEUDOREPLICATED bin-level p-values - NOT FOR INFERENCE ===\n")
        fh.write("# (250kb bins are autocorrelated; use the block-bootstrap CIs above.)\n")
        fh.write(f"PSEUDOREP_specificity_fisher_OR\t{ps_or:.4f}\n")
        fh.write(f"PSEUDOREP_specificity_fisher_p\t{ps_p:.3e}\n")
        fh.write(f"PSEUDOREP_transition_fisher_p\t{fisher_p:.3e}\n")
        fh.write(f"PSEUDOREP_dPC1_spearman_rho_B10treat\t{rho_b10:.4f}\n")
        fh.write(f"PSEUDOREP_dPC1_spearman_rho_parental\t{rho_par:.4f}\n")
    print(f"[INFO] Stats -> {os.path.join(args.outdir, 'Reversal_stats.txt')}")

    if args.gene_list:
        print("[INFO] Exporting enrichment-ready gene lists...")
        export_gene_lists(df, args, args.outdir)
    else:
        print("[INFO] --gene_list not given; skipping gene-list export.")

    print(f"\n[INFO] Done -> {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()