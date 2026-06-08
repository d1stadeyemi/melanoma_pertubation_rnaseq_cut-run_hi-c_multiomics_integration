"""
=============================================================================
compartment_reversal_analysis.py
  Does UNC (EZH2i) treatment REVERSE migration-acquired compartment switches?
=============================================================================

THE CAUSAL LOGIC (read first)
  Migration (Parental -> Bottom10) drives compartment switching. If H3K27me3
  activity is REQUIRED to maintain those switches, then inhibiting it (UNC) in
  Bottom10 cells should push switched bins BACK toward their parental state.
  A clean reversal claim needs TWO things, both built here:
    (1) Migration switches revert under UNC-in-B10  (the reversal), AND
    (2) UNC-in-PARENTAL does NOT move those same bins the same way (specificity
        control — otherwise UNC just perturbs compartments in any cell).
  We also extract treatment-UNRESPONSIVE bins (negative set).

INPUTS (from clean_hic_pc1.py output)
  --pc1   A375-CompartmentsAcrossConditions.txt   needs:
            coordinate, chr, start, end
            numeric PC1:  Parental_orig, Bottom10_orig,
                          UNC-DMSO_Parental, UNC-DMSO_Bottom10
            switch calls (BtoA/AtoB/TowardA/TowardB/BothtoA/BothtoB tokens):
              migration : P->B10_Switch_Orig / _Shift_Orig / _Both_original
              P+UNC     : P_DMSO->UNC_switch / _shift / _Both
              B10+UNC   : B10_DMSO->UNC_switch / _shift / _Both
  --gene_list (optional) hg19 250kb-binned gene list (coordinate + gene_name)
              -> enables enrichment-ready gene-list export.

OUTPUTS (into --outdir)
  Reversal_classified.txt        per-bin: all three switch classes + reversal flag
  Fig_alluvial_states.pdf/png    Parental->B10->B10+UNC compartment flow (Sankey/alluvial)
  Fig_reversal_bars.pdf/png      of migration switchers: %% reverted vs stay vs further,
                                 with the PARENTAL-treatment control beside it
  Fig_transition_heatmap.pdf/png migration-class x treatment-class bin counts (+ Fisher)
  Fig_dPC1_scatter.pdf/png       ΔPC1(migration) vs ΔPC1(B10 treatment); Spearman
  Reversal_stats.txt             counts, reversal rate, specificity, Fisher/Spearman
  genelists/<set>.txt            enrichment-ready symbol lists (if --gene_list)
  genelists/background_universe.txt   the correct enrichment background

USAGE
  python compartment_reversal_analysis.py \
      --pc1 data/processed/A375-CompartmentsAcrossConditions.txt \
      --outdir results/exploratory/reversal \
      --gene_list data/processed/hg19_gene_list_binned.txt
  # add stable + unresponsive gene sets:
  python compartment_reversal_analysis.py ... --gene_list ... --include_stable --include_unresponsive
"""

import argparse
import os
import sys
from itertools import product

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from scipy.stats import fisher_exact, spearmanr, rankdata


# ── palette ──────────────────────────────────────────────────────────────────
COL = {"A": "#C0392B", "B": "#1A5276", "switch": "#E67E22",
       "revert": "#27AE60", "stay": "#95A5A6", "further": "#8E44AD",
       "AtoB": "#1A5276", "BtoA": "#C0392B", "Stable": "#D5D8DC"}


# =============================================================================
# Switch-token parsing  (vocabulary: BtoA/AtoB, TowardA/TowardB, BothtoA/BothtoB)
# =============================================================================
def parse_switch_token(both, switch, shift):
    """Return 'B_to_A' / 'A_to_B' / 'Stable' from a row's three switch columns,
    read in priority order Both -> Switch -> Shift. Same vocabulary as the rest
    of the study, extended with BothtoA/BothtoB."""
    for v in (both, switch, shift):
        s = str(v).strip()
        if s in ("-", "", "nan", "None", "NaN"):
            continue
        # toward-A / activated tokens
        if "BtoA" in s or "BothtoA" in s or "TowardA" in s or s == "toA":
            return "B_to_A"
        # toward-B / silenced tokens
        if "AtoB" in s or "BothtoB" in s or "TowardB" in s or s == "toB":
            return "A_to_B"
    return "Stable"


def classify_axis(df, both_c, switch_c, shift_c):
    """Vectorized-ish classification for one axis (returns a Series)."""
    return df.apply(lambda r: parse_switch_token(
        r.get(both_c, "-"), r.get(switch_c, "-"), r.get(shift_c, "-")), axis=1)


# =============================================================================
# Compartment-state sign helper (for the alluvial axes and ΔPC1 scatter)
# =============================================================================
def pc1_state(series):
    """Map a numeric PC1 column to 'A' (>0) / 'B' (<0) / NA.
    Build as an object Series so string labels and missing values coexist
    (np.where mixing str with np.nan raises a DType error in modern numpy)."""
    s = pd.to_numeric(series, errors="coerce")
    out = pd.Series(pd.NA, index=s.index, dtype="object")
    out[s > 0] = "A"
    out[s < 0] = "B"
    return out.to_numpy()


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pc1", required=True)
    p.add_argument("--outdir", default="results/exploratory/reversal")
    p.add_argument("--gene_list", default=None,
                   help="hg19 250kb-binned gene list (coordinate + gene_name) -> "
                        "enables enrichment-ready gene export.")
    p.add_argument("--coord_col", default="coordinate")
    p.add_argument("--gene_name_col", default="gene_name")
    # which numeric PC1 columns to use for the continuous (ΔPC1) view
    p.add_argument("--pc1_parental", default="Parental_orig")
    p.add_argument("--pc1_b10", default="Bottom10_orig")
    p.add_argument("--pc1_parental_dmso", default="Parental_DMSO")
    p.add_argument("--pc1_b10_dmso", default="Bottom10-DMSO")
    p.add_argument("--pc1_parental_unc", default="Parental_UNC")
    p.add_argument("--pc1_b10_unc", default="Bottom10-UNC")
    # switch-call column groups
    p.add_argument("--mig_cols", nargs=3,
                   default=["P->B10_Both_original", "P->B10_Switch_Orig", "P->B10_Shift_Orig"],
                   metavar=("BOTH", "SWITCH", "SHIFT"))
    p.add_argument("--punc_cols", nargs=3,
                   default=["P_DMSO->UNC_Both", "P_DMSO->UNC_switch", "P_DMSO->UNC_shift"],
                   metavar=("BOTH", "SWITCH", "SHIFT"))
    p.add_argument("--b10unc_cols", nargs=3,
                   default=["B10_DMSO->UNC_Both", "B10_DMSO->UNC_switch", "B10_DMSO->UNC_shift"],
                   metavar=("BOTH", "SWITCH", "SHIFT"))
    # enrichment gene-set flags
    p.add_argument("--include_stable", action="store_true",
                   help="Also export Stable_A / Stable_B gene lists.")
    p.add_argument("--include_unresponsive", action="store_true",
                   help="Also export treatment-UNRESPONSIVE switcher gene list "
                        "(migration switchers that do NOT move under UNC-in-B10).")
    p.add_argument("--unresp_thresh", type=float, default=None,
                   help="|ΔPC1| magnitude threshold defining the migration-SPECIFIC "
                        "reversal subset: parental treatment moves the bin LESS "
                        "than this (unresponsive in parental) while B10 treatment "
                        "moves it MORE than this in the reverting direction. "
                        "Default: auto = 33rd percentile of |ΔPC1 parental-treat| "
                        "among migration switchers. A sensitivity sweep around the "
                        "chosen value is always reported.")
    return p.parse_args()

ARROW = "\u2192"

# =============================================================================
# Core classification
# =============================================================================
def build_reversal_table(df, args):
    cc = args.coord_col

    # --- categorical switch classes per axis -------------------------------
    df["mig_switch"]    = classify_axis(df, *args.mig_cols)      # P -> B10
    df["punc_switch"]   = classify_axis(df, *args.punc_cols)     # P -> P+UNC  (control)
    df["b10unc_switch"] = classify_axis(df, *args.b10unc_cols)   # B10 -> B10+UNC (reversal test)

    # --- baseline stable classes (for gene sets / alluvial start state) ----
    par = pd.to_numeric(df[args.pc1_parental], errors="coerce")
    b10 = pd.to_numeric(df[args.pc1_b10], errors="coerce")
    df["baseline_class"] = np.where(
        (par > 0) & (b10 > 0), "Stable_A",
        np.where((par < 0) & (b10 < 0), "Stable_B",
                 np.where(df["mig_switch"] == "A_to_B", "A_to_B",
                          np.where(df["mig_switch"] == "B_to_A", "B_to_A", "Other"))))

    # --- REVERSAL logic (categorical) --------------------------------------
    # A migration switcher REVERTS if the B10+UNC treatment moves it in the
    # OPPOSITE direction to the migration switch:
    #   migration A_to_B  (went toward B)  -> revert needs B10+UNC = B_to_A (toward A)
    #   migration B_to_A  (went toward A)  -> revert needs B10+UNC = A_to_B (toward B)
    def reversal_status(r):
        m, t = r["mig_switch"], r["b10unc_switch"]
        if m == "Stable":
            return "not_switcher"
        if t == "Stable":
            return "unresponsive"          # switched in migration, untouched by UNC
        if (m == "A_to_B" and t == "B_to_A") or (m == "B_to_A" and t == "A_to_B"):
            return "reverted"
        if m == t:
            return "reinforced"            # UNC pushed further the same way
        return "other"
    df["reversal_status"] = df.apply(reversal_status, axis=1)

    # --- specificity check: did UNC move the SAME bin in PARENTAL cells? ----
    # If a bin reverts in B10+UNC but parental+UNC leaves it alone, the effect
    # is migration-state-specific (the strong claim).
    df["parental_unc_moved"] = df["punc_switch"] != "Stable"

    # --- continuous ΔPC1 (support evidence) --------------------------------
    pu  = pd.to_numeric(df[args.pc1_parental_unc], errors="coerce")
    bu  = pd.to_numeric(df[args.pc1_b10_unc], errors="coerce")
    df["dPC1_migration"]   = b10 - par      # P -> B10
    df["dPC1_b10_treat"]   = bu  - b10      # B10 -> B10+UNC
    df["dPC1_par_treat"]   = pu  - par      # P  -> P+UNC  (control)
    return df


# =============================================================================
# Figure 1 — alluvial (Parental -> B10 -> B10+UNC compartment state)
# =============================================================================
def fig_alluvial(df, outdir, args):
    """Three-axis compartment-state flow. Each bin has an A/B state at each of:
    Parental, Bottom10, Bottom10+UNC. Reversal = a ribbon that leaves and returns."""
    par = pc1_state(pd.to_numeric(df[args.pc1_parental], errors="coerce"))
    b10 = pc1_state(pd.to_numeric(df[args.pc1_b10], errors="coerce"))
    bu  = pc1_state(pd.to_numeric(df[args.pc1_b10_unc], errors="coerce"))
    sdf = pd.DataFrame({"Parental": par, "Bottom10": b10, "Bottom10+UNC": bu}).dropna()

    axes_names = ["Parental", "Bottom10", "Bottom10+UNC"]
    states = ["A", "B"]
    x_pos = {n: i for i, n in enumerate(axes_names)}

    # node y-extents (stacked A on top, B below) per axis
    fig, ax = plt.subplots(figsize=(8, 6))
    gap = 0.04
    node_y = {}
    for n in axes_names:
        counts = sdf[n].value_counts()
        total = counts.sum()
        y0 = 0.0
        for s in states:
            h = counts.get(s, 0) / total
            node_y[(n, s)] = (y0, y0 + h)
            ax.add_patch(patches.Rectangle((x_pos[n] - 0.03, y0), 0.06, h,
                         facecolor=COL[s], edgecolor="white", lw=0.5, zorder=3))
            ax.text(x_pos[n], y0 + h / 2, f"{s}", ha="center", va="center",
                    color="white", fontweight="bold", zorder=4, fontsize=9)
            y0 += h + gap

    # ribbons between consecutive axes
    def ribbons(a_name, b_name):
        sub = sdf.groupby([a_name, b_name]).size().reset_index(name="n")
        total = sub["n"].sum()
        # track running offset within each source/target node
        off_src = {s: node_y[(a_name, s)][0] for s in states}
        off_tgt = {s: node_y[(b_name, s)][0] for s in states}
        for _, row in sub.sort_values([a_name, b_name]).iterrows():
            s_s, s_t, n = row[a_name], row[b_name], row["n"]
            h = n / total
            y0s = off_src[s_s]; y0t = off_tgt[s_t]
            off_src[s_s] += h; off_tgt[s_t] += h
            x0 = x_pos[a_name] + 0.03; x1 = x_pos[b_name] - 0.03
            xm = (x0 + x1) / 2
            verts = [(x0, y0s), (xm, y0s), (xm, y0t), (x1, y0t),
                     (x1, y0t + h), (xm, y0t + h), (xm, y0s + h), (x0, y0s + h), (x0, y0s)]
            codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
                     Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
            # colour ribbon by whether the state CHANGED (highlight switches)
            rib_col = COL["switch"] if s_s != s_t else COL[s_s]
            alpha = 0.75 if s_s != s_t else 0.28
            ax.add_patch(patches.PathPatch(Path(verts, codes), facecolor=rib_col,
                         edgecolor="none", alpha=alpha, zorder=2))
    ribbons("Parental", "Bottom10")
    ribbons("Bottom10", "Bottom10+UNC")

    ax.set_xlim(-0.3, len(axes_names) - 1 + 0.3); ax.set_ylim(-0.05, 1.15)
    ax.set_xticks(list(x_pos.values())); ax.set_xticklabels(axes_names, fontsize=11)
    ax.set_yticks([]); 
    for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
    ax.set_title("Compartment state flow: migration then UNC treatment\n"
                 "(orange ribbons = state change; a leave-then-return = reversal)",
                 fontweight="bold", fontsize=11, pad=10)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"Fig_alluvial_states.{ext}"),
                    bbox_inches="tight", dpi=600)
    plt.close(fig)
    print("  Saved -> Fig_alluvial_states.[pdf|png]")


# =============================================================================
# Figure 2 — reversal bars (B10 treatment) with parental-treatment control
# =============================================================================
def fig_reversal_bars(df, outdir):
    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])]
    cats = ["reverted", "unresponsive", "reinforced", "other"]
    colmap = {"reverted": COL["revert"], "unresponsive": COL["stay"],
              "reinforced": COL["further"], "other": "#BDC3C7"}

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharey=True)
    for ax, (grp, title) in zip(axes, [
            (sw, "B10 + UNC  (reversal test)"),
            (sw, "Parental + UNC  (specificity control)")]):
        if "control" in title.lower():
            # control: classify parental-treatment direction relative to migration
            def ctrl_status(r):
                m, t = r["mig_switch"], r["punc_switch"]
                if t == "Stable": return "unresponsive"
                if (m == "A_to_B" and t == "B_to_A") or (m == "B_to_A" and t == "A_to_B"):
                    return "reverted"
                if m == t: return "reinforced"
                return "other"
            status = grp.apply(ctrl_status, axis=1)
        else:
            status = grp["reversal_status"]
        for direction in ["A_to_B", "B_to_A"]:
            mask = grp["mig_switch"] == direction
            sub = status[mask]
            total = len(sub) if len(sub) else 1
            bottom = 0
            xlab = "A\u2192B\nmigration" if direction == "A_to_B" else "B\u2192A\nmigration"
            xi = 0 if direction == "A_to_B" else 1
            for c in cats:
                frac = (sub == c).sum() / total * 100
                ax.bar(xi, frac, bottom=bottom, color=colmap[c], width=0.6,
                       edgecolor="white", linewidth=0.6,
                       label=c if (xi == 0 and ax is axes[0]) else None)
                if frac > 6:
                    ax.text(xi, bottom + frac / 2, f"{frac:.0f}%", ha="center",
                            va="center", fontsize=8, color="white", fontweight="bold")
                bottom += frac
        ax.set_xticks([0, 1]); ax.set_xticklabels(["A\u2192B", "B\u2192A"])
        ax.set_title(title, fontweight="bold", fontsize=10)
        ax.set_ylim(0, 100)
    axes[0].set_ylabel("% of migration switchers")
    axes[0].legend(frameon=False, fontsize=8.5, loc="upper center",
                   bbox_to_anchor=(0.5, -0.12), ncol=4)
    fig.suptitle("Do migration-acquired switches revert under UNC?\n"
                 "Reversal in B10 (left) should EXCEED the parental control (right) "
                 "for a migration-specific effect",
                 fontweight="bold", fontsize=11, y=1.02)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"Fig_reversal_bars.{ext}"),
                    bbox_inches="tight", dpi=600)
    plt.close(fig)
    print("  Saved -> Fig_reversal_bars.[pdf|png]")


# =============================================================================
# Figure 3 — transition heatmap (migration class x B10-treatment class) + Fisher
# =============================================================================
def fig_transition_heatmap(df, outdir):
    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])]
    mig = ["A_to_B", "B_to_A"]
    trt = ["A_to_B", "B_to_A", "Stable"]
    M = np.zeros((len(mig), len(trt)), dtype=int)
    for i, m in enumerate(mig):
        for j, t in enumerate(trt):
            M[i, j] = ((sw["mig_switch"] == m) & (sw["b10unc_switch"] == t)).sum()

    # Fisher on the 2x2 of reversal vs not (anti-diagonal enrichment)
    rev = M[0, 1] + M[1, 0]                      # A_to_B->B_to_A  +  B_to_A->A_to_B
    same = M[0, 0] + M[1, 1]
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
                 f"(reverted={rev:,}, same-direction={same:,}; Fisher p={fisher_p:.1e})",
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
# Figure 4 — continuous ΔPC1 scatter (migration vs B10 treatment)
# =============================================================================
def fig_dpc1_scatter(df, outdir):
    """Two panels, SHARED axes:
      LEFT  (reversal test)       x=ΔPC1 migration, y=ΔPC1 B10-treatment
      RIGHT (specificity control) x=ΔPC1 migration, y=ΔPC1 PARENTAL-treatment
    The causal logic: the LEFT panel should show a strong NEGATIVE correlation
    (treatment undoes the migration shift). The RIGHT panel is the control — if
    UNC moved these same bins regardless of cell state, it would ALSO be
    negative. Strong-negative LEFT vs near-zero RIGHT is the migration-SPECIFIC
    reversal signal. Reporting only the left panel cannot distinguish the two."""
    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])].copy()
    needed = ["dPC1_migration", "dPC1_b10_treat", "dPC1_par_treat"]
    sw = sw.dropna(subset=needed)
    if len(sw) < 10:
        print("  [WARN] too few bins for ΔPC1 scatter; skipping.")
        return np.nan, np.nan, np.nan, np.nan

    rho_b10, p_b10 = spearmanr(sw["dPC1_migration"], sw["dPC1_b10_treat"])
    rho_par, p_par = spearmanr(sw["dPC1_migration"], sw["dPC1_par_treat"])

    # shared symmetric limits so the two panels are visually comparable
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
        ax.set_title(f"{title}\nSpearman \u03c1={rho:+.2f}  (p={p:.1e})",
                     fontweight="bold", fontsize=10, pad=8)
    axes[0].legend(frameon=False, fontsize=9, loc="best")
    fig.suptitle("Continuous reversal vs specificity control\n"
                 "Strong-negative LEFT + near-zero RIGHT = migration-specific reversal",
                 fontweight="bold", fontsize=11, y=1.04)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"Fig_dPC1_scatter.{ext}"),
                    bbox_inches="tight", dpi=600)
    plt.close(fig)
    print(f"  Saved -> Fig_dPC1_scatter.[pdf|png]  (B10 \u03c1={rho_b10:+.2f}, "
          f"parental-control \u03c1={rho_par:+.2f})")
    return rho_b10, p_b10, rho_par, p_par


# =============================================================================
# Figure 5 — ARTIFACT DIAGNOSTIC: is the ±0.6 real, or shared-variable / RTM?
# =============================================================================
def _partial_spearman(x, y, z):
    """Partial Spearman correlation of x and y, controlling for z.
    Rank-transform all three, then apply the first-order partial-correlation
    formula. Tells us how much of the x~y rank relationship SURVIVES once the
    confounder z is removed."""
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
    """Quantify how much of the reversal (\u03c1\u2248\u22120.6) and control (\u03c1\u2248+0.6)
    correlations is a SHARED-VARIABLE / regression-to-the-mean (attractor)
    artifact rather than migration-specific biology.

    THE CONCERN
      x = \u0394PC1_migration = (B10 \u2212 Par) shares a term with each y:
        reversal y = (B10_UNC \u2212 B10) shares B10 (opposite sign) -> spurious NEG
        control  y = (Par_UNC \u2212 Par) shares Par (opposite sign) -> spurious POS
      The SAME \u2212/+ pattern is also produced if UNC simply pulls PC1 toward a
      common attractor level (regression to the mean), independent of migration.

    THREE DIAGNOSTICS (this figure + printed/saved numbers)
      (A) RTM / attractor: treatment \u0394PC1 vs OWN baseline level. A strong
          NEGATIVE slope means the treatment change is governed by where the bin
          STARTED (pull toward a common level), not by migration history.
      (B) Partial Spearman: does \u0394PC1_migration still predict treatment \u0394PC1
          after controlling for the baseline level? If it COLLAPSES toward 0,
          the apparent reversal was baseline/shared-variable, not migration.
      (C) Convergence: are the two TREATED states (B10_UNC, Par_UNC) closer to
          each other than the two UNTREATED states (B10, Par)? If yes, UNC drives
          both toward a common attractor -> explains \u2212/+ without reversal.
    """
    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])].copy()
    par_b = pd.to_numeric(sw[args.pc1_parental], errors="coerce")
    b10_b = pd.to_numeric(sw[args.pc1_b10], errors="coerce")
    paru  = pd.to_numeric(sw[args.pc1_parental_unc], errors="coerce")
    b10u  = pd.to_numeric(sw[args.pc1_b10_unc], errors="coerce")
    dMig  = sw["dPC1_migration"]
    dB10  = sw["dPC1_b10_treat"]
    dPar  = sw["dPC1_par_treat"]

    def sp(a, b):
        m = np.isfinite(a) & np.isfinite(b)
        return spearmanr(a[m], b[m])[0] if m.sum() > 10 else np.nan

    # (A) RTM: treatment change vs own baseline
    rho_rtm_b10 = sp(dB10, b10_b)
    rho_rtm_par = sp(dPar, par_b)

    # raw correlations (the ones from Fig 4)
    raw_rev  = sp(dMig, dB10)
    raw_ctrl = sp(dMig, dPar)

    # (B) partial correlations controlling for the SHARED baseline term
    part_rev_b10 = _partial_spearman(dMig, dB10, b10_b)   # control the shared B10
    part_rev_par = _partial_spearman(dMig, dB10, par_b)   # control the other baseline
    part_ctrl_par = _partial_spearman(dMig, dPar, par_b)  # control the shared Par

    # (C) convergence to a common attractor
    d_untreated = np.nanmedian(np.abs(b10_b - par_b))
    d_treated   = np.nanmedian(np.abs(b10u - paru))

    # ---- figure: 3 panels ----------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.9))

    # Panel A1: RTM B10
    for d, c in [("A_to_B", COL["AtoB"]), ("B_to_A", COL["BtoA"])]:
        m = sw["mig_switch"] == d
        axes[0].scatter(b10_b[m], dB10[m], s=10, alpha=0.45, color=c, edgecolors="none")
    axes[0].axhline(0, color="#888", lw=0.7)
    axes[0].set_xlabel("Bottom10 baseline PC1")
    axes[0].set_ylabel("\u0394PC1 B10 treatment (B10+UNC \u2212 B10)")
    axes[0].set_title(f"RTM / attractor (B10)\n\u03c1={rho_rtm_b10:+.2f}  "
                      "(neg = pull toward common level)", fontsize=9.5, fontweight="bold")

    # Panel A2: RTM Parental
    for d, c in [("A_to_B", COL["AtoB"]), ("B_to_A", COL["BtoA"])]:
        m = sw["mig_switch"] == d
        axes[1].scatter(par_b[m], dPar[m], s=10, alpha=0.45, color=c, edgecolors="none")
    axes[1].axhline(0, color="#888", lw=0.7)
    axes[1].set_xlabel("Parental baseline PC1")
    axes[1].set_ylabel("\u0394PC1 parental treatment (Par+UNC \u2212 Par)")
    axes[1].set_title(f"RTM / attractor (Parental)\n\u03c1={rho_rtm_par:+.2f}  "
                      "(neg = pull toward common level)", fontsize=9.5, fontweight="bold")

    # Panel B: raw vs partial correlations
    labels = ["Reversal\nraw", "Reversal\n| B10 base", "Reversal\n| Par base",
              "Control\nraw", "Control\n| Par base"]
    vals = [raw_rev, part_rev_b10, part_rev_par, raw_ctrl, part_ctrl_par]
    cols = [COL["BtoA"] if v < 0 else COL["AtoB"] for v in vals]  # neg red, pos blue
    bars = axes[2].bar(range(len(vals)), vals, color=cols, edgecolor="white", width=0.7)
    axes[2].axhline(0, color="#333", lw=1.0)
    for i, v in enumerate(vals):
        axes[2].text(i, v + (0.03 if v >= 0 else -0.03), f"{v:+.2f}",
                     ha="center", va="bottom" if v >= 0 else "top", fontsize=8.5)
    axes[2].set_xticks(range(len(labels)))
    axes[2].set_xticklabels(labels, fontsize=8)
    axes[2].set_ylabel("Spearman \u03c1")
    axes[2].set_ylim(-1, 1)
    axes[2].set_title("Raw vs PARTIAL correlation\n(collapse toward 0 = was artifact)",
                      fontsize=9.5, fontweight="bold")

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
    print(f"  [DIAG] RTM \u03c1: B10={rho_rtm_b10:+.2f}, parental={rho_rtm_par:+.2f} "
          "(strong-neg both = attractor dominates)")
    print(f"  [DIAG] reversal \u03c1 raw={raw_rev:+.2f} -> partial|B10base={part_rev_b10:+.2f}, "
          f"partial|Parbase={part_rev_par:+.2f}")
    print(f"  [DIAG] control  \u03c1 raw={raw_ctrl:+.2f} -> partial|Parbase={part_ctrl_par:+.2f}")
    print(f"  [DIAG] convergence: treated dist {d_treated:.4f} vs untreated "
          f"{d_untreated:.4f} -> {conv_txt}")
    return dict(rho_rtm_b10=rho_rtm_b10, rho_rtm_par=rho_rtm_par,
                raw_rev=raw_rev, part_rev_b10=part_rev_b10, part_rev_par=part_rev_par,
                raw_ctrl=raw_ctrl, part_ctrl_par=part_ctrl_par,
                dist_treated=d_treated, dist_untreated=d_untreated)


# =============================================================================
# Migration-SPECIFIC reversal subset (the clean causal set)
# =============================================================================
def add_specific_reversal_flag(df, args):
    """Isolate MIGRATION-SPECIFIC reversal — reversal that is NOT explained by the
    global attractor (pull toward a common PC1 level).

    WHY A NAIVE CRITERION FAILS
      "parental-unresponsive AND B10-reverting" does NOT work: a bin whose
      parental baseline already sits near the attractor target looks parentally
      unresponsive (already there) yet its B10 is still pulled toward that target,
      mimicking reversal. Those ambiguous bins would contaminate the set.

    THE FIX — attractor residuals
      Fit the attractor directly: regress (ΔPC1 B10-treat) on (B10 baseline level)
      over switchers. The fitted line IS the attractor (RTM slope \u2248 \u22120.93 here).
      The RESIDUAL = the part of the B10 treatment move NOT explained by baseline
      level. A genuine migration-specific reversal appears as a residual that
      OPPOSES the migration shift (moves back toward parental beyond what the
      attractor predicts). We additionally require parental to be unresponsive,
      as a second independent filter."""
    sw_mask = df["mig_switch"].isin(["A_to_B", "B_to_A"])
    dmig = df["dPC1_migration"]; db10 = df["dPC1_b10_treat"]; dpar = df["dPC1_par_treat"]
    b10_base = pd.to_numeric(df[args.pc1_b10], errors="coerce")

    # Fit the attractor (RTM) on switchers and subtract it -> residual B10 move.
    m = sw_mask & np.isfinite(db10) & np.isfinite(b10_base)
    if m.sum() > 10:
        slope, intercept = np.polyfit(b10_base[m].values, db10[m].values, 1)
        pred = slope * b10_base + intercept
    else:
        pred = pd.Series(0.0, index=df.index)
    resid_b10 = db10 - pred
    df["resid_b10_treat"] = resid_b10

    # thresholds (data-driven): residual must be a real, reverting move; parental
    # must be quiet. tau on the residual is the knob the sensitivity sweep varies.
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
    """Two panels: (left) ATTRACTOR-REMOVED view — residual B10 move vs migration
    shift, with the migration-specific subset highlighted (real reversal sits on
    the anti-diagonal HERE, after the attractor is regressed out); (right) a
    threshold sensitivity sweep so the subset size isn't pinned to one tau."""
    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])].copy()
    spec = sw[sw["specific_reversal"]]
    n_spec = len(spec); n_sw = len(sw)

    med_par_move = float(np.nanmedian(np.abs(spec["dPC1_par_treat"]))) if n_spec else np.nan
    med_b10_move = float(np.nanmedian(np.abs(spec["dPC1_b10_treat"]))) if n_spec else np.nan
    # correlation in residual space (attractor removed) — the honest effect
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

    # sensitivity sweep over residual threshold
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
    axes[1].set_ylabel("# migration-specific reversal bins")
    axes[1].set_title("Sensitivity of subset size to tau\n(smooth trend = not cherry-picked)",
                      fontsize=9.5, fontweight="bold")
    axes[1].legend(frameon=False, fontsize=8.5)

    fig.suptitle("Clean migration-specific reversal (attractor regressed out)\n"
                 "residual-reverting \u2227 parental-unresponsive bins",
                 fontweight="bold", fontsize=11, y=1.04)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"Fig_specific_reversal.{ext}"),
                    bbox_inches="tight", dpi=600)
    plt.close(fig)
    print(f"  Saved -> Fig_specific_reversal.[pdf|png]")
    print(f"  [SPECIFIC] migration-specific reversal: {n_spec}/{n_sw} switchers "
          f"({n_spec/max(n_sw,1):.1%}); tau={tau:.4f}; residual \u03c1={rho_resid:+.2f}; "
          f"parental |ΔPC1| median={med_par_move:.4f} (should be small)")
    return dict(n_specific=n_spec, n_switchers=n_sw, tau=tau, rho_resid=rho_resid,
                med_par_move=med_par_move, med_b10_move=med_b10_move)


# =============================================================================
# Gene-list export (enrichment-ready) + correct background universe
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

    # map each gene to its bin's classes
    keep = [cc, "mig_switch", "reversal_status", "baseline_class"]
    if "specific_reversal" in df.columns:
        keep.append("specific_reversal")
    g = genes.merge(df[keep], on=cc, how="inner")

    def write_set(symbols, name):
        syms = sorted(set(s for s in symbols if s and s.lower() != "nan"))
        path = os.path.join(gdir, f"{name}.txt")
        with open(path, "w") as fh:
            fh.write("\n".join(syms) + "\n")
        print(f"  gene set [{name}]: {len(syms)} genes -> {path}")

    # --- BACKGROUND UNIVERSE (critical): all genes assignable to ANY bin ----
    # Switchers are gene-poor, so the enrichment background MUST be the set of
    # genes that could have been assigned — NOT the whole genome — or every
    # term is spuriously enriched.
    write_set(g[gn], "background_universe")

    # --- default sets: switchers + reverted ---------------------------------
    write_set(g.loc[g["mig_switch"] == "A_to_B", gn], "switcher_A_to_B")
    write_set(g.loc[g["mig_switch"] == "B_to_A", gn], "switcher_B_to_A")
    write_set(g.loc[g["reversal_status"] == "reverted", gn], "reverted_under_UNC")

    # --- the CLEAN causal set: migration-specific (attractor-removed) reversal -
    # This is the defensible reversal gene set — bins that revert in B10 AND are
    # unresponsive in parental, so the reversal is NOT the global attractor.
    if "specific_reversal" in g.columns:
        write_set(g.loc[g["specific_reversal"] == True, gn],
                  "migration_specific_reversal")

    # --- flagged sets -------------------------------------------------------
    if args.include_unresponsive:
        write_set(g.loc[g["reversal_status"] == "unresponsive", gn],
                  "switcher_unresponsive_to_UNC")
    if args.include_stable:
        write_set(g.loc[g["baseline_class"] == "Stable_A", gn], "stable_A")
        write_set(g.loc[g["baseline_class"] == "Stable_B", gn], "stable_B")

    # auditable per-gene table
    audit_cols = [gn, cc, "mig_switch", "reversal_status", "baseline_class"]
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

    # verify required columns exist (loud failure beats silent misclassification)
    need = ([args.pc1_parental, args.pc1_b10, args.pc1_parental_unc, args.pc1_b10_unc]
            + args.mig_cols + args.punc_cols + args.b10unc_cols)
    missing = [c for c in need if c not in df.columns]
    if missing:
        sys.exit("[ERROR] Missing required columns:\n  " + "\n  ".join(missing)
                 + "\n\nAvailable:\n  " + "\n  ".join(map(str, df.columns)))

    df = build_reversal_table(df, args)
    tau_spec = add_specific_reversal_flag(df, args)

    # counts
    sw = df[df["mig_switch"].isin(["A_to_B", "B_to_A"])]
    n_sw = len(sw)
    rev_rate = (sw["reversal_status"] == "reverted").mean() if n_sw else np.nan
    unresp_rate = (sw["reversal_status"] == "unresponsive").mean() if n_sw else np.nan
    print(f"[INFO] Migration switchers: {n_sw:,}")
    print(f"[INFO] Reverted under UNC-in-B10: {rev_rate:.1%}")
    print(f"[INFO] Unresponsive to UNC-in-B10: {unresp_rate:.1%}")

    # save classified table
    df.to_csv(os.path.join(args.outdir, "Reversal_classified.txt"),
              sep="\t", index=False)

    # figures
    print("[INFO] Figures...")
    fig_alluvial(df, args.outdir, args)
    fig_reversal_bars(df, args.outdir)
    rev, same, fisher_p = fig_transition_heatmap(df, args.outdir)
    rho_b10, sp_p_b10, rho_par, sp_p_par = fig_dpc1_scatter(df, args.outdir)
    diag = fig_artifact_diagnostic(df, args.outdir, args)
    spec_summary = fig_specific_reversal(df, args.outdir, tau_spec, args)

    # specificity: reversal rate in B10 vs apparent-reversal in parental control
    def ctrl_rev(r):
        m, t = r["mig_switch"], r["punc_switch"]
        return ((m == "A_to_B" and t == "B_to_A") or (m == "B_to_A" and t == "A_to_B"))
    sw_rev_b10 = (sw["reversal_status"] == "reverted")
    sw_rev_par = sw.apply(ctrl_rev, axis=1)
    rev_rate_b10 = sw_rev_b10.mean() if n_sw else np.nan
    ctrl_rate = sw_rev_par.mean() if n_sw else np.nan

    # FORMAL test of EXCESS reversal: is reversal MORE common in B10 than in the
    # parental control? 2x2 = (reverted / not) x (B10 / parental). This is the
    # number that quantifies the migration-SPECIFIC effect — the raw B10 rate
    # alone conflates specific reversal with non-specific UNC movement.
    a = int(sw_rev_b10.sum()); b = int((~sw_rev_b10).sum())
    c = int(sw_rev_par.sum()); d_ = int((~sw_rev_par).sum())
    try:
        excess_OR, excess_p = fisher_exact([[a, b], [c, d_]])
    except Exception:
        excess_OR, excess_p = np.nan, np.nan
    excess_pts = (rev_rate_b10 - ctrl_rate)
    print(f"[INFO] EXCESS reversal (B10 \u2212 parental): "
          f"{excess_pts:+.1%}  (Fisher OR={excess_OR:.2f}, p={excess_p:.1e})")

    with open(os.path.join(args.outdir, "Reversal_stats.txt"), "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"migration_switchers\t{n_sw}\n")
        fh.write(f"reverted_rate_B10UNC\t{rev_rate_b10:.4f}\n")
        fh.write(f"apparent_reverted_rate_parentalUNC_control\t{ctrl_rate:.4f}\n")
        fh.write(f"excess_reversal_points_B10_minus_parental\t{excess_pts:.4f}\n")
        fh.write(f"excess_reversal_fisher_OR\t{excess_OR:.4f}\n")
        fh.write(f"excess_reversal_fisher_p\t{excess_p:.3e}\n")
        fh.write(f"unresponsive_rate_B10UNC\t{unresp_rate:.4f}\n")
        fh.write(f"transition_reverted_bins\t{rev}\n")
        fh.write(f"transition_same_direction_bins\t{same}\n")
        fh.write(f"transition_fisher_p\t{fisher_p:.3e}\n")
        fh.write(f"dPC1_spearman_rho_B10treat\t{rho_b10:.4f}\n")
        fh.write(f"dPC1_spearman_p_B10treat\t{sp_p_b10:.3e}\n")
        fh.write(f"dPC1_spearman_rho_parentalControl\t{rho_par:.4f}\n")
        fh.write(f"dPC1_spearman_p_parentalControl\t{sp_p_par:.3e}\n")
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
    print(f"[INFO] Stats -> {os.path.join(args.outdir, 'Reversal_stats.txt')}")
    print(f"[INFO] Specificity (categorical): B10 reversal {rev_rate_b10:.1%} vs "
          f"parental-control {ctrl_rate:.1%}  -> excess {excess_pts:+.1%}.")
    print(f"[INFO] Specificity (continuous): B10 \u03c1={rho_b10:+.2f} vs "
          f"parental-control \u03c1={rho_par:+.2f}  "
          f"(near-zero control = migration-specific reversal).")

    if args.gene_list:
        print("[INFO] Exporting enrichment-ready gene lists...")
        export_gene_lists(df, args, args.outdir)
    else:
        print("[INFO] --gene_list not given; skipping gene-list export.")

    print(f"\n[INFO] Done -> {os.path.abspath(args.outdir)}")
    print("[NOTE] Compartment/PC1 analysis is spike-in-INDEPENDENT (Hi-C derived). "
          "Any H3K27me3 *amount* comparison is a SEPARATE analysis that REQUIRES "
          "spike-in normalization (see compartment_reversal docs).")


if __name__ == "__main__":
    main()