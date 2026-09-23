import numpy as np
import pandas as pd
import matplotlib


def _in_notebook():
    """True inside a Jupyter/IPython kernel."""
    try:
        from IPython import get_ipython
        ip = get_ipython()
        return ip is not None and "IPKernelApp" in ip.config
    except Exception:
        return False


# Force the headless backend ONLY when running as a script (e.g. via SLURM on
# ISAAC, where there is no display). Doing this unconditionally would stop
# figures rendering inline in a notebook.
if not _in_notebook():
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def _norm(s):
    """Normalise a category label for matching.

    Underscores, hyphens and spaces are interchangeable, and case is ignored,
    so 'No_change', 'No change' and 'no-change' all resolve to the same key.
    This exists because the legend renders labels with underscores replaced by
    spaces: without normalisation a mismatch between the data's convention and
    this module's constants is INVISIBLE -- the legend looks right while the
    colours silently fall back and `exclude` silently fails.
    """
    return str(s).replace("_", " ").replace("-", " ").strip().lower()


# Canonical stacking order: most B-ward at one end, most A-ward at the other,
# with No_change in the middle. Matching is normalisation-insensitive, so the
# underscores here are cosmetic. Any category not listed is appended at the end
# so unexpected labels still appear rather than being silently dropped.
CANONICAL_ORDER = [
    "Strong_Switch_AtoB", "Weak_Switch_AtoB",
    "Strong_Shift_in_A_towardB", "Weak_Shift_in_A_towardB",
    "Strong_Shift_in_B_towardB", "Weak_Shift_in_B_towardB",
    "No_change",
    "Weak_Shift_in_B_towardA", "Strong_Shift_in_B_towardA",
    "Weak_Shift_in_A_towardA", "Strong_Shift_in_A_towardA",
    "Weak_Switch_BtoA", "Strong_Switch_BtoA",
    "Undefined",
]

# Diverging palette. Hue encodes DIRECTION (cool = toward B, warm = toward A);
# lightness encodes MAGNITUDE, darkest at the extremes. Blue/orange rather than
# blue/red: it survives deuteranopia and protanopia, stays legible in greyscale
# because the lightness ramp is monotone, and reads as editorial rather than
# default-matplotlib. Neutral warm grey for No_change keeps the eye on the tails.
CATEGORY_COLORS = {
    # toward B -- cool arm. Deep navy through to teal: the ramp uses BOTH
    # lightness and a hue rotation (blue -> teal) so six steps stay separable
    # without any of them going pale. Nothing here is lighter than mid-tone.
    "Strong_Switch_AtoB":        "#0A2E4D",
    "Weak_Switch_AtoB":          "#14567E",
    "Strong_Shift_in_A_towardB": "#1E7BA6",
    "Weak_Shift_in_A_towardB":   "#3D9CBC",
    "Strong_Shift_in_B_towardB": "#63BAC4",
    "Weak_Shift_in_B_towardB":   "#93CFCB",
    # neutral -- warm mid grey, dark enough to sit in the same tonal family
    "No_change":                 "#CFCBC4",
    # toward A -- warm arm, mirroring the cool one: amber through to dark maroon.
    "Weak_Shift_in_B_towardA":   "#F2C078",
    "Strong_Shift_in_B_towardA": "#E79B49",
    "Weak_Shift_in_A_towardA":   "#D4712F",
    "Strong_Shift_in_A_towardA": "#B54B26",
    "Weak_Switch_BtoA":          "#8C2B1E",
    "Strong_Switch_BtoA":        "#551210",
    # unassignable -- deliberately the flattest colour so it recedes
    "Undefined":                 "#B8B4AD",
}

_ORDER_LOOKUP = {_norm(c): i for i, c in enumerate(CANONICAL_ORDER)}
_COLOR_LOOKUP = {_norm(k): v for k, v in CATEGORY_COLORS.items()}
_FALLBACK = plt.get_cmap("tab20").colors


def _order_categories(found):
    """Canonical order first (matched loosely), then anything unrecognised."""
    known = sorted((c for c in found if _norm(c) in _ORDER_LOOKUP),
                   key=lambda c: _ORDER_LOOKUP[_norm(c)])
    extra = sorted(c for c in found if _norm(c) not in _ORDER_LOOKUP)
    return known + extra

def plot_category_stack(
    df, columns, labels=None, normalize="percent", orientation="vertical",
    title=None, subtitle=None, min_label=4.0, show_n=True, bar_width=0.62,
    figsize=None, legend_ncol=1, palette=None,
    exclude=("No_change", "Undefined"), report="print", label_as="auto",
    category_order=None, ax=None, save=None, dpi=600,
):
    """Stacked bar of category composition, one bar per column.

    df              : DataFrame holding the category columns.
    columns         : list of column names to plot (any number).
    labels          : x tick labels; defaults to the column names.
    normalize       : 'percent' (composition, bars sum to 100) or 'count'.
    orientation     : 'vertical' or 'horizontal'.
    min_label       : minimum size for a segment to get an inline label, as a
                      percent. In 'percent' mode this is a share of the segment's
                      own bar; in 'count' mode it is a share of the TALLEST bar,
                      i.e. physical size on the page, so short bars are not
                      crowded with unreadable text. Set None to disable.
    label_as        : text of the inline labels -- 'percent', 'count', or
                      'auto' (default), which follows `normalize` so a count
                      axis gets count labels and a percent axis gets percents.
    show_n          : annotate each bar with its n.
    exclude         : categories held OUT of the figure and reported separately
                      instead (default: No_change and Undefined). The remaining
                      bars renormalise over the changed categories only, so the
                      y axis reads '% of changed bins'. Pass () to plot
                      everything.
    report          : 'print' writes a plain-text summary of the excluded
                      categories to stdout; 'none' stays silent. The summary is
                      returned either way.
    palette         : dict overriding CATEGORY_COLORS.
    category_order  : explicit stacking order, overriding the canonical one.
    save            : path to write (extension decides format); also returns fig.

    Returns (fig, ax, table, summary):
      table   : counts of the PLOTTED categories, one column per input column
      summary : per-column counts and percentages for the excluded categories
                plus the changed/total totals -- report these in the legend or
                the text rather than in the bars.
    """
    if isinstance(columns, str):
        columns = [columns]
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise KeyError(f"columns not in df: {missing}")
    if normalize not in ("percent", "count"):
        raise ValueError("normalize must be 'percent' or 'count'")
    if orientation not in ("vertical", "horizontal"):
        raise ValueError("orientation must be 'vertical' or 'horizontal'")

    full = pd.DataFrame({c: df[c].value_counts() for c in columns}).fillna(0)
    grand = full.sum(axis=0)

    # Match loosely so 'No_change' in the data is caught by 'No change' here.
    exclude_keys = {_norm(e) for e in (exclude or ())}
    held = [c for c in full.index if _norm(c) in exclude_keys]
    tab = full.drop(index=held) if held else full.copy()

    # Per-column summary of what was held out, plus the totals.
    rows = {}
    for c in held:
        rows[f"{c} (n)"] = full.loc[c]
        rows[f"{c} (% of all bins)"] = (full.loc[c] / grand * 100).round(2)
    rows["changed (n)"] = tab.sum(axis=0)
    rows["changed (% of all bins)"] = (tab.sum(axis=0) / grand * 100).round(2)
    rows["all bins (n)"] = grand
    summary = pd.DataFrame(rows).T[columns]

    if report == "print" and held:
        print(f"Categories held out of the figure: {', '.join(held)}")
        for c in columns:
            parts = ", ".join(
                f"{h.replace('_', ' ')} {int(full.loc[h, c]):,} "
                f"({full.loc[h, c] / grand[c] * 100:.1f}%)" for h in held)
            print(f"  {c}: {parts}; changed {int(tab[c].sum()):,} "
                  f"({tab[c].sum() / grand[c] * 100:.1f}%) of {int(grand[c]):,} bins")

    cats = category_order or _order_categories(set(tab.index))
    cats = [c for c in cats if c in tab.index]
    tab = tab.loc[cats, columns]

    totals = tab.sum(axis=0)
    if (totals == 0).any():
        empty = list(totals[totals == 0].index)
        raise ValueError(f"no plottable categories left in {empty} after "
                         f"excluding {held}")
    pct = tab.div(totals, axis=1) * 100.0
    vals = pct if normalize == "percent" else tab

    lookup = dict(_COLOR_LOOKUP)
    if palette:
        lookup.update({_norm(k): v for k, v in palette.items()})
    colors, unknown = {}, []
    for i, c in enumerate(cats):
        if _norm(c) in lookup:
            colors[c] = lookup[_norm(c)]
        else:
            colors[c] = _FALLBACK[i % len(_FALLBACK)]
            unknown.append(c)
    if unknown:
        print(f"[WARN] no palette entry for {unknown} -- using fallback colours. "
              f"Check the category spelling against CANONICAL_ORDER.")

    n = len(columns)
    if figsize is None:
        figsize = ((max(4.2, 1.35 * n + 3.0), 5.6) if orientation == "vertical"
                   else (8.4, max(3.0, 0.72 * n + 1.8)))

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    pos = np.arange(n)
    base = np.zeros(n)
    for cat in cats:
        v = vals.loc[cat].to_numpy(dtype=float)
        if orientation == "vertical":
            ax.bar(pos, v, bottom=base, width=bar_width, color=colors[cat],
                   edgecolor="white", linewidth=0.7, zorder=3)
        else:
            ax.barh(pos, v, left=base, height=bar_width, color=colors[cat],
                    edgecolor="white", linewidth=0.7, zorder=3)

        if min_label is not None:
            share = pct.loc[cat].to_numpy(dtype=float)
            raw = tab.loc[cat].to_numpy(dtype=float)
            mode = normalize if label_as == "auto" else label_as
            # Gate on PHYSICAL size in count mode (share of the tallest bar) so
            # short bars are not crowded; on within-bar share in percent mode,
            # where every bar is the same height anyway.
            gate = (raw / float(totals.max()) * 100.0 if normalize == "count"
                    else share)
            for i in range(n):
                if gate[i] >= min_label:
                    cx, cy = pos[i], base[i] + v[i] / 2.0
                    if orientation == "horizontal":
                        cx, cy = cy, pos[i]
                    _c = colors[cat]
                    lum = (sum(int(str(_c).lstrip("#")[j:j+2], 16)
                               for j in (0, 2, 4)) / 3
                           if isinstance(_c, str) and _c.startswith("#") else 200)
                    txt_i = (f"{int(round(raw[i])):,}" if mode == "count"
                             else f"{share[i]:.0f}%")
                    ax.text(cx, cy, txt_i, ha="center", va="center",
                            fontsize=8, zorder=4,
                            color="white" if lum < 140 else "#222222")
        base = base + v

    scope = "changed bins" if held else "bins"
    unit = f"% of {scope}" if normalize == "percent" else f"number of {scope}"
    tick_labels = labels if labels is not None else list(columns)
    if orientation == "vertical":
        longest = max(len(str(t)) for t in tick_labels)
        # Auto-rotate: estimate how many characters fit under one bar and only
        # tilt when they would collide. The axes take roughly 60% of the figure
        # width once the legend is placed, and a character at this font size is
        # about 0.075 in wide. Horizontal labels stay horizontal while there is
        # room; past that they go to 45 degrees.
        avail_in = 0.60 * fig.get_size_inches()[0] / max(n, 1)
        rot = 0 if longest <= (avail_in / 0.075) else 45
        ax.set_xticks(pos)
        ax.set_xticklabels(tick_labels, rotation=rot,
                           ha="center" if rot == 0 else "right",
                           rotation_mode="anchor" if rot else None,
                           fontsize=9.5)
        ax.set_ylabel(unit, labelpad=6)
        if normalize == "percent":
            ax.set_ylim(0, 106 if show_n else 100); ax.set_yticks(range(0, 101, 20))
        ax.grid(axis="y", color="#DDDDDD", lw=0.6, zorder=0)
        ax.set_axisbelow(True)
    else:
        ax.set_yticks(pos); ax.set_yticklabels(tick_labels)
        ax.invert_yaxis()
        ax.set_xlabel(unit, labelpad=6)
        if normalize == "percent":
            ax.set_xlim(0, 112 if show_n else 100); ax.set_xticks(range(0, 101, 20))
        ax.grid(axis="x", color="#DDDDDD", lw=0.6, zorder=0)
        ax.set_axisbelow(True)

    if show_n:
        for i, c in enumerate(columns):
            txt = (f"n = {int(totals[c]):,} / {int(grand[c]):,}" if held
                   else f"{int(totals[c]):,}")
            if orientation == "vertical":
                ax.text(pos[i], (100.8 if normalize == "percent" else totals[c] * 1.01),
                        txt, ha="center", va="bottom", fontsize=8, color="#444444")
            else:
                ax.text((101 if normalize == "percent" else totals[c] * 1.01), pos[i],
                        txt, ha="left", va="center", fontsize=8, color="#444444")

    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    handles = [Patch(facecolor=colors[c], edgecolor="#B0B0B0", linewidth=0.5,
                     label=str(c).replace("_", " ")) for c in cats]
    ax.legend(handles=handles[::-1], frameon=False, fontsize=8.5,
              loc="center left", bbox_to_anchor=(1.02, 0.5), ncol=legend_ncol,
              handlelength=1.1, handleheight=1.1, labelspacing=0.45)

    if title:
        ax.set_title(title + ("\n" + subtitle if subtitle else ""),
                     fontweight="bold", fontsize=11.5,
                     pad=20 if (show_n and orientation == "vertical") else 12)

    fig.tight_layout()
    if save:
        fig.savefig(save, bbox_inches="tight", pad_inches=0.06, dpi=dpi,
                    facecolor="white")
    return fig, ax, tab, summary