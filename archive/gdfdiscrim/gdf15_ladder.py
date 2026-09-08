"""
GDF15 cutoff charts — all 322 patients against the two separation boundaries.

Renders the raw plasma GDF15 distribution for all 322 individuals in the source
cohort against the two empirical boundaries:

    lower = min(PVOD)        no PVOD patient falls below it
    upper = max(non-PVOD)    only PVOD patients fall above it

Both are raw order statistics of the observed data, not fitted thresholds. No
interval is drawn and none is implied.

Five chart styles, each rendered on a log and a linear axis, as PNG + PDF + SVG
into gdf15_charts/:

    Fig1   grouped   one swarm per diagnosis, names on the x axis
    Fig2   pooled    all patients in one swarm, symbols keyed by a legend
    Fig3   ranked    x = rank order, so the points trace a rising curve
    Fig4   scatter   x = random jitter, an unstructured cloud
    FigS1  labelled  every patient carries its diagnosis in text

Fig1 and Fig2 use a beeswarm and so take on a violin-like silhouette; Fig3 and
Fig4 deliberately avoid that shape.

By default everything is drawn at screen proportions. Pass --print to also emit
the journal-column sizes (183 mm and 89 mm).

Source data: Prabhakar et al., J Clin Invest 2026;136(2):e199159, Supplementary
Source Data, sheet `Fig_1a` (panel Fig. 1a), extracted to data/df_all.csv.

Usage:  venv/bin/python gdf15_ladder.py [--print]
"""

from __future__ import annotations

import csv
import shutil
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy import stats

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
OUT = ROOT / "gdf15_charts"

# ── expected values, asserted against the input so drift fails loudly ─────────
EXPECTED_N = 322
EXPECTED_COUNTS = {
    "Ctrl": 68, "PVOD": 26, "IPAH": 98, "HPAH": 42, "Other PH": 52, "COPD": 36,
}
EXPECTED_LOWER = 860.513
EXPECTED_UPPER = 1809.46

# Cohorts ordered by median GDF15, so every ordered display (legend, tables,
# lightness ramp) tells the same story.
COHORTS = ["Ctrl", "Other PH", "HPAH", "COPD", "IPAH", "PVOD"]

# ── palette ──────────────────────────────────────────────────────────────────
# Six slots taken from the dataviz reference palette (blue/aqua/yellow/magenta/
# green/violet — orange and red dropped). Validated as a set for this chart form:
#   node scripts/validate_palette.js "#2a78d6,#1baf7a,#eda100,#e87ba4,#008300,#4a3aa7" \
#        --mode light --surface "#ffffff" --pairs all
#   -> PASS lightness band, PASS chroma floor, PASS normal-vision floor (15.6),
#      WARN CVD 6.1 (6-8 band), WARN contrast (3 slots under 3:1)
# Both WARNs are discharged by mandatory secondary encoding, which this figure
# carries twice over: a distinct marker shape per cohort, and -- in FigS1 -- a
# visible direct label on every single point. data/gdf15_ladder_points.csv is
# the accompanying table view.
#
# Hues are assigned for grayscale separation, measured with the ITU-R 601 luma
# actually used by print conversion (0.299R + 0.587G + 0.114B on sRGB), not by
# eye and not by WCAG relative luminance -- those two disagree sharply here.
#
# Within these six hues exactly two pairs collapse in grayscale: violet/green
# (75 vs 77) and magenta/yellow (160 vs 165). They cannot be avoided while using
# all six, so they are *placed* where they cost nothing:
#
#   violet/green  -> PVOD and Ctrl. These two cohorts never occupy the same
#                    region of the axis at all (Ctrl max 684 < PVOD min 861,
#                    asserted below), so they can never be confused for one
#                    another in position.
#   magenta/yellow-> Other PH and HPAH, the two cohorts contributing fewest
#                    points to the ambiguous band (13 and 10), and which carry
#                    clearly different marker silhouettes.
#
# That frees real separation for the pairs the figure exists to distinguish:
# PVOD vs IPAH is Δ50 in grayscale, PVOD vs COPD Δ32, IPAH vs COPD Δ17.
COLOR = {
    "Ctrl":     "#008300",  # gray  77
    "Other PH": "#e87ba4",  # gray 160
    "HPAH":     "#eda100",  # gray 165
    "COPD":     "#2a78d6",  # gray 107
    "IPAH":     "#1baf7a",  # gray 125
    "PVOD":     "#4a3aa7",  # gray  75
}
# Redundant channel. Shapes are chosen to stay separable at 300 dpi small sizes,
# and to disambiguate the two grayscale-adjacent pairs (HPAH/Other PH,
# IPAH/COPD) by silhouette alone.
MARKER = {
    "Ctrl": "x", "Other PH": "D", "HPAH": "v", "COPD": "s", "IPAH": "^", "PVOD": "o",
}

# Project chrome, matching figures/02_-06_ (see plot_fagan, gdf15_lr.py:412).
INK, INK_SOFT, RULE, BAND_FILL = "#1a1a19", "#52514e", "#8a8a86", "#c9c9c6"


# ── data ─────────────────────────────────────────────────────────────────────
def load_points():
    """Read the verified Fig_1a extraction. Asserts n and per-cohort counts."""
    with open(DATA / "df_all.csv", newline="") as fh:
        rows = [
            {"gdf15": float(r["gdf15_raw"]),
             "cohort": r["cohort"].replace("\xa0", " ").strip(),
             "sex": r["sex"]}
            for r in csv.DictReader(fh)
        ]

    assert len(rows) == EXPECTED_N, f"expected {EXPECTED_N} rows, got {len(rows)}"
    counts = {c: sum(r["cohort"] == c for r in rows) for c in COHORTS}
    assert counts == EXPECTED_COUNTS, f"cohort counts drifted: {counts}"

    # The palette gives PVOD and Ctrl the one grayscale-degenerate hue pair, on
    # the grounds that their ranges are disjoint. If that ever stops being true
    # the colour assignment is no longer safe -- so check it, don't assume it.
    ctrl_max = max(r["gdf15"] for r in rows if r["cohort"] == "Ctrl")
    pvod_min = min(r["gdf15"] for r in rows if r["cohort"] == "PVOD")
    assert ctrl_max < pvod_min, (
        f"Ctrl and PVOD ranges now overlap ({ctrl_max:.1f} >= {pvod_min:.1f}); "
        "they can no longer share a grayscale-degenerate hue pair")
    return rows


def compute_lines(rows):
    """The two boundaries, derived from the data -- never hardcoded."""
    lower = min(r["gdf15"] for r in rows if r["cohort"] == "PVOD")
    upper = max(r["gdf15"] for r in rows if r["cohort"] != "PVOD")
    assert abs(lower - EXPECTED_LOWER) < 0.01, f"lower line drifted: {lower}"
    assert abs(upper - EXPECTED_UPPER) < 0.01, f"upper line drifted: {upper}"
    assert lower < upper
    return lower, upper


def zone_of(v, lower, upper):
    if v < lower:
        return "below"
    return "above" if v > upper else "band"


def zone_counts(rows, lower, upper):
    out = {z: {c: 0 for c in COHORTS} for z in ("below", "band", "above")}
    for r in rows:
        out[zone_of(r["gdf15"], lower, upper)][r["cohort"]] += 1
    return out


def cohort_stats(rows):
    """n / median / IQR per cohort, plus the tests quoted in the caption."""
    by = {c: np.array([r["gdf15"] for r in rows if r["cohort"] == c]) for c in COHORTS}
    desc = {
        c: {"n": len(v), "median": float(np.median(v)),
            "q1": float(np.percentile(v, 25)), "q3": float(np.percentile(v, 75))}
        for c, v in by.items()
    }

    H, p_kw = stats.kruskal(*by.values())

    others = [c for c in COHORTS if c != "PVOD"]
    raw = [stats.mannwhitneyu(by["PVOD"], by[c], alternative="two-sided")[1] for c in others]
    # Holm step-down
    order, running, adj = np.argsort(raw), 0.0, np.empty(len(raw))
    for i, j in enumerate(order):
        running = max(running, (len(raw) - i) * raw[j])
        adj[j] = min(running, 1.0)

    return {
        "desc": desc,
        "kruskal": {"H": float(H), "df": len(COHORTS) - 1, "p": float(p_kw)},
        "pairwise": {c: {"p": float(pr), "p_holm": float(pa)}
                     for c, pr, pa in zip(others, raw, adj)},
    }


# ── layout helpers ───────────────────────────────────────────────────────────
def assign_columns(y_frac, n_cols):
    """Deal points into columns round-robin by rank.

    Each column then holds every n_cols-th patient by value, so all columns span
    the full axis and carry equal density -- which is what buys each label
    n_cols times the vertical room and lets all 322 fit without dropping any.

    A point and its own label live in the *same* column. That is the property
    that matters: leader lines stay short and near-horizontal, and no leader
    ever crosses into a neighbouring column.
    """
    col = np.empty(len(y_frac), dtype=int)
    for rank, idx in enumerate(np.argsort(y_frac)):
        col[idx] = rank % n_cols
    return col


def declutter(ys, min_gap):
    """Spread values so no two are closer than min_gap, staying near the input.

    Forward pass pushes each label up off its predecessor; backward pass pulls
    the whole run back down off the ceiling. Order is preserved, so a label
    never crosses its neighbours and leader lines cannot tangle.
    """
    order = np.argsort(ys)
    v = np.array(ys, dtype=float)[order]

    for i in range(1, len(v)):
        v[i] = max(v[i], v[i - 1] + min_gap)
    for i in range(len(v) - 2, -1, -1):
        v[i] = min(v[i], v[i + 1] - min_gap)

    out = np.empty_like(v)
    out[order] = v
    return out


def beeswarm(y_frac, r_x, r_y):
    """Symmetric offsets so points touch but never overlap (clean panel)."""
    xs = np.zeros(len(y_frac))
    placed = []
    for idx in np.argsort(y_frac):
        y = y_frac[idx]
        near = [px for px, py in placed if abs(py - y) < 2 * r_y]
        step = 0
        while True:
            for cand in ((0.0,) if step == 0 else (step * 2 * r_x, -step * 2 * r_x)):
                if all(abs(cand - px) >= 2 * r_x * 0.98 for px in near):
                    xs[idx] = cand
                    placed.append((cand, y))
                    break
            else:
                step += 1
                continue
            break
    return xs


# ── style ────────────────────────────────────────────────────────────────────
def apply_style(base_pt):
    plt.rcParams.update({
        "font.family": "sans-serif",
        # Arial first, not Helvetica: macOS ships Helvetica only as a .dfont,
        # whose glyph metrics fail to load under Agg (and which embeds poorly in
        # PDF). Both are acceptable to Nature-family journals.
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": base_pt,
        "axes.labelsize": base_pt,
        "xtick.labelsize": base_pt, "ytick.labelsize": base_pt,
        "legend.fontsize": base_pt,
        # Real, selectable text in both vector formats -- production desks need
        # to restyle captions and axis text without re-rendering the figure.
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        # Exponents are drawn with mathtext (Arial has no Unicode
        # superscript glyphs); keep it in the body face, not DejaVu.
        "mathtext.fontset": "custom", "mathtext.default": "regular",
        "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "ytick.major.size": 2.5, "ytick.minor.size": 1.4,
        "ytick.minor.width": 0.4,
        "axes.edgecolor": INK_SOFT, "text.color": INK,
        "axes.labelcolor": INK, "ytick.color": INK_SOFT,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })


MM = 1 / 25.4


def frac_to_data(f, scale, lo, hi):
    """Axes-fraction -> data coordinate, for the active y scale."""
    if scale == "log":
        return 10 ** (np.log10(lo) + np.asarray(f) * (np.log10(hi) - np.log10(lo)))
    return lo + np.asarray(f) * (hi - lo)


def data_to_frac(v, scale, lo, hi):
    if scale == "log":
        return (np.log10(v) - np.log10(lo)) / (np.log10(hi) - np.log10(lo))
    return (np.asarray(v, dtype=float) - lo) / (hi - lo)


# ── the figure ───────────────────────────────────────────────────────────────
# Column count and canvas height are sized from the measured worst-case density
# of the data, not guessed: the tightest 10% window of the log axis holds 98 of
# the 322 points. At 183 mm, 12 columns x 247 mm clears that window with no
# crowding, so labels sit essentially at their true value. At 89 mm the label
# text cannot shrink below ~4.6 pt and only 7 columns fit across the page, so
# that window overflows by ~1.3x and labels there are displaced further from
# their point -- still complete, still non-overlapping, with longer leaders.
# 183 mm is therefore the version to submit; 89 mm is the fallback.
#
# "screen" is the default: ordinary landscape proportions for reading on a
# monitor, with type sized to be legible there rather than under a loupe. The
# two millimetre presets exist only for journal submission and are off by
# default -- pass --print to emit them.
#
# All heights are inches. The labelled form is the one that cannot be short:
# it needs enough vertical room for 322 lines of text, so it stays tall in
# every preset, and its linear variant taller still (a linear axis pushes 217
# of the 322 into the bottom 14%).
IN, MM_IN = 1.0, 1 / 25.4
LAYOUT = {
    "screen": dict(w=10.0, dpi=200, label_pt=7.0, base_pt=9.0, cols=16,
                   gutter=0.130, legend_ncol=6, zone_stack=False,
                   left_lab=0.085, left_other=0.105,
                   h_labelled={"log": 10.0, "linear": 13.0},
                   h_clean=6.0, h_pooled=6.8, h_ranked=6.0, h_jitter=6.4,
                   mk_lab=3.2, mk_clean=4.5, mk_pool=4.5, pool_spread=3.0,
                   mk_rank=4.0, mk_jit=4.2, jitter_w=0.78),
    "183mm": dict(w=183 * MM_IN, dpi=300, label_pt=5.8, base_pt=7.0, cols=12,
                  gutter=0.140, legend_ncol=6, zone_stack=False,
                  left_lab=0.105, left_other=0.135,
                  h_labelled={"log": 247 * MM_IN, "linear": 330 * MM_IN},
                  h_clean=132 * MM_IN, h_pooled=170 * MM_IN,
                  h_ranked=120 * MM_IN, h_jitter=150 * MM_IN,
                  mk_lab=2.4, mk_clean=3.0, mk_pool=3.2, pool_spread=3.0,
                  mk_rank=2.6, mk_jit=2.8, jitter_w=0.78),
    "89mm":  dict(w=89 * MM_IN, dpi=300, label_pt=4.6, base_pt=5.8, cols=7,
                  gutter=0.250, legend_ncol=3, zone_stack=True,
                  left_lab=0.145, left_other=0.205,
                  h_labelled={"log": 240 * MM_IN, "linear": 320 * MM_IN},
                  h_clean=112 * MM_IN, h_pooled=150 * MM_IN,
                  h_ranked=100 * MM_IN, h_jitter=130 * MM_IN,
                  mk_lab=1.8, mk_clean=2.1, mk_pool=2.4, pool_spread=2.0,
                  mk_rank=1.9, mk_jit=2.0, jitter_w=0.72),
}
PRINT_SIZES = ("183mm", "89mm")

# The panel forms this script can draw.
#   labelled  every patient carries its diagnosis in text (columns make room)
#   grouped   one swarm per diagnosis, cohort names on the x axis, no legend
#   pooled    all patients in one swarm, identified by symbol via a legend
#   ranked    x = rank order, so the points trace a rising curve
#   jitter    x = random scatter within a band, an unstructured cloud
#
# `grouped` and `pooled` both use a beeswarm, which stacks points symmetrically
# about a centre line and therefore takes on a violin-like silhouette. `ranked`
# and `jitter` exist to avoid that shape entirely: neither builds a symmetric
# profile, so neither can be mistaken for a density curve.
MODES = ("labelled", "grouped", "pooled", "ranked", "jitter")


def add_caption(fig, fig_w, fig_h, box, lead, body, label, cap_pt=7.0, margin=0.030):
    """Extend the page downward and typeset a journal caption beneath the axes.

    Call before adding the axes: the chart's own proportions are preserved by
    scaling its box into the taller page rather than by redesigning it.
    Returns the new (fig_h, box).
    """
    usable_pt = (1 - 2 * margin) * fig_w * 72
    wrap_at = int(usable_pt / (cap_pt * 0.52))   # Arial mean advance, measured
    head = textwrap.wrap(f"{label} | {lead}", wrap_at)
    lines = head + textwrap.wrap(body, wrap_at)

    lead_in = (cap_pt * 1.45) / 72
    cap_in = len(lines) * lead_in + 0.16
    new_h = fig_h + cap_in
    k = fig_h / new_h
    box = [box[0], (cap_in / new_h) + box[1] * k, box[2], box[3] * k]
    fig.set_size_inches(fig_w, new_h)

    y = (0.16 * 0.5 + len(lines) * lead_in) / new_h
    drawn = []
    for i, line in enumerate(lines):
        y -= lead_in / new_h
        drawn.append(fig.text(margin, y, line, ha="left", va="baseline",
                              fontsize=cap_pt, color=INK,
                              fontweight="bold" if i < len(head) else "normal"))

    # The wrap width is estimated from mean glyph advance; confirm against the
    # real renderer that no line actually overruns the text column.
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    widest = max(t.get_window_extent(r).width for t in drawn)
    limit = (1 - 2 * margin) * fig_w * fig.dpi
    assert widest <= limit, (
        f"caption line overruns the text column by {widest - limit:.1f}px")
    return new_h, box


def plot_figure(rows, stats_, lower, upper, *, mode, scale, size, bold=False,
                caption=None, stem_override=None, formats=None):
    assert mode in MODES, mode
    # `bold` keeps the plotted data identical and only enlarges the reading
    # furniture -- legend, zone summary, boundary captions -- for a version
    # meant to be seen from further away than a desk.
    emph = 1.35 if bold else 1.0
    labelled, grouped = mode == "labelled", mode == "grouped"

    cfg = LAYOUT[size]
    label_pt, n_cols = cfg["label_pt"], cfg["cols"]
    apply_style(cfg["base_pt"])

    fig_w = cfg["w"]
    fig_h = (cfg["h_labelled"][scale] if labelled else
             {"grouped": cfg["h_clean"], "pooled": cfg["h_pooled"],
              "ranked": cfg["h_ranked"], "jitter": cfg["h_jitter"]}[mode])
    if bold:
        fig_h *= 1.12          # room for the larger header block
    fig = plt.figure(figsize=(fig_w, fig_h))

    # Explicit axes box, so axes height in inches is known exactly and the
    # decluttering can be done in real typographic units.
    if labelled:
        left = cfg["left_lab"]
        box = [left, 0.035, 0.995 - left, 0.885]
    elif grouped:
        left = cfg["left_other"]
        box = [left, 0.130, 0.985 - left, 0.800]
    elif mode == "ranked":
        # Room at the bottom for the "patients, ranked" caption line.
        left = cfg["left_other"]
        box = [left, 0.105, 0.985 - left, 0.815]
    else:
        # No x-axis categories to make room for, so the plot can run lower.
        left = cfg["left_other"]
        box = [left, 0.055, 0.985 - left, 0.715 if bold else 0.875]
    if caption is not None:
        fig_h, box = add_caption(fig, fig_w, fig_h, box, *caption, "Figure 4")

    ax = fig.add_axes(box)
    ax_h_in = box[3] * fig_h
    ax_w_in = box[2] * fig_w

    vals = np.array([r["gdf15"] for r in rows])
    cohorts = [r["cohort"] for r in rows]

    # Axis limits: pad in the active space so nothing sits on a spine.
    if scale == "log":
        lo, hi = 10 ** (np.log10(vals.min()) - 0.09), 10 ** (np.log10(vals.max()) + 0.09)
        ax.set_yscale("log")
    else:
        pad = 0.045 * (vals.max() - vals.min())
        lo, hi = max(0.0, vals.min() - pad), vals.max() + pad
    ax.set_ylim(lo, hi)
    ax.set_xlim(0, 1)

    y_frac = data_to_frac(vals, scale, lo, hi)
    blend = ax.get_yaxis_transform()  # x in axes fraction, y in data

    # ── shaded ambiguity band + the two boundaries ──
    ax.axhspan(lower, upper, color=BAND_FILL, alpha=0.30, lw=0, zorder=0)
    for y in (lower, upper):
        ax.axhline(y, color=INK, lw=1.5 if bold else 0.8, zorder=4)

    # ── points ──
    gutter = cfg["gutter"]          # reserved strip carrying the line annotations
    if labelled:
        mk = cfg["mk_lab"]
        col = assign_columns(y_frac, n_cols)
        pitch = (1.0 - gutter) / n_cols
        x_pt = gutter + pitch * (col + 0.06)
    elif grouped:
        # One beeswarm per cohort across the width. A single 1-D swarm cannot
        # fill a 183 mm column; grouping by diagnosis uses the space to say
        # something -- each cohort's whole distribution read against the two
        # boundaries -- and lets the cohort names act as direct labels, so the
        # main-text panel needs no legend box at all.
        mk = cfg["mk_clean"]
        r_x = ((mk / 2 + 0.35) / 72) / ax_w_in
        r_y = ((mk / 2 + 0.35) / 72) / ax_h_in
        slot = (1.0 - gutter) / len(COHORTS)
        x_pt = np.zeros(len(vals))
        for k, c in enumerate(COHORTS):
            m = np.array([i for i, cc in enumerate(cohorts) if cc == c])
            x_pt[m] = gutter + slot * (k + 0.5) + beeswarm(y_frac[m], r_x, r_y)
    elif mode == "pooled":
        # All 322 in one swarm, cohort carried by symbol alone. Every patient is
        # read against every other on a single axis -- which is the comparison
        # the two boundary lines are actually about. The horizontal spread is
        # widened past the marker radius so the swarm uses the column width
        # rather than collapsing into a thin spine.
        mk = cfg["mk_pool"]
        r_y = ((mk / 2 + 0.35) / 72) / ax_h_in
        r_x = ((mk / 2 + 0.35) * cfg["pool_spread"] / 72) / ax_w_in
        centre = gutter + (1.0 - gutter) / 2
        x_pt = centre + beeswarm(y_frac, r_x, r_y)
    elif mode == "ranked":
        # x is simply each patient's rank, lowest GDF15 to highest, so the
        # points trace one rising curve. Nothing is stacked and nothing is
        # mirrored, so there is no density profile to mistake for a violin --
        # and where the two boundaries cut the curve is directly readable.
        mk = cfg["mk_rank"]
        order = np.argsort(vals)
        rank = np.empty(len(vals))
        rank[order] = np.arange(len(vals))
        pad = 0.012
        x_pt = gutter + pad + (rank / (len(vals) - 1)) * (1 - gutter - 2 * pad)
    else:
        # x is random within a fixed band -- an unstructured cloud. Seeded, so
        # the figure is reproducible run to run.
        mk = cfg["mk_jit"]
        rng = np.random.default_rng(20260821)
        w = cfg["jitter_w"] * (1.0 - gutter)
        centre = gutter + (1.0 - gutter) / 2
        x_pt = centre + rng.uniform(-w / 2, w / 2, len(vals))

    for c in COHORTS:
        m = [i for i, cc in enumerate(cohorts) if cc == c]
        ax.plot(x_pt[m], vals[m], linestyle="none", marker=MARKER[c],
                markersize=mk, color=COLOR[c],
                markeredgecolor=COLOR[c] if c != "Ctrl" else COLOR[c],
                markeredgewidth=0.45 if c == "Ctrl" else 0.3,
                markerfacecolor="none" if c == "Ctrl" else COLOR[c],
                transform=blend, zorder=5, clip_on=False)

    # ── every point labelled ──
    n_labels = 0
    if labelled:
        min_gap = ((label_pt * 1.28) / 72) / ax_h_in   # leading -> axes fraction
        text_dx = ((mk / 2 + 1.6) / 72) / ax_w_in    # marker edge -> text start

        for k in range(n_cols):
            idx = np.where(col == k)[0]
            if not len(idx):
                continue
            placed = np.clip(declutter(y_frac[idx], min_gap), 0.003, 0.997)

            for i, yf in zip(idx, placed):
                c = cohorts[i]
                anchor = x_pt[i] + text_dx
                # Leader: point -> its own label, always within one column.
                ax.plot([x_pt[i], anchor - text_dx * 0.25],
                        frac_to_data([y_frac[i], yf], scale, lo, hi),
                        color=COLOR[c], lw=0.3, alpha=0.7,
                        transform=blend, zorder=2, clip_on=False,
                        solid_capstyle="round")
                ax.text(anchor, frac_to_data(yf, scale, lo, hi), c,
                        transform=blend, fontsize=label_pt, color=INK,
                        ha="left", va="center", zorder=6, clip_on=False)
                n_labels += 1

        # The one thing the brief marked non-negotiable: no label may be lost.
        assert n_labels == EXPECTED_N, f"expected {EXPECTED_N} labels, drew {n_labels}"
        label_artists = [t for t in ax.texts if t.get_text() in COHORTS]

    # ── boundary annotations, in the reserved left gutter ──
    ann = cfg["base_pt"] * emph
    for y, above, txt in ((upper, True, f"{upper:,.0f} pg/mL\nonly PVOD above"),
                          (lower, False, f"{lower:,.0f} pg/mL\nno PVOD below")):
        if scale == "log":
            y_txt = y * (1.012 if above else 0.988)
        else:
            y_txt = y + (0.004 if above else -0.004) * (hi - lo)
        ax.text(0.005, y_txt, txt, transform=blend, ha="left",
                va="bottom" if above else "top", linespacing=1.3,
                fontsize=ann, color=INK, zorder=7,
                fontweight="bold" if bold else "normal",
                bbox=dict(boxstyle="square,pad=0.32", fc="white",
                          ec="none", alpha=0.85) if bold else None)

    # ── zone counts, one compact line under the legend ──
    zc = zone_counts(rows, lower, upper)
    n_below = sum(zc["below"].values())
    n_band = sum(zc["band"].values())
    n_above = sum(zc["above"].values())
    parts = [f"PVOD excluded below {lower:,.0f}  ·  n = {n_below}",
             f"Ambiguous {lower:,.0f}–{upper:,.0f}  ·  n = {n_band} "
             f"(PVOD {zc['band']['PVOD']})",
             f"PVOD only above {upper:,.0f}  ·  n = {n_above}"]
    # One line is far too wide for a single column; stack it there instead.
    if bold:
        # Three separate chips read as three zones far better than one grey
        # run-on line does, and they sit directly above the bands they name.
        for frac, txt in zip((0.0, 0.355, 0.715), parts):
            ax.text(frac, 1.012, txt, transform=ax.transAxes,
                    ha="left", va="bottom", fontsize=cfg["base_pt"] * 1.2,
                    color=INK, fontweight="bold", linespacing=1.35,
                    bbox=dict(boxstyle="round,pad=0.42", fc="#f2f1ee",
                              ec=BAND_FILL, lw=0.8))
    else:
        joiner = "\n" if cfg["zone_stack"] else "          "
        ax.text(0, 1.006, joiner.join(parts), transform=ax.transAxes,
                ha="left", va="bottom", linespacing=1.45,
                fontsize=ann, color=INK_SOFT)

    # ── cohort medians, drawn as a short rule across each swarm ──
    if grouped:
        slot = (1.0 - gutter) / len(COHORTS)
        for k, c in enumerate(COHORTS):
            cx = gutter + slot * (k + 0.5)
            med = stats_["desc"][c]["median"]
            ax.plot([cx - slot * 0.30, cx + slot * 0.30], [med, med],
                    color=INK, lw=0.9, transform=blend, zorder=6,
                    solid_capstyle="butt")

    # ── axis furniture ──
    ax.set_ylabel("Plasma GDF15 (pg/mL)", fontsize=cfg["base_pt"] * emph, labelpad=3)
    if mode == "ranked":
        # Ticks off but the axis itself left visible -- hiding the whole axis
        # would take the caption with it.
        ax.set_xticks([])
        ax.set_xlabel(f"{len(rows)} patients, ordered by GDF15 (lowest to highest)",
                      fontsize=cfg["base_pt"], color=INK_SOFT, labelpad=6)
    elif not grouped:
        ax.xaxis.set_visible(False)
    else:
        slot = (1.0 - gutter) / len(COHORTS)
        ax.set_xticks([gutter + slot * (k + 0.5) for k in range(len(COHORTS))])
        ax.set_xticklabels([f"{c}\nn = {stats_['desc'][c]['n']}" for c in COHORTS],
                           fontsize=cfg["base_pt"], linespacing=1.4)
        ax.tick_params(axis="x", length=0, pad=3, colors=INK)
    for s in ("top", "right", "bottom"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color(RULE)

    if scale == "log":
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
            lambda v, _: f"{v:,.0f}"))
        ax.set_yticks([25, 50, 100, 250, 500, 1000, 2500, 5000])
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    else:
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
            lambda v, _: f"{v:,.0f}"))

    # ── legend, with n per cohort ──
    # The grouped panel is the one that skips it: there the cohort names sit
    # under each swarm, which is a stronger direct label than a legend box. The
    # pooled panel has no other way to name a symbol, so the legend is the whole
    # key to it.
    handles = [] if grouped else [
        Line2D([], [], linestyle="none", marker=MARKER[c],
               markersize=max(mk, 3.2) * (1.5 if bold else 1.0), color=COLOR[c],
               markerfacecolor="none" if c == "Ctrl" else COLOR[c],
               markeredgecolor=COLOR[c], markeredgewidth=0.45,
               label=f"{c} (n = {stats_['desc'][c]['n']})")
        for c in reversed(COHORTS)
    ]
    if handles:
        ax.legend(handles=handles, loc="lower left",
                  bbox_to_anchor=(0, 1.10 if bold else 1.032),
                  ncol=max(2, cfg["legend_ncol"] // 2) if bold
                       else cfg["legend_ncol"],
                  frameon=bold,
                  edgecolor=BAND_FILL, facecolor="white", framealpha=1.0,
                  handletextpad=0.5 if bold else 0.35,
                  columnspacing=1.5 if bold else 1.1,
                  borderpad=0.6 if bold else 0,
                  labelspacing=0.35,
                  fontsize=cfg["base_pt"] * emph)

    # ── geometric gate: labels must not collide, and must stay on the canvas ──
    # Cheaper and far more reliable than eyeballing 322 labels per variant.
    worst_overlap = 0.0
    if labelled:
        fig.canvas.draw()
        boxes = [t.get_window_extent(fig.canvas.get_renderer()) for t in label_artists]
        fw, fh = fig.get_size_inches() * fig.dpi

        off = [b for b in boxes if b.x0 < -0.5 or b.x1 > fw + 0.5
               or b.y0 < -0.5 or b.y1 > fh + 0.5]
        assert not off, f"{len(off)} label(s) fall outside the canvas"

        # A label wider than the column pitch would collide with the next
        # column, which the per-column scan below cannot see.
        pitch_px = ((1.0 - gutter) / n_cols) * ax_w_in * fig.dpi
        widest = max(b.x1 - b.x0 for b in boxes)
        assert widest + text_dx * ax_w_in * fig.dpi < pitch_px, (
            f"widest label {widest:.1f}px + offset exceeds column pitch "
            f"{pitch_px:.1f}px -- labels would collide across columns")

        # Only same-column labels can collide, and within a column they are
        # y-ordered, so comparing each to the next by y is sufficient.
        for k in range(n_cols):
            bs = sorted((boxes[i] for i in np.where(col == k)[0]), key=lambda b: b.y0)
            for a, b in zip(bs, bs[1:]):
                gap = b.y0 - a.y1
                if gap < 0 and (min(a.x1, b.x1) - max(a.x0, b.x0)) > 0:
                    worst_overlap = max(worst_overlap, -gap)
        assert worst_overlap < 0.5, f"labels overlap by up to {worst_overlap:.2f} px"

    # ── write ──
    stem = {"labelled": "FigS1_gdf15_cutoffs_labelled",
            "grouped": "Fig1_gdf15_cutoffs",
            "pooled": "Fig2_gdf15_cutoffs_symbols",
            "ranked": "Fig3_gdf15_cutoffs_ranked",
            "jitter": "Fig4_gdf15_cutoffs_scatter"}[mode]
    stem = (f"{stem}_{scale}" + ("" if size == "screen" else f"_{size}")
            + ("_bold" if bold else ""))
    dest = OUT
    if stem_override is not None:
        dest, stem = stem_override
    dest.mkdir(parents=True, exist_ok=True)
    for ext, kw in (formats or (("pdf", {}), ("svg", {}), ("png", {"dpi": cfg["dpi"]}))):
        fig.savefig(dest / f"{stem}.{ext}", format=ext, **kw)
    plt.close(fig)
    return stem, n_labels


# ── captions ─────────────────────────────────────────────────────────────────
SOURCE = ("Prabhakar et al., *J Clin Invest* 2026;136(2):e199159, Supplementary "
          "Source Data, sheet `Fig_1a` (panel Fig. 1a).")


_SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def _sci(p):
    e = int(np.floor(np.log10(p)))
    return f"{p / 10 ** e:.1f} × 10{str(e).translate(_SUP)}"


def _sci_tex(p):
    """Same value as _sci, as one unbreakable mathtext token."""
    e = int(np.floor(np.log10(p)))
    return f"${p / 10 ** e:.1f}{{\\times}}10^{{{e}}}$"


def caption_parts(rows, stats_, lower, upper, math=False):
    """The Figure 4 caption as (bold lead, body), in plain text.

    Returned rather than inlined so the rendered proof and CAPTIONS.md are
    generated from one string and cannot drift apart.
    """
    d, zc = stats_["desc"], zone_counts(rows, lower, upper)
    kw = stats_["kruskal"]
    ns = ", ".join(f"{'control' if c == 'Ctrl' else c} (n = {d[c]['n']})"
                   for c in COHORTS)
    band = ", ".join(f"{c} {zc['band'][c]}" for c in reversed(COHORTS) if zc["band"][c])
    worst = max(v["p_holm"] for v in stats_["pairwise"].values())
    n_band = sum(zc["band"].values())
    sci = _sci_tex if math else _sci

    lead = ("Plasma GDF15 in all 322 individuals, shown as an unstructured "
            "scatter against the two separation boundaries.")
    body = (
        f"Raw plasma GDF15 (pg/mL) for {len(rows)} individuals: {ns}. Each point "
        "is one patient. Horizontal position is random jitter applied solely to "
        "separate overlapping points; it carries no information about the "
        "patient, and the jitter is generated from a fixed seed so the figure is "
        "reproducible. No density estimate, kernel or distributional envelope is "
        "fitted or drawn — the horizontal extent of the cloud is uniform by "
        "construction and must not be read as a distribution shape. Cohort is "
        "encoded redundantly by both colour and marker shape, and hues were "
        "assigned so that the cohorts the figure exists to separate stay distinct "
        "after grayscale conversion and under simulated protan/deutan "
        "colour-vision deficiency. Horizontal lines mark the two empirical "
        f"separation boundaries: {lower:,.1f} pg/mL, below which no PVOD patient "
        f"falls, and {upper:,.1f} pg/mL, above which only PVOD patients fall. Both "
        "lines are raw order statistics of these observed data — respectively the "
        "minimum of the PVOD cohort and the maximum of the pooled non-PVOD "
        "cohorts — not fitted or optimised thresholds; no confidence interval is "
        "shown and none is implied. The shaded band between them contains "
        f"{n_band} patients ({band}) whose GDF15 value alone does not resolve "
        "their diagnosis. GDF15 differs across cohorts (Kruskal–Wallis "
        f"H = {kw['H']:.1f}, df = {kw['df']}, P = {sci(kw['p'])}); PVOD exceeds "
        "every other cohort pairwise (Mann–Whitney U, Holm-adjusted "
        f"P ≤ {sci(worst)} for all {len(stats_['pairwise'])} comparisons). "
        "Cohort medians (IQR), pg/mL: "
        + "; ".join(f"{c} {d[c]['median']:,.1f} ({d[c]['q1']:,.1f}–{d[c]['q3']:,.1f})"
                    for c in COHORTS)
        + ". y axis is linear. Source data: Prabhakar et al., J Clin Invest "
        "2026;136(2):e199159, Supplementary Source Data, sheet Fig_1a "
        "(panel Fig. 1a).")
    return lead, body


def write_captions(rows, stats_, lower, upper):
    d, zc = stats_["desc"], zone_counts(rows, lower, upper)
    kw = stats_["kruskal"]
    ns = ", ".join(f"{'control' if c == 'Ctrl' else c} (n = {d[c]['n']})"
                   for c in COHORTS)
    band = ", ".join(f"{c} {zc['band'][c]}" for c in reversed(COHORTS) if zc['band'][c])
    worst = max(v["p_holm"] for v in stats_["pairwise"].values())
    n_band = sum(zc["band"].values())

    lines = [
        "# Draft figure captions",
        "",
        f"Source data for every figure: {SOURCE}",
        "",
        "Generated by `gdf15_ladder.py` from the computed values — do not retype "
        "the numbers by hand, or they will drift from the figures.",
        "",
        "---",
        "",
        "## Figure 1",
        "",
        "**Plasma GDF15 separates PVOD from other pulmonary hypertension only at "
        "the extremes.**",
        "",
        f"Raw plasma GDF15 (pg/mL) for {len(rows)} individuals: {ns}. Each point "
        "is one patient, plotted as a within-cohort swarm; the black rule across "
        "each swarm is that cohort's median. Cohort is encoded redundantly by "
        "both colour and marker shape, and hues were assigned so that the "
        "cohorts the figure exists to separate stay distinct after grayscale "
        "conversion and under simulated protan/deutan colour-vision deficiency. "
        "Horizontal lines mark the two empirical separation "
        f"boundaries: {lower:,.1f} pg/mL, below which no PVOD patient falls, and "
        f"{upper:,.1f} pg/mL, above which only PVOD patients fall. **Both lines "
        "are raw order statistics of these observed data — respectively the "
        "minimum of the PVOD cohort and the maximum of the pooled non-PVOD "
        "cohorts — not fitted or optimised thresholds; no confidence interval is "
        "shown and none is implied.** The shaded band between them contains "
        f"{n_band} patients ({band}) whose GDF15 value alone does not resolve "
        "their diagnosis. GDF15 differs across cohorts (Kruskal–Wallis "
        f"H = {kw['H']:.1f}, df = {kw['df']}, P = {_sci(kw['p'])}); PVOD exceeds "
        "every other cohort pairwise (Mann–Whitney U, Holm-adjusted "
        f"P ≤ {_sci(worst)} for all {len(stats_['pairwise'])} comparisons). "
        "Cohort medians (IQR), pg/mL: "
        + "; ".join(f"{c} {d[c]['median']:,.1f} ({d[c]['q1']:,.1f}–{d[c]['q3']:,.1f})"
                    for c in COHORTS)
        + ". *y* axis is log₁₀-scaled (linear variant: *y* axis is linear). "
        f"Source data: {SOURCE}",
        "",
        "---",
        "",
        "## Supplementary Figure S1",
        "",
        "**As Figure 1, with every patient individually labelled by diagnosis.**",
        "",
        f"Identical data and boundaries to Figure 1 (n = {len(rows)}; group n as "
        f"above), with all {len(rows)} patients labelled individually. Patients "
        "are dealt by rank into evenly spaced columns purely to make room for the "
        "labels; horizontal position carries no meaning, and every point is drawn "
        "at its true value on the shared vertical axis. Where a label has been "
        "displaced vertically to avoid overlapping its neighbour, a leader line "
        "connects it to its point. Lines as in Figure 1 — raw order statistics of "
        "the observed data, no interval implied. Note the single IPAH patient at "
        f"{upper:,.1f} pg/mL lying exactly on the upper boundary: that one patient "
        "is what prevents a cleaner threshold. Per-point values underlying every "
        "figure are tabulated in `gdf15_ladder_points.csv`. "
        f"Source data: {SOURCE}",
        "",
    ]

    # Shared closing sentences, so every caption states the boundaries and the
    # tests in exactly the same words.
    boundaries = (
        f"Horizontal lines mark the two empirical separation boundaries: "
        f"{lower:,.1f} pg/mL, below which no PVOD patient falls, and "
        f"{upper:,.1f} pg/mL, above which only PVOD patients fall. **Both lines "
        "are raw order statistics of these observed data — respectively the "
        "minimum of the PVOD cohort and the maximum of the pooled non-PVOD "
        "cohorts — not fitted or optimised thresholds; no confidence interval is "
        "shown and none is implied.** The shaded band between them contains "
        f"{n_band} patients ({band}) whose GDF15 value alone does not resolve "
        "their diagnosis.")
    tests = (
        "GDF15 differs across cohorts (Kruskal–Wallis "
        f"H = {kw['H']:.1f}, df = {kw['df']}, P = {_sci(kw['p'])}); PVOD exceeds "
        "every other cohort pairwise (Mann–Whitney U, Holm-adjusted "
        f"P ≤ {_sci(worst)} for all {len(stats_['pairwise'])} comparisons).")
    encoding = (
        "Cohort is encoded redundantly by both colour and marker shape, and hues "
        "were assigned so that the cohorts the figure exists to separate stay "
        "distinct after grayscale conversion and under simulated protan/deutan "
        "colour-vision deficiency.")

    lines += [
        "---",
        "",
        "## Figure 2",
        "",
        "**Plasma GDF15 in all 322 individuals against the two separation "
        "boundaries.**",
        "",
        f"Raw plasma GDF15 (pg/mL) for {len(rows)} individuals: {ns}. Each point "
        "is one patient; all cohorts are pooled into a single swarm so that every "
        "patient is read against every other on one axis. Horizontal position "
        "reflects only the packing needed to keep points from overlapping and "
        f"carries no information about the patient. {encoding} {boundaries} "
        f"{tests} *y* axis is log₁₀-scaled (linear variant: *y* axis is linear). "
        f"Source data: {SOURCE}",
        "",
        "---",
        "",
        "## Figure 3",
        "",
        "**Plasma GDF15 across the cohort in rank order, showing where the two "
        "boundaries fall.**",
        "",
        f"Raw plasma GDF15 (pg/mL) for {len(rows)} individuals: {ns}. Each point "
        "is one patient, ordered along the *x* axis by ascending GDF15; the *x* "
        "axis is therefore rank, not a measured quantity, and the points trace "
        f"the empirical cumulative distribution of the cohort. {encoding} "
        f"{boundaries} {tests} *y* axis is log₁₀-scaled (linear variant: *y* axis "
        f"is linear). Source data: {SOURCE}",
        "",
        "---",
        "",
        "## Figure 4",
        "",
        "**Plasma GDF15 in all 322 individuals, shown as an unstructured "
        "scatter against the two separation boundaries.**",
        "",
        f"Raw plasma GDF15 (pg/mL) for {len(rows)} individuals: {ns}. Each point "
        "is one patient. Horizontal position is random jitter applied solely to "
        "separate overlapping points; it carries no information about the "
        "patient, and the jitter is generated from a fixed seed so the figure is "
        "reproducible. No density estimate, kernel or distributional envelope is "
        "fitted or drawn — the horizontal extent of the cloud is uniform by "
        f"construction and must not be read as a distribution shape. {encoding} "
        f"{boundaries} {tests} Cohort medians (IQR), pg/mL: "
        + "; ".join(f"{c} {d[c]['median']:,.1f} ({d[c]['q1']:,.1f}–{d[c]['q3']:,.1f})"
                    for c in COHORTS)
        + ". *y* axis is linear (log variant: *y* axis is log₁₀-scaled). "
        f"Per-point values are tabulated in `gdf15_ladder_points.csv`. "
        f"Source data: {SOURCE}",
        "",
        "*A second copy of the linear panel, "
        "`Fig4_gdf15_cutoffs_scatter_linear_bold`, is supplied for presentation "
        "use. It plots identical data and identical boundaries; only the legend, "
        "zone labels and boundary captions are enlarged. The caption above "
        "applies to it unchanged.*",
        "",
    ]
    (OUT / "CAPTIONS.md").write_text("\n".join(lines))


def write_readme(rows, stats_, lower, upper):
    """A plain-language guide to the folder, generated so the numbers stay true."""
    zc = zone_counts(rows, lower, upper)
    n_below, n_band, n_above = (sum(zc[z].values()) for z in ("below", "band", "above"))
    band = ", ".join(f"{c} {zc['band'][c]}" for c in reversed(COHORTS) if zc["band"][c])

    lines = [
        "# GDF15 cutoff charts",
        "",
        f"All {len(rows)} patients plotted by their plasma GDF15 level, with two",
        "horizontal lines marking where the diagnosis becomes certain.",
        "",
        "## The two lines",
        "",
        f"- **{lower:,.0f} pg/mL** — the lowest PVOD patient. Nobody below this line is",
        f"  PVOD. That clears **{n_below} patients**, about "
        f"{round(100 * n_below / len(rows))}% of everyone.",
        f"- **{upper:,.0f} pg/mL** — the highest patient who is *not* PVOD (an IPAH case).",
        f"  Everyone above this line is PVOD: **{n_above} patients**, which is "
        f"{round(100 * n_above / stats_['desc']['PVOD']['n'])}% of all PVOD cases.",
        f"- **In between** sit **{n_band} patients** ({band}) where GDF15 alone does not",
        "  settle the diagnosis. That grey band is the real finding.",
        "",
        "Both lines are just the highest and lowest actual patients — not thresholds",
        "fitted by a model. Nothing here implies they are precise.",
        "",
        "## The five charts",
        "",
        "| File | What it shows | Points labelled | Violin-like shape |",
        "|---|---|---|---|",
        "| `Fig1_gdf15_cutoffs` | one group per diagnosis, with medians | no | yes |",
        "| `Fig2_gdf15_cutoffs_symbols` | everyone pooled, symbols + legend | no | yes |",
        "| `Fig3_gdf15_cutoffs_ranked` | everyone in order, lowest to highest | no | no |",
        "| `Fig4_gdf15_cutoffs_scatter` | everyone scattered at random | no | no |",
        "| `FigS1_gdf15_cutoffs_labelled` | every patient named individually | yes | no |",
        "",
        "Each comes in `_log` and `_linear`, as PNG, PDF and SVG.",
        "",
        "`Fig4_gdf15_cutoffs_scatter_linear_bold` is a second copy of the linear",
        "scatter with the legend, the zone labels and the line captions enlarged,",
        "for showing on a screen or a slide. The points and the lines are identical",
        "to the plain version — only the labelling is louder.",
        "",
        "- **log** spreads out the crowded 200–900 region so the grey band is readable.",
        "- **linear** is true spacing, but squashes most patients into the bottom strip.",
        "",
        "On every chart the **only real measurement is vertical** — height is the GDF15",
        "level. Horizontal position is layout only, except in `Fig1` (where it is the",
        "diagnosis group) and `Fig3` (where it is rank order).",
        "",
        "## Also here",
        "",
        "- `CAPTIONS.md` — formal figure captions with group sizes and statistics.",
        f"- `gdf15_ladder_points.csv` — the {len(rows)} plotted values, so any dot can be",
        "  traced back to its number.",
        "",
        "## Regenerating",
        "",
        "```",
        "venv/bin/python gdf15_ladder.py            # these screen-sized charts",
        "venv/bin/python gdf15_ladder.py --print    # adds 183 mm / 89 mm journal sizes",
        "```",
        "",
        f"Source data: {SOURCE}",
        "",
    ]
    OUT.mkdir(exist_ok=True)
    (OUT / "README.md").write_text("\n".join(lines))


def write_points_table(rows, lower, upper):
    OUT.mkdir(exist_ok=True)
    with open(OUT / "gdf15_ladder_points.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["gdf15_raw", "cohort", "sex", "zone"])
        for r in sorted(rows, key=lambda r: -r["gdf15"]):
            w.writerow([f"{r['gdf15']:.3f}", r["cohort"], r["sex"],
                        zone_of(r["gdf15"], lower, upper)])


SUBMIT = ROOT / "figure4_submission"


def build_submission(rows, stats_, lower, upper):
    """Assemble the single-figure submission folder.

    Two existing artefacts are reorganised into it -- the chart itself and its
    caption -- and one new file is produced: the two combined into a single
    journal-format page, 183 mm wide (Nature double-column), vector PDF with
    embedded Type 42 fonts, caption set at 7 pt sans-serif beneath the panel.
    """
    lead, body = caption_parts(rows, stats_, lower, upper)
    lead_m, body_m = caption_parts(rows, stats_, lower, upper, math=True)
    SUBMIT.mkdir(exist_ok=True)

    # (1) reorganised: the chart, exactly as it already exists
    src = OUT / "Fig4_gdf15_cutoffs_scatter_linear.svg"
    shutil.copy2(src, SUBMIT / src.name)

    # (2) reorganised: this figure's caption, lifted out of the combined file
    (SUBMIT / "Figure4_caption.md").write_text(
        f"**Figure 4 | {lead}**\n\n{body}\n")

    # (3) new: chart and caption composed as one submission-ready page
    stem, _ = plot_figure(rows, stats_, lower, upper, mode="jitter",
                          scale="linear", size="183mm",
                          caption=(lead_m, body_m),
                          stem_override=(SUBMIT, "Figure4_submission_proof"),
                          # Vector PDF only: one file, and the format
                          # journals actually want for line art.
                          formats=(("pdf", {}),))
    return stem


# ── main ─────────────────────────────────────────────────────────────────────
def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    # Journal-page sizes are opt-in: they are tall and narrow because they are
    # cut to a print column, which is the wrong shape for reading on a screen.
    sizes = ["screen"] + (list(PRINT_SIZES) if "--print" in argv else [])

    rows = load_points()
    lower, upper = compute_lines(rows)
    st = cohort_stats(rows)
    zc = zone_counts(rows, lower, upper)

    print(f"n = {len(rows)}    lower = {lower:,.1f}    upper = {upper:,.1f}\n")
    print(f"{'cohort':<10}{'n':>5}{'median':>10}{'IQR':>22}")
    for c in COHORTS:
        d = st["desc"][c]
        iqr = f"{d['q1']:,.1f} – {d['q3']:,.1f}"
        print(f"{c:<10}{d['n']:>5}{d['median']:>10,.1f}{iqr:>22}")
    print()
    for z, name in (("below", "below lower"), ("band", "between lines"), ("above", "above upper")):
        tot = sum(zc[z].values())
        detail = ", ".join(f"{c} {zc[z][c]}" for c in reversed(COHORTS) if zc[z][c])
        print(f"{name:<15} n = {tot:>3}   {detail}")
    print()

    write_points_table(rows, lower, upper)
    write_captions(rows, st, lower, upper)
    write_readme(rows, st, lower, upper)

    for size in sizes:
        for mode in MODES:
            for scale in ("log", "linear"):
                stem, n = plot_figure(rows, st, lower, upper,
                                      mode=mode, scale=scale, size=size)
                print(f"wrote {stem}.[pdf|svg|png]" + (f"   labels = {n}" if n else ""))
    # A louder copy of the linear scatter -- same data and same points as
    # Fig4_..._scatter_linear, only the reading furniture is enlarged.
    stem, _ = plot_figure(rows, st, lower, upper, mode="jitter",
                          scale="linear", size="screen", bold=True)
    print(f"wrote {stem}.[pdf|svg|png]")

    build_submission(rows, st, lower, upper)
    print(f"\nsubmission folder -> {SUBMIT}")

    if len(sizes) == 1:
        print("\n(journal 183 mm / 89 mm versions not written; pass --print for those)")

    print(f"\nreadme   -> {OUT / 'README.md'}")
    print(f"captions -> {OUT / 'CAPTIONS.md'}")
    print(f"table    -> {OUT / 'gdf15_ladder_points.csv'}")


if __name__ == "__main__":
    main()
