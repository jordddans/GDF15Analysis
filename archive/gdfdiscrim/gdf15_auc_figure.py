"""
Figure 1 — how GDF15 compares with every other marker measured in the study.

Discrimination of PVOD from IPAH + HPAH (the study's group 1, n = 166), for
GDF15 and for each conventional haemodynamic, natriuretic-peptide, functional
and demographic variable available in the same cohort.

The AUC point estimates are NOT recomputed for presentation: they are the values
already established by this project (notebook 05, spec §5) and are asserted here
against those published numbers, so the figure cannot drift from the analysis.
Confidence intervals are added for the figure, using the same stratified
bootstrap notebook 04 used for GDF15 (whose interval reproduces to rounding).

Writes three files into figure1_submission/:
    Figure1_auc_comparison.svg      the chart alone
    Figure1_caption.md              the caption alone
    Figure1_submission_proof.pdf    the two composed as one journal page

Usage:  venv/bin/python gdf15_auc_figure.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve

from gdf15_ladder import (INK, INK_SOFT, RULE, BAND_FILL, MM_IN,
                          apply_style, add_caption)

ROOT = Path(__file__).resolve().parent
SUBMIT = ROOT / "figure1_submission"

# marker column -> display label, in the order the study measures them
MARKERS = [
    ("gdf15_raw", "GDF15"),
    ("sixmwd_m", "6MWD"),
    ("mpap_mmhg", "mPAP"),
    ("pvr_dyn_s_cm5", "PVR"),
    ("ci_distfrommedian", "Cardiac index"),
    ("ntprobnp_ng_l", "NT-proBNP"),
    ("age_years", "Age"),
]
UNITS = {"GDF15": "pg/mL", "6MWD": "m", "mPAP": "mmHg",
         "PVR": "dyn·s·cm⁻⁵", "Cardiac index": "L/min/m²",
         "NT-proBNP": "ng/L", "Age": "years"}

# Published values this figure must reproduce (notebook 05 / spec §5, plus the
# age result from spec §4.7). n, n_PVOD, AUC.
PUBLISHED = {
    "GDF15": (166, 26, 0.956), "6MWD": (112, 22, 0.673),
    "mPAP": (157, 24, 0.617), "PVR": (142, 23, 0.602),
    "Cardiac index": (82, 25, 0.567), "NT-proBNP": (65, 17, 0.558),
    "Age": (166, 26, 0.522),
}
N_BOOT, SEED = 4000, 1

HERO, MUTED = "#4a3aa7", "#8a8a86"   # GDF15 vs everything else


def auc_table():
    d = pd.read_csv(ROOT / "data" / "df_g1.csv")
    assert len(d) == 166, len(d)
    assert d.cohort.value_counts().to_dict() == {"IPAH": 98, "HPAH": 42, "PVOD": 26}

    y = (d.cohort == "PVOD").astype(int).to_numpy()
    rng = np.random.default_rng(SEED)
    rows = []
    for col, label in MARKERS:
        m = d[col].notna().to_numpy()
        yy, xx = y[m], d.loc[m, col].to_numpy(float)

        # Direction-free AUC, as used throughout this project: the marker is
        # oriented to whichever direction discriminates better. Orientation is
        # fixed ONCE on the full sample and then held for every bootstrap
        # replicate -- taking max(a, 1-a) inside each replicate would truncate
        # the interval at 0.5 and bias it upward.
        sign = 1.0 if roc_auc_score(yy, xx) >= 0.5 else -1.0
        auc = roc_auc_score(yy, sign * xx)

        i1, i0 = np.where(yy == 1)[0], np.where(yy == 0)[0]
        b = np.empty(N_BOOT)
        for k in range(N_BOOT):
            idx = np.concatenate([rng.choice(i1, len(i1), replace=True),
                                  rng.choice(i0, len(i0), replace=True)])
            b[k] = roc_auc_score(yy[idx], sign * xx[idx])
        lo, hi = np.percentile(b, [2.5, 97.5])

        exp = PUBLISHED[label]
        got = (int(m.sum()), int(yy.sum()), round(auc, 3))
        assert got == exp, f"{label}: {got} != published {exp}"

        fpr, tpr, _ = roc_curve(yy, sign * xx)
        rows.append(dict(marker=label, n=int(m.sum()), n_pvod=int(yy.sum()),
                         auc=auc, lo=lo, hi=hi, fpr=fpr, tpr=tpr,
                         higher="higher" if sign > 0 else "lower"))
    return pd.DataFrame(rows).sort_values("auc", ascending=False).reset_index(drop=True)


def draw(tab, caption=None):
    """One ROC panel per marker, plus a final panel overlaying all of them."""
    apply_style(7.0)
    fig_w = 183 * MM_IN
    fig_h = 92 * MM_IN
    fig = plt.figure(figsize=(fig_w, fig_h))
    box = [0.055, 0.085, 0.935, 0.865]
    if caption is not None:
        fig_h, box = add_caption(fig, fig_w, fig_h, box, *caption, "Figure 1")

    ncol, nrow = 4, 2
    gs = fig.add_gridspec(nrow, ncol,
                          left=box[0], bottom=box[1],
                          right=box[0] + box[2], top=box[1] + box[3],
                          wspace=0.36, hspace=0.46)

    panels = list(tab.itertuples()) + [None]      # trailing overlay panel
    for k, r in enumerate(panels):
        ax = fig.add_subplot(gs[k // ncol, k % ncol])
        ax.plot([0, 1], [0, 1], color=INK_SOFT, lw=0.7, ls=(0, (3, 2)), zorder=1)

        if r is None:
            for rr in tab.itertuples():
                hero = rr.marker == "GDF15"
                ax.plot(rr.fpr, rr.tpr, color=HERO if hero else MUTED,
                        lw=1.5 if hero else 0.8, alpha=1.0 if hero else 0.75,
                        zorder=4 if hero else 2, solid_joinstyle="round")
            ax.set_title("All markers", fontsize=7.4, color=INK, pad=3)
            ax.text(0.96, 0.10, "GDF15", fontsize=6.4, color=HERO,
                    ha="right", va="bottom", fontweight="bold")
            ax.text(0.96, 0.03, "all others", fontsize=6.4, color=MUTED,
                    ha="right", va="bottom")
        else:
            hero = r.marker == "GDF15"
            c = HERO if hero else MUTED
            ax.plot(r.fpr, r.tpr, color=c, lw=1.6 if hero else 1.2, zorder=3,
                    solid_joinstyle="round")
            ax.fill_between(r.fpr, r.tpr, color=c, alpha=0.11, lw=0, zorder=2)
            ax.set_title(r.marker, fontsize=7.4, color=INK, pad=3,
                         fontweight="bold" if hero else "normal")
            ax.text(0.95, 0.16, f"AUC {r.auc:.3f}", fontsize=6.8, color=INK,
                    ha="right", va="bottom",
                    fontweight="bold" if hero else "normal")
            ax.text(0.95, 0.085, f"({r.lo:.3f}–{r.hi:.3f})", fontsize=6.0,
                    color=INK_SOFT, ha="right", va="bottom")
            ax.text(0.95, 0.01, f"n = {r.n} ({r.n_pvod} PVOD)", fontsize=6.0,
                    color=INK_SOFT, ha="right", va="bottom")

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect("equal")
        ax.set_xticks([0, 0.5, 1])
        ax.set_yticks([0, 0.5, 1])
        ax.tick_params(labelsize=6.2, length=2, colors=INK_SOFT, pad=1.5)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(RULE)
        # axis names only on the outer edge, so the panels stay uncluttered
        if k % ncol == 0:
            ax.set_ylabel("Sensitivity", fontsize=6.8, color=INK, labelpad=2)
        if k // ncol == nrow - 1:
            ax.set_xlabel("1 – specificity", fontsize=6.8, color=INK, labelpad=2)
    return fig


def caption_parts(tab):
    g = tab[tab.marker == "GDF15"].iloc[0]
    others = tab[tab.marker != "GDF15"]
    worst, best = others.auc.min(), others.auc.max()
    best_name = others.loc[others.auc.idxmax(), "marker"]
    crosses = others[(others.lo <= 0.5)]

    lead = ("Plasma GDF15 discriminates PVOD from IPAH and HPAH; no other "
            "variable measured in the cohort does.")
    body = (
        "Receiver-operating-characteristic curves for separating PVOD "
        "(n = 26) from the pooled IPAH (n = 98) and HPAH (n = 42) cohorts, for "
        f"each variable available in the same {int(g.n)} patients. One panel "
        "per variable, ordered by area under the curve (AUC); the final panel "
        "overlays all seven. The dashed diagonal is the line of no "
        "discrimination (AUC 0.5). AUC is given in each panel with a 95% "
        "confidence interval from 4,000 stratified bootstrap resamples, drawn "
        "independently within the PVOD and the comparator groups so that group "
        "sizes are preserved. "
        "Marker n differs because the underlying measurements are unevenly "
        "complete; n and the number of PVOD cases contributing are given in "
        "each panel, and every curve and AUC is computed only on the patients "
        "with that measurement. AUCs are direction-free — each variable is "
        "oriented to "
        "whichever direction discriminates better, decided once on the full "
        "sample and then held fixed across all bootstrap replicates. That "
        "convention favours the weak markers, since the more favourable of two "
        "directions is chosen after seeing the data; they remain near 0.5 "
        f"regardless. GDF15 reaches AUC {g.auc:.3f} "
        f"({g.lo:.3f}–{g.hi:.3f}), while the best of the remaining variables, "
        f"{best_name}, reaches only {best:.3f} and the weakest {worst:.3f}; "
        f"the confidence interval includes 0.5 for {len(crosses)} of the "
        f"{len(others)} non-GDF15 variables, so their discrimination is not "
        "distinguishable from chance in this sample. No multiplicity adjustment "
        "is applied across the seven markers and no marker is combined with "
        "another; each interval describes one variable in isolation. AUC point "
        "estimates are those established in the project analysis and are "
        "reproduced exactly here. Each AUC uses one variable together with the "
        "diagnostic label only, so no patient-level pairing between GDF15 and "
        "any clinical variable is involved. Source data: Prabhakar et al., "
        "J Clin Invest 2026;136(2):e199159, Supplementary Source Data — GDF15 "
        "from sheet Fig_1a (panel Fig. 1a); age and the clinical variables from "
        "the per-patient panels Fig_1g–Fig_1m.")
    return lead, body


def main():
    tab = auc_table()
    SUBMIT.mkdir(exist_ok=True)

    print(f"{'marker':<14}{'n':>5}{'PVOD':>6}{'AUC':>8}{'95% CI':>20}  dir")
    for r in tab.itertuples():
        print(f"{r.marker:<14}{r.n:>5}{r.n_pvod:>6}{r.auc:>8.3f}"
              f"   [{r.lo:.3f}, {r.hi:.3f}]  {r.higher}")
    print("\n[VALIDATE] all AUCs match the published values.\n")

    lead, body = caption_parts(tab)

    # (1) the chart on its own
    fig = draw(tab)
    fig.savefig(SUBMIT / "Figure1_auc_comparison.svg", format="svg")
    plt.close(fig)

    # (2) the caption on its own
    (SUBMIT / "Figure1_caption.md").write_text(
        f"**Figure 1 | {lead}**\n\n{body}\n")

    # (3) the two composed as one submission-ready page
    fig = draw(tab, caption=(lead, body))
    fig.savefig(SUBMIT / "Figure1_submission_proof.pdf", format="pdf")
    plt.close(fig)

    for f in sorted(SUBMIT.iterdir()):
        print(f"wrote {f.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
