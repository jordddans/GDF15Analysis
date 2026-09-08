"""
Figure 3 — what a clinical classifier built on this data can currently do.

A plain classifier is fitted to the data exactly as it stands and evaluated
fairly. Nothing is added to make it behave: none of the capping, Jeffreys
flooring or extrapolation-refusal machinery in gdf15_lr.py is used, and that
module is deliberately not imported. Equally nothing is tilted to make it look
worse — the thresholds, the fitting and the evaluation are all the ordinary
defaults.

Classifier : logistic regression on log10(GDF15), single predictor, library
             defaults, no tuning, no capping, no flooring.
Evaluation : leave-one-out cross-validation; every probability reported is
             out-of-fold. Intervals are 4,000 stratified bootstrap resamples
             drawn within PVOD and within the comparator so group sizes hold.
Population : PVOD (n = 26) vs IPAH (n = 98) + HPAH (n = 42); n = 166.

Writes three files into JournalFigures/figure3/:
    Figure3_classifier.svg          the chart alone
    Figure3_caption.md              the caption alone
    Figure3_submission_proof.pdf    the two composed as one journal page

Usage:  venv/bin/python gdf15_classifier_figure.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score, roc_curve

from gdf15_ladder import (INK, INK_SOFT, RULE, BAND_FILL, MM_IN,
                          apply_style, add_caption)

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "JournalFigures" / "figure3"

COMPARATOR = ["IPAH", "HPAH"]
N_BOOT, SEED = 4000, 1
HERO, MUTED = "#4a3aa7", "#8a8a86"

# Values established during analysis; asserted so the figure cannot drift.
EXPECT = dict(n=166, n_pvod=26, auc_insample=0.956, auc_loo=0.951,
              cut=1659.0, cut_lo=978.0, cut_hi=2551.0, n_distinct_cuts=13,
              p_ceiling=0.80)


def load():
    d = pd.read_csv(ROOT / "data" / "df_all.csv")
    d["cohort"] = d.cohort.str.replace("\xa0", " ").str.strip()
    sub = d[d.cohort.isin(["PVOD"] + COMPARATOR)].reset_index(drop=True)
    y = (sub.cohort == "PVOD").astype(int).to_numpy()
    x = sub.gdf15_raw.to_numpy(float)
    assert len(y) == EXPECT["n"] and y.sum() == EXPECT["n_pvod"], (len(y), y.sum())
    return y, x


def fit_loo(y, x):
    """Out-of-fold probabilities. Every patient is scored by a model that never saw them."""
    X = np.log10(x).reshape(-1, 1)
    p = np.empty(len(y))
    for tr, te in LeaveOneOut().split(X):
        p[te] = LogisticRegression().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
    p_ins = LogisticRegression().fit(X, y).predict_proba(X)[:, 1]

    auc_loo = roc_auc_score(y, p)
    auc_ins = roc_auc_score(y, p_ins)
    assert round(auc_loo, 3) == EXPECT["auc_loo"], auc_loo
    assert round(auc_ins, 3) == EXPECT["auc_insample"], auc_ins
    assert round(p[y == 1].max(), 2) == EXPECT["p_ceiling"], p[y == 1].max()
    return p, auc_loo, auc_ins


def youden_cut(yy, xx):
    fpr, tpr, thr = roc_curve(yy, xx)
    return float(thr[np.argmax(tpr - fpr)])


def bootstrap(y, x, p):
    """Resample within each group; refit the cutoff and re-measure, 4,000 times.

    The cutoff is re-selected inside every replicate, so its interval carries the
    uncertainty of *choosing* the threshold as well as of estimating it.
    """
    i1, i0 = np.where(y == 1)[0], np.where(y == 0)[0]
    rng = np.random.default_rng(SEED)
    cuts = np.empty(N_BOOT)
    sens = np.empty(N_BOOT)
    spec = np.empty(N_BOOT)
    aucs = np.empty(N_BOOT)
    for b in range(N_BOOT):
        s1 = rng.choice(i1, len(i1), True)
        s0 = rng.choice(i0, len(i0), True)
        idx = np.concatenate([s1, s0])
        c = youden_cut(y[idx], x[idx])
        cuts[b] = c
        sens[b] = (x[y == 1] >= c).mean()
        spec[b] = (x[y == 0] < c).mean()
        aucs[b] = roc_auc_score(y[idx], p[idx])

    cut = youden_cut(y, x)
    lo, hi = np.percentile(cuts, [2.5, 97.5])
    assert round(cut) == EXPECT["cut"], cut
    assert round(lo) == EXPECT["cut_lo"] and round(hi) == EXPECT["cut_hi"], (lo, hi)
    n_distinct = len(np.unique(np.round(cuts)))
    assert n_distinct == EXPECT["n_distinct_cuts"], n_distinct

    return dict(cut=cut, cuts=cuts, cut_lo=lo, cut_hi=hi, n_distinct=n_distinct,
                sens=float((x[y == 1] >= cut).mean()),
                spec=float((x[y == 0] < cut).mean()),
                sens_ci=tuple(np.percentile(sens, [2.5, 97.5])),
                spec_ci=tuple(np.percentile(spec, [2.5, 97.5])),
                auc_ci=tuple(np.percentile(aucs, [2.5, 97.5])))


def ppv_npv(prev, sens, spec):
    ppv = sens * prev / (sens * prev + (1 - spec) * (1 - prev))
    npv = spec * (1 - prev) / (spec * (1 - prev) + (1 - sens) * prev)
    return ppv, npv


def draw(y, x, p, r, auc_loo, auc_ins, caption=None):
    apply_style(7.0)
    fig_w, fig_h = 183 * MM_IN, 74 * MM_IN
    fig = plt.figure(figsize=(fig_w, fig_h))
    box = [0.075, 0.145, 0.910, 0.755]
    if caption is not None:
        fig_h, box = add_caption(fig, fig_w, fig_h, box, *caption, "Figure 3")
    gs = fig.add_gridspec(1, 3, left=box[0], bottom=box[1],
                          right=box[0] + box[2], top=box[1] + box[3], wspace=0.30)

    def tidy(ax):
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(RULE)
        ax.tick_params(labelsize=6.4, length=2, colors=INK_SOFT, pad=2)

    # ── A: what the tool outputs for one patient ────────────────────────────
    ax = fig.add_subplot(gs[0, 0])
    rng = np.random.default_rng(7)
    for k, (lab, mask, col) in enumerate((("IPAH + HPAH", y == 0, MUTED),
                                          ("PVOD", y == 1, HERO))):
        xs = k + rng.uniform(-0.19, 0.19, mask.sum())
        ax.plot(xs, p[mask], linestyle="none", marker="o", ms=2.8, color=col,
                mec="white", mew=0.35, alpha=0.9, zorder=3)
        med = np.median(p[mask])
        ax.plot([k - 0.3, k + 0.3], [med, med], color=INK, lw=1.1, zorder=4)
    ceiling = p[y == 1].max()
    ax.axhline(ceiling, color=INK, lw=0.8, ls=(0, (3, 2)), zorder=2)
    ax.text(1.42, ceiling, f"highest value\nreached: {ceiling:.2f}", fontsize=6.0,
            color=INK, ha="right", va="bottom", linespacing=1.3)
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["IPAH + HPAH\nn = 140", "PVOD\nn = 26"], fontsize=6.8,
                       color=INK, linespacing=1.4)
    ax.tick_params(axis="x", length=0)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_ylabel("Predicted probability of PVOD\n(out-of-fold)", fontsize=7.0,
                  labelpad=3, linespacing=1.4)
    ax.set_title("A   Output for each patient", fontsize=7.4, color=INK,
                 loc="left", pad=5)
    tidy(ax)

    # ── B: where the cutoff lands ───────────────────────────────────────────
    ax = fig.add_subplot(gs[0, 1])
    vals, counts = np.unique(np.round(r["cuts"]), return_counts=True)
    pct = 100 * counts / counts.sum()
    ax.axvspan(r["cut_lo"], r["cut_hi"], color=HERO, alpha=0.10, lw=0, zorder=1)
    ax.vlines(vals, 0, pct, color=HERO, lw=2.6, zorder=3)
    ax.axvline(1658, color=INK, lw=0.9, ls=(0, (3, 2)), zorder=4)
    top = pct.max()
    # 1,658 coincides with the modal cutoff, so the label goes beside the spike
    ax.text(1560, top * 0.97, "published\n1,658", fontsize=6.0, color=INK,
            ha="right", va="top", linespacing=1.3)
    # the interval spans the tallest bars, so it is drawn clear above them
    ax.annotate("", xy=(r["cut_lo"], top * 1.12), xytext=(r["cut_hi"], top * 1.12),
                arrowprops=dict(arrowstyle="<->", color=HERO, lw=0.9))
    ax.text((r["cut_lo"] + r["cut_hi"]) / 2, top * 1.19,
            f"95%: {r['cut_lo']:,.0f}–{r['cut_hi']:,.0f}", fontsize=6.2,
            color=HERO, ha="center", va="bottom", fontweight="bold")
    ax.text(0.02, 0.99, f"only {r['n_distinct']} distinct\ncutoffs reachable",
            transform=ax.transAxes, fontsize=6.0, color=INK_SOFT,
            ha="left", va="top", linespacing=1.3)
    ax.set_xlim(600, 3100)
    ax.set_ylim(0, top * 1.40)
    ax.set_xlabel("Selected cutoff (pg/mL)", fontsize=7.0, labelpad=2)
    ax.set_ylabel("% of 4,000 resamples", fontsize=7.0, labelpad=3)
    ax.set_title("B   Where the cutoff lands", fontsize=7.4, color=INK,
                 loc="left", pad=5)
    tidy(ax)

    # ── C: what that cutoff would mean clinically ───────────────────────────
    ax = fig.add_subplot(gs[0, 2])
    prev = np.linspace(0.01, 0.30, 200)
    at = np.argmin(np.abs(prev - 0.05))     # the curves are far apart here
    for lab, i, col, va in (("PPV", 0, HERO, "top"), ("NPV", 1, MUTED, "bottom")):
        mid = np.array([ppv_npv(q, r["sens"], r["spec"])[i] for q in prev])
        b1 = np.array([ppv_npv(q, r["sens_ci"][0], r["spec_ci"][0])[i] for q in prev])
        b2 = np.array([ppv_npv(q, r["sens_ci"][1], r["spec_ci"][1])[i] for q in prev])
        ax.fill_between(prev * 100, np.minimum(b1, b2) * 100,
                        np.maximum(b1, b2) * 100, color=col, alpha=0.16, lw=0, zorder=1)
        ax.plot(prev * 100, mid * 100, color=col, lw=1.6, zorder=3)
        ax.text(prev[at] * 100 + 1.0, mid[at] * 100, lab, fontsize=7.0, color=col,
                ha="left", va=va, fontweight="bold", zorder=5,
                bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none",
                          alpha=0.85))
    ax.axvline(15.7, color=INK, lw=0.8, ls=(0, (3, 2)), zorder=2)
    ax.text(16.6, 3, "this cohort: 15.7%\n(reflects how it\nwas assembled)",
            fontsize=5.8, color=INK_SOFT, ha="left", va="bottom", linespacing=1.35,
            zorder=5, bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none",
                                alpha=0.85))
    ax.set_xlim(0, 34)
    ax.set_ylim(0, 103)
    ax.set_xticks([0, 10, 20, 30])
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Assumed prevalence of PVOD (%)", fontsize=7.0, labelpad=2)
    ax.set_ylabel("Predictive value (%)", fontsize=7.0, labelpad=3)
    ax.set_title("C   What a positive or negative means", fontsize=7.4, color=INK,
                 loc="left", pad=5)
    tidy(ax)
    return fig


SOURCE = ("Prabhakar et al., J Clin Invest 2026;136(2):e199159, Supplementary "
          "Source Data, sheet Fig_1a (panel Fig. 1a)")


def caption_parts(y, p, r, auc_loo, auc_ins, math=False):
    """Returns (lead, body). One source for both the caption file and the proof."""
    n, n_pv = len(y), int(y.sum())
    lead = ("A diagnostic classifier for PVOD built from plasma GDF15 alone, "
            "trained and evaluated on the present cohort.")
    body = (
        f"Logistic regression on log10(plasma GDF15), a single predictor, fitted "
        f"to {n} patients — PVOD (n = {n_pv}) versus pooled IPAH (n = 98) and "
        "HPAH (n = 42) — with library default settings and no tuning, penalty, "
        "capping or flooring applied to the fitted model or its output. "
        "(a) Predicted probability of PVOD for every patient, grouped by true "
        "diagnosis; horizontal rules are group medians. Every probability shown "
        "is out-of-fold, obtained by leave-one-out cross-validation, so no "
        "patient contributed to the model that scored them. The dashed line "
        f"marks the highest probability the classifier assigned to any patient "
        f"({p[y == 1].max():.2f}). Cross-validated area under the ROC curve was "
        f"{auc_loo:.3f} (95% CI {r['auc_ci'][0]:.3f}–{r['auc_ci'][1]:.3f}); the "
        f"in-sample value was {auc_ins:.3f}. "
        "(b) Distribution of the decision cutoff across 4,000 stratified "
        "bootstrap resamples, drawn independently within the PVOD and the "
        "comparator groups so that group sizes are preserved. The cutoff was "
        "re-selected by the Youden criterion inside each resample, so the "
        "interval shown reflects the uncertainty of choosing the threshold as "
        "well as of estimating it. Bars are the distinct cutoff values the "
        f"procedure selected — {r['n_distinct']} in total, since with "
        f"{n_pv} PVOD cases a cutoff can only fall between adjacent observed "
        f"values. The shaded band and arrow span the 95% interval "
        f"({r['cut_lo']:,.0f}–{r['cut_hi']:,.0f} pg/mL); the point estimate was "
        f"{r['cut']:,.0f} pg/mL. The dashed line marks the 1,658 pg/mL cutoff "
        "published in the source study, which falls at the "
        f"{100 * (r['cuts'] <= 1658).mean():.0f}th percentile of this "
        "distribution. "
        "(c) Positive and negative predictive value at the point-estimate "
        f"cutoff (sensitivity {r['sens']:.2f}, 95% CI {r['sens_ci'][0]:.2f}–"
        f"{r['sens_ci'][1]:.2f}; specificity {r['spec']:.2f}, 95% CI "
        f"{r['spec_ci'][0]:.2f}–{r['spec_ci'][1]:.2f}), plotted against assumed "
        "disease prevalence, with bands propagated from those intervals. "
        f"PVOD comprises {100 * y.mean():.1f}% of this cohort, a consequence of "
        "how the cohort was assembled rather than an estimate of prevalence in "
        "any clinical population; predictive values are therefore shown across a "
        "range of assumed prevalences rather than evaluated at a single one. All "
        "intervals are percentile bootstrap intervals from 4,000 resamples. "
        f"Source data: {SOURCE}.")
    return lead, body


def main():
    y, x = load()
    p, auc_loo, auc_ins = fit_loo(y, x)
    r = bootstrap(y, x, p)
    OUT.mkdir(parents=True, exist_ok=True)

    print(f"  n = {len(y)}   PVOD = {int(y.sum())}   prevalence {y.mean():.3f}")
    print(f"  AUC   in-sample {auc_ins:.3f}   LOO-CV {auc_loo:.3f} "
          f"[{r['auc_ci'][0]:.3f}, {r['auc_ci'][1]:.3f}]")
    print(f"  out-of-fold p   PVOD median {np.median(p[y==1]):.2f}  max {p[y==1].max():.2f}"
          f"   |  comparator median {np.median(p[y==0]):.2f}  max {p[y==0].max():.2f}")
    print(f"  cutoff  {r['cut']:,.0f} pg/mL   95% [{r['cut_lo']:,.0f}, {r['cut_hi']:,.0f}]"
          f"   distinct values {r['n_distinct']}")
    print(f"  sens {r['sens']:.2f} [{r['sens_ci'][0]:.2f}, {r['sens_ci'][1]:.2f}]"
          f"   spec {r['spec']:.2f} [{r['spec_ci'][0]:.2f}, {r['spec_ci'][1]:.2f}]")
    for q in (0.05, 0.157):
        a, b = ppv_npv(q, r["sens"], r["spec"])
        print(f"  at prevalence {q*100:>4.1f}%  PPV {a*100:.1f}%  NPV {b*100:.1f}%")
    print("\n  [VALIDATE] all assertions passed.\n")

    lead, body = caption_parts(y, p, r, auc_loo, auc_ins)

    fig = draw(y, x, p, r, auc_loo, auc_ins)
    fig.savefig(OUT / "Figure3_classifier.svg", format="svg")
    plt.close(fig)

    (OUT / "Figure3_caption.md").write_text(f"**Figure 3 | {lead}**\n\n{body}\n")

    fig = draw(y, x, p, r, auc_loo, auc_ins, caption=(lead, body))
    fig.savefig(OUT / "Figure3_submission_proof.pdf", format="pdf")
    plt.close(fig)

    for f in sorted(OUT.iterdir()):
        print(f"  wrote {f.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
