"""Draw the figures from the values written by ``run_analysis.py``.

Every number shown is read from ``results/`` -- nothing is recomputed here, so the
figures cannot drift from the analysis.

Usage:  python analysis/make_figures.py
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np
import pandas as pd

import data

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

MM = 1 / 25.4
HERO = "#4a3aa7"          # the fitted curve
COMPARATOR_COLOURS = {"IPAH": "#2a78d6", "HPAH": "#eb6834",
                      "Other PH": "#1baf7a", "COPD": "#4a3aa7"}
INK, INK_SOFT, RULE = "#1a1a19", "#52514e", "#d8d7d2"
BAND, SHADE = "#4a3aa7", "#f0efe9"

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 7,
    "axes.labelsize": 7.5, "axes.titlesize": 8, "legend.fontsize": 6.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.edgecolor": INK_SOFT, "axes.linewidth": 0.6,
    "xtick.color": INK_SOFT, "ytick.color": INK_SOFT,
    "text.color": INK, "axes.labelcolor": INK,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def round10(value):
    """Round a concentration to the nearest 10 pg/mL, for reader-facing text only.

    Bootstrap endpoints are quoted to this precision in figures, captions and the README;
    ``results/`` keeps them at full precision.
    """
    return int(round(value / 10.0) * 10)


def load():
    results = json.loads((RESULTS / "results.json").read_text())
    curves = pd.read_csv(RESULTS / "lr_curves.csv")
    return results, curves


def save(fig, name):
    FIGURES.mkdir(exist_ok=True)
    for suffix in ("pdf", "png"):
        # no creation timestamp, so a rerun reproduces the file byte for byte
        fig.savefig(FIGURES / f"{name}.{suffix}", dpi=400,
                    metadata={"CreationDate": None} if suffix == "pdf" else None)
    plt.close(fig)
    print(f"wrote figures/{name}.pdf and .png")


def style_lr_axes(ax, xlim=(300, 6000), ylim=(0.004, 2000)):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xticks([500, 1000, 2000, 5000])
    ax.set_xticklabels(["500", "1,000", "2,000", "5,000"])
    ax.set_yticks([0.01, 0.1, 1, 10, 100, 1000])
    ax.set_yticklabels(["0.01", "0.1", "1", "10", "100", "1,000"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(which="minor", length=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def figure_lr_curve(results, curves):
    """Primary figure: the continuous likelihood-ratio curve for PVOD vs IPAH+HPAH."""
    entry = results["models"]["IPAH+HPAH"]
    support = entry["support"]
    lr10 = entry["gdf15_at_lr"]["10"]
    threshold = entry["threshold_comparison"]
    curve = curves[curves.model == "IPAH+HPAH"]
    gdf15 = data.load_gdf15()

    fig, (ax, rug) = plt.subplots(
        2, 1, figsize=(120 * MM, 96 * MM), sharex=True,
        gridspec_kw=dict(height_ratios=[6, 1], hspace=0.08))
    style_lr_axes(ax)

    lower, upper = support["pvod_min_gdf15"], support["comparator_max_gdf15"]
    for axis in (ax, rug):
        axis.axvspan(300, lower, color=SHADE, zorder=0)
        axis.axvspan(upper, 6000, color=SHADE, zorder=0)
        for boundary in (lower, upper):
            axis.axvline(boundary, color=INK_SOFT, lw=0.7, ls=(0, (3, 2)), zorder=1)

    for level in (1, 2, 5, 10):
        ax.axhline(level, color=RULE, lw=0.6, zorder=1)
        ax.text(312, level, f"LR {level}", color=INK_SOFT, fontsize=6,
                va="center", ha="left",
                bbox=dict(fc=SHADE, ec="none", pad=0.8), zorder=4)

    ax.fill_between(curve.gdf15, curve.lr_ci_low, curve.lr_ci_high,
                    color=BAND, alpha=0.16, lw=0, zorder=2)
    ax.plot(curve.gdf15, curve.lr, color=HERO, lw=1.6, zorder=3)

    ax.plot([lr10["ci_low"], lr10["ci_high"]], [10, 10], color=HERO, lw=1.3,
            solid_capstyle="butt", zorder=5)
    ax.plot([lr10["gdf15"]], [10], "o", ms=4.5, color=HERO, mec="white", mew=0.8, zorder=6)
    ax.annotate(f"LR = 10 at {round10(lr10['gdf15']):,} pg/mL\n"
                f"(95% CI {round10(lr10['ci_low']):,}–{round10(lr10['ci_high']):,})",
                xy=(lr10["ci_high"], 10), xytext=(10, -22), textcoords="offset points",
                color=INK, fontsize=6.5, ha="left", va="top", zorder=6)

    published = threshold["published_cutoff"]
    lr_published = threshold["lr_at_published_cutoff"]
    ax.plot([published], [lr_published], "D", ms=3.4, color="white",
            mec=INK_SOFT, mew=0.9, zorder=6)
    ax.annotate(f"1,658 pg/mL: LR {lr_published:.1f}",
                xy=(published, lr_published), xytext=(-9, 0),
                textcoords="offset points", color=INK_SOFT, fontsize=6.5,
                ha="right", va="center", zorder=5)

    ax.text(lower * 0.93, 1300,
            f"below ~{lower:,.0f} pg/mL\nmodel-supported beyond\nobserved PVOD data",
            fontsize=6, color=INK_SOFT, ha="right", va="top")
    ax.text(upper * 1.07, 1300,
            f"above ~{upper:,.1f} pg/mL\nmodel-supported beyond observed\n"
            f"IPAH/HPAH comparator data",
            fontsize=6, color=INK_SOFT, ha="left", va="top")

    ax.set_ylabel("Likelihood ratio for PVOD")
    ax.set_title("PVOD versus IPAH + HPAH  ·  "
                 f"n = {entry['n']}, {entry['n_pvod']} PVOD", loc="left", color=INK)

    # observed values, so the reader can see where the data actually are
    pvod = gdf15.loc[gdf15.cohort == "PVOD", "gdf15"]
    comparator = gdf15.loc[gdf15.cohort.isin(["IPAH", "HPAH"]), "gdf15"]
    rug.plot(pvod, np.full(len(pvod), 1.0), "|", ms=6, mew=0.8, color=HERO, zorder=3)
    rug.plot(comparator, np.full(len(comparator), 0.0), "|", ms=6, mew=0.8,
             color=INK_SOFT, alpha=0.8, zorder=3)
    rug.set_ylim(-0.7, 1.7)
    rug.set_yticks([1.0, 0.0])
    rug.set_yticklabels([f"PVOD ({len(pvod)})", f"IPAH + HPAH ({len(comparator)})"],
                        fontsize=6)
    rug.tick_params(axis="y", length=0)
    rug.tick_params(which="minor", length=0)
    for spine in ("top", "right", "left"):
        rug.spines[spine].set_visible(False)
    rug.set_xlabel("Plasma GDF15 (pg/mL)")
    save(fig, "figure1_lr_curve")


def figure_comparators(results, curves):
    """Supporting figure: the same model fitted against each comparator separately."""
    order = ["IPAH", "HPAH", "Other PH", "COPD"]
    main = curves[curves.model == "IPAH+HPAH"]

    fig, axes = plt.subplots(2, 2, figsize=(150 * MM, 105 * MM), sharex=True, sharey=True)
    for ax, name in zip(axes.ravel(), order):
        entry = results["models"][name]
        curve = curves[curves.model == name]
        lr10 = entry["gdf15_at_lr"]["10"]
        colour = COMPARATOR_COLOURS[name]
        style_lr_axes(ax)
        ax.axhline(10, color=RULE, lw=0.6, zorder=1)
        ax.plot(main.gdf15, main.lr, color=RULE, lw=1.0, zorder=2)
        ax.fill_between(curve.gdf15, curve.lr_ci_low, curve.lr_ci_high,
                        color=colour, alpha=0.15, lw=0, zorder=3)
        ax.plot(curve.gdf15, curve.lr, color=colour, lw=1.4, zorder=4)
        ax.plot([lr10["ci_low"], lr10["ci_high"]], [10, 10], color=colour, lw=1.1, zorder=5)
        ax.plot([lr10["gdf15"]], [10], "o", ms=4, color=colour, mec="white", mew=0.7, zorder=6)
        ax.text(0.03, 0.95, f"PVOD vs {name}", transform=ax.transAxes, fontsize=7.5,
                color=colour, va="top", fontweight="bold")
        ax.text(0.03, 0.855,
                f"n = {entry['n']} ({entry['n_pvod']} PVOD)\n"
                f"LR = 10 at {lr10['gdf15']:,.0f} pg/mL\n"
                f"95% CI {lr10['ci_low']:,.0f}–{lr10['ci_high']:,.0f}",
                transform=ax.transAxes, fontsize=6.2, color=INK_SOFT, va="top")
        ax.axvspan(entry["comparator_max_gdf15"], 6000, color=SHADE, zorder=0)

    for ax in axes[1]:
        ax.set_xlabel("Plasma GDF15 (pg/mL)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Likelihood ratio for PVOD")
    controls = results["healthy_controls"]
    fig.text(0.0, -0.035,
             f"Grey curve in each panel is the PVOD vs IPAH + HPAH model. Shaded region: "
             f"no comparator patient observed above that concentration.\n"
             f"PVOD versus healthy controls ({controls['n_controls']}) cannot be fitted: "
             f"the highest control value ({controls['control_max_gdf15']:,.0f} pg/mL) lies below "
             f"the lowest PVOD value ({controls['pvod_min_gdf15']:,.0f} pg/mL), so the groups are "
             f"perfectly separated and\nthe maximum-likelihood coefficient is unbounded.",
             fontsize=6.2, color=INK_SOFT, ha="left", va="top")
    fig.tight_layout()
    save(fig, "figure2_comparators")


def figure_calibration(results):
    """Observed versus predicted PVOD frequency, by fifth of predicted probability."""
    entry = results["models"]["IPAH+HPAH"]
    rows = pd.DataFrame(entry["calibration_quintiles"])
    position = np.arange(1, 6)

    fig, ax = plt.subplots(figsize=(108 * MM, 74 * MM))
    ax.vlines(position, rows.observed_ci_low, rows.observed_ci_high,
              color=HERO, lw=1.0, zorder=2)
    ax.plot(position, rows.observed_proportion, "o", ms=5.5, color=HERO,
            mec="white", mew=0.8, zorder=3, label="observed proportion with PVOD")
    ax.plot(position, rows.predicted_mean, "D", ms=4.5, color="white", mec=INK_SOFT,
            mew=1.0, zorder=4, label="mean predicted probability")

    ax.set_xticks(position)
    ax.set_xticklabels([f"Q{int(row.quintile)}\n{row.gdf15_min:,.0f}–{row.gdf15_max:,.0f}"
                        for _, row in rows.iterrows()])
    ax.set_xlim(0.5, 5.5)
    ax.set_ylim(-0.13, 1.0)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.spines["left"].set_bounds(0, 1.0)
    ax.set_xlabel("Quintile of predicted probability (GDF15 range, pg/mL)")
    ax.set_ylabel("Probability of PVOD")
    ax.set_title("Calibration by quintile of predicted probability", loc="left", color=INK)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for x, row in zip(position, rows.itertuples()):
        ax.text(x, -0.09, f"{row.observed_events}/{row.n}", fontsize=6.2,
                color=INK_SOFT, ha="center", va="center")
    ax.legend(loc="upper left", frameon=False, handletextpad=0.4)
    ax.text(0.02, 0.72, "PVOD cases / n shown below each quintile",
            transform=ax.transAxes, fontsize=6.2, color=INK_SOFT, va="top")
    save(fig, "figure3_calibration")


SOURCE = ("Source data: Prabhakar et al., J Clin Invest 2026;136(2):e199159, "
          "Supporting data values, sheet GDF15-1.")


def write_captions(results):
    """Write the figure captions, with every number taken from results.json."""
    main = results["models"]["IPAH+HPAH"]
    counts = results["cohort_counts"]
    support = main["support"]
    lr10 = main["gdf15_at_lr"]["10"]
    threshold = main["threshold_comparison"]
    example = main["published_cutoff_example"]
    controls = results["healthy_controls"]
    calibration = main["calibration_quintiles"]
    validation = main["internal_validation"]

    captions = {
        "figure1_caption.md": f"""# Figure 1. Likelihood ratio for PVOD as a continuous function of plasma GDF15.

Logistic regression on log10(plasma GDF15), a single predictor, fitted to
{main['n']} patients \u2014 PVOD (n = {main['n_pvod']}) versus pooled IPAH
(n = {counts['IPAH']}) and HPAH (n = {counts['HPAH']}). The fitted log-odds are converted to a likelihood ratio by subtracting the
log-odds of the sampling prevalence of this cohort
({main['sampling_prevalence']:.1%}, {main['sampling_log_odds']:.4f} on the log-odds scale),
which removes the cohort's case mix from the estimate; the resulting curve is a property of
the two conditional GDF15 distributions and carries no assumption about prevalence in any
clinical population. A clinician supplies the pre-test probability and multiplies its odds
by the likelihood ratio read from this curve. Intercept {main['intercept']:.3f}
(SE {main['intercept_se']:.3f}); log10(GDF15) coefficient {main['slope_log10_gdf15']:.3f}
(SE {main['slope_se']:.3f}, 95% CI {main['slope_ci'][0]:.2f} to {main['slope_ci'][1]:.2f});
{main['n_pvod']} events for one predictor.

The band is a pointwise 95% interval from a parametric bootstrap: outcomes were simulated
from the fitted model at the observed GDF15 values and the model refitted
{results['n_bootstrap']:,} times (seed {results['seed']}), so the interval reflects
uncertainty in the fitted parameters and is defined across the whole measured range rather
than only where cases were observed. Horizontal rules mark likelihood ratios of 1, 2, 5 and
10. The likelihood ratio reaches 10 at {round10(lr10['gdf15']):,} pg/mL
(95% CI {round10(lr10['ci_low']):,}–{round10(lr10['ci_high']):,}). That interval is a stratified
case bootstrap — PVOD and comparator patients resampled separately at their observed sample sizes,
{results['n_bootstrap']:,} replicates, seed {results['seed']} — which is the scheme under which this
landmark and a selected threshold can be compared on identical footing (see below). The parametric
bootstrap gives {round10(lr10['ci_low_parametric']):,}–{round10(lr10['ci_high_parametric']):,} for the
same landmark and is reported as a supplementary sensitivity analysis. Bootstrap endpoints are
quoted to the nearest 10 pg/mL; `results/results.json` holds them at full precision.

Dashed lines mark the limits of the observed data. The region below
~{support['pvod_min_gdf15']:,.0f} pg/mL is labelled on the figure as model-supported beyond observed
PVOD data: no PVOD patient was measured below {support['pvod_min_gdf15']:,.3f} pg/mL, so the curve
there is the fitted model form rather than observation. The region above
~{support['comparator_max_gdf15']:,.1f} pg/mL is labelled as model-supported beyond observed
IPAH/HPAH comparator data: no comparator patient was measured above that value, and
{support['pvod_above_comparator_max']} of {main['n_pvod']} PVOD patients
({support['pvod_above_comparator_max_fraction']:.0%}) lie above that value; the
concentration at which the likelihood ratio reaches 10 therefore falls in a region
containing no comparator observations, which is why the band widens there. Rug marks below
the axis show every measured value in each group. The open diamond marks
{threshold['published_cutoff']:,.0f} pg/mL, the cutoff reported for these data, which
carries a likelihood ratio of {threshold['lr_at_published_cutoff']:.1f}
({example['continuous_lr_ci'][0]:.2f}–{example['continuous_lr_ci'][1]:.2f}). Read as a point on
this curve, that measurement takes a {example['by_pre_test'][1]['pre_test_probability']:.0%}
pre-test probability to {example['by_pre_test'][1]['continuous_post_test']:.1%}
({example['by_pre_test'][1]['continuous_post_test_ci'][0]:.1%}–{example['by_pre_test'][1]['continuous_post_test_ci'][1]:.1%});
read as a positive result from the dichotomy, whose likelihood ratio is
{example['cutoff_positive_lr']:.1f}, it would take the same patient to
{example['by_pre_test'][1]['cutoff_post_test']:.1%}. Both readings use the same patients and the
same outcome and differ only in what is done with the measurement.

{SOURCE}
""",

        "figure2_caption.md": f"""# Figure 2. The same model fitted against each comparator separately.

Logistic regression on log10(plasma GDF15) fitted to PVOD (n = {main['n_pvod']}) against
each comparator cohort in turn, with the same prevalence adjustment and the same parametric
bootstrap ({results['n_bootstrap']:,} replicates, seed {results['seed']}) as Figure 1. Each
panel shows the fitted likelihood-ratio curve with its pointwise 95% band; the grey curve
repeated in every panel is the PVOD versus IPAH + HPAH model, for reference. Shaded regions
mark concentrations above the highest value observed in that comparator cohort. The
concentration at which the likelihood ratio reaches 10 differs by comparator; the intervals in this
figure are parametric-bootstrap throughout, unlike the stratified case-bootstrap interval quoted for
the main model in Figure 1:
{", ".join(f"{name} {results['models'][name]['gdf15_at_lr']['10']['gdf15']:,.0f} pg/mL (95% CI {results['models'][name]['gdf15_at_lr']['10']['ci_low']:,.0f}–{results['models'][name]['gdf15_at_lr']['10']['ci_high']:,.0f})" for name in ["IPAH", "HPAH", "Other PH", "COPD"])};
the same measurement therefore carries different evidential weight depending on which
diagnosis it is being weighed against.

PVOD versus healthy controls cannot be fitted. The highest control value
({controls['control_max_gdf15']:,.1f} pg/mL, n = {controls['n_controls']}) lies below the
lowest PVOD value ({controls['pvod_min_gdf15']:,.1f} pg/mL), so the two groups are perfectly
separated, the log10(GDF15) coefficient has no finite maximum-likelihood estimate, and no
curve exists to plot.

{SOURCE}
""",

        "figure3_caption.md": f"""# Figure 3. Calibration of the fitted model by quintile of predicted probability.

The {main['n']} patients of the PVOD versus IPAH + HPAH model, ordered by predicted
probability and split into fifths. Open diamonds are the mean predicted probability within
each fifth; filled circles are the observed proportion with PVOD, with exact
(Clopper–Pearson) 95% intervals. Counts below each fifth are PVOD cases over patients;
the GDF15 range spanned by each fifth is given on the axis. Observed proportions:
{", ".join(f"Q{row['quintile']} {row['observed_events']}/{row['n']}" for row in calibration)}.
Predicted probabilities are those of the fitted model and therefore apply at the sampling
prevalence of this cohort ({main['sampling_prevalence']:.1%}); the prevalence-independent
form of the model is the likelihood-ratio curve in Figure 1.

Internal validation is a bootstrap optimism correction ({validation['n_replicates']:,} replicates,
seed {results['seed']}), reported apparent against optimism-corrected. C statistic:
{validation['apparent_auc']:.4f} apparent, {validation['corrected_auc']:.4f} corrected
(optimism {validation['auc_optimism']:+.4f}). Calibration slope:
{validation['apparent_calibration_slope']:.4f} apparent, {validation['corrected_calibration_slope']:.4f}
corrected (optimism {validation['calibration_slope_optimism']:+.4f}). The apparent slope is 1 by
construction, so the corrected value is what the resampling has to say about it.

{SOURCE}
""",
    }
    FIGURES.mkdir(exist_ok=True)
    for name, text in captions.items():
        (FIGURES / name).write_text(text)
        print(f"wrote figures/{name}")


def main():
    results, curves = load()
    figure_lr_curve(results, curves)
    figure_comparators(results, curves)
    figure_calibration(results)
    write_captions(results)


if __name__ == "__main__":
    main()
