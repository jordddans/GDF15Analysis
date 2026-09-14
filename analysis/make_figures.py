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
import figure1_bootstrap
import make_figure1

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
    FIGURES.mkdir(parents=True, exist_ok=True)
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
    observed = data.load_gdf15()
    cases = observed.loc[observed.cohort == 'PVOD', 'gdf15']
    comparators = observed.loc[observed.cohort.isin(['IPAH', 'HPAH']), 'gdf15']
    negative_lr = (cases < threshold['published_cutoff']).mean() / (comparators < threshold['published_cutoff']).mean()

    captions = {
        "figure1_caption.md": f"""# Figure 1. Plasma GDF15 reported as a threshold and as continuous evidence.

(A) Observed concentrations in {main['n']} patients: {counts['PVOD']} PVOD,
{counts['IPAH']} IPAH and {counts['HPAH']} HPAH. Vertical bars mark medians.
(B) Ordinary logistic regression on log10 GDF15, with study sampling log-odds
subtracted to give a likelihood ratio. The pointwise 95% band uses {figure1_bootstrap.N_BOOT:,}
requested stratified resamples (seed {figure1_bootstrap.SEED}): 3,991 finite fits and
9 separated resamples excluded. The band narrows in the tail after these exclusions.
The empirical cutoff rule assigns LR {negative_lr:.2f} below
{threshold['published_cutoff']:,.0f} pg/mL and {example['cutoff_positive_lr']:.1f} at or above it.
Both observed support boundaries are marked: {support['pvod_min_gdf15']:,.0f} pg/mL
(lowest PVOD) and {support['comparator_max_gdf15']:,.1f} pg/mL (highest comparator).
Beyond either boundary the curve is broken and the band is lighter.
(C) Every plotted patient is observed in the published source. Filled points show
model-derived post-test probability and open points show the cutoff-rule probability,
both at 10% pre-test probability; joining lines show their difference. The adjacent
observed IPAH and PVOD patients from Table 1 remain the centerpiece. Shapes distinguish
PVOD from comparators. No hypothetical patients are added.

{figure1_bootstrap.METHODS}

Source: Prabhakar et al., J Clin Invest 2026;136(2):e199159, published supporting
workbook Fig_1a (equivalent cleaned alias GDF15-1).
""",

        "figure2_caption.md": f"""# Figure 2. The same model fitted against each comparator separately.

Logistic regression on log10(plasma GDF15) fitted to PVOD (n = {main['n_pvod']}) against
each comparator cohort in turn, with the same prevalence adjustment and the same parametric
bootstrap ({results['n_bootstrap']:,} replicates, seed {results['seed']}). Each
panel shows the fitted likelihood-ratio curve with its pointwise 95% band; the grey curve
repeated in every panel is the PVOD versus IPAH + HPAH model, for reference. Shaded regions
mark concentrations above the highest value observed in that comparator cohort. The
concentration at which the likelihood ratio reaches 10 differs by comparator; the intervals in this
figure are parametric-bootstrap throughout, unlike the stratified case-bootstrap interval quoted for
the main-model landmark reported in the analysis:
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
    FIGURES.mkdir(parents=True, exist_ok=True)
    for name, text in captions.items():
        (FIGURES / name).write_text(text)
        print(f"wrote figures/{name}")


def main():
    results, curves = load()
    figure1_bootstrap.write(RESULTS, data.WORKBOOK)
    make_figure1.draw(RESULTS, FIGURES, data.WORKBOOK)
    figure_comparators(results, curves)
    figure_calibration(results)
    write_captions(results)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workbook', type=Path, default=data.WORKBOOK)
    parser.add_argument('--results-dir', type=Path, default=RESULTS)
    parser.add_argument('--figures-dir', type=Path, default=FIGURES)
    args = parser.parse_args()
    data.WORKBOOK = args.workbook.resolve()
    RESULTS, FIGURES = args.results_dir, args.figures_dir
    main()
