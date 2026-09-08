"""Run the GDF15 -> PVOD likelihood-ratio analysis and write every reported number.

Outputs (all under ``results/``):

    results.json            every value the figures, tables and manuscript use
    model_coefficients.csv  fitted coefficients for the main and comparator models
    comparator_models.csv   per-comparator LR = 10 concentrations
    calibration.csv         observed vs predicted PVOD frequency by quintile
    lr_curves.csv           the plotted LR curves with their confidence bands
    marker_auc.csv          GDF15 against the other markers measured in the same cohort
    auc_reconciliation.csv  every AUC in circulation for these data, with its calculation
    lr_reference.csv        likelihood ratio at reference GDF15 concentrations
    internal_validation.csv apparent and optimism-corrected AUC and calibration slope
    leave_one_case_out.csv  the 26 refits, each with one PVOD case removed
    functional_form.csv     linear against linear-plus-quadratic log10(GDF15)
    marker_auc_complete_case.csv  the clinical comparators on the complete-case subset
    complete_case_patients.csv    the Patient_IDs making up that subset
    clinical_completeness.csv     which measurements exist, by diagnosis
    ntprobnp_sources.csv    NT-proBNP per PVOD patient in each of the two source sheets

Figures and their captions are drawn from these outputs by ``make_figures``, which this
script calls once the tables are written.

Patient-level outputs carry ``Patient_ID``.  Outputs derived from sheet ``GDF15-1`` do not
and cannot: that sheet is anonymous group-level concentrations, so ``leave_one_case_out.csv``
identifies its cases by source-sheet position and GDF15 value, which is not a patient.

Usage:  python analysis/run_analysis.py
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

import data
import make_figures
import model

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"

COMPARATORS = {
    "IPAH+HPAH": ["IPAH", "HPAH"],
    "IPAH": ["IPAH"],
    "HPAH": ["HPAH"],
    "Other PH": ["Other PH"],
    "COPD": ["COPD"],
}
MAIN = "IPAH+HPAH"
REFERENCE_GDF15 = [500, 861, 1000, 1200, 1658, 1700, 2000, 2500, 3000]
REFERENCE_LR = [1, 2, 5, 10]
CURVE_GRID = np.geomspace(300, 6000, 400)
PUBLISHED_CUTOFF = 1658.0

# AUC values printed in the source study (Prabhakar et al., J Clin Invest 2026;136(2):e199159,
# Figure 1, panels D and E), as percentages with their 95% confidence intervals.  The legend
# states the method: "C statistics to calculate the AUC, followed by Z tests".
#
# Provenance differs between the two panels and is recorded per value.  The panel D values were
# read from the published figure directly.  The panel E values are those a co-author of the source
# study reports from the 600-dpi source render; this analysis has not inspected that render, so
# they are carried as reported and their local verification is pending.
PANEL_D = "read from the published figure"
PANEL_E = ("reported by a co-author of the source study from the 600-dpi source render; "
           "local verification pending")
PUBLISHED_AUC = {
    "IPAH": (94.94, 90.35, 99.52, "Figure 1D", PANEL_D),
    "Other PH": (97.26, 94.14, 100.00, "Figure 1D", PANEL_D),
    "COPD": (95.19, 90.19, 100.00, "Figure 1D", PANEL_D),
    "IPAH+HPAH": (95.60, 91.64, 99.57, "Figure 1E", PANEL_E),
    "IPAH+HPAH+Other PH": (96.05, 92.45, 99.66, "Figure 1E", PANEL_E),
}
AUC_COHORTS = {
    "IPAH": ["IPAH"], "HPAH": ["HPAH"], "Other PH": ["Other PH"], "COPD": ["COPD"],
    "IPAH+HPAH": ["IPAH", "HPAH"], "IPAH+HPAH+Other PH": ["IPAH", "HPAH", "Other PH"],
}
MARKERS = [("gdf15_ageadj", "Age-adjusted GDF15", "ageadj_gdf15"),
           ("sixmwd", "6MWD", "sixmwd_m"),
           ("mpap", "mPAP", "mpap"),
           ("pvr", "PVR", "pvr"),
           ("ci", "Cardiac index", "ci"),
           ("ntprobnp", "NT-proBNP", "ntprobnp"),
           ("age", "Age", "age_years")]

# The complete-case table carries the clinical comparators only.  Age-adjusted GDF15 is
# excluded deliberately: on this subset the two groups do not overlap, so its AUC is 1.00,
# which is the arithmetic of complete separation and not a measure of performance.  Printed
# beside Table 1's 0.956 it would read as the marker improving on the subset, when what
# changed is both the patients and the variable.  The separation is reported as two ranges
# in ``ageadj_gdf15_separation`` instead.
CLINICAL_MARKERS = [entry for entry in MARKERS if entry[0] != "gdf15_ageadj"]

# Availability of these five defines the complete-case subset.  Stated explicitly rather
# than taken as "every marker column": age and age-adjusted GDF15 are present for all 166,
# so including them would not change n, but the subset is a statement about which clinical
# measurements a patient has and should say so.
COMPLETE_CASE_COLUMNS = ["sixmwd_m", "mpap", "pvr", "ci", "ntprobnp"]
PRE_TEST = (0.10, 0.20)   # pre-test probabilities used in the worked examples
N_SEEDS = 10          # seeds over which bootstrap-interval sensitivity is recorded


def build_model(gdf15, comparator):
    subset = gdf15[gdf15.cohort.isin(["PVOD"] + COMPARATORS[comparator])]
    return model.fit_lr_model(subset.gdf15, subset.cohort == "PVOD",
                              f"PVOD vs {comparator}")


def auc_ci_bootstrap(scores, is_case, n_boot=model.N_BOOT, seed=model.SEED):
    """Stratified percentile bootstrap interval for an AUC."""
    scores, is_case = np.asarray(scores, float), np.asarray(is_case).astype(bool)
    rng = np.random.default_rng(seed)
    cases, controls = np.where(is_case)[0], np.where(~is_case)[0]
    values = []
    for _ in range(n_boot):
        take = np.concatenate([rng.choice(cases, len(cases), replace=True),
                               rng.choice(controls, len(controls), replace=True)])
        values.append(model.auc(scores[take], is_case[take]))
    return float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))


def auc_ci_provenance(asymptotic):
    """Exactly how the asymptotic AUC interval in this repository is calculated.

    Recorded because that interval is very close to -- but not the same as -- the one
    printed in the source study's Figure 1E, and the difference between computing a value
    and copying one has to be visible rather than asserted.
    """
    published_low, published_high = PUBLISHED_AUC[MAIN][1] / 100, PUBLISHED_AUC[MAIN][2] / 100
    return {
        "function": "model.auc_ci_asymptotic, in analysis/model.py",
        "statistic": ("the C statistic computed from the Mann-Whitney statistic with "
                      "midranks, so tied GDF15 values contribute one half each"),
        "variance": ("the structural-component estimate. v10 is, for each case, the "
                     "probability that it outranks a random comparator; v01 is the "
                     "complement for each comparator; the variance of the C statistic is "
                     "v10.var(ddof=1) / m + v01.var(ddof=1) / n"),
        "interval": ("area +/- 1.959964 * sqrt(variance): a normal interval formed on "
                     "the AUC scale, with no transformation"),
        "bounded": ("not by construction. The normal interval is formed first and then "
                    "clipped to [0, 1], so an endpoint reported as 1.000 has been "
                    "truncated rather than held below 1 by the method itself. For this "
                    f"comparison no clipping occurs -- the upper endpoint is "
                    f"{asymptotic[1]:.5f} -- but three comparator rows in "
                    "auc_reconciliation.csv do reach 1.0 and are truncated there"),
        "independently_computed": True,
        "computed_interval": [asymptotic[0], asymptotic[1]],
        "published_panel_e_interval": [published_low, published_high],
        "difference_from_published": [asymptotic[0] - published_low,
                                      asymptotic[1] - published_high],
        "evidence_it_was_computed": (
            "the two intervals disagree in the fourth decimal place "
            f"({asymptotic[0]:.5f} against {published_low:.4f}, {asymptotic[1]:.5f} "
            f"against {published_high:.4f}). A transcribed interval would agree exactly; "
            "these agree only as closely as the same data and the same family of "
            "calculation force them to"),
        "verification": ("verify.py recomputes the variance from the definition of the "
                         "two structural components and requires it to match"),
        "other_calculations": (
            "a different asymptotic calculation, reported during review, returns 0.900 "
            "to 1.012 for the same comparison. That interval is unbounded and crosses 1, "
            "which an AUC cannot exceed. It is neither the calculation used here nor "
            "reproduced here, and is recorded only to note that asymptotic AUC intervals "
            "differ by estimator"),
        "reported_primary": ("neither. The reported interval is the stratified "
                             "percentile bootstrap, which cannot leave [0, 1]"),
    }


def auc_reconciliation(gdf15, n_seeds=N_SEEDS):
    """Recompute every AUC printed for these data, using the source study's stated method.

    The study reports the AUC as a C statistic with an asymptotic (Z-test) interval, which
    is what :func:`model.auc_ci_asymptotic` computes.  Recomputing each published comparison
    from the supporting data shows which printed values the data reproduce and which they
    do not.  The stratified percentile bootstrap is reported alongside for the main model,
    with its sensitivity to the random seed, because a bootstrap interval on 26 cases moves
    by more than its printed precision from one seed to the next.
    """
    rows = []
    for name, cohorts in AUC_COHORTS.items():
        subset = gdf15[gdf15.cohort.isin(["PVOD"] + cohorts)]
        is_case = (subset.cohort == "PVOD").to_numpy()
        area, ci = model.auc_ci_asymptotic(subset.gdf15.to_numpy(), is_case)
        row = {"comparison": f"PVOD vs {name}", "n_comparator": int((~is_case).sum()),
               "auc": 100 * area, "ci_low": 100 * ci[0], "ci_high": 100 * ci[1],
               "method": "C statistic, asymptotic 95% CI",
               # the asymptotic interval is the one that pairs with the published values,
               # which is why these rows come first; the reported interval is the bootstrap
               "reported_primary": False}
        if name in PUBLISHED_AUC:
            published = PUBLISHED_AUC[name]
            row.update({"published_auc": published[0], "published_ci_low": published[1],
                        "published_ci_high": published[2], "published_source": published[3],
                        "published_provenance": published[4],
                        "auc_difference": 100 * area - published[0],
                        "ci_low_difference": 100 * ci[0] - published[1]})
        rows.append(row)

    main = gdf15[gdf15.cohort.isin(["PVOD"] + AUC_COHORTS["IPAH+HPAH"])]
    is_case = (main.cohort == "PVOD").to_numpy()
    low, high = auc_ci_bootstrap(main.gdf15.to_numpy(), is_case)
    rows.append({"comparison": "PVOD vs IPAH+HPAH", "n_comparator": int((~is_case).sum()),
                 "auc": 100 * model.auc(main.gdf15.to_numpy(), is_case),
                 "ci_low": 100 * low, "ci_high": 100 * high,
                 "method": f"stratified percentile bootstrap, {model.N_BOOT} replicates, "
                           f"seed {model.SEED}",
                 "reported_primary": True})

    # how far the bootstrap interval moves with the random seed alone
    endpoints = [auc_ci_bootstrap(main.gdf15.to_numpy(), is_case, seed=seed)
                 for seed in range(n_seeds)]
    seed_sensitivity = {
        "n_seeds": n_seeds,
        "ci_low_range": [100 * min(e[0] for e in endpoints), 100 * max(e[0] for e in endpoints)],
        "ci_high_range": [100 * min(e[1] for e in endpoints), 100 * max(e[1] for e in endpoints)],
    }
    return rows, seed_sensitivity


def marker_row(subset, key, label, column):
    """One marker's discrimination of PVOD from IPAH+HPAH, on the patients given."""
    is_case = (subset.group == "PVOD").to_numpy()
    area = model.auc(subset[column], is_case)
    # each marker is oriented to whichever direction discriminates better
    oriented = max(area, 1 - area)
    values = subset[column] if area >= 0.5 else -subset[column]
    _, ci = model.auc_ci_asymptotic(values, is_case)
    return {"marker": label, "key": key, "n": len(subset),
            "n_pvod": int(is_case.sum()), "auc": oriented,
            "auc_ci_low": ci[0], "auc_ci_high": ci[1]}


def marker_table(clinical):
    """Discrimination of PVOD from IPAH+HPAH for every marker measured in the cohort.

    Each marker uses every patient who has it, so the denominators differ by marker.
    """
    return [marker_row(clinical.dropna(subset=[column]), key, label, column)
            for key, label, column in MARKERS]


def complete_case_marker_table(clinical):
    """The clinical comparators on the patients complete on all five measurements.

    The per-marker table above lets each marker use its own denominator, which runs from
    63 to 166, so no two rows describe the same patients.  This restricts to the patients
    who have all five clinical measurements and recomputes every AUC on that one set.  No
    value is imputed and nothing is modelled about why a value is missing.

    What the restriction buys is a common denominator.  What it costs is representativeness,
    and the cost is the larger of the two: whether a patient has a measurement depends on
    their diagnosis, most strongly for the variables missing most often, so the complete
    cases are roughly three times as enriched for PVOD as the full cohort.  The subset is
    therefore not a smaller version of the study and cannot confirm the full-cohort
    ranking; see ``clinical_completeness_by_diagnosis`` for the mechanism.

    NT-proBNP here is the distribution sheet, the source used throughout the analysis.
    The completeness table counts a patient as measured if either source records them,
    which is the right question for availability and a different one from this.
    """
    complete = clinical.dropna(subset=COMPLETE_CASE_COLUMNS)
    is_case = (complete.group == "PVOD").to_numpy()
    full_case_fraction = float((clinical.group == "PVOD").mean())
    patient_ids = sorted(int(value) for value in complete.patient_id)
    return {
        "n": len(complete),
        "n_pvod": int(is_case.sum()),
        "n_by_group": {group: int(count) for group, count
                       in complete.group.value_counts().items()},
        "pvod_fraction": float(is_case.mean()),
        "full_cohort_pvod_fraction": full_case_fraction,
        "case_enrichment": float(is_case.mean() / full_case_fraction),
        "required_columns": list(COMPLETE_CASE_COLUMNS),
        "ntprobnp_source": "NT-proBNP distribution sheet",
        "patient_ids": patient_ids,
        "patient_ids_by_group": {
            group: sorted(int(value) for value
                          in complete.patient_id[complete.group == group])
            for group in sorted(complete.group.unique())},
        "markers": [marker_row(complete, key, label, column)
                    for key, label, column in CLINICAL_MARKERS],
        "ageadj_gdf15_separation": ageadj_separation(complete),
        "raw_gdf15": ("raw GDF15 is not in this table and cannot be. It is published in "
                      "sheet GDF15-1 as six anonymous cohort columns of unequal length, "
                      "with no patient identifier, no age and no row order shared with "
                      "the per-patient sheets, so there is no key to attach it to a "
                      "clinical row and none can be constructed without guessing which "
                      "patient goes with which value"),
        "note": ("the patients complete on all five clinical measurements are a small and "
                 "strongly case-enriched subset: PVOD is 43.5% of them against 15.7% of "
                 "the full cohort, about three times enriched. This is not a smaller "
                 "version of the study but a different group of patients, because whether "
                 "a patient has a measurement depends on their diagnosis and depends on it "
                 "most for the variables missing most often. A complete-case comparison "
                 "therefore cannot confirm the full-cohort comparator ranking, and the "
                 "rankings among the clinical variables do move substantially at this "
                 "sample size. These rows are not an unbiased re-estimate of any AUC: they "
                 "are conditional on the measurement being available, and availability is "
                 "associated with diagnosis"),
    }


def ageadj_separation(complete):
    """Where age-adjusted GDF15 puts the two groups on the complete-case subset.

    Reported as two ranges rather than an AUC.  The groups do not overlap here, so the AUC
    is exactly 1.00 -- which is what the arithmetic gives when two small groups happen not
    to overlap, and with 10 patients against 13 that is not rare.  Quoted as a figure of
    merit it invites the reader to conclude overfitting, so the ranges are given instead.
    """
    case = complete.ageadj_gdf15[complete.group == "PVOD"]
    other = complete.ageadj_gdf15[complete.group != "PVOD"]
    return {
        "pvod_n": int(len(case)),
        "pvod_min": float(case.min()),
        "pvod_max": float(case.max()),
        "comparator_n": int(len(other)),
        "comparator_min": float(other.min()),
        "comparator_max": float(other.max()),
        "gap": float(case.min() - other.max()),
        "complete_separation": bool(case.min() > other.max()),
        "note": ("age-adjusted GDF15 separates the two groups completely on this subset: "
                 "every PVOD patient sits above every IPAH and HPAH patient, with a clean "
                 "gap between them. This is reported as complete separation and as the two "
                 "ranges, deliberately not as an AUC. The AUC would be 1.00, which measures "
                 "the absence of overlap in 10 cases against 13 comparators rather than the "
                 "performance of the marker, and beside the full-cohort 0.956 of Table 1 it "
                 "would read as the marker improving on the subset when what changed is "
                 "both the patients and the variable"),
    }


def completeness_by_diagnosis(patients):
    """Which clinical measurements exist, per diagnosis, from the per-patient reconstruction.

    This is the reason the complete-case subset is not representative: availability is
    associated with diagnosis, most strongly for cardiac index and NT-proBNP, the two
    variables missing most often. PVOD patients are investigated more thoroughly and so
    arrive with fuller data.

    Only presence is read, never the value.  The reconstruction recovers raw numbers by
    adding a median back to a published deviation, which leaves two NT-proBNP cells
    negative; nothing here interprets these cells as concentrations and no AUC uses them.

    NT-proBNP counts a patient as measured if either source sheet records them, which is
    17 of 26 PVOD. The complete-case subset uses the distribution sheet alone, 15 of 26.
    The two answer different questions -- was the patient measured at all, against which
    patients share one consistent analysable set -- and both are reported as such.
    """
    totals = patients.Group.value_counts().to_dict()
    rows = []
    for column, label in data.RECONSTRUCTED_MEASUREMENTS.items():
        present = patients.groupby("Group")[column].count().to_dict()
        row = {"variable": label, "column": column,
               "n_present": int(patients[column].notna().sum()),
               "n_total": int(len(patients))}
        for group in data.GROUP1:
            row[f"{group}_present"] = int(present.get(group, 0))
            row[f"{group}_total"] = int(totals[group])
            row[f"{group}_percent"] = 100.0 * present.get(group, 0) / totals[group]
        rows.append(row)
    return rows


def main():
    RESULTS.mkdir(exist_ok=True)
    gdf15 = data.load_gdf15()
    clinical = data.load_clinical()
    patients = data.load_patients_166()

    pvod = gdf15.loc[gdf15.cohort == "PVOD", "gdf15"]
    out = {
        "source_workbook": data.WORKBOOK.name,
        "patient_identifier": {
            "source": data.RECONSTRUCTED.name,
            "column": "Patient_ID",
            "range": "1-166, unique",
            "joined_on": ("age-adjusted GDF15, which is distinct across all 166 patients "
                          "in both files; the join is one-to-one and exact, and the "
                          "diagnosis label agrees on every row"),
            "provenance": ("Patient_ID was created for this per-patient reconstruction. It "
                           "is not a clinical identifier from the contributing centres, and "
                           "anything reaching the supplement that uses it needs a line "
                           "saying so"),
            "caution": ("the ids are blocked by diagnosis (1-26 PVOD, 27-124 IPAH, 125-166 "
                        "HPAH), so an id implies a group; ordering a figure by Patient_ID "
                        "would manufacture a pattern that is only that blocking"),
            "not_available_for": ("raw GDF15. Sheet GDF15-1 is anonymous group-level "
                                  "concentrations with no patient identifier, so no output "
                                  "derived from it carries a Patient_ID"),
        },
        "seed": model.SEED,
        "n_bootstrap": model.N_BOOT,
        "cohort_counts": gdf15.cohort.value_counts().to_dict(),
        "models": {},
    }

    # ---- main and comparator models -------------------------------------------------
    fitted_models, draws = {}, {}
    for comparator in COMPARATORS:
        fitted = build_model(gdf15, comparator)
        replicates = model.parametric_bootstrap(fitted)
        fitted_models[comparator], draws[comparator] = fitted, replicates
        comparator_values = gdf15.loc[gdf15.cohort.isin(COMPARATORS[comparator]), "gdf15"]
        coef_ci = fitted.coef_ci()
        entry = {
            "comparator": comparator,
            "n": fitted.n,
            "n_pvod": fitted.n_events,
            "n_comparator": fitted.n - fitted.n_events,
            "sampling_prevalence": fitted.prevalence,
            "sampling_log_odds": fitted.sampling_log_odds,
            "intercept": fitted.intercept,
            "intercept_se": fitted.se_intercept,
            "intercept_ci": list(coef_ci["intercept"]),
            "slope_log10_gdf15": fitted.slope,
            "slope_se": fitted.se_slope,
            "slope_ci": list(coef_ci["slope"]),
            "n_predictors": 1,
            "events_per_variable": fitted.n_events,
            "auc": model.auc(fitted.gdf15, fitted.is_case),
            "bootstrap_replicates_kept": replicates["n_kept"],
            "bootstrap_replicates_discarded": replicates["n_discarded"],
            "gdf15_at_lr": {},
            "comparator_max_gdf15": float(comparator_values.max()),
        }
        for target in REFERENCE_LR:
            low, high = model.gdf15_at_lr_ci(fitted, target, replicates)
            entry["gdf15_at_lr"][str(target)] = {
                "gdf15": float(fitted.gdf15_at_lr(target)),
                "ci_low": low, "ci_high": high,
            }
        out["models"][comparator] = entry

    main_model, main_draws = fitted_models[MAIN], draws[MAIN]
    main_entry = out["models"][MAIN]

    # ---- likelihood ratio at reference concentrations --------------------------------
    lr_rows = []
    for value in REFERENCE_GDF15:
        low, high = model.lr_curve_ci(main_model, [value], main_draws)
        lr_rows.append({
            "gdf15": value,
            "lr": float(main_model.lr(value)),
            "lr_ci_low": float(low[0]), "lr_ci_high": float(high[0]),
            **{f"post_test_from_{pre:.0%}".replace("%", "pct"): model.post_test(
                pre, float(main_model.lr(value)))["post_test_probability"]
               for pre in PRE_TEST},
        })
    main_entry["lr_at_reference_gdf15"] = lr_rows
    main_entry["worked_example"] = {
        str(value): {f"{pre:.0%}": model.post_test(pre, float(main_model.lr(value)))
                     for pre in PRE_TEST}
        for value in (1200, 1700)
    }

    # ---- observed support and the extrapolation boundary ------------------------------
    comparator_values = gdf15.loc[gdf15.cohort.isin(COMPARATORS[MAIN]), "gdf15"]
    n_above = int((pvod > comparator_values.max()).sum())
    main_entry["support"] = {
        "pvod_min_gdf15": float(pvod.min()),
        "pvod_max_gdf15": float(pvod.max()),
        "comparator_min_gdf15": float(comparator_values.min()),
        "comparator_max_gdf15": float(comparator_values.max()),
        "pvod_above_comparator_max": n_above,
        "pvod_above_comparator_max_fraction": n_above / len(pvod),
        "note": ("below the lowest observed PVOD value the curve is the fitted model "
                 "form rather than observed PVOD data; above the highest comparator "
                 "value no comparator patient was observed"),
    }

    # ---- discrimination --------------------------------------------------------------
    # The reported interval is the stratified percentile bootstrap.  It is bounded by
    # construction -- every replicate's AUC is itself an AUC, so no endpoint can fall
    # outside [0, 1] -- and it makes no distributional assumption about a statistic
    # estimated from 26 cases.  The asymptotic interval is kept beside it, and the source
    # study's published interval is quoted separately and never merged with either.
    _, asymptotic = model.auc_ci_asymptotic(main_model.gdf15, main_model.is_case)
    main_entry["auc_ci"] = list(auc_ci_bootstrap(main_model.gdf15, main_model.is_case))
    main_entry["auc_ci_method"] = (
        f"stratified percentile bootstrap, {model.N_BOOT:,} replicates, seed "
        f"{model.SEED}; PVOD and comparator patients resampled separately at their "
        "observed sample sizes")
    main_entry["auc_ci_asymptotic_secondary"] = list(asymptotic)
    main_entry["auc_ci_asymptotic_method"] = (
        "asymptotic C-statistic interval from the structural-component variance; "
        "secondary output, not combined with the bootstrap interval")
    main_entry["auc_ci_published_panel_e"] = [PUBLISHED_AUC[MAIN][1] / 100,
                                              PUBLISHED_AUC[MAIN][2] / 100]
    main_entry["auc_ci_published_panel_e_source"] = (
        f"{PUBLISHED_AUC[MAIN][3]} of the source study, as reported by a co-author; "
        "quoted for comparison and not computed here")
    out["auc_ci_provenance"] = auc_ci_provenance(asymptotic)

    # ---- functional form of the predictor ---------------------------------------------
    main_entry["functional_form"] = model.functional_form_test(main_model.gdf15,
                                                               main_model.is_case)

    # ---- stability of the curve against a selected threshold --------------------------
    # Both quantities come out of one stratified case bootstrap, recomputed inside every
    # replicate, so the selected cutoff and the curve landmark are compared like for like.
    stability = model.stability_bootstrap(main_model)
    # A percentile bootstrap on 26 cases moves between seeds, so how far is recorded rather
    # than left for a reader to discover; the same is done for the AUC below.
    across_seeds = [model.stability_bootstrap(main_model, seed=seed)
                    for seed in range(N_SEEDS)]
    parametric_lr10 = main_entry["gdf15_at_lr"]["10"]
    main_entry["threshold_comparison"] = {
        "method": stability["method"],
        "n_replicates": stability["n_replicates"],
        "n_replicates_dropped": stability["n_dropped"],
        "youden_cutoff": stability["youden"]["point"],
        "youden_ci_low": stability["youden"]["ci_low"],
        "youden_ci_high": stability["youden"]["ci_high"],
        "youden_ci_width": stability["youden"]["ci_width"],
        "lr10_gdf15": stability["landmark"]["point"],
        "lr10_ci_low": stability["landmark"]["ci_low"],
        "lr10_ci_high": stability["landmark"]["ci_high"],
        "lr10_ci_width": stability["landmark"]["ci_width"],
        "width_ratio": stability["width_ratio"],
        "seed_sensitivity": {
            "n_seeds": N_SEEDS,
            "lr10_ci_low_range": [min(r["landmark"]["ci_low"] for r in across_seeds),
                                  max(r["landmark"]["ci_low"] for r in across_seeds)],
            "lr10_ci_high_range": [min(r["landmark"]["ci_high"] for r in across_seeds),
                                   max(r["landmark"]["ci_high"] for r in across_seeds)],
            "youden_ci_low_range": [min(r["youden"]["ci_low"] for r in across_seeds),
                                    max(r["youden"]["ci_low"] for r in across_seeds)],
            "youden_ci_high_range": [min(r["youden"]["ci_high"] for r in across_seeds),
                                     max(r["youden"]["ci_high"] for r in across_seeds)],
            "width_ratio_range": [min(r["width_ratio"] for r in across_seeds),
                                  max(r["width_ratio"] for r in across_seeds)],
            "note": ("the LR = 10 endpoints move by a few pg/mL between seeds while the "
                     "Youden endpoints do not move at all, because a selected cutoff is a "
                     "discrete choice among observed values and the percentile lands on the "
                     "same two of them in every seed"),
        },
        "published_cutoff": PUBLISHED_CUTOFF,
        "lr_at_published_cutoff": float(main_model.lr(PUBLISHED_CUTOFF)),
        "lr10_ci_parametric_supplementary": [parametric_lr10["ci_low"],
                                             parametric_lr10["ci_high"]],
        "lr10_ci_width_parametric_supplementary": (parametric_lr10["ci_high"]
                                                   - parametric_lr10["ci_low"]),
        "note": ("the stratified case bootstrap is the reported uncertainty for the LR = 10 "
                 "landmark, because it is the scheme under which the landmark and a selected "
                 "threshold can be compared on identical footing; the parametric-bootstrap "
                 "interval is retained as a supplementary sensitivity analysis, and the "
                 "pointwise band on the fitted curve itself remains parametric"),
    }
    # the landmark's reported interval is the stratified one; the parametric is kept beside it
    main_entry["gdf15_at_lr"]["10"] = {
        "gdf15": stability["landmark"]["point"],
        "ci_low": stability["landmark"]["ci_low"],
        "ci_high": stability["landmark"]["ci_high"],
        "ci_method": stability["method"],
        "ci_low_parametric": parametric_lr10["ci_low"],
        "ci_high_parametric": parametric_lr10["ci_high"],
        "ci_method_parametric": ("parametric bootstrap, supplementary sensitivity analysis"),
    }

    # ---- internal validation ----------------------------------------------------------
    main_entry["internal_validation"] = model.bootstrap_optimism(main_model)

    positive = main_model.gdf15 >= PUBLISHED_CUTOFF
    cases = main_model.is_case == 1
    sensitivity = int(positive[cases].sum()), int(cases.sum())
    specificity = int((~positive[~cases]).sum()), int((~cases).sum())
    main_entry["published_cutoff_performance"] = {
        "cutoff": PUBLISHED_CUTOFF,
        "sensitivity": sensitivity[0] / sensitivity[1],
        "sensitivity_ci": list(model.clopper_pearson(*sensitivity)),
        "sensitivity_counts": sensitivity,
        "specificity": specificity[0] / specificity[1],
        "specificity_ci": list(model.clopper_pearson(*specificity)),
        "specificity_counts": specificity,
    }

    # ---- leave one PVOD case out --------------------------------------------------------
    # How far one patient moves the answer depends entirely on where on the curve the
    # answer is read, so both locations are named and fixed across all 26 refits: the
    # published cutoff, and the highest concentration anyone in this cohort was measured
    # at.  A percentage quoted without its concentration cannot be checked by a reader.
    top_of_range = float(main_model.gdf15.max())
    loo = model.leave_one_case_out(main_model, [PUBLISHED_CUTOFF, top_of_range])
    loo_rows = []
    for refit in loo["refits"]:
        row = {"case_id": refit["case_id"], "excluded_gdf15": refit["excluded_gdf15"],
               "n": refit["n"], "n_pvod": refit["n_events"],
               "intercept": refit["intercept"], "slope_log10_gdf15": refit["slope"]}
        for index, point in enumerate(loo["evaluation_points"]):
            row[f"lr_at_{point:.0f}"] = float(refit["lr"][index])
            row[f"percent_change_at_{point:.0f}"] = float(refit["percent_change"][index])
        loo_rows.append(row)
    main_entry["leave_one_case_out"] = {
        "n_refits": loo["n_refits"],
        "method": loo["method"],
        "case_id_definition": ("position 1-26 in the PVOD column of sheet GDF15-1, the "
                               "sheet the model reads. This is a position in a source "
                               "sheet, not a patient: GDF15-1 is anonymous group-level "
                               "concentrations, so no Patient_ID exists for these cases "
                               "and none is reported. Each case is quoted with its own "
                               "GDF15 value, which distinguishes it within this table "
                               "only. A case_id here has no relation to the Patient_ID "
                               "used everywhere else, and the two must never be equated"),
        "patient_id": ("not available and not constructed. Attaching one would mean "
                       "guessing which patient goes with which raw value, which this "
                       "project refuses to do here as everywhere else"),
        "top_of_observed_range": top_of_range,
        "at_evaluation_point": loo["at_evaluation_point"],
        "note": ("a sensitivity analysis, not a validation procedure: the reported model "
                 "remains the fit on all 26 cases and nothing here corrects it"),
    }

    # ---- the published cutoff, read two ways ---------------------------------------------
    # The same patient at the same concentration, once as a point on the fitted curve and
    # once as a positive result from the published dichotomy, combined with the same
    # pre-test probabilities so that the reading is the only thing that differs.
    cutoff = main_entry["published_cutoff_performance"]
    lr_at_cutoff = float(main_model.lr(PUBLISHED_CUTOFF))
    band_low, band_high = model.lr_curve_ci(main_model, [PUBLISHED_CUTOFF], main_draws)
    positive_lr = cutoff["sensitivity"] / (1 - cutoff["specificity"])
    main_entry["published_cutoff_example"] = {
        "gdf15": PUBLISHED_CUTOFF,
        "continuous_lr": lr_at_cutoff,
        "continuous_lr_ci": [float(band_low[0]), float(band_high[0])],
        "continuous_lr_ci_method": (f"parametric bootstrap, {model.N_BOOT:,} replicates, "
                                    f"seed {model.SEED}"),
        "cutoff_positive_lr": positive_lr,
        "cutoff_positive_lr_derivation": (
            f"sensitivity / (1 - specificity) = "
            f"({cutoff['sensitivity_counts'][0]}/{cutoff['sensitivity_counts'][1]}) / "
            f"({cutoff['specificity_counts'][1] - cutoff['specificity_counts'][0]}/"
            f"{cutoff['specificity_counts'][1]})"),
        "by_pre_test": [
            {"pre_test_probability": pre,
             "continuous_post_test": model.post_test(
                 pre, lr_at_cutoff)["post_test_probability"],
             "continuous_post_test_ci": [
                 model.post_test(pre, float(band_low[0]))["post_test_probability"],
                 model.post_test(pre, float(band_high[0]))["post_test_probability"]],
             "cutoff_post_test": model.post_test(
                 pre, positive_lr)["post_test_probability"]}
            for pre in PRE_TEST],
        "note": ("the dichotomy assigns every value above the cutoff the same weight, so "
                 "a patient measured at the cutoff itself is handed the average evidence "
                 "of everyone above it. The curve reads that patient's own value, which "
                 "at the cutoff carries far less weight. Both readings use the same "
                 "patients and the same outcome; they differ only in what is done with "
                 "the measurement"),
    }

    # ---- calibration ------------------------------------------------------------------
    calibration = model.calibration_by_quintile(main_model)
    main_entry["calibration_quintiles"] = calibration

    # ---- healthy controls --------------------------------------------------------------
    controls = gdf15.loc[gdf15.cohort == "Ctrl", "gdf15"]
    out["healthy_controls"] = {
        "n_controls": len(controls),
        "n_pvod": len(pvod),
        "control_max_gdf15": float(controls.max()),
        "pvod_min_gdf15": float(pvod.min()),
        "perfectly_separated": bool(controls.max() < pvod.min()),
        "gap_pg_ml": float(pvod.min() - controls.max()),
        "note": ("GDF15 separates PVOD from healthy controls completely, so the "
                 "maximum-likelihood estimate of the log10(GDF15) coefficient is "
                 "unbounded and no finite curve can be fitted"),
    }

    # ---- every AUC in circulation for these data ------------------------------------------
    reconciliation, seed_sensitivity = auc_reconciliation(gdf15)
    out["auc_reconciliation"] = reconciliation
    out["auc_bootstrap_seed_sensitivity"] = seed_sensitivity

    out["source_figure_verification"] = {
        "status": "pending local verification",
        "panel_d": {
            "provenance": PANEL_D,
            "values": ["PVOD vs IPAH 94.94", "PVOD vs Other PH 97.26", "PVOD vs COPD 95.19"],
            "reproduced": "all three reproduce from the supporting data to two decimal places",
        },
        "panel_e": {
            "provenance": PANEL_E,
            "values": ["PVOD vs IPAH+HPAH 95.60 (91.64-99.57)",
                       "PVOD vs IPAH+HPAH+Other PH 96.05 (92.45-99.66)"],
            "reproduced": ("both reproduce from the supporting data: 95.60 and 96.05 by "
                           "C statistic"),
            "note": ("these panel E values are those a co-author of the source study reports "
                     "from the 600-dpi source render. This analysis has not inspected that "
                     "render and records them as reported, pending local verification. No "
                     "claim is made here about any other value for panel E"),
        },
    }

    # ---- NT-proBNP: the two source sheets and their denominators ---------------------------
    ntprobnp = data.load_ntprobnp_correlation()
    extra = ntprobnp[ntprobnp.in_correlation_sheet & ~ntprobnp.in_distribution_sheet]
    out["ntprobnp_denominator"] = {
        "primary_source": "NT-proBNP",
        "primary_n": int(clinical.ntprobnp.notna().sum()),
        "primary_n_by_group": {group: int(count) for group, count
                               in clinical.groupby("group").ntprobnp.count().items()},
        "correlation_source": "6MWD PVR mPAP NT-proBNP",
        "correlation_n_pvod": int(ntprobnp.in_correlation_sheet.sum()),
        "pvod_ids_in_correlation_sheet_only": [int(i) for i in extra.pvod_id],
        "pvod_values_in_correlation_sheet_only": [float(v) for v
                                                  in extra.ntprobnp_correlation_sheet],
        "median_offset_ng_l": ntprobnp.attrs["median_offset"],
        "note": ("the analysis uses the NT-proBNP distribution sheet, which records 63 "
                 "measurements across the 166 group-1 patients (PVOD 15, IPAH 33, HPAH 15). "
                 "The PVOD-only correlation sheet records 17 NT-proBNP measurements for the "
                 "same 26 PVOD patients: the 15 of the distribution sheet, which agree "
                 "exactly once the distribution sheet's median centring is added back, and "
                 "two further patients, PVOD ids 24 and 26, which the distribution sheet "
                 "records as NA. The two sheets therefore support different denominators "
                 "for PVOD. This analysis states both and uses the distribution sheet, whose "
                 "denominator matches the one published for this variable"),
    }

    # ---- other markers measured in the same cohort ---------------------------------------
    out["marker_auc"] = marker_table(clinical)
    out["marker_auc_complete_case"] = complete_case_marker_table(clinical)
    out["clinical_completeness"] = {
        column: {"n": int(clinical[column].notna().sum()),
                 **{group: int(count) for group, count
                    in clinical.groupby("group")[column].count().items()}}
        for column in ["ntprobnp", "pvr", "mpap", "ci", "sixmwd_m"]
    }
    out["clinical_completeness_by_diagnosis"] = {
        "source": data.RECONSTRUCTED.name,
        "rows": completeness_by_diagnosis(patients),
        "ntprobnp_source": ("either source sheet counts as measured here, 17 of 26 PVOD. "
                            "The complete-case subset uses the distribution sheet alone, "
                            "15 of 26. The two answer different questions -- whether the "
                            "patient was measured at all, against which patients share one "
                            "consistent analysable set -- and are reported separately"),
        "values_not_used": ("only the presence or absence of a measurement is read from "
                            "the reconstruction. Its numeric values are recovered by adding "
                            "a median back to a published deviation, which leaves NT-proBNP "
                            "negative for Patient_ID 38 and 58; they are not interpreted as "
                            "concentrations and no AUC is estimated from them. Every AUC in "
                            "this project is computed from the original workbook"),
        "note": ("whether a patient has a measurement depends on their diagnosis, and it "
                 "depends on it most for cardiac index and NT-proBNP, the two variables "
                 "missing most often. PVOD patients are investigated more thoroughly and "
                 "arrive with fuller data. This is why the complete-case subset is "
                 "case-enriched and why comparator AUCs on it are conditional on the "
                 "measurement being available"),
    }

    # ---- write ---------------------------------------------------------------------------
    (RESULTS / "results.json").write_text(json.dumps(out, indent=2, default=float))

    pd.DataFrame([
        {"model": name, "n": entry["n"], "n_pvod": entry["n_pvod"],
         "sampling_prevalence": entry["sampling_prevalence"],
         "intercept": entry["intercept"], "intercept_se": entry["intercept_se"],
         "slope_log10_gdf15": entry["slope_log10_gdf15"], "slope_se": entry["slope_se"],
         "slope_ci_low": entry["slope_ci"][0], "slope_ci_high": entry["slope_ci"][1],
         "events_per_variable": entry["events_per_variable"], "auc": entry["auc"]}
        for name, entry in out["models"].items()
    ]).to_csv(RESULTS / "model_coefficients.csv", index=False)

    pd.DataFrame([
        {"comparator": name, "n": entry["n"], "n_pvod": entry["n_pvod"],
         "gdf15_at_lr10": entry["gdf15_at_lr"]["10"]["gdf15"],
         "ci_low": entry["gdf15_at_lr"]["10"]["ci_low"],
         "ci_high": entry["gdf15_at_lr"]["10"]["ci_high"],
         "ci_method": entry["gdf15_at_lr"]["10"].get(
             "ci_method", "parametric bootstrap"),
         "comparator_max_gdf15": entry["comparator_max_gdf15"]}
        for name, entry in out["models"].items()
    ]).to_csv(RESULTS / "comparator_models.csv", index=False)

    curves = []
    for name, fitted in fitted_models.items():
        low, high = model.lr_curve_ci(fitted, CURVE_GRID, draws[name])
        curves.append(pd.DataFrame({"model": name, "gdf15": CURVE_GRID,
                                    "lr": fitted.lr(CURVE_GRID),
                                    "lr_ci_low": low, "lr_ci_high": high}))
    pd.concat(curves).to_csv(RESULTS / "lr_curves.csv", index=False)

    pd.DataFrame(loo_rows).to_csv(RESULTS / "leave_one_case_out.csv", index=False)

    functional = main_entry["functional_form"]
    pd.DataFrame([
        {"model": "log10(GDF15)", "parameters": 2,
         "loglikelihood": functional["linear_loglikelihood"],
         "aic": functional["linear_aic"]},
        {"model": "log10(GDF15) + log10(GDF15)^2", "parameters": 3,
         "loglikelihood": functional["quadratic_loglikelihood"],
         "aic": functional["quadratic_aic"]},
    ]).assign(lr_statistic=functional["lr_statistic"], df=functional["df"],
              p_value=functional["p_value"], preferred=functional["preferred"]
              ).to_csv(RESULTS / "functional_form.csv", index=False)

    complete_case = out["marker_auc_complete_case"]
    pd.DataFrame(complete_case["markers"]).assign(
        common_n=complete_case["n"], n_pvod=complete_case["n_pvod"]
        ).to_csv(RESULTS / "marker_auc_complete_case.csv", index=False)

    # the exact subset, by identifier, so anyone can reproduce it without re-deriving it
    complete = clinical[clinical.patient_id.isin(complete_case["patient_ids"])]
    complete[["patient_id", "group"]].sort_values("patient_id").to_csv(
        RESULTS / "complete_case_patients.csv", index=False)

    pd.DataFrame(out["clinical_completeness_by_diagnosis"]["rows"]).to_csv(
        RESULTS / "clinical_completeness.csv", index=False)

    pd.DataFrame(calibration).to_csv(RESULTS / "calibration.csv", index=False)
    pd.DataFrame(out["marker_auc"]).to_csv(RESULTS / "marker_auc.csv", index=False)
    pd.DataFrame(reconciliation).to_csv(RESULTS / "auc_reconciliation.csv", index=False)
    pd.DataFrame(lr_rows).to_csv(RESULTS / "lr_reference.csv", index=False)

    validation = main_entry["internal_validation"]
    pd.DataFrame([
        {"quantity": "C statistic (AUC)",
         "apparent": validation["apparent_auc"],
         "optimism": validation["auc_optimism"],
         "optimism_corrected": validation["corrected_auc"]},
        {"quantity": "calibration slope",
         "apparent": validation["apparent_calibration_slope"],
         "optimism": validation["calibration_slope_optimism"],
         "optimism_corrected": validation["corrected_calibration_slope"]},
    ]).assign(method=validation["method"],
              replicates=validation["n_replicates"],
              replicates_dropped=validation["n_dropped"]
              ).to_csv(RESULTS / "internal_validation.csv", index=False)

    ntprobnp.to_csv(RESULTS / "ntprobnp_sources.csv", index=False)

    print(f"n = {main_entry['n']}, PVOD = {main_entry['n_pvod']}, "
          f"prevalence = {main_entry['sampling_prevalence']:.3f}")
    print(f"intercept {main_entry['intercept']:.3f}, "
          f"log10(GDF15) coefficient {main_entry['slope_log10_gdf15']:.3f} "
          f"(SE {main_entry['slope_se']:.3f})")
    landmark = main_entry["gdf15_at_lr"]["10"]
    threshold = main_entry["threshold_comparison"]
    print(f"LR = 10 at {landmark['gdf15']:.0f} pg/mL "
          f"(95% CI {landmark['ci_low']:.0f}-{landmark['ci_high']:.0f}, "
          f"stratified case bootstrap)")
    print(f"Youden cutoff {threshold['youden_cutoff']:.0f} pg/mL "
          f"(95% CI {threshold['youden_ci_low']:.0f}-{threshold['youden_ci_high']:.0f}); "
          f"width ratio {threshold['width_ratio']:.2f}")
    print(f"AUC {main_entry['auc']:.4f} "
          f"(95% CI {main_entry['auc_ci'][0]:.3f}-{main_entry['auc_ci'][1]:.3f}, "
          f"stratified percentile bootstrap)")
    print(f"AUC apparent {validation['apparent_auc']:.4f}, "
          f"optimism-corrected {validation['corrected_auc']:.4f}; "
          f"calibration slope apparent 1.0000, "
          f"optimism-corrected {validation['corrected_calibration_slope']:.4f}")
    functional = main_entry["functional_form"]
    print(f"quadratic term: LR statistic {functional['lr_statistic']:.3f}, "
          f"p = {functional['p_value']:.3f}; AIC {functional['linear_aic']:.2f} linear "
          f"against {functional['quadratic_aic']:.2f} quadratic -- "
          f"{functional['preferred']} model kept")
    for point in main_entry["leave_one_case_out"]["at_evaluation_point"]:
        print(f"leave-one-case-out at {point['gdf15']:,.0f} pg/mL: LR "
              f"{point['full_cohort_lr']:,.3f}, largest change "
              f"{point['signed_percent_change']:+.1f}% on removing PVOD case "
              f"{point['case_id']} ({point['excluded_gdf15']:,.0f} pg/mL)")
    example = main_entry["published_cutoff_example"]
    for entry in example["by_pre_test"]:
        print(f"at {example['gdf15']:,.0f} pg/mL, {entry['pre_test_probability']:.0%} "
              f"pre-test: curve {entry['continuous_post_test']:.1%} "
              f"({entry['continuous_post_test_ci'][0]:.1%}-"
              f"{entry['continuous_post_test_ci'][1]:.1%}), "
              f"cutoff {entry['cutoff_post_test']:.1%}")
    complete_case = out["marker_auc_complete_case"]
    print(f"complete on all {len(COMPLETE_CASE_COLUMNS)} clinical measurements: "
          f"n = {complete_case['n']} "
          f"({', '.join(f'{g} {c}' for g, c in sorted(complete_case['n_by_group'].items()))}), "
          f"PVOD {complete_case['pvod_fraction']:.1%} against "
          f"{complete_case['full_cohort_pvod_fraction']:.1%} in the full cohort "
          f"({complete_case['case_enrichment']:.1f}x enriched)")
    print(f"  Patient_IDs: {', '.join(str(i) for i in complete_case['patient_ids'])}")
    separation = complete_case["ageadj_gdf15_separation"]
    print(f"  age-adjusted GDF15 separates completely: PVOD "
          f"{separation['pvod_min']:.2f}-{separation['pvod_max']:.2f} against IPAH/HPAH "
          f"{separation['comparator_min']:.2f}-{separation['comparator_max']:.2f} "
          f"(reported as ranges, not as an AUC)")
    for row in out["clinical_completeness_by_diagnosis"]["rows"]:
        print(f"  {row['variable']:<14}"
              + "  ".join(f"{g} {row[f'{g}_percent']:.0f}%" for g in data.GROUP1))
    print(f"wrote {RESULTS}/results.json and 14 tables")
    make_figures.main()


if __name__ == "__main__":
    main()
