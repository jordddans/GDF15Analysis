"""Numerical checks on the implementation.

Each check derives its expected value independently -- by brute force, by numerical
root-finding, or from the definition of the quantity -- rather than comparing against a
stored number, so a check can fail on a real error rather than on a stale constant.

Usage:  python analysis/verify.py
"""

import contextlib
import hashlib
import io
import tempfile
from pathlib import Path

import numpy as np
import statsmodels.api as sm
from scipy.optimize import brentq
from scipy.stats import chi2

import data
import make_figures
import model
import run_analysis
import json
import reconstruct
import verify_centerpiece
import audit_numbers
import verify_figure1

CHECKS = []


def check(name, condition, detail=""):
    CHECKS.append({"name": name, "passed": bool(condition), "detail": detail})
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {name}{'  ' + detail if detail else ''}")
    if not condition:
        raise AssertionError(name)


def main():
    gdf15 = data.load_gdf15()
    subset = gdf15[gdf15.cohort.isin(["PVOD", "IPAH", "HPAH"])]
    values = subset.gdf15.to_numpy()
    is_case = (subset.cohort == "PVOD").to_numpy()
    fitted = model.fit_lr_model(values, is_case, "PVOD vs IPAH+HPAH")

    # 1. the likelihood ratio is the fitted odds divided by the sampling odds
    for gdf in (500.0, 1200.0, 1700.0, 3000.0):
        probability = fitted.probability(gdf)
        fitted_odds = probability / (1 - probability)
        sampling_odds = fitted.prevalence / (1 - fitted.prevalence)
        check(f"LR({gdf:,.0f}) equals fitted odds / sampling odds",
              np.isclose(fitted.lr(gdf), fitted_odds / sampling_odds, rtol=1e-12))

    # 2. the fitted probabilities follow the sample's case mix; the likelihood ratio does not.
    #    Duplicating every comparator patient leaves both conditional GDF15 distributions
    #    exactly as they were and halves the sampling odds.
    comparator = np.where(~is_case)[0]
    doubled = model.fit_lr_model(np.concatenate([values, values[comparator]]),
                                 np.concatenate([is_case, is_case[comparator]]), "doubled")
    sampling_odds_ratio = ((fitted.prevalence / (1 - fitted.prevalence))
                           / (doubled.prevalence / (1 - doubled.prevalence)))
    probability_odds_ratio = ((fitted.probability(1000) / (1 - fitted.probability(1000)))
                              / (doubled.probability(1000) / (1 - doubled.probability(1000))))
    lr_ratio = fitted.lr(1000) / doubled.lr(1000)
    check("doubling the comparator group halves the sampling odds",
          np.isclose(sampling_odds_ratio, 2.0, rtol=0.02),
          f"odds ratio {sampling_odds_ratio:.3f}")
    check("the fitted odds at 1,000 pg/mL move with the sampling odds",
          abs(np.log(probability_odds_ratio) - np.log(sampling_odds_ratio))
          < abs(np.log(probability_odds_ratio)),
          f"fitted odds ratio {probability_odds_ratio:.3f}")
    check("the likelihood ratio at 1,000 pg/mL does not",
          abs(np.log(lr_ratio)) < 0.5 * abs(np.log(probability_odds_ratio)),
          f"LR ratio {lr_ratio:.3f}")

    # 3. the pre-test -> post-test chain inverts back to the likelihood ratio
    chain = model.post_test(0.10, fitted.lr(1700))
    posterior = chain["post_test_probability"]
    recovered = (posterior / (1 - posterior)) / (0.10 / 0.90)
    check("post-test probability inverts to the likelihood ratio",
          np.isclose(recovered, fitted.lr(1700), rtol=1e-12),
          f"post-test {posterior:.4f} from a 10% pre-test probability")

    # 4. the closed-form solution for LR = 10 agrees with a numerical root
    root = brentq(lambda g: fitted.lr(g) - 10, 1000, 6000, xtol=1e-8)
    check("closed-form GDF15 at LR = 10 matches a numerical root",
          np.isclose(fitted.gdf15_at_lr(10), root, rtol=1e-8),
          f"{fitted.gdf15_at_lr(10):,.2f} pg/mL")

    # 5. the AUC and its asymptotic variance match direct computation from the definition
    cases, controls = values[is_case], values[~is_case]
    kernel = lambda a, b: (a > b) * 1.0 + (a == b) * 0.5
    case_component = np.array([kernel(c, controls).mean() for c in cases])
    control_component = np.array([kernel(cases, k).mean() for k in controls])
    variance = (case_component.var(ddof=1) / len(cases)
                + control_component.var(ddof=1) / len(controls))
    area, interval = model.auc_ci_asymptotic(values, is_case)
    half_width = (interval[1] - interval[0]) / 2
    check("AUC matches the direct Mann-Whitney computation",
          np.isclose(area, case_component.mean(), atol=1e-12), f"{area:.6f}")
    check("asymptotic variance matches the structural components computed directly",
          np.isclose((half_width / 1.959963984540054) ** 2, variance, rtol=1e-9),
          f"95% CI {interval[0]:.4f}-{interval[1]:.4f}")

    # 6. exact binomial limits bracket the observed proportion
    low, high = model.clopper_pearson(21, 26)
    check("exact binomial interval brackets the observed proportion",
          low < 21 / 26 < high, f"{low:.4f}-{high:.4f}")
    check("exact binomial interval for zero events starts at zero",
          model.clopper_pearson(0, 34)[0] == 0.0)

    # 7. the parametric bootstrap is reproducible from its seed
    first = model.parametric_bootstrap(fitted, n_boot=200, seed=model.SEED)
    second = model.parametric_bootstrap(fitted, n_boot=200, seed=model.SEED)
    check("the parametric bootstrap repeats exactly for a fixed seed",
          np.array_equal(first["slope"], second["slope"]))

    # 8. the optimism correction satisfies its own identity, and an independent
    #    reimplementation of the estimator from its definition reproduces it
    validation = model.bootstrap_optimism(fitted, n_boot=200, seed=model.SEED)
    check("apparent AUC in the optimism correction is the C statistic itself",
          np.isclose(validation["apparent_auc"], model.auc(values, is_case), rtol=1e-12),
          f"{validation['apparent_auc']:.6f}")
    check("corrected AUC = apparent - optimism",
          np.isclose(validation["corrected_auc"],
                     validation["apparent_auc"] - validation["auc_optimism"], rtol=1e-12),
          f"{validation['apparent_auc']:.4f} - ({validation['auc_optimism']:+.4f}) "
          f"= {validation['corrected_auc']:.4f}")
    check("corrected calibration slope = apparent - optimism",
          np.isclose(validation["corrected_calibration_slope"],
                     validation["apparent_calibration_slope"]
                     - validation["calibration_slope_optimism"], rtol=1e-12),
          f"{validation['corrected_calibration_slope']:.4f}")

    # the same estimator, written out here from its definition against its own fits
    x = np.log10(values)
    design = sm.add_constant(x)
    outcome = is_case.astype(int)
    rng = np.random.default_rng(model.SEED)
    auc_optimism, slope_optimism = [], []
    for _ in range(200):
        take = rng.integers(0, len(outcome), len(outcome))
        resampled = outcome[take]
        if resampled.sum() < 2 or resampled.sum() > len(outcome) - 2:
            continue
        fit = sm.Logit(resampled, design[take]).fit(disp=0)
        if not fit.mle_retvals.get("converged", False) or not 0 < fit.params[1] < 100:
            continue
        in_replicate = fit.params[0] + fit.params[1] * x[take]
        in_original = fit.params[0] + fit.params[1] * x
        recalibrated = sm.Logit(outcome, sm.add_constant(in_original)).fit(disp=0)
        if not recalibrated.mle_retvals.get("converged", False):
            continue
        # a model's C statistic where it was fitted, minus where it was not
        auc_optimism.append(model.auc(in_replicate, resampled.astype(bool))
                            - model.auc(in_original, is_case))
        # the calibration slope is 1 in the sample a model was fitted to
        slope_optimism.append(1.0 - float(recalibrated.params[1]))
    check("optimism recomputed independently from the definition matches",
          np.isclose(validation["auc_optimism"], np.mean(auc_optimism), rtol=1e-10)
          and np.isclose(validation["calibration_slope_optimism"],
                         np.mean(slope_optimism), rtol=1e-10),
          f"AUC optimism {np.mean(auc_optimism):+.6f}, "
          f"slope optimism {np.mean(slope_optimism):+.6f}")

    # scoring a model on the very sample it was fitted to leaves no optimism to remove
    self_fit = sm.Logit(outcome, design).fit(disp=0)
    self_lp = self_fit.params[0] + self_fit.params[1] * x
    self_slope = sm.Logit(outcome, sm.add_constant(self_lp)).fit(disp=0).params[1]
    check("optimism of a model scored on its own fitting sample is zero",
          np.isclose(self_slope, 1.0, atol=1e-6),
          f"apparent calibration slope {self_slope:.8f}")

    # 9. the stability bootstrap computes both quantities inside one shared replicate.
    #    With a single replicate both percentile limits collapse onto that replicate's
    #    value, so recomputing the replicate here from its own resample proves that the
    #    cutoff and the curve landmark came from the same patients.
    single = model.stability_bootstrap(fitted, n_boot=1, seed=7)
    take = model.stratified_indices(is_case, np.random.default_rng(7))
    replicate_fit = sm.Logit(outcome[take], sm.add_constant(x[take])).fit(disp=0)
    offset = np.log(outcome[take].mean() / (1 - outcome[take].mean()))
    landmark = 10 ** ((np.log(10) + offset - replicate_fit.params[0]) / replicate_fit.params[1])
    cutoff = model.youden_cutoff(values[take], is_case[take])
    check("both stability quantities come from one shared stratified replicate",
          np.isclose(single["landmark"]["ci_low"], landmark, rtol=1e-9)
          and np.isclose(single["landmark"]["ci_high"], landmark, rtol=1e-9)
          and np.isclose(single["youden"]["ci_low"], cutoff, rtol=1e-12)
          and np.isclose(single["youden"]["ci_high"], cutoff, rtol=1e-12),
          f"LR = 10 at {landmark:,.1f}, Youden {cutoff:,.2f} in that replicate")
    check("the stability bootstrap repeats exactly for a fixed seed",
          model.stability_bootstrap(fitted, n_boot=100, seed=3)
          == model.stability_bootstrap(fitted, n_boot=100, seed=3))

    # 10. the two NT-proBNP source sheets, and what separates their denominators
    ntprobnp = data.load_ntprobnp_correlation()
    only_correlation = ntprobnp[ntprobnp.in_correlation_sheet
                                & ~ntprobnp.in_distribution_sheet]
    check("the NT-proBNP denominators differ by exactly two PVOD patients",
          ntprobnp.in_correlation_sheet.sum() - ntprobnp.in_distribution_sheet.sum() == 2
          and len(only_correlation) == 2,
          f"PVOD ids {sorted(only_correlation.pvod_id)} "
          f"({', '.join(f'{v:,.0f}' for v in only_correlation.ntprobnp_correlation_sheet)} ng/L)")
    check("no distribution-sheet value is absent from the correlation sheet",
          not (ntprobnp.in_distribution_sheet & ~ntprobnp.in_correlation_sheet).any())
    shared = ntprobnp[ntprobnp.in_correlation_sheet & ntprobnp.in_distribution_sheet]
    difference = (shared.ntprobnp_correlation_sheet - shared.ntprobnp_distribution_sheet
                  - ntprobnp.attrs["median_offset"]).abs().max()
    check("the sheets agree on every patient they share, to the recorded rounding",
          difference < 0.5,
          f"largest difference {difference:.1f} ng/L beyond the "
          f"{ntprobnp.attrs['median_offset']:,.0f} ng/L median centring")


    # 11. leave-one-case-out removes each PVOD case exactly once and nothing else.
    #     The 26 refits are identified here by which patient is missing from each, so a
    #     refit that silently dropped the wrong patient, or the same one twice, fails.
    top_of_range = float(values.max())
    loo = model.leave_one_case_out(fitted, [1658.0, top_of_range])
    excluded = [r["excluded_gdf15"] for r in loo["refits"]]
    case_values = sorted(values[is_case])
    check("leave-one-case-out refits the model once per PVOD case",
          loo["n_refits"] == int(is_case.sum()) == 26, f"{loo['n_refits']} refits")
    check("each refit drops exactly one PVOD case, and the 26 are all of them",
          sorted(excluded) == case_values,
          f"{len(set(excluded))} distinct excluded values")
    check("every refit keeps one fewer patient and one fewer case",
          all(r["n"] == fitted.n - 1 and r["n_events"] == fitted.n_events - 1
              for r in loo["refits"]),
          f"n {loo['refits'][0]['n']}, events {loo['refits'][0]['n_events']}")
    # the reported change is recomputed here from the two fits it claims to compare
    worst = loo["at_evaluation_point"][0]
    without = next(r for r in loo["refits"] if r["case_id"] == worst["case_id"])
    independent = 100 * (float(without["lr"][0]) - float(fitted.lr(1658.0))) / float(
        fitted.lr(1658.0))
    check("the reported percentage change is the change between those two fits",
          np.isclose(worst["signed_percent_change"], independent, rtol=1e-12),
          f"{worst['signed_percent_change']:+.2f}% at 1,658 pg/mL, "
          f"PVOD case {worst['case_id']}")
    check("a case's influence differs by where on the curve it is read",
          (loo["at_evaluation_point"][1]["max_absolute_percent_change"]
           > 5 * loo["at_evaluation_point"][0]["max_absolute_percent_change"]),
          f"{loo['at_evaluation_point'][0]['max_absolute_percent_change']:.1f}% at 1,658 "
          f"against {loo['at_evaluation_point'][1]['max_absolute_percent_change']:.1f}% "
          f"at {top_of_range:,.0f} pg/mL")

    # 12. the reported AUC interval is a bootstrap, so it repeats for its seed and cannot
    #     leave [0, 1]; the asymptotic interval beside it is formed and then truncated
    first = run_analysis.auc_ci_bootstrap(values, is_case)
    second = run_analysis.auc_ci_bootstrap(values, is_case)
    area = model.auc(values, is_case)
    check("the bootstrap AUC interval repeats exactly for a fixed seed",
          first == second, f"{first[0]:.4f}-{first[1]:.4f}")
    check("the bootstrap AUC interval brackets the point estimate and stays in [0, 1]",
          first[0] < area < first[1] and 0 <= first[0] and first[1] <= 1,
          f"AUC {area:.4f}")
    _, asymptotic = model.auc_ci_asymptotic(values, is_case)
    check("the asymptotic interval is the normal interval, clipped to [0, 1]",
          np.isclose(asymptotic[1], min(1.0, area + (asymptotic[1] - asymptotic[0]) / 2),
                     rtol=1e-12),
          f"{asymptotic[0]:.5f}-{asymptotic[1]:.5f}, no clipping applied here")

    # 13. the quadratic comparison is a likelihood-ratio test, recomputed here from the
    #     two fits and from the definitions of the statistic and of the AIC
    functional = model.functional_form_test(values, is_case)
    linear_fit = sm.Logit(is_case.astype(int), sm.add_constant(np.log10(values))).fit(disp=0)
    quadratic_fit = sm.Logit(
        is_case.astype(int),
        sm.add_constant(np.column_stack([np.log10(values), np.log10(values) ** 2]))
    ).fit(disp=0)
    statistic = 2 * (quadratic_fit.llf - linear_fit.llf)
    check("the quadratic LR statistic is twice the log-likelihood difference",
          np.isclose(functional["lr_statistic"], statistic, rtol=1e-10),
          f"{statistic:.4f} on 1 df")
    check("its p-value is the chi-square upper tail on one degree of freedom",
          np.isclose(functional["p_value"], chi2.sf(statistic, 1), rtol=1e-12),
          f"p = {functional['p_value']:.4f}")
    check("both AICs are -2 x log-likelihood + 2 x parameters",
          np.isclose(functional["linear_aic"], -2 * linear_fit.llf + 4, rtol=1e-12)
          and np.isclose(functional["quadratic_aic"], -2 * quadratic_fit.llf + 6,
                         rtol=1e-12),
          f"{functional['linear_aic']:.3f} linear, "
          f"{functional['quadratic_aic']:.3f} quadratic")
    check("the linear model is the one the AIC prefers",
          functional["preferred"] == "linear"
          and functional["linear_aic"] < functional["quadratic_aic"])

    # 14. the worked example at the published cutoff: the curve's chain inverts, and the
    #     dichotomy's likelihood ratio is the one its own 2x2 counts imply
    lr_at_cutoff = float(fitted.lr(1658.0))
    for pre in (0.10, 0.20):
        posterior = model.post_test(pre, lr_at_cutoff)["post_test_probability"]
        recovered = (posterior / (1 - posterior)) / (pre / (1 - pre))
        check(f"the 1,658 pg/mL chain inverts to the likelihood ratio at a {pre:.0%} "
              "pre-test probability",
              np.isclose(recovered, lr_at_cutoff, rtol=1e-12),
              f"LR {lr_at_cutoff:.3f} -> {posterior:.1%}")
    positive = values >= 1658.0
    true_positive, false_positive = int(positive[is_case].sum()), int(positive[~is_case].sum())
    positive_lr = (true_positive / is_case.sum()) / (false_positive / (~is_case).sum())
    check("the cutoff's likelihood ratio is sensitivity / (1 - specificity) from its counts",
          np.isclose(positive_lr,
                     (21 / 26) / (4 / 140), rtol=1e-12)
          and (true_positive, false_positive) == (21, 4),
          f"({true_positive}/26) / ({false_positive}/140) = {positive_lr:.2f}")
    check("the dichotomy claims more at the cutoff than the curve does there",
          positive_lr > 10 * lr_at_cutoff,
          f"LR {positive_lr:.1f} against {lr_at_cutoff:.2f} at the same concentration")

    # 15. the complete-case table is computed on one set of patients, the same for
    #     every marker in it
    clinical = data.load_clinical()
    complete_case = run_analysis.complete_case_marker_table(clinical)
    columns = list(run_analysis.COMPLETE_CASE_COLUMNS)
    expected = clinical.dropna(subset=columns)
    check("the complete-case subset is the patients with every clinical measurement present",
          complete_case["n"] == len(expected) == int(clinical[columns].notna().all(axis=1).sum()),
          f"n = {complete_case['n']} of {len(clinical)}")
    check("every marker in the complete-case table uses that same subset",
          all(row["n"] == complete_case["n"] for row in complete_case["markers"])
          and all(row["n_pvod"] == complete_case["n_pvod"]
                  for row in complete_case["markers"]),
          f"{len(complete_case['markers'])} markers, all at n = {complete_case['n']} "
          f"({complete_case['n_pvod']} PVOD)")
    check("the per-marker table does not use one subset, which is why the check is needed",
          len({row["n"] for row in run_analysis.marker_table(clinical)}) > 1,
          f"denominators {sorted({row['n'] for row in run_analysis.marker_table(clinical)})}")

    # 17. Patient_ID: present, valid and unique wherever a row is a patient
    patients = data.load_patients_166()
    check("Patients_166 gives one row per patient, Patient_ID 1-166 unique",
          patients.Patient_ID.tolist() == list(range(1, 167)),
          f"{len(patients)} patients")
    check("the clinical frame carries a valid unique Patient_ID on every row",
          "patient_id" in clinical.columns
          and clinical.patient_id.notna().all()
          and clinical.patient_id.is_unique
          and clinical.patient_id.between(1, 166).all()
          and len(clinical) == 166,
          f"{clinical.patient_id.nunique()} distinct ids over {len(clinical)} rows")

    # the join is by key, never by position: the clinical frame is in distribution-sheet
    # order and Patient_ID is in 6MWD-1 order, so a positional join would mislabel rows
    by_id = patients.set_index("Patient_ID")
    check("Patient_ID is joined on age-adjusted GDF15 exactly, group agreeing on every row",
          all(by_id.AgeAdj_GDF15[pid] == value and by_id.Group[pid] == group
              for pid, value, group
              in zip(clinical.patient_id, clinical.ageadj_gdf15, clinical.group)),
          "166 of 166 exact, zero residual")
    check("that join is not positional, so the key is doing real work",
          clinical.patient_id.tolist() != list(range(1, 167)),
          f"frame order starts {clinical.patient_id.head(3).tolist()}")

    ntprobnp_rows = data.load_ntprobnp_correlation()
    check("the NT-proBNP per-patient table carries Patient_ID",
          "patient_id" in ntprobnp_rows.columns
          and ntprobnp_rows.patient_id.is_unique
          and ntprobnp_rows.patient_id.between(1, 26).all(),
          f"{len(ntprobnp_rows)} PVOD patients")

    # 18. the complete-case subset is exactly the approved one, and its ids are the
    #     recomputed intersection rather than a stored list
    check("the complete-case subset is n = 23, PVOD 10, IPAH 7, HPAH 6",
          complete_case["n"] == 23
          and complete_case["n_by_group"] == {"PVOD": 10, "IPAH": 7, "HPAH": 6},
          f"n = {complete_case['n']} {complete_case['n_by_group']}")
    check("the reported Patient_IDs are the actual intersection of the five measurements",
          complete_case["patient_ids"] == sorted(int(i) for i in expected.patient_id),
          f"{len(complete_case['patient_ids'])} ids: "
          f"{complete_case['patient_ids'][:4]}...{complete_case['patient_ids'][-2:]}")
    check("those ids are unique and inside the cohort",
          len(set(complete_case["patient_ids"])) == 23
          and all(1 <= i <= 166 for i in complete_case["patient_ids"]),
          "23 distinct ids in 1-166")
    check("the subset is case-enriched against the full cohort, computed not asserted",
          complete_case["pvod_fraction"] > 2 * complete_case["full_cohort_pvod_fraction"],
          f"PVOD {complete_case['pvod_fraction']:.1%} against "
          f"{complete_case['full_cohort_pvod_fraction']:.1%}, "
          f"{complete_case['case_enrichment']:.1f}x")

    # 19. raw GDF15 is not in the complete-case table and no AUC of 1.00 is printed there
    check("the complete-case table carries the clinical comparators only, no GDF15 row",
          not any("GDF15" in row["marker"] or "gdf15" in row["key"]
                  for row in complete_case["markers"]),
          f"{[row['marker'] for row in complete_case['markers']]}")
    check("no AUC of 1.00 is reported in the complete-case table",
          all(row["auc"] < 1.0 for row in complete_case["markers"]),
          f"largest {max(row['auc'] for row in complete_case['markers']):.4f}")

    # the separation is real, and is what the removed row would have been reporting
    separation = complete_case["ageadj_gdf15_separation"]
    case = expected.ageadj_gdf15[expected.group == "PVOD"]
    other = expected.ageadj_gdf15[expected.group != "PVOD"]
    check("age-adjusted GDF15 separates the subset completely, ranges recomputed",
          separation["complete_separation"]
          and separation["pvod_min"] == float(case.min())
          and separation["pvod_max"] == float(case.max())
          and separation["comparator_min"] == float(other.min())
          and separation["comparator_max"] == float(other.max())
          and float(case.min()) > float(other.max()),
          f"PVOD {separation['pvod_min']:.2f}-{separation['pvod_max']:.2f} against "
          f"{separation['comparator_min']:.2f}-{separation['comparator_max']:.2f}")
    check("that separation is exactly what an AUC of 1.00 would have been reporting",
          model.auc(expected.ageadj_gdf15, (expected.group == "PVOD").to_numpy()) == 1.0,
          "complete separation, reported as two ranges instead")

    # 20. raw GDF15 is never joined to a patient by an invented key
    check("no Patient_ID is attached to the raw-GDF15 leave-one-out cases",
          all("patient_id" not in refit for refit in loo["refits"]),
          f"{len(loo['refits'])} refits, identified by source-sheet position only")
    gdf15_frame = data.load_gdf15()
    check("the raw-GDF15 frame carries no patient identifier to join on",
          "patient_id" not in gdf15_frame.columns
          and set(gdf15_frame.columns) == {"cohort", "gdf15"},
          f"columns {sorted(gdf15_frame.columns)}")
    check("raw source covers six cohorts; clinical reconstruction covers group 1 only",
          len(gdf15_frame) == 322 and len(clinical) == 166,
          "322 concentrations over six cohorts against 166 group-1 patients")

    # 21. completeness by diagnosis recomputes from the reconstruction
    completeness = run_analysis.completeness_by_diagnosis(patients)
    totals = patients.Group.value_counts().to_dict()
    recomputed = all(
        row[f"{group}_present"] == int(patients[patients.Group == group][row["column"]].notna().sum())
        and row[f"{group}_total"] == totals[group]
        and abs(row[f"{group}_percent"]
                - 100.0 * row[f"{group}_present"] / totals[group]) < 1e-9
        for row in completeness for group in data.GROUP1)
    check("completeness percentages recompute from the reconstructed workbook",
          recomputed and len(completeness) == 5,
          f"{len(completeness)} variables over {len(totals)} diagnoses")
    check("measurement availability is associated with diagnosis, most for cardiac index",
          next(r for r in completeness if r["variable"] == "Cardiac index")["PVOD_percent"]
          > 2 * next(r for r in completeness if r["variable"] == "Cardiac index")["IPAH_percent"],
          "cardiac index PVOD "
          f"{next(r for r in completeness if r['variable'] == 'Cardiac index')['PVOD_percent']:.0f}%"
          " against IPAH "
          f"{next(r for r in completeness if r['variable'] == 'Cardiac index')['IPAH_percent']:.0f}%")
    # the two NT-proBNP source definitions, each right for its own question
    check("NT-proBNP availability counts both sources; the subset uses the distribution sheet",
          next(r for r in completeness if r["variable"] == "NT-proBNP")["PVOD_present"] == 17
          and int(clinical[clinical.group == "PVOD"].ntprobnp.notna().sum()) == 15,
          "17 of 26 PVOD measured in either sheet, 15 of 26 in the distribution sheet")

    # 16. the whole workflow is reproducible: run twice into two empty directories, every
    #     output file must be byte for byte identical, figures included
    digests = []
    for _ in range(2):
        with tempfile.TemporaryDirectory() as scratch:
            root = Path(scratch)
            (root / "results").mkdir()
            (root / "figures").mkdir()
            run_results, figure_results, figure_figures = (
                run_analysis.RESULTS, make_figures.RESULTS, make_figures.FIGURES)
            try:
                run_analysis.RESULTS = root / "results"
                make_figures.RESULTS = root / "results"
                make_figures.FIGURES = root / "figures"
                with contextlib.redirect_stdout(io.StringIO()):
                    run_analysis.main()
                digests.append({
                    path.relative_to(root).as_posix():
                        hashlib.sha256(path.read_bytes()).hexdigest()
                    for path in sorted(root.rglob("*")) if path.is_file()})
            finally:
                run_analysis.RESULTS = run_results
                make_figures.RESULTS = figure_results
                make_figures.FIGURES = figure_figures
    check("two clean runs produce every output file byte for byte identical",
          digests[0] == digests[1] and len(digests[0]) > 0,
          f"{len(digests[0])} files")

    independent = verify_centerpiece.verify(data.WORKBOOK, run_analysis.RESULTS)
    for name, passed in independent['checks'].items():
        check('centerpiece: ' + name, passed)
    for name, passed in audit_numbers.audit(data.WORKBOOK, run_analysis.RESULTS)['checks'].items():
        check('manuscript: ' + name, passed)
    patients, _, provenance = reconstruct.reconstruct(data.WORKBOOK)
    check('source-only reconstruction has 166 unique keys',
          len(patients) == 166 and patients.AgeAdj_GDF15.is_unique)
    check('clinical key residuals below recorded source precision',
          provenance.key_residual.abs().max() < 1e-6)
    try:
        reconstruct.unique_match([1., 1.], 1.)
        rejected = False
    except ValueError:
        rejected = True
    check('ambiguous clinical keys are rejected', rejected)
    try:
        reconstruct.unique_match([1., 2.], 3.)
        rejected = False
    except ValueError:
        rejected = True
    check('missing clinical keys are rejected', rejected)
    verify_figure1.verify(run_analysis.RESULTS, check)
    print(f"\nall {len(CHECKS)} checks passed")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workbook', type=Path, default=data.WORKBOOK)
    parser.add_argument('--results-dir', type=Path, default=run_analysis.RESULTS)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    data.WORKBOOK = args.workbook.resolve()
    run_analysis.RESULTS = make_figures.RESULTS = args.results_dir
    try:
        main()
    finally:
        if args.report:
            args.report.parent.mkdir(parents=True, exist_ok=True)
            args.report.write_text(json.dumps({'checks': CHECKS, 'count':len(CHECKS),
                'failures':sum(not c['passed'] for c in CHECKS)}, indent=2) + '\n')
