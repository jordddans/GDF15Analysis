"""Continuous log10(GDF15) logistic model, converted to a likelihood-ratio curve.

The model is a single-predictor logistic regression

    logit P(PVOD | x) = b0 + b1 * log10(GDF15)

fitted on a cohort in which PVOD cases are over-represented by design.  The fitted
intercept therefore carries the sampling prevalence of that cohort, and the predicted
probabilities apply only to a population with that same prevalence.

Subtracting the sampling log-odds removes it:

    log LR(x) = [b0 + b1 * log10(x)] - log(pi / (1 - pi))

where pi is the proportion of PVOD cases in the fitted sample.  What remains is the
likelihood ratio, a property of the two conditional distributions of GDF15 given
disease status, which does not depend on how the sample was assembled.  A clinician
supplies the pre-test probability; the model supplies only the update:

    pre-test odds x LR = post-test odds

Uncertainty for the fitted curve comes from a parametric bootstrap: outcomes are
simulated from the fitted model at the observed GDF15 values and the model is refitted,
so the interval reflects uncertainty in the fitted parameters and is defined across the
whole GDF15 range rather than only where cases were observed.  Count-based quantities
use exact (Clopper-Pearson) binomial limits.

Three further procedures sit alongside that band and answer different questions.
``stability_bootstrap`` resamples patients, stratified by class, and recomputes both a
selected threshold and a landmark on the fitted curve inside every replicate, so the two
are compared under one scheme.  ``bootstrap_optimism`` is the internal validation: it
corrects the apparent C statistic and calibration slope for the optimism of fitting and
evaluating on the same patients.  ``leave_one_case_out`` refits the model once per case
with that case removed and reads the likelihood ratio at fixed concentrations, which
measures how much of the curve rests on any single patient; it is a sensitivity analysis
and corrects nothing.

``functional_form_test`` justifies the one modelling choice the curve makes: that
log10(GDF15) enters linearly.
"""

from dataclasses import dataclass

import numpy as np
import statsmodels.api as sm
from scipy.stats import beta, chi2, norm

SEED = 1
N_BOOT = 4000


def _logit(p):
    return np.log(p / (1 - p))


@dataclass
class LRModel:
    """A fitted log10(GDF15) model and the likelihood-ratio curve derived from it."""

    name: str
    n: int
    n_events: int
    prevalence: float
    intercept: float
    slope: float
    se_intercept: float
    se_slope: float
    cov: np.ndarray
    gdf15: np.ndarray
    is_case: np.ndarray

    @property
    def sampling_log_odds(self):
        """The log-odds of being a case in the fitted sample, removed to obtain the LR."""
        return _logit(self.prevalence)

    def log_lr(self, gdf15):
        """Natural log of the likelihood ratio at a GDF15 concentration (pg/mL)."""
        x = np.log10(np.asarray(gdf15, dtype=float))
        return self.intercept + self.slope * x - self.sampling_log_odds

    def lr(self, gdf15):
        """Likelihood ratio for PVOD carried by a GDF15 concentration (pg/mL)."""
        return np.exp(self.log_lr(gdf15))

    def probability(self, gdf15):
        """Fitted probability of PVOD -- valid only at the sampling prevalence."""
        x = np.log10(np.asarray(gdf15, dtype=float))
        return 1 / (1 + np.exp(-(self.intercept + self.slope * x)))

    def gdf15_at_lr(self, target):
        """The GDF15 concentration whose likelihood ratio equals ``target``."""
        return 10 ** ((np.log(target) + self.sampling_log_odds - self.intercept) / self.slope)

    def coef_ci(self, level=0.95):
        z = norm.ppf(0.5 + level / 2)
        return {
            "intercept": (self.intercept - z * self.se_intercept,
                          self.intercept + z * self.se_intercept),
            "slope": (self.slope - z * self.se_slope, self.slope + z * self.se_slope),
        }


def fit_lr_model(gdf15, is_case, name):
    """Fit the single-predictor model and return an :class:`LRModel`."""
    gdf15 = np.asarray(gdf15, dtype=float)
    is_case = np.asarray(is_case).astype(int)
    design = sm.add_constant(np.log10(gdf15))
    fit = sm.Logit(is_case, design).fit(disp=0)
    cov = np.asarray(fit.cov_params())
    return LRModel(
        name=name,
        n=len(is_case),
        n_events=int(is_case.sum()),
        prevalence=float(is_case.mean()),
        intercept=float(fit.params[0]),
        slope=float(fit.params[1]),
        se_intercept=float(fit.bse[0]),
        se_slope=float(fit.bse[1]),
        cov=cov,
        gdf15=gdf15,
        is_case=is_case,
    )



def functional_form_test(gdf15, is_case):
    """Is a quadratic term in log10(GDF15) needed?  Likelihood-ratio test and AIC.

    The reported model enters log10(GDF15) linearly.  This refits it with a squared term
    added and compares the two by the usual likelihood-ratio test on one degree of
    freedom, alongside each model's AIC.  It justifies the functional form actually used;
    it is not a model-selection step, and the linear model is reported either way.
    """
    x = np.log10(np.asarray(gdf15, dtype=float))
    y = np.asarray(is_case).astype(int)
    linear = sm.Logit(y, sm.add_constant(x)).fit(disp=0)
    quadratic = sm.Logit(y, sm.add_constant(np.column_stack([x, x ** 2]))).fit(disp=0)
    statistic = 2 * (float(quadratic.llf) - float(linear.llf))
    aic = lambda fit, k: -2 * float(fit.llf) + 2 * k
    return {
        "linear_loglikelihood": float(linear.llf),
        "quadratic_loglikelihood": float(quadratic.llf),
        "lr_statistic": statistic,
        "df": 1,
        "p_value": float(chi2.sf(statistic, 1)),
        "linear_aic": aic(linear, 2),
        "quadratic_aic": aic(quadratic, 3),
        "aic_difference": aic(quadratic, 3) - aic(linear, 2),
        "quadratic_coefficient": float(quadratic.params[2]),
        "quadratic_coefficient_se": float(quadratic.bse[2]),
        "preferred": "linear" if aic(linear, 2) <= aic(quadratic, 3) else "quadratic",
        "method": ("likelihood-ratio test of log10(GDF15) against "
                   "log10(GDF15) + log10(GDF15)^2, one degree of freedom, with the AIC "
                   "of each model"),
    }

def post_test(pre_test_probability, likelihood_ratio):
    """pre-test probability -> pre-test odds -> x LR -> post-test odds -> probability."""
    pre_odds = pre_test_probability / (1 - pre_test_probability)
    post_odds = pre_odds * likelihood_ratio
    return {
        "pre_test_probability": pre_test_probability,
        "pre_test_odds": pre_odds,
        "likelihood_ratio": likelihood_ratio,
        "post_test_odds": post_odds,
        "post_test_probability": post_odds / (1 + post_odds),
    }


def parametric_bootstrap(model, n_boot=N_BOOT, seed=SEED):
    """Refit the model on outcomes simulated from it; return the replicate parameters.

    GDF15 values are held at their observed values and only the outcomes are simulated,
    so no replicate invents a concentration that was never measured.  Replicates that
    fail to converge, or that separate (a non-positive or unbounded slope), are dropped
    and counted.
    """
    import warnings

    rng = np.random.default_rng(seed)
    design = sm.add_constant(np.log10(model.gdf15))
    fitted = model.probability(model.gdf15)
    intercepts, slopes, offsets = [], [], []
    discarded = 0
    # a replicate that fails to converge is detected below and counted, so the optimiser's
    # own warning about it carries no information the discard count does not already give
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for _ in range(n_boot):
            outcome = rng.binomial(1, fitted)
            if outcome.sum() < 2 or outcome.sum() > len(outcome) - 2:
                discarded += 1
                continue
            try:
                fit = sm.Logit(outcome, design).fit(disp=0)
            except Exception:
                discarded += 1
                continue
            if not fit.mle_retvals.get("converged", False) or not 0 < fit.params[1] < 100:
                discarded += 1
                continue
            intercepts.append(float(fit.params[0]))
            slopes.append(float(fit.params[1]))
            offsets.append(_logit(outcome.mean()))
    return {
        "intercept": np.array(intercepts),
        "slope": np.array(slopes),
        "sampling_log_odds": np.array(offsets),
        "n_kept": len(intercepts),
        "n_discarded": discarded,
    }


def lr_curve_ci(model, grid, draws, level=0.95):
    """Pointwise percentile interval for the LR curve over ``grid`` (pg/mL)."""
    x = np.log10(np.asarray(grid, dtype=float))
    log_lr = (draws["intercept"][:, None]
              + draws["slope"][:, None] * x[None, :]
              - draws["sampling_log_odds"][:, None])
    lo, hi = (1 - level) / 2 * 100, (1 + level) / 2 * 100
    return np.exp(np.percentile(log_lr, lo, axis=0)), np.exp(np.percentile(log_lr, hi, axis=0))


def gdf15_at_lr_ci(model, target, draws, level=0.95):
    """Percentile interval for the GDF15 value at which the LR equals ``target``."""
    values = 10 ** ((np.log(target) + draws["sampling_log_odds"] - draws["intercept"])
                    / draws["slope"])
    lo, hi = (1 - level) / 2 * 100, (1 + level) / 2 * 100
    return float(np.percentile(values, lo)), float(np.percentile(values, hi))


def clopper_pearson(successes, total, level=0.95):
    """Exact binomial confidence interval for a proportion."""
    alpha = 1 - level
    lower = 0.0 if successes == 0 else beta.ppf(alpha / 2, successes, total - successes + 1)
    upper = 1.0 if successes == total else beta.ppf(1 - alpha / 2, successes + 1, total - successes)
    return float(lower), float(upper)


def calibration_by_quintile(model):
    """Observed vs predicted PVOD frequency in fifths of predicted probability."""
    predicted = model.probability(model.gdf15)
    order = np.argsort(predicted)
    rows = []
    for quintile, chunk in enumerate(np.array_split(order, 5), start=1):
        observed = model.is_case[chunk]
        successes, total = int(observed.sum()), len(chunk)
        low, high = clopper_pearson(successes, total)
        rows.append({
            "quintile": quintile,
            "n": total,
            "gdf15_min": float(model.gdf15[chunk].min()),
            "gdf15_max": float(model.gdf15[chunk].max()),
            "predicted_mean": float(predicted[chunk].mean()),
            "predicted_min": float(predicted[chunk].min()),
            "predicted_max": float(predicted[chunk].max()),
            "observed_events": successes,
            "observed_proportion": successes / total,
            "observed_ci_low": low,
            "observed_ci_high": high,
        })
    return rows


def auc(scores, is_case):
    """Area under the ROC curve, computed from the Mann-Whitney statistic with midranks."""
    scores = np.asarray(scores, dtype=float)
    is_case = np.asarray(is_case).astype(bool)
    case, control = scores[is_case], scores[~is_case]
    ranks = _midrank(np.concatenate([case, control]))
    return float((ranks[:len(case)].sum() - len(case) * (len(case) + 1) / 2)
                 / (len(case) * len(control)))


def _midrank(values):
    order = np.argsort(values)
    sorted_values = values[order]
    ranks = np.empty(len(values), dtype=float)
    i = 0
    while i < len(values):
        j = i
        while j < len(values) - 1 and sorted_values[j + 1] == sorted_values[i]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1
        i = j + 1
    return ranks


def auc_ci_asymptotic(scores, is_case, level=0.95):
    """Asymptotic confidence interval for a single C statistic (AUC).

    The variance is the structural-component estimate: the AUC is the mean over cases of
    each case's probability of outranking a random control, and the variance of that mean
    is the sum of the between-case and between-control component variances.  The components
    are computed here from midranks; ``verify.py`` checks them against the same quantities
    computed directly from the definition.
    """
    scores = np.asarray(scores, dtype=float)
    is_case = np.asarray(is_case).astype(bool)
    case, control = scores[is_case], scores[~is_case]
    m, n = len(case), len(control)
    area = auc(scores, is_case)
    v10 = (_midrank(np.concatenate([case, control]))[:m] - _midrank(case)) / n
    v01 = 1 - (_midrank(np.concatenate([control, case]))[:n] - _midrank(control)) / m
    variance = v10.var(ddof=1) / m + v01.var(ddof=1) / n
    z = norm.ppf(0.5 + level / 2)
    half = z * np.sqrt(variance)
    return area, (max(0.0, area - half), min(1.0, area + half))


def youden_cutoff(gdf15, is_case):
    """Cutoff maximising sensitivity + specificity - 1; ties resolved to the highest."""
    gdf15 = np.asarray(gdf15, dtype=float)
    is_case = np.asarray(is_case).astype(bool)
    best_j, best_cut = -np.inf, np.nan
    for candidate in np.unique(gdf15):
        positive = gdf15 >= candidate
        sensitivity = positive[is_case].mean()
        specificity = 1 - positive[~is_case].mean()
        j = sensitivity + specificity - 1
        if j >= best_j - 1e-12:
            best_j, best_cut = max(j, best_j), candidate
    return float(best_cut)


def stratified_indices(is_case, rng):
    """One stratified case-resample: cases and comparators drawn separately, at their sizes."""
    is_case = np.asarray(is_case).astype(bool)
    cases, controls = np.where(is_case)[0], np.where(~is_case)[0]
    return np.concatenate([rng.choice(cases, len(cases), replace=True),
                           rng.choice(controls, len(controls), replace=True)])


def stability_bootstrap(fitted, n_boot=N_BOOT, seed=SEED, lr_target=10, level=0.95):
    """Youden cutoff and the LR = ``lr_target`` concentration, from one shared resample.

    PVOD and comparator patients are resampled separately at their observed sample sizes,
    and *both* quantities are recomputed inside every replicate from the same patients, so
    a selected threshold and a landmark on the fitted curve are compared on identical
    footing.  A replicate whose two classes no longer overlap, or whose fit does not
    converge, has no finite likelihood-ratio curve; it is dropped from both quantities
    rather than from one, so the paired comparison is never made on unequal replicate sets.
    """
    import warnings

    gdf15 = np.asarray(fitted.gdf15, dtype=float)
    is_case = np.asarray(fitted.is_case).astype(bool)
    x = np.log10(gdf15)
    rng = np.random.default_rng(seed)
    cutoffs, landmarks = [], []
    dropped = 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for _ in range(n_boot):
            take = stratified_indices(is_case, rng)
            outcome, values = is_case[take], gdf15[take]
            if values[outcome].min() > values[~outcome].max():
                dropped += 1
                continue
            try:
                fit = sm.Logit(outcome.astype(int), sm.add_constant(x[take])).fit(disp=0)
            except Exception:
                dropped += 1
                continue
            if not fit.mle_retvals.get("converged", False) or not 0 < fit.params[1] < 100:
                dropped += 1
                continue
            offset = _logit(outcome.mean())
            landmarks.append(10 ** ((np.log(lr_target) + offset - fit.params[0]) / fit.params[1]))
            cutoffs.append(youden_cutoff(values, outcome))

    lo, hi = (1 - level) / 2 * 100, (1 + level) / 2 * 100
    cutoffs, landmarks = np.array(cutoffs), np.array(landmarks)
    youden = {"point": youden_cutoff(gdf15, is_case),
              "ci_low": float(np.percentile(cutoffs, lo)),
              "ci_high": float(np.percentile(cutoffs, hi))}
    landmark = {"point": float(fitted.gdf15_at_lr(lr_target)),
                "ci_low": float(np.percentile(landmarks, lo)),
                "ci_high": float(np.percentile(landmarks, hi))}
    for entry in (youden, landmark):
        entry["ci_width"] = entry["ci_high"] - entry["ci_low"]
    return {
        "youden": youden,
        "landmark": landmark,
        "lr_target": lr_target,
        # undefined when a run is too short for the landmark interval to have any width
        "width_ratio": (youden["ci_width"] / landmark["ci_width"]
                        if landmark["ci_width"] > 0 else float("nan")),
        "n_replicates": len(landmarks),
        "n_dropped": dropped,
        "method": (f"stratified case bootstrap, {n_boot:,} replicates, seed {seed}; PVOD and "
                   "comparator patients resampled separately at their observed sample sizes, "
                   "both quantities recomputed in every replicate"),
    }


def bootstrap_optimism(fitted, n_boot=N_BOOT, seed=SEED):
    """Bootstrap optimism correction (Harrell) for the fitted single-predictor model.

    In each replicate the model is refitted on a bootstrap sample of the patients and
    scored twice: on the sample it was fitted to, and on the original data.  The
    difference is that replicate's optimism, and the corrected value is the apparent
    performance minus the mean optimism.

    Two quantities are corrected.  The C statistic measures discrimination.  The
    calibration slope is the coefficient obtained by regressing the outcome on the linear
    predictor; in the sample a model was fitted to it is exactly 1 by construction, so the
    apparent slope is 1 and the correction is what the resampling has to say about it.

    Optimism is not constrained in sign: for a strongly separating single predictor it can
    come out slightly negative, so the corrected value may sit marginally above the
    apparent one.  Replicates that fail to converge, or that separate, are dropped and
    counted.
    """
    import warnings

    gdf15 = np.asarray(fitted.gdf15, dtype=float)
    outcome = np.asarray(fitted.is_case).astype(int)
    x = np.log10(gdf15)
    design = sm.add_constant(x)
    n = len(outcome)
    rng = np.random.default_rng(seed)
    auc_optimism, slope_optimism = [], []
    dropped = 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for _ in range(n_boot):
            take = rng.integers(0, n, n)
            resampled = outcome[take]
            if resampled.sum() < 2 or resampled.sum() > n - 2:
                dropped += 1
                continue
            try:
                fit = sm.Logit(resampled, design[take]).fit(disp=0)
                if (not fit.mle_retvals.get("converged", False)
                        or not 0 < fit.params[1] < 100):
                    dropped += 1
                    continue
                # the replicate model's linear predictor, in each sample in turn
                lp_replicate = fit.params[0] + fit.params[1] * x[take]
                lp_original = fit.params[0] + fit.params[1] * x
                recalibrated = sm.Logit(outcome, sm.add_constant(lp_original)).fit(disp=0)
                if not recalibrated.mle_retvals.get("converged", False):
                    dropped += 1
                    continue
            except Exception:
                dropped += 1
                continue
            auc_optimism.append(auc(lp_replicate, resampled.astype(bool))
                                - auc(lp_original, outcome.astype(bool)))
            slope_optimism.append(1.0 - float(recalibrated.params[1]))

    apparent_auc = auc(gdf15, outcome.astype(bool))
    mean_auc_optimism = float(np.mean(auc_optimism))
    mean_slope_optimism = float(np.mean(slope_optimism))
    return {
        "apparent_auc": float(apparent_auc),
        "auc_optimism": mean_auc_optimism,
        "corrected_auc": float(apparent_auc - mean_auc_optimism),
        "apparent_calibration_slope": 1.0,
        "calibration_slope_optimism": mean_slope_optimism,
        "corrected_calibration_slope": 1.0 - mean_slope_optimism,
        "n_replicates": len(auc_optimism),
        "n_dropped": dropped,
        "method": (f"bootstrap optimism correction, {n_boot:,} replicates, seed {seed}; "
                   "patients resampled with replacement, the model refitted in each "
                   "replicate and scored on both the replicate and the original sample"),
    }


def leave_one_case_out(fitted, evaluation_points):
    """Refit the model once per case, dropping that case, and read the LR at fixed points.

    One PVOD case is removed, the model is refitted on the remaining patients, and the
    likelihood ratio is read at each of ``evaluation_points`` -- the same concentrations
    in every refit.  Holding the evaluation points fixed is what makes the answer
    interpretable: the influence of a single case on the likelihood ratio depends on where
    on the curve it is measured, so a change quoted without its concentration cannot be
    checked.  Each case is identified by its position in the PVOD column of the source
    sheet and by its own GDF15 value, which is what makes it distinguishable at all.

    These cases carry no ``Patient_ID`` and cannot be given one.  The model reads sheet
    ``GDF15-1``, which is anonymous group-level concentrations -- no identifier, no age,
    no row order shared with the per-patient sheets -- so there is nothing to join on.
    A ``case_id`` here is a position in a source sheet, not a patient, and it has no
    relation to the ``Patient_ID`` used in the patient-level tables.

    This is a sensitivity analysis and not a validation procedure: nothing here corrects
    the reported model, which remains the fit on all cases.
    """
    gdf15 = np.asarray(fitted.gdf15, dtype=float)
    is_case = np.asarray(fitted.is_case).astype(int)
    points = np.asarray(evaluation_points, dtype=float)
    full_cohort_lr = np.asarray(fitted.lr(points), dtype=float)

    refits = []
    for case_id, position in enumerate(np.where(is_case == 1)[0], start=1):
        keep = np.ones(len(is_case), dtype=bool)
        keep[position] = False
        refit = fit_lr_model(gdf15[keep], is_case[keep],
                             f"{fitted.name} without PVOD case {case_id}")
        lr = np.asarray(refit.lr(points), dtype=float)
        refits.append({
            "case_id": case_id,
            "excluded_gdf15": float(gdf15[position]),
            "n": refit.n,
            "n_events": refit.n_events,
            "intercept": refit.intercept,
            "slope": refit.slope,
            "lr": lr,
            "percent_change": 100 * (lr - full_cohort_lr) / full_cohort_lr,
        })

    at_point = []
    for index, point in enumerate(points):
        changes = np.array([abs(r["percent_change"][index]) for r in refits])
        worst = refits[int(np.argmax(changes))]
        at_point.append({
            "gdf15": float(point),
            "full_cohort_lr": float(full_cohort_lr[index]),
            "max_absolute_percent_change": float(changes.max()),
            "signed_percent_change": float(worst["percent_change"][index]),
            "lr_without_that_case": float(worst["lr"][index]),
            "case_id": worst["case_id"],
            "excluded_gdf15": worst["excluded_gdf15"],
            "median_absolute_percent_change": float(np.median(changes)),
        })

    return {
        "n_refits": len(refits),
        "evaluation_points": [float(p) for p in points],
        "at_evaluation_point": at_point,
        "refits": refits,
        "method": ("the model refitted once per PVOD case with that case removed, the "
                   "likelihood ratio read at fixed concentrations in every refit; cases "
                   "identified by position in the PVOD column of sheet GDF15-1 and by "
                   "their own GDF15 value. That sheet is anonymous, so these cases have "
                   "no Patient_ID and a case_id is not a patient"),
    }
