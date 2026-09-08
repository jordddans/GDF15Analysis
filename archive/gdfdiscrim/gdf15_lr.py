"""GDF15 -> PVOD likelihood-ratio engine.

Not validated for clinical use; requires prospective validation in an independent cohort.

Built incrementally by the notebooks in this project:
  03_lr_curves.ipynb  - fit_lr_model, lr, lr_ci, banded_lr
  04_honesty.ipynb    - rule-of-three ceiling, extrapolation-zone flags
  05_bayes_engine.ipynb - post_test, Fagan-nomogram helpers

Likelihood ratios only. This module never fits or returns P(PVOD | GDF15) directly -
see spec section 0 for why the 26/322 sample fraction must never be treated as a prevalence.
"""
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy import stats


@dataclass
class LRModel:
    """A fitted two-parameter lognormal LR model for one comparator definition."""
    name: str
    mu1: float
    sigma1: float
    mu0: float
    sigma0: float
    n1: int
    n0: int
    pvod_vals: np.ndarray
    comp_vals: np.ndarray
    # filled in by 04_honesty.ipynb
    ceiling: Optional[float] = None
    support_pvod: Optional[tuple] = None
    support_comp: Optional[tuple] = None


def _fit_lognormal(x):
    y = np.log(np.asarray(x, dtype=float))
    return float(y.mean()), float(y.std(ddof=1))


def fit_lr_model(pvod_vals, comparator_vals, name="model") -> LRModel:
    """Fit independent lognormal distributions to PVOD and comparator raw GDF15 values."""
    pvod_vals = np.asarray(pvod_vals, dtype=float)
    comparator_vals = np.asarray(comparator_vals, dtype=float)
    mu1, sigma1 = _fit_lognormal(pvod_vals)
    mu0, sigma0 = _fit_lognormal(comparator_vals)
    return LRModel(
        name=name, mu1=mu1, sigma1=sigma1, mu0=mu0, sigma0=sigma0,
        n1=len(pvod_vals), n0=len(comparator_vals),
        pvod_vals=pvod_vals, comp_vals=comparator_vals,
    )


def lr(x, model: LRModel):
    """Point estimate of the parametric LR at raw GDF15 value(s) x. No cap, no CI."""
    lx = np.log(np.asarray(x, dtype=float))
    num = stats.norm.pdf(lx, model.mu1, model.sigma1)
    den = stats.norm.pdf(lx, model.mu0, model.sigma0)
    return num / den


def lr_ci(x, model: LRModel, n_boot=4000, seed=1, ci=0.95):
    """Stratified bootstrap CI for the parametric LR at x.

    Resamples within PVOD and within the comparator separately (preserving group
    sizes), refits mu/sigma on each replicate, and evaluates LR(x) each time.
    Percentile method. Returns (point, lo, hi); x may be a scalar or array.
    """
    rng = np.random.RandomState(seed)
    x_arr = np.atleast_1d(np.asarray(x, dtype=float))
    lx = np.log(x_arr)

    boot = np.empty((n_boot, len(x_arr)))
    for i in range(n_boot):
        bs1 = rng.choice(model.pvod_vals, size=model.n1, replace=True)
        bs0 = rng.choice(model.comp_vals, size=model.n0, replace=True)
        m1, s1 = _fit_lognormal(bs1)
        m0, s0 = _fit_lognormal(bs0)
        boot[i] = stats.norm.pdf(lx, m1, s1) / stats.norm.pdf(lx, m0, s0)

    point = lr(x_arr, model)
    lo = np.percentile(boot, 100 * (1 - ci) / 2, axis=0)
    hi = np.percentile(boot, 100 * (1 + ci) / 2, axis=0)

    if np.isscalar(x) or np.ndim(x) == 0:
        return float(point[0]), float(lo[0]), float(hi[0])
    return point, lo, hi


def banded_lr(pvod_vals, comparator_vals, edges, n_boot=4000, seed=1, ci=0.95):
    """Model-free banded LR with Haldane-Anscombe continuity correction and a
    stratified bootstrap CI per band. Returns arrays (a_counts, b_counts, lr, ci_lo, ci_hi).
    """
    pvod_vals = np.asarray(pvod_vals, dtype=float)
    comparator_vals = np.asarray(comparator_vals, dtype=float)
    n1, n0 = len(pvod_vals), len(comparator_vals)

    def counts(vals):
        return np.array([np.sum((vals >= edges[i]) & (vals < edges[i + 1])) for i in range(len(edges) - 1)])

    a = counts(pvod_vals)
    b = counts(comparator_vals)
    point = ((a + 0.5) / n1) / ((b + 0.5) / n0)

    rng = np.random.RandomState(seed)
    boot = np.empty((n_boot, len(edges) - 1))
    for i in range(n_boot):
        bs1 = rng.choice(pvod_vals, size=n1, replace=True)
        bs0 = rng.choice(comparator_vals, size=n0, replace=True)
        ab, bb = counts(bs1), counts(bs0)
        boot[i] = ((ab + 0.5) / n1) / ((bb + 0.5) / n0)

    lo = np.percentile(boot, 100 * (1 - ci) / 2, axis=0)
    hi = np.percentile(boot, 100 * (1 + ci) / 2, axis=0)
    return a, b, point, lo, hi


# ---------------------------------------------------------------------------
# §4.1 / §4.2 honesty layer: rule-of-three ceiling and extrapolation zones.
#
# RULE_OF_THREE[key]['thresholds'] holds the Jeffreys-bound "defensible LR >="
# figure computed at each threshold in spec Table 4.1 (sens lower bound / FPR
# upper bound). The *reported* defensible floor at a given x is the running
# maximum of all thresholds <= x: a bound established at threshold t remains a
# valid conservative claim for any x >= t (higher GDF15 is only more extreme),
# so the strongest still-valid claim is the max, not the bound at the nearest
# threshold. For Model A this plateaus at 32.9 for all x >= 2000 pg/mL, which
# is exactly the flat "LR >= 33" figure the spec's display message uses.
RULE_OF_THREE = {
    'A': dict(thresholds={1658: 9.5, 2000: 32.9, 2500: 28.3, 3000: 23.9}, cap=50.0),
    # Model B's cap is not specified by the spec (only Model A's is, at 50).
    # We deliberately reuse 50 here rather than raising it to match B's higher
    # defensible ceiling (~45) -- one shared, conservative cap keeps A and B
    # directly comparable and never overstates B relative to A. See notebook 04.
    'B': dict(thresholds={1658: 12.9, 2000: 45.1, 2500: 38.7, 3000: 32.7}, cap=50.0),
}

# Observed support (min, max raw GDF15), spec §4.2.
SUPPORT = {
    'A': dict(pvod=(860.5, 5986.4), comp=(91.0, 1809.5)),
    'B': dict(pvod=(860.5, 5986.4), comp=(21.0, 1809.5)),
}


def defensible_floor(x, model_key):
    """Running-max Jeffreys/rule-of-three LOWER bound on LR, for the RULE-IN (high GDF15) tail.
    None if x is below the smallest evaluated threshold (1658).

    "The data firmly supports LR >= this" -- protects against UNDER-claiming caution when the
    raw point estimate has run away upward with no comparator data to constrain it."""
    thresholds = RULE_OF_THREE[model_key]['thresholds']
    applicable = [bound for t, bound in thresholds.items() if t <= x]
    return max(applicable) if applicable else None


def rule_out_floor(model: LRModel):
    """Mirror image of defensible_floor(), for the RULE-OUT (low GDF15) tail.

    Below the smallest observed PVOD value there is zero PVOD data to anchor how small the
    true LR could really be -- the raw parametric point estimate there is pure extrapolation of
    the fitted curve's tail, unconstrained by any actual PVOD observation. This computes the
    Jeffreys UPPER bound on P(GDF15 <= pvod_min | PVOD) (0 PVOD observations that low, by
    definition of pvod_min) divided by the Jeffreys LOWER bound on P(GDF15 <= pvod_min |
    comparator) (well-estimated: comparator has plenty of data down there). That ratio is a
    defensible UPPER bound on the true region LR -- i.e. "the data cannot rule out that the true
    LR down there is as weak (close to 1, uninformative) as this value" -- so it is the LEAST
    extreme (safest, closest to 1) claim the data can defend anywhere below pvod_min. The
    reported LR must never be MORE extreme (further from 1, smaller) than this floor:
    reported = max(point, floor).

    Deliberately evaluated ONCE, anchored at the model's own PVOD minimum (well-populated on the
    comparator side -- typically ~60% of the comparator sample), rather than re-evaluated
    pointwise at each x. Evaluating this Jeffreys ratio near very small comparator counts (close
    to the comparator's own minimum, where the denominator count can be 0-2) is numerically
    unstable -- the Jeffreys lower bound on a near-zero count is itself near zero, so the ratio
    can explode. This single, stably-anchored value is then used as a flat floor for the entire
    region below pvod_min, exactly mirroring how RULE_OF_THREE's high-tail plateau (a small set
    of well-populated thresholds, not a continuous evaluation) avoids the same instability on
    the other tail. Unlike RULE_OF_THREE, which mirrors the spec's own published threshold
    numbers, this side of the honesty layer has no published spec table to match -- it is
    derived from first principles using the same Jeffreys machinery, applied symmetrically.
    """
    pvod_min = float(model.pvod_vals.min())
    n1, n0 = model.n1, model.n0
    x0 = int(np.sum(model.comp_vals <= pvod_min))
    jeff_upper_pvod = stats.beta.ppf(0.975, 0.5, n1 + 0.5)  # 0 PVOD observations <= their own min, by definition
    jeff_lower_comp = stats.beta.ppf(0.025, x0 + 0.5, n0 - x0 + 0.5)
    return jeff_upper_pvod / jeff_lower_comp


def rule_in_floor(model: LRModel):
    """Generic mirror of rule_out_floor(), for the RULE-IN (high-value) tail of ANY fitted
    LRModel -- not tied to the GDF15-specific RULE_OF_THREE threshold table (which mirrors the
    spec's own published, validated numbers and remains the primary tool for the GDF15 models).
    For markers with no curated clinical thresholds to anchor on, this gives a defensible LOWER
    bound on LR for the entire region beyond the comparator's observed maximum, anchored ONCE at
    that maximum (well-populated on the PVOD side) rather than evaluated pointwise -- same
    numerical-stability rationale as rule_out_floor. Symmetric construction: Jeffreys LOWER bound
    on P(marker >= comp_max | PVOD) divided by Jeffreys UPPER bound on P(marker >= comp_max |
    comparator) (0 comparator observations that high, by definition of comp_max).
    """
    comp_max = float(model.comp_vals.max())
    n1, n0 = model.n1, model.n0
    x1 = int(np.sum(model.pvod_vals >= comp_max))
    jeff_lower_pvod = stats.beta.ppf(0.025, x1 + 0.5, n1 - x1 + 0.5)
    jeff_upper_comp = stats.beta.ppf(0.975, 0.5, n0 + 0.5)  # 0 comparator observations >= their own max
    return jeff_lower_pvod / jeff_upper_comp


def extrapolation_zone(x, model_key):
    """'interpolation' | 'single_class_extrapolation' | 'full_extrapolation', per spec §4.2.

    Symmetric on both tails: full_extrapolation applies above the highest PVOD value ever
    observed (no data of either class up there) AND below the lowest comparator value ever
    observed (no data of either class down there either)."""
    pvod_lo, pvod_hi = SUPPORT[model_key]['pvod']
    comp_lo, comp_hi = SUPPORT[model_key]['comp']
    if x > pvod_hi or x < comp_lo:
        return 'full_extrapolation'
    if pvod_lo <= x <= comp_hi:
        return 'interpolation'
    return 'single_class_extrapolation'


def reported_lr(x, model: LRModel, model_key):
    """The honest, user-facing LR: capped/floored and zone-flagged per §4.1/§4.2, symmetric on
    both tails.

    `lr()` above is the raw, uncapped parametric point estimate used for model fitting and
    diagnostics. This function is what any user-facing output (including the Bayesian engine in
    05) should call.
    """
    x = float(x)
    zone = extrapolation_zone(x, model_key)
    point = float(lr(x, model))
    pvod_lo, pvod_hi = SUPPORT[model_key]['pvod']
    comp_lo, comp_hi = SUPPORT[model_key]['comp']

    if zone == 'full_extrapolation':
        if x > pvod_hi:
            floor = defensible_floor(pvod_hi, model_key)
            return dict(
                x=x, zone=zone, point=None, reported=floor, capped=True,
                message=(f"No GDF15 data of any kind (PVOD or comparator) exists above {pvod_hi:.0f} pg/mL. "
                          f"Refusing to return a point estimate; the data firmly supports LR >= {floor:.0f}, "
                          f"no more can be claimed."),
            )
        else:  # x < comp_lo
            ceiling = rule_out_floor(model)
            return dict(
                x=x, zone=zone, point=None, reported=ceiling, capped=True,
                message=(f"No GDF15 data of any kind (PVOD or comparator) exists below {comp_lo:.0f} pg/mL. "
                          f"Refusing to return a point estimate; the data cannot defend a rule-out claim "
                          f"stronger than LR <= {ceiling:.2f} here, no more extreme claim can be made."),
            )

    if zone == 'interpolation':
        return dict(x=x, zone=zone, point=point, reported=point, capped=False, message=None)

    # single_class_extrapolation -- two very different sub-cases, handled symmetrically
    if x > comp_hi:
        # rule-IN tail: only PVOD data exists; raw point estimate can run away UPWARD
        cap = RULE_OF_THREE[model_key]['cap']
        floor = defensible_floor(x, model_key)
        reported = min(point, cap)
        message = None
        if point > cap:
            message = (f"Point estimate ({point:.0f}) exceeds what {model.n0} comparator patients can resolve. "
                        f"The data firmly supports LR >= {floor:.0f}; the true value may be higher but this "
                        f"sample cannot demonstrate it. Reporting the capped value.")
        return dict(x=x, zone=zone, point=point, reported=reported, capped=point > cap,
                    floor=floor, message=message)
    else:
        # rule-OUT tail (x < pvod_lo): only comparator data exists; raw point estimate can run away DOWNWARD
        floor = rule_out_floor(model)
        reported = max(point, floor)
        message = None
        if point < floor:
            message = (f"Point estimate ({point:.3f}) is more extreme than {model.n1} PVOD patients (0 observed "
                        f"this low) can defend. The data cannot demonstrate rule-out evidence stronger than "
                        f"LR <= {floor:.2f} in this region; reporting the more conservative floor instead.")
        return dict(x=x, zone=zone, point=point, reported=reported, capped=point < floor,
                    floor=floor, message=message)


# ---------------------------------------------------------------------------
# §5 Bayesian chaining engine.
#
# No claim of clinical validity anywhere below: header disclaimer applies to
# every output this section produces.
DISCLAIMER = "Not validated for clinical use; requires prospective validation in an independent cohort."

INDEPENDENCE_DISCLAIMER = (
    "Chaining assumes conditional independence given disease status. Within group-1, GDF15 is "
    "statistically near-orthogonal to NT-proBNP, mPAP, PVR, and CI (spec Table 5) -- chaining GDF15 "
    "onto any of those is reasonably defensible. GDF15 correlates with 6MWD at rho=-0.42: chaining "
    "onto a functional-capacity measure will overstate the combined evidence. No data exists in this "
    "cohort on CT findings, DLCO, or EIF2AK4 genotype -- the most important pre-test elements in "
    "practice -- so independence from those is ASSUMED and UNTESTABLE with this dataset."
)

# Ships empty by design (spec §5): populate with fit_lr_model()-external, literature-sourced LRs.
# Every entry MUST carry a citation. Example shape (not populated):
#   EVIDENCE_LIBRARY['hrct_pvod_pattern'] = dict(
#       lr=4.0, citation="Author et al., Journal Year;Vol:pp — replace with a real source",
#       note="Illustrative shape only; this LR has not been sourced.",
#   )
EVIDENCE_LIBRARY = {}

_FUNCTIONAL_CAPACITY_MARKERS = ('6mwd', 'six_minute', 'walk', 'functional_capacity')


def chain_pretest(base_rate, evidence_names, library=None):
    """Multiply a base rate through a chain of named, cited literature LRs into a pre-test probability.

    Every entry used must exist in `library` (default: the module-level EVIDENCE_LIBRARY, which
    ships empty) and carry a non-empty 'citation' field, or this raises rather than silently
    proceeding with an unsourced number.
    """
    if library is None:
        library = EVIDENCE_LIBRARY
    if not (0 < base_rate < 1):
        raise ValueError("base_rate must be strictly between 0 and 1")

    odds = base_rate / (1 - base_rate)
    warnings, used = [], []
    for name in evidence_names:
        if name not in library:
            raise KeyError(f"'{name}' is not in the evidence library. Add it with lr= and citation= "
                            f"before chaining on it -- this tool never chains on an unsourced LR.")
        entry = library[name]
        if not entry.get('citation'):
            raise ValueError(f"Evidence entry '{name}' has no citation. Every chain entry must be sourced.")
        odds *= entry['lr']
        used.append(name)
        if any(term in name.lower() for term in _FUNCTIONAL_CAPACITY_MARKERS):
            warnings.append(
                f"'{name}' looks like a functional-capacity measure. GDF15 correlates with 6MWD at "
                f"rho=-0.42 in this cohort (spec §5) -- chaining GDF15 onto this will overstate the "
                f"combined evidence. Treat the resulting post-test probability as optimistic."
            )
    pre_prob = odds / (1 + odds)
    return dict(pre_prob=pre_prob, pre_odds=odds, used=used, warnings=warnings)


def post_test(pre_prob, x, model: LRModel, model_key, n_boot=4000, seed=1):
    """The headline function: pre-test odds x LR(GDF15) = post-test odds.

    Always uses the capped/floored, zone-flagged LR (never the raw uncapped lr()), and always
    propagates a bootstrap CI through to the post-test probability -- never a bare number.
    Symmetric on both tails: a rule-in cap above the data (§4.1) and a mirror-image rule-out
    floor below the data (no comparably extreme claim is defensible in either direction).
    """
    if not (0 < pre_prob < 1):
        raise ValueError("pre_prob must be strictly between 0 and 1")

    zone = extrapolation_zone(x, model_key)
    pre_odds = pre_prob / (1 - pre_prob)
    pvod_lo, pvod_hi = SUPPORT[model_key]['pvod']
    comp_lo, comp_hi = SUPPORT[model_key]['comp']

    def odds_to_prob(o):
        return o / (1 + o) if o is not None else None

    if zone == 'full_extrapolation':
        rep = reported_lr(x, model, model_key)
        bound = rep['reported']
        point_raw = None
        if x > pvod_hi:
            lr_point = lr_lo = bound
            lr_hi = None
        else:  # x < comp_lo
            lr_point = lr_hi = bound
            lr_lo = None
        capped = True
    elif zone == 'interpolation':
        point_raw, lr_lo, lr_hi = lr_ci(x, model, n_boot=n_boot, seed=seed)
        lr_point = point_raw
        capped = False
    elif x > comp_hi:
        # rule-in tail: clamp from ABOVE at the cap
        cap = RULE_OF_THREE[model_key]['cap']
        point_raw, lo, hi = lr_ci(x, model, n_boot=n_boot, seed=seed)
        lr_point, lr_lo, lr_hi = min(point_raw, cap), min(lo, cap), min(hi, cap)
        capped = point_raw > cap
    else:
        # rule-out tail: clamp from BELOW at the rule_out_floor
        floor = rule_out_floor(model)
        point_raw, lo, hi = lr_ci(x, model, n_boot=n_boot, seed=seed)
        lr_point, lr_lo, lr_hi = max(point_raw, floor), max(lo, floor), max(hi, floor)
        capped = point_raw < floor

    post_prob = odds_to_prob(pre_odds * lr_point)
    post_prob_lo = odds_to_prob(pre_odds * lr_lo) if lr_lo is not None else None
    post_prob_hi = odds_to_prob(pre_odds * lr_hi) if lr_hi is not None else None

    return dict(
        pre_prob=pre_prob, x=x, zone=zone, model_key=model_key,
        lr_point=lr_point, lr_ci=(lr_lo, lr_hi), lr_raw=point_raw,
        capped=capped,
        post_prob=post_prob, post_prob_ci=(post_prob_lo, post_prob_hi),
        delta_pp=100 * (post_prob - pre_prob),
    )


def _logit(p):
    return np.log(p / (1 - p))


def plot_fagan(scenarios, ax=None):
    """Fagan-nomogram-style plot: 3 parallel axes (pre-test %, LR, post-test %).

    Each scenario's connecting line is drawn from its own exact pre/post-test
    probabilities -- the line is mathematically exact by construction, not read off
    a shared middle scale. (A single shared, pretest-independent LR scale on the
    middle axis is not achievable with 3 equally-spaced straight parallel lines
    for an additive log-odds relationship -- the correct crossing point at x=1
    depends on the specific pretest value -- so rather than print a scale that is
    only exactly correct at one reference pretest, we label each line's own LR
    directly at its crossing point.)
    scenarios: list of dicts with pre_prob, lr, post_prob, label, color.
    """
    if ax is None:
        import matplotlib.pyplot as plt
        _, ax = plt.subplots(figsize=(7.5, 7))

    prob_ticks = [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999]
    for xpos in [0, 1, 2]:
        ax.axvline(xpos, color='#c9c9c6', linewidth=1, zorder=0)
    for p in prob_ticks:
        y = _logit(p)
        ax.text(-0.08, y, f"{p*100:g}", ha='right', va='center', fontsize=7.5, color='#52514e')
        ax.text(2.08, y, f"{p*100:g}", ha='left', va='center', fontsize=7.5, color='#52514e')
        ax.plot([-0.03, 0.03], [y, y], color='#8a8a86', linewidth=0.8)
        ax.plot([1.97, 2.03], [y, y], color='#8a8a86', linewidth=0.8)

    y_span = [_logit(0.001), _logit(0.999)]
    ax.set_xlim(-0.75, 2.75)
    ax.set_ylim(y_span[0] - 0.4, y_span[1] + 0.4)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['pre-test\nprobability (%)', 'likelihood\nratio', 'post-test\nprobability (%)'])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    for sc in scenarios:
        y1, y3 = _logit(sc['pre_prob']), _logit(sc['post_prob'])
        y2 = (y1 + y3) / 2  # exact crossing for THIS line; not a shared universal scale
        color = sc.get('color', '#2a78d6')
        ax.plot([0, 1, 2], [y1, y2, y3], color=color, linewidth=2, marker='o', markersize=5, zorder=2)
        ax.annotate(f"LR={sc['lr']:.1f}", xy=(1, y2), xytext=(1.08, y2 + 0.12),
                    fontsize=8, color=color)
        if sc.get('label'):
            ax.annotate(sc['label'], xy=(0, y1), xytext=(-0.72, y1),
                        fontsize=8, ha='left', va='center', color=color)
    return ax
