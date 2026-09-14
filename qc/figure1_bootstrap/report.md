# Figure 1 bootstrap reconciliation — integration stopped

> Historical audit, superseded by Giorgio's approval of 11 September 2026. The finite-fit band is now integrated into production. The earlier stop decision and numerical evidence below are retained as provenance; they are not current blockers. See [current integration](../figure1_integration.md) and [freeze report](../published_data_freeze.md). The audit command now exits successfully when the approved 3,991/9 accounting reproduces.

The supplied-script gate is resolved. The uncertainty-band gate remains open: nine strictly separated bootstrap samples have no finite unpenalized logistic maximum-likelihood estimate. Removing these from the percentile calculation changes the band materially. Production code, results, figures and both manuscripts remain unchanged.

## Exact resamples and solver evidence

The original recipe is archived unchanged in `giorgio_original.py`. `analysis/audit_figure1_bootstrap.py` reproduces its BFGS objective, initialization and default settings. It uses the published workbook through the verified shared data reader, ordering PVOD, IPAH and HPAH exactly as Giorgio does. Seed 20260907 and 4,000 stratified draws reproduce all 59 precision-loss warnings.

`warning_fits.json` contains all 59 warning fits, four additional separated fits with original success status, and one additional inaccurate successful fit. Each record includes its zero-based replicate index, RNG state before sampling, all 166 source-array indices, sample SHA-256, original coefficients, objective, numerical and analytic gradient norms, optimizer status, and independent solver diagnostics. Indices refer to anonymous public raw observations, not reconstructed clinical patient IDs. `all_fits.json` records every draw and identifies the exact accepted set.

## Deterministic fit policy

Every draw uses standardized log10 concentrations and independent statsmodels Newton logistic maximum likelihood, with 100 iterations and parameter tolerance 1e-11. Acceptance requires overlapping classes, full-rank design, solver convergence, finite coefficients, positive-definite observed information, and maximum absolute standardized analytic score below 1e-7. The resulting unique stationary point is the finite maximum of the concave logistic likelihood. No penalization, slope cap, replacement draws or forced estimates are used.

Separated draws are rejected with an explicit gap between the lowest case and highest comparator. For all nine, an additional independent statsmodels GLM/IRLS attempt is saved as diagnostic evidence; finite terminal iterates, even with a convergence flag, cannot overcome the mathematical nonexistence of a finite MLE. No genuine Newton failures or other non-identifiability occurred among the 3,991 overlapping draws.

The numerical comparison tolerance was set before inspecting outcomes: a 0.1% maximum relative LR difference across the original 400-point 700–6,500 grid plus the eight requested concentrations distinguishes A from B. This is an audit tolerance, not a clinical significance threshold.

## Disposition and counts

- 59 warning fits: **54 A (valid despite warning), 0 B, 5 C (strictly separated)**.
- All 4,000 samples: **3,990 A, 1 B, 9 C**. Four of the separated samples had a misleading original success flag.
- Requested: **4,000**. Accepted finite refits: **3,991**. Rejected: **9**.
- Independent refits/diagnostic attempts: **4,000**. All **3,991** accepted fits use Newton coefficients consistently; **54** warning cases are resolved and **1** additional successful BFGS fit is numerically repaired.
- Zero-based replicate 2296 is the additional B fit: maximum relative LR difference 0.1913%, while the negative log likelihood improves from 3.248630491618 to 3.248630489954. Maximum warning-case relative LR difference among finite fits is 0.03783%.
- Separated replicate indices: **149, 769, 1111, 1134, 1392, 1654, 1873, 3464, 3536**.

## Candidate finite-fit percentile band versus original

The candidate uses precisely the 3,991 accepted Newton fits, with 2.5th/97.5th percentiles on the LR scale and the fixed sampling log-odds log(26/140). It is conditional on a finite fit and has not been adopted as production uncertainty. Original percentiles include all 4,000 BFGS outputs, including the nine invalid finite approximations to separated fits.

| GDF15 (pg/mL) | Original LR interval | Candidate finite-fit interval | Lower change | Upper change |
|---:|---:|---:|---:|---:|
| 861 | 0.0029188–0.254636 | 0.00346551–0.254768 | +18.7308% | +0.0516% |
| 1658 | 1.23895–4.13077 | 1.26584–4.13096 | +2.1704% | +0.0046% |
| 1809.5 | 2.40747–6.92314 | 2.44344–6.92478 | +1.4942% | +0.0236% |
| 2156 | 6.60283–25.9717 | 6.60431–25.5306 | +0.0224% | -1.6983% |
| 2500 | 12.6908–93.6201 | 12.6871–87.6967 | -0.0293% | -6.3270% |
| 3000 | 26.5606–502.753 | 26.5597–469.617 | -0.0034% | -6.5909% |
| 4000 | 82.929–7578.81 | 82.8999–6473.98 | -0.0351% | -14.5778% |
| 5986 | 406.606–316161 | 405.896–270854 | -0.1747% | -14.3304% |

Maximum absolute endpoint difference: **85,364.923722 LR units**, upper endpoint at 6,500 pg/mL. Maximum relative endpoint difference: **24.808277%**, lower endpoint at 711.827639 pg/mL. The CSV retains both absolute and relative differences at all grid/key concentrations.

The change is material to the plotted uncertainty, including manuscript-visible interval rounding (at 1,658, lower LR endpoint 1.24 becomes 1.27). At the eight requested concentrations, whether the interval lies above/below LR=1 or LR=10 is unchanged; this does not justify claiming the band itself is unchanged. The main-model coefficients, pointwise LR curve, patient probabilities, classification counts, AUC, calibration, LOO and missingness calculations were not modified.

## Verification and withheld work

**81 focused audit checks, 0 failures.** Checks replay each warning sample from its saved RNG state, independently confirm all nine separation gaps, reconstruct the accepted-set percentile band, and repeat every validated fit with byte-identical coefficients and identical accepted/rejected sample hashes. The shared reader gives identical 322 values and diagnoses for published Fig_1a and cleaned GDF15-1. This verifies audit-reader compatibility, not an integration of the supplied standalone plotting script.

Reproduce the evidence with `venv/bin/python analysis/audit_figure1_bootstrap.py` (expected nonzero exit at the explicit material-band-change guard), then `venv/bin/python analysis/verify_figure1_audit.py` (expected 81 checks passed). The guard is a scientific integration blocker, distinct from failures of the audit verification checks.

Per the requested stop condition, Figure 1 production integration, annotation changes, support-marker changes, QC DOCX edits, full-workflow regression and release tagging are withheld. No fresh-environment full-workflow claim is made. No commit, tag, push, publication or private-workbook tracking occurred. No modification of the original manuscript occurred.

The next decision is how to document the nine separated resamples and the resulting finite-fit-conditioned interval. The original all-4,000 finite-parameter band cannot be certified as a band of valid logistic MLEs. No replacement scientific procedure has been silently adopted.
