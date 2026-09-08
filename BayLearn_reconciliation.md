# BayLearn abstract — reconciliation against the final analysis

Every statement or number in the submitted abstract that the final analysis changes,
restates differently, or no longer supports.

**Source:** `BayLearn_Submission.pdf` — *"The Clinical Utility of a Promising Biomarker for Two
Overlapping Lung Diseases: What 26 Cases Support"*.

**The submitted abstract is not modified.** This is a record for the authors, prepared so the
extent of the change is visible before any decision about the abstract is taken.

Final-analysis values are from `results/results.json` unless another file is named. Bootstrap
endpoints are quoted to the nearest 10 pg/mL, as in the figures and captions.

---

## 1. Numbers that change

| § | As submitted | Final analysis | Why it changes |
|---|---|---|---|
| Findings | GDF15 AUC 0.956, **95% CI 0.912–0.990** | 0.956, **95% CI 0.912–0.989** | Same method, same primary status: the stratified percentile bootstrap, 4,000 replicates, is the reported interval. The upper endpoint differs by 0.001, which is inside the seed-to-seed range of the interval itself (0.9890 to 0.9904 over ten seeds, `results.json` → `auc_bootstrap_seed_sensitivity`). The final build fixes the seed at 1 and quotes what that seed gives. |
| §3 | "We refit the model 26 times, each time removing one PVOD case. Removing the patient with the highest GDF15 value changed the estimate at that extreme by **15%**." | **Reproduced, and larger than stated** — but only once the location is named. At 1,658 pg/mL the largest change is **−8.1%**; at 5,986 pg/mL, the highest concentration measured in the cohort, it is **+177.3%**. | The claim reproduces; what made it uncheckable was quoting one percentage without saying where on the curve it was measured. Both locations are now reported, with the case behind each. One correction to the abstract's wording: the case with the most influence at the top of the range is the **lowest** PVOD value (861 pg/mL), not the highest — removing it steepens the slope, which has its greatest leverage at the extremes. `results/leave_one_case_out.csv`. |
| §1 | "…logistic regression on log GDF15 with a quadratic term." | Linear only, **now justified rather than assumed**: likelihood-ratio test 1.84 on 1 df, **p = 0.175**; AIC 54.32 linear against 54.48 quadratic. | The quadratic term does not improve the fit, and the abstract gave no test either way. `results/functional_form.csv`. |
| Findings | At 1,658 pg/mL, 20% pre-test → **37%**, interval **19–46%** | LR **2.56 (1.40–4.74)** → **39.0% (25.9–54.2%)** | Different estimator. The abstract's figure came from the lognormal primary estimator; the final model is the single logistic on log10(GDF15). The qualitative point does not survive: 25.9% now sits above the 20% starting point, so the interval no longer contains where the clinician began. |
| Findings | At 3,000 pg/mL, 20% → **86%** (78–93%) | **93.3%** (84.0–99.7%) | Same cause. |
| Findings | Below 500 pg/mL, model-free LR **0.06** (0.05–0.08) | LR **0.005** (0.00002–0.028) | Same cause. Direction and conclusion unchanged; magnitude far stronger and the interval far wider, because the banded estimator's continuity correction floored the estimate. |
| Table 1 | NT-proBNP: **65** patients, AUC **0.558** | **63** patients, AUC **0.528** | The analysis uses the `NT-proBNP` distribution sheet (63; PVOD 15). The extra two are `Patient_ID` 24 and 26, which appear only in the PVOD-only correlation sheet. See §3. |
| Findings | Threshold reports **81% sensitivity and 98% specificity** | **80.8%** (21/26) and **97.1%** (136/140) | Recomputed against this comparator set. The abstract quotes the source letter's rounding; the exact counts are in `results.json` → `published_cutoff_performance`. |
| §3 | 0.632+ bootstrap changed the aggregate estimate by **0.0001** | Bootstrap optimism correction: AUC **0.9560 apparent → 0.9562 corrected** (optimism −0.0001); calibration slope **1.0000 apparent → 0.9676 corrected** | Different correction, and now two quantities rather than one. The abstract's headline claim — that the standard diagnostic barely moves — holds for the C statistic and is *contradicted* for the calibration slope, which moves by 3.2%. |

## 2. The 1,658 pg/mL example, in full

The abstract's sentence at this concentration is superseded. What replaces it is a
stronger version of the same argument, at the exact value the source letter recommends.
The same measurement is read twice — as a point on the fitted curve, and as a positive
result from the published dichotomy — and combined with the same pre-test probability, so
the reading is the only thing that differs.

| Pre-test probability | Curve: LR 2.56 (1.40–4.74) | Dichotomy: LR+ 28.3 |
|---|---|---|
| 10% | **22.1%** (13.4–34.5%) | **75.9%** |
| 20% | **39.0%** (25.9–54.2%) | **87.6%** |

The dichotomy's likelihood ratio is sensitivity / (1 − specificity) = (21/26) / (4/140) =
28.3, so a patient measured *at* the cutoff is handed the average evidence of everyone
above it, when their own value carries a factor of 2.56.

Both rows are needed because the gap depends on the clinician's starting suspicion, and a
figure quoted without its pre-test probability is not checkable — the same point the
leave-one-case-out result makes about location on the curve. `results.json` →
`models.IPAH+HPAH.published_cutoff_example`.

## 3. Statements the final analysis no longer supports

| § | As submitted | Position now |
|---|---|---|
| §1 | "We estimate the ratio three ways … a two-parameter lognormal … a model-free banded estimate … and logistic regression on log GDF15 with a quadratic term." | The final model is **one** estimator: logistic regression on log10(GDF15), entered linearly. Agreement-across-estimators is not available as evidence. The choice of the linear form is now justified by test rather than left implicit (§1). |
| §2 | "We can support only a floor: values in this region are **at least 33 times** more common in PVOD." | No rule-of-three floor is computed. The fitted curve gives LR 21.6 (10.1–248) at 2,500 pg/mL and 55.6 (20.9–1,563) at 3,000 — a point estimate with a wide band, not a bound. |
| §2 | "an inconclusive band of roughly **1,000 to 2,500** pg/mL" | No banded or zone construct exists. The comparable landmarks are LR = 1 at 1,384 pg/mL and LR = 10 at 2,160 pg/mL (95% CI 1,920–2,360). |
| §5 | The complete-case subset confirms the Table 1 comparator ranking (an internal claim of an earlier build, not of the submitted abstract) | **Withdrawn.** The 23 complete cases are 43.5% PVOD against 15.7% in the full cohort, because measurement availability is associated with diagnosis. A complete-case comparison on a differently-composed group cannot confirm the full-cohort ranking, and the comparator AUCs are conditional on the measurement being available (§5). |
| §1 | "Intervals come from 4,000 stratified bootstrap resamples." | True of the AUC, the LR = 10 landmark and the Youden comparison. The curve's own band is a **parametric** bootstrap; count-based quantities use exact binomial limits. |

## 4. NT-proBNP, stated in full

The abstract's 65 / 17 and the published 63 / 15 are both readings of the same workbook.

- The `NT-proBNP` distribution sheet covers all 166 patients and holds **63** measurements
  (PVOD 15, IPAH 33, HPAH 15). This is the analysis source and matches the published denominator.
- The PVOD-only correlation sheet `6MWD PVR mPAP NT-proBNP` holds **17** for the same 26 PVOD
  patients: the same 15, which agree exactly once the distribution sheet's median centring of
  1,264 ng/L is added back, plus **`Patient_ID` 24 and 26** (2,974 and 4,296 ng/L), which the
  distribution sheet records as NA.

Neither sheet is treated as the erroneous one. Per-patient detail is in
`results/ntprobnp_sources.csv`.

## 5. Added since the abstract

Two analyses in the final build have no counterpart in the submitted abstract.

- **Markers on a common set of patients, and why the subset does not settle anything.**
  Table 1's denominators run from 63 to 166, so no two rows describe the same patients. The
  clinical comparators are therefore recomputed on the 23 patients complete on all five
  clinical measurements, listed by `Patient_ID` in `results/complete_case_patients.csv`.
  This was expected to confirm the Table 1 ranking. **It does not, and is not reported as
  though it does.** The complete cases are 43.5% PVOD against 15.7% in the full cohort,
  about three times enriched, because measurement availability is associated with diagnosis:
  cardiac index is present for 96% of PVOD against 43% of IPAH and 36% of HPAH
  (`results/clinical_completeness.csv`). The subset is a different group of patients, not a
  smaller version of the study, and the rankings among the clinical variables move
  substantially on it (cardiac index rises from fourth to first, six-minute walk distance
  falls from first to third). Those movements are the symptom; the completeness table is the
  cause. The limitation carried is that comparator AUCs are conditional on the measurement
  being available, and availability is associated with diagnosis. This costs the abstract
  nothing, because the comparator table was never its argument.
- **Age-adjusted GDF15 on that subset, and raw GDF15's absence from it.** On the 23 complete
  cases age-adjusted GDF15 separates the two groups completely, running 1.05 to 2.64 in the
  10 PVOD patients against −0.06 to 1.00 in the 13 others. This is reported as complete
  separation and as two ranges, deliberately not as an AUC: the AUC would be 1.00, which is
  what the arithmetic gives when two small groups do not overlap, and printed beside the
  0.956 of Table 1 it would read as the marker improving on the subset when what changed is
  both the patients and the variable. Raw GDF15 is not in the table at all and cannot be:
  sheet `GDF15-1` is anonymous group-level concentrations with no patient identifier, so it
  cannot be keyed to a clinical row. `results/marker_auc_complete_case.csv`.
- **Blinding of outcome assessment.** Resolved in both directions and stated in the TRIPOD+AI
  checklist: diagnoses were established at the contributing centres under ESC/ERS criteria
  before and independently of GDF15 measurement, the assay was run blinded to group, and GDF15
  forms no part of those criteria, so the index test could not have contributed to the
  reference standard.

## 6. Unchanged — no action needed

| As submitted | Status |
|---|---|
| 166 of 322 patients; 26 PVOD, 140 idiopathic or heritable PAH | Unchanged |
| GDF15 AUC point estimate 0.956 | Unchanged, and the primary AUC |
| Table 1: GDF15 166 / 0.956, 6MWD 112 / 0.673, mPAP 157 / 0.617, PVR 142 / 0.602, cardiac index 82 / 0.567, age 166 / 0.522 | All six reproduce exactly |
| "Eighteen of the 26 PVOD patients exceed 2,500 pg/mL" | Confirmed, 18 of 26 |
| Healthy controls excluded from the diagnostic comparison | Unchanged |
| Threshold framing: above 1,658 suggests PVOD, below suggests IPAH | Unchanged |
| Complete-case, no imputation, denominators differing by variable | Unchanged as a method. The *interpretation* of the common-subset version has changed — see §3 |
| "26 PVOD cases from three registries with no independent test set; nothing here is validated for clinical use" | Unchanged |

Two near-misses worth a line rather than a row:

- "no comparison patient appears above roughly **2,000** pg/mL" — the exact comparator maximum is
  **1,809.5** pg/mL, and 20 of 26 PVOD patients lie above it. Figure 1 now labels that boundary
  directly, so the abstract's rounding and the figure would disagree in print.
- "the interval **19–46%** … still contains where they started" — see §1. Under the final model the
  interval no longer contains the 20% starting point, so the sentence's conclusion changes, not
  just its numbers.

## 7. The abstract's central argument, on the final build

The abstract pairs a near-zero optimism correction with a large leave-one-case-out shift, and
argues that aggregate validation and individual support are different claims. **Both halves now
hold on the final analysis**, and the pairing is sharper than the abstract stated:

- The C statistic moves by 0.0001 under bootstrap optimism correction — the aggregate
  diagnostic is as stable as the abstract said.
- One PVOD case moves the likelihood ratio by 8.1% at the published cutoff and by 177.3% at the
  top of the observed range. Where the two groups were both observed, one patient is worth a few
  per cent; where no comparator patient was ever measured, one patient is worth a factor of
  nearly three.

The argument is therefore available in the form the abstract intended, with the location on the
curve supplied. One qualification the abstract does not anticipate: the calibration slope moves
by 3.2% under the same correction, so "the standard diagnostic barely moves" is true of
discrimination and not of calibration.
