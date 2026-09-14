# GDF15 as continuous evidence for pulmonary veno-occlusive disease

Plasma GDF15 measured in a patient with pulmonary hypertension is a continuous quantity, and
how strongly a given concentration argues for pulmonary veno-occlusive disease (PVOD) depends
on the value measured. This analysis fits a single-predictor logistic regression on
log10(GDF15) and converts it into a **likelihood-ratio curve**: for any measured concentration
it returns the factor by which that measurement multiplies the odds of PVOD, with a confidence
band, an explicitly marked range over which the curve is supported by observation, and the same
curve fitted separately against each comparator diagnosis.

The likelihood ratio is used rather than a predicted probability because the cohort was
assembled with PVOD cases deliberately over-represented (15.7% of the analysis sample). A
predicted probability from such a sample applies only to a population with that same case mix.
A likelihood ratio does not: it is a property of the two conditional distributions of GDF15
given disease status, so a clinician can combine it with their own pre-test probability.

## QC status

The published-data analysis uses Giorgio's integrated three-panel Figure 1 with 4,000 requested stratified resamples: 3,991 finite ordinary logistic fits and 9 separated resamples excluded. The upper edge is effectively unchanged near/below the cutoff and narrows in the tail. See [freeze report](qc/published_data_freeze.md), [clean reproduction](qc/published_data_clean_run.json), and [integration record](qc/figure1_integration.md). The freeze tag is `analysis-v1.0-published-data`; no remote is configured and nothing has been published.

## Quick start

Use Python 3.14 (tested with 3.14.6 on macOS arm64). The only data input is the published supporting workbook, bundled unchanged at `archive/gdfdiscrim/jci-136-199159-s068.xlsx`. The private reconstructed workbook is neither required nor distributed.

```sh
python3.14 -m venv venv
venv/bin/python -m pip install -r requirements-lock.txt
venv/bin/python analysis/workflow.py --verify
```

This reconstructs clinical records into ignored `generated/`, fits the existing models, writes `results/`, generates Table 1 and existing figures/captions, and runs independent checks plus two deterministic temporary reruns. No network is required after environment installation. The new Figure 1 remains gated; this command does not claim to regenerate Giorgio's PNG.

A separate source workbook and output directory can be supplied:

```sh
venv/bin/python analysis/workflow.py --workbook path/to/published.xlsx --output-dir output --verify
```

Individual steps are also available:

```sh
venv/bin/python analysis/reconstruct.py --workbook path/to/published.xlsx --output-dir generated
venv/bin/python analysis/run_analysis.py --workbook path/to/published.xlsx --output-dir output
venv/bin/python analysis/make_figures.py --workbook path/to/published.xlsx --results-dir output/results --figures-dir output/figures
venv/bin/python analysis/verify.py --workbook path/to/published.xlsx --results-dir output/results --report output/results/verification.json
```

The `requirements.txt` compatibility ranges remain available; the lock file is the tested environment for exact reproduction. PDF/PNG byte determinism is tested within that environment; different fonts, platforms, or library versions may change figure bytes without changing numbers.

## Repository layout

| Path | Contents |
|---|---|
| `archive/gdfdiscrim/jci-136-199159-s068.xlsx` | Unmodified published input; active source exception within archive |
| `Supporting data values_...xlsx` | Historical cleaned copy; numeric agreement verified, original untouched |
| `analysis/` | Reconstruction, analysis, figure/table generation and independent verification |
| `generated/` | Ignored local clinical records and provenance; do not publish |
| `results/` | Model outputs, generated Table 1, anonymous source-cell verification |
| `figures/` | Reproducible existing figures and captions; new Figure 1 pending |
| `qc/` | Baseline, audits, source determinations, QC report and supplied figure |
| `TRIPOD-AI_mapping.md` | Reporting checklist with author-owned items still pending |
| `BayLearn_reconciliation.md` | Submitted abstract versus current analysis and post-submission notes |
| `archive/` | Historical exploratory material; clinical exports/notebooks stay local |

## Intended use

The model is intended for patients under evaluation for pulmonary hypertension in whom PVOD is
already part of the differential diagnosis, before PAH-targeted therapy is started. Its role is as
an add-on test that updates pre-test suspicion alongside HRCT, DLCO, EIF2AK4 genotyping and other
clinical information. It is not a screening test, a stand-alone diagnostic, or an independent
rule-in / rule-out test. The decision it informs is whether the evidence for PVOD is sufficient to
affect the safety assessment before starting pulmonary vasodilator therapy. The likelihood ratio is
the model output; a post-test probability requires the clinician's own pre-test probability.
`TRIPOD-AI_mapping.md` carries this in full.

## Input and reconstruction

The published supporting data accompany Prabhakar et al., *J Clin Invest* 2026;136(2):e199159, DOI [10.1172/JCI199159](https://www.jci.org/articles/view/199159). `analysis/reconstruct.py` maps original `Fig_1*`/`Fig_S*` sheet names to the historical cleaned aliases used below. It accepts either workbook. All 4,609 numeric cells agree; the cleaned workbook adds diagnosis labels absent from original `Fig_1g`.

Raw GDF15 comes from `Fig_1a` (`GDF15-1`): Ctrl 68, PVOD 26, IPAH 98, HPAH 42, Other PH 52, COPD 36. No patient key exists for these concentrations.

Clinical IDs, age, age-adjusted GDF15 and measurements are rebuilt from `Fig_1g`–`Fig_1m` and `Fig_S2b`. The diagnosis block convention and every join, median offset, source rounding rule and missingness distinction are documented in [reconstruction provenance](qc/reconstruction.md). The generated records are local artifacts. A local comparison against the historical private workbook passed all 22 columns across both data tabs.

No model quantity is taken from the private workbook. The reconstructed clinical numeric values are also not used for marker AUCs: those use the published distribution values. Reconstructed NT-proBNP has two negative artifacts, which are preserved for provenance and never interpreted as concentrations.

### Patient_ID

`Patient_ID` runs 1–166, every value unique, and it is a stable local identifier for the reconstructed clinical rows. Age cannot — 117 distinct ages over 166 patients, nine
of them 69. Group cannot. Age-adjusted GDF15 happens to be unique, but it is a measurement,
and measurements get transformed, rounded and filtered.

It is carried through every patient-level table and every intermediate frame. Report a
subset, an outlier, or a dropped case by `Patient_ID`, never by row number: row numbers
change the moment anyone sorts or filters, and then two people are quietly comparing
different patients.

Two cautions. The ids are blocked by diagnosis — 1–26 PVOD, 27–124 IPAH, 125–166 HPAH — so
an id implies a group, and a figure ordered by `Patient_ID` would show a pattern that is
only that blocking. And it is an identifier created for this reconstruction, not a clinical
one from the contributing centres, so anything reaching the supplement that uses it needs a
line saying so.

`Patient_ID` does **not** reach anything derived from `GDF15-1`. That sheet is anonymous, so
`leave_one_case_out.csv` identifies its 26 cases by source-sheet position and GDF15 value.
A `case_id` there is not a patient and has no relation to a `Patient_ID`.

## Dependencies

Python 3.9 or newer and the packages in `requirements.txt`:

```
python -m venv venv
venv/bin/pip install -r requirements.txt
```

## Running it

```
venv/bin/python analysis/run_analysis.py  # the analysis, its tables, and the figures
venv/bin/python analysis/verify.py        # optional: numerical checks on the implementation
```

`run_analysis.py` is the whole workflow: it reads the workbook, fits every model, computes
every reported value, writes them to `results/`, and then calls `make_figures.py`, which draws
the figures and writes their captions. It takes about a minute, most of it the bootstraps and
the seed-sensitivity records described below. `make_figures.py` can also
be run on its own to redraw the figures from the values already in `results/`.

## What is generated

Everything below is written by the workflow; nothing else in the repository is modified.

`results/`

| File | Contents |
|---|---|
| `results.json` | every value the figures, captions and manuscript use — the authoritative source |
| `model_coefficients.csv` | intercept, log10(GDF15) coefficient, standard errors, events per variable, AUC, for every model |
| `comparator_models.csv` | the GDF15 concentration at which the likelihood ratio reaches 10, per comparator |
| `lr_reference.csv` | likelihood ratio and post-test probability at reference concentrations |
| `lr_curves.csv` | the plotted curves with their confidence bands |
| `calibration.csv` | observed versus predicted PVOD frequency by quintile |
| `marker_auc.csv` | GDF15 against every other marker measured in the same 166 patients |
| `auc_reconciliation.csv` | every AUC printed for these data, recomputed, with the calculation used |
| `internal_validation.csv` | apparent and optimism-corrected C statistic and calibration slope |
| `leave_one_case_out.csv` | the 26 refits, one per PVOD case removed, with the likelihood ratio each gives |
| `functional_form.csv` | log10(GDF15) against log10(GDF15) plus a quadratic term |
| `marker_auc_complete_case.csv` | the clinical comparators recomputed on the complete-case subset |
| `complete_case_patients.csv` | the 23 `Patient_ID`s making up that subset |
| `clinical_completeness.csv` | which clinical measurements exist, by diagnosis |
| `ntprobnp_sources.csv` | NT-proBNP per PVOD patient, in each of the two source sheets |

`figures/` — `figure1_lr_curve`, `figure2_comparators`, `figure3_calibration`, each as PDF and
PNG, each with a caption in Markdown whose numbers are read from `results.json`.

## How the likelihood ratio is calculated

The model is

```
logit P(PVOD | x) = b0 + b1 * log10(GDF15)
```

fitted by maximum likelihood to PVOD (26) versus pooled IPAH and HPAH (140). The fitted
intercept absorbs the case mix of that sample, so it is removed:

```
log LR(x) = [b0 + b1 * log10(x)] - log(pi / (1 - pi))
```

where `pi` is the proportion of PVOD cases in the fitted sample (26/166). `LRModel.lr` in
`analysis/model.py` implements exactly this. A clinician's own suspicion then enters through

```
pre-test odds  x  LR  =  post-test odds
```

implemented as `model.post_test`. From a 10% pre-test probability, 1,200 pg/mL carries LR 0.48
and gives a post-test probability of 5.0%; 1,700 pg/mL carries LR 2.9 and gives 24.4%.

### The published cutoff, read two ways

At 1,658 pg/mL — the cutoff reported for these data — the two readings of the same
measurement diverge, and by how much depends on the clinician's own suspicion:

| Pre-test probability | Curve: LR 2.56 (1.40–4.74) | Dichotomy: LR+ 28.3 |
|---|---|---|
| 10% | 22.1% (13.4–34.5%) | 75.9% |
| 20% | 39.0% (25.9–54.2%) | 87.6% |

The dichotomy's likelihood ratio is sensitivity / (1 − specificity) = (21/26) / (4/140),
so a patient measured *at* the cutoff is handed the average evidence of everyone above
it. The curve reads that patient's own value, which at the cutoff is close to
uninformative. Both readings use the same patients and the same outcome; they differ
only in what is done with the measurement. `results.json` →
`models.IPAH+HPAH.published_cutoff_example`.

## Functional form of the predictor

log10(GDF15) is entered linearly. Adding a quadratic term does not improve the fit:
likelihood-ratio test 1.84 on 1 degree of freedom, p = 0.175, with the AIC slightly
favouring the linear model (54.32 against 54.48). The linear model is the one reported.
`results/functional_form.csv`.

## Influence of a single case

The model is refitted 26 times, each time with one PVOD case removed, and the likelihood
ratio is read at the same two concentrations in every refit. How far one patient moves
the answer depends entirely on where on the curve the answer is read:

| Read at | LR on all 26 cases | Largest change on removing one case |
|---|---|---|
| 1,658 pg/mL (the published cutoff) | 2.56 | −8.1%, removing the case at 1,046 pg/mL |
| 5,986 pg/mL (the highest value measured) | 2,012 | +177.3%, removing the case at 861 pg/mL |

Within the range where both groups were observed, one patient is worth a few per cent.
At the top of the range, where no comparator patient was ever measured, one patient is
worth a factor of nearly three. This is why a change of this kind cannot be quoted
without saying where on the curve it was measured. Note that the case with the most
influence at the top of the range is the *lowest* PVOD value, not the highest: removing
it steepens the fitted slope, and slope has its greatest leverage at the extremes.

This is a sensitivity analysis, not a validation procedure — the reported model remains
the fit on all 26 cases. `results/leave_one_case_out.csv` holds all 26 refits.

## Markers on a common set of patients

`marker_auc.csv` lets every marker use all the patients who have it, so its denominators
run from 63 to 166 and no two rows describe the same patients.
`marker_auc_complete_case.csv` recomputes the clinical comparators on the 23 patients who
have all five clinical measurements. `complete_case_patients.csv` lists those 23 by
`Patient_ID`, so the exact subset can be reproduced without re-deriving it:

> 1, 3, 4, 5, 6, 7, 10, 16, 22, 25 (PVOD) · 47, 58, 99, 107, 110, 111, 117 (IPAH) ·
> 130, 159, 160, 161, 164, 165 (HPAH)

NT-proBNP here is the distribution sheet, the source used throughout the analysis. Taking
the union with the PVOD-only correlation sheet instead would give 25 patients rather than
23. No value is imputed and nothing is modelled about why a value is missing.

**This subset does not confirm the full-cohort ranking, and is not reported as though it
does.** It is 43.5% PVOD against 15.7% in the full cohort — about three times enriched. It
is not a smaller version of the study but a different group of patients, and the rankings
among the clinical variables move substantially at this size: cardiac index rises from
fourth to first (0.567 to 0.792) and six-minute walk distance falls from first to third
(0.673 to 0.712). Those movements are the symptom. The cause is in the next table.

### Why the complete cases are not representative

Whether a patient has a measurement depends on their diagnosis, and it depends on it most
for the two variables missing most often. PVOD patients are investigated more thoroughly
and arrive with fuller data. From `clinical_completeness.csv`:

| Variable | PVOD | IPAH | HPAH |
|---|---|---|---|
| Cardiac index | 96% | 43% | 36% |
| NT-proBNP | 65% | 34% | 36% |
| 6MWD | 85% | 64% | 64% |
| PVR | 88% | 86% | 83% |
| mPAP | 92% | 95% | 95% |

So the limitation to carry is this one:

> Comparator AUCs are conditional on the measurement being available, and measurement
> availability is associated with diagnosis.

A complete-case comparison cannot be treated as an unbiased re-estimate of the full-cohort
comparator ranking. That is a better limitation than the one the analysis set out to check,
and it costs nothing, because the comparator table was never the argument of the paper.

The completeness table counts a patient as measured if *either* NT-proBNP sheet records
them — 17 of 26 PVOD — because the question there is whether the patient was measured at
all. The complete-case subset uses the distribution sheet alone, 15 of 26, because the
question there is which patients share one consistent analysable set. The two are
different questions and are reported separately.

### GDF15 on this subset

The table holds the clinical comparators only. Age-adjusted GDF15 is reported here in a
sentence instead, deliberately:

> On the 23 complete cases, age-adjusted GDF15 separates the two groups completely. It
> runs from **1.05 to 2.64** in the 10 PVOD patients and from **−0.06 to 1.00** in the 13
> IPAH and HPAH patients, with a clean gap between them.

That is stated as complete separation and as two ranges, not as an AUC. The AUC would be
1.00, which is what the arithmetic gives when two small groups happen not to overlap — with
10 patients against 13 that is not rare — and it measures the absence of overlap rather
than the performance of the marker. Printed beside Table 1's 0.956 it would read as the
marker improving on the subset, when what changed is both the patients *and* the variable.

**Raw GDF15 cannot be in this table at all.** Sheet `GDF15-1` is six anonymous cohort
columns of unequal length, with no patient identifier, no age and no row order shared with
the per-patient sheets. Neither join mechanism reaches it, and no key can be constructed
without guessing which patient goes with which value. This is the same limitation that
prevents a combined GDF15-plus-clinical model and confines the main analysis to GDF15 with
the diagnosis label only.

## Internal validation

The model is fitted and evaluated on the same 166 patients, so both reported performance measures
carry an optimism correction beside them (`results/internal_validation.csv`). In each of 4,000
bootstrap replicates the model is refitted on a resample of the patients and scored twice, on the
resample and on the original data; the mean difference is the optimism, and the corrected value is
the apparent value minus it. Two quantities are corrected: the C statistic, and the calibration
slope, which is 1 by construction in the sample a model was fitted to, so the corrected value is
what the resampling has to say about it. Optimism is not constrained in sign — with a single
strongly separating predictor the C-statistic optimism comes out very slightly negative, so the
corrected C statistic sits marginally above the apparent one. Nothing beyond this correction is
applied as validation: no cross-validation and no shrinkage. The leave-one-case-out refitting
described above is a sensitivity analysis reported beside it, not a second validation procedure,
and it does not alter the reported model.

Figure 1 uses a stratified case-bootstrap band (seed 20260907), with 3,991 finite fits from 4,000 requested draws and 9 complete-separation exclusions. These are ordinary logistic MLEs; no second estimator is introduced. `results/figure1_band.csv` and `results/figure1_bootstrap.json` are the production band and accounting. The original parametric intervals in `lr_curves.csv` and `results.json` remain supplementary outputs, and Figure 2 retains its existing parametric bands.

For those retained **parametric bootstrap** intervals: outcomes are simulated
from the fitted model at the observed GDF15 concentrations and the model is refitted 4,000
times. Holding the concentrations fixed and simulating only the outcomes lets the interval be
defined across the whole measured range, including concentrations at which one class was never
observed. Count-based quantities — sensitivity, specificity, the observed proportion in each
calibration quintile — use exact (Clopper–Pearson) binomial limits instead.

## Reproducibility notes

- **Verified environments.** The workflow was run from a clean checkout on Python 3.14.6
  (numpy 2.5.1, pandas 3.0.5, scipy 1.18.0, statsmodels 0.14.6, matplotlib 3.11.1) and, from a
  fresh virtual environment built from `requirements.txt`, on Python 3.9.6 (numpy 2.0.2,
  pandas 2.3.3, scipy 1.13.1, statsmodels 0.14.6, matplotlib 3.9.4). Tables, captions and
  figures agree between them, with one exception: a single exact-binomial limit differs in the
  17th significant digit between scipy versions. On Python 3.14.6 with statsmodels 0.15.0
  (numpy 2.5.2, scipy 1.18.1, pandas 3.0.5, matplotlib 3.11.1) the figures are byte-identical
  and every reported value agrees, with fitted coefficients and their standard errors
  differing in the 12th significant digit and beyond — the optimiser's tolerance, not the
  analysis.
- **Seeds.** Existing model/validation calculations use seed 1 (`model.SEED`) and 4,000 replicates
  (`model.N_BOOT`); Figure 1 uses seed 20260907 and the approved 3,991/9 accounting.
  Two clean runs in the pinned environment reproduce all generated output files
  byte for byte identical, figures included; `verify.py` asserts this by running the whole
  workflow twice into temporary directories and comparing every file. A percentile bootstrap on 26 cases still
  moves in the third decimal between seeds, so any bootstrap interval should be quoted with its
  seed. `results.json` records that seed sensitivity for the AUC under
  `auc_bootstrap_seed_sensitivity`, and for the threshold comparison under
  `threshold_comparison.seed_sensitivity`. Over ten seeds the LR = 10 endpoints move by about
  7 pg/mL each and the width ratio between 3.5 and 3.6, while the Youden endpoints do not move at
  all: a selected cutoff is a discrete choice among observed values, and the percentile lands on
  the same two of them every time.
- **Paths.** All paths are resolved relative to the repository root; none is absolute.
- **NT-proBNP.** Two sheets of the workbook carry this variable and they support different
  denominators for PVOD. The `NT-proBNP` distribution sheet covers all 166 patients and holds 63
  measurements (PVOD 15, IPAH 33, HPAH 15); it is the source for the analysis, and its denominator
  is the one published for this variable. The PVOD-only correlation sheet
  `6MWD PVR mPAP NT-proBNP` holds 17 measurements for the same 26 PVOD patients: the same 15,
  which agree exactly once the distribution sheet's median centring of 1,264 ng/L is added back,
  and two further patients — `Patient_ID` 24 and 26, at 2,974 and 4,296 ng/L — that the distribution
  sheet records as NA. Neither sheet is treated here as the erroneous one. Both counts are
  asserted in `analysis/data.py`, the per-patient comparison is written to
  `results/ntprobnp_sources.csv`, and `results.json` records the difference under
  `ntprobnp_denominator`.
- **Row keying and `Patient_ID`.** The clinical sheets of the source workbook carry no patient
  identifier. Each row is keyed to its cohort label through the age-adjusted GDF15 value, and the
  keying is checked end-to-end by requiring raw and median-centred 6MWD to differ by exactly
  294 m on all 112 patients who have that measurement. `Patient_ID` is then joined on that same
  age-adjusted GDF15 value, which is distinct across all 166 patients in both files; the join is
  asserted to be one-to-one, to match exactly rather than to a tolerance, and to agree with the
  `6MWD-1` diagnosis label on every row. **The join is by key and never by position.** The
  clinical frame is in the distribution sheets' order, which runs HPAH, then IPAH, then PVOD;
  `Patient_ID` is in the `6MWD-1` order, which runs PVOD first. Lining the two up by row would
  silently mislabel all 166 patients, and `verify.py` asserts that the frame is *not* in
  `Patient_ID` order precisely so that this cannot pass unnoticed.
- **AUC.** The reported interval is a **stratified percentile bootstrap** — PVOD and
  comparator patients resampled separately at their observed sample sizes, 4,000
  replicates, seed 1 — giving 0.956 (95% CI 0.912–0.989). It is preferred because it is
  bounded by construction: every replicate's AUC is itself an AUC, so no endpoint can
  fall outside [0, 1]. Two other intervals are recorded beside it and never combined with
  it. The first is this repository's asymptotic C-statistic interval, 0.91605–0.99604,
  computed by `model.auc_ci_asymptotic` from the structural-component variance
  `v10.var(ddof=1)/m + v01.var(ddof=1)/n` and checked in `verify.py` against the same
  quantities computed directly from the definition; it is a normal interval on the AUC
  scale that is **clipped** to [0, 1] afterwards rather than bounded by construction, so
  an endpoint printed as 1.000 anywhere in `auc_reconciliation.csv` has been truncated. The
  same applies to the age-adjusted GDF15 row of `marker_auc.csv`, whose upper endpoint prints
  as 1.0 for the same reason. (`marker_auc_complete_case.csv` no longer carries a GDF15 row at
  all; the complete separation it would have reported as an AUC of 1.00 is given as two ranges
  under "Markers on a common set of patients".)
  The second is the source study's published Figure 1E interval, 0.9164–0.9957, quoted
  for comparison only. The first two are computed here; the third is not. The full
  derivation is in `results.json` → `auc_ci_provenance`.
- **The LR = 10 landmark.** Its reported interval is a stratified case bootstrap: PVOD and
  comparator patients resampled separately at their observed sample sizes, 4,000 replicates,
  seed 1, with the landmark and the Youden cutoff both recomputed inside every replicate so the
  two are compared on identical footing. The parametric-bootstrap interval for the same landmark
  is kept beside it as a supplementary sensitivity analysis. Figure 1's pointwise band uses the separately recorded finite-fit stratified bootstrap. Reader-facing text rounds bootstrap endpoints to the nearest
  10 pg/mL; `results/` keeps full precision.
- **Source-figure values.** The three Figure 1D values (94.94, 97.26, 95.19) were read from the
  published figure and reproduce from these supporting data to two decimal places. The Figure 1E
  values — 95.60 (91.64–99.57) for PVOD versus IPAH + HPAH and 96.05 (92.45–99.66) for PVOD versus
  IPAH + HPAH + Other PH — are those a co-author of the source study reports from the 600-dpi
  source render; both reproduce from the supporting data. This analysis has not inspected that
  render at that stage; the first-author QC brief now confirms independent verification of both Figure 1E values, recorded under
  `source_figure_verification` in `results.json`. Nothing here claims to correct the source study.
- **Perfect separation.** PVOD and healthy controls do not overlap in GDF15 (highest control
  684.1 pg/mL, lowest PVOD 860.5 pg/mL), so that comparison has no finite maximum-likelihood
  coefficient and no curve is fitted for it.
- **Not validated for clinical use.** Single cohort, no external validation set; the analysis
  reports the evidential weight available in these data and nothing beyond it.

## Confirmed declarations and source checks

Jordan's secondary work was supported through the **2026 Responsible AI Summer Research Program at Santa Clara University**, under **“Transparent AI for Rare Disease Diagnosis”**. The programme had no role in the analysis. Original-cohort funding remains separate; see `qc/sources.md`. The study is a post hoc secondary analysis, with no prospective protocol and no registration; the model form was selected after seeing the data. There was no patient/public involvement, and the proposed presentation has not been tested with clinicians or patients.

A1, A2, A4–A7, the separate SCU determination, and A9 publication details remain author-owned. A3 source text and A8 are verified; neither resolves the other author items. No remote or Zenodo release has been created.

## Negative likelihood ratios

At the fixed published cutoff, LR− = (5/26)/(136/140) = 0.1979638009, displayed as 0.20. Two historical cross-validated classifier PDFs display 0.24 correctly for a different procedure: (6/26)/(136/140) = 0.2375565611. Do not substitute one for the other. The complete occurrence ledger is `qc/audit_024.csv`.

## Published-data freeze and future validation

Run a fresh environment with `requirements-lock.txt`, then `python analysis/freeze_reproduce.py --python /absolute/path/to/fresh/venv/bin/python`. This copies only analysis code and the original published workbook into an isolated directory, runs the full workflow with verification, repeats it, compares all artifact bytes, and promotes successful outputs. The exact supplied Figure 1 script is preserved in `qc/figure1_bootstrap/giorgio_original.py`; root `make_figure1_v2.py` now delegates to production and accepts either workbook.

Patient-level reconstructions and exports remain ignored local artifacts. Five older tracked exports have been removed from the frozen tree; they still exist in earlier local history. The private reconstructed workbook was never tracked. Nothing has been pushed. Jordan must create or approve a private remote before any push; no public repository or Zenodo release is authorized.

`data.py` loads development inputs separately from `model.py` fitting/scoring. Its cohort assertions intentionally lock this published development cohort. A future validation source can have its own loader and be scored using frozen coefficients, without altering or pooling into this development fit. No future samples were analyzed or simulated and no v2.0 was created. The participant-count disagreement and assay/blinding/author-order confirmations remain author-side submission items; they do not change this published-workbook analysis.
