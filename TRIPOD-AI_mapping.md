# TRIPOD+AI checklist — QC revision; author items pending

Checklist: TRIPOD+AI (Collins et al., *BMJ* 2024;385:e078378), 27 items. This study
**develops** a diagnostic prediction model and does not evaluate it on separate data, so
items that apply only to evaluation are marked N/A with the reason given.

This is the checklist as it will be submitted. Each item states what the study actually
did, and where in the manuscript it is reported. **Manuscript locations are section
placeholders**, because the manuscript is in draft and has no page numbers yet; replace
each with a page number at submission. The section scheme is the one used in the draft:

| Placeholder | Section |
|---|---|
| Methods §2.1 | Data source and participants |
| Methods §2.2 | Outcome and predictor |
| Methods §2.3 | Model specification |
| Methods §2.4 | Statistical analysis |
| Methods §2.5 | Missing data |
| Results §3.1 | Participants |
| Results §3.2 | Fitted model |
| Results §3.3 | Performance and internal validation |
| Results §3.4 | Sensitivity analyses |
| Discussion §4.1–§4.3 | Interpretation, limitations, implications |

Funding, protocol, registration and patient/public involvement are confirmed below. This checklist retains the manuscript's author-review status; the current author-item summary is in `qc/published_data_freeze.md`. The historical issue tracker remains at `analysis-v1.0-published-data:qc/hostile_review.md`.

---

## Intended use

The framing used throughout this checklist and in the manuscript.

**Target population.** Patients under evaluation for pulmonary hypertension in whom PVOD
is already part of the differential diagnosis, before PAH-targeted therapy is started.

**Intended role.** An add-on test that updates pre-test suspicion alongside HRCT, DLCO,
EIF2AK4 genotyping and other clinical information. It is not a screening test, a
stand-alone diagnostic, or an independent rule-in / rule-out test.

**Clinical decision informed.** Whether the evidence for PVOD is sufficient to affect the
safety assessment before starting pulmonary vasodilator therapy.

**Model output.** The likelihood ratio is the model output; a post-test probability
requires the clinician's own pre-test probability.

---

## Title and abstract

| Item | Checklist item | How this study addresses it | Location |
|---|---|---|---|
| 1 | Identify the study as developing or evaluating a multivariable prediction model, the target population, and the outcome | Development of a single-predictor logistic regression converting plasma GDF15 into a likelihood ratio for PVOD against IPAH and HPAH. Target population: patients under evaluation for pulmonary hypertension in whom PVOD is already in the differential, before PAH-targeted therapy. Outcome: the diagnosis of PVOD. | Title page |
| 2 | See TRIPOD+AI for Abstracts checklist | Completed on the separate TRIPOD+AI for Abstracts checklist with the final abstract. | Abstract |

## Introduction

| Item | Checklist item | How this study addresses it | Location |
|---|---|---|---|
| 3a | Healthcare context and rationale, with references to existing models | Diagnostic. The existing model for these data is the published 1,658 pg/mL cutoff. This study reports what that concentration carries as continuous evidence: a likelihood ratio of 2.56 (95% CI 1.40–4.74), against the 28.3 implied by treating the cutoff as a positive test. | Introduction §1 |
| 3b | Target population and intended use in the care pathway, including intended users | As stated under **Intended use**: an add-on test alongside HRCT, DLCO and EIF2AK4 genotyping, informing the pre-treatment safety assessment. Intended users are the clinicians making that assessment. | Introduction §1; Discussion §4.3 |
| 3c | Known health inequalities between sociodemographic groups | Not addressed. No analysis of health inequalities between sociodemographic groups was performed, and the source study states that sex was not considered as a biological variable. Recorded as a limitation. | Discussion §4.2 |
| 4 | Study objectives, including whether development or validation | Development only. No validation dataset exists for this model, and none is claimed. | Introduction §1 |

## Methods

| Item | Checklist item | How this study addresses it | Location |
|---|---|---|---|
| 5a | Sources of data, rationale, representativeness | Secondary analysis of the published supporting data values of Prabhakar et al., *J Clin Invest* 2026;136(2):e199159 (DOI 10.1172/JCI199159). 322 individuals in six cohorts; the 166 group-1 patients (PVOD 26, IPAH 98, HPAH 42) form the analysis sample. PVOD cases are over-represented by design, which is why the model reports a likelihood ratio rather than a probability. | Methods §2.1 |
| 5b | Dates of participant accrual and follow-up | UK National Cohort Study: recruited 19 February 2014 to December 2021. UCSF PH Biobank: February 2023 onwards. UCSF Lung Transplant Registry: March 2013 onwards. COPD cases diagnosed March 1994 to December 2021. Cross-sectional; no follow-up. | Methods §2.1 |
| 6a | Study setting, number and location of centres | Three recruitment sources: the UK National Cohort Study of Idiopathic and Heritable PAH, eight participating centres (ClinicalTrials.gov NCT01907295); the UCSF PH Biobank; the UCSF Lung Transplant Registry. | Methods §2.1 |
| 6b | Eligibility criteria | IPAH, HPAH and PVOD diagnosed by contemporary international guidelines, with resting mPAP and pulmonary arterial wedge pressure from routine clinical visits and the diagnosis date taken as confirmatory right heart catheterisation; COPD by GOLD criteria, excluding COPD patients with clinical evidence of PH; healthy volunteers self-declared free of cardiovascular and respiratory disease. This analysis uses those cohort labels as published and applies no further exclusion. | Methods §2.1 |
| 6c | Treatments received and how handled | All PVOD patients were treated with either vasodilators or calcium channel blockers. Treatment is not a variable in the model and is not adjusted for; GDF15 was measured on samples collected in the course of routine care. | Methods §2.1 |
| 7 | Data pre-processing and quality checking | Implemented in `analysis/data.py`: sheet parsing, normalisation of non-breaking spaces in headers, `"NA"` to missing. The clinical sheets of the source workbook carry no patient identifier, so each row is keyed to its cohort label through the age-adjusted GDF15 value, and the keying is checked end to end by requiring raw and median-centred 6MWD to differ by exactly 294 m on all 112 patients with that measurement. `Patient_ID` is then joined on that same value from the per-patient reconstruction the source-only generated reconstruction (`analysis/reconstruct.py`), one-to-one and exactly, with the diagnosis label required to agree on every row; it is carried through every patient-level output. The join is by key and never by position, since the two files order their patients differently. Sheet `GDF15-1` stands outside this: it is anonymous group-level concentrations, so no output derived from it carries a `Patient_ID`. Cohort counts and both NT-proBNP denominators are asserted at load time. Independent QA of the source workbook is in `QA_note_Prabhakar_et_al_2025_supporting_data.md`. | Methods §2.1; Supplement S2 |
| 8a | Outcome definition, time horizon, how and when assessed | The outcome is the diagnosis of PVOD versus the comparator diagnosis. Diagnosis was made by a multidisciplinary specialist team and dated to confirmatory right heart catheterisation. Cross-sectional; no time horizon. | Methods §2.2 |
| 8b | Qualifications of outcome assessors, if subjective | A multidisciplinary team of radiologists, cardiologists and cardiothoracic surgeons at the contributing centres. | Methods §2.2 |
| 8c | Actions to blind outcome assessment | **PROVISIONAL AUTHOR WORDING; A4/A5 pending.** **Diagnoses of PVOD, IPAH and HPAH were established at the contributing centres according to ESC/ERS diagnostic criteria for pulmonary hypertension, before and independently of GDF15 measurement. Diagnostic labels accompanied the samples but were withheld from the investigator performing the GDF15 assay, who was blinded to group throughout. GDF15 does not form part of the ESC/ERS criteria, so the index test could not have contributed to the reference standard.** Amit must confirm the blinding statement (A5), and Akiko/Amit must confirm diagnostic timing (A4). Neither direction is treated as closed by this QC. | Methods §2.2 |
| 9a | Choice of initial predictors and any pre-selection | One pre-specified predictor, log10(GDF15). No selection, screening, or comparison of candidate predictors was performed, so no selection bias enters the reported model. | Methods §2.3 |
| 9b | Definition and measurement of all predictors | Plasma GDF15 in pg/mL by ELISA (R&D Systems Quantikine DGD150 and DuoSet DY957), plasma diluted 1:10 in PBS, absorbance read at 450 nm, total protein by Bradford assay used to adjust sample volume. This analysis applies a base-10 logarithm and nothing else. | Methods §2.2 |
| 9c | Qualifications of predictor assessors, if subjective | N/A — a quantitative laboratory assay with no subjective interpretation. | — |
| 10 | How study size was arrived at; sample size justification | The sample is the published cohort in full: n = 166 with 26 events and one predictor, that is 26 events per variable. No prospective sample-size calculation was performed, and none was possible for a secondary analysis of a fixed dataset. The consequences of 26 events are reported rather than left implicit: see items 12c and 23a. | Methods §2.1; Discussion §4.2 |
| 11 | How missing data were handled | GDF15 is complete for all 322 individuals, so the model itself uses no incomplete records. The comparator markers are used complete-case, never imputed, with n reported per marker (63 to 166). Because those denominators differ, the comparator table is also recomputed on the 23 patients complete on all five clinical measurements, listed by `Patient_ID` in `results/complete_case_patients.csv`. **That subset cannot confirm the full-cohort ranking and is not reported as though it does.** It is 43.5% PVOD against 15.7% in the full cohort, about three times enriched, because measurement availability is itself associated with diagnosis: cardiac index is present for 96% of PVOD patients against 43% of IPAH and 36% of HPAH, NT-proBNP for 65% against 34% and 36% (`results/clinical_completeness.csv`). PVOD patients are investigated more thoroughly and arrive with fuller data. The complete cases are therefore a different group of patients rather than a smaller version of the study, the rankings among the clinical variables move substantially at that size, and those rows are not an unbiased re-estimate of any AUC. The limitation carried is that **comparator AUCs are conditional on the measurement being available, and measurement availability is associated with diagnosis.** For NT-proBNP the workbook supports two denominators — the distribution sheet used here gives 63 (PVOD 15), the PVOD-only correlation sheet gives 17 for PVOD — and both are stated per patient. | Methods §2.5; Results §3.4 |
| 12a | How the data were used; partitioning; sample size requirements | All 166 patients contribute to one fit; no partitioning into development and test sets, because 26 events cannot support a split that leaves either part informative. Apparent performance is reported beside its optimism-corrected value. | Methods §2.4 |
| 12b | How predictors were handled (functional form, rescaling, transformation) | log10(GDF15) entered linearly; no rescaling, standardisation, spline or categorisation. The linear form is justified rather than assumed: adding a quadratic term gives a likelihood-ratio statistic of 1.84 on 1 degree of freedom, p = 0.175, with the AIC slightly favouring the linear model (54.32 against 54.48). | Methods §2.3; Results §3.2 |
| 12c | Type of model, rationale, model building steps, hyperparameter tuning, internal validation | Binary logistic regression by maximum likelihood on a single predictor (`model.fit_lr_model`). No tuning, no penalty, no variable selection, no stepwise procedure. The fitted log-odds are converted to a likelihood ratio by subtracting the sampling log-odds of the cohort, which removes its deliberate case enrichment. Internal validation is a bootstrap optimism correction (`model.bootstrap_optimism`, 4,000 replicates, seed 1) applied to the C statistic and the calibration slope. No cross-validation and no shrinkage. Leave-one-case-out refitting is reported as a sensitivity analysis (item 12e), not as validation, and does not alter the reported model. | Methods §2.3, §2.4 |
| 12d | Heterogeneity across clusters | N/A. The cohort pools three recruitment sources and their per-cohort counts are published (PVOD 19 from the National Cohort Study, 4 from the PH Biobank, 3 from the Lung Transplant Registry), but the supporting data carry no patient-level *source* identifier, so heterogeneity across sources cannot be estimated. The `Patient_ID` used elsewhere does not help here: it identifies a patient across sheets but carries no information about which centre recruited them. Recorded as a limitation. | Methods §2.1; Discussion §4.2 |
| 12e | Measures and plots used to evaluate performance | C statistic with a stratified percentile bootstrap interval; calibration by quintile of predicted probability with exact (Clopper–Pearson) intervals; the Figure 1 likelihood-ratio curve with a pointwise stratified case-bootstrap percentile band from 3,991 finite ordinary logistic fits out of 4,000 requested resamples, excluding 9 complete-separation resamples, and the range over which each group was actually observed marked on it. The upper edge is effectively unchanged near/below the published cutoff and narrows in the tail. Figure 2 and supplementary parametric intervals retain their existing procedure. Two sensitivity analyses accompany them: the model refitted 26 times with one PVOD case removed each time, with the likelihood ratio read at two fixed concentrations; and the comparator table recomputed on the patients complete on every clinical measurement, reported as a missingness diagnostic rather than as a re-estimate of the ranking (item 11). | Methods §2.4; Results §3.3, §3.4 |
| 12f | Model updating arising from evaluation | N/A — no separate evaluation was performed, so no updating arises. | — |
| 12g | How predictions were calculated for evaluation | N/A — no separate evaluation. The formula allowing prediction in new individuals is item 22. | — |
| 13 | Class imbalance methods | No resampling, class weighting or threshold shifting was applied. The case fraction (15.7%) is a property of how the cohort was assembled, not a nuisance to be corrected: it is removed analytically by subtracting the sampling log-odds, which is what makes the reported likelihood ratio independent of the case mix. | Methods §2.3 |
| 14 | Approaches used to address model fairness | No fairness analysis was performed and no fairness-related approach was applied. With 26 outcome events, subgroup estimates could not be made with useful precision. Recorded rather than left unstated. | Methods §2.4; Discussion §4.2 |
| 15 | Model output, and thresholds if it classifies | The output is a likelihood ratio: continuous, with no threshold and no classification. A post-test probability requires the clinician's own pre-test probability; the model supplies only the update. The fitted probability is also available but applies only at the sampling prevalence of this cohort. | Methods §2.3; Results §3.2 |
| 16 | Differences between development and evaluation data | N/A — no evaluation data. | — |
| 17 | Ethics approval and consent | National Cohort Study: UK REC Ref. 13/EE/0203 (East of England Research Ethics Committee); gene identification studies 08/H0802/32, with written informed consent for whole-blood collection. UCSF PH Biobank: IRB #22-36931. UCSF Lung Transplant Registry: IRB #13-10738. This secondary analysis uses only published, de-identified supporting data. | Declarations |

## Open science

| Item | Checklist item | How this study addresses it | Location |
|---|---|---|---|
| 18a | Funding and role of funders | Jordan’s secondary work: 2026 Responsible AI Summer Research Program at Santa Clara University, award “Transparent AI for Rare Disease Diagnosis”. The programme had no role in analysis. Original-cohort funding is separate (qc/sources.md). | Declarations |
| 18b | Conflicts of interest | *Outstanding — author to supply.* | Declarations |
| 18c | Protocol availability | No prospective protocol; post hoc secondary analysis; no registration; model form selected after data were seen. | Declarations |
| 18d | Registration | The source cohort is registered as ClinicalTrials.gov NCT01907295. This secondary analysis is not itself registered. | Declarations |
| 18e | Data availability | The analysis reads only the published supporting data values of Prabhakar et al. 2026 (DOI 10.1172/JCI199159). No additional participant data were used, and no data beyond the published file are required to reproduce any reported number. | Data and code availability |
| 18f | Code availability | The full analysis code is in the accompanying repository. One command reproduces every table, figure and caption from the published workbook; `README.md` documents dependencies, the entry point, the fixed seeds and every generated output, and `analysis/verify.py` runs numerical checks on the implementation. | Data and code availability |

## Patient and public involvement

| Item | Checklist item | How this study addresses it | Location |
|---|---|---|---|
| 19 | Patient and public involvement in design, conduct, reporting, interpretation or dissemination | None. The proposed presentation has not been tested with clinicians or patients; this remains a limitation. | Declarations |

## Results

| Item | Checklist item | How this study addresses it | Location |
|---|---|---|---|
| 20a | Flow of participants, with and without the outcome | 322 individuals in six published cohorts. 166 enter the main model: 26 PVOD and 140 pooled IPAH (98) and HPAH (42). Other PH (52) and COPD (36) enter only their own comparator models; the 68 healthy controls enter none, because GDF15 separates them from PVOD completely and no finite coefficient exists for that comparison. No patient was excluded by this analysis. | Results §3.1 |
| 20b | Characteristics overall and by data source, including key dates, predictors, sample size, events, missing data | Age, the five clinical measurements and their completeness are reported per group, with the denominator for each; the completeness breakdown by diagnosis is `results/clinical_completeness.csv`, and it shows availability differing sharply by group for the two variables missing most often. Characteristics by recruitment source cannot be given: the supporting data carry no patient-level source identifier (item 12d). | Results §3.1 |
| 20c | Distribution comparison for evaluation data | N/A — no evaluation data. | — |
| 21 | Number of participants and outcome events in each analysis | Reported for every model: main (166, 26 events), and each comparator model separately — IPAH (124, 26), HPAH (68, 26), Other PH (78, 26), COPD (62, 26) — together with the per-quintile counts underlying the calibration figure. | Results §3.1, §3.3 |
| 22 | Full model specification allowing prediction in new individuals | logit P(PVOD) = −39.247 + 11.959 × log10(GDF15 in pg/mL), fitted on a sample with 15.7% PVOD, sampling log-odds −1.6835. The likelihood ratio is log LR(x) = [−39.247 + 11.959 × log10(x)] + 1.6835, and post-test odds = pre-test odds × LR. Everything needed to compute a likelihood ratio for a new patient is in this row. | Results §3.2 |
| 23a | Model performance with confidence intervals, including key subgroups | C statistic 0.956 (95% CI 0.912–0.989; stratified percentile bootstrap, 4,000 replicates, seed 1). An asymptotic C-statistic interval computed here, 0.9160–0.9960, is reported beside it and never combined with it; the source study's published Figure 1E interval, 0.9164–0.9957, is quoted separately as a published value and not as ours. Quoted to four decimals because at three they coincide: the two are close but not the same interval, and the difference is what distinguishes a value computed from these data from one transcribed from the published figure. Calibration by quintile with exact intervals. Internal validation: C statistic 0.9560 apparent against 0.9562 optimism-corrected; calibration slope 1.0000 apparent against 0.9676 corrected. Performance is not reported in any subgroup: with 26 events, subgroup estimates would not be informative. | Results §3.3 |
| 23b | Heterogeneity across clusters | N/A — see item 12d. | — |
| 24 | Results of model updating | N/A — no model updating was performed. | — |

## Discussion

| Item | Checklist item | How this study addresses it | Location |
|---|---|---|---|
| 25 | Overall interpretation, including fairness, in the context of objectives and previous studies | The model reproduces the source study's discrimination and reinterprets its cutoff: at 1,658 pg/mL the measurement carries a likelihood ratio of 2.56, taking a 20% pre-test probability to 39.0% (25.9–54.2%), where treating the cutoff as a positive test would give 87.6%. Aggregate discrimination is high and stable, and that is not the same claim as support for an individual patient's estimate. No fairness analysis was performed (item 14). | Discussion §4.1 |
| 26 | Limitations and their effect on bias, statistical uncertainty and generalisability | 26 events, one cohort, no external validation set. No PVOD patient was measured below 860.5 pg/mL and no IPAH or HPAH patient above 1,809.5 pg/mL, while 20 of 26 PVOD patients lie above that value, so the curve's upper region rests on the model's form rather than on observed comparators. Leave-one-case-out refitting quantifies that directly: removing one PVOD case moves the likelihood ratio by at most 8.1% at 1,658 pg/mL but by 177.3% at 5,986 pg/mL. Recruitment-source heterogeneity is not estimable. No fairness or subgroup analysis. The two source Figure 1E values are independently verified from the supplied full figure, as confirmed in the first-author QC brief; this item is closed. | Discussion §4.2 |
| 27a | How poor quality or unavailable input data should be handled in implementation | The model takes one input and returns a likelihood ratio; with no GDF15 measurement there is no output and pre-test suspicion is unchanged. Above 1,809.5 pg/mL and below 860.5 pg/mL the returned value reflects the fitted model form rather than observation of both groups, and both regions are marked on Figure 1 and should be treated as such. | Discussion §4.3 |
| 27b | Whether users interact with the input data or model, and the expertise required | The intended users are clinicians making the pre-treatment safety assessment. They supply their own pre-test probability and combine the returned likelihood ratio with HRCT, DLCO, EIF2AK4 genotyping and the rest of the clinical picture. It is an add-on test, not a stand-alone one, and it requires no interaction with the underlying data. | Discussion §4.3 |
| 27c | Next steps for future research, applicability and generalisability | External validation in an independent cohort with comparator patients measured across the full GDF15 range is the necessary next step; the upper region of the curve cannot be settled by these data at any sample size, because the comparator observations are absent rather than sparse. Joint modelling of GDF15 with the clinical variables is blocked twice over. Raw GDF15 cannot enter a patient-level model at all: the source sheet `GDF15-1` contains anonymous group-level concentrations without a patient identifier, so it cannot be joined to a clinical row and no combined model can be built. Age-adjusted GDF15 can be joined, but only 23 patients have every clinical measurement and they are strongly case-enriched, with availability associated with diagnosis, so a complete-case multivariable fit would estimate a different population from the one of interest. | Discussion §4.3 |

---

## Outstanding items

A1 author list/order/corresponding author; A2 assay details; A4 accrual/guideline timing; A5 final blinding confirmation; A6 COI (Giorgio and counsel); A7 Akiko's framing approval; separate SCU determination; A9 repository/Zenodo details. No COI wording is finalized.

Funding (18a), protocol (18c), registration and patient/public involvement (19) are confirmed. A3 has verified source wording, with consent scoped to what the supplement states. A8 is verified (qc/sources.md). Abstract checklist and final manuscript page mapping still need submission assembly.

## Closed source Figure 1E item

Independent verification of the supplied full figure was confirmed in the first-author QC brief: PVOD vs IPAH+HPAH+Other PH **96.05% (92.45–99.66)**; PVOD vs IPAH+HPAH **95.60% (91.64–99.57)**. These are published values; our computed intervals are separate.

## Positions recorded rather than left open

- **Blinding (8c)**: A5 remains pending Amit’s final confirmation. Prior project wording is not a substitute for that confirmation.
- **Fairness (14) and health inequalities (3c).** No fairness analysis was performed;
  with 26 events, subgroup estimates could not be made with useful precision.
- **Subgroup performance (23a).** None reported, for the same reason.
- **Clustering by recruitment source (12d, 23b).** Not estimable: the supporting data
  carry no patient-level identifier for the three recruitment sources. This is a separate
  question from patient identity, which is settled — `Patient_ID` names a patient across
  sheets but says nothing about which centre recruited them.
- **Raw GDF15 at patient level.** Not available and not constructible. Sheet `GDF15-1` is
  six anonymous cohort columns of unequal length with no identifier, no age and no shared
  row order, so neither the positional join nor the age-adjusted GDF15 key reaches it. This
  is why no combined GDF15-plus-clinical model was built, why the main analysis uses GDF15
  with the diagnosis label only, and why the complete-case table cannot include raw values.
  Constructing a key would mean guessing which patient goes with which value.
