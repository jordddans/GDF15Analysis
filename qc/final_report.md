# First-author QC completion report

7 September 2026. **The numerical centerpiece passes; the draft is not submission-ready.** The agreed missing-script gate prevents new Figure 1 integration and the `analysis-v1.0` release. The original manuscript, published workbooks and submitted BayLearn PDFs are unchanged. Mechanical corrections are in the separate [QC DOCX](../GDF15_evidence_curve_QC.docx), with a [rendered PDF](manuscript_QC.pdf) and [32-edit log](manuscript_changes.json). Scientific framing concerns are recommendations, not silent rewrites.

## A. Independent Table 1 — PASS

Direct original-workbook extraction, independent of Giorgio's script:

| Diagnosis | Original cell | GDF15, pg/mL | Continuous LR | Post-test probability, 10% pre-test | Cutoff LR | Cutoff post-test probability |
|---|---|---:|---:|---:|---:|---:|
| IPAH | Fig_1a!C36 | 1657.8125 | 2.5554087466205235 | 22.1143951084% | 0.19796380090497737 | 2.1522567950% |
| PVOD | Fig_1a!B27 | 1658.81 | 2.5634044012059505 | 22.1682500435% | 28.26923076923076 | 75.8513931889% |

The exact decimal gap is **0.9975 pg/mL** (about 0.0602% of the concentration). Exhaustive enumeration of **26 × 140 = 3,640** cross-group pairs establishes a **unique global minimum**. Both values exist in the correct columns, are distinct at available precision and occur nowhere else in the six relevant diagnosis columns. The published 1658 cutoff lies strictly between them and classifies both correctly. No clinical Patient_ID is attached to either raw measurement.

The final coefficients and prevalence offset independently reproduce both LRs and probabilities. [Machine-readable verification](../results/centerpiece_verification.json), [all pairwise differences](../results/cross_group_gaps.csv), and [programmatically generated Table 1](../results/table1.csv) preserve full values and source coordinates.

## B. Figure 1B/C point estimates — PASS

At 1658 pg/mL: TP 21, FN **5/26**, FP **4/140**, TN 136; **157/166** correct. Sensitivity **80.7692%**, specificity **97.1429%**, LR− **0.1979638** (0.20), LR+ **28.2692** (28.3). Continuous LR near the cutoff is approximately 2.56; the adjacent measurements both give approximately 22%, versus approximately 2% and 76% under the cutoff rule.

Across the entire cohort, **34/166** differ by strictly more than 10 percentage points and **11/166** by strictly more than 20. [Observation-level comparison](../results/cutoff_probability_comparison.csv) uses original raw-data cell references, not unsupported clinical IDs. Numerical verification does not certify the supplied image's bootstrap band or generation code.

## C. 0.24 audit

**394 occurrences in 38 existing files**, plus OCR of 24 raster figures/previews with zero additional hits. Every occurrence and classification is listed in [the ledger](audit_024.csv); scope and per-file counts are in [the audit summary](audit_024_summary.json). Environments/caches and recursive matches in new audit explanations are excluded.

**Zero erroneous fixed-cutoff LR− occurrences; zero numerical replacements.** The two displayed LR− 0.24 values in archived classifier PDFs correctly describe nested leave-one-out/reselected-cutoff performance: six false negatives yield 0.2375566. They are not the fixed-cutoff five-false-negative LR−. These PDFs remain unchanged. README, archive README and BayLearn reconciliation now explain this distinction. Other matches are measurements, model output, code or vector geometry. Private full contexts remain local and ignored.

## D. Ethics source — VERIFIED, with scope limits

Original JCI supplement `archive/gdfdiscrim/jci-136-199159-s067.pdf`, printed pp. 5–6, Study Approval; consent in the preceding participant section on p. 5. Independently corroborated against the [published supplement](https://dm5migu4zj3pb.cloudfront.net/manuscripts/199000/199159/199159-JCI-RL-RV-3_sd_949922.pdf).

Verified manuscript wording: “The original study reported approval of the UK National Cohort Study protocol (13/EE/0203) by the East of England Research Ethics Committee and of the UCSF PH Biobank and Lung Transplant Registry protocols (22-36931 and 13-10738, respectively) by the UCSF Institutional Review Board. The NCS gene-identification protocol (08/H0802/32) was approved by the UK Research Ethics Committee; the supplement reports written informed consent for whole-blood collection for those gene-identification studies.”

The source does not establish universal consent wording for every plasma sample. The separate SCU determination remains pending. [Exact local-source passages and provenance](sources.md).

## E. 2026 PVOD/PCH score — VERIFIED

Andruska AM, Cruz-Utrilla A, Nossent EJ, et al. *Identification of pulmonary arterial hypertension patients with venous or capillary involvement.* AJRCCM 2026;212(9):2199–2210. DOI **10.1093/ajrccm/aamag289**; PMID **42275164**. The [published Discussion](https://academic.oup.com/ajrccm/article/212/9/2199/8706424), paragraph beginning “The approach for developing…”, explicitly identifies “information loss from dichotomization at thresholds”. A8 is closed with this source. It does not substantiate the separate assertion that external validation data are unobtainable.

## F. Funding and BayLearn

The QC manuscript and project documentation identify the **2026 Responsible AI Summer Research Program at Santa Clara University**, award **“Transparent AI for Rare Disease Diagnosis”**, and that the programme had no role in the analysis. Original-cohort funding remains separate. An acknowledgement was added to the editable archived BayLearn TeX source, explicitly documented as a post-submission update. The submitted PDFs are unchanged. No summer final report existed; [a standalone course-correction note](summer_report_note.md) describes the scientifically motivated refinement from an explainable classifier to a continuous evidence model.

## G. Reconstruction and privacy — PASS

The private workbook is ignored, untracked and absent from all local git history. Historical private clinical CSVs, notebooks/chat with potentially embedded records, and generated clinical output are also excluded. No remote exists.

`analysis/reconstruct.py` rebuilds Patient_ID 1–166, grouping, clinical values and provenance from the published workbook alone. Local comparison against the known-good workbook passed **all 22 fields** across Patients_166 and Source_values, including missingness. All values match exactly except a maximum NT-proBNP floating-point difference of 8.53e−14 (absolute tolerance 1e−12). See [comparison](reconstruction_verification.json) and [join/provenance documentation](reconstruction.md).

The original age-adjusted sheet lacks diagnostic labels: grouping uses the established, independently cross-checked block convention, not invented labels claimed to exist in that sheet. Joins reject ambiguous or missing keys. **Raw GDF15 remains unjoined to clinical Patient_ID.** The published and cleaned source sheets match at all **4,609 numeric cells**, while text provenance differences are explicitly documented.

## H. New Figure 1 — BLOCKED

Giorgio's script was not present in the project or DOCX; only its PNG was supplied. Per the user's explicit decision, the historical reproducible Figure 1 is retained and the new image remains in the QC manuscript for review. No substitute implementation has been invented.

There is a material methods conflict: the draft describes a stratified case-bootstrap LR band, whereas the active repository uses a parametric-bootstrap LR band. Existing AUC/stability procedures use stratified resampling, and optimism uses ordinary resampling. Existing 4,000-replicate requests and deterministic seeds are preserved. The new image also does not visibly mark the upper comparator support limit despite the text claiming it does; the exact accessibility metric is unverified. [Integration gate and reconciliation requirements](figure1_integration.md).

## I. Hostile-review findings

[The complete audit](hostile_review.md) and [machine-readable table](hostile_review.csv) contain **64 entries** with original section/body-block locations, exact claims, category, severity, independent verification status and recommended action. There are **34 OPEN entries: 10 blocking, 20 important and 4 minor**, plus two partly resolved important entries. These counts include existing author placeholders, not ten newly discovered numerical failures.

Blocking issues include the Figure 1/bootstrap conflict, an unsupported decision-curve-analysis claim, unresolved author/assay/timing/blinding/COI/framing items, SCU determination and repository availability. Important scientific concerns include population transport claims, treating fitted probabilities as individual truth, overstating what negligible optimism proves, the claimed conventional events-per-variable limit, and unsupported assertions about where information loss occurs or why validation is unobtainable.

Mechanical corrections include the AUC CI upper endpoint **0.989**, corrected calibration slope **0.968**, stale calibration-bin values, and Table 2 LR=10 bootstrap interval/width rounding (**1920–2360**, width **430** when each quantity is rounded from full precision). Main analysis numbers were not changed. Table 3, LOO landmarks, 43.5%/2.8-fold enrichment and comparator support reproduce. Minor final assembly work includes a blank cover-separator page 2 in the QC PDF; remove it when preparing the submission file.

## J. Reference audit

All **eight bibliographic identities**, including DOI/PMID, are verified in [the reference audit](reference_audit.csv), with sentence-level support limits and recommended actions. Prabhakar, Altman/Royston, ESC/ERS, TRIPOD+AI, STARD and the diagnostic-score source are accounted for. Full-text access was limited for the Aguado review and Steyerberg methods paper: their identities and accessible abstract-level support are verified, but a blanket full-text claim audit is not asserted. Bootstrap/calibration calculations were checked independently in code. References do not validate stronger manuscript rhetoric merely by being relevant.

## K. Git state

Baseline commit: **c4a5a20dd2024f6a35bb0b898829cdd18e347936**, created after privacy exclusions and before analysis edits. The completion commit contains this report; obtain its immutable identifier with `git log -1 --format=%H` (also supplied in the final handoff). It is a **QC completion commit**, not a final scientific release.

**No remote, no push, no Zenodo release, no `analysis-v1.0` tag.** The tag is deliberately withheld under the agreed Figure 1 script gate. Once that gate is resolved, the next repository publication action requires creating/approving a **PRIVATE remote**; no public remote has been created.

## L. Verification — PASS for available workflow

A fresh Python 3.14.6 environment installed `requirements-lock.txt`; an isolated code tree used only the original published workbook. The complete available workflow passed **98 checks, zero failures**. These include the independent centerpiece and manuscript-number checks and deterministic reruns. **32 analysis/figure/report artifacts** and **four generated reconstruction artifacts** were byte-identical across the final independent runs. The separate QC DOCX generation was also byte-identical on repeat.

Baseline regression found **zero changes across 685 numerical/Boolean fields** (675 numbers and 10 Booleans). All existing scientific CSV values, model coefficients, main AUC, LR landmarks, calibration, internal validation, missingness/complete-case and LOO outputs are preserved. All nine existing figure/caption artifacts, including Figures 2 and 3, are byte-identical. Original workbooks, original manuscript and submitted PDFs match baseline hashes. [Clean run](clean_run.json), [regression](baseline_regression.json), [workflow checks](../results/verification.json), [document checks](document_verification.json).

This is not a PASS for new Figure 1 reproducibility, which cannot be tested without its script. README documents the simple source-only workflow and this limitation.

## M. Remaining external-author blockers

A1 author list/order/corresponding author; A2 assay details; A4 accrual and guideline timing; A5 blinding confirmation; A6 Giorgio/counsel COI wording; A7 Akiko's Discussion approval remain untouched as author decisions. A3 published source is verified, but the independent SCU determination remains pending. A8 is verified. A9 repository URL/Zenodo remains pending. Giorgio's Figure 1 script and bootstrap reconciliation are required before release. No author confirmation has been inferred from silence.

## N. Files changed

[Complete relative-path inventory](files_changed.txt). Principal changes: source-only reconstruction and independent verification code; workflow/README and environment lock; generated Table 1 and machine-readable audit outputs; QC manuscript copy and PDF; TRIPOD/STARD mappings; reference/ethics and hostile-review audits; funding, BayLearn reconciliation and editable acknowledgement; summer-report note; privacy exclusions. The original source workbooks, manuscript and submitted PDFs were not edited.
