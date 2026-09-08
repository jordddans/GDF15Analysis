# GDF15 → PVOD Likelihood Ratio Engine — Build Specification

**For:** Claude Code, implementing in Jupyter Notebook
**Source study:** Prabhakar et al., *J Clin Invest* 2026;136(2):e199159 (DOI 10.1172/JCI199159)
**Status of this document:** authoritative build spec. Numbers marked **[VALIDATE]** are verified against the source data and must be reproduced exactly by the implementation, or the pipeline is wrong.

---

## 0. Purpose and framing

### What we are building

A tool that accepts a raw plasma GDF15 concentration (pg/mL) and returns **the likelihood ratio that value carries for pulmonary veno-occlusive disease**, with honest uncertainty bounds, plus a Bayesian chaining engine that takes a clinician's pre-test probability (assembled from evidence external to this study) and shows **how far GDF15 moves it**.

### What we are replacing

The source paper reduces GDF15 to a single Youden-derived hard cutoff of **1,658 pg/mL** (81% sensitivity, 98% specificity vs IPAH+HPAH). A hard cutoff:

- Discards magnitude. A patient at 1,700 and a patient at 5,000 are treated identically, when the second is enormously more informative.
- Reports point estimates with no uncertainty, from **26 PVOD cases**.
- Produces a binary "positive/negative" that clinicians correctly distrust, because it does not compose with the rest of the diagnostic picture.

### The thesis this tool exists to support

> From only 26 PVOD patients, here is precisely how much diagnostic persuasion a GDF15 measurement can honestly supply — no more, no less — and here is how it composes with the evidence a clinician already has.

The tool must be *more* conservative than the source paper, not less. Every design decision below resolves in favour of understatement.

### Critical framing decision: likelihood ratios, not probabilities

**Do not fit a model that outputs P(PVOD | GDF15) directly.**

The study sample is 26 PVOD out of 322 (8%). That figure is an artifact of how the UK NCS, UCSF PH Biobank, and UCSF Lung Transplant Registry were assembled — it is not a prevalence in any real population. A logistic regression fitted on this data bakes that artificial base rate into its intercept, and every probability it emits is wrong in a way that is invisible to the end user.

The **likelihood ratio** is a property of the two conditional distributions `f(x | PVOD)` and `f(x | comparator)`. It is invariant to sampling prevalence, so it transfers to any clinical setting. The clinician supplies the pre-test probability; the model supplies only the update.

This is also exactly the quantity the project is about. An LR *is* the unit of persuasion.

---

## 1. Data assets and reconstruction

All data derives from `jci136199159s068.xlsx` (JCI source data file). **This file contains individual-level values, not summary statistics.** That is the whole reason this project is feasible.

### 1.1 Primary extraction — sheet `Fig_1a`

Row 3 is the header. Six columns, ragged length, one raw plasma GDF15 value (pg/mL) per cell.

> **Note:** column headers contain non-breaking spaces (`\xa0`). Normalise with `str.replace('\xa0', ' ')` before matching, or `Other PH` will silently fail to parse.

**[VALIDATE]** After extraction:

| Cohort | n | Median (pg/mL) | IQR |
|---|---|---|---|
| Ctrl | 68 | 255.7 | 157.3 – 339.4 |
| PVOD | 26 | 3170.7 | 2151.2 – 3634.1 |
| IPAH | 98 | 777.0 | 484.9 – 1115.0 |
| HPAH | 42 | 606.1 | 398.6 – 716.1 |
| Other PH | 52 | 548.1 | 302.9 – 805.3 |
| COPD | 36 | 620.4 | 474.9 – 1087.1 |

These reproduce Figure 1A of the paper exactly. Total n = 322.

### 1.2 Sex labels — sheet `Fig_S2a`

Same 322 values, split into M/F columns per cohort. Join back to `Fig_1a` values to attach a sex label to every observation. Counts: Ctrl 27/41, PVOD 10/16, IPAH 36/62, HPAH 9/33, Other PH 21/31, COPD 21/15.

Used only for the subgroup-stability check (§4.6). The paper found no significant sex difference; we verify rather than assume.

### 1.3 Age linkage — non-obvious reconstruction, implement carefully

`cohort_166_template.csv` holds 166 group-1 patients (PVOD 26, IPAH 98, HPAH 42) with `Age_years`, `AgeAdj_GDF15`, and assorted clinical variables — **but not raw GDF15**. The raw values live in `Fig_1a` and **are not row-aligned** with the template (verified: Spearman ≈ 0.07 under a naive row join).

The linkage is recoverable. The paper's `AgeAdj_GDF15` is the residual from a GAM, `gam(log_GDF ~ s(log_age))`, per their R script. Therefore:

```
log(GDF15_i) = f(L_i) + residual_i,    where L_i = log(age_i + 1)
```

We know `age_i` and `residual_i` for all 166, and we know the *multiset* of raw GDF15 values per group. Fit `f` by minimising the squared distance between the sorted predicted log values and the sorted observed log values. A quadratic in `L` suffices:

```python
f(L) = 7.77930799 - 1.97578884*L + 0.40243762*L**2
```

Then, **within each group separately**, rank patients by `f(L_i) + residual_i` and assign the correspondingly ranked observed raw GDF15 value.

**[VALIDATE]** Reconstruction quality: RMSE ≈ 0.010 in log units overall (PVOD 0.0087, HPAH 0.0070, IPAH 0.0176). Recovered group medians must equal 3170.7 / 777.0 / 606.1 exactly. If RMSE exceeds 0.05, stop and report — the assignment is unreliable.

**Caveat to carry forward:** IPAH max deviation is 0.113 log units, so one or two adjacent IPAH patients may be transposed. Immaterial for distributional work; do **not** use this reconstruction for any per-patient claim.

### 1.4 Final analysis frames

- **`df_all`** — 322 rows: `gdf15_raw`, `cohort`, `sex`. Complete, no imputation.
- **`df_g1`** — 166 rows: adds `age_years`, `ageadj_gdf15`, and clinical variables (`6MWD_m`, `mPAP_mmHg`, `PVR_dyn_s_cm5`, `NTproBNP_ng_L`, `CI_distFromMedian`). Clinical variables are heavily missing (NT-proBNP 101/166 missing, CI 84/166, 6MWD 54/166) — **complete-case only, never impute**, and report the n for every analysis that uses them.

---

## 2. Comparator groups

Build **two parallel models**, reported side by side throughout. Never collapse them into one headline number.

| Model | Comparator | n | Clinical question |
|---|---|---|---|
| **A** | IPAH + HPAH | 140 | "This patient has precapillary PAH and I have labelled them idiopathic or heritable. Is it actually PVOD?" |
| **B** | IPAH + HPAH + Other PH | 192 | "This patient has pulmonary hypertension of some kind. Is it PVOD?" |

Model A is the decision that matters — it is the one where getting it wrong means giving a vasodilator to a PVOD patient and precipitating pulmonary oedema. Model A also matches the paper's 1,658 cutoff derivation. Treat A as primary, B as the generalisability check.

### Exclusions — mandatory, and state the reasoning in the notebook

**Healthy controls (n=68) must be excluded from all diagnostic modelling.** Median age 33 vs 64 in PVOD, and no cardiopulmonary disease. The 100%/100% sensitivity/specificity in Figure 1B is substantially an age-and-health artifact. No clinician is ever distinguishing PVOD from a healthy 33-year-old. Including controls would inflate every performance metric and is the single largest source of optimism in the source paper.

**COPD (n=36) excluded from primary models.** Different clinical pathway; not a differential for PVOD. Retain as an out-of-scope specificity probe only.

---

## 3. The models

Fit all three. Agreement between them is evidence; disagreement is a finding to report.

### 3.1 Primary — two-parameter lognormal per group

Work throughout on `y = log(gdf15_raw)`.

Fit `N(μ₁, σ₁)` to PVOD and `N(μ₀, σ₀)` to comparator, **allowing unequal variances**. Then

```
LR(x) = φ((log x − μ₁)/σ₁)/σ₁  ÷  φ((log x − μ₀)/σ₀)/σ₀
```

**[VALIDATE]** Fitted parameters (log scale):

| Group | n | μ | σ |
|---|---|---|---|
| PVOD | 26 | 7.887 | 0.509 |
| Comparator A (IPAH+HPAH) | 140 | 6.504 | 0.607 |
| Comparator B (+Other PH) | 192 | 6.406 | 0.687 |

**[VALIDATE]** Resulting LR curve:

| GDF15 (pg/mL) | LR, Model A | LR, Model B |
|---|---|---|
| 400 | 0.004 | 0.004 |
| 800 | 0.08 | 0.09 |
| 1200 | 0.56 | 0.65 |
| **1658** (paper cutoff) | **2.37** | **2.56** |
| 2000 | 5.21 | 5.22 |
| 3000 | 24.79 | 19.77 |
| 4000 | 66.90 | 42.71 |
| 6000 | 230.75 | 99.03 |

Note the paper's cutoff of 1,658 corresponds to an LR of only ~2.4. That is weak evidence. Surface this prominently — it is one of the strongest arguments for replacing the cutoff with a continuous function.

**Assumption check:** Shapiro–Wilk on log scale — PVOD p = 0.023, HPAH p = 0.223, IPAH p < 0.001, Other PH p = 0.001. Lognormality is **imperfect**, particularly in the comparator. Report these p-values in the notebook. This is precisely why the empirical model below is not optional.

Levene's test for equal variance: A p = 0.180, B p = 0.098. Unequal variance is not strongly supported but is permitted; fit both and report whether the LR curve differs materially.

### 3.2 Secondary — empirical banded LR (model-free)

Bands: `<500, 500–1000, 1000–1658, 1658–2500, 2500–4000, ≥4000`.

`LR_band = (a + 0.5)/n₁ ÷ (b + 0.5)/n₀`, Haldane–Anscombe correction. Bootstrap CI, 4000 stratified resamples.

**[VALIDATE]** Model A:

| Band (pg/mL) | PVOD | IPAH+HPAH | LR | 95% CI |
|---|---|---|---|---|
| <500 | 0 | 42 | 0.06 | 0.05 – 0.10 |
| 500–1000 | 2 | 62 | 0.22 | 0.04 – 0.50 |
| 1000–1658 | 3 | 32 | 0.58 | 0.09 – 1.3 |
| 1658–2500 | 3 | 4 | 4.19 | 0.8 – 20 |
| 2500–4000 | 14 | 0 | 156 | 102 – 210 |
| ≥4000 | 4 | 0 | 48 | 16 – 92 |

**Two things the notebook must say explicitly about this table:**

1. **The non-monotonicity (156 then 48) is noise, not biology.** It reflects 18 patients split across two bins. This is the empirical argument for preferring the smooth parametric LR as the primary estimate.
2. **The top two rows are not real numbers.** Zero comparator patients in those bands means the point estimate is driven entirely by the continuity correction, and the bootstrap CI is fake — it cannot express uncertainty about a cell that never varies. See §4.1.

### 3.3 Tertiary — logistic regression on log GDF15

`isPVOD ~ log(gdf15)` and `isPVOD ~ log(gdf15) + log(gdf15)²`. The LR is recoverable from the fitted odds after removing the intercept's prevalence contribution. The linear model is mathematically equivalent to the equal-variance lognormal case; the quadratic term tests curvature. Use as a cross-check on §3.1 and to formally test whether unequal variance is warranted (likelihood ratio test on the quadratic term).

---

## 4. The honesty layer

This is the part that distinguishes the tool. Implement all seven.

### 4.1 Rule-of-three ceiling on the upper LR — the single most important measure

Above ~2,000 pg/mL there are **zero** comparator patients. The point estimates (156, 231) are model extrapolation into a region with no negative-class data, and they are not defensible.

What the data supports is a **lower confidence bound** on LR, computed from the Jeffreys lower bound on sensitivity and the Jeffreys upper bound on the comparator's false-positive rate:

**[VALIDATE]** Model A:

| Threshold | Sens [95% CI] | FP / n | FPR upper bound | **Defensible LR ≥** |
|---|---|---|---|---|
| 1658 | 0.81 [0.63–0.92] | 4/140 | 0.0665 | **9.5** |
| 2000 | 0.77 [0.58–0.90] | 0/140 | 0.0178 | **32.9** |
| 2500 | 0.69 [0.50–0.84] | 0/140 | 0.0178 | **28.3** |
| 3000 | 0.62 [0.42–0.78] | 0/140 | 0.0178 | **23.9** |

Model B (n=192) gives correspondingly tighter bounds: LR ≥ 12.9 / 45.1 / 38.7 / 32.7 at the same thresholds.

**Implementation rule:** the tool's reported LR is `min(parametric_LR, ceiling)` where the ceiling is set by the rule-of-three bound at the relevant threshold. Recommend a **hard cap of LR = 50 for Model A**, with the displayed message:

> "Point estimate exceeds what 140 comparator patients can resolve. The data firmly supports LR ≥ 33; the true value may be higher but this sample cannot demonstrate it."

This is the honest answer to "how much persuasion do 26 patients buy," and it is the headline result of the whole project.

### 4.2 Extrapolation flagging

Observed PVOD support is **860.5 – 5,986.4 pg/mL**. Comparator A support is **91.0 – 1,809.5 pg/mL**.

Three zones, displayed to the user:

- **Interpolation** (860 – 1,810): both classes observed. Estimates are data-driven.
- **Single-class extrapolation** (< 860 or 1,810 – 5,986): only one class observed. Direction is trustworthy, magnitude is model assumption.
- **Full extrapolation** (> 5,986): no data of any kind. Refuse to return a point estimate; return the capped bound and a warning.

### 4.3 Optimism correction

In-sample AUC is optimistic. Report:

- Apparent AUC — **[VALIDATE]** Model A: **0.956, bootstrap 95% CI 0.913 – 0.991**
- Leave-one-out cross-validated LR estimates (refit the lognormal excluding each patient, evaluate LR at that patient's value)
- 0.632+ bootstrap optimism-corrected AUC

Expect a small correction — AUC is nearly unbiased — but the *LR estimates* at the tails will degrade noticeably under LOO. Show that.

### 4.4 Full CI propagation

Every LR reported must carry a stratified bootstrap CI (resample within PVOD and within comparator, preserving group sizes, ≥4000 reps, percentile method). The CI must propagate all the way through to the post-test probability in §5 — the user sees an interval, never a bare number.

### 4.5 The n=26 resolution limit

With 26 cases you cannot resolve an LR curve finely. Demonstrate this concretely: report the number of PVOD patients contributing to each region of the curve, and show a "one patient added / one patient removed" sensitivity analysis — refit with each PVOD patient dropped and plot the resulting envelope of LR curves. The width of that envelope *is* the answer to how much a 26-patient sample can tell you.

### 4.6 Subgroup stability

Refit separately by sex (PVOD 10 M / 16 F) and report whether the LR curve shifts. With 10 male PVOD patients this analysis is underpowered by construction — state that as the finding rather than reporting a spuriously precise comparison.

### 4.7 Age handling

Per the decision on scope: **raw GDF15 is the primary input; age is optional.**

**[VALIDATE] The key age result — report this prominently, it is good news:**

| Marker | AUC for PVOD vs IPAH+HPAH |
|---|---|
| Raw GDF15 | 0.956 |
| Age-adjusted GDF15 | 0.967 |
| **Age alone** | **0.522** |

Age alone is coin-flip. The discrimination is **not** an age artifact, despite PVOD patients (median 64.9) being older than HPAH (47.7). Note that PVOD and IPAH have near-identical median ages (64.9 vs 69.0), which is why age carries no signal for the distinction that matters.

Within group-1, log GDF15 correlates with age at Spearman rho = +0.334 (p < 0.001), so an age-adjusted mode gains a little (0.956 → 0.967). Offer it as an optional second input using the recovered GAM function from §1.3, clearly labelled as reconstructed rather than original.

---

## 5. The Bayesian chaining engine

### Design

```
pre-test odds  ×  LR_GDF15  =  post-test odds
```

The pre-test probability comes from **outside this study** — the clinician's assessment based on published PVOD evidence: HRCT findings (septal lines, centrilobular ground-glass, mediastinal lymphadenopathy), markedly reduced DLCO, biallelic *EIF2AK4* mutation, occupational exposure, pulmonary oedema on vasodilator challenge, and the published base rate of PVOD among cases labelled IPAH.

### Required interface

Three ways to supply pre-test probability:

1. **Direct entry** — user types a percentage.
2. **Evidence chain** — user enters a base rate plus a series of named LRs from external literature (a small editable YAML/dict library, user-populated), which multiply into a pre-test odds. **The tool ships with this library empty**; it is populated from sources the user supplies, and every entry must carry a citation field.
3. **Preset scenarios** — 2–3 worked examples with explicitly labelled assumed inputs.

### Output — the movement is the point

The headline is not the post-test probability. It is the **delta**:

```
Pre-test probability of PVOD:        12%   (from CT + DLCO, user-supplied)
GDF15 measured:                      3,140 pg/mL
Likelihood ratio (Model A):          24.8  [capped at 50; defensible floor 28.3]
Post-test probability:               77%   (95% CI 58% – 89%)
──────────────────────────────────────────
GDF15 moved this diagnosis by:       +65 percentage points
```

Include an interactive Fagan nomogram. Include a plot of post-test probability against GDF15 across the full input range, at several fixed pre-test probabilities — this single figure is the clearest statement of the project's thesis.

### Conditional independence — state the assumption, then test it

Chaining LRs assumes conditional independence given disease status. This is an assumption, and it must be declared every time the chain is used.

**We can partially test it, and the results support chaining.** **[VALIDATE]** Within group-1, Spearman correlation of log GDF15 with:

| Marker | n | rho | p |
|---|---|---|---|
| NT-proBNP | 65 | +0.080 | 0.528 |
| mPAP | 157 | −0.117 | 0.143 |
| PVR | 142 | −0.140 | 0.095 |
| 6MWD | 112 | **−0.421** | **<0.001** |
| Age | 166 | +0.334 | <0.001 |

And, discriminating PVOD from IPAH+HPAH, every conventional marker is near-useless:

| Marker | n (PVOD) | AUC |
|---|---|---|
| GDF15 | 166 (26) | **0.956** |
| 6MWD | 112 (22) | 0.673 |
| mPAP | 157 (24) | 0.617 |
| PVR | 142 (23) | 0.602 |
| CI | 82 (25) | 0.567 |
| NT-proBNP | 65 (17) | 0.558 |

**This is the core quantitative argument of the project.** GDF15 is essentially uncorrelated with the haemodynamic and natriuretic-peptide data a clinician already has, and it discriminates PVOD from IPAH/HPAH where none of those markers can. It is supplying genuinely new information, not restating the right-heart catheterisation. Chaining GDF15 onto NT-proBNP, mPAP, PVR, or CI is therefore reasonably defensible.

**Exception:** 6MWD correlates at rho = −0.42. Chaining GDF15 with 6MWD will overstate the combined evidence. Emit an explicit warning if the user's evidence chain includes a functional-capacity term.

No data is available here on CT findings, DLCO, or *EIF2AK4* — the most important pre-test elements. Independence from those is **assumed and untestable with this dataset**. Say so in the output, every time.

---

## 6. Notebook structure

```
01_extract.ipynb          §1  — parse xlsx, reconstruct age linkage, build df_all / df_g1,
                                run all [VALIDATE] assertions as hard `assert` statements
02_distributions.ipynb    §3  — lognormality checks, fit models A and B, all three estimators
03_lr_curves.ipynb        §3  — LR functions, bootstrap CIs, banded table, curve comparison
04_honesty.ipynb          §4  — ceilings, extrapolation zones, LOO, leave-one-out envelope,
                                subgroup stability, age analysis
05_bayes_engine.ipynb     §5  — chaining engine, Fagan nomogram, scenario plots
06_report.ipynb           —     assembled figures + a plain-language summary of exactly what
                                26 PVOD patients can and cannot support
gdf15_lr.py               —     importable module: fit_lr_model(), lr(x), lr_ci(x),
                                post_test(pre, x), with the cap and flags built in
```

Dependencies: `numpy pandas scipy scikit-learn statsmodels matplotlib openpyxl`. Set seeds (`np.random.seed(1)`) everywhere; bootstrap results must be reproducible.

---

## 7. Hard constraints

1. **No claim of clinical validity.** This is a demonstration of the evidential weight available in a 26-case sample, built for methodological argument. Header on every output: *not validated for clinical use; requires prospective validation in an independent cohort.*
2. **Never report an uncapped LR above the §4.1 ceiling.**
3. **Never report a probability without its interval.**
4. **Never use healthy controls in a diagnostic metric.**
5. **Never impute the clinical variables.** Report n for every analysis touching them.
6. **Every [VALIDATE] number becomes an `assert`.** If the pipeline cannot reproduce them, it fails loudly rather than producing a plausible wrong answer.
7. **Label reconstructed data as reconstructed.** The age↔GDF15 linkage (§1.3) is inferred, not published. Flag it wherever it is used.
8. **Single-cohort limitation.** All 26 PVOD patients come from UK NCS (19), UCSF PHB (4), and UCSF LTR (3). There is no external validation set and none can be constructed from this data. The tool cannot estimate its own generalisation error, and must say so.

---

## 8. Expected headline findings

The notebook should build to these conclusions:

1. The paper's 1,658 pg/mL cutoff carries **LR ≈ 2.4** — weak evidence, far weaker than "81% sensitivity, 98% specificity" implies to a reader.
2. GDF15 is a **much stronger rule-out than rule-in** test. Below 500 pg/mL, LR ≈ 0.06 with a tight CI — the most reliable single finding in the dataset, and one the paper does not emphasise.
3. Above ~2,500 pg/mL, the honest claim is **LR ≥ ~30**, not the 156 the point estimate suggests. Sample size, not biology, is the binding constraint.
4. GDF15 is **statistically orthogonal** to the haemodynamic and natriuretic-peptide markers already in use, and outperforms all of them for this specific distinction by a wide margin.
5. Therefore: 26 PVOD patients are enough to establish GDF15 as a **genuinely additive** diagnostic input worth prospective study — and not enough to establish a deployable threshold.
