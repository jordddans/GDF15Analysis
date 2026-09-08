# Figure 1. Likelihood ratio for PVOD as a continuous function of plasma GDF15.

Logistic regression on log10(plasma GDF15), a single predictor, fitted to
166 patients — PVOD (n = 26) versus pooled IPAH
(n = 98) and HPAH (n = 42). The fitted log-odds are converted to a likelihood ratio by subtracting the
log-odds of the sampling prevalence of this cohort
(15.7%, -1.6835 on the log-odds scale),
which removes the cohort's case mix from the estimate; the resulting curve is a property of
the two conditional GDF15 distributions and carries no assumption about prevalence in any
clinical population. A clinician supplies the pre-test probability and multiplies its odds
by the likelihood ratio read from this curve. Intercept -39.247
(SE 7.757); log10(GDF15) coefficient 11.959
(SE 2.422, 95% CI 7.21 to 16.71);
26 events for one predictor.

The band is a pointwise 95% interval from a parametric bootstrap: outcomes were simulated
from the fitted model at the observed GDF15 values and the model refitted
4,000 times (seed 1), so the interval reflects
uncertainty in the fitted parameters and is defined across the whole measured range rather
than only where cases were observed. Horizontal rules mark likelihood ratios of 1, 2, 5 and
10. The likelihood ratio reaches 10 at 2,160 pg/mL
(95% CI 1,920–2,360). That interval is a stratified
case bootstrap — PVOD and comparator patients resampled separately at their observed sample sizes,
4,000 replicates, seed 1 — which is the scheme under which this
landmark and a selected threshold can be compared on identical footing (see below). The parametric
bootstrap gives 1,800–2,490 for the
same landmark and is reported as a supplementary sensitivity analysis. Bootstrap endpoints are
quoted to the nearest 10 pg/mL; `results/results.json` holds them at full precision.

Dashed lines mark the limits of the observed data. The region below
~861 pg/mL is labelled on the figure as model-supported beyond observed
PVOD data: no PVOD patient was measured below 860.513 pg/mL, so the curve
there is the fitted model form rather than observation. The region above
~1,809.5 pg/mL is labelled as model-supported beyond observed
IPAH/HPAH comparator data: no comparator patient was measured above that value, and
20 of 26 PVOD patients
(77%) lie above that value; the
concentration at which the likelihood ratio reaches 10 therefore falls in a region
containing no comparator observations, which is why the band widens there. Rug marks below
the axis show every measured value in each group. The open diamond marks
1,658 pg/mL, the cutoff reported for these data, which
carries a likelihood ratio of 2.6
(1.40–4.74). Read as a point on
this curve, that measurement takes a 20%
pre-test probability to 39.0%
(25.9%–54.2%);
read as a positive result from the dichotomy, whose likelihood ratio is
28.3, it would take the same patient to
87.6%. Both readings use the same patients and the
same outcome and differ only in what is done with the measurement.

Source data: Prabhakar et al., J Clin Invest 2026;136(2):e199159, Supporting data values, sheet GDF15-1.
