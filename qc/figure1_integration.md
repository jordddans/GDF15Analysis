# Figure 1 integration — approved finite-fit production

Giorgio approved the nine complete-separation exclusions on 11 September 2026. Production now uses 4,000 requested stratified draws (seed 20260907), 3,991 validated finite ordinary logistic MLEs and 9 separated exclusions. No replacement draws or second estimator are used; any failure other than separation stops production.

`analysis/make_figure1.py` adapts Giorgio's three-panel script. The exact supplied source is preserved at `figure1_bootstrap/giorgio_original.py` (SHA-256 e01247c941dbaf64f26fed177411007eb8bfc58610a2237163ee886e088ece45). The root compatibility entry point now runs production; the archived all-4,000 implementation is only audit provenance.

`analysis/figure1_bootstrap.py` generates `results/figure1_band.csv`, `figure1_bootstrap_coefficients.csv`, and `figure1_bootstrap.json`. Plotting uses the existing verified main-model coefficients and computed counts, support limits, medians and adjacent observed pair. Original Fig_1a and cleaned GDF15-1 inputs are supported by the shared reader. Both extrapolated regions are marked. No hypothetical patients or private reconstruction are used by Figure 1.

The QC DOCX embeds the new production PNG and discloses the approved handling in Methods and caption. The supplied original DOCX remains unchanged. The upper edge narrows in the extrapolated tail, approximately as in Giorgio's note; exact lower-edge and landmark qualifications are in [the freeze report](published_data_freeze.md). Historical audit conclusions below the dated notice in `figure1_bootstrap/report.md` record the earlier gate and are superseded by this approval.
