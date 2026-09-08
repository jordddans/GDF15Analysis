# Figure 1 integration gate

BLOCKED: Giorgio's script is absent. The supplied DOCX contains one PNG and no embedded script, code attachment, or script hyperlink. Its image is preserved as `qc/figure1_supplied.png`. The original DOCX remains unchanged.

The user explicitly selected: retain the script gate; do not independently replace Giorgio's implementation; withhold analysis-v1.0 until this is resolved. The current reproducible historical Figure 1 remains in `figures/` and is preserved in baseline commit c4a5a20. The QC DOCX retains Giorgio's supplied image and is not a submission-ready replacement.

Independent source-cell checks reproduce panels B/C's point estimates and counts. This establishes numerical agreement with those claims, not reproducibility of the image or its uncertainty band.

Method conflict: the DOCX labels the band a 4,000-replicate stratified case bootstrap; current repository LR bands use a 4,000-replicate parametric bootstrap. Stratified case bootstrapping is used for AUC and the paired Youden/LR=10 stability comparison; ordinary case resampling is used for optimism correction. Existing numerical results and methods have been preserved. The image does not visibly mark the 1,809.46 pg/mL upper comparator support boundary despite the manuscript saying both support limits are marked. The colour-distance claim also requires its calculation or source.

When supplied: inspect the script; compare source selection, all point estimates, both bands/methods, seed and replicate handling against the independent verification. STOP on disagreement. Only after reconciliation integrate source/result-driven panels A–C, observed-point counts, medians, pair annotations, accessible shapes/lines and support boundaries, then rerun deterministic checks. Do not equate a raw observation's cell reference with Patient_ID.
