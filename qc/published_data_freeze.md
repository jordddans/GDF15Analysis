# Published-data analysis freeze

The permanent full historical record is the unchanged `analysis-v1.0-published-data` tag. Current `main` retains the audit's required code, public sources, manuscript inputs and verification/reporting evidence; superseded development material is available through that tag. This report's scientific conclusions are unchanged by the branch cleanup.

The approved Figure 1 bootstrap is integrated. All 211 verification checks pass, with zero failures: 109 full-workflow checks, 84 independent bootstrap-audit checks, and 18 freeze/document/privacy checks. Separately, all 685 baseline numerical/Boolean fields are unchanged.

## Bootstrap and Figure 1

4,000 stratified resamples requested, seed 20260907; 3,991 finite ordinary logistic maximum-likelihood fits; 9 complete-separation exclusions. The nine are not replaced. Non-separation fitting failures stop the pipeline. No penalization or second estimator was added.

Giorgio's three-panel script is adapted in `analysis/make_figure1.py`, called by the production workflow. The exact supplied script is archived unchanged. The retired root compatibility entry point remains only in frozen history. The percentile band is generated deterministically from all and only the accepted fits; counts and coefficients are saved separately. The point curve uses the pre-existing verified coefficients. Group counts, medians, both support boundaries and the real adjacent-patient annotations are computed from public observations. There are no hypothetical patients.

## Exact Methods wording used

“Nine of 4,000 stratified bootstrap resamples produced complete separation, for which the maximum likelihood estimate is not finite. Percentile intervals for the Figure 1 band were computed from the 3,991 resamples with finite fits. Excluding the nine leaves the upper edge effectively unchanged at and below the published cutoff and narrows it within the extrapolated tail, by approximately 2% at 2,156 pg/mL and 17% at 5,000 pg/mL.”

This wording appears in the QC DOCX and generated Figure 1 caption. The Methods paragraph also distinguishes the existing AUC/paired-landmark stratified procedures and retained parametric supplementary/comparator intervals. Only QC status, bootstrap Methods and Figure 1 caption text changed relative to the preceding QC manuscript; the reproduced image replaced the supplied preview. Giorgio's original DOCX is unchanged.

## Agreement with Giorgio's note

The upper edge is effectively unchanged at 1,000 pg/mL and the cutoff and narrows in the extrapolated tail. Exact signed upper-edge differences against all 4,000 original BFGS outputs are:

| GDF15 (pg/mL) | Upper-edge difference |
|---:|---:|
| 1,000 | +0.00430% |
| 1,658 | +0.00457% |
| 2,156 | -1.69832% |
| 3,000 | -6.59088% |
| 5,000 | -16.65908% |

These support the approximate 2%, 6% and 17% tail reductions. “Unchanged” is not literal for the entire interval: the lower LR endpoint moves from 0.011615 to 0.014001 at 1,000 pg/mL and from 1.238952 to 1.265842 at the cutoff. The Methods wording therefore specifies the **upper edge effectively unchanged**.

The note's LR=10 endpoints **1,913–2,365 pg/mL do not reproduce exactly** with the supplied script/seed. Its all-4,000 interval is 1914.109618–2363.412551; the finite-fit interval is 1913.854906–2363.420522. The exclusion shifts endpoints by less than 1 pg/mL. The previously verified production/manuscript landmark uses the existing seed-1 paired comparison and remains 2155.886177 pg/mL, interval 1922.636169–2356.969036 (manuscript rounds to 1,920–2,360). No landmark was overwritten to force agreement with the email. This is a documented numerical qualification, not a new statistical design.

## Reproduction and invariants

A fresh Python environment installed `requirements-lock.txt`. An isolated directory received only the analysis code and original published supporting workbook, with no private workbook, cached results, or archived patient datasets. The full workflow and its internal repeats passed. A second independent full run reproduced all 38 analysis/figure/reconstruction artifacts byte for byte. The QC DOCX also reproduces byte for byte. The local Word-rendered manuscript PDF is a review derivative.

Table 1, main coefficients, AUC, LR landmarks, calibration, internal validation, complete-case/missingness and LOO results are unchanged. Figures 2 and 3 remain byte-identical; Figure 2's caption only corrects its outdated Figure 1 bootstrap cross-reference. Both original and cleaned public workbook readers agree. The protected original DOCX and published workbooks retain their hashes. Historical submitted PDFs remain unchanged in the frozen tag and are intentionally absent from current `main`; the regression checker checks protected files present in the current tree. Historical paths in `baseline_manifest.json` describe the baseline, not additional current reproduction inputs.

Evidence: `published_data_clean_run.json`, `published_data_verification.json`, `baseline_regression.json`, `../results/verification.json`, `figure1_bootstrap/verification.json`, and `figure1_bootstrap/giorgio_note_comparison.json`.

## Privacy, version and remote

The private reconstructed workbook remains ignored and untracked and is absent from all local Git history. Five previously tracked public-source patient-level CSV exports have been removed from the frozen tree and added to ignore rules; local copies remain recoverable. After the subsequent authorized privacy cleanup, those exports are absent from every reachable commit and tag; obsolete backup refs, reflogs and objects were removed. See `history_rewrite.json`. No generated patient-level dataset is added in this freeze; the published source workbook and scientific tables/audit summaries remain reproducible.

The final local commit is to be identified by `git rev-parse analysis-v1.0-published-data^{commit}` after the verified commit/tag step. Only `analysis-v1.0-published-data` is authorized; no plain `analysis-v1.0`, v2.0 or Zenodo release. Jordan has approved `https://github.com/jordddans/GDF15Analysis`, independently confirmed private. Push is conditional on passing the repeated scientific and history checks. No public publication or Zenodo release is authorized.

## Future validation and remaining author items

Development loading (`data.py`) is separate from fitting and scoring (`model.py`). A future validation dataset can have a separate input loader and be scored with frozen coefficients; the current cohort-count assertions prevent silently changing the development cohort. No future UCSF samples were analyzed or simulated.

Giorgio's note closes COI, corresponding author (Akiko Hata), cohort-date/guideline and ethics/declaration decisions. They were not reopened or rewritten in this narrow pass. Assay methods, blinding confirmation, Amit's author-list position and the participant-count discrepancy remain author-side submission items. The published workbook's internally consistent counts are preserved. The email landmark discrepancy above remains explicitly recorded for authors. The authorized history-cleanup and private-push verification are documented separately in `history_rewrite.json`.

The original freeze changed-file inventory remains at `analysis-v1.0-published-data:qc/published_data_files_changed.txt`. Earlier QC reports are historical and superseded for Figure 1/freeze status. Historical paths embedded in retained verification/provenance records describe their original runs, rather than current-tree links.
