"""Write the narrow published-data handoff from completed verification artifacts."""
import json
from pathlib import Path
import subprocess
from figure1_bootstrap import METHODS

ROOT = Path(__file__).resolve().parents[1]


def main():
    verified = json.loads((ROOT/'qc/published_data_verification.json').read_text())
    clean = json.loads((ROOT/'qc/published_data_clean_run.json').read_text())
    note = json.loads((ROOT/'qc/figure1_bootstrap/giorgio_note_comparison.json').read_text())
    main_model = json.loads((ROOT/'results/results.json').read_text())['models']['IPAH+HPAH']
    assert verified['status'] == 'PASS' and clean['failures'] == 0
    changed = subprocess.check_output(['git','diff','--name-only','HEAD'],cwd=ROOT,text=True).splitlines()
    added = subprocess.check_output(['git','ls-files','--others','--exclude-standard'],cwd=ROOT,text=True).splitlines()
    files = sorted(set(changed + added + ['qc/published_data_freeze.md','qc/published_data_files_changed.txt']))
    (ROOT/'qc/published_data_files_changed.txt').write_text('\n'.join(files)+'\n')
    rows = '\n'.join(f"| {g:,.0f} | {note['signed_percent_difference'][1][i]:+.5f}% |"
                     for i,g in enumerate(note['concentrations']))
    landmark = main_model['gdf15_at_lr']['10']
    report = f'''# Published-data analysis freeze

The approved Figure 1 bootstrap is integrated. All {verified['total_checks']} verification checks pass, with zero failures: {verified['workflow_checks']} full-workflow checks, {verified['bootstrap_audit_checks']} independent bootstrap-audit checks, and {verified['count']} freeze/document/privacy checks. Separately, all 685 baseline numerical/Boolean fields are unchanged.

## Bootstrap and Figure 1

4,000 stratified resamples requested, seed 20260907; 3,991 finite ordinary logistic maximum-likelihood fits; 9 complete-separation exclusions. The nine are not replaced. Non-separation fitting failures stop the pipeline. No penalization or second estimator was added.

Giorgio's three-panel script is adapted in `analysis/make_figure1.py`, called by the production workflow. The exact supplied script is archived unchanged. The root compatibility entry point now delegates to production. The percentile band is generated deterministically from all and only the accepted fits; counts and coefficients are saved separately. The point curve uses the pre-existing verified coefficients. Group counts, medians, both support boundaries and the real adjacent-patient annotations are computed from public observations. There are no hypothetical patients.

## Exact Methods wording used

“{METHODS}”

This wording appears in the QC DOCX and generated Figure 1 caption. The Methods paragraph also distinguishes the existing AUC/paired-landmark stratified procedures and retained parametric supplementary/comparator intervals. Only QC status, bootstrap Methods and Figure 1 caption text changed relative to the preceding QC manuscript; the reproduced image replaced the supplied preview. Giorgio's original DOCX is unchanged.

## Agreement with Giorgio's note

The upper edge is effectively unchanged at 1,000 pg/mL and the cutoff and narrows in the extrapolated tail. Exact signed upper-edge differences against all 4,000 original BFGS outputs are:

| GDF15 (pg/mL) | Upper-edge difference |
|---:|---:|
{rows}

These support the approximate 2%, 6% and 17% tail reductions. “Unchanged” is not literal for the entire interval: the lower LR endpoint moves from 0.011615 to 0.014001 at 1,000 pg/mL and from 1.238952 to 1.265842 at the cutoff. The Methods wording therefore specifies the **upper edge effectively unchanged**.

The note's LR=10 endpoints **1,913–2,365 pg/mL do not reproduce exactly** with the supplied script/seed. Its all-4,000 interval is {note['original_lr10_interval'][0]:.6f}–{note['original_lr10_interval'][1]:.6f}; the finite-fit interval is {note['finite_lr10_interval'][0]:.6f}–{note['finite_lr10_interval'][1]:.6f}. The exclusion shifts endpoints by less than 1 pg/mL. The previously verified production/manuscript landmark uses the existing seed-1 paired comparison and remains {landmark['gdf15']:.6f} pg/mL, interval {landmark['ci_low']:.6f}–{landmark['ci_high']:.6f} (manuscript rounds to 1,920–2,360). No landmark was overwritten to force agreement with the email. This is a documented numerical qualification, not a new statistical design.

## Reproduction and invariants

A fresh Python environment installed `requirements-lock.txt`. An isolated directory received only the analysis code and original published supporting workbook, with no private workbook, cached results, or archived patient datasets. The full workflow and its internal repeats passed. A second independent full run reproduced all {clean['repeated_artifact_count']} analysis/figure/reconstruction artifacts byte for byte. The QC DOCX also reproduces byte for byte. The local Word-rendered manuscript PDF is a review derivative.

Table 1, main coefficients, AUC, LR landmarks, calibration, internal validation, complete-case/missingness and LOO results are unchanged. Figures 2 and 3 remain byte-identical; Figure 2's caption only corrects its outdated Figure 1 bootstrap cross-reference. Both original and cleaned public workbook readers agree. The protected original DOCX, published workbooks and submitted PDFs retain their hashes.

Evidence: `published_data_clean_run.json`, `published_data_verification.json`, `baseline_regression.json`, `../results/verification.json`, `figure1_bootstrap/verification.json`, and `figure1_bootstrap/giorgio_note_comparison.json`.

## Privacy, version and remote

The private reconstructed workbook remains ignored and untracked and is absent from all local Git history. Five previously tracked public-source patient-level CSV exports have been removed from the frozen tree and added to ignore rules; local copies remain recoverable. After the subsequent authorized privacy cleanup, those exports are absent from every reachable commit and tag; obsolete backup refs, reflogs and objects were removed. See `history_rewrite.json`. No generated patient-level dataset is added in this freeze; the published source workbook and scientific tables/audit summaries remain reproducible.

The final local commit is to be identified by `git rev-parse analysis-v1.0-published-data^{{commit}}` after the verified commit/tag step. Only `analysis-v1.0-published-data` is authorized; no plain `analysis-v1.0`, v2.0 or Zenodo release. Jordan has approved `https://github.com/jordddans/GDF15Analysis`, independently confirmed private. Push is conditional on passing the repeated scientific and history checks. No public publication or Zenodo release is authorized.

## Future validation and remaining author items

Development loading (`data.py`) is separate from fitting and scoring (`model.py`). A future validation dataset can have a separate input loader and be scored with frozen coefficients; the current cohort-count assertions prevent silently changing the development cohort. No future UCSF samples were analyzed or simulated.

Giorgio's note closes COI, corresponding author (Akiko Hata), cohort-date/guideline and ethics/declaration decisions. They were not reopened or rewritten in this narrow pass. Assay methods, blinding confirmation, Amit's author-list position and the participant-count discrepancy remain author-side submission items. The published workbook's internally consistent counts are preserved. The email landmark discrepancy above remains explicitly recorded for authors. The authorized history-cleanup and private-push verification are documented separately in `history_rewrite.json`.

Full changed-file inventory: [published_data_files_changed.txt](published_data_files_changed.txt). Earlier QC reports are historical and superseded for Figure 1/freeze status.
'''
    (ROOT/'qc/published_data_freeze.md').write_text(report)


if __name__ == '__main__':
    main()
