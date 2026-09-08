# Audit of 0.24 and variants

Found **394 literal matches in 38 files**, including decimal substrings in measurements, coefficients and SVG geometry. All occurrences are listed in audit_024.csv. Full contexts involving private clinical rows are retained only in ignored qc/private/.

**No incorrect fixed-cutoff LR− occurrence was found; no numeric replacement was made.** Two historical PDF labels show 0.24 for a different, cross-validated procedure, correctly rounded from (6/26)/(136/140) = 0.2375565611. The fixed published cutoff instead gives (5/26)/(136/140) = 0.1979638009, displayed as 0.20. The archived caption identifies nested leave-one-out threshold reselection and six false negatives. Neither PDF was altered.

Text extraction covered PDFs, DOCX components, spreadsheet cells, source code, notes, notebook source/outputs, chat, and decompressed Pages IWA content; no extraction errors. Local OCR inspected 24 raster images/previews, with 0 additional matches. The original draft, README, BayLearn reconciliation, submitted PDFs and relevant captions were included even where no match was found.

Clarifying notes added to README.md, archive/README.md and BayLearn_reconciliation.md prevent reuse of the historical CV estimand as the fixed-cutoff LR−. The ledger includes legitimate post-test probabilities around 0.244 at 1700 pg/mL, source clinical values, log-likelihoods, coefficients, continuous LR/CI values and vector coordinates. These were left unchanged.

| File | Occurrences |
|---|---:|
| Reconstructed_per_patient_Prabhakar_2025_2.xlsx | 15 |
| Supporting data values_Prabhakar et al 2025_CLEANED.xlsx | 39 |
| TRIPOD-AI_mapping.md | 2 |
| archive/figures 2/figure1/Figure1_auc_comparison.svg | 11 |
| archive/figures 2/figure2/Fig4_gdf15_cutoffs_scatter_linear_revised.svg | 5 |
| archive/figures 2/figure3/figure_gdf15_classifier.pdf | 1 |
| archive/figures 2/figure3/figure_gdf15_classifier_captioned.pdf | 1 |
| archive/gdfdiscrim/.ipynb_checkpoints/04_honesty-checkpoint.ipynb | 1 |
| archive/gdfdiscrim/04_honesty.ipynb | 1 |
| archive/gdfdiscrim/GitHub_scripts_PVOD_project.R | 1 |
| archive/gdfdiscrim/JournalFigures/figure1/Figure1_auc_comparison.svg | 11 |
| archive/gdfdiscrim/JournalFigures/figure2/Fig4_gdf15_cutoffs_scatter_linear.svg | 5 |
| archive/gdfdiscrim/JournalFigures/figure3/Figure3_classifier.svg | 21 |
| archive/gdfdiscrim/data/df_all.csv | 2 |
| archive/gdfdiscrim/data/df_g1.csv | 4 |
| archive/gdfdiscrim/gdf15_charts/Fig1_gdf15_cutoffs_linear.svg | 16 |
| archive/gdfdiscrim/gdf15_charts/Fig1_gdf15_cutoffs_log.svg | 14 |
| archive/gdfdiscrim/gdf15_charts/Fig2_gdf15_cutoffs_symbols_linear.svg | 4 |
| archive/gdfdiscrim/gdf15_charts/Fig2_gdf15_cutoffs_symbols_log.svg | 7 |
| archive/gdfdiscrim/gdf15_charts/Fig3_gdf15_cutoffs_ranked_linear.svg | 8 |
| archive/gdfdiscrim/gdf15_charts/Fig3_gdf15_cutoffs_ranked_log.svg | 5 |
| archive/gdfdiscrim/gdf15_charts/Fig4_gdf15_cutoffs_scatter_linear.svg | 5 |
| archive/gdfdiscrim/gdf15_charts/Fig4_gdf15_cutoffs_scatter_linear_bold.svg | 11 |
| archive/gdfdiscrim/gdf15_charts/Fig4_gdf15_cutoffs_scatter_log.svg | 6 |
| archive/gdfdiscrim/gdf15_charts/FigS1_gdf15_cutoffs_labelled_linear.svg | 17 |
| archive/gdfdiscrim/gdf15_charts/FigS1_gdf15_cutoffs_labelled_log.svg | 14 |
| archive/gdfdiscrim/gdf15_charts/gdf15_ladder_points.csv | 2 |
| archive/gdfdiscrim/jci-136-199159-s068.xlsx | 39 |
| archive/gdfdiscrim/pvod_166.csv | 4 |
| archive/gdfdiscrim/untitled.chat | 10 |
| figures/figure1_caption.md | 1 |
| results/comparator_models.csv | 1 |
| results/functional_form.csv | 1 |
| results/lr_curves.csv | 99 |
| results/lr_reference.csv | 1 |
| results/model_coefficients.csv | 2 |
| results/ntprobnp_sources.csv | 1 |
| results/results.json | 6 |
