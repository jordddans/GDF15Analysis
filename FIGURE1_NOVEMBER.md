# November Figure 1 update

Run `python analysis/make_figure1.py` from the repository root to render the November manuscript layout from the archived bootstrap coefficients. This plotting command performs no new fitting or resampling.

Outputs in `figures/`: `GDF15_Figure1.pdf`, `GDF15_Figure1.tif` (600 dpi, LZW), the existing `figure1_lr_curve.pdf/png` names, and `figure1_validation.json`.

The renderer checks seed 20260907, 4,000 requested resamples, 3,991 finite fits and 9 separation exclusions. It checks 34 differences greater than 10 percentage points: 17 within observed overlap and 17 above it. Curves are dashed outside the exact observed overlap, whose labels are rounded to 861–1,809.5 pg/ml. Labels say model-estimated; panel titles are neutral.

The existing full pipeline still regenerates analyses and bootstrap fits. Use the standalone command above when reproducing the figure from the archived run.

Validation: rendered successfully from the archived inputs; inputs match GitHub main by Git blob hashes. TIFF is 4320 × 5580 pixels at 600 dpi. The full historical verification suite was not rerun for this layout update.

Public release is pending. This repository contains manuscript drafts and internal review materials, including in its history; review these before changing visibility. No Zenodo DOI has been created by this update.
