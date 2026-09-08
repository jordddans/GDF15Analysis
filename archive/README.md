# Archive

Superseded exploratory work. **Exception: `gdfdiscrim/jci-136-199159-s068.xlsx` is the unmodified published workbook and is now the active source input.** Other material is historical. It is kept because it is the record of how the final
model was arrived at, not because anything depends on it.

| Path | What it is |
|---|---|
| `gdfdiscrim/` | The first pass: exploratory notebooks, three earlier figure scripts, the BayLearn abstract source, and correspondence. Superseded by `analysis/`. |
| `figures 2/` | Earlier figure drafts, superseded by `figures/`. |
| `Professor notes.pages` | Supervisor notes from the exploratory phase. |

To reproduce the reported analysis, ignore this directory entirely and follow the
quick start in the repository's top-level `README.md`.

## QC distinctions

Patient-level reconstructed CSVs, exploratory notebooks, and chat logs are local-only and excluded from git. In particular `gdfdiscrim/data/df_g1.csv` contains an unsupported association of raw GDF15 with patient IDs; it must not be used as a patient-level source. See ../qc/reconstruction.md.

The two PDFs in `figures 2/figure3/` show LR− 0.24 for a cross-validated, reselected cutoff procedure (six false negatives). This is distinct from LR− 0.20 for the fixed published cutoff (five false negatives). Historical PDFs remain unchanged; their values must not be imported into current fixed-cutoff results.
