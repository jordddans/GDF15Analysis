# Source-only reconstruction and privacy

`analysis/reconstruct.py` reads the published `jci-136-199159-s068.xlsx`. A checked sheet-name mapping also accepts the cleaned workbook. Neither workbook is modified. All 4,609 numeric cells agree at the same coordinates across corresponding sheets; text labels are not identical.

The original `Fig_1g` has no diagnosis text in column A. The historical reconstruction assigns rows 4–29 to PVOD, 30–127 to IPAH, and 128–169 to HPAH. This explicit block convention creates Patient_ID 1–166. The cleaned workbook supplies those labels, and the private workbook independently agrees on every diagnosis and ID. This is verification-based reconstruction metadata, not a claim that the published sheet contains those strings. The PVOD block additionally matches all 26 age-adjusted keys in `Fig_S2b`; assignment of the two comparator blocks retains the established reconstruction convention.

`Fig_1h`–`Fig_1m` share exactly the same age vector. Age alone is not unique. Their order is matched to `Fig_1g` using age-adjusted GDF15 after restoring the median, 0.19355683, derived from `Fig_1g`. Every match must be unique within 1e-6, and the match must be one-to-one. Missing or ambiguous matches fail. Independently, paired 6MWD values differ by exactly 294 m across all diagnosis groups.

The four correlation panels of `Fig_S2b` are matched by their respective age-adjusted keys. Raw values are used when tabulated there. Otherwise mPAP, PVR and NT-proBNP are restored from deviations using offsets 51, 871.925 and 1264, derived and checked from overlapping values. PVR deviations have 0.0001-level rounding; NT-proBNP has one 0.2-level rounding difference. Historical reconstructed PVR values are rounded to three decimals; measured PVR retains its source precision. Cardiac index remains a deviation, not a reconstructed absolute value.

The two negative reconstructed NT-proBNP values are preserved as reconstruction artifacts, never interpreted as concentrations. Model comparison uses the source distribution values, not these reconstructed magnitudes. Missingness counts from either source (65 NT-proBNP measurements overall; 17 PVOD) remain separate from distribution-only counts (63 overall; 15 PVOD) and the 23-patient complete-case subset.

Local comparison passed all 22 columns across `Patients_166` and `Source_values`. IDs, diagnoses, missingness, age, age-adjusted GDF15, and all but one numerical field agree exactly; maximum NT-proBNP floating residual is 8.53e-14. The tolerance is 1e-12 with no relative tolerance. Full comparison details and generated data are in ignored `generated/`.

No raw GDF15 is joined to a clinical identifier. Raw-observation outputs use source cells, and LOO case numbers refer only to source-column positions. The historical `archive/gdfdiscrim/data/df_g1.csv` contains an unsupported raw-GDF15 patient linkage and is excluded from git, as is `pvod_166.csv`. Exploratory notebooks and chat logs are conservatively excluded because they may embed these exports. The local originals remain intact.

To reconstruct: `python analysis/reconstruct.py --workbook path/to/published.xlsx --output-dir generated`.
The optional `--compare-private` switch is only for the local verification custodian. Do not obtain, distribute, commit, or upload the private workbook to reproduce the project.
