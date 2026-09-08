# QA Review — Supporting Data, Prabhakar et al. 2025

**File reviewed:** *Supporting data values_199159JCIRLRV2_Prabhakar et al 2025.xlsx* (18 sheets)
**Date:** July 18, 2026
**Accompanying file:** *Supporting data values_Prabhakar et al 2025_CLEANED.xlsx*

## Bottom line

The dataset is structurally sound and internally consistent. There is **one item that should be verified against the original lab data before use**, plus a few cosmetic fixes that have already been applied in the cleaned copy. Nothing else in the workbook was altered.

## 1. Needs verification — possible copy/paste artifact (GDF15-7)

On the **GDF15-7** tab, second table (*Plasma GDF15 amount, ng/ml*, MMC time course, rows 17–28):

- The **Day 4** row and the **Day 8** row share an exact run of **10 identical values**: cells **M24:V24** are identical to **H25:Q25** (17.47, 11.685, 18.05, 17.46, 16.13, 16.205, 16.16, 18.295, 18.345, 18.335).
- Within the **Day 16** row, a block of 5 values repeats itself: **M26:Q26** equals **R26:V26** (22.23, 22.2, 22.92, 22.97, 22.95).

Ten measured values matching to three decimals across two different timepoints does not occur by chance; this is the signature of an autofill or copy/paste slip. These cells are **highlighted in yellow and commented** in the cleaned copy. The values were **left unchanged** — they should be checked against the source records and corrected there. A full workbook scan found **no other duplicated runs** of this kind.

## 2. Cosmetic fixes already applied in the cleaned copy

- **GDF15-6**, header: "Specificity at 100% Sans­itivity" → **"Sensitivity"** (GDF15-3 already spells it correctly).
- **GDF15-7**, header: "RelativeGDF15 mRNA amount" → **"Relative GDF15 mRNA amount"** (missing space).
- **22 text labels** used non-breaking spaces instead of normal spaces (e.g., "Ctrl vs PVOD", "MMC Duration (days)"). Invisible on screen, but they break text searches and downstream parsing — now normalized.

## 3. Context worth noting (not errors)

- **Pediatric patients are present.** Ages include 1, ~14.7, ~15.7, and 16 years. This is consistent with a heritable/pediatric PVOD-PAH cohort and explains the very low 6-minute-walk distances (down to 15 m, mostly PVOD/HPAH). If pediatric inclusion is intended, these values are internally consistent; if not, the age-1 row and the 15–30 m walk distances are the first to re-check.
- **Missing-data completeness varies by measure** (coded consistently as "NA"): NT-proBNP ~62% missing, CI ~51%, 6MWD ~33%, PVR ~14%, mPAP ~5%. Normal for retrospective clinical data, but worth stating the denominators alongside any summary statistics.

## 4. Verified clean

- No formulas, no error cells (#REF!/#DIV/0!), no hidden rows or columns, no embedded charts or macros.
- Strong internal consistency:
  - The Ctrl and PVOD reference columns are **identical everywhere they are reused** (GDF15-1, GDF15-2, GDF15-4, GDF15-5).
  - The male + female split in **GDF15-10 reproduces the pooled GDF15-1 values** (differences ≤ 0.1, from rounding only), and group counts add up exactly: Ctrl 68, PVOD 26, IPAH 98, HPAH 42, Other PH 52, COPD 36.
  - The age vector is **identical across all six distribution sheets** (Age-GDF15, NT-proBNP, PVR, mPAP, CI, 6MWD).
- The repeated Ctrl / PVOD / "Age-adjusted GDF15" columns are **by design** (paired ROC and correlation analyses), not accidental duplicates.

---
*Prepared as a pre-share data check. The only substantive action item is verifying the highlighted GDF15-7 cells against the original data; all other changes are cosmetic and already applied in the cleaned copy.*
