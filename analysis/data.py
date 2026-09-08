"""Build the analysis dataset from the JCI supporting-data workbook.

Source: Prabhakar et al., J Clin Invest 2026;136(2):e199159, Supporting data values
(``Supporting data values_Prabhakar et al 2025_CLEANED.xlsx``).  Every value analysed here
is read from that workbook.

A second input, ``Reconstructed_per_patient_Prabhakar_2025_2.xlsx``, supplies the patient
identifier and nothing else.  It is a per-patient reconstruction of the same study,
assembled by matching the workbook's figure-panel sheets, and its ``Patients_166`` tab is
the only place in this project where one row is one patient under a stable name.  It is
used for exactly two things: to attach ``Patient_ID`` to the clinical frame, and to count
which measurements exist per diagnosis.  No analysed quantity is taken from it -- see
:func:`load_patients_166` for why.

Frames produced:

``load_gdf15``    322 raw plasma GDF15 values with cohort labels, from sheet ``GDF15-1``.
                  This is the only input the likelihood-ratio models use.  The sheet is
                  six anonymous cohort columns: it carries no patient identifier, no age
                  and no row order shared with anything else, so nothing in it can be
                  attached to a named patient.

``load_clinical`` the 166 group-1 patients (PVOD, IPAH, HPAH) with ``patient_id``, age,
                  age-adjusted GDF15 and the five clinical measurements, from sheets
                  ``6MWD-1``, ``Age-GDF15``, ``NT-proBNP``, ``PVR``, ``mPAP``, ``CI`` and
                  ``6MWD-2``.  Used for the marker-comparison table only.

``load_patients_166``
                  the reconstruction's per-patient tab, carrying ``Patient_ID`` 1-166.

``load_ntprobnp_correlation``
                  the 26 PVOD rows of the PVOD-only correlation sheet
                  ``6MWD PVR mPAP NT-proBNP``, keyed to their ``NT-proBNP``
                  distribution-sheet rows, so that the NT-proBNP denominators of the two
                  sheets can be stated per patient.

Sheet layout notes that matter for parsing:

* Column headers contain non-breaking spaces (``Other\\xa0PH``); they are normalised.
* Missing values are the string ``"NA"``.
* ``6MWD-1`` carries the group label for each of the 166 patients in column A
  (PVOD rows 4-29, IPAH rows 30-127, HPAH rows 128-169) and the age-adjusted
  GDF15 value in column B.
* The six distribution sheets share one patient order, which is *not* the order of
  ``6MWD-1``.  ``Age-GDF15`` column C holds the age-adjusted GDF15 value expressed as
  a distance from the cohort median (0.19355683); adding that constant back recovers
  the ``6MWD-1`` values to 7e-8 and keys every distribution row to its group label.
* The five clinical sheets likewise store each measurement as a distance from its
  cohort median.  Every statistic computed here is rank-based, so the centring
  constant cancels; the values are kept as published rather than un-centred.
* ``Patient_ID`` is an identifier created for the reconstruction, not a clinical
  identifier from the contributing centres.  It is blocked by diagnosis (1-26 PVOD,
  27-124 IPAH, 125-166 HPAH), so an id implies a group and no figure should be ordered
  by it.  It is nonetheless the only unambiguous handle on a patient in this project:
  age repeats, group repeats, and age-adjusted GDF15, though unique, is a measurement
  and so gets transformed, rounded and filtered.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook

ROOT = Path(__file__).resolve().parents[1]

WORKBOOK = ROOT / "Supporting data values_Prabhakar et al 2025_CLEANED.xlsx"

# the per-patient reconstruction: source of Patient_ID and of the completeness counts only
RECONSTRUCTED = ROOT / "Reconstructed_per_patient_Prabhakar_2025_2.xlsx"

COHORTS = ["Ctrl", "PVOD", "IPAH", "HPAH", "Other PH", "COPD"]
GROUP1 = ["PVOD", "IPAH", "HPAH"]

# Patients_166 column -> the measurement it records.  Availability of these five defines
# the complete-case subset; their per-group counts are the completeness table.
RECONSTRUCTED_MEASUREMENTS = {
    "CI_distFromMedian": "Cardiac index",
    "NTproBNP_ng_L": "NT-proBNP",
    "6MWD_m": "6MWD",
    "PVR_dyn_s_cm5": "PVR",
    "mPAP_mmHg": "mPAP",
}

# Patient_ID blocks, as published in the reconstruction's READ ME tab
ID_BLOCKS = {"PVOD": (1, 26), "IPAH": (27, 124), "HPAH": (125, 166)}

# age-adjusted GDF15 median used by the workbook's distribution sheets
AGEADJ_MEDIAN = 0.19355683

# distribution sheet -> analysis column name
CLINICAL_SHEETS = {
    "NT-proBNP": "ntprobnp",
    "PVR": "pvr",
    "mPAP": "mpap",
    "CI": "ci",
    "6MWD-2": "sixmwd",
}


def _num(value):
    """Workbook cell -> float, with 'NA' and blanks as NaN."""
    if value is None:
        return np.nan
    text = str(value).strip()
    if text in ("", "NA", "None"):
        return np.nan
    try:
        return float(text)
    except ValueError:
        return np.nan


def load_gdf15(workbook=WORKBOOK):
    """322 raw plasma GDF15 values (pg/mL) with cohort labels, from sheet GDF15-1."""
    sheet = load_workbook(workbook, data_only=True)["GDF15-1"]
    headers = {
        cell.column - 1: str(cell.value).replace("\xa0", " ")
        for cell in sheet[3]
        if cell.value is not None
    }
    rows = []
    for row in sheet.iter_rows(min_row=4, max_row=sheet.max_row,
                               min_col=1, max_col=6, values_only=True):
        for index, value in enumerate(row):
            if index in headers and value is not None:
                rows.append({"cohort": headers[index], "gdf15": float(value)})
    df = pd.DataFrame(rows)
    counts = df.cohort.value_counts().to_dict()
    expected = {"Ctrl": 68, "PVOD": 26, "IPAH": 98, "HPAH": 42, "Other PH": 52, "COPD": 36}
    assert counts == expected, f"GDF15-1 cohort counts {counts} != {expected}"
    return df


def load_patients_166(reconstructed=RECONSTRUCTED):
    """The reconstruction's per-patient tab: one row per patient, carrying ``Patient_ID``.

    Two things this frame is used for, and nothing else.

    ``Patient_ID`` is the identifier.  It is the only stable handle on a patient in this
    project, and :func:`load_clinical` joins it onto the analysis frame.

    The presence or absence of each measurement is the completeness table.  Only presence
    is read: the reconstruction recovers raw values by adding a median back to a published
    deviation, which is exact for the join but leaves two NT-proBNP cells negative
    (``Patient_ID`` 38 and 58) where the deviation sat below the median by more than the
    median itself.  A negative concentration is not a concentration, so **no numeric value
    from this file is analysed**.  Every AUC in this project is computed from the original
    workbook's own distribution values; this file answers only "was this patient measured".
    """
    frame = pd.read_excel(reconstructed, sheet_name="Patients_166")

    expected = ["Patient_ID", "Group", "Age_years", "AgeAdj_GDF15"]
    missing = [column for column in expected + list(RECONSTRUCTED_MEASUREMENTS)
               if column not in frame.columns]
    assert not missing, f"Patients_166 is missing columns {missing}"

    assert len(frame) == 166, f"Patients_166 has {len(frame)} rows, expected 166"
    assert frame.Patient_ID.is_unique, "Patient_ID is not unique"
    assert frame.Patient_ID.tolist() == list(range(1, 167)), "Patient_ID is not 1-166"

    counts = frame.Group.value_counts().to_dict()
    assert counts == {"IPAH": 98, "HPAH": 42, "PVOD": 26}, counts
    for group, (first, last) in ID_BLOCKS.items():
        block = frame.Patient_ID[frame.Group == group]
        assert block.min() == first and block.max() == last and len(block) == last - first + 1, \
            f"{group} is not the contiguous Patient_ID block {first}-{last}"
    return frame


def load_clinical(workbook=WORKBOOK, reconstructed=RECONSTRUCTED):
    """The 166 group-1 patients with ``Patient_ID``, age, age-adjusted GDF15 and measurements."""
    book = load_workbook(workbook, data_only=True)

    labelled = book["6MWD-1"]
    group, ageadj, sixmwd_raw = [], [], []
    current = None
    for row in range(4, 170):
        label = labelled.cell(row, 1).value
        if label:
            current = str(label).strip()
        group.append(current)
        ageadj.append(_num(labelled.cell(row, 2).value))
        sixmwd_raw.append(_num(labelled.cell(row, 3).value))
    labels = pd.DataFrame({"group": group, "ageadj_gdf15": ageadj,
                           "sixmwd_m": sixmwd_raw})
    counts = labels.group.value_counts().to_dict()
    assert counts == {"IPAH": 98, "HPAH": 42, "PVOD": 26}, counts

    ages = book["Age-GDF15"]
    age = np.array([_num(ages.cell(row, 2).value) for row in range(5, 171)])
    key = np.array([_num(ages.cell(row, 3).value) for row in range(5, 171)]) + AGEADJ_MEDIAN

    # key each distribution row to its labelled row by age-adjusted GDF15 value
    order = []
    for value in key:
        match = np.argmin(np.abs(labels.ageadj_gdf15.to_numpy() - value))
        assert abs(labels.ageadj_gdf15.iloc[match] - value) < 1e-6, "age-adjusted GDF15 key failed"
        order.append(match)
    assert len(set(order)) == 166, "age-adjusted GDF15 key is not one-to-one"

    df = labels.iloc[order].reset_index(drop=True)
    df["age_years"] = age
    for sheet_name, column in CLINICAL_SHEETS.items():
        sheet = book[sheet_name]
        df[column] = [_num(sheet.cell(row, 3).value) for row in range(5, 171)]

    # end-to-end check on the row keying: 6MWD-1 holds raw 6MWD, 6MWD-2 holds the same
    # measurement centred on the cohort median (294 m), so the two must differ by exactly
    # that constant on every patient once the two sheet orders are aligned.
    offsets = (df.sixmwd_m - df.sixmwd).dropna().unique()
    assert len(offsets) == 1 and abs(offsets[0] - 294.0) < 1e-9, f"6MWD keying failed: {offsets}"

    present = df.groupby("group")[list(CLINICAL_SHEETS.values())].count()
    assert present.loc["PVOD", "ntprobnp"] == 15, "NT-proBNP PVOD count != 15"
    assert df.ntprobnp.notna().sum() == 63, "NT-proBNP n != 63"

    df.insert(0, "patient_id", _patient_ids_for(df, reconstructed))
    return df


def _patient_ids_for(clinical, reconstructed=RECONSTRUCTED):
    """``Patient_ID`` for each row of the clinical frame, joined on age-adjusted GDF15.

    Age-adjusted GDF15 is distinct across all 166 patients in both files, so it keys the
    two to each other exactly -- the join is asserted below to be one-to-one, to match to
    the last bit rather than to a tolerance, and to agree with the ``6MWD-1`` group label
    on every row.

    The join is by key and never by position.  The clinical frame is in the distribution
    sheets' order, which runs HPAH, then IPAH, then PVOD; ``Patient_ID`` is in the
    ``6MWD-1`` order, which runs PVOD first.  Lining the two up by row would silently
    mislabel all 166 patients, which is the whole reason the key is used.
    """
    patients = load_patients_166(reconstructed)
    ids = patients.Patient_ID.to_numpy()
    groups = patients.Group.to_numpy()
    keys = patients.AgeAdj_GDF15.to_numpy()

    matched, seen = [], set()
    for value, group in zip(clinical.ageadj_gdf15, clinical.group):
        position = int(np.argmin(np.abs(keys - value)))
        assert keys[position] == value, (
            f"age-adjusted GDF15 {value!r} is not in Patients_166 exactly "
            f"(closest {keys[position]!r})")
        assert groups[position] == group, (
            f"Patient_ID {ids[position]} is {groups[position]} in Patients_166 "
            f"but {group} in 6MWD-1")
        assert position not in seen, f"Patient_ID {ids[position]} matched twice"
        seen.add(position)
        matched.append(int(ids[position]))

    assert len(seen) == 166, f"the join reached {len(seen)} of 166 patients"
    return matched


def load_ntprobnp_correlation(workbook=WORKBOOK):
    """The 26 PVOD rows of the correlation sheet, keyed to the ``NT-proBNP`` sheet.

    ``6MWD PVR mPAP NT-proBNP`` is a PVOD-only sheet holding the pairs used for the
    published correlation analyses.  Its 26 rows are the PVOD patients in the order of the
    ``6MWD-1`` PVOD block; the workbook carries no patient identifier of its own, so each
    row is keyed to its ``NT-proBNP`` distribution-sheet row through the age-adjusted
    GDF15 value, the same key :func:`load_clinical` uses, and reported under the
    ``patient_id`` that key resolves to.  ``pvod_id``, the position in the PVOD block, is
    kept beside it because earlier write-ups quote it, but ``patient_id`` is the
    identifier to use.

    The two sheets do not carry the same number of NT-proBNP measurements for PVOD.  This
    function reports both, per patient, so the denominators can be stated rather than
    reconciled; it takes no position on which sheet a downstream analysis should use.
    """
    book = load_workbook(workbook, data_only=True)
    clinical = load_clinical(workbook)
    pvod = clinical[clinical.group == "PVOD"].reset_index(drop=True)

    sheet = book["6MWD PVR mPAP NT-proBNP"]
    rows = []
    for offset in range(26):
        ageadj = _num(sheet.cell(5 + offset, 10).value)
        match = np.argmin(np.abs(pvod.ageadj_gdf15.to_numpy() - ageadj))
        assert abs(pvod.ageadj_gdf15.iloc[match] - ageadj) < 1e-6, "correlation-sheet key failed"
        rows.append({
            "patient_id": int(pvod.patient_id.iloc[match]),
            "pvod_id": match + 1,
            "ageadj_gdf15": ageadj,
            "ntprobnp_correlation_sheet": _num(sheet.cell(5 + offset, 11).value),
            "ntprobnp_distribution_sheet": pvod.ntprobnp.iloc[match],
        })
    df = pd.DataFrame(rows).sort_values("patient_id").reset_index(drop=True)
    assert len(set(df.pvod_id)) == 26, "correlation-sheet key is not one-to-one"
    assert len(set(df.patient_id)) == 26, "patient_id is not one-to-one on the PVOD block"
    # the PVOD Patient_ID block is 1-26 in the 6MWD-1 order, which is the order pvod_id
    # counts, so the two coincide here; asserted rather than assumed
    assert (df.patient_id == df.pvod_id).all(), "patient_id and pvod_id disagree on PVOD"

    # The correlation sheet carries raw NT-proBNP and the distribution sheet the same
    # measurement centred on the cohort median, so wherever both record a patient the two
    # must differ by that one constant.  The distribution sheet is written to whole ng/L,
    # so the single non-integer measurement (396.8) differs by 1,263.8 rather than 1,264;
    # the tolerance below is that rounding and nothing more.
    offsets = (df.ntprobnp_correlation_sheet - df.ntprobnp_distribution_sheet).dropna()
    median_offset = float(offsets.median())
    assert (offsets - median_offset).abs().max() < 0.5, \
        f"NT-proBNP sheets disagree by more than rounding: {sorted(offsets.unique())}"
    df.attrs["median_offset"] = median_offset
    df.attrs["max_rounding_difference"] = float((offsets - median_offset).abs().max())

    df["in_correlation_sheet"] = df.ntprobnp_correlation_sheet.notna()
    df["in_distribution_sheet"] = df.ntprobnp_distribution_sheet.notna()
    assert df.in_distribution_sheet.sum() == 15, "distribution-sheet PVOD count != 15"
    assert df.in_correlation_sheet.sum() == 17, "correlation-sheet PVOD count != 17"
    assert not (df.in_distribution_sheet & ~df.in_correlation_sheet).any(), \
        "a distribution-sheet value has no correlation-sheet counterpart"
    return df
