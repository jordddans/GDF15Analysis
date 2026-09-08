"""Deterministic local clinical reconstruction; never joins anonymous raw GDF15."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'archive/gdfdiscrim/jci-136-199159-s068.xlsx'
ALIASES = dict(zip(
    ['GDF15-1','GDF15-2','GDF15-3','GDF15-4','GDF15-5','GDF15-6','6MWD-1',
     'Age-GDF15','NT-proBNP','PVR','mPAP','CI','6MWD-2','GDF15-7','GDF15-8',
     'GDF15-9','GDF15-10','6MWD PVR mPAP NT-proBNP'],
    ['Fig_1a','Fig_1b','Fig_1c','Fig_1d','Fig_1e','Fig_1f','Fig_1g','Fig_1h',
     'Fig_1i','Fig_1j','Fig_1k','Fig_1l','Fig_1m','Fig_S1a','Fig_S1b','Fig_S1c',
     'Fig_S2a','Fig_S2b']))


def book(path):
    w = load_workbook(path, data_only=True)
    return {alias: w[alias if alias in w.sheetnames else original]
            for alias, original in ALIASES.items()}


def num(v):
    return float(v) if isinstance(v, (int, float)) else np.nan


def unique_match(keys, value, tolerance=1e-6):
    matches = np.flatnonzero(np.abs(np.asarray(keys)-value) < tolerance)
    if len(matches) != 1:
        raise ValueError(f'Clinical key must match exactly one row; got {len(matches)}')
    return int(matches[0])


def reconstruct(path=SOURCE):
    w = book(path)
    labelled, ages, corr = w['6MWD-1'], w['Age-GDF15'], w['6MWD PVR mPAP NT-proBNP']
    keys = np.array([num(labelled.cell(r,2).value) for r in range(4,170)])
    assert np.isfinite(keys).all() and len(set(keys)) == 166
    # Existing reconstruction convention, verified against the known-good workbook.
    # Original Fig_1g has unlabeled blocks; do not claim these strings occur there.
    groups = ['PVOD']*26 + ['IPAH']*98 + ['HPAH']*42
    age_values = np.array([num(ages.cell(r,2).value) for r in range(5,171)])
    distkeys = np.array([num(ages.cell(r,3).value) for r in range(5,171)])
    median = float(np.median(keys))
    assert abs(median - 0.19355683) < 1e-12
    positions = [unique_match(keys, v+median) for v in distkeys]
    assert len(set(positions)) == 166
    for name in ['NT-proBNP','PVR','mPAP','CI','6MWD-2']:
        assert np.array_equal(age_values, [num(w[name].cell(r,2).value) for r in range(5,171)])
    distribution_row = {pos: i+5 for i,pos in enumerate(positions)}
    correlation_rows = {}
    for r in range(5,31):
        idx = unique_match(keys, num(corr.cell(r,1).value))
        assert idx < 26 and idx not in correlation_rows
        for c in (4,7,10):
            assert unique_match(keys,num(corr.cell(r,c).value)) == idx
        correlation_rows[idx] = r
    # Derive offsets from paired source values; source PVR deviations are rounded.
    definitions = [('mPAP',7,8,51.0,1e-9),('PVR',4,5,871.925,0.00011),
                   ('NT-proBNP',10,11,1264.0,0.5)]
    offsets = {}
    for name,_,c,expected,tol in definitions:
        differences=[]
        for idx,r in correlation_rows.items():
            a=num(corr.cell(r,c).value); b=num(w[name].cell(distribution_row[idx],3).value)
            if np.isfinite(a) and np.isfinite(b): differences.append(a-b)
        observed=float(np.median(differences))
        # Reconstructed workbook specifies PVR offset to three decimal places.
        offset=round(observed,3) if name=='PVR' else observed
        assert abs(offset-expected)<1e-9 and max(abs(v-offset) for v in differences)<tol
        offsets[name]=offset
    rows=[]; provenance=[]; source_rows=[]
    for idx,key in enumerate(keys):
        r=idx+4; d=distribution_row[idx]; cr=correlation_rows.get(idx)
        raw6=num(labelled.cell(r,3).value)
        dist={name:num(w[name].cell(d,3).value) for name in ['mPAP','PVR','NT-proBNP','CI','6MWD-2']}
        if np.isfinite(raw6): assert raw6-dist['6MWD-2']==294
        measured={name:num(corr.cell(cr,c).value) if cr else np.nan
                  for name,c in [('mPAP',8),('PVR',5),('NT-proBNP',11)]}
        if cr:
            assert (np.isnan(raw6) and np.isnan(num(corr.cell(cr,2).value))) or raw6==num(corr.cell(cr,2).value)
        row={'Patient_ID':idx+1,'Group':groups[idx],'Age_years':num(ages.cell(d,2).value),
             'AgeAdj_GDF15':key,'6MWD_m':raw6}
        for name,field in [('mPAP','mPAP_mmHg'),('PVR','PVR_dyn_s_cm5'),('NT-proBNP','NTproBNP_ng_L')]:
            row[field]=measured[name] if np.isfinite(measured[name]) else dist[name]+offsets[name]
            if name == 'PVR' and not np.isfinite(measured[name]):
                row[field] = round(row[field], 3)  # historical reconstruction precision
        row['CI_distFromMedian']=dist['CI'];rows.append(row)
        provenance.append({'Patient_ID':idx+1,'clinical_key_cell':f'{labelled.title}!B{r}',
            'distribution_row':d,'correlation_row':cr,'diagnosis_provenance':'verified clinical block convention; original Fig_1g labels absent',
            'key_residual':float(key-(num(ages.cell(d,3).value)+median))})
        source_rows.append([idx+1,groups[idx],row['Age_years'],key,raw6,measured['mPAP'],dist['mPAP'],
                            measured['PVR'],dist['PVR'],measured['NT-proBNP'],dist['NT-proBNP'],dist['CI'],dist['6MWD-2']])
    source_columns=['Patient_ID','Group','Age_years','AgeAdj_GDF15','6MWD_m\n(tab7 raw)',
                    'mPAP meas\n(tab18)','mPAP distFromMed\n(tab11)','PVR meas\n(tab18)',
                    'PVR distFromMed\n(tab10)','NTproBNP meas\n(tab18)','NTproBNP distFromMed\n(tab9)',
                    'CI distFromMed\n(tab12)','6MWD distFromMed\n(tab13)']
    return pd.DataFrame(rows),pd.DataFrame(source_rows,columns=source_columns),pd.DataFrame(provenance)


def write(path, output, compare=None):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    patients,sources,provenance=reconstruct(path)
    patients.to_csv(output/'patients_166.csv',index=False)
    sources.to_csv(output/'source_values.csv',index=False)
    provenance.to_csv(output/'provenance.csv',index=False)
    report={'source_sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),'rows':len(patients),
            'raw_gdf15_joined':False,'private_workbook_required':False,'comparison':None}
    if compare:
        checks=[]
        for name,actual in [('Patients_166',patients),('Source_values',sources)]:
            expected=pd.read_excel(compare,sheet_name=name)
            assert list(actual.columns)==list(expected.columns)
            for col in actual:
                a,b=actual[col],expected[col]
                assert a.isna().equals(b.isna()),(name,col,'missingness')
                if pd.api.types.is_numeric_dtype(a):
                    equal=np.allclose(a,b,rtol=0,atol=1e-12,equal_nan=True)
                    delta=float(np.nanmax(np.abs(a-b)))
                else:equal=a.equals(b);delta=None
                checks.append({'sheet':name,'field':col,'equal':bool(equal),'max_abs_difference':delta})
        report['comparison']={'status':'PASS' if all(c['equal'] for c in checks) else 'FAIL',
                              'numeric_tolerance':1e-12,'fields':checks}
        (output/'reconstruction_verification.json').write_text(json.dumps(report,indent=2)+'\n')
        assert report['comparison']['status']=='PASS',report['comparison']
    else:(output/'reconstruction_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workbook',type=Path,default=SOURCE)
    p.add_argument('--output-dir',type=Path,default=ROOT/'generated')
    p.add_argument('--compare-private',type=Path)
    a=p.parse_args();print(json.dumps(write(a.workbook,a.output_dir,a.compare_private),indent=2))
