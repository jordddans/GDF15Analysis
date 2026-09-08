"""Independent source-cell verification; does not import production model/data helpers."""
import argparse
from decimal import Decimal
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from scipy.optimize import minimize
from scipy.special import expit

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'archive/gdfdiscrim/jci-136-199159-s068.xlsx'


def verify(workbook=SOURCE, results=ROOT/'results'):
    results=Path(results);results.mkdir(parents=True,exist_ok=True)
    w=load_workbook(workbook,data_only=True)
    sheet=w['Fig_1a' if 'Fig_1a' in w.sheetnames else 'GDF15-1']
    # Source first sheet: retain XLSX numeric tokens rather than rounded display text.
    ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with ZipFile(workbook) as z:
        xml=ET.fromstring(z.read('xl/worksheets/sheet1.xml'))
    tokens={c.attrib['r']:c.find('s:v',ns).text for c in xml.findall('.//s:c',ns)
            if c.find('s:v',ns) is not None and c.attrib.get('t','n')=='n'}
    observations=[]
    for col in range(1,7):
        diagnosis=str(sheet.cell(3,col).value).replace('\xa0',' ')
        for row in range(4,sheet.max_row+1):
            c=sheet.cell(row,col)
            if isinstance(c.value,(int,float)):
                observations.append({'diagnosis':diagnosis,'raw_concentration':tokens[c.coordinate],
                                     'source_cell':f'{sheet.title}!{c.coordinate}',
                                     'source_column':c.column_letter})
    pv=[r for r in observations if r['diagnosis']=='PVOD']
    comp=[r for r in observations if r['diagnosis'] in ('IPAH','HPAH')]
    dec=lambda r:Decimal(r['raw_concentration'])
    pairs=[{'pvod_source_cell':a['source_cell'],'comparator_source_cell':b['source_cell'],
            'absolute_gap':str(abs(dec(a)-dec(b)))} for a in pv for b in comp]
    minimum=min(Decimal(r['absolute_gap']) for r in pairs)
    minima=[r for r in pairs if Decimal(r['absolute_gap'])==minimum]
    selected={minima[0]['pvod_source_cell'],minima[0]['comparator_source_cell']}
    cohort=comp+pv; x=np.array([float(dec(r)) for r in cohort]);y=np.array([r['diagnosis']=='PVOD' for r in cohort],int)
    predictor=np.log10(x);mean=predictor.mean();sd=predictor.std()
    design=np.column_stack([np.ones(len(x)),(predictor-mean)/sd])
    fit=minimize(lambda b:float(np.sum(np.logaddexp(0,design@b)-y*(design@b))),
                 [0.,0.],jac=lambda b:design.T@(expit(design@b)-y),method='BFGS',options={'gtol':1e-9})
    assert np.max(np.abs(design.T@(expit(design@fit.x)-y)))<1e-6
    fitted=np.array([fit.x[0]-fit.x[1]*mean/sd,fit.x[1]/sd])
    saved=pd.read_csv(results/'model_coefficients.csv').set_index('model').loc['IPAH+HPAH']
    coefficients=np.array([saved['intercept'],saved['slope_log10_gdf15']])
    cutoff=Decimal('1658');positive=x>=float(cutoff)
    tp=int(sum(positive&(y==1)));fn=int(sum(~positive&(y==1)))
    fp=int(sum(positive&(y==0)));tn=int(sum(~positive&(y==0)))
    sensitivity=tp/(tp+fn);specificity=tn/(tn+fp)
    lrplus=sensitivity/(1-specificity);lrminus=(1-sensitivity)/specificity
    offset=np.log(len(pv)/len(comp))
    rows=[]
    for r in cohort:
        v=float(dec(r));lr=float(np.exp(coefficients[0]+coefficients[1]*np.log10(v)-offset))
        rule=lrplus if dec(r)>=cutoff else lrminus
        p=lr/(9+lr);pr=rule/(9+rule)
        rows.append(dict(r,continuous_lr=lr,continuous_post_test_probability=p,rule_lr=rule,
                         rule_post_test_probability=pr,absolute_percentage_point_difference=100*abs(p-pr),
                         cutoff_positive=dec(r)>=cutoff))
    adjacent=[dict(r,duplicate_source_cells=[o['source_cell'] for o in observations
                                            if dec(o)==dec(r) and o['source_cell']!=r['source_cell']])
              for r in rows if r['source_cell'] in selected]
    adjacent.sort(key=lambda r:dec(r))
    checks={
      'cohort_counts':len(pv)==26 and len(comp)==140,
      'all_cross_group_pairs_evaluated':len(pairs)==len(pv)*len(comp),
      'unique_global_minimum':len(minima)==1,
      'distinct_concentrations':minimum>0,
      'no_cross_column_duplicates_for_pair':all(not r['duplicate_source_cells'] for r in adjacent),
      'cutoff_between_observations':dec(adjacent[0])<cutoff<=dec(adjacent[1]),
      'pair_classified_correctly':adjacent[0]['diagnosis']=='IPAH' and adjacent[1]['diagnosis']=='PVOD',
      'coefficients_independently_reproduce':bool(np.allclose(fitted,coefficients,rtol=0,atol=1e-6)),
      'counts_match_manuscript':(tp,fn,fp,tn)==(21,5,4,136),
      'correct_157':tp+tn==157,
      'difference_over_10_is_34':sum(r['absolute_percentage_point_difference']>10 for r in rows)==34,
      'difference_over_20_is_11':sum(r['absolute_percentage_point_difference']>20 for r in rows)==11,
      'both_continuous_probabilities_round_to_22_percent':all(round(r['continuous_post_test_probability']*100)==22 for r in adjacent),
      'rule_lr_minus_rounds_to_020':round(lrminus,2)==0.20,
      'rule_lr_plus_rounds_to_283':round(lrplus,1)==28.3,
      'rule_probabilities_round_to_2_and_76':round(lrminus/(9+lrminus)*100)==2 and round(lrplus/(9+lrplus)*100)==76,
      'raw_observations_have_no_patient_id':all('patient_id' not in r and 'Patient_ID' not in r for r in rows),
    }
    out={'status':'PASS' if all(checks.values()) else 'FAIL','source_sha256':hashlib.sha256(Path(workbook).read_bytes()).hexdigest(),
         'source_sheet':sheet.title,'cutoff':str(cutoff),'pre_test_probability':0.1,
         'sampling_prevalence':len(pv)/len(cohort),'sampling_log_odds':float(offset),
         'final_coefficients':coefficients.tolist(),'independent_coefficients':fitted.tolist(),
         'absolute_gap':str(minimum),'gap_percent_of_cutoff':float(minimum/cutoff*100),
         'number_of_pairs_evaluated':len(pairs),'global_minimum_pairs':minima,'adjacent_observations':adjacent,
         'confusion_matrix':{'tp':tp,'fn':fn,'fp':fp,'tn':tn},'correctly_classified':tp+tn,
         'sensitivity':sensitivity,'specificity':specificity,'lr_plus':lrplus,'lr_minus':lrminus,
         'continuous_lr_at_cutoff':float(np.exp(coefficients[0]+coefficients[1]*np.log10(float(cutoff))-offset)),
         'difference_over_10_pp':sum(r['absolute_percentage_point_difference']>10 for r in rows),
         'difference_over_20_pp':sum(r['absolute_percentage_point_difference']>20 for r in rows),'checks':checks}
    (results/'centerpiece_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    pd.DataFrame(rows).to_csv(results/'cutoff_probability_comparison.csv',index=False)
    pd.DataFrame(pairs).to_csv(results/'cross_group_gaps.csv',index=False)
    pd.DataFrame(adjacent).drop(columns=['duplicate_source_cells']).to_csv(results/'table1.csv',index=False)
    lines=['| Diagnosis | Measured GDF15 (pg/mL) | Continuous LR | Post-test at 10% | Cutoff-rule LR | Rule post-test |',
           '|---|---:|---:|---:|---:|---:|']
    for r in adjacent:
        lines.append(f"| {r['diagnosis']} | {float(dec(r)):,.2f} | {r['continuous_lr']:.2f} | {r['continuous_post_test_probability']:.0%} | {r['rule_lr']:.2f} | {r['rule_post_test_probability']:.0%} |")
    (results/'table1.md').write_text('\n'.join(lines)+'\n\nRows identify anonymous source observations, not reconstructed Patient_IDs. Full precision and source cells are in table1.csv.\n')
    assert all(checks.values()),checks
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--workbook',type=Path,default=SOURCE)
    p.add_argument('--results-dir',type=Path,default=ROOT/'results');a=p.parse_args()
    r=verify(a.workbook,a.results_dir);print(r['status'],len(r['checks']),'checks; gap',r['absolute_gap'])
