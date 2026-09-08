"""Independent calibration, AUC, and leave-one-case-out manuscript checks."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from openpyxl import load_workbook
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import beta

ROOT=Path(__file__).resolve().parents[1]

def fit(x,y):
    l=np.log10(x);mean=l.mean();sd=l.std();d=np.column_stack([np.ones(len(x)),(l-mean)/sd])
    r=minimize(lambda b:np.sum(np.logaddexp(0,d@b)-y*(d@b)),[0.,0.],
               jac=lambda b:d.T@(expit(d@b)-y),method='BFGS',options={'gtol':1e-9})
    assert max(abs(d.T@(expit(d@r.x)-y)))<1e-6
    return np.array([r.x[0]-r.x[1]*mean/sd,r.x[1]/sd])

def audit(workbook,results):
    results=Path(results);w=load_workbook(workbook,data_only=True);s=w['Fig_1a' if 'Fig_1a' in w.sheetnames else 'GDF15-1']
    raw=[(float(s.cell(r,c).value),c==2,s.cell(r,c).coordinate) for r in range(4,102) for c in (2,3,4) if isinstance(s.cell(r,c).value,(int,float))]
    x=np.array([r[0] for r in raw]);y=np.array([r[1] for r in raw],int)
    co=fit(x,y);pred=expit(co[0]+co[1]*np.log10(x));main=json.loads((results/'results.json').read_text())['models']['IPAH+HPAH']
    checks={};quintiles=[]
    for i,chunk in enumerate(np.array_split(np.argsort(x),5)):
        n=len(chunk);events=int(sum(y[chunk]));p=float(pred[chunk].mean());obs=events/n
        saved=main['calibration_quintiles'][i]
        low=0 if events==0 else beta.ppf(.025,events,n-events+1)
        high=1 if events==n else beta.ppf(.975,events+1,n-events)
        checks[f'calibration_quintile_{i+1}']=bool(n==saved['n'] and events==saved['observed_events'] and
            np.allclose([p,obs,low,high],[saved['predicted_mean'],saved['observed_proportion'],saved['observed_ci_low'],saved['observed_ci_high']],rtol=0,atol=1e-8))
        quintiles.append({'quintile':i+1,'n':n,'predicted':p,'observed':obs})
    cases=x[y==1];controls=x[y==0]
    kernel=(cases[:,None]>controls).astype(float)+.5*(cases[:,None]==controls)
    checks['auc_pairwise']=abs(kernel.mean()-main['auc'])<1e-14
    # Independent implementation of stratified AUC bootstrap in the same source-row order.
    rng=np.random.default_rng(1);ix1=np.flatnonzero(y);ix0=np.flatnonzero(1-y);areas=[]
    for _ in range(4000):
        a=x[rng.choice(ix1,len(ix1),replace=True)];b=x[rng.choice(ix0,len(ix0),replace=True)]
        areas.append(np.mean((a[:,None]>b)+.5*(a[:,None]==b)))
    ci=np.percentile(areas,[2.5,97.5]);checks['auc_bootstrap_4000']=bool(np.allclose(ci,main['auc_ci'],rtol=0,atol=1e-14))
    points=np.array([1658.,2500.,3000.,4000.,5986.41]);offset=np.log(26/140)
    # Use independently fitted coefficients, not rounded manuscript coefficients.
    lr=np.exp(co[0]+co[1]*np.log10(points)-offset);deleted=[]
    for idx in ix1:
        take=np.arange(len(x))!=idx;b=fit(x[take],y[take])
        deleted.append(np.exp(b[0]+b[1]*np.log10(points)-np.log(25/140)))
    deleted=np.array(deleted);change=100*(deleted/lr-1);rows=[]
    manuscript=[8.1,31.3,53.5,96.3,177.]
    for j,p in enumerate(points):
        worst=int(np.argmax(abs(change[:,j])))
        value=float(abs(change[worst,j]));shown=round(value,1) if j<4 else round(value)
        checks[f'loo_manuscript_{p:g}']=shown==manuscript[j]
        rows.append({'gdf15':p,'independent_full_lr':lr[j],'max_absolute_percent_change':value,
                     'signed_percent_change':change[worst,j],'excluded_source_cell':s.title+'!'+raw[ix1[worst]][2]})
    checks['comparator_support_max']=round(float(controls.max()),1)==1809.5
    cc=json.loads((results/'results.json').read_text())['marker_auc_complete_case']
    checks['enrichment_435pct_28fold']=round(100*10/23,1)==43.5 and round((10/23)/(26/166),1)==2.8 and cc['n']==23 and cc['n_pvod']==10
    out={'status':'PASS' if all(checks.values()) else 'FAIL','checks':{k:bool(v) for k,v in checks.items()},
         'calibration':quintiles,'auc':float(kernel.mean()),'auc_bootstrap_ci':ci.tolist(),
         'loo':rows,'maximum_observed_concentration':float(x.max()),
         'note':'5986 in the manuscript labels the observed maximum 5986.41; no clinical Patient_ID is assigned.'}
    (results/'manuscript_number_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    pd.DataFrame(rows).to_csv(results/'loo_manuscript_landmarks.csv',index=False)
    assert all(checks.values()),checks
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--workbook',type=Path,default=ROOT/'archive/gdfdiscrim/jci-136-199159-s068.xlsx')
    p.add_argument('--results-dir',type=Path,default=ROOT/'results');a=p.parse_args()
    out=audit(a.workbook,a.results_dir);print(out['status'],len(out['checks']),'checks')
