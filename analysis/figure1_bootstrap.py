"""Figure 1 stratified case bootstrap; independent of the historical parametric bands.

Every sample is fit with standardized-predictor Newton maximum likelihood. Acceptance
requires overlapping classes, full rank, solver convergence, finite parameters,
positive-definite observed information and max absolute analytic score < 1e-7.
No resampling to replace rejected draws, no slope cap, no penalization.
"""
import hashlib
import json
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.special import expit
import data

SEED = 20260907
N_BOOT = 4000
GRID = np.logspace(np.log10(700), np.log10(6500), 400)
KEYS = np.array([861.,1658.,1809.5,2156.,2500.,3000.,4000.,5986.])

METHODS = (
    'Nine of 4,000 stratified bootstrap resamples produced complete separation, '
    'for which the maximum likelihood estimate is not finite. Percentile intervals '
    'for the Figure 1 band were computed from the 3,991 resamples with finite fits. '
    'Excluding the nine leaves the upper edge effectively unchanged at and below '
    'the published cutoff and narrows it within the extrapolated tail, by '
    'approximately 2% at 2,156 pg/mL and 17% at 5,000 pg/mL.'
)

def cohort(workbook=None):
    frame=data.load_gdf15(workbook)
    groups=[frame.loc[frame.cohort==g,'gdf15'].to_numpy() for g in data.GROUP1]
    return np.concatenate(groups),np.r_[np.ones(len(groups[0]),int),np.zeros(len(groups[1])+len(groups[2]),int)],groups

def samples(y):
    rng=np.random.default_rng(SEED)
    case,comp=np.flatnonzero(y==1),np.flatnonzero(y==0)
    for i in range(N_BOOT):
        state=rng.bit_generator.state
        take=np.r_[rng.choice(case,len(case),True),rng.choice(comp,len(comp),True)]
        yield i,state,take

def validated_fit(x,y):
    lx=np.log10(x); center=lx.mean();scale=lx.std()
    if not scale>0:raise ValueError('constant predictor / non-identifiability')
    if x[y==1].min()>=x[y==0].max() or x[y==1].max()<=x[y==0].min():
        raise ValueError('complete or quasi-complete separation')
    design=np.column_stack([np.ones(len(x)),(lx-center)/scale])
    if np.linalg.matrix_rank(design)!=2:raise ValueError('rank-deficient design')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fit=sm.Logit(y,design).fit(method='newton',disp=0,maxiter=100,tol=1e-11)
    eta=design@fit.params;p=expit(eta)
    score=design.T@(p-y);info=design.T@((p*(1-p))[:,None]*design)
    eig=np.linalg.eigvalsh(info)
    if not fit.mle_retvals['converged']:raise ValueError('Newton nonconvergence')
    if not np.all(np.isfinite(fit.params)) or not np.all(np.isfinite(eig)):raise ValueError('nonfinite fit')
    if eig.min()<=0 or np.linalg.norm(score,np.inf)>=1e-7:raise ValueError('score/information validation failed')
    coefficients=np.array([fit.params[0]-fit.params[1]*center/scale,fit.params[1]/scale])
    return coefficients,{'objective':float(np.sum(np.logaddexp(0,eta)-y*eta)),
        'standardized_score_max_abs':float(np.linalg.norm(score,np.inf)),
        'information_min_eigenvalue':float(eig.min()),'converged':True}

def lr(coefficients,grid,offset):
    b=np.asarray(coefficients)
    return np.exp(b[...,0,None]+b[...,1,None]*np.log10(grid)-offset)

def bootstrap(x,y):
    coefficients=[];records=[]
    for i,state,take in samples(y):
        rec={'replicate':i,'sample_sha256':hashlib.sha256(take.astype('<i8').tobytes()).hexdigest()}
        try:
            b,diagnostic=validated_fit(x[take],y[take])
            coefficients.append(b);rec.update(accepted=True,coefficients=b.tolist(),**diagnostic)
        except ValueError as e:
            # Only separation is approved for exclusion; any other failure is fatal.
            if str(e) != 'complete or quasi-complete separation':
                raise
            rec.update(accepted=False,reason=str(e),rng_state=state,indices=take.tolist())
        records.append(rec)
    if not coefficients:raise ValueError('No valid bootstrap fits')
    return np.array(coefficients),records


def write(results_dir, workbook=None):
    """Production band from public inputs; retain only aggregate/fit diagnostics."""
    dest = Path(results_dir)
    result = json.loads((dest / 'results.json').read_text())['models']['IPAH+HPAH']
    x, y, _ = cohort(workbook)
    coefficients, records = bootstrap(x, y)
    separated = [r['replicate'] for r in records if not r['accepted']]
    assert len(records) == N_BOOT and len(coefficients) == 3991 and len(separated) == 9
    offset = np.log(y.sum() / (len(y) - y.sum()))
    bounds = np.percentile(lr(coefficients, GRID, offset), [2.5, 97.5], axis=0)
    point = [result['intercept'], result['slope_log10_gdf15']]
    pd.DataFrame({'gdf15': GRID, 'lr': lr(point, GRID, offset),
                  'lr_ci_low': bounds[0], 'lr_ci_high': bounds[1]}).to_csv(
                      dest / 'figure1_band.csv', index=False)
    pd.DataFrame(coefficients, columns=['intercept', 'slope_log10_gdf15']).to_csv(
        dest / 'figure1_bootstrap_coefficients.csv', index=False)
    report = {'seed': SEED, 'requested': N_BOOT, 'finite': len(coefficients),
              'separated': len(separated), 'separated_replicates_zero_based': separated,
              'estimator': 'ordinary unpenalized logistic maximum likelihood',
              'percentiles': [2.5, 97.5], 'percentile_scale': 'LR',
              'methods_wording': METHODS,
              'sample_hashes_sha256': hashlib.sha256(''.join(
                  r['sample_sha256'] for r in records).encode()).hexdigest(),
              'max_standardized_score': max(r['standardized_score_max_abs']
                  for r in records if r['accepted']),
              'min_information_eigenvalue': min(r['information_min_eigenvalue']
                  for r in records if r['accepted'])}
    (dest / 'figure1_bootstrap.json').write_text(json.dumps(report, indent=2) + '\n')
    return report
