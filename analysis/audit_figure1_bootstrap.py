"""Reproduce Giorgio's exact BFGS recipe and audit against independent Newton fits."""
from pathlib import Path
import json
import hashlib
import warnings
import statsmodels.api as sm
import pandas as pd
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
import figure1_bootstrap as fb

ROOT=Path(__file__).resolve().parents[1]
def main():
    dest=ROOT/'qc/figure1_bootstrap';dest.mkdir(exist_ok=True)
    x,y,_=fb.cohort();offset=np.log(y.sum()/(len(y)-y.sum()))
    grid=np.r_[fb.GRID,fb.KEYS];original=[];validated=[];records=[];all_records=[]
    for i,state,take in fb.samples(y):
        xx=x[take];yy=y[take];lx=np.log10(xx)
        f=lambda b:-np.sum(yy*(b[0]+b[1]*lx)-np.log1p(np.exp(b[0]+b[1]*lx)))
        fit=minimize(f,[0,1],method='BFGS')
        original.append(fit.x)
        try:
            b,diagnostic=fb.validated_fit(xx,yy);validated.append(b)
            delta=np.abs(fb.lr(b,grid,offset)/fb.lr(fit.x,grid,offset)-1)
            disposition='A' if delta.max()<1e-3 else 'B'
            record={'validated_coefficients':b.tolist(),'validated_diagnostics':diagnostic,
                    'max_absolute_coefficient_difference':float(np.max(np.abs(b-fit.x))),
                    'max_relative_lr_difference':float(delta.max()),'disposition':disposition}
        except (ValueError,np.linalg.LinAlgError) as e:
            # Diagnostic independent IRLS attempt; a separated sample has no finite MLE,
            # even when the solver returns a large finite approximation with success.
            design=np.column_stack([np.ones(len(xx)),(lx-lx.mean())/lx.std()])
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                independent=sm.GLM(yy,design,family=sm.families.Binomial()).fit(maxiter=100,tol=1e-11)
            ib=np.array([independent.params[0]-independent.params[1]*lx.mean()/lx.std(),independent.params[1]/lx.std()])
            record={'disposition':'C','reason':str(e),
                'case_min':float(xx[yy==1].min()),'comparator_max':float(xx[yy==0].max()),
                'separation_gap':float(xx[yy==1].min()-xx[yy==0].max()),
                'independent_IRLS_diagnostic_only':{'coefficients':ib.tolist(),
                  'converged':bool(independent.converged),'iterations':int(independent.fit_history['iteration']),
                  'warnings':sorted(set(str(w.message) for w in caught)),
                  'accepted':False,'reason':'Strict separation proves no finite unpenalized MLE; finite solver iterates are not estimates'}}
        all_records.append({'replicate':i,'sample_sha256':hashlib.sha256(take.astype('<i8').tobytes()).hexdigest(),
            'original_success':bool(fit.success),'original_coefficients':fit.x.tolist(),**record})
        if not fit.success or record['disposition']!='A':
            score=np.column_stack([np.ones(len(xx)),lx]).T@(expit(fit.x[0]+fit.x[1]*lx)-yy)
            record.update(replicate=i,seed=fb.SEED,rng_state=state,indices=take.tolist(),
                sample_sha256=hashlib.sha256(take.astype('<i8').tobytes()).hexdigest(),
                original_coefficients=fit.x.tolist(),original_objective=float(fit.fun),
                original_gradient_norm=float(np.linalg.norm(fit.jac)),analytic_gradient_norm=float(np.linalg.norm(score)),
                original_status=int(fit.status),original_success=bool(fit.success),original_message=str(fit.message),
                original_lr_at_keys=fb.lr(fit.x,fb.KEYS,offset).tolist(),
                validated_lr_at_keys=fb.lr(b,fb.KEYS,offset).tolist() if record['disposition']!='C' else None)
            records.append(record)
    old=np.percentile(fb.lr(np.array(original),grid,offset),[2.5,97.5],axis=0)
    new=np.percentile(fb.lr(np.array(validated),grid,offset),[2.5,97.5],axis=0)
    absolute=np.abs(new-old);relative=absolute/old
    summary={'seed':fb.SEED,'requested':fb.N_BOOT,'accepted':len(validated),'rejected':fb.N_BOOT-len(validated),
        'original_warning_count':sum(not r['original_success'] for r in records),
        'warning_dispositions':{c:sum(not r['original_success'] and r['disposition']==c for r in records) for c in 'ABC'},
        'all_samples_refit':fb.N_BOOT,'warning_fits_repaired':sum(not r['original_success'] and r['disposition']!='C' for r in records),
        'materiality_relative_lr_tolerance':1e-3,'max_absolute_band_difference':float(absolute.max()),
        'max_relative_band_difference':float(relative.max()),
        'max_warning_relative_lr_difference':max(r.get('max_relative_lr_difference',0) for r in records if not r['original_success']),
        'band_materially_changed':bool(relative.max()>=1e-3),
        'key_concentrations':[{'gdf15':float(v),'original':old[:,-len(fb.KEYS)+j].tolist(),
            'validated':new[:,-len(fb.KEYS)+j].tolist(),'absolute_difference':absolute[:,-len(fb.KEYS)+j].tolist(),
            'relative_difference':relative[:,-len(fb.KEYS)+j].tolist()} for j,v in enumerate(fb.KEYS)]}
    extrema={}
    for label,arr in [('absolute',absolute),('relative',relative)]:
        endpoint,j=np.unravel_index(arr.argmax(),arr.shape)
        extrema[label]={'endpoint':['lower','upper'][endpoint],'gdf15':float(grid[j]),'difference':float(arr[endpoint,j])}
    summary['band_difference_maxima']=extrema
    summary['all_sample_dispositions']={c:sum(r['disposition']==c for r in all_records) for c in 'ABC'}
    summary['nonwarning_material_refits']=[r['replicate'] for r in all_records if r['original_success'] and r['disposition']=='B']
    (dest/'all_fits.json').write_text(json.dumps(all_records,indent=2)+'\n')
    pd.DataFrame({'gdf15':grid,'original_low':old[0],'original_high':old[1],
        'validated_low':new[0],'validated_high':new[1],
        'absolute_low_difference':absolute[0],'absolute_high_difference':absolute[1],
        'relative_low_difference':relative[0],'relative_high_difference':relative[1]}).to_csv(dest/'band_comparison.csv',index=False)
    (dest/'warning_fits.json').write_text(json.dumps(records,indent=2)+'\n')
    (dest/'band_comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    assert summary['accepted'] == 3991 and summary['rejected'] == 9
if __name__=='__main__':main()
