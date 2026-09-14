"""Check bootstrap audit evidence and the approved note against exact finite-fit results."""
from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
import figure1_bootstrap as fb
import data

ROOT=Path(__file__).resolve().parents[1]
def main():
    dest=ROOT/'qc/figure1_bootstrap'
    summary=json.loads((dest/'band_comparison.json').read_text())
    all_fits=json.loads((dest/'all_fits.json').read_text())
    exceptional=json.loads((dest/'warning_fits.json').read_text())
    x,y,_=fb.cohort();checks=[]
    def check(name,condition):
        checks.append({'name':name,'passed':bool(condition)})
    warnings=[r for r in exceptional if not r['original_success']]
    check('59 original warning fits reproduced',len(warnings)==59)
    check('all 4000 draws accounted for',len(all_fits)==4000 and [r['replicate'] for r in all_fits]==list(range(4000)))
    for r in warnings:
        rng=np.random.default_rng();rng.bit_generator.state=r['rng_state']
        case,comp=np.flatnonzero(y==1),np.flatnonzero(y==0)
        take=np.r_[rng.choice(case,len(case),True),rng.choice(comp,len(comp),True)]
        check(f"warning {r['replicate']}: exact sample replays from RNG state",
            take.tolist()==r['indices'] and hashlib.sha256(take.astype('<i8').tobytes()).hexdigest()==r['sample_sha256'])
    for r in exceptional:
        if r['disposition']=='C':
            take=np.array(r['indices']);xx=x[take];yy=y[take]
            check(f"rejected {r['replicate']}: strict separation independently confirmed",
                xx[yy==1].min()>xx[yy==0].max() and not r['independent_IRLS_diagnostic_only']['accepted'])
    accepted=[r for r in all_fits if r['disposition']!='C']
    check('accepted fits meet analytic score/information criteria',all(
        r['validated_diagnostics']['standardized_score_max_abs']<1e-7 and
        r['validated_diagnostics']['information_min_eigenvalue']>0 for r in accepted))
    check('warning dispositions 54 valid, zero materially different, five unusable',summary['warning_dispositions']=={'A':54,'B':0,'C':5})
    check('3991 accepted, nine rejected',len(accepted)==3991 and summary['rejected']==9)
    check('four separated fits incorrectly reported original success',sum(r['disposition']=='C' and r['original_success'] for r in all_fits)==4)
    check('successful BFGS fit 2296 requires numerical replacement',summary['nonwarning_material_refits']==[2296])
    grid=np.r_[fb.GRID,fb.KEYS];offset=np.log(y.sum()/(len(y)-y.sum()))
    expected=np.percentile(fb.lr(np.array([r['validated_coefficients'] for r in accepted]),grid,offset),[2.5,97.5],axis=0)
    saved=pd.read_csv(dest/'band_comparison.csv',float_precision='round_trip')
    check('band uses exactly accepted fits, LR-scale percentiles',np.allclose(expected,saved[['validated_low','validated_high']].to_numpy().T,rtol=1e-14,atol=0))
    repeated,records=fb.bootstrap(x,y)
    check('repeat independent fits byte-identical coefficients',np.array_equal(repeated,np.array([r['validated_coefficients'] for r in accepted])))
    check('repeat accepted/rejected sample identities identical',all(
        a['sample_sha256']==b['sample_sha256'] and a['accepted']==(b['disposition']!='C') for a,b in zip(records,all_fits)))
    original=data.load_gdf15(ROOT/'archive/gdfdiscrim/jci-136-199159-s068.xlsx')
    cleaned=data.load_gdf15(ROOT/'Supporting data values_Prabhakar et al 2025_CLEANED.xlsx')
    check('shared reader returns identical 322 source values and diagnoses',len(original)==322 and original.equals(cleaned))
    check('archived Giorgio script unchanged',hashlib.sha256((dest/'giorgio_original.py').read_bytes()).hexdigest()=='e01247c941dbaf64f26fed177411007eb8bfc58610a2237163ee886e088ece45')
    check('historical band difference retained transparently after approval',summary['band_materially_changed'] and summary['max_relative_band_difference']>.24)
    # Check Giorgio's rounded upper-edge note against the full-precision audit fits.
    note_grid = np.array([1000., 1658., 2156., 3000., 5000.])
    old = np.percentile(fb.lr(np.array([r['original_coefficients'] for r in all_fits]), note_grid, offset), [2.5,97.5], axis=0)
    new = np.percentile(fb.lr(np.array([r['validated_coefficients'] for r in accepted]), note_grid, offset), [2.5,97.5], axis=0)
    percent = 100 * (new / old - 1)
    check('upper edge effectively unchanged at 1000 and cutoff', np.max(np.abs(percent[1,:2])) < .01)
    check('upper edge narrows approximately 2%, 6%, 17% in tail', np.all(np.abs(percent[1,2:] + [2,6,17]) < .7))
    landmark = lambda b: 10 ** ((np.log(10) + offset - b[:,0]) / b[:,1])
    old_ci = np.percentile(landmark(np.array([r['original_coefficients'] for r in all_fits])), [2.5,97.5])
    new_ci = np.percentile(landmark(np.array([r['validated_coefficients'] for r in accepted])), [2.5,97.5])
    check('same-seed LR10 interval effectively stable within 1 pg/mL', np.max(np.abs(new_ci-old_ci)) < 1)
    (dest/'giorgio_note_comparison.json').write_text(json.dumps({
        'concentrations': note_grid.tolist(), 'original_band':old.tolist(), 'finite_band':new.tolist(),
        'signed_percent_difference':percent.tolist(), 'original_lr10_interval':old_ci.tolist(),
        'finite_lr10_interval':new_ci.tolist(), 'note_lr10_interval':[1913,2365],
        'qualification':'Upper-edge percentages match approximately; lower edge is not literally unchanged. The note landmark endpoints are not reproduced exactly. Existing seed-1 manuscript landmarks remain unchanged.'}, indent=2)+'\n')
    report={'scope':'bootstrap reconciliation evidence; finite-fit exclusion approved 2026-09-11',
        'count':len(checks),'failures':sum(not r['passed'] for r in checks),'checks':checks}
    (dest/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(report['count'],'checks;',report['failures'],'failures; approved finite-fit audit')
    assert report['failures']==0
if __name__=='__main__':main()
