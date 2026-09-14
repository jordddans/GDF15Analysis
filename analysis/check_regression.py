"""Verify baseline scientific invariants and protected-file hashes, without git mutation."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
BASELINE='c4a5a20dd2024f6a35bb0b898829cdd18e347936'

def check(output=ROOT):
    output=Path(output)
    manifest=json.loads((ROOT/'qc/baseline_manifest.json').read_text())['files']
    old=json.loads(subprocess.check_output(['git','show',BASELINE+':results/results.json'],cwd=ROOT))
    new=json.loads((output/'results/results.json').read_text())
    numeric_count=0;differences=[];text_changes=[]
    def walk(a,b,path):
        nonlocal numeric_count
        if isinstance(a,dict):
            assert a.keys()==b.keys(),path
            for k in a:walk(a[k],b[k],path+'/'+k)
        elif isinstance(a,list):
            assert len(a)==len(b),path
            for i,(x,y) in enumerate(zip(a,b)):walk(x,y,path+f'/{i}')
        elif isinstance(a,(int,float)):
            numeric_count+=1
            if a!=b:differences.append(path)
        elif a!=b:text_changes.append({'path':path,'before':a,'after':b})
    walk(old,new,'results.json')
    protected=[];output_files=[]
    for entry in manifest:
        path=entry['path'];p=ROOT/path
        if path.endswith('.xlsx') or path in ['BayLearn_Submission.pdf','archive/gdfdiscrim/BayLearn_Submission.pdf','GDF15_evidence_curve_DRAFT_v0.1_2-1.docx']:
            if p.exists():protected.append({'path':path,'unchanged':hashlib.sha256(p.read_bytes()).hexdigest()==entry['sha256']})
        if path.startswith('results/') and path.endswith('.csv'):
            p=output/path
            if path=='results/auc_reconciliation.csv':
                import io
                a=pd.read_csv(io.BytesIO(subprocess.check_output(['git','show',BASELINE+':'+path],cwd=ROOT)))
                b=pd.read_csv(p)
                equal=a.drop(columns='published_provenance').equals(b.drop(columns='published_provenance'))
            else:equal=hashlib.sha256(p.read_bytes()).hexdigest()==entry['sha256']
            output_files.append({'path':path,'science_unchanged':equal})
        if path.startswith('figures/') and not path.startswith('figures/figure1') and path != 'figures/figure2_caption.md':
            output_files.append({'path':path,'science_unchanged':hashlib.sha256((output/path).read_bytes()).hexdigest()==entry['sha256']})
    forbidden=['Reconstructed_per_patient_Prabhakar_2025_2.xlsx','archive/gdfdiscrim/pvod_166.csv','archive/gdfdiscrim/data/df_g1.csv']
    history=subprocess.check_output(['git','log','--all','--format=','--name-only'],cwd=ROOT,text=True).splitlines()
    tracked=subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines()
    private_absent=all(p not in history and p not in tracked for p in forbidden)
    r={'numeric_fields_compared':numeric_count,'numeric_differences':differences,'text_changes':text_changes,
       'protected_files':protected,'baseline_artifacts':output_files,'private_exports_absent_from_git_history':private_absent,
       'status':'PASS' if not differences and all(x['unchanged'] for x in protected) and all(x['science_unchanged'] for x in output_files) and private_absent else 'FAIL'}
    assert r['status']=='PASS',r
    return r
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,default=ROOT);p.add_argument('--report',type=Path,default=ROOT/'qc/baseline_regression.json');a=p.parse_args()
    result=check(a.output_dir);a.report.write_text(json.dumps(result,indent=2)+'\n');print(result['status'],result['numeric_fields_compared'],'numeric fields')
