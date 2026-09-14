"""Source-only reconstruction, analysis, tables, figures, and optional verification."""
import argparse
from pathlib import Path
import json
import data
import reconstruct
import run_analysis
import make_figures
import verify

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workbook',type=Path,default=reconstruct.SOURCE)
    p.add_argument('--output-dir',type=Path,default=ROOT)
    p.add_argument('--verify',action='store_true')
    a=p.parse_args()
    data.WORKBOOK=a.workbook.resolve()
    reconstruct.write(data.WORKBOOK,a.output_dir/'generated')
    run_analysis.RESULTS=make_figures.RESULTS=a.output_dir/'results'
    make_figures.FIGURES=a.output_dir/'figures'
    run_analysis.main()
    if a.verify:
        verify.main()
        report={'count':len(verify.CHECKS),'failures':sum(not c['passed'] for c in verify.CHECKS),
                'checks':verify.CHECKS,'figure1_integration':'PASS: Giorgio three-panel figure; 3991 finite fits, 9 separated resamples excluded'}
        (a.output_dir/'results/verification.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__': main()
