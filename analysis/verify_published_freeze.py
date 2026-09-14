"""Final regression, QC-document scope, privacy and loader-boundary checks."""
from pathlib import Path
from zipfile import ZipFile
import hashlib
import json
import re
import subprocess
import sys
import inspect
import io
import check_regression
import data
import model
from make_qc_manuscript import plain
from figure1_bootstrap import METHODS

ROOT = Path(__file__).resolve().parents[1]


def main():
    checks = []
    def check(name, condition):
        checks.append({'name': name, 'passed': bool(condition)})
        assert condition, name
    regression = check_regression.check()
    (ROOT / 'qc/baseline_regression.json').write_text(json.dumps(regression, indent=2) + '\n')
    check('685 baseline numerical/Boolean fields unchanged', regression['numeric_fields_compared'] == 685 and not regression['numeric_differences'])
    check('baseline scientific CSV outputs and Figures 2/3 unchanged', all(r['science_unchanged'] for r in regression['baseline_artifacts']))
    check('Table 1 CSV and Markdown byte-identical to prior verified QC', all(
          (ROOT / path).read_bytes() == subprocess.check_output(['git','show','ac8339e:' + path], cwd=ROOT)
          for path in ['results/table1.csv','results/table1.md']))
    old_caption = subprocess.check_output(['git','show','ac8339e:figures/figure2_caption.md'], cwd=ROOT, text=True)
    expected_caption = old_caption.replace(') as Figure 1. Each', '). Each').replace('the main model in Figure 1:', 'the main-model landmark reported in the analysis:')
    check('Figure 2 caption changes only the outdated Figure 1 bootstrap cross-reference', (ROOT/'figures/figure2_caption.md').read_text() == expected_caption)
    check('original DOCX, workbooks and protected PDFs unchanged', all(r['unchanged'] for r in regression['protected_files']))
    check('private reconstructed workbooks absent from all Git history', regression['private_exports_absent_from_git_history'])
    before = subprocess.check_output(['git', 'show', 'ac8339e:GDF15_evidence_curve_QC.docx'], cwd=ROOT)
    docx = ROOT / 'GDF15_evidence_curve_QC.docx'
    with ZipFile(io.BytesIO(before)) as z:
        original = {n: z.read(n) for n in z.namelist()}
    with ZipFile(docx) as z:
        revised = {n: z.read(n) for n in z.namelist()}
    blocks = lambda x: [plain(b) for b in re.findall(r'<w:tbl(?:\s[^>]*)?>.*?</w:tbl>|<w:p(?:\s[^>]*)?>.*?</w:p>', x.decode(), re.S)]
    old, new = blocks(original['word/document.xml']), blocks(revised['word/document.xml'])
    changed = [i for i,(a,b) in enumerate(zip(old,new)) if a != b]
    check('QC manuscript changes confined to status, bootstrap Methods and Figure 1 caption', len(old) == len(new) and changed == [2,33,98])
    check('approved Methods wording present exactly', METHODS in new[33])
    image = [n for n in revised if n.startswith('word/media/') and n.endswith('.png')]
    check('QC DOCX embeds exactly the reproduced Figure 1', len(image) == 1 and revised[image[0]] == (ROOT/'figures/figure1_lr_curve.png').read_bytes())
    check('other DOCX package entries unchanged', original.keys() == revised.keys() and all(original[n] == revised[n] for n in original if n not in ['word/document.xml', *image]))
    digest = hashlib.sha256(docx.read_bytes()).hexdigest()
    subprocess.run([sys.executable, str(ROOT/'analysis/make_qc_manuscript.py')], cwd=ROOT, check=True)
    check('QC DOCX repeated generation byte-identical', hashlib.sha256(docx.read_bytes()).hexdigest() == digest)
    tracked = subprocess.check_output(['git','ls-files'], cwd=ROOT, text=True).splitlines()
    excluded = ['Reconstructed_per_patient_Prabhakar_2025_2.xlsx', 'archive/gdfdiscrim/data/df_all.csv',
                'archive/gdfdiscrim/gdf15_charts/gdf15_ladder_points.csv',
                'results/ntprobnp_sources.csv', 'results/complete_case_patients.csv', 'results/cutoff_probability_comparison.csv']
    check('private workbook and patient-level exports untracked', all(p not in tracked for p in excluded) and not any(p.startswith('generated/') for p in tracked))
    ignored = subprocess.check_output(['git','check-ignore',*excluded,'generated/patients_166.csv'], cwd=ROOT, text=True).splitlines()
    check('private workbook and generated patient exports ignored', set(ignored) == set(excluded + ['generated/patients_166.csv']))
    original_source = ROOT/'archive/gdfdiscrim/jci-136-199159-s068.xlsx'
    cleaned = ROOT/'Supporting data values_Prabhakar et al 2025_CLEANED.xlsx'
    check('published and cleaned workbook readers agree', data.load_gdf15(original_source).equals(data.load_gdf15(cleaned)))
    check('loader input is explicit and fitting independent of data loading',
          'workbook' in inspect.signature(data.load_gdf15).parameters and
          list(inspect.signature(model.fit_lr_model).parameters) == ['gdf15','is_case','name'] and
          not re.search(r'^import data|^from data ', inspect.getsource(model), re.M))
    check('scoring uses frozen coefficients without fitting', 'fit(' not in inspect.getsource(model.LRModel.probability))
    clean = json.loads((ROOT/'qc/published_data_clean_run.json').read_text())
    check('clean public-input reproduction and repeat artifact comparisons pass', clean['failures'] == 0 and clean['no_private_workbook_or_cached_results_in_inputs'] and clean['repeated_artifact_count'] > 30)
    audit = json.loads((ROOT/'qc/figure1_bootstrap/verification.json').read_text())
    check('independent bootstrap audit passes', audit['failures'] == 0)
    report = {'status':'PASS', 'count':len(checks), 'failures':0, 'checks':checks,
              'workflow_checks':clean['workflow_checks'], 'bootstrap_audit_checks':audit['count'],
              'total_checks':len(checks) + clean['workflow_checks'] + audit['count'],
              'numeric_fields_compared_separately':685, 'qc_docx_sha256':digest,
              'qc_changed_body_blocks':changed,
              'history_note':'Private reconstruction absent throughout. Five published-source patient-level exports were tracked in earlier local commits; removed from frozen tree, retained locally and ignored. History was not rewritten.'}
    (ROOT/'qc/published_data_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    document = {'qc_docx':docx.name, 'docx_sha256':digest,
                'logged_mechanical_edits':len(json.loads((ROOT/'qc/manuscript_changes.json').read_text())),
                'docx_repeat_byte_identical':True, 'original_docx_unchanged':True,
                'figure1_image':'integrated production PNG; 3991 finite fits, 9 separated exclusions',
                'changed_body_blocks_from_previous_qc':changed,
                'pdf':'qc/manuscript_QC.pdf',
                'pdf_note':'Local Microsoft Word rendering is a review derivative; analysis figure PDF/PNG and QC DOCX byte determinism are checked separately.'}
    (ROOT/'qc/document_verification.json').write_text(json.dumps(document, indent=2)+'\n')
    print(report['total_checks'], 'total checks; zero failures')


if __name__ == '__main__':
    main()
