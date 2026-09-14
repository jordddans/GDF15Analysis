"""Isolate code plus the published workbook; reproduce twice with a supplied fresh Python.

Usage: python analysis/freeze_reproduce.py --python /path/to/fresh/venv/bin/python
Install requirements-lock.txt into that new environment first. No private inputs or
existing result artifacts are copied into the isolated tree.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python', required=True, type=Path)
    args = parser.parse_args()
    scratch = Path(tempfile.mkdtemp(prefix='gdf15-published-reproduce-'))
    code = scratch / 'project'
    (code / 'analysis').mkdir(parents=True)
    for path in (ROOT / 'analysis').glob('*.py'):
        shutil.copy2(path, code / 'analysis' / path.name)
    source = Path('archive/gdfdiscrim/jci-136-199159-s068.xlsx')
    (code / source).parent.mkdir(parents=True)
    shutil.copy2(ROOT / source, code / source)
    shutil.copy2(ROOT / 'requirements-lock.txt', code / 'requirements-lock.txt')
    inputs = sorted(str(p.relative_to(code)) for p in code.rglob('*') if p.is_file())
    env = __import__('os').environ.copy()
    env['MPLCONFIGDIR'] = str(scratch / 'matplotlib')
    for name in ('run1', 'run2'):
        with (scratch / f'{name}.log').open('w') as log:
            command = [str(args.python), str(code / 'analysis/workflow.py'),
                       '--output-dir', str(scratch / name)]
            if name == 'run1':
                command.append('--verify')
            subprocess.run(command, cwd=code, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    def hashes(path):
        return {str(p.relative_to(path)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(path.rglob('*')) if p.is_file() and p.name != 'verification.json'}
    first, second = hashes(scratch / 'run1'), hashes(scratch / 'run2')
    assert first == second, 'Clean repeated outputs differ'
    verification = json.loads((scratch / 'run1/results/verification.json').read_text())
    assert verification['failures'] == 0
    # Only after success, promote the reproduced results/figures and local reconstruction.
    for directory in ('results', 'figures', 'generated'):
        shutil.copytree(scratch / 'run1' / directory, ROOT / directory, dirs_exist_ok=True)
    report = {'status': 'PASS', 'python': subprocess.check_output([str(args.python), '--version'], text=True).strip(),
              'fresh_environment_python': str(args.python), 'isolated_directory': str(scratch),
              'input_files': inputs, 'public_workbook_sha256': hashlib.sha256((code / source).read_bytes()).hexdigest(),
              'workflow_checks': verification['count'], 'failures': 0,
              'repeated_artifact_count': len(first), 'repeated_artifact_sha256': first,
              'no_private_workbook_or_cached_results_in_inputs': True}
    (ROOT / 'qc/published_data_clean_run.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('input_files', 'repeated_artifact_sha256')}, indent=2))


if __name__ == '__main__':
    main()
