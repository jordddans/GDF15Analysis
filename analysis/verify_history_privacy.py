"""Scan every reachable commit/tag and stored object for excluded patient exports.

Supports a normal checkout or a fresh remote mirror. Patient file contents are never
read or printed. Removed blob IDs are provenance hashes, not reconstructed records.
"""
import argparse
import json
from pathlib import Path
import re
import subprocess

EXCLUDED = [
    'archive/gdfdiscrim/data/df_all.csv',
    'archive/gdfdiscrim/gdf15_charts/gdf15_ladder_points.csv',
    'results/complete_case_patients.csv',
    'results/ntprobnp_sources.csv',
    'results/cutoff_probability_comparison.csv',
]


def verify(repo, removed_blobs=(), private_workbook_blobs=()):
    def git(*args):
        return subprocess.check_output(['git', *args], cwd=repo)
    commits = git('rev-list', '--all').decode().splitlines()
    tags = git('tag', '--list').decode().splitlines()
    refs = git('for-each-ref', '--format=%(refname)').decode().splitlines()
    bad = []
    names = {Path(p).name for p in EXCLUDED}
    for commit in commits:
        for entry in git('ls-tree', '-r', '-z', commit).split(b'\0'):
            if not entry:
                continue
            metadata, path_bytes = entry.split(b'\t', 1)
            path = path_bytes.decode()
            blob = metadata.decode().split()[2]
            if (Path(path).name in names or re.search(r'(^|/)Reconstructed_per_patient[^/]*\.xlsx$', path)
                    or path.startswith('generated/') or path in ('archive/gdfdiscrim/pvod_166.csv', 'archive/gdfdiscrim/data/df_g1.csv')
                    or blob in removed_blobs or blob in private_workbook_blobs):
                bad.append({'commit':commit, 'path':path, 'object':blob})
    tag_commits = {t:git('rev-parse', t+'^{commit}').decode().strip() for t in tags}
    stored = set(git('cat-file', '--batch-all-objects', '--batch-check=%(objectname)').decode().splitlines())
    bare = git('rev-parse','--is-bare-repository').strip() == b'true'
    checks = {
        'no_excluded_paths_or_blobs_in_any_reachable_commit': not bad,
        'every_tag_target_in_scanned_commits': all(c in commits for c in tag_commits.values()),
        'no_rewrite_backup_refs': not any(r.startswith(('refs/original/', 'refs/replace/')) for r in refs),
        'removed_patient_blobs_absent_from_object_database': not (set(removed_blobs) & stored),
        'private_workbook_blobs_absent_from_object_database': not (set(private_workbook_blobs) & stored),
    }
    if not bare:
        ignored = EXCLUDED + ['Reconstructed_per_patient_Prabhakar_2025_2.xlsx', 'generated/patients_166.csv']
        checks['gitignore_excludes_every_patient_export_and_private_workbook'] = set(
            git('check-ignore', '--no-index', *ignored).decode().splitlines()) == set(ignored)
        checks['working_tree_clean'] = not git('status', '--porcelain').strip()
    report = {'status':'PASS' if all(checks.values()) else 'FAIL', 'checks':checks,
              'commit_count':len(commits), 'commits':commits, 'tag_targets':tag_commits,
              'refs':refs, 'excluded_paths':EXCLUDED, 'forbidden_occurrences':bad,
              'removed_blob_ids':list(removed_blobs), 'private_workbook_blob_ids':list(private_workbook_blobs),
              'removed_blobs_still_stored':sorted(set(removed_blobs)&stored)}
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--removed-blobs-from', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    evidence = json.loads(args.removed_blobs_from.read_text())
    report = verify(args.repo, evidence['removed_blob_ids'], evidence.get('private_workbook_blob_ids', []))
    args.report.write_text(json.dumps(report, indent=2)+'\n')
    print(report['status'], report['commit_count'], 'commits;', len(report['tag_targets']), 'tags;',
          len(report['forbidden_occurrences']), 'forbidden occurrences')
    assert report['status'] == 'PASS', report['checks']
