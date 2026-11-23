#!/usr/bin/env python3
"""Read-only audit of the temperament-psychopathology genetics release repository.

Usage: python3 scripts/final_release_audit.py [--json /tmp/audit.json]
No files are modified unless --json is explicitly provided.
"""
import argparse
import csv
import hashlib
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ISSUES = []


def add(level, code, detail):
    ISSUES.append({'severity': level, 'check': code, 'detail': detail})


def git(*args):
    p = subprocess.run(['git', '-C', str(ROOT), *args], capture_output=True, text=True)
    return p.returncode, p.stdout.strip(), p.stderr.strip()


def tracked():
    p = subprocess.run(['git', '-C', str(ROOT), 'ls-files', '-z'], capture_output=True)
    if p.returncode:
        add('FAIL', 'git', 'Not a Git repository or git unavailable')
        return []
    return [x.decode('utf-8', errors='replace') for x in p.stdout.split(b'\0') if x]


def required_files(files):
    expected = [
        'README.md',
        'analysis/genomicsem/04_joint_ldsc_developmental.R',
        'analysis/genomicsem/05_test_developmental_contrasts.R',
        'analysis/genomicsem/07_age_specific_contrasts.R',
        'analysis/gwas_qc/07_prepare_age_specific.R',
        'analysis/gwas_qc/finalise_arisk.py',
        'scripts/reproduce_cross_age_temperament.py',
        'scripts/reproduce_psychiatric_gwas.py',
        'scripts/reproduce_arisk.py',
        'scripts/check_results.py',
        'config/source_script_manifest.tsv',
        'docs/CROSS_AGE_REPRODUCTION.md',
        'docs/TEMPERAMENT_DATA_PROVENANCE.md',
    ]
    for name in expected:
        if not (ROOT / name).is_file():
            add('FAIL', 'required-file', f'Missing file: {name}')
        elif name not in files:
            add('FAIL', 'untracked-required', f'Not committed to Git: {name}')


def reports(files):
    required = {
        'cross_age_temperament_20261123.json': ('results', 4),
        'age_specific_temperament_20261123.json': ('datasets', 6),
        'psychiatric_gwas_20261123.json': ('results', 5),
        'arisk_20261123.json': ('results', 4),
        'cncr_20261123.json': ('results', 2),
        'ldsc_statistical_validation_20261123.json': (None, None),
    }
    base = ROOT / 'results/reproducibility'
    for name, (field, expected_n) in required.items():
        p = base / name
        rel = str(p.relative_to(ROOT))
        if not p.exists():
            add('FAIL', 'missing-validation', rel)
            continue
        if rel not in files:
            add('FAIL', 'untracked-validation', rel)
        try:
            d = json.loads(p.read_text())
        except Exception as exc:
            add('FAIL', 'invalid-report', f'{rel}: {exc}')
            continue
        if field is not None:
            entries = d.get(field)
            if not isinstance(entries, dict) or len(entries) != expected_n:
                add('FAIL', 'report-incomplete', f'{rel}: expected {expected_n} {field} entries')
            elif not all(isinstance(v, dict) and v.get('status') == 'IDENTICAL' for v in entries.values()):
                add('FAIL', 'report-nonidentical', f'{rel}: at least one entry not IDENTICAL')
        if name.startswith('ldsc_'):
            if d.get('status') != 'PASS':
                add('FAIL', 'ldsc-validation', 'LDSC statistical validation report not PASS')
            elif d.get('age_specific', {}).get('models') != 18:
                add('FAIL', 'ldsc-validation', 'LDSC report does not identify 18 age-specific models')


def scan_tracked(files):
    home_paths = re.compile(r'(/Users/[^/\s"\x27]+/|/home/[^/\s"\x27]+/|[A-Za-z]:\\Users\\)', re.I)
    secret_patterns = [
        ('aws-key-id', re.compile(r'AKIA[0-9A-Z]{16}|ASIA[0-9A-Z]{16}')),
        ('private-key', re.compile(r'-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----')),
        ('github-token', re.compile(r'gh[pousr]_[A-Za-z0-9_]{20,}')),
    ]
    allowed_binary = {'.rds', '.gz', '.png', '.pdf', '.xlsx', '.zip', '.jpg', '.jpeg', '.webp'}
    historical = ('analysis/source_project/',)
    for rel in files:
        p = ROOT / rel
        if p.is_symlink():
            resolved = p.resolve()
            if not resolved.is_relative_to(ROOT):
                add('FAIL', 'external-symlink', f'{rel} points outside repository')
        if not p.is_file():
            continue
        size = p.stat().st_size
        if size > 50 * 1024 * 1024:
            add('FAIL', 'large-tracked-file', f'{rel}: {size / 1024**2:.1f} MiB')
        elif size > 10 * 1024 * 1024:
            add('WARN', 'large-tracked-file', f'{rel}: {size / 1024**2:.1f} MiB')
        if rel.endswith(('.sumstats.gz', '.vcf', '.vcf.gz')):
            add('FAIL', 'tracked-analysis-input', f'Dataset tracked in Git: {rel}')
        if p.suffix.lower() in allowed_binary or size > 2 * 1024 * 1024:
            continue
        try:
            content = p.read_text(encoding='utf8')
        except (UnicodeError, OSError):
            continue
        for key, pattern in secret_patterns:
            if pattern.search(content):
                add('FAIL', key, f'Potential credential in {rel}; inspect privately (no value printed)')
        path_content = content
        if rel in ('scripts/audit_repository.py', 'scripts/final_release_audit.py'):
            path_content = '\n'.join(
                line for line in content.splitlines()
                if not ('re.compile(' in line and 'Users' in line)
            )
        if home_paths.search(path_content):
            level = 'WARN' if rel.startswith(historical) else 'FAIL'
            add(level, 'local-absolute-path', f'Local account path in {rel}')
        if rel.startswith('scripts/reproduce_'):
            required_options = ('--source-root', '--analysis-root', '--workdir')
            if not all(option in content for option in required_options):
                add('WARN', 'workflow-portability', f'{rel} lacks explicit directory CLI options')
    for filename in ('README.md', 'CITATION.cff', 'LICENSE', '.gitignore'):
        if not (ROOT / filename).is_file():
            add('WARN', 'release-metadata', f'Missing {filename}')
    if (ROOT / '.gitignore').is_file():
        ignores = (ROOT / '.gitignore').read_text()
        if not any(x in ignores for x in ('*.sumstats.gz', 'data/', '*.gz')):
            add('WARN', 'gitignore', 'No obvious exclusion for large input GWAS data')


def source_manifest(files):
    p = ROOT / 'config/source_script_manifest.tsv'
    if not p.is_file():
        return
    try:
        with p.open(newline='') as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        if not rows:
            add('FAIL', 'source-manifest', 'Manifest empty')
            return
        for row in rows:
            rel = row.get('repository_path', '')
            if row.get('status') != 'COPIED':
                add('WARN', 'source-manifest', f'{row.get("script", "unknown")}: status {row.get("status")}')
                continue
            f = ROOT / rel
            if not f.is_file():
                add('FAIL', 'source-manifest', f'Missing imported source script: {rel}')
                continue
            actual = hashlib.sha256(f.read_bytes()).hexdigest()
            if row.get('sha256_source') != actual:
                add('FAIL', 'source-manifest-hash', f'Original source hash mismatch: {rel}')
            if rel not in files:
                add('FAIL', 'untracked-source', f'Imported source script not tracked: {rel}')
    except Exception as exc:
        add('FAIL', 'source-manifest', f'Unable to audit manifest: {exc}')


def git_checks(files):
    rc, out, err = git('status', '--porcelain')
    if rc:
        add('FAIL', 'git-status', err)
    elif out:
        add('FAIL', 'dirty-working-tree', 'Repository contains uncommitted changes; run git status --short')
    rc, out, _ = git('remote', '-v')
    if rc == 0 and not out:
        add('WARN', 'no-remote', 'No Git remote configured; GitHub publication not yet possible')
    rc, _, _ = git('diff', '--check', 'HEAD')
    if rc:
        add('WARN', 'git-whitespace', 'Check git diff --check for whitespace errors')
    rc, out, _ = git('branch', '--show-current')
    if rc == 0 and out:
        print('Branch:', out)
    print('Tracked files:', len(files))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json', type=Path, help='Write portable findings to this file')
    args = parser.parse_args()
    print('Repository:', ROOT)
    files = tracked()
    git_checks(files)
    required_files(files)
    source_manifest(files)
    reports(files)
    scan_tracked(files)
    print('\nFINAL REPOSITORY AUDIT')
    for severity in ('FAIL', 'WARN'):
        entries = [v for v in ISSUES if v['severity'] == severity]
        print(f'\n{severity}: {len(entries)}')
        for item in entries:
            print(f'  [{item["check"]}] {item["detail"]}')
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps({'checks': ISSUES}, indent=2) + '\n')
        print('\nSaved:', args.json)
    if not ISSUES:
        print('\nPASS: No issues detected by this static audit.')
    else:
        print('\nNOTE: Static checks cannot certify scientific assumptions or independent execution.')
    return 2 if any(v['severity'] == 'FAIL' for v in ISSUES) else 0


if __name__ == '__main__':
    raise SystemExit(main())
