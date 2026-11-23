#!/usr/bin/env python3
"""Reproduce historical Hollowell cross-age munging in an isolated workspace.

Default: preflight only. --run executes the original 13/14/16/17 scripts, without
modifying the original GWAS inputs or already validated munged summary statistics.
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime, timezone

TRAITS = ('EMO', 'ACT', 'SHY', 'SOC')
SCRIPT_NAMES = (
    '13_prepare_emotionality.py',
    '14_munge_emotionality.R',
    '16_prepare_remaining_temperament.py',
    '17_munge_remaining_temperament.R',
)
EXPECTED = {
    'EMO': 'a8960ccd465efd85884d7730f0160dc8f6b3a99e81f11902c422f06f225a42c2',
    'ACT': '88511e90222fa87e95e3ae41856b54c7540bb81adaa17409b82331b251a27c79',
    'SHY': '16ffd868213ddcb578c20b98687f86763f347e55d26d7b417ec661a92180cc13',
    'SOC': 'd9a0fb009580376944165f609e7395e45179b01ac169f0dcc48621a7930137d8',
}


def sha256_decompressed(path):
    h = hashlib.sha256()
    with gzip.open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def require_file(path, failures):
    if not path.is_file():
        failures.append(f'Missing file: {path}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, default=Path.home() / 'temperament_project')
    parser.add_argument('--analysis-root', type=Path, default=Path.home() / 'temperament_psychopathology')
    parser.add_argument('--workdir', type=Path, default=Path.home() / 'temperament_cross_age_reproduction')
    parser.add_argument('--run', action='store_true', help='Execute historical preprocessing; otherwise preflight only')
    args = parser.parse_args()
    root = args.source_root.expanduser().resolve()
    analysis = args.analysis_root.expanduser().resolve()
    workdir = args.workdir.expanduser().resolve()
    if workdir == root or workdir == analysis or root in workdir.parents or analysis in workdir.parents:
        parser.error('--workdir must be outside both original projects')

    failures = []
    raw = root / 'data/raw/temperament'
    hm3 = root / 'references/ldsc/w_hm3.snplist'
    for t, name in (('EMO','Emotionality'),('ACT','Activity'),('SHY','Shyness'),('SOC','Sociability')):
        require_file(raw / f'EASgwas_sumstats_{name}_Hollowell2025.txt', failures)
        require_file(root / 'gwas/processed/munged/temperament' / f'{t}.sumstats.gz', failures)
        require_file(analysis / 'data/ldsc_munged/temperament' / f'{t}.sumstats.gz', failures)
    require_file(hm3, failures)
    for name in SCRIPT_NAMES:
        require_file(root / 'scripts' / name, failures)
    for tool in ('python3','Rscript'):
        if not shutil.which(tool):
            failures.append(f'Executable not found on PATH: {tool}')
    if failures:
        for issue in failures:
            print('ERROR:', issue, file=sys.stderr)
        return 2

    print('Preflight: all historical scripts, original GWAS, reference and validated files exist.')
    verified = True
    for trait in TRAITS:
        original = root / 'gwas/processed/munged/temperament' / f'{trait}.sumstats.gz'
        current = analysis / 'data/ldsc_munged/temperament' / f'{trait}.sumstats.gz'
        o_hash, c_hash = sha256_decompressed(original), sha256_decompressed(current)
        ok = o_hash == c_hash == EXPECTED[trait]
        print(f'{trait}: original={o_hash} current={c_hash} verified={ok}')
        verified = verified and ok
    if not verified:
        print('ERROR: baseline file content differs from recorded verified hashes.', file=sys.stderr)
        return 3
    if not args.run:
        print('DRY RUN ONLY. To regenerate in an isolated directory, add --run.')
        print('Workdir (not created in preflight):', workdir)
        return 0

    if workdir.exists() and any(workdir.iterdir()):
        print('ERROR: workdir already exists and is nonempty; choose a fresh directory.', file=sys.stderr)
        return 4
    (workdir / 'scripts').mkdir(parents=True)
    (workdir / 'data/raw').mkdir(parents=True)
    (workdir / 'data/processed/temperament').mkdir(parents=True)
    (workdir / 'data/munged').mkdir(parents=True)
    (workdir / 'references/ldsc').mkdir(parents=True)
    (workdir / 'logs').mkdir(parents=True)
    (workdir / 'data/raw/temperament').symlink_to(raw, target_is_directory=True)
    (workdir / 'references/ldsc/w_hm3.snplist').symlink_to(hm3)
    for name in SCRIPT_NAMES:
        shutil.copy2(root / 'scripts' / name, workdir / 'scripts' / name)

    manifest = {'timestamp_utc': datetime.now(timezone.utc).isoformat(),
                'source_root': str(root), 'analysis_root': str(analysis),
                'workdir': str(workdir),
                'scripts_sha256': {name: hashlib.sha256((workdir/'scripts'/name).read_bytes()).hexdigest() for name in SCRIPT_NAMES},
                'results': {}}
    commands = [
        ['python3', 'scripts/13_prepare_emotionality.py'],
        ['Rscript', 'scripts/14_munge_emotionality.R'],
        ['python3', 'scripts/16_prepare_remaining_temperament.py'],
        ['Rscript', 'scripts/17_munge_remaining_temperament.R'],
    ]
    for n, command in enumerate(commands, 1):
        log = workdir / 'logs' / f'{n:02d}_{Path(command[-1]).stem}.log'
        print('Executing:', ' '.join(command), '| log:', log, flush=True)
        with log.open('w') as out:
            ret = subprocess.run(command, cwd=workdir, stdout=out, stderr=subprocess.STDOUT, check=False).returncode
        if ret:
            print(f'ERROR: command exited {ret}; inspect {log}', file=sys.stderr)
            manifest['failed_command'] = command
            manifest['exit_code'] = ret
            (workdir / 'reproduction_report.json').write_text(json.dumps(manifest, indent=2) + '\n')
            return 5

    equal_all = True
    for trait in TRAITS:
        generated = workdir / 'data/munged' / f'{trait}.sumstats.gz'
        if not generated.exists():
            print(f'ERROR: no generated file for {trait} at {generated}', file=sys.stderr)
            manifest['results'][trait] = {'status': 'MISSING'}
            equal_all = False
            continue
        digest = sha256_decompressed(generated)
        equal = digest == EXPECTED[trait]
        manifest['results'][trait] = {'status': 'IDENTICAL' if equal else 'DIFFERENT', 'decompressed_sha256': digest, 'expected_sha256': EXPECTED[trait]}
        print(f'{trait}: {manifest["results"][trait]["status"]}, decompressed SHA256={digest}')
        equal_all = equal_all and equal
    (workdir / 'reproduction_report.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print('Report:', workdir / 'reproduction_report.json')
    if not equal_all:
        print('WARNING: differences require investigation; this does not automatically imply scientific disagreement.')
        return 6
    print('SUCCESS: all regenerated cross-age files are text-identical to the verified final analytical inputs.')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
