#!/usr/bin/env python3
"""Read-only reproducibility/security audit of imported manuscript source code.

Usage: python3 scripts/audit_repository.py [--strict]
Writes a shareable Markdown audit with file names/line numbers but never secret values.
"""
from __future__ import annotations
import argparse
import datetime as dt
import json
import re
import subprocess
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CODE_EXT = {'.py', '.R', '.r', '.sh', '.bash', '.yml', '.yaml'}
SKIP_DIRS = {'.git', '__pycache__', '.venv', 'node_modules', '.Rproj.user'}
# Patterns are indicators for manual review, not guaranteed security flaws.
RULES = [
 ('absolute_user_path','portability', re.compile(r'/' r'Users/[^/\s\"\']+|/home/[^/\s\"\']+|/mnt/data/[^\s\"\']+')),
 ('home_project_path','portability', re.compile(r'~/temperament_(?:project|psychopathology)|\$HOME/temperament_(?:project|psychopathology)|Path\.home\(\)\s*/\s*[\"\']temperament_(?:project|psychopathology)')),
 ('private_cloud_location','review', re.compile(r's3://|gs://|s3a://|globus://', re.I)),
 ('token_assignment','security', re.compile(r'(?i)\b(?:aws_secret_access_key|github_token|gh_token|api[_-]?key|access[_-]?token|secret[_-]?key|password)\b\s*(?:=|:)\s*[\"\']?[^\s,;}\)\]]+')),
 ('private_key_marker','security', re.compile(r'-----BEGIN (?:RSA |OPENSSH |EC |DSA )?PRIVATE KEY-----')),
 ('absolute_setwd','portability', re.compile(r'\bsetwd\s*\(')),
 ('os_specific_invocation','review', re.compile(r'\b(?:brew install|conda activate|module load)\b')),
 ('uncontrolled_install','reproducibility', re.compile(r'\b(?:pip install|install\.packages|remotes::install_github)\b')),
]

def collect():
    candidates = [p for folder in ('analysis','scripts','.github') for p in (ROOT/folder).rglob('*') if p.is_file() and p.suffix in CODE_EXT and p.name != 'audit_repository.py' and not any(x in p.parts for x in SKIP_DIRS)]
    findings=[]
    by_area=Counter()
    for p in sorted(candidates):
        rel=p.relative_to(ROOT).as_posix()
        try: lines=p.read_text(encoding='utf-8').splitlines()
        except (UnicodeError, OSError): continue
        by_area[rel.split('/')[0]] += 1
        for i,line in enumerate(lines,1):
            if line.lstrip().startswith('#') and not line.lstrip().startswith('#!'):
                # Documented commands are informative but not execution issues.
                continue
            for name,level,pat in RULES:
                if pat.search(line):
                    findings.append({'path':rel,'line':i,'rule':name,'category':level})
    return candidates,findings,by_area

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--strict',action='store_true',help='return nonzero if warnings exist')
    ap.add_argument('--output',default='docs/LOCAL_CODE_AUDIT.md')
    args=ap.parse_args()
    candidates, findings, areas=collect()
    out=ROOT/args.output
    out.parent.mkdir(parents=True, exist_ok=True)
    imported=sorted(p for p in (ROOT/'analysis').rglob('*') if p.is_file() and p.suffix in CODE_EXT)
    expected=['01_joint_ldsc.R','02_joint_ldsc_7traits.R','03_test_rg_contrasts.R','04_joint_ldsc_developmental.R','05_test_developmental_contrasts.R','06_validate_rg_se.R','07_age_specific_contrasts.R','07_prepare_age_specific.R','clean_arisk.py','harmonise_shy_arisk.py','run_developmental_rg.sh']
    existing={p.name for p in imported}
    missing=[n for n in expected if n not in existing]
    counts=Counter(f['category'] for f in findings)
    lines=[
      '# Local manuscript code audit', '',
      f'Generated: {dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")}', '',
      'This is a **static, read-only scan**, not an execution test. It deliberately omits matched values, which could include secrets.', '',
      f'- Imported analysis scripts found: **{len(imported)}**',
      f'- Code files scanned (including repository utilities): **{len(candidates)}**',
      f'- Findings: **{len(findings)}** (security: {counts["security"]}, portability: {counts["portability"]}, other: {counts["review"]+counts["reproducibility"]})',
      f'- Expected analysis scripts missing: **{len(missing)}**', '',
      '## Missing analysis scripts', '',
    ]
    lines.extend([f'- `{x}`' for x in missing] or ['None.'])
    lines += ['', '## Findings', '', '| Category | File | Line | Rule |', '|---|---|---:|---|']
    for x in findings:
        lines.append(f'| {x["category"]} | `{x["path"]}` | {x["line"]} | `{x["rule"]}` |')
    if not findings: lines.append('| — | — | — | None detected |')
    lines += ['', '## What to do before public release', '',
        '1. Review each flagged line locally; do not paste actual credentials into reports or chats.',
        '2. Parameterise source and output roots, then rerun each affected script on a clean working directory. Keep the original analysis revisions in Git history.',
        '3. Import missing **essential** source preparation and LDSC scripts from `~/temperament_project` into `analysis/`; record the original Git commit or SHA256 where possible.',
        '4. Verify ordering and parameterisation of GenomicSEM sampling covariance `V` before marking contrast P-values final.',
        '5. Capture `R --version`, `sessionInfo()`, the LDSC Python environment, and full source GWAS metadata; pin dependency versions.',
        '6. Check `git status --short`, `git diff --cached --stat`, and `git diff --cached --check` before pushing; never commit bulk GWAS or private data.', '',
        '## Limitations', '',
        'Regex-based detection cannot prove that a repository is free of secrets or that all research results reproduce. False positives and false negatives are possible.', '',
    ]
    out.write_text('\n'.join(lines),encoding='utf-8')
    print(f'Audit: {len(imported)} imported scripts, {len(findings)} flagged locations, {len(missing)} missing expected scripts')
    print(f'Report: {out.relative_to(ROOT)}')
    if missing: print('Missing:', ', '.join(missing))
    if findings: print('Categories:', dict(counts))
    if args.strict and (findings or missing): raise SystemExit(1)

if __name__=='__main__': main()
