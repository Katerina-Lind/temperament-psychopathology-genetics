#!/usr/bin/env python3
"""Print which expected original analysis scripts/data products are on local investigator's Mac."""
from pathlib import Path
import argparse
p = argparse.ArgumentParser()
p.add_argument('--project', type=Path, default=Path.home()/'temperament_psychopathology')
a = p.parse_args(); root=a.project.expanduser()
paths = [
  'scripts/04_genomicsem/04_joint_ldsc_developmental.R',
  'scripts/04_genomicsem/05_test_developmental_contrasts.R',
  'scripts/04_genomicsem/06_validate_rg_se.R',
  'scripts/04_genomicsem/07_age_specific_contrasts.R',
  'scripts/01_gwas_qc/07_prepare_age_specific.R',
  'results/genomicsem/developmental_7traits.rds',
  'results/genomicsem/developmental_rg_contrasts.csv',
  'results/genomicsem/age_specific/age_specific_contrasts.csv',
]
for rel in paths:
    print(('OK     ' if (root/rel).is_file() else 'MISSING'), rel)
