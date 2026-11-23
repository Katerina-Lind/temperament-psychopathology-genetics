#!/usr/bin/env python3
"""Cross-check internal consistency of published snapshot tables; no external dependencies."""
from pathlib import Path
import csv, math
root = Path(__file__).resolve().parents[1]/'results'/'tables'
def rows(name):
    with (root/name).open(newline='') as f: return list(csv.DictReader(f))
rg=rows('Table_2_joint_LDSC_genetic_correlations.csv')
cross=rows('Table_3_cross_phenotype_contrasts.csv')
age=rows('Table_4_age_specific_contrasts.csv')
h2=rows('Table_5_age_specific_heritability.csv')
assert (len(rg),len(cross),len(age),len(h2))==(21,4,18,6)
assert len({(r['trait_1'],r['trait_2']) for r in rg}) == 21
assert len({(r['temperament'],r['outcome']) for r in age})==18
for r in cross:
  d=float(r['delta_rg']); s=float(r['SE_delta'])
  assert math.isclose(float(r['CI_low']),d-1.96*s,abs_tol=1e-8)
  assert math.isclose(float(r['CI_high']),d+1.96*s,abs_tol=1e-8)
for r in age:
  d=float(r['delta_rg_2Y_minus_3Y']); s=float(r['SE_delta'])
  assert math.isclose(d,float(r['rg_2Y'])-float(r['rg_3Y']),abs_tol=1e-5)
  assert math.isclose(float(r['CI_low']),d-1.96*s,abs_tol=1e-8)
  assert math.isclose(float(r['CI_high']),d+1.96*s,abs_tol=1e-8)
  assert float(r['min_eigenvalue_S'])>0 and float(r['min_eigenvalue_V'])>0
assert sum(float(r['FDR_q'])<.05 for r in cross)==3
assert sum(float(r['FDR_q'])<.05 for r in age)==3
print('PASS: 21 rg, 4 cross-phenotype, 18 age-specific, 6 h2; CI and FDR snapshot checks.')
print('NOTE: does not independently reproduce source GWAS or covariance parameter ordering.')
