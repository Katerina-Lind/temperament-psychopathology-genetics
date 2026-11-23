#!/usr/bin/env bash
set -euo pipefail
SRC="${1:-$HOME/temperament_psychopathology}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ ! -d "$SRC" ]]; then echo "No such local analysis project: $SRC" >&2; exit 2; fi
mkdir -p "$ROOT/analysis/genomicsem" "$ROOT/analysis/gwas_qc" "$ROOT/analysis/ldsc"
copied=0
copy_one() {
  local from="$1" dest="$2"
  if [[ -f "$from" ]]; then
    cp -n "$from" "$dest/" || true
    printf 'Copied (or already exists): %s\n' "$from"
    copied=$((copied+1))
  else
    printf 'NOT FOUND: %s\n' "$from" >&2
  fi
}
for name in 01_joint_ldsc.R 02_joint_ldsc_7traits.R 03_test_rg_contrasts.R 04_joint_ldsc_developmental.R 05_test_developmental_contrasts.R 06_validate_rg_se.R 07_age_specific_contrasts.R 08_age_specific_manuscript_outputs.R; do
  copy_one "$SRC/scripts/04_genomicsem/$name" "$ROOT/analysis/genomicsem"
done
copy_one "$SRC/scripts/01_gwas_qc/07_prepare_age_specific.R" "$ROOT/analysis/gwas_qc"
copy_one "$SRC/scripts/01_gwas_qc/clean_arisk.py" "$ROOT/analysis/gwas_qc"
copy_one "$SRC/scripts/01_gwas_qc/harmonise_shy_arisk.py" "$ROOT/analysis/gwas_qc"
copy_one "$SRC/scripts/03_ldsc/run_developmental_rg.sh" "$ROOT/analysis/ldsc"
printf '\nImported or located %d scripts. REVIEW paths, source code, disclosure and provenance before git add.\n' "$copied"
printf 'Original GWAS and numerical RDS products are intentionally not copied.\n'
