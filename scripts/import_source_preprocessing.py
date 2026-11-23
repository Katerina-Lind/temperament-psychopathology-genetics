#!/usr/bin/env python3
"""Non-destructive, allowlisted import of original source GWAS processing scripts.

Usage: python3 scripts/import_source_preprocessing.py --source ~/temperament_project
Run at repository root. Never copies raw GWAS or derived patient-level data.
"""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

ALLOWLIST = {
    "core_preparation": [
        "03_munge.R", "30_munge_ARISK.R", "33_prepare_remaining_gwas.R",
        "40_make_strict_allele_inputs.R",
    ],
    "quality_control": [
        "26_audit_psychopathology_gwas.R", "27_validate_gwas_statistics.R",
        "34_check_ldsc_coverage.R", "39_audit_cross_gwas_alleles.R",
        "47_audit_sample_size_ldscore.R", "52_gwas_provenance_audit.R",
        "53_annotate_gwas_provenance.R", "54_inventory_gwas_build_metadata.py",
        "55_inspect_original_gwas_metadata.py", "56_validate_gwas_effect_directions.R",
    ],
    "ldsc_reference": [
        "37_ldsc_childhood.R", "38_ldsc_temperament_childhood.R",
        "41_ldsc_full_ten_traits.R", "42_test_developmental_contrasts.R",
        "49_ldsc_full_individual_highN.R", "50_validate_rg_gradients.R",
        "51_finalise_rg_results.R",
    ],
}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True, help="Original temperament_project directory")
    parser.add_argument("--repository", type=Path, default=Path.cwd())
    args = parser.parse_args()
    source = args.source.expanduser().resolve()
    repo = args.repository.expanduser().resolve()
    if not (source / "scripts").is_dir():
        parser.error(f"Missing source scripts directory: {source / 'scripts'}")
    if not (repo / "analysis").is_dir():
        parser.error(f"Missing repository analysis directory: {repo / 'analysis'}")

    dest_base = repo / "analysis" / "source_project"
    dest_base.mkdir(parents=True, exist_ok=True)
    rows = []
    for group, filenames in ALLOWLIST.items():
        target_dir = dest_base / group
        target_dir.mkdir(exist_ok=True)
        for filename in filenames:
            src = source / "scripts" / filename
            dst = target_dir / filename
            if not src.is_file():
                status, digest = "NOT_FOUND", ""
            elif dst.exists():
                digest = sha256(src)
                status = "EXISTS_IDENTICAL" if sha256(dst) == digest else "EXISTS_DIFFERENT_REVIEW"
            else:
                data = src.read_bytes()
                dst.write_bytes(data)
                digest = hashlib.sha256(data).hexdigest()
                status = "COPIED"
            rows.append({"group": group, "script": filename, "status": status,
                         "sha256_source": digest, "repository_path": str(dst.relative_to(repo))})
            print(f"{status:25} {group}/{filename}")

    cfg = repo / "config"
    cfg.mkdir(exist_ok=True)
    manifest = cfg / "source_script_manifest.tsv"
    with manifest.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nManifest: {manifest.relative_to(repo)}")
    print("This command did not execute the scripts, copy data, or modify previously imported code.")
    print("Review all copied scripts before staging files with git.")

if __name__ == "__main__":
    main()
