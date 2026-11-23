#!/usr/bin/env python3
import gzip
import hashlib
import json
import os
import subprocess
import argparse
from pathlib import Path

home = Path.home()
repo = Path(__file__).resolve().parents[1]

parser = argparse.ArgumentParser()
parser.add_argument(
    "--source-root", type=Path,
    default=home / "temperament_project"
)
parser.add_argument(
    "--analysis-root", type=Path,
    default=home / "temperament_psychopathology"
)
parser.add_argument(
    "--workdir", type=Path,
    default=home / "temperament_arisk_reproduction_20261123"
)
args = parser.parse_args()

source_root = args.source_root.expanduser().resolve()
analysis_root = args.analysis_root.expanduser().resolve()
workdir = args.workdir.expanduser().resolve()

if (
    workdir in (source_root, analysis_root)
    or source_root in workdir.parents
    or analysis_root in workdir.parents
):
    parser.error("Workdir must be outside both source projects")


raw_rel = Path(
    "gwas/adult/risk_taking_Clifton2018/"
    "Risk_SumStats_Clifton_2018.csv.gz"
)
hm3_rel = Path("references/ldsc/w_hm3.snplist")

scripts = [
    source_root / "scripts/29_prepare_ARISK_unambiguous.R",
    source_root / "scripts/30_munge_ARISK.R",
    repo / "analysis/gwas_qc/clean_arisk.py",
    repo / "analysis/gwas_qc/finalise_arisk.py",
]

for filename in [
    source_root / raw_rel,
    source_root / hm3_rel,
    *scripts,
]:
    if not filename.is_file():
        raise SystemExit(f"Missing dependency: {filename}")

if workdir.exists():
    raise SystemExit(
        f"Refusing to use existing workspace: {workdir}"
    )

# All generated outputs stay in this fresh workspace.
(workdir / "scripts").mkdir(parents=True)
(workdir / "logs").mkdir()
(workdir / "gwas/processed/qc").mkdir(parents=True)
(workdir / "gwas/processed/munged/adult").mkdir(parents=True)
(workdir / "data/ldsc_munged/cleaned").mkdir(parents=True)

for relative in (raw_rel, hm3_rel):
    destination = workdir / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.symlink_to(source_root / relative)

def run(command, name, env=None):
    logfile = workdir / "logs" / f"{name}.log"
    print(f"Running {name}...", flush=True)

    with logfile.open("w") as out:
        result = subprocess.run(
            command,
            cwd=workdir,
            env=env,
            stdout=out,
            stderr=subprocess.STDOUT,
            check=False
        )

    if result.returncode:
        raise SystemExit(
            f"{name} failed. Inspect {logfile}"
        )

def sha256_uncompressed(filename):
    digest = hashlib.sha256()
    with gzip.open(filename, "rb") as stream:
        for chunk in iter(
            lambda: stream.read(1024 * 1024), b""
        ):
            digest.update(chunk)
    return digest.hexdigest()

# Patch only the hardcoded project-root assignment.
for original in scripts[:2]:
    code = original.read_text()
    old = 'root <- path.expand("~/temperament_project")'
    if code.count(old) != 1:
        raise SystemExit(
            f"Unexpected root assignment: {original}"
        )

    modified = code.replace(
        old, f"root <- {json.dumps(str(workdir))}"
    )

    target = workdir / "scripts" / original.name
    target.write_text(modified)

    run(["Rscript", str(target)], original.stem)

# Execute the previously validated cleaning script.
env = os.environ.copy()
env["TEMPERAMENT_SOURCE_ROOT"] = str(workdir)
env["TEMPERAMENT_ANALYSIS_ROOT"] = str(workdir)

run(
    ["python3", str(scripts[2])],
    "clean_arisk",
    env
)

# Apply the single-SNP exclusion.
cleaned = (
    workdir /
    "data/ldsc_munged/cleaned/ARISK.sumstats.gz"
)
final = (
    workdir /
    "data/ldsc_munged/cleaned/ARISK_final.sumstats.gz"
)

run(
    [
        "python3", str(scripts[3]),
        "--input", str(cleaned),
        "--output", str(final)
    ],
    "finalise_arisk"
)

# Compare every available intermediate and final output.
comparisons = {
    "intermediate": (
        workdir / "gwas/processed/qc/ARISK_unambiguous_hm3.tsv.gz",
        source_root / "gwas/processed/qc/ARISK_unambiguous_hm3.tsv.gz"
    ),
    "munged": (
        workdir / "gwas/processed/munged/adult/ARISK.sumstats.gz",
        source_root / "gwas/processed/munged/adult/ARISK.sumstats.gz"
    ),
    "cleaned": (
        cleaned,
        analysis_root / "data/ldsc_munged/cleaned/ARISK.sumstats.gz"
    ),
    "final": (
        final,
        analysis_root / "data/ldsc_munged/cleaned/ARISK_final.sumstats.gz"
    )
}

results = {}

for label, (actual, expected) in comparisons.items():
    if not actual.is_file() or not expected.is_file():
        raise SystemExit(f"Missing comparison file: {label}")

    observed_hash = sha256_uncompressed(actual)
    expected_hash = sha256_uncompressed(expected)

    status = (
        "IDENTICAL"
        if observed_hash == expected_hash
        else "DIFFERENT"
    )

    results[label] = {
        "status": status,
        "reproduced_sha256": observed_hash,
        "expected_sha256": expected_hash
    }
    print(f"{label}: {status}", flush=True)

report = {
    "workflow": "Clifton 2018 ARISK preprocessing",
    "method": "Original scripts 29 and 30 with isolated root",
    "scripts": {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in scripts
    },
    "results": results
}

(workdir / "reproduction_report.json").write_text(
    json.dumps(report, indent=2) + "\n"
)

passed = sum(
    item["status"] == "IDENTICAL"
    for item in results.values()
)

print(f"\nRESULT: {passed}/4 stages identical")

if passed != 4:
    raise SystemExit("ARISK reproduction requires investigation")

print("COMPLETE ARISK REPRODUCTION PASSED")
