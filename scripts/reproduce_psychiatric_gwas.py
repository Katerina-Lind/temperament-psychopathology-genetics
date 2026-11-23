#!/usr/bin/env python3
import csv
import gzip
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

home = Path.home()
repo = Path(__file__).resolve().parents[1]

source_root = home / "temperament_project"
analysis_root = home / "temperament_psychopathology"
workdir = home / "temperament_psychiatric_reproduction_20261123"

traits = {
    "CINT": "childhood",
    "CADHD": "childhood",
    "CAGG": "childhood",
    "AINT": "adult",
    "ASUD": "adult",
}

original_script = (
    repo / "analysis/source_project/core_preparation/"
    "33_prepare_remaining_gwas.R"
)

manifest_file = (
    source_root / "config/"
    "temperament_psychopathology_datasets.tsv"
)

hm3_file = (
    source_root / "references/ldsc/w_hm3.snplist"
)

def sha256_file(filename):
    h = hashlib.sha256()
    with gzip.open(filename, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

# Validate all source dependencies before creating anything.
with manifest_file.open() as stream:
    manifest = {
        row["trait"]: row
        for row in csv.DictReader(stream, delimiter="\t")
    }

for trait in traits:
    if trait not in manifest:
        raise SystemExit(f"Missing manifest entry: {trait}")

    source = source_root / manifest[trait]["file"]
    if not source.is_file():
        raise SystemExit(f"Missing GWAS source: {source}")

if not hm3_file.is_file():
    raise SystemExit(f"Missing HapMap3 reference: {hm3_file}")

if not original_script.is_file():
    raise SystemExit(f"Missing R script: {original_script}")

if workdir.exists():
    raise SystemExit(
        f"Workdir exists; refusing to overwrite: {workdir}"
    )

# Create isolated workspace.
(workdir / "config").mkdir(parents=True)
(workdir / "references/ldsc").mkdir(parents=True)
(workdir / "scripts").mkdir(parents=True)
(workdir / "logs").mkdir(parents=True)

(workdir / "config/temperament_psychopathology_datasets.tsv").symlink_to(
    manifest_file
)

(workdir / "references/ldsc/w_hm3.snplist").symlink_to(
    hm3_file
)

for trait in traits:
    relative = Path(manifest[trait]["file"])
    source = source_root / relative
    destination = workdir / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.symlink_to(source)

# Change only the hardcoded project-root assignment.
original_code = original_script.read_text()

old_assignment = 'root <- path.expand("~/temperament_project")'
new_assignment = f'root <- {json.dumps(str(workdir))}'

if original_code.count(old_assignment) != 1:
    raise SystemExit("Unexpected R script structure")

reproduction_code = original_code.replace(
    old_assignment, new_assignment
)

temporary_script = workdir / "scripts/33_prepare_remaining_gwas.R"
temporary_script.write_text(reproduction_code)

# Execute historical preprocessing in the isolated workspace.
log_file = workdir / "logs/psychiatric_preparation.log"

print("Running original psychiatric GWAS preparation...", flush=True)

with log_file.open("w") as log:
    result = subprocess.run(
        ["Rscript", str(temporary_script)],
        cwd=workdir,
        stdout=log,
        stderr=subprocess.STDOUT,
        check=False
    )

if result.returncode != 0:
    raise SystemExit(
        f"R preprocessing failed. Inspect: {log_file}"
    )

# Compare complete decompressed contents.
results = {}
passed = 0

for trait, domain in traits.items():
    relative = Path(domain) / f"{trait}.sumstats.gz"

    reproduced = (
        workdir / "gwas/processed/munged" / relative
    )

    expected = (
        analysis_root / "data/ldsc_munged" / relative
    )

    if not reproduced.is_file() or not expected.is_file():
        status = "MISSING"
        h1 = h2 = None
    else:
        h1 = sha256_file(reproduced)
        h2 = sha256_file(expected)
        status = "IDENTICAL" if h1 == h2 else "DIFFERENT"

    results[trait] = {
        "status": status,
        "reproduced_sha256": h1,
        "expected_sha256": h2,
    }

    print(f"{trait}: {status}", flush=True)
    if status == "IDENTICAL":
        passed += 1

report = {
    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    "historical_script": str(original_script.relative_to(repo)),
    "historical_script_sha256": hashlib.sha256(
        original_script.read_bytes()
    ).hexdigest(),
    "method": "Isolated execution; only project-root assignment changed",
    "results": results,
    "passed": passed,
    "total": len(traits),
}

report_file = workdir / "reproduction_report.json"
report_file.write_text(json.dumps(report, indent=2) + "\n")

print(f"\nRESULT: {passed}/{len(traits)} identical")
print(f"Report: {report_file}")
print(f"Log: {log_file}")

if passed != len(traits):
    raise SystemExit("Some GWAS files require investigation")

print("PSYCHIATRIC GWAS REPRODUCTION PASSED")
