from pathlib import Path
import gzip
import re

root = Path.home() / "temperament_project"
gwas = root / "gwas"

sources = {
    "EMO": "temperament/raw/EASgwas_sumstats_Emotionality_Hollowell2025.txt",
    "ACT": "temperament/raw/EASgwas_sumstats_Activity_Hollowell2025.txt",
    "SHY": "temperament/raw/EASgwas_sumstats_Shyness_Hollowell2025.txt",
    "SOC": "temperament/raw/EASgwas_sumstats_Sociability_Hollowell2025.txt",
    "CINT": "childhood/internalising_GCST90054778/GCST90054778_buildGRCh37.tsv.gz",
    "CADHD": "childhood/ADHD_symptoms_GCST90568440/GCST90568440.tsv.gz",
    "CAGG": "childhood/aggression_Pappa2015_GCST002986/pappa 150810 EAGLE MA aggression_METAL RESULTS.txt",
    "AINT": "adult/internalising_PGC2025/F4_Internalizing_2025.tsv.gz",
    "ASUD": "adult/substance_use_PGC2025/F5_SubstanceUse_2025.tsv.gz",
    "ARISK": "adult/risk_taking_Clifton2018/Risk_SumStats_Clifton_2018.csv.gz",
}

outdir = root / "results/qc/gwas_provenance"
outdir.mkdir(parents=True, exist_ok=True)

report = []

def emit(*items):
    line = " ".join(map(str, items))
    print(line)
    report.append(line)

emit("=== SOURCE METADATA FILES ===")

metadata_suffixes = (
    ".yaml", ".yml", ".json", ".md",
    ".readme", ".html"
)

for p in sorted(gwas.rglob("*")):
    if not p.is_file():
        continue
    if "processed" in p.relative_to(gwas).parts:
        continue

    name = p.name.lower()

    if (
        name.startswith("readme")
        or name.startswith("meta")
        or name.endswith(metadata_suffixes)
    ):
        emit("\nFILE:", p.relative_to(root))

        if p.stat().st_size < 100_000:
            try:
                content = p.read_text(errors="replace")
                for line in content.splitlines():
                    if re.search(
                        r"build|assembl|grch|hg19|hg38|"
                        r"ancestry|sample|effect|allele|"
                        r"reference|genom|doi|ftp",
                        line,
                        flags=re.I
                    ):
                        emit("  ", line[:250])
            except Exception as exc:
                emit("  Cannot read:", str(exc)[:150])

emit("\n=== ORIGINAL GWAS VARIANT EXAMPLES ===")

for trait, relative in sources.items():
    p = gwas / relative

    emit("\nTRAIT:", trait)

    if not p.exists():
        emit("MISSING:", p)
        continue

    opener = gzip.open if p.suffix == ".gz" else open

    try:
        with opener(p, "rt", errors="replace") as f:
            lines = []

            for line in f:
                line = line.strip()
                if not line or line.startswith("##"):
                    continue

                lines.append(line[:350])

                if len(lines) >= 4:
                    break

        for i, line in enumerate(lines):
            emit(
                "HEADER:" if i == 0 else f"ROW{i}:",
                line
            )

        if trait == "CAGG":
            emit(
                "NOTE: Inspect MarkerName to determine "
                "whether it contains rsIDs or coordinates."
            )

    except Exception as exc:
        emit("READ ERROR:", str(exc)[:200])

emit("\n=== LOCAL REFERENCE FILES ===")

ref = root / "references"

if ref.exists():
    for p in sorted(ref.rglob("*")):
        if p.is_file() and (
            "hm3" in p.name.lower()
            or "ldscore" in p.name.lower()
            or "bim" in p.name.lower()
            or "snp" in p.name.lower()
        ):
            emit(p.relative_to(root))

emit("\n=== COMPLETE ===")

output = outdir / "original_source_metadata_report.txt"
output.write_text("\n".join(report) + "\n")

print("\nSaved:", output)
