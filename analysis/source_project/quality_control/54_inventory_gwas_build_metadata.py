from pathlib import Path
import gzip

root = Path.home() / "temperament_project"

gwas_dir = root / "gwas"
ref_dir = root / "references"

extensions = (
    ".tsv", ".txt", ".csv", ".gz",
    ".bgz", ".ssf", ".sumstats"
)

print("\n=== ORIGINAL GWAS FILE INVENTORY ===")

files = sorted(
    p for p in gwas_dir.rglob("*")
    if p.is_file()
    and p.name.lower().endswith(extensions)
    and "processed" not in p.relative_to(gwas_dir).parts
)

print("Original candidate files:", len(files))

for p in files[:60]:
    relative = p.relative_to(root)
    size_mb = p.stat().st_size / 1024**2

    print(f"\nFILE: {relative}")
    print(f"SIZE_MB: {size_mb:.1f}")

    try:
        opener = gzip.open if p.name.endswith((".gz", ".bgz")) else open

        with opener(p, "rt", errors="replace") as f:
            header = ""
            for line in f:
                if line.strip() and not line.startswith("##"):
                    header = line.strip()[:500]
                    break

        print("FIRST NON-METADATA LINE:", header)

    except Exception as e:
        print("HEADER UNAVAILABLE:", str(e)[:160])

if len(files) > 60:
    print(f"\nAdditional files not shown: {len(files)-60}")

print("\n=== REFERENCE GENOME / VARIANT FILES ===")

keywords = (
    "fasta", ".fa", ".fai", ".bim",
    ".pvar", ".vcf", "dbsnp",
    "grch37", "grch38", "hg19", "hg38"
)

found = []

if ref_dir.exists():
    for p in ref_dir.rglob("*"):
        if p.is_file() and any(
            k in p.name.lower() for k in keywords
        ):
            found.append(p)

for p in sorted(found)[:60]:
    print(p.relative_to(root))

print("\nReference candidates:", len(found))

print("\n=== COMPLETE ===")
