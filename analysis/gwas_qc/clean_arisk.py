import gzip
from pathlib import Path

src = (
    Path.home()
    / "temperament_project/gwas/processed/munged/adult/ARISK.sumstats.gz"
)

dst = Path(
    "data/ldsc_munged/cleaned/ARISK.sumstats.gz"
)

allowed = {"A", "C", "G", "T"}
ambiguous = {
    ("A", "T"), ("T", "A"),
    ("C", "G"), ("G", "C")
}

counts = {
    "total": 0,
    "kept": 0,
    "invalid_alleles": 0,
    "strand_ambiguous": 0
}

with gzip.open(src, "rt") as fin, gzip.open(dst, "wt") as fout:
    header = fin.readline()
    fout.write(header)

    columns = header.strip().split()
    a1_idx = columns.index("A1")
    a2_idx = columns.index("A2")

    for line in fin:
        counts["total"] += 1

        parts = line.strip().split()

        if len(parts) != len(columns):
            counts["invalid_alleles"] += 1
            continue

        a1 = parts[a1_idx].upper()
        a2 = parts[a2_idx].upper()

        if a1 not in allowed or a2 not in allowed:
            counts["invalid_alleles"] += 1
            continue

        if (a1, a2) in ambiguous:
            counts["strand_ambiguous"] += 1
            continue

        fout.write(line)
        counts["kept"] += 1

print("ARISK cleaning completed")
for name, value in counts.items():
    print(f"{name}: {value:,}")

print("Output:", dst)
