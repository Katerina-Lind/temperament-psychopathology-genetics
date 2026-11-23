import gzip
from pathlib import Path

root = Path.home() / "temperament_psychopathology"

shy = root / "data/ldsc_munged/temperament/SHY.sumstats.gz"
arisk = root / "data/ldsc_munged/cleaned/ARISK.sumstats.gz"

out = root / "data/ldsc_munged/cleaned/ARISK_SHY_compatible.sumstats.gz"
report = root / "results/ldsc/developmental_rg/SHY_ARISK_allele_mismatches.tsv"

complement = str.maketrans("ACGT", "TGCA")

def comp(allele):
    return allele.translate(complement)

def compatible(pair1, pair2):
    a, b = pair1

    valid = {
        (a, b),
        (b, a),
        (comp(a), comp(b)),
        (comp(b), comp(a))
    }

    return pair2 in valid

# Read SHY allele definitions.
shy_alleles = {}

with gzip.open(shy, "rt") as f:
    header = f.readline().split()

    snp_idx = header.index("SNP")
    a1_idx = header.index("A1")
    a2_idx = header.index("A2")

    for line in f:
        fields = line.split()

        shy_alleles[fields[snp_idx]] = (
            fields[a1_idx].upper(),
            fields[a2_idx].upper()
        )

counts = {
    "total": 0,
    "kept": 0,
    "not_in_SHY": 0,
    "incompatible": 0
}

with (
    gzip.open(arisk, "rt") as fin,
    gzip.open(out, "wt") as fout,
    open(report, "w") as freport
):
    header_line = fin.readline()
    fout.write(header_line)

    header = header_line.split()

    snp_idx = header.index("SNP")
    a1_idx = header.index("A1")
    a2_idx = header.index("A2")

    freport.write(
        "SNP\tSHY_A1\tSHY_A2\tARISK_A1\tARISK_A2\n"
    )

    for line in fin:
        counts["total"] += 1

        fields = line.split()
        snp = fields[snp_idx]

        if snp not in shy_alleles:
            counts["not_in_SHY"] += 1
            continue

        pair_shy = shy_alleles[snp]

        pair_arisk = (
            fields[a1_idx].upper(),
            fields[a2_idx].upper()
        )

        if not compatible(pair_shy, pair_arisk):
            counts["incompatible"] += 1

            freport.write(
                f"{snp}\t"
                f"{pair_shy[0]}\t{pair_shy[1]}\t"
                f"{pair_arisk[0]}\t{pair_arisk[1]}\n"
            )
            continue

        fout.write(line)
        counts["kept"] += 1

print("\nSHY–ARISK allele compatibility audit")
for key, value in counts.items():
    print(f"{key}: {value:,}")

print("\nCleaned file:", out)
print("Mismatch report:", report)
