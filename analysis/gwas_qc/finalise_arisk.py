#!/usr/bin/env python3

import argparse
import gzip
import hashlib
from pathlib import Path

EXCLUDED_SNPS = {"rs10038745"}

EXPECTED_SHA256 = (
    "54633122986bc096734e48275f890"
    "f928fe1e76ca8da7bd438cc43aae7543bb1"
)


def decompressed_sha256(filename):
    digest = hashlib.sha256()

    with gzip.open(filename, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    source = Path(args.input)
    destination = Path(args.output)

    if not source.is_file():
        raise SystemExit(f"Missing input: {source}")

    if destination.exists():
        raise SystemExit(f"Refusing to overwrite: {destination}")

    if source.resolve() == destination.resolve():
        raise SystemExit("Input and output must differ")

    destination.parent.mkdir(parents=True, exist_ok=True)

    total = 0
    removed = 0

    with gzip.open(source, "rt") as fin, gzip.open(destination, "wt") as fout:
        header_line = fin.readline()
        columns = header_line.strip().split()
        snp_index = columns.index("SNP")

        fout.write(header_line)

        for line in fin:
            fields = line.split()
            total += 1

            if fields[snp_index] in EXCLUDED_SNPS:
                removed += 1
                continue

            fout.write(line)

    if removed != 1:
        raise SystemExit(
            f"Expected one excluded record; observed {removed}"
        )

    observed_sha256 = decompressed_sha256(destination)

    print("Input records:", total)
    print("Excluded records:", removed)
    print("Output records:", total - removed)
    print("SHA256:", observed_sha256)

    if observed_sha256 != EXPECTED_SHA256:
        raise SystemExit("Output differs from the validated final ARISK")

    print("ARISK FINAL REPRODUCTION PASSED")


if __name__ == "__main__":
    main()
