# Cross-age temperament GWAS provenance

## Dataset

Hollowell cross-age temperament GWAS summary statistics:

- EMO: Emotionality
- ACT: Activity
- SHY: Shyness
- SOC: Sociability

## Historical processing

The original project contains:

- 13_prepare_emotionality.py
- 14_munge_emotionality.R
- 16_prepare_remaining_temperament.py
- 17_munge_remaining_temperament.R

These scripts prepare HapMap3-filtered input files and
run GenomicSEM::munge with SNP-specific sample sizes.

The historical pilot outputs under data/munged/
are no longer present.

## Final analytical files

Original project:
temperament_project/gwas/processed/munged/temperament/

Current analysis:
temperament_psychopathology/data/ldsc_munged/temperament/

## Verification

On 2026-11-23, the uncompressed SHA256 hashes
were compared across both projects.

| Trait | SHA256 | Result |
|---|---|---|
| EMO | a8960ccd465efd85884d7730f0160dc8f6b3a99e81f11902c422f06f225a42c2 | Identical |
| ACT | 88511e90222fa87e95e3ae41856b54c7540bb81adaa17409b82331b251a27c79 | Identical |
| SHY | 16ffd868213ddcb578c20b98687f86763f347e55d26d7b417ec661a92180cc13 | Identical |
| SOC | d9a0fb009580376944165f609e7395e45179b01ac169f0dcc48621a7930137d8 | Identical |

## Interpretation

The final temperament summary statistics were transferred
between analysis projects without changes to their
uncompressed contents.

This does not constitute an independent reproduction
of the original GWAS preprocessing.

The original preparation workflow must still be tested
against the final analytical inputs.
