# Data availability and source verification

The manuscript uses publicly available or externally distributed summary-level GWAS data. Bulk data are not part of this Git repository. Dataset-specific accession/identifiers are in `config/data_sources.tsv`.

**Hollowell EAS GWAS:** European-ancestry, cross-age and age-specific files, with source metadata and checksum information from the investigator's `hollowell_figshare.json`. Download using the official Figshare links. Retain Figshare article metadata, DOI, file ID, download date and hash in the provenance record. Relevant sample sizes are *reported study/meta-analysis counts* and may not be unique individuals.

**Childhood traits:** GCST90568440 (ADHD), GCST90054778 (internalising); verify differing effective and unique sample sizes. 

**Adult traits:** PGC 2025 F4 and F5 GWAS; Clifton 2018 risk-taking. Worry and depressed-affect GWAS require completed DOI/source citation entry before release.

**Reference data:** European-ancestry LDSC LD reference (`eur_w_ld_chr`) and HapMap3 merge SNP list (`w_hm3.snplist`) acquired from legitimate source. Do not redistribute in this repository without checking licences.

**Required checks before release:** original source publication references/DOIs, source file hashes, GWAS allele/variant build provenance, exact software versions, and confirmation that open GWAS licences permit derived summary results and citations.
