# Checklist before publishing and citing the GitHub repository

- [ ] Import author-written analysis scripts from both local projects; inspect and preserve actual final versions.
- [ ] Remove user-specific absolute paths, usernames, embedded tokens, credentials, copyrighted source attachments or unlicensed GWAS data.
- [ ] Regenerate figures and tables directly from exact numerical analysis outputs; avoid rounded transcriptions.
- [ ] Confirm GenomicSEM sampling covariance ordering and delta-method contrast standard errors against package implementation.
- [ ] Confirm GWAS data sources, sample-size conventions, participant overlap, ancestry and reference genome builds.
- [ ] Record tested hypotheses and the scope of multiplicity correction accurately.
- [ ] Add complete authors, affiliations, acknowledgments, a final repository title, manuscript DOI, `CITATION.cff` metadata and an agreed software licence.
- [ ] Pin software versions in a reproducible environment (`renv.lock` or equivalent) after verifying the original Mac installations.
- [ ] Choose an immutable public GitHub release tag (e.g. `v1.0.0`) and archive it in Zenodo to obtain a version-specific DOI.
- [ ] Include the repository release URL and DOI in manuscript Code Availability and Supplementary Methods.
