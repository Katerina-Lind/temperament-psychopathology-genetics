# Publish after code/provenance review

1. Create a new empty public GitHub repository under your own account or lab organisation. Suggested name: `early-temperament-psychiatric-genetics`. Do not initialise it with a README because this archive already has one.
2. Extract this archive into a folder, open a terminal there, and import the original analysis scripts on your Mac.
3. Verify and redact imported scripts, source hashes, citation and licence. Run `python3 scripts/check_results.py`.
4. Create the remote repo with GitHub CLI if authenticated:

```bash
git init
git add .
git commit -m "Prepare manuscript reproducibility repository"
gh auth status
gh repo create early-temperament-psychiatric-genetics --public --source=. --remote=origin --push
```

If the remote repository already exists, use `git remote add origin <YOUR-GITHUB-REPOSITORY-URL>` and `git push -u origin main` after checking branch name. Never paste tokens into scripts or chat.

5. When the analysis is final, tag a version and archive the release via Zenodo. Use the immutable archived DOI in Supplementary Materials, not just an unversioned `main` branch URL.

**Do not present the current provisional snapshot release as sufficient for end-to-end reproduction.**
