# Publish the reviewed source

The source is ready for a new `williamtbarker/sigtrellis` repository. This work did not create or publish a remote repository. The available repository lookup returned 404; do not overwrite a different repository if the name is already in use or inaccessible.

The release archive includes code, installable distributions and executed evidence. `.gitignore` excludes raw downloads, generated reports, distribution files and the release checksum inventory. `release_evidence/` remains available in the download while staying out of a normal source commit. The source and original documentation are MIT; derived public reports retain their separate dataset terms.

After unpacking, run the checks in `README.md`. To create and publish the source using the authenticated GitHub CLI:

```bash
cd sigtrellis
git init -b main
git add .
git diff --cached --stat
git commit -m "Release SigTrellis 0.3.0 with native single-cell validation"
gh repo create williamtbarker/sigtrellis --public --source=. --remote=origin --push
gh run list --repo williamtbarker/sigtrellis
```

The `git diff` step shows precisely what will be published. The commit should contain no raw datasets or `release_evidence/`. GitHub Actions checks Python 3.12 and 3.13, including a separate wheel-install smoke test. Resolve remote CI failures before tagging a public release. Verify the public repository URL before adding it to package metadata or reserving a PyPI name.

The package can be used directly from source or wheel; PyPI publication is a separate action. Add the repository URL and later a DOI to `CITATION.cff` only when they exist.
