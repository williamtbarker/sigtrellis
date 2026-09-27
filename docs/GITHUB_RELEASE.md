# Maintaining GitHub releases

The public source repository is [williamtbarker/sigtrellis](https://github.com/williamtbarker/sigtrellis). Clone it normally to preserve its history:

```bash
git clone https://github.com/williamtbarker/sigtrellis.git
cd sigtrellis
```

Before preparing a release, run the checks in [CONTRIBUTING.md](../CONTRIBUTING.md), build the wheel/source distribution, and test the installed wheel outside the source tree. Inspect [VALIDATION.md](../VALIDATION.md) and retain both successful and failed scientific gates. Verify the GitHub Actions results for the exact commit being released; local checks do not establish hosted CI success.

Review the diff and package contents before committing or tagging. The source repository includes compact validation evidence and excludes raw datasets, virtual environments, and full generated reports. A separately prepared evidence archive may include those reports subject to the dataset reuse terms; check its inventory rather than assuming they are present in a source checkout.

Keep release tags attached to the reviewed source commit. Check the tag author and target, package version, changelog, citation metadata, distributions, and release notes for consistency. Do not move an existing release tag to a newer source revision as part of ordinary documentation maintenance.

The package can be installed from source or a built wheel. PyPI publication and an archived DOI are separate publication steps; do not advertise either unless it exists.
