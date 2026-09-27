# Public-release checklist

The code is intended for publication as a general-purpose research package for
candidate transcriptomic signature discovery. It learns a model for the supplied
study; it is not a pretrained universal cell-line biomarker classifier.

The release archive contains the source tree, wheel/source distribution, Git
history bundle, executed reports and verification evidence. Large raw datasets,
virtual environments and caches are excluded. The Git repository ignores raw
data and full generated evidence; public preparation scripts and compact
validation records are tracked. Preserve all dataset attribution/reuse notices.

The public source is maintained at
[williamtbarker/sigtrellis](https://github.com/williamtbarker/sigtrellis).
Before preparing a release, inspect the exact source commit, package contents,
GitHub Actions results, and `VALIDATION.md`. See [GITHUB_RELEASE.md](GITHUB_RELEASE.md)
for the maintenance workflow. Source availability does not establish PyPI publication. Do not advertise clinical
qualification, causal discovery, formal gene-level stability error control, or
an ability to extract a valid marker from every possible experimental design.

Remaining research questions include conditional inference with covariates,
better treatment of technical noise in cell variability, calibration under
cohort/prevalence shift, formal correlated-group error control, learned cell-state
representations with frozen external assignment, targeted-assay normalization,
and independent cross-laboratory transfer. These questions do not prevent the
implemented predictive discovery workflow from being useful; they delimit what
its evidence establishes.
