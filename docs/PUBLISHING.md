# Public-release checklist

The code is intended for publication as a general-purpose research package for
candidate transcriptomic signature discovery. It learns a model for the supplied
study; it is not a pretrained universal cell-line biomarker classifier.

The release archive contains the source tree, wheel/source distribution, Git
history bundle, executed reports and verification evidence. Large raw datasets,
virtual environments and caches are excluded. The Git repository ignores raw
data and full generated evidence; public preparation scripts and compact
validation records are tracked. Preserve all dataset attribution/reuse notices.

The package has not been uploaded to PyPI or pushed to a remote repository by
this preparation task. Before making a public release, choose the desired GitHub
repository and visibility, inspect the tagged source and `VALIDATION.md`, and use
the supplied bundle if preserving history is desired. Do not advertise clinical
qualification, causal discovery, formal gene-level stability error control, or
an ability to extract a valid marker from every possible experimental design.

Remaining research questions include conditional inference with covariates,
better treatment of technical noise in cell variability, calibration under
cohort/prevalence shift, formal correlated-group error control, learned cell-state
representations with frozen external assignment, targeted-assay normalization,
and independent cross-laboratory transfer. These questions do not prevent the
implemented predictive discovery workflow from being useful; they delimit what
its evidence establishes.
