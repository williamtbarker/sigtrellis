# Validation evidence archive

These files preserve historical software checks and scientific-control results.
Use the versioned directories and [current validation report](../../VALIDATION.md)
to distinguish release evidence; an earlier run is not evidence for later code.
Failed baseline regressions and warnings are retained alongside successful checks.

Machine-specific absolute path prefixes in the archived text logs have been
replaced with `/path/to/python-runtime` and `/path/to/validation-workspace`.
These are descriptive placeholders, not runnable installation locations. Relative
paths, package versions, source and wheel hashes, test outcomes, warning text,
coverage, timings, and scientific results are unchanged. This normalization does
not rerun or upgrade any historical validation claim.

The source hashes in adjacent JSON manifests identify the recorded source
snapshots; they do not certify the byte contents of these normalized text logs.
