# Changelog

## 0.1.1 — 2026-09-26

- Preserve sample, donor, batch, stratum and categorical outcome strings during CSV/TSV parsing, including leading zeros and literal `NA`; reject blank required metadata.
- Canonicalize signed zeros before fingerprinting so numerically identical specimens cannot evade exact-duplicate safeguards.
- Restrict paired count DE to one observation per biological group/condition; reject unsupported repeated-condition designs before backend fitting.
- Verify saved model, configuration and training metadata against the run manifest before frozen external evaluation. Export explicitly labeled probability columns and model fingerprints.
- Keep correlation-group and substitution frequencies separate for every multinomial contrast.
- Block nonfinite statistics and invalid permutation p-values; report missing-class discrimination and constant-outcome R² as undefined.
- Validate H5AD matrix dimensions against cell/gene metadata before aggregation; reject boolean numeric settings and duplicated hyperparameter-grid entries.
- Preserve batch tuning choices, complete tuning tables, coefficients and batch-axis labels; correct filtering-stage counts and apply coefficient tolerance consistently in manifests.
- Add 27 targeted review tests, including a real nested PyDESeq2 screen. The original release fails 26 of the new adversarial cases; the original fold-local DE integration passes.
- Refresh public/synthetic reports and package verification for the revised release. See `docs/REVIEW_0.1.1.md`.

## 0.1.0 — 2026-09-26

- Original bulk and chunked H5AD pseudobulk adapters with biological-group contracts.
- Binary, multinomial, and continuous elastic net with nested grouped validation.
- Training-only normalization, filtering, optional association/DE screening and scaling.
- Retuned biological-group subsampling, perturbation-specific frequencies, sign/rank evidence, correlated-feature substitution diagnostics.
- Full-procedure permutation controls and purged batch holdouts.
- Optional binary PyDESeq2 with explicit paired/covariate/batch design checks.
- JSON/CSV/NPZ evidence, portable model state, frozen external evaluation and self-contained scientific reports.
- Synthetic adversarial tests, pinned public bulk/single-cell preparation, release documentation and CI configuration.
