# Changelog

## 0.3.0 — 2026-09-26

- Add fixed expression-tail fractions and within-cell gene/program correlation predictors without count pseudobulk or cell pseudoreplication.
- Make cell perturbations invariant to row order, chunking and unrelated-donor removal; detect renamed/reordered exact copies of cell collections.
- Enforce the complete RNA normalization-universe contract for frozen native single-cell evaluation; legacy distribution models require refitting.
- Reject misaligned perturbation matrices and ambiguous joint-feature configuration; preserve canonical program-pair identity across serialization.
- Stream cellular moments without retaining score arrays unless exact quantiles are requested; stream public AnnData subset preparation and support all cells per specimen.
- Report specimen-level feature distributions, technical cell-quality summaries, explicit partner genes/programs and thresholds; do not attach mean-count DE to correlation features.
- Support a single informative predictor and missing low-information variance/correlation estimates.
- Preserve reconstructable filtering decisions using a shared feature universe instead of repeating excluded features in every fit record.
- Add 23 adversarial/scientific tests and rerun public bulk and native all-cell workflows; retain earlier release evidence with its original labels.

## 0.2.0 — 2026-09-26

- Add fixed cell-state abundance, gene detection/mean/variance, and predefined program mean/variance/quantile/activation features, with explicit units and biological sample identity.
- Add within-sample cell perturbations, per-perturbation stability gates, sparse chunk processing and bounded sample-level accumulators.
- Fix capture-depth dependence of abundance pseudocounts using a constant offset on proportions; test unequal-cell-count negative controls and reject protected metadata as cell-state labels.
- Select/refit compact panels inside outer training sets, permutations and batch holdouts. Save panel-specific coefficients, gates, numeric states and external evaluations.
- Add fold-local median imputation for transformed features; preserve missing-state meaning and reject missing raw counts.
- Extend count DE to continuous slopes and multiclass reference contrasts, with explicit effects and appropriate joint BH families.
- Add forward temporal group splitting and conservative exchangeability reporting.
- Add strict Matrix Market/Seurat-export import and frozen native single-cell external evaluation.
- Add pinned replicated yeast cultures and a larger donor-level single-cell processing-cohort holdout. Correct repeated-processing specimen identity before purging held-out donors.
- Expand adversarial tests and documentation for new cell lines, feature interpretation and independent verification; preserve earlier release evidence as historical.

## 0.1.1 — 2026-09-26

- Preserve sample, donor, batch, stratum and categorical outcome strings during CSV/TSV parsing, including leading zeros and literal `NA`; reject blank required metadata.
- Canonicalize signed zeros before fingerprinting so numerically identical specimens cannot evade exact-duplicate safeguards.
- Restrict paired count DE to one observation per biological group/condition; reject unsupported repeated-condition designs before backend fitting.
- Verify saved model, configuration and training metadata against the run manifest before frozen external evaluation. Export explicitly labeled probability columns and model fingerprints.
- Keep correlation-group and substitution frequencies separate for every multinomial contrast.
- Block nonfinite statistics and invalid permutation p-values; report missing-class discrimination and constant-outcome R² as undefined.
- Validate H5AD matrix dimensions against cell/gene metadata before aggregation; reject boolean numeric settings and duplicated hyperparameter-grid entries.
- Preserve batch tuning choices, complete tuning tables, coefficients and batch-axis labels; correct filtering-stage counts and apply coefficient tolerance consistently in manifests.
- Add 27 targeted regression tests, including a real nested PyDESeq2 screen, covering identifier handling, artifact verification, multinomial evidence, nonfinite statistics, and H5AD validation.
- Refresh public/synthetic reports and package verification for the revised release.

## 0.1.0 — 2026-09-26

- Original bulk and chunked H5AD pseudobulk adapters with biological-group contracts.
- Binary, multinomial, and continuous elastic net with nested grouped validation.
- Training-only normalization, filtering, optional association/DE screening and scaling.
- Retuned biological-group subsampling, perturbation-specific frequencies, sign/rank evidence, correlated-feature substitution diagnostics.
- Full-procedure permutation controls and purged batch holdouts.
- Optional binary PyDESeq2 with explicit paired/covariate/batch design checks.
- JSON/CSV/NPZ evidence, portable model state, frozen external evaluation and self-contained scientific reports.
- Synthetic adversarial tests, pinned public bulk/single-cell preparation, release documentation and CI configuration.
