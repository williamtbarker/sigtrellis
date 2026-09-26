# Configuration and study contracts

YAML, TOML and JSON contain the same flat keys. CLI values override file values. Unknown keys, incompatible scales, nonsparse mixing ratios, and unsupported contrasts fail validation. Every run writes its resolved `configuration.json`; that file is itself a valid configuration input.

The complete typed contract and defaults are in `src/sigtrellis/config.py`. Use `examples/biological_study.yaml` as a study-design starting point, not as a universal validated protocol. Choose all scientific settings before inspecting held-out performance.

| Area | Main keys | Contract |
|---|---|---|
| Outcome | `outcome`, `outcome_type`, `positive_class` | Binary, multiclass or continuous; reference order is saved. Binary positive class can be explicit. |
| Replicates | `sample_id`, `group` | Sample identifies a specimen; group identifies the independent donor/biological unit. All group rows move together. |
| Nuisance | `batch`, `covariates` | Batch is a diagnostic/DE term. Numeric or categorical covariates are jointly penalized predictors. Unseen categories fail. |
| Scale | `input_scale`, `normalization` | Counts use logCPM or frozen median ratio. Declared transformed expression uses `none`; upstream validity remains the user's responsibility. |
| Screening | `candidate_method`, `max_features`, `min_count`, `min_prevalence` | All fitted decisions occur in each training fold. `association` ranks with F statistics; it is not count-DE inference. |
| DE | `candidate_method: deseq2`, `supporting_de`, `de_fdr`, `de_pair_group` | Binary raw-count NB contrasts only. Optional paired donor fixed effects must be identifiable. Supporting DE is post-evaluation same-cohort evidence. |
| Validation | `outer_folds`, `inner_folds`, `repeats`, `cv_strategy` | Grouped nesting; optional leave-group/batch-out outer evaluation. Batch holdouts purge overlapping donors. |
| Model | `strengths`, `l1_ratios`, `tuning_rule`, `max_iter`, `tolerance` | Grid over mean-loss lambda and L1 mixing. One-SE is a sparsity heuristic, not a confidence bound. Nonconvergence is recorded. |
| Prediction | `decision_threshold`, `coefficient_tolerance` | Predeclared binary threshold and numerical definition of nonzero. No threshold optimization on outer tests. |
| Stability | `stability_resamples`, `stability_fraction`, `stability_normalizations` | Whole-group subsampling with retuning; optional alternating normalization perturbations. Frequencies are also reported within each normalization. |
| Gates | `selection_threshold`, `sign_threshold`, `outer_selection_threshold`, `min_groups_gate` | Operational candidate-association gates; no feature-level false-discovery guarantee. Final and resampled coefficient directions must agree. |
| Null control | `permutations`, `permutation_scheme`, `permutation_strata`, `permutation_alpha` | Full nested reruns. Group shuffling for constant donor outcomes; within-group shuffling for exchangeable repeated conditions. Strata are only for group shuffling. |
| Correlation | `correlation_threshold`, `correlation_max_features` | Bounded complete-linkage groups from absolute marginal correlation, after evaluation. |
| Cells | `cell_type`, `cell_type_value`, `layer`, `min_cells` | Sum raw counts within specimen/type; omit a type value to run every supplied type separately. |
| Cell QC | `cell_min_counts`, `cell_min_genes`, `mitochondrial_prefix`, `max_mito_fraction` | Fixed thresholds. No organism-specific mitochondrial prefix is assumed. |
| Memory | `chunk_size`, `max_dense_mb` | Cell read chunk and dense pseudobulk-accumulator guard. These do not cap every solver/report allocation. |
| Reproducibility | `seed` | Seeds all stochastic splits/subsamples/solvers; actual groups, fits and source hashes are recorded. |

## Paired experiments

Two specimens from the same donor need distinct `sample_id` values and one shared `group`. A condition contrast varying within donor can use `de_pair_group: true` and `permutation_scheme: within_group` if treatment labels are exchangeable. Do not use within-donor label permutation for a time series or a systematically ordered intervention without justifying exchangeability. Technical replicates should be combined according to the assay design before modeling; two rows do not create two donors.

## Compact panels and independent performance

The candidate table is a discovery product assembled after procedure validation. To evaluate a chosen assay panel, freeze the panel and protocol and evaluate on new donors/cohorts; do not report the earlier procedure CV score as that panel's measured accuracy. `sigtrellis external` evaluates the already fitted final model, which can include more genes than the gated candidate list. It does not silently refit a compact panel.

## Multiple types and normalization perturbations

When every cell type is requested, the CLI records the tested family and divides the global permutation eligibility threshold by the number of requested types. At least `ceil(1/alpha)-1` permutations are needed even to resolve a p-value as small as the resulting threshold. The package will not manufacture a passing result with a coarse permutation grid. This family-level protection is not gene-level FDR control. Comparing many configurations after looking at results requires additional validation.

When alternating normalizations, coefficients are expressed per fit-specific SD; frequency and direction are more comparable than raw coefficient magnitudes. Normalization-specific frequencies must each pass the declared threshold. A normalization that fails for a sparse subsample is an explicit error rather than an omitted inconvenient resample.
