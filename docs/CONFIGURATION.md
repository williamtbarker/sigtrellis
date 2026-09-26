# Configuration and study contracts

YAML, TOML and JSON contain the same flat keys. CLI values override file values. Unknown keys, incompatible scales, nonsparse mixing ratios, and unsupported contrasts fail validation. Every run writes its resolved `configuration.json`; that file is itself a valid configuration input.

The complete typed contract and defaults are in `src/sigtrellis/config.py`. Use `examples/biological_study.yaml` as a study-design starting point, not as a universal validated protocol. Choose all scientific settings before inspecting held-out performance.

Compact-panel permutation validation is expensive: every null run repeats the
panel-selection policy, including subsample tuning. Its dominant fit count grows
approximately as `(permutations + 1) × outer_folds × repeats × stability_resamples
× inner_folds × grid_size`. The larger example configurations can require hours
on one CPU; cell-line study templates with larger budgets can take longer. Use a
small exploratory run to verify inputs, then fix the scientific budget before
evaluating the reserved cohort. The core limits numerical threads for reproducibility;
it does not silently substitute a cheaper, differently defined null test.

| Area | Main keys | Contract |
|---|---|---|
| Outcome | `outcome`, `outcome_type`, `positive_class` | Binary, multiclass or continuous; reference order is saved. Binary positive class can be explicit. |
| Replicates | `sample_id`, `group` | Sample identifies a specimen; group identifies the independent donor/biological unit. All group rows move together. |
| Nuisance | `batch`, `covariates` | Batch is a diagnostic/DE term. Numeric or categorical covariates are jointly penalized predictors. Unseen categories fail. |
| Scale | `input_scale`, `normalization` | Counts use logCPM or frozen median ratio. Declared transformed expression uses `none`; upstream validity remains the user's responsibility. |
| Screening | `candidate_method`, `max_features`, `min_count`, `min_prevalence` | All fitted decisions occur in each training fold. `association` ranks with F statistics; it is not count-DE inference. |
| DE | `candidate_method: deseq2`, `supporting_de`, `de_fdr`, `de_pair_group` | Binary, continuous or reference-class multiclass raw-count NB contrasts. Paired donor fixed effects must be identifiable. Supporting DE is same-cohort evidence. |
| Validation | `outer_folds`, `inner_folds`, `repeats`, `cv_strategy` | Grouped nesting; optional leave-group/batch-out outer evaluation. Batch holdouts purge overlapping donors. |
| Model | `strengths`, `l1_ratios`, `tuning_rule`, `max_iter`, `tolerance` | Grid over mean-loss lambda and L1 mixing. One-SE is a sparsity heuristic, not a confidence bound. Nonconvergence is recorded. |
| Prediction | `decision_threshold`, `coefficient_tolerance` | Predeclared binary threshold and numerical definition of nonzero. No threshold optimization on outer tests. |
| Stability | `stability_resamples`, `stability_fraction`, `stability_normalizations` | Whole-group subsampling with retuning; optional alternating normalization perturbations. Frequencies are also reported within each normalization. |
| Gates | `selection_threshold`, `sign_threshold`, `outer_selection_threshold`, `min_groups_gate` | Operational candidate-association gates; no feature-level false-discovery guarantee. Final and resampled coefficient directions must agree. |
| Null control | `permutations`, `permutation_scheme`, `permutation_strata`, `permutation_alpha` | Full nested reruns. Group shuffling for constant donor outcomes; within-group shuffling for exchangeable repeated conditions. Strata are only for group shuffling. |
| Correlation | `correlation_threshold`, `correlation_max_features` | Bounded complete-linkage groups from absolute marginal correlation, after evaluation. |
| Cells | `cell_type`, `cell_type_value`, `layer`, `min_cells` | Pseudobulk uses type/value; native distribution uses declared `cell_states` and rejects `cell_type_value`. Raw counts are required by both adapters. |
| Cell QC | `cell_min_counts`, `cell_min_genes`, `mitochondrial_prefix`, `max_mito_fraction` | Fixed thresholds. No organism-specific mitochondrial prefix is assumed. |
| Memory | `chunk_size`, `max_dense_mb` | Cell read chunk and dense pseudobulk-accumulator guard. These do not cap every solver/report allocation. |
| Reproducibility | `seed` | Seeds all stochastic splits/subsamples/solvers; actual groups, fits and source hashes are recorded. |
| Imputation | `imputation`, `max_missing_fraction` | Optional training-only medians for transformed expression/features; missing raw counts are rejected. |
| Compact panel | `panel_validation`, `panel_max_features` | Stability-based selection/refit inside every outer training, permutation and batch boundary. |
| Temporal | `cv_strategy: temporal`, `time`, `temporal_gap`, `temporal_train_fraction` | Numeric group-aware forward splitting; training ends before testing begins, including inner CV. One forward repeat only. |
| Cell representation | `single_cell_mode`, `cell_states`, `feature_blocks`, `feature_genes` | Pseudobulk or predefined sample-local distributions. State/gene identities must be explicit. |
| Programs | `programs`, `program_thresholds` | Maps names to fixed gene lists; activation fractions require a predeclared threshold per program. |
| Joint/tail cell features | `gene_pairs`, `gene_thresholds` | Prespecified pairs for within-cell Pearson correlation and fixed positive log1p(CP10K) thresholds for high-expression fractions. Duplicate/self pairs fail. Program correlations use all pairs of declared programs. |
| Cellular perturbations | `cell_resamples`, `cell_fraction` | Distribution mode only. Perturbed matrices must come from the adapter; every perturbation requires at least one biological-group stability resample. |

## Paired experiments

Two specimens from the same donor need distinct `sample_id` values and one shared `group`. A condition contrast varying within donor can use `de_pair_group: true` and `permutation_scheme: within_group` if treatment labels are exchangeable. Do not use within-donor label permutation for a time series or a systematically ordered intervention without justifying exchangeability. Technical replicates should be combined according to the assay design before modeling; two rows do not create two donors.

The paired DE backend requires at most one specimen per group/condition. Repeated visits in the same condition cannot be treated as additional independent observations by setting `de_pair_group`. Use a design-specific repeated-measures model for that inference. Grouped predictive validation can still keep such rows together, subject to its stated estimand and exchangeability assumptions.

Identifiers are exact strings: `001`, `01` and `1` denote different units, and a literal identifier `NA` is preserved. Blank required values fail. Numeric settings cannot be booleans; hyperparameter grids cannot contain duplicate values, which would otherwise duplicate fold evidence in the tuning summary.

## Compact panels and independent performance

The candidate table is a discovery product. With `panel_validation: true`, panel discovery is independently evaluated inside outer folds, and the final panel is saved in `panel_state.json`. Use `sigtrellis external --model-kind panel` on untouched specimens to evaluate that fixed panel. Without the flag, external evaluation uses the full model. Neither path silently reselects features. `final_panel_size` is the retained predictor set; `final_nonzero_feature_count` records actual nonzero weights. Count normalization still needs the full gene universe, so a targeted assay needs additional measurement validation.

## Multiple types and normalization perturbations

When every cell type is requested, the CLI records the tested family and divides the global permutation eligibility threshold by the number of requested types. At least `ceil(1/alpha)-1` permutations are needed even to resolve a p-value as small as the resulting threshold. The package will not manufacture a passing result with a coarse permutation grid. This family-level protection is not gene-level FDR control. Comparing many configurations after looking at results requires additional validation.

When alternating normalizations, coefficients are expressed per fit-specific SD; frequency and direction are more comparable than raw coefficient magnitudes. Normalization-specific frequencies must each pass the declared threshold. A normalization that fails for a sparse subsample is an explicit error rather than an omitted inconvenient resample.
