# Architecture

## Scientific data boundary

`Dataset` contains a sample-by-feature frame, exactly aligned metadata, typed `Feature` identities, optional cell-type count tables and cell-perturbed feature frames, gene-symbol mappings, input hashes and upstream QC. Metadata cannot silently become features. Bulk aligns either orientation. Single-cell adapters read sparse HDF5 chunks and use fixed cell QC. Pseudobulk produces count sums; distribution mode produces sample-local moments/program scores and relative abundance. Both share the same downstream machinery without treating their upstream estimands as identical.

## Components

| Module | Responsibility |
|---|---|
| `config` | Frozen configuration, YAML/TOML/JSON loading, unknown-key rejection |
| `domain` | Dataset, split/audit types, hashes and strict JSON serialization |
| `io`, `singlecell`, `cell_features`, `matrix_import` | Bulk alignment, count sums, typed sample-local distributions, strict sparse export import |
| `qc` | Integrity, replicate, target-copy, duplicate and confounding checks |
| `preprocessing` | Training-only normalizer, gene screen, scaler and covariate encoder |
| `de` | Rank-checked binary, continuous and multiclass-reference PyDESeq2 contrasts |
| `splits` | Stratified group splitting, mixed-label group splitting, purged batch splits, subsampling, permutation |
| `modeling` | Solver wrappers, parameter mapping, inner tuning, coefficient extraction, frozen numeric state |
| `validation` | Nested out-of-group prediction and full-procedure permutation control |
| `stability` | Retuned group subsampling, coefficient/sign/rank/frequency/Jaccard evidence |
| `panel` | Compact feature selection/refit policy inside each outer training boundary |
| `evidence` | Scientific feature identities and post-validation count-DE joins |
| `correlation` | Bounded descriptive clustering and substitution diagnostics |
| `metrics` | Proper-loss selection, group weights, complementary metrics/calibration |
| `workflow` | Orchestration, positive-gate eligibility, manifests, errors and artifacts |
| `prediction` | Frozen external-cohort evaluation without transform refitting |
| `integrity` | Verify required frozen-run artifacts against their recorded hashes |
| `reporting` | Purpose-labeled plots and self-contained HTML/Markdown |
| `simulate` | Deterministic fixtures with separately recorded planted truth |
| `cli` | Explicit user-facing contracts and modality routing |

## Validation state ownership

Each training boundary owns a fresh `Prepared` object. It holds the normalization reference, prevalence/variance/candidate decisions, training scales, and covariate categories. Transforming validation rows cannot modify it. Hyperparameter candidates within the **same** training fold share that fitted preprocessing safely; no candidate can access the validation labels except through its scalar scoring step.

Outer folds are never used to choose their own hyperparameters. Final full-data tuning, stability selection and exploratory analyses run after the outer evaluation and do not replace its inputs. The final model is a new fit; its exact performance on an untouched cohort remains unknown until the external command is used.

Panel validation owns a fresh stability engine in each outer training set. It emits separate predictions, coefficients, selected features and tuning records. Batch and permutation evaluations run the same panel policy. Cell perturbations travel inside `Dataset.subset`, so a training subset cannot acquire held-out sample rows. Their definitions use no cross-sample fitted statistics. Imputation and scaling remain within `Prepared`, including when native external H5AD is converted using the frozen configuration.

`Audit` records every observed-data fit and every training/test boundary. Permutation controls reconstruct splits from permuted labels and preserve separate audit records. Fingerprints describe normalization/scaling state; row IDs and gene inclusion/exclusion lists make the boundary inspectable. Output manifests are marked `running`, then `complete` or `failed`. Nonempty output directories are never overwritten.

## Statistical scope

- Independent samples: sample ID is the default group, with an explicit assumption warning.
- Multiple specimens from a donor: an explicit shared group controls splitting, weights and subsampling.
- Paired DE: donor fixed effects are included only when requested and identifiable; repeated groups without a pairing contract are rejected.
- Paired DE additionally requires one observation per group/condition. Technical replicates must be combined appropriately; general repeated-measures designs require another model.
- Covariates: typed numeric or categorical, jointly penalized in prediction and unpenalized design terms in count DE. These are different statistical roles and are labeled accordingly.
- Batch: diagnostics, optional purged holdouts, and count-DE nuisance terms. No global batch correction or batch predictor.
- Multiclass: multinomial logistic elastic net; coefficient contrasts against class 0.
- Correlated genes: descriptive groups, not fitted latent features or biological pathway assignments.
- Multinomial correlation-group and substitution frequencies retain the explicit coefficient-contrast axis. Group membership can be shared while selection evidence remains contrast specific.

CSV identifiers are parsed as strings before constructing the index; a numeric-looking identifier is never silently normalized. Input metadata blanks and H5AD matrix/metadata shape disagreement fail before modeling. Numeric fingerprints canonicalize signed zero. External evaluation verifies model state, configuration and training metadata against the saved manifest before reading the new cohort. This detects changed artifacts against a trusted manifest; it is not cryptographic authentication of a manifest an attacker can also replace.

## Performance and scale

HDF5 cell matrices remain sparse/backed and are read in `chunk_size` blocks. Dense memory scales primarily with sample × cell type × genes. `max_dense_mb` bounds the accumulator, and a single cell type can be selected. Model fitting uses dense sample-level arrays and a training-local feature cap; this is not an out-of-core solver for millions of biological donors. Correlation analysis is capped separately to avoid quadratic all-gene memory. Repeated DE and nested permutations are intentionally expensive; manifests expose the actual work.

BLAS/OpenMP threads are limited to one within `run_analysis` for reproducibility and to avoid oversubscription. CV is sequential in this release. No model is loaded from pickle and no model-generated code is executed. Analysis commands make no network calls; public downloads are confined to an explicit example script.

## Extension rules

New supervised screens must implement the training-only `Prepared.fit` boundary and pass the boundary-spy and global-screen attack tests. New normalization must separate fitted reference state from sample-local transform. New validators must preserve biological groups and declare their estimand. New DE backends must expose design rank, replicate assumptions, effect contrasts and missing/filtered tests rather than mapping all methods onto a misleading universal p-value schema.

Deferred extensions include conditional/random-effect inference, validated continuous/multiclass count contrasts, temporal splits, exchangeability diagnostics, feature-level error control, train-frozen batch harmonization, pathway-score feature adapters, and panel-specific nested stability selection. They require statistical design and validation, not merely more configuration switches.
