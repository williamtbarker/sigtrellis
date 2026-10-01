# SigTrellis

**Stable candidate transcriptomic signatures, with the biological replicate kept inside the validation boundary.**

SigTrellis is a local Python CLI and library for discovering candidate phenotype-associated genes and transcriptomic features from bulk or single-cell RNA-seq. It accepts user-defined organisms, cell lines, gene identifiers and outcomes. It combines training-only preprocessing and optional differential-expression screening, nested grouped elastic-net modeling, compact-panel validation, stability analysis, negative controls, and an auditable report.

**Version 0.3.0, for research.** A selected feature is a candidate association. Statistical selection does not establish biological validation, mechanism, causality, clinical utility or biomarker qualification. Both the full modeling procedure and an optional compact-panel discovery procedure can be evaluated by nested CV; their final fixed models need untouched external evaluation. Read [LIMITATIONS.md](LIMITATIONS.md).

**Elastic net in plain language:** give the model gene measurements and the phenotype you want to predict. It learns a weight for each gene, while charging a penalty for complexity. One part of that penalty pushes weak weights to zero; the other restrains large weights and helps handle genes that carry overlapping information. The remaining genes form a candidate signature. SigTrellis tests whether that signature predicts unseen biological samples and whether the same genes or correlated groups recur when the data change.

For a new cell line, start with [the practical cell-line guide](docs/CELL_LINE_GUIDE.md). You need independent cultures or experiments and a measurable phenotype with variation. A single culture, an unspecified phenotype, or perfect treatment/batch confounding cannot support a reproducible phenotype signature, regardless of the number of sequenced cells.

The methodological rationale, assumptions, and literature basis are documented in [RESEARCH.md](RESEARCH.md).

## Install

Python **3.12 or newer**. The verified runtime uses Python 3.12; CI is configured for 3.12 and 3.13. No cloud service or R installation is required.

From this source directory:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install '.[de,demo]'
sigtrellis --help
```

`de` adds PyDESeq2. `demo` adds the parser needed to convert the pinned public single-cell deposit. Neither is required for the core elastic-net analysis. A release archive also includes a built wheel; this project is not assumed to be published on PyPI.

## A reproducible first run

```bash
sigtrellis simulate --output data/toy --samples 80 --features 200 --seed 7
sigtrellis bulk \
  --expression data/toy/counts.tsv \
  --metadata data/toy/metadata.tsv \
  --config data/toy/config.yaml \
  --output results/toy
```

Open `results/toy/report.html`. The fixture's truth is recorded separately in `data/toy/truth.json`; the modeling code never reads it. Output directories must be new or empty. A run with no credible candidates still succeeds and writes a report. Invalid data, nonidentifiable designs, or unrecoverable fit errors return exit code 1 and preserve a failed-run manifest where a run has started.

## Bulk RNA-seq

```bash
sigtrellis bulk \
  --expression counts.tsv --metadata metadata.tsv \
  --outcome infection_status --outcome-type binary \
  --positive-class infected --sample-id sample_id \
  --group donor --batch sequencing_batch \
  --config examples/biological_study.yaml \
  --output results/bulk
```

Expression is CSV/TSV, optionally gzipped, with a first identifier column. Both orientations are supported. Auto-detection requires an exact match to metadata sample IDs; ambiguous axes are rejected. Use `--orientation genes_by_samples` or `samples_by_genes` when needed. The first metadata column must match `sample_id`; metadata are aligned by identifier, never guessed from row position.

Raw counts must be finite nonnegative integers; missing counts are rejected. Set `input_scale: log_expression` and `normalization: none` for deliberately supplied transformed expression. `input_scale: features` supports explicitly defined sample-level features. For those two transformed input scales, optional `imputation: median` learns fills independently inside each training fold; all-missing training features are excluded. Upstream transformation safety cannot be inferred from a matrix alone. Gene symbols may be joined with `--annotations genes.csv` containing unique `gene_id,gene_symbol`; they never influence selection.

Binary logistic, multinomial logistic, and continuous elastic-net regression are supported. Class order and positive/reference class are saved. Continuous outcomes must be numeric. Multiclass training folds must contain every class. Covariates are explicit, jointly regularized predictors; donor IDs and batch are excluded from ordinary predictors. Batch enters count-DE designs and robustness diagnostics. Optional paired DE uses `de_pair_group: true`; this is a fixed-effect paired contrast, not a general mixed model.

## Single-cell RNA-seq

```bash
sigtrellis single-cell \
  --input experiment.h5ad --layer counts \
  --outcome infection_status --sample-id specimen \
  --group donor --single-cell-mode distribution \
  --batch sequencing_batch --config examples/single_cell_native.yaml \
  --output results/cell_distributions
```

The chosen matrix or layer must contain **raw counts**. `X` is used if `--layer` is absent. Two representations share the donor/culture-aware modeling engine:

| Mode | Biological question and representation |
|---|---|
| `pseudobulk` | Which genes change in sample-level expression within a declared cell type? Raw counts are summed per specimen/type for count DE and gene modeling. |
| `distribution` | Does phenotype relate to cell-state abundance, detection, variability, or a predefined gene program's distribution? Cells contribute to specimen-level features; they never become independent phenotype replicates. |

The command above uses **no summed-count pseudobulk**. Distribution mode supports `abundance`, `gene_mean`, `gene_detection`, `gene_variance`, fixed-threshold `gene_tail`, declared `gene_correlation` pairs, `program_mean`, `program_variance`, `program_q90`, `program_correlation`, and threshold-defined `program_fraction`. It preserves rare-cell and within-cell coupling information that count sums can erase. These are interpretable specimen-level predictors; cells still do not become independent biological replicates. Add `cell_type` and prespecified `cell_states` to distinguish cell populations. State names, program members and activation thresholds are declared before validation. Every cell's normalization is sample-local; no cohort-trained embedding is hidden in the adapter. Within-sample cell subsampling probes sensitivity to cellular sampling separately from biological replication. Missing state expression is unavailable, not zero. See [the single-cell methods guide](docs/SINGLE_CELL.md) and `examples/replicated_single_cell_native.yaml`. Without an explicit mode, the backward-compatible CLI default remains `pseudobulk`.

Use distinct specimen IDs for conditions/processing aliquots from one donor, and put the shared donor ID in `group`. A donor column can itself be `sample_id` only when one specimen/condition exists per donor. Inconsistent outcome, group or batch metadata within a sample are rejected. Low-cell-count pseudobulks are logged and excluded using the predeclared `min_cells` threshold. Input annotations may carry upstream bias.

Seurat/Matrix Market exports can be converted explicitly:

```bash
sigtrellis import-mtx --matrix matrix.mtx.gz --genes features.tsv.gz \
  --barcodes barcodes.tsv.gz --metadata cell_metadata.tsv --output experiment.h5ad
```

The matrix is genes × cells; feature/barcode files are headerless. Metadata's first column contains exact cell IDs. Other metadata remain strings unless declared using repeated `--numeric-column` arguments. Cells are aligned by barcode and raw counts are checked.

For pseudobulk, select one cell type in advance or run types separately. If all types are run, each is analyzed separately; the CLI applies a Bonferroni threshold across the requested cell-type family for the global permutation gate. This does not provide gene-level FDR control.

## Public data demonstrations

The expanded demonstrations use **86 independent yeast cultures** (66 training, 20 untouched) and a **261-donor single-cell source** with a predefined 137-donor training subset and 96-donor processing-cohort holdout. The current single-cell demonstration uses cell-state distributions and within-cell gene/program correlations from all cells of the selected specimens. Every holdout donor is purged from training, including their other processing aliquots. Exact commands, source terms and measured results are in [PUBLIC_DEMOS.md](docs/PUBLIC_DEMOS.md) and [VALIDATION.md](VALIDATION.md). The large H5AD download is explicit (12.2 GB); the current preparation retains **526,899 training cells and 375,261 held-out cells**, with no per-specimen cap. An explicit cap remains available for smaller practice runs. The original 400-cell and pseudobulk results are retained as historical 0.2.0 evidence.

Smaller examples remain useful for testing paired designs and confounding:

```bash
python examples/prepare_public.py --dataset all --output data/public

sigtrellis bulk \
  --expression data/public/bulk/counts.tsv \
  --metadata data/public/bulk/metadata.tsv \
  --config examples/public_bulk.yaml \
  --output results/public_bulk

sigtrellis single-cell \
  --input data/public/single_cell/kang.h5ad \
  --config examples/public_single_cell.yaml \
  --output results/public_single_cell
```

The bulk airway example has eight samples from four paired donors. The Kang single-cell example starts with 24,673 cells from eight donors, analyzes the supplied CD14+ monocyte annotation, and holds donors together. **Neither small study qualifies a biomarker panel.** Kang condition is tied to the pooled processing library; a count-DE contrast adjusted for that library is nonidentifiable and is reported as unavailable. Strong prediction does not resolve this confounding.

Downloads are checksum-pinned, licenses and changes are documented, and datasets are not committed. See [DATA_SOURCES.md](docs/DATA_SOURCES.md) and the executed results in [VALIDATION.md](VALIDATION.md).

## What is fitted where?

1. Split biological groups into outer training and test folds.
2. Within each outer training set, run inner grouped CV. Each inner training fit learns its own normalization reference, gene filters/screen, scales, covariate encoding and model.
3. Select lambda and L1 mixing ratio from inner proper loss; optionally apply a sparsity-oriented one-SE heuristic.
4. Refit the selected procedure on outer training rows and predict only held-out groups.
5. Repeat the complete nested procedure under label permutation as a negative control.
6. Independently subsample whole biological groups and retune to describe feature stability.
7. Fit a final full-data model, calculate exploratory full-data evidence, and report operational robustness gates.

With `panel_validation: true`, each outer training set performs its own stability selection, applies the predeclared frequency/sign thresholds and `panel_max_features` cap, and refits that panel before predicting the outer test set. Its regularization choice is the modal subsample-tuning result; it is not retuned against inner labels that already selected the panel. The full procedure is repeated under permutation. The final panel and full model have separate numeric states and separate external evaluations. A panel may legitimately contain no genes and predict only the baseline.

No globally selected gene set feeds nested CV. Full-cohort DE, PCA, correlation groups, and coefficient paths are descriptive. The optional `candidate_method: deseq2` performs real count-aware DE **inside every training fold**; `association` is a faster ranking statistic, not a substitute for count-aware inference.

## Elastic-net parameters and outputs

| Parameter/evidence | Meaning |
|---|---|
| `l1_ratios` | L1/L2 mixing; corresponds to glmnet's alpha; values in (0,1] |
| `strengths` | Lambda on mean loss; larger means stronger regularization |
| sklearn regression `alpha` | Equal to SigTrellis lambda |
| sklearn classification `C` | `1/(n_training_samples * lambda)` with weights summing to n |
| Coefficient | Change per training-fold SD of normalized expression; binary log odds, multinomial logit contrast, or continuous response units |
| Selection frequency | Fraction of retuned group subsamples with nonzero coefficient |
| Sign consistency | Dominant sign among selecting fits; zero when never selected |
| Rank/coef variation | Descriptive distributions, not confidence intervals |
| Correlation-group frequency | At least one member selected; no pathway or gene-level error-control claim |
| Robustness gate | Operational evidence threshold; never biological validation |

Tuning minimizes log loss or squared error. Evaluation additionally reports discrimination, average precision, balanced accuracy, sensitivity, specificity, MCC, Brier/calibration, or RMSE/MAE/R². Samples are weighted so each biological group has equal total weight. Classification thresholds are predeclared, not chosen from outer test outcomes.

Default stability uses 50 subsamples and 99 permutations. These are starting settings. `examples/biological_study.yaml` uses more repetitions and permutations. Permutations must be exchangeable for the actual study: use `group` for donor-constant outcomes or `within_group` only for exchangeable repeated conditions. Do not permute an observational time series as independent rows. With covariates, global permutations do not establish conditional gene significance and positive robustness gates are withheld.

`cv_strategy: temporal` with a numeric `time` column supplies forward validation. Training groups must end before the earliest test time minus `temporal_gap`. Inner splits obey the same rule. Unrestricted label permutations are not assumed valid for a time series, so temporal runs retain exploratory evidence status. Survival/censoring models are outside the supported outcome contract.

## Results and audit trail

```text
results/
  run_manifest.json             configuration, versions, hashes, warnings, outcomes
  configuration.json            reusable resolved configuration
  qc_report.json / .html        input contract and biological replicate counts
  model_metrics.json            nested performance, baselines, permutation and batch checks
  biomarkers.csv                ranked contrast-specific evidence and gates
  biomarker_stability.csv        resampling frequencies, signs, coefficients and ranks
  coefficients.csv              final refit coefficients, with explicit directions
  cv_results.csv                inner losses, hyperparameters and chosen settings
  out_of_fold_predictions.csv    every held-out prediction
  audit.json                    training IDs, splits, filtering and fitted-transform hashes
  feature_universe.json          one shared feature universe; exclusions are reconstructable
  permutation_audit.json         permuted labels, splits, fit evidence
  stability_resamples.json       rows and tuned settings for each perturbation
  resample_coefficients.npz      numeric coefficients, no pickled objects
  correlation_groups.csv        group stability
  feature_substitution.csv      pairwise substitution/co-selection
  model_state.json               portable frozen numeric model, no pickle
  report.md / .html              scientific report; HTML embeds all images
  figures/                      PCA, path, tuning, ROC/PR/calibration, stability,
                                correlations, permutations, and DE when identifiable
```

Some additional artifacts are conditional, including count-DE results and pseudobulk counts. Exact output hashes are saved in the final manifest. No timestamps enter the numerical evidence. Identical versions, inputs, configuration and hardware produce deterministic seeded fits; portable numerical identity across BLAS/platform versions is not guaranteed. HDF5/gzip container bytes may include serialization metadata, so compare scientific values as well as raw file hashes.

Distribution runs add `feature_schema.json`, `sample_features.tsv.gz`, `cell_distribution_evidence.csv`, `cell_quality_by_sample.csv`, and a specimen-level distribution figure; `feature_id`, source `gene_id`, gene/program partner, threshold, cell type, program, kind and unit remain distinct. Compact-panel runs add `panel.csv`, `panel_state.json`, panel OOF predictions and panel metrics. DE supports binary contrasts, continuous outcome slopes, and multinomial reference contrasts. Multiclass DE controls the gene-by-contrast family; multiple distribution cell types share a DE family. DE effect units differ from predictive coefficient units and are explicitly recorded.

## Frozen external validation

```bash
sigtrellis external --run results/bulk \
  --expression external_counts.tsv --metadata external_metadata.tsv \
  --output results/external
```

Add `--model-kind panel` to evaluate the frozen compact panel. For single-cell inputs, use `--input external.h5ad` in place of the expression/metadata pair. The stored representation and training transforms are reused.

The external cohort must have the exact training gene/feature universe and declared metadata. A compact predictor panel is not yet a small targeted assay: count normalization still requires the original gene universe. Sample/group overlap and exact copied expression profiles are rejected, including copies under renamed IDs. Distribution inputs use original cell-profile fingerprints when available; matching low-dimensional summaries alone are not treated as copied specimens. Exact reordered/renamed copies of entire cell collections are detected. Partial overlaps, near-duplicates, changed annotations and biological relatives still require study-specific review. Native distribution evaluation also checks the full RNA normalization gene universe. Native distribution models from 0.2.0 must be refitted to supply this contract. Covariate categories absent from training are rejected.

Before prediction, the model, configuration and training metadata must match their saved hashes. Modified artifacts fail verification. External probability columns identify their classes explicitly. Hash verification assumes a trusted manifest; it is not a digital signature.

## Develop and verify

```bash
python -m pip install -e '.[dev,de,demo]'
python -m ruff check .
python -m ruff format --check .
python -m mypy src
python -m pytest --cov=sigtrellis --cov-report=term-missing
python -m build
```

Tests include planted counts, correlated replacements, pure noise, imbalance, outliers, batch-only and batch-specific effects, paired cells, too few donors, a deliberately unsafe global screen, exact duplicates, target copies, and train-boundary spies. See [VALIDATION.md](VALIDATION.md) for measured results and [ARCHITECTURE.md](ARCHITECTURE.md) for extension boundaries.

Version 0.3.0 adds within-cell correlations, rare-tail features, order-invariant
cell perturbations/duplicate checks, strict external normalization compatibility,
streamed all-cell preparation, and specimen-level distribution/QC reports.
See [REVIEW_0.3.0.md](docs/REVIEW_0.3.0.md) for verified counterexamples and fixes,
[VALIDATION.md](VALIDATION.md) for executed tests and public results, and
[RELEASE_GUIDE.md](RELEASE_GUIDE.md) for the packaged evidence and reproduction.
A richer representation can expose different signal; it does not guarantee better
cross-cohort prediction or a biomarker that passes the robustness gates.

MIT software license. Third-party public datasets retain their own licenses. Cite the software with `CITATION.cff` and the methods and datasets used in your analysis.
