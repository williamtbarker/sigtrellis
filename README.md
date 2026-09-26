# SigTrellis

**Stable candidate transcriptomic signatures, with the biological replicate kept inside the validation boundary.**

SigTrellis is a local Python CLI and library for bulk RNA-seq and donor-aware single-cell pseudobulk analysis. It combines fold-local preprocessing and optional differential-expression screening, nested grouped elastic-net modeling, retuned group subsampling, correlated-feature diagnostics, negative controls, and an auditable scientific report.

**Research software, version 0.1.0.** A selected gene is a candidate association—not a validated biomarker, mechanism, causal effect, diagnostic test, or qualified clinical measurement. The final consensus panel has not inherited the nested CV performance of the training procedure. Read [LIMITATIONS.md](LIMITATIONS.md) before interpreting results.

This is an original, general-purpose implementation of public methods. It does not reproduce any previous employer's code, data, or confidential methodology. The [research review](RESEARCH.md) explains the design and the limited conceptual relationship to Stabilomics.

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

Raw counts must be finite nonnegative integers. Set `input_scale: log_expression` and `normalization: none` only for a deliberately supplied transformed matrix; its upstream leakage safety cannot be verified. Missing values are rejected. Gene symbols may be joined with `--annotations genes.csv` containing unique `gene_id,gene_symbol`; they never influence selection.

Binary logistic, multinomial logistic, and continuous elastic-net regression are supported. Class order and positive/reference class are saved. Continuous outcomes must be numeric. Multiclass training folds must contain every class. Covariates are explicit, jointly regularized predictors; donor IDs and batch are excluded from ordinary predictors. Batch enters count-DE designs and robustness diagnostics. Optional paired DE uses `de_pair_group: true`; this is a fixed-effect paired contrast, not a general mixed model.

## Single-cell RNA-seq

```bash
sigtrellis single-cell \
  --input experiment.h5ad --layer counts \
  --outcome infection_status --sample-id specimen \
  --group donor --cell-type cell_type --cell-type-value Monocytes \
  --batch sequencing_batch --config examples/biological_study.yaml \
  --output results/monocytes
```

The chosen matrix or layer must contain **raw counts**. `X` is used if `--layer` is absent. Cells are summed per experimental sample and cell type, in chunks. The result is a sample-level model. Cells are never independent phenotype replicates.

Use distinct specimen IDs for two conditions from one donor, and put the shared donor ID in `group`. A donor column can itself be `sample_id` only when one specimen/condition exists per donor. Inconsistent outcome or donor metadata within a sample are rejected. Low-cell-count pseudobulks are logged and excluded using the predeclared `min_cells` threshold. Input cell-type annotations are treated as supplied; learned annotations or integrated data may carry upstream bias. Export Seurat raw counts and metadata to H5AD or use already aggregated sample counts through `bulk`.

Selecting one cell type in advance is recommended. If all types are run, each is analyzed separately; the CLI applies a Bonferroni threshold across the requested cell-type family for the global permutation gate. This does not provide gene-level FDR control.

## Public data demonstrations

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

## Frozen external validation

```bash
sigtrellis external --run results/bulk \
  --expression external_counts.tsv --metadata external_metadata.tsv \
  --output results/external
```

The external cohort must have the exact training gene universe and declared metadata. Training normalization and scaling are reused without fitting. Sample/group overlap and exact copied expression profiles are rejected, including copies under renamed IDs. This evaluates the frozen final model, not a new compact panel refitted from the stability table. Near-duplicates and biological relatives still require study-specific review. Covariate categories absent from training are rejected in both CV and external evaluation.

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

MIT software license. Third-party public datasets retain their own licenses. Cite the software with `CITATION.cff` and the methods and datasets used in your analysis.
