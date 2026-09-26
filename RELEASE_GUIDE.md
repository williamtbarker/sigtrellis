# SigTrellis 0.2.0 release guide

The package is ready for public review as general-purpose **research software**.
It discovers candidate associations for a supplied phenotype across independent
cultures or donors. It cannot guarantee a valid biomarker from every experimental
design. Nothing has been pushed to GitHub or uploaded to PyPI by this preparation.

## How elastic net finds a candidate

Imagine measuring thousands of genes in independently grown cultures, some resistant
and some sensitive. The model predicts that label by assigning each gene a weight.
It pays a penalty for complexity: one part pushes unnecessary weights to zero;
another restrains large weights and helps handle overlapping information. Genes
with remaining weights become candidate predictors. The weight's direction describes
an association within this model, not a causal effect.

SigTrellis repeatedly changes the training cultures, retunes the model, and tests
predictions on cultures it did not train on. A useful candidate should recur with
consistent direction, help predict unseen samples and survive appropriate batch
and negative-control checks. Correlated genes may substitute for each other, so
the report also tracks their group. Laboratory experiments and independent cohorts
provide additional evidence that a coefficient alone cannot establish.

## Repository contents

| Path | Purpose |
|---|---|
| `README.md`, `RESEARCH.md`, `ARCHITECTURE.md` | Entry point, primary-source rationale, design |
| `VALIDATION.md`, `LIMITATIONS.md` | Measured execution and scientific boundaries |
| `src/sigtrellis/` | 25 typed modules plus `py.typed`; no notebook-only logic |
| `tests/` | Unit, integration, scientific and leakage-adversarial tests |
| `examples/` | Public/synthetic preparation, culture configurations, runnable studies |
| `docs/` | Cell-line/single-cell guides, exact commands, source terms, review and verification |
| `.github/workflows/`, `Dockerfile` | CI and local container build configuration |
| `pyproject.toml`, `LICENSE`, `CITATION.cff`, `CHANGELOG.md` | Installation, MIT, citation and history |
| `dist/` | Built 0.2.0 wheel and source distribution in the archive |
| `release_evidence/` | Reports and exact execution source snapshots in the archive; ignored by Git |
| `sigtrellis-history.bundle` | Local Git history in the archive |

Raw public downloads, virtual environments and caches are excluded.

## Architecture and methods

Modality-specific adapters produce a typed specimen-level dataset. Bulk keeps the
supplied gene matrix. Single-cell either sums counts per specimen/state or computes
fixed state abundance, gene detection/variability and predefined program distribution
features. Cells do not become phenotype replicates. Biological groups control splits
and group weighting.

The shared engine nests learned preprocessing and optional DE/association screening
inside grouped CV. Binary/multinomial logistic elastic net uses SAGA; continuous
outcomes use elastic-net regression. Both regularization and L1 mixing are tuned
on inner proper loss. Group resampling, measurement perturbations and retuning
characterize selection/sign/rank stability. Compact-panel selection repeats inside
each outer training set and has its own predictive, permutation and batch evaluation.

PyDESeq2 supplies binary, continuous and reference-class multiclass count contrasts
with identifiable fixed covariates/pairing. Full-data DE is supporting/exploratory;
it never prefilters reported CV globally. Correlation/substitution diagnostics track
shared signal. Manifests and frozen numeric states preserve the fitted contract.

## Measured validation

The clean installed-wheel suite passes **152 tests**, with 87.10% branch-aware
coverage. Ruff and strict mypy pass. The installed toy recovers all four planted
genes. Controls include pure noise, confounding, duplicates, supervised-screen
leakage, annotation leakage and unequal single-cell capture depth.

| Public compact panel | Held-out data | AUC | Specificity at 0.5 | Interpretation |
|---|---|---:|---:|---|
| Yeast, 10 genes | 20 cultures | 1.000 | 1.000 | Strong same-study genotype contrast; no cross-lab claim |
| Monocytes, 3 genes | 96 donors | 0.747 | 0.227 | Poor specificity and probability transfer |
| Cell-state/program distributions, 5 features | 96 donors | 0.804 | 0.500 | Cross-batch gates fail; no fully supported candidates |

The distribution holdout was reused after a synthetic correctness fix and is
engineering verification. All holdouts are within their source studies. See
`VALIDATION.md` for full/panel metrics, gates, source versions and limitations.

## Exact commands

Run inside the unpacked repository. Python 3.12 is the verified runtime.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -c docs/validation/v0.2.0/runtime_constraints.txt '.[de,demo]'

sigtrellis simulate --output data/toy --samples 80 --features 200 --seed 7
sigtrellis bulk --expression data/toy/counts.tsv --metadata data/toy/metadata.tsv \
  --config data/toy/config.yaml --panel-validation --output results/toy

python examples/prepare_replicated.py --dataset yeast --output data/public/yeast
sigtrellis bulk --expression data/public/yeast/train/counts.tsv.gz \
  --metadata data/public/yeast/train/metadata.tsv \
  --config examples/replicated_bulk.yaml --output results/replicated_bulk
sigtrellis external --run results/replicated_bulk --model-kind panel \
  --expression data/public/yeast/external/counts.tsv.gz \
  --metadata data/public/yeast/external/metadata.tsv \
  --output results/replicated_bulk_panel_external

# Explicit 12.2 GB download; alternatively pass --source /path/to/source.h5ad.
python examples/prepare_replicated.py --dataset lupus --download --output data/public/lupus
sigtrellis single-cell --input data/public/lupus/train.h5ad \
  --config examples/replicated_single_cell.yaml --output results/replicated_single_cell
sigtrellis external --run results/replicated_single_cell --model-kind panel \
  --input data/public/lupus/external.h5ad --output results/replicated_single_cell_panel_external
sigtrellis single-cell --input data/public/lupus/train.h5ad \
  --config examples/replicated_cell_distributions.yaml --output results/replicated_cell_distributions
sigtrellis external --run results/replicated_cell_distributions --model-kind panel \
  --input data/public/lupus/external.h5ad --output results/replicated_cell_distributions_panel_external
```

Full-model external commands are in `docs/PUBLIC_DEMOS.md`. Use fresh output
directories. Nested/permutation analyses take substantial CPU time; the configuration
guide explains the fitting budget. Open each `report.html`. For a new cell line,
edit `examples/cell_line.yaml` and follow `docs/CELL_LINE_GUIDE.md`.

## Limits and remaining research

Operational gates are not formal feature-level error control. Penalized covariates
do not provide conditional gene significance. Cell variability includes technical
noise, state abundance is compositional, and annotations can carry upstream bias.
A compact RNA-seq predictor still needs its full gene universe for normalization;
it is not automatically a small targeted assay. Dataset terms remain separate from MIT.

Research priorities are independent cross-laboratory cohorts, calibration under
prevalence shift, conditional/correlated-group inference, deconvolved cell variability,
frozen assignment of learned states and targeted-assay measurement. Cell-level mixed
models, neighborhood inference and neural bag models are distinct future backends.

Publication as transparent research software is justified by the implemented
boundaries and executed evidence. Clinical qualification, causal interpretation
and guaranteed biomarker recovery from every arbitrary cell line are not justified.
