# SigTrellis 0.3.0 release guide

This release supports bulk RNA-seq and a native single-cell distribution workflow
that does not sum cells into pseudobulk. It is a general-purpose research tool for
candidate phenotype-associated signatures. It is not a pretrained universal
biomarker classifier, a clinical assay, or a causal discovery method.

## What the single-cell model learns

The model can ask whether phenotype relates to a cell state's abundance, a gene's
mean/detection/variability, rare high-expression cells, predefined program
activity, or genes/programs varying together within the same cells. Each specimen
contributes one interpretable feature vector. Related specimens stay together in
all validation splits. Individual cells improve the measurement of that vector;
they do not become extra independent donors.

Elastic net learns a weight for each feature. Its penalty shrinks weak weights
toward zero while accommodating overlapping predictors. A nonzero weight is only
an initial association. Nested testing, biological-group resampling, cell
perturbations, negative controls and batch checks determine how much confidence
the evidence supports. Correlated groups and substitutions are reported alongside
individual feature frequencies.

## Repository and archive structure

The table distinguishes tracked source from artifacts prepared in the original
release archive. A Git clone does not include ignored distributions, full evidence,
or the historical Git bundle; inspect any separately obtained archive inventory.

| Path | Purpose |
|---|---|
| `src/sigtrellis/` | 28 typed production modules: adapters, cell moments/identity, preprocessing, models, splits, stability, panels, DE, external prediction and reporting |
| `tests/` | Unit, integration, synthetic and adversarial controls; exact coupling benchmark recorder |
| `examples/` | Generic configurations, pinned public preparation and scientific benchmark commands |
| `docs/SINGLE_CELL.md` | Exact feature definitions, units, assumptions and missing-value behavior |
| `docs/REVIEW_0.3.0.md` | Five verified previous-release defects, fixes, new capabilities and remaining risks |
| `docs/validation/v0.3.0/` | Compact current execution logs, metrics, dependencies, coverage and provenance checks |
| `README.md`, `RESEARCH.md`, `ARCHITECTURE.md` | Installation, methodological sources and implementation boundaries |
| `VALIDATION.md`, `LIMITATIONS.md` | Measured results and interpretation limits |
| `.github/workflows/`, `Dockerfile`, `pyproject.toml` | Local packaging and configured CI; remote CI/Docker execution is not claimed |
| `LICENSE`, `CITATION.cff`, `CHANGELOG.md` | MIT source license, attribution and release history |
| `release_evidence/` | Full generated reports and numerical/audit artifacts, included in the archive but ignored by Git |
| `dist/` | Built 0.3.0 wheel and source distribution |
| `sigtrellis-history.bundle` | Historical archive copy of local commit/tag history |
| `RELEASE_INVENTORY.json`, `SHA256SUMS.txt` | Archive inventory and per-file hashes |

Large raw data, environments and caches are excluded. Earlier reports retain their
original version labels. Source commits include compact evidence, not raw studies
or large generated results. The archive can be unpacked directly, or its Git
bundle can be cloned to preserve history before copying desired evidence files.

## Validation boundary

The adapters perform only fixed, specimen-local transformations. Inner training
folds own candidate screening, normalization references where applicable,
imputation, scaling and model tuning. Outer folds hold out entire biological
groups. Compact-panel selection repeats inside those boundaries and has its own
permutation/batch evaluation. Frozen external models never reselect features or
refit transformations. The native adapter checks both the predictor dictionary
and full RNA normalization universe.

Count DE uses optional PyDESeq2 with explicit design/rank checks. It is either a
training-fold screen for raw-count models or descriptive supporting evidence.
Native distribution features are not interpreted as DESeq2 counts. Observed
within-cell correlations/variances are not technical-noise-deconvolved statistics.

## Reproduce

From the unpacked source directory, with Python 3.12:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -c docs/validation/v0.3.0/runtime_constraints.txt '.[dev,de,demo]'
python -m pytest
python -m ruff check .
python -m ruff format --check .
python -m mypy src
python tests/record_cell_benchmarks.py --output results/cell_benchmarks.json
python examples/validate_science.py --output results/scientific_benchmarks.json

python examples/prepare_replicated.py --dataset yeast --output data/public/yeast
sigtrellis bulk --expression data/public/yeast/train/counts.tsv.gz \
  --metadata data/public/yeast/train/metadata.tsv \
  --config examples/replicated_bulk.yaml --output results/bulk
sigtrellis external --run results/bulk --model-kind panel \
  --expression data/public/yeast/external/counts.tsv.gz \
  --metadata data/public/yeast/external/metadata.tsv --output results/bulk_panel_external

python examples/prepare_replicated.py --dataset lupus --download \
  --cells-per-sample 0 --output data/public/lupus_all_cells
sigtrellis single-cell --input data/public/lupus_all_cells/train.h5ad \
  --config examples/replicated_single_cell_native.yaml --output results/cell_native
sigtrellis external --run results/cell_native --model-kind panel \
  --input data/public/lupus_all_cells/external.h5ad --output results/cell_panel_external
```

Omit `--model-kind panel` and use a new output directory to evaluate a full model.
Use `--source /path/to/source.h5ad` instead of `--download` to reuse the pinned
12.2 GB single-cell source. Both public workflows need substantial CPU time for
nested selection and permutations. Do not choose scientific settings by repeatedly
optimizing performance on the reserved cohort. See `docs/PUBLIC_DEMOS.md` for
source terms, specimen selection and complete commands.

## Evidence and publication recommendation

See `VALIDATION.md` and its machine-readable summary for current tests, all public
metrics, gate failures and exact execution scope. The public single-cell rerun
uses 526,899 training cells and 375,261 reserved cells from 137 and 96 donors.
It is development validation on a previously examined same-study holdout, not
independent prospective confirmation. No universal improvement over pseudobulk
is claimed from these demonstrations.

Publication as transparent research software is appropriate when accompanied by
these evidence records and limitations. Claims of clinical qualification,
causality, guaranteed transfer across arbitrary systems, or formal gene-level
stability error control are not supported. A valid input study can correctly yield
no candidate that passes all robustness gates.

Remaining research questions include technical-noise-aware cell variability,
conditional inference with nuisance covariates, learned cell encoders with frozen
external assignment, comparison with correctly nested neural multiple-instance
models, independent cross-laboratory transfer, calibration under prevalence shift,
and targeted-assay measurement/normalization. These are distinct research tasks;
the current fixed distribution embedding does not claim to implement them.
