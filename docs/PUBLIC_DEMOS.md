# Reproduce the public demonstrations

Run from the repository after `python -m pip install '.[de,demo]'`. Use fresh output
directories. No public raw data are committed. SHA-256 source checksums and fixed
random seeds make the preparation explicit. The distributions retain their own
reuse terms; see [DATA_SOURCES.md](DATA_SOURCES.md).

## Replicated nonhuman cultures

```bash
python examples/prepare_replicated.py --dataset yeast --output data/public/yeast
sigtrellis bulk --expression data/public/yeast/train/counts.tsv.gz \
  --metadata data/public/yeast/train/metadata.tsv \
  --config examples/replicated_bulk.yaml --output results/replicated_bulk
sigtrellis external --run results/replicated_bulk \
  --expression data/public/yeast/external/counts.tsv.gz \
  --metadata data/public/yeast/external/metadata.tsv \
  --output results/replicated_bulk_external
sigtrellis external --run results/replicated_bulk --model-kind panel \
  --expression data/public/yeast/external/counts.tsv.gz \
  --metadata data/public/yeast/external/metadata.tsv \
  --output results/replicated_bulk_panel_external
```

Schurch et al. provide 48 wild-type and 48 SNF2-knockout cultures. Published QC
excludes ten, leaving 86. Seven sequencing lanes per culture have already been
combined; lanes are not extra biological replicates. Before model fitting, seed
2026 reserves ten cultures per genotype. The other 66 train the models.
This tests genotype-associated expression within one experiment, not universal
biomarker transfer or cross-laboratory validation. Batch metadata sufficient for
cross-laboratory testing are unavailable; they are not invented.

## Donor-level single-cell cohort holdout

The source download is **12.2 GB**, from a 261-donor, 1,263,676-cell study. Keep
at least 15 GB free for the source and prepared subsets. The preparation retains
up to 400 randomly selected cells per chosen specimen, all 30,172 source genes,
and original raw counts. It produces approximately 165 MB of prepared H5AD files.
This validates sparse preparation; the model demonstration is not a million-cell
runtime benchmark.

```bash
python examples/prepare_replicated.py --dataset lupus --download --output data/public/lupus
```

If already downloaded, use `--source /path/to/source.h5ad` instead of `--download`.
The complete source SHA-256 must match in either case. Specimen/processing-cohort
aliquots are distinct. Cohort 4 supplies 96 held-out donors. **All measurements of
those donors are removed from training**; one specimen per remaining donor in
cohorts 2/3 gives 137 training donors. Cohort 1 contains controls only and was
excluded by the predeclared case-control design. This is a processing-cohort
holdout within the same observational study, not an independent clinical study.

```bash
sigtrellis single-cell --input data/public/lupus/train.h5ad \
  --config examples/replicated_single_cell.yaml --output results/replicated_single_cell
sigtrellis external --run results/replicated_single_cell \
  --input data/public/lupus/external.h5ad --output results/replicated_single_cell_external
sigtrellis external --run results/replicated_single_cell --model-kind panel \
  --input data/public/lupus/external.h5ad --output results/replicated_single_cell_panel_external

sigtrellis single-cell --input data/public/lupus/train.h5ad \
  --config examples/replicated_cell_distributions.yaml --output results/replicated_cell_distributions
sigtrellis external --run results/replicated_cell_distributions \
  --input data/public/lupus/external.h5ad --output results/replicated_cell_distributions_external
sigtrellis external --run results/replicated_cell_distributions --model-kind panel \
  --input data/public/lupus/external.h5ad --output results/replicated_cell_distributions_panel_external
```

The first analysis prespecifies classical monocytes and discovers gene predictors.
Low-cell-count specimens are explicitly excluded according to its minimum of 20.
The second prespecifies five common cell states and two transparent immune-response
programs, retaining abundance, mean, variability and upper-tail scores. It asks a
different feature-level question; its programs were fixed before external
performance was inspected. Neither representation is chosen retrospectively as a
winner. Disease, treatment and other observational factors may be entangled.

The reserved donors never enter fitting or tuning. During adversarial review, a
synthetic capture-depth counterexample motivated a correction to abundance
smoothing after the first external evaluation. The corrected representation is
rerun on the same reserved cohort, with both versions preserved in release
evidence. Its final holdout numbers are development verification, not a pristine
preregistered confirmation. A genuinely new cohort is needed for that next claim.

Nested prediction, batch holdouts, panel performance and external prediction are
all reported, including failures and empty panels. Nineteen permutations provide
only 0.05 p-value resolution. The small grids and 8/9 stability resamples keep these
demonstrations practical; definitive research needs finer uncertainty assessment,
larger stability budgets and an independently justified analysis plan.

The smaller airway and Kang paired examples remain reproducible using
`examples/prepare_public.py` and the commands in the README. Historical 0.1.1
reports are retained separately and never relabeled as current-release executions.
