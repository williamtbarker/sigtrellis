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

## Native single-cell distributions, with every selected cell

The source download is **12.2 GB**, from a 261-donor, 1,263,676-cell study. Keep
at least 18 GB free for the source, prepared subsets and run artifacts. The current
preparation retains **526,899 training cells from 137 donors** and **375,261
held-out cells from 96 donors**, with all 30,172 source genes and original raw
counts. `--cells-per-sample 0` disables the demonstration cap. The analysis itself
never imposes a cell cap. Sparse preparation writes chunks in AnnData format.

```bash
python examples/prepare_replicated.py --dataset lupus --download \
  --cells-per-sample 0 --output data/public/lupus_all_cells

sigtrellis single-cell --input data/public/lupus_all_cells/train.h5ad \
  --config examples/replicated_single_cell_native.yaml \
  --output results/single_cell_native
sigtrellis external --run results/single_cell_native \
  --input data/public/lupus_all_cells/external.h5ad \
  --output results/single_cell_native_external
sigtrellis external --run results/single_cell_native --model-kind panel \
  --input data/public/lupus_all_cells/external.h5ad \
  --output results/single_cell_native_panel_external
```

If already downloaded, replace `--download` with `--source /path/to/source.h5ad`.
The complete source SHA-256 must match in either case. The supplied native
configuration explicitly sets `single_cell_mode: distribution`; it does not
construct summed-count pseudobulk and leaves optional count-DE support disabled.
The 445 candidate features describe five fixed cell states, ten fixed program
member genes, three expression-tail thresholds, 21 within-program gene pairs,
and two fixed programs. The feature dictionary uses the same previously declared
program genes; external outcomes were not used to choose additional genes or tune
thresholds. A feature cap is fitted independently inside each training fold.

Specimen/processing-cohort aliquots are distinct. Cohort 4 supplies 96 held-out
donors. **Every aliquot of those donors is purged from training**; one specimen
per remaining donor in cohorts 2/3 gives 137 training donors. Cohort 1 contains
controls only and is excluded by the case-control design. We retain all cells of
these chosen specimens, not every cell of the original study. Donors remain the
independent validation units. Disease, treatment and processing factors can still
be entangled.

These cohorts were examined during earlier releases. Reserved donors never enter
training or tuning, but current reruns are **development validation, not a new
preregistered confirmation**. This processing-cohort holdout is from the same
observational study, not an independent clinical cohort. A new study is needed
for stronger transfer claims. More distribution information does not guarantee
better AUC, calibration, specificity or batch robustness.

Nested prediction, batch holdouts, panel performance and frozen prediction are
all reported, including failures and empty panels. Nineteen permutations provide
only 0.05 p-value resolution. The small grids and 8/9 stability resamples keep these
demonstrations practical; definitive research needs a larger prespecified budget.
See [VALIDATION.md](../VALIDATION.md) for measured current results.

## Historical comparisons and smaller practice runs

To reproduce the previous capped preparation, use `--cells-per-sample 400`
(the preparation script's backward-compatible default) and a different output
path. `examples/replicated_single_cell.yaml` supplies classical-monocyte
pseudobulk; `examples/replicated_cell_distributions.yaml` supplies the earlier
35-feature representation. Their executed 0.2.0 results remain in
`docs/validation/v0.2.0/VALIDATION.md`. They are not relabeled as 0.3.0 runs or
used as a controlled claim of improvement: cell counts and representations differ.

The smaller airway and Kang paired examples remain reproducible using
`examples/prepare_public.py` and the commands in the README. Historical 0.1.1
reports retain their original version labels.
