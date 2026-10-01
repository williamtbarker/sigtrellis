# Executed validation — SigTrellis 0.3.0

Executed locally on Linux / Python 3.12.14, 2026-09-26. These are measured
software and scientific-control results. The public single-cell workflow uses
native cell distributions and within-cell coupling, without summed-count
pseudobulk. Earlier release evidence remains versioned; 0.2.0 results are in
`docs/validation/v0.2.0/VALIDATION.md`.

## Software verification

| Check | Observed result |
|---|---|
| Final installed-wheel suite, outside source tree | **175 passed**, one fixture warning, 226.95 seconds |
| Branch-aware coverage | 87.66% overall; 90.57% statements; 79.84% branches |
| Ruff lint/formatting and strict mypy | Pass; all 28 production modules type checked |
| Wheel/source build and clean installation | Pass; dependency versions recorded |
| Native and bulk CLI smoke/quickstart | Complete reports from installed package |
| Public bulk/native-cell analysis | Complete, including nested full/panel selection, permutations and batch checks where identifiable |
| Frozen external prediction | Full and panel models evaluated for both public modalities |
| Provenance | All current installed modules match delivered source; run output hashes verified |
| Remote CI / Docker | Configurations supplied; not executed here |

Compact logs, exact metrics, coverage and source/output verification are under
`docs/validation/v0.3.0/`. Full reports are in the archive's
`release_evidence/v0.3.0/`. The single pytest warning arises while constructing a
fixture with temporarily duplicated AnnData observation names; it is not a public
run warning. An initial coverage run was interrupted by exhausted temporary disk
space; after removing incomplete intermediates, the clean installed suite was
rerun successfully. Its interruption log is retained separately.

The large native run exposed repeated per-group sparse product work. Chunk-level
product accumulation reduced a fixed local benchmark's median update time from
3.408 to 0.727 seconds. The final 175-test suite was repeated after this numerical
optimization. Bulk and the earlier smoke executions retain their exact
pre-optimization module snapshot in `execution_sources/`; the changed code concerns
only within-cell products. Their source hashes are resolved explicitly rather
than relabeled. This is not a general end-to-end speed guarantee.

## Public model results

Every outer split holds out biological groups. Hyperparameters, screens and
compact panels are selected only from the corresponding training data. The full
model and compact-panel discovery procedure have separate scores.

| Representation | Training groups | Features | Full nested AUC | Panel nested AUC | Nonzero final panel | Full features passing gates |
|---|---:|---:|---:|---:|---:|---:|
| Yeast bulk | 66 | 7126 | 1.000 | 1.000 | 10 | 146 |
| Native single-cell distributions | 137 | 445 | 0.973 | 0.955 | 8 | 0 |
| Synthetic bulk quickstart | 80 | 200 | 1.000 | 1.000 | 4 | 4 |

| Frozen model | Held-out units | AUC | Average precision | Log loss | Brier | Sensitivity | Specificity |
|---|---:|---:|---:|---:|---:|---:|---:|
| Yeast full | 20 | 1.000 | 1.000 | 0.01085 | 0.00012 | 1.000 | 1.000 |
| Yeast panel | 20 | 1.000 | 1.000 | 0.01573 | 0.00025 | 1.000 | 1.000 |
| Native single-cell full | 96 | 0.801 | 0.866 | 0.56667 | 0.18949 | 0.750 | 0.659 |
| Native single-cell panel | 96 | 0.810 | 0.877 | 0.54720 | 0.18442 | 0.750 | 0.659 |

Threshold-based metrics use the predeclared 0.5 threshold. Reports also retain
balanced accuracy, MCC and calibration bins. Nineteen permutations provide a
minimum p-value of 0.05; this is a coarse demonstration budget. Stability uses
8/9 subsamples and a small grid. Neither budget is offered as a universal
recommendation for definitive biomarker research.

**Yeast:** 86 independent cultures remain after published QC; 66 train and 20 are
reserved. Sequencing lanes are combined within culture. This is a strong genotype
contrast in one experiment, without cross-laboratory batch evidence.

**Single-cell:** source data contain 1,263,676 cells from 261 donors. Preparation
retains all cells from the selected specimens: 526,899 training cells from 137
donors and 375,261 held-out cells from 96 donors, with all 30,172 RNA genes.
Every aliquot of a held-out donor is purged from training. Training uses one
specimen per donor in processing cohorts 2/3; cohort 4 supplies the holdout.
The all-control first cohort is excluded by the declared design. The native
representation uses five states and 445 features, including three tail thresholds,
21 gene pairs and two fixed programs. Raw counts are never summed into pseudobulk
for this predictive analysis.

The holdouts were examined during prior releases. No reserved donor enters model
fitting or tuning, but current reruns are **development validation**, not fresh
preregistered confirmation. Cell counts and feature dictionaries differ from
0.2.0: differences in scores are not a controlled claim that the new method is
superior. Disease, treatment, depth and processing may remain confounded.

Yeast: full/panel permutation p-values are 0.05 / 0.05; 9 final panel features pass panel-specific gates. Full-model gate blockers: `[]`.


Native single-cell: full/panel permutation p-values are 0.05 / 0.05; 0 final panel features pass panel-specific gates. Full-model gate blockers: `['cross_batch_performance_failure', 'compact_panel_cross_batch_performance_failure']`.

- Batch 2.0: tested; full loss improvement -2.375721738729307; panel loss improvement -0.5087516318967151.
- Batch 3.0: tested; full loss improvement 0.1736787395797308; panel loss improvement 0.09960027448298864.

## Controlled scientific evidence

The coupling fixture contains 48 donors and 80 cells per donor. It changes which
genes are expressed together in a cell while keeping the marginal mechanisms
uninformative. Identical donor-level nested CV is applied to all representations.

| Representation | Nested ROC-AUC |
|---|---:|
| count_pseudobulk | 0.385417 |
| marginal_summaries | 0.354167 |
| within_cell_coupling | 1.000000 |

Reproduce this exact table with:

```bash
python tests/record_cell_benchmarks.py --output results/cell_benchmarks.json
python examples/validate_science.py --output results/scientific_benchmarks.json
```

A separate passing test plants four rare high-expression cells among 80 while
keeping the raw gene count total at 400 in every specimen; tail features recover
the phenotype. Direct-array checks verify tail and correlation calculations.
These fixtures show available information, not a universal predictive advantage.

The bulk benchmark recovers all four planted strong genes. It also reports a
stable nonplanted gene, demonstrating why selection frequency is not biological
validation. Its pure-noise run has AUC 0.472 and permutation p=0.5; three features
can look stable by frequency alone, while predictive/global gates prevent a
qualified panel. Separate full-pipeline null tests verify these gates. Perfect
batch confounding blocks eligibility; batch-specific signal fails transfer in
held-out batches. Correlated predictors, imbalance, outliers, multinomial and
continuous outcomes are also exercised. Finite simulations do not establish a
zero false-positive rate.

## Leakage and regression evidence

Regression tests cover cell ordering, unrelated-donor perturbation effects, renamed/reordered duplicates, changed normalization gene universes, and misaligned perturbation rows. The suite additionally checks supervised DE inside training folds, unsafe global screening as a positive control, training-only transforms, donor splits, duplicated bulk profiles, target/identifier encodings, cell-state confounding, paired DE, external artifact integrity, and canonical feature serialization.

The single-cell regression tests are in `tests/test_cell_adversarial.py`. Exact full-cell source hashes, resolved configurations, and model outputs are preserved. Hashes assume a trusted manifest; they are not digital signatures.

## Interpretation and remaining research

Passing software tests establishes that these tested contracts work. Favorable
cross-validation does not establish transportability, biological mechanism,
causality or clinical utility. A run that finds no features passing all gates is
a valid outcome. Observed variance/correlation includes technical and compositional
variation, annotations can encode upstream bias, and fixed summaries do not
preserve every cell interaction. Independent cohorts, technical-noise modeling,
conditional inference and correctly nested learned cell representations remain
research priorities. Publication as transparent research software is justified
within this scope; universal or clinical biomarker claims are not.
