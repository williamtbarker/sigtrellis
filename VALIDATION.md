# Executed validation — SigTrellis 0.1.1

Executed locally on Linux / Python 3.12.14, 2026-09-26. Results below are observed execution results, not anticipated CI outcomes. Source hashes in all final demonstration manifests match the delivered Python modules. Every recorded output hash was checked after each completed run.

## Software verification

| Check | Result | Evidence |
|---|---|---|
| Full development-environment pytest suite | **116 passed**, 59.61 s | `docs/validation/pytest.log` |
| Full suite against installed wheel in a fresh venv | **116 passed**, 62.53 s | `docs/validation/clean_pytest.log` |
| Branch-aware coverage | **89%**; 1,838 statements, 650 branches | `docs/validation/coverage.json` and test logs |
| Ruff lint and format | Pass | `docs/validation/ruff.log` |
| Strict mypy | Pass, 21 source files | `docs/validation/mypy.log` |
| Build and clean wheel installation | Pass | build/install logs under `docs/validation/` |
| Installed dependency consistency | Pass | `docs/validation/pip_check.log` |
| Documented synthetic quickstart, outside source tree | Complete; four planted candidates, no extras | `release_evidence/toy/` |
| Public bulk and single-cell preparation | Complete; checksums and H5AD count round-trip verified | preparation log and JSON under `release_evidence/preparation/` |
| Public end-to-end runs from installed wheel | Both complete | `release_evidence/public_bulk/`, `release_evidence/public_single_cell/` |
| Report visual inspection | Volcano, coefficient IQR and correlation figures inspected; readable axes/labels | PNGs and self-contained HTML in evidence directories |
| Remote GitHub Actions / Docker execution | **Not executed here**; configuration supplied | `.github/workflows/ci.yml`, `Dockerfile` |

The fresh environment independently resolved newer versions: NumPy 2.5.3, pandas 2.3.3, SciPy 1.18.1, scikit-learn 1.9.1, AnnData 0.12.19, h5py 3.16.0, matplotlib 3.11.2, and PyDESeq2 0.5.4. The development environment used scikit-learn 1.8.0 and a different NumPy/SciPy/matplotlib combination. Both passed the same suite. Exact fresh-environment pins are in `docs/validation/runtime_constraints.txt`; this is evidence for these combinations, not a guarantee over every allowed dependency version.

## Public demonstrations and exact interpretation

| Run | Biological units | Input features | Nested ROC-AUC | Log loss | Brier | Permutation p | Genes passing all gates |
|---|---|---:|---:|---:|---:|---:|---:|
| Synthetic quickstart, seed 7 | 80 samples / 80 groups | 200 | 1.000 | 0.012016 | 0.000197 | 0.05 | 4 |
| Public airway bulk | 8 samples / **4 paired groups** | 64,102 | 1.000 | 0.032057 | 0.002020 | 0.25 | **0** |
| Public Kang monocyte pseudobulk | 16 samples / **8 donors** | 15,706 | 1.000 | 0.016379 | 0.000451 | 0.05 | **0** |

All three runs also had average precision, balanced accuracy, sensitivity, specificity and MCC equal to 1.0. These are deliberately easy signal/protocol examples; they are not evidence of universal accuracy. Nineteen permutations only resolve p-values in increments of 0.05. The small public examples are software demonstrations, not biomarker qualification studies.

The synthetic gated panel was exactly `gene_0000`, `gene_0001`, `gene_0002`, `gene_0003`. The modeling code never reads the separately generated truth file.

The bulk run used repeated 2×2 nested grouped CV, 12 retuned group subsamples and both logCPM/frozen-median-ratio perturbations. It failed both the biological-group-count gate and the permutation gate. Paired PyDESeq2 completed; its dispersion uncertainty warning is preserved. Eleven analytical figures include the exploratory count-DE volcano.

The single-cell input contains **24,673 cells**. The predeclared CD14+ monocyte annotation contributes **5,697 cells**, summed into 16 donor-condition samples. Outer CV has four donor-held-out folds and two inner folds. Stability uses 12 group subsamples and two normalizations. Condition and its known pooled processing library are perfectly confounded. Batch holdouts cannot identify the comparison; library-adjusted supporting DE is rank deficient and reported as unavailable. The run fails group-count, perfect-confounding and cross-batch-evidence gates. Eleven figures document the remaining predictive/descriptive evidence without presenting unavailable DE as a result.

Both examples preserve source-specific reuse terms. See `docs/DATA_SOURCES.md`; MIT software does not relicense source data or included derived reports.

## Controlled synthetic benchmark

`examples/validate_science.py` regenerates nine declared scenarios (80 samples, 200 features, seed 17), with common nested 3×2 validation. These are controlled checks of identifiable behavior, not an exhaustive benchmark. Machine-readable results and their resolved common configuration are in `docs/validation/scientific_benchmarks.json`.

| Scenario | Held-out result | Interpretation |
|---|---|---|
| Strong planted signal | AUC 1.000; log loss 0.024912; permutation p 0.05 | All four planted genes meet raw frequency/sign thresholds. One nonplanted gene also meets those thresholds. |
| Pure noise | AUC 0.471875; log loss 0.739316; permutation p 0.50 | Three noise genes look stable by frequency/sign alone; failed global evidence blocks confident claims. |
| Correlated signal | AUC 1.000; log loss 0.048454 | Predictable redundant signal; group/substitution evidence is reported separately from individual stability. |
| Perfect batch confounding | AUC 1.000; log loss 0.013916 | High discrimination does not rescue perfect confounding; QC and unavailable batch contrasts block confidence. |
| Signal in only one batch | AUC 0.544375; log loss 0.696342 | Generalization weakens; batch holdouts include failures of improvement. |
| Class imbalance | AUC 1.000; log loss 0.021866 | Strong planted signal remains detectable; class-sensitive metrics remain explicit. |
| Library-depth outlier | AUC 1.000; log loss 0.024912 | Depth outlier is flagged. Uniform count multiplication cancels under logCPM; this is not general robustness to arbitrary outliers. |
| Multiclass | Macro OvR AUC 1.000; log loss 0.036025 | Multinomial prediction, contrasts and report generation work. |
| Continuous | RMSE 0.125972; MAE 0.094106; R² 0.983383 | Regression and coefficient units behave as declared. |

Additional automated noise fixtures use independent seeds 3, 29 and 101. In each, even a hypothetical perfectly stable feature cannot bypass the observed failed global evidence gates. This is a finite negative-control battery, not proof of a zero false-positive rate.

A constructed correlated substitution test alternates two identical-expression genes across ten fits: each individual frequency is 0.5 while the group's any-member frequency is 1.0 and exclusive-selection fraction is 1.0. Empty/empty feature sets receive undefined Jaccard similarity, preventing false reassurance from repeatedly selecting nothing.

## Leakage and pseudoreplication attacks

- Every supervised training fit is checked against its recorded training/test sample sets. A simulated DE backend must be called separately inside every training boundary.
- A deliberately unsafe global association screen on 2,000 independent noise features inflates AUC relative to the fold-local procedure; the test fails if the safe path behaves like that leaky positive control.
- Frozen normalization/scaling tests alter held-out distributions without changing training state. IDs containing phenotype text cannot enter the model.
- Donor overlap and duplicated expression profiles across groups fail. External validation also rejects renamed exact training-profile copies.
- Group permutations preserve donor units; paired permutations preserve within-donor labels; stratified permutations preserve declared strata.
- Batch holdouts purge any training specimen sharing a held-out donor.
- H5AD dense/sparse layer aggregation matches manual raw-count sums. Fractional/integrated values and inconsistent sample metadata fail. Many cells from only two donors cannot support nested validation.
- Real PyDESeq2 tests check planted effect directions, BH-adjusted significance and paired design behavior; nonidentifiable/repeated-group misuse fails explicitly.
- Binary, multinomial and continuous frozen numeric predictions match their fitted sklearn estimators. The lambda-to-C convention is tested directly.

## Reproduction

From the unpacked repository, optionally pin the verified dependency combination:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -c docs/validation/runtime_constraints.txt '.[de,demo]'

sigtrellis simulate --output data/toy --samples 80 --features 200 --seed 7
sigtrellis bulk --expression data/toy/counts.tsv --metadata data/toy/metadata.tsv \
  --config data/toy/config.yaml --output results/toy

python examples/prepare_public.py --dataset all --output data/public
sigtrellis bulk --expression data/public/bulk/counts.tsv \
  --metadata data/public/bulk/metadata.tsv --config examples/public_bulk.yaml \
  --output results/public_bulk
sigtrellis single-cell --input data/public/single_cell/kang.h5ad \
  --config examples/public_single_cell.yaml --output results/public_single_cell

python examples/validate_science.py --output results/scientific_benchmarks.json
```

Use fresh output directories. Public demonstration reports in the archive were generated with the same delivered source from the built wheel outside the source tree. The synthetic benchmark table was generated in the recorded fresh installed-wheel environment; small numerical differences under other dependency versions are not hidden. Model-manifest input hashes and source hashes enable a precise comparison.

## Acceptance and remaining work

Both modalities, all three outcome families, elastic-net tuning, empirical stability, nested/grouped validation, negative controls, required reports, reproducibility records, public downloads/preparation, clean install and documented quickstart have executed successfully. The requested scientific limitations are enforced or explicitly documented rather than represented as supported inference.

The release is scientifically credible as **research software with tested safeguards**. It is not a completed clinical qualification, independent reproduction of every DE backend, a validated compact assay panel, a formal feature-level error-control method, or an atlas-scale performance benchmark. Remaining research priorities are conditional inference with covariates, panel-specific nested discovery, larger independent cohorts, hierarchical/random-effect designs, and gene-program-level validation. See `docs/ADVERSARIAL_REVIEW.md` for the final risk review.

## Follow-up adversarial source review

The 0.1.1 source review adds 27 tests. Twenty-six adversarial cases reproduce failures in 0.1.0; a real nested PyDESeq2 integration check also passes on the original as a positive control. All 116 tests pass after correction. `docs/REVIEW_0.1.1.md` records severity, fixes and limits. Original records are preserved under `docs/validation/v0.1.0/`. Both public input expression hashes, primary metrics (within 1e-10) and candidate decisions match the original valid-case demonstrations. Source and output hashes now describe the corrected release.
