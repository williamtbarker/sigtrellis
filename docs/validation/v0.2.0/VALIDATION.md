# Executed validation — SigTrellis 0.2.0

Executed locally on Linux / Python 3.12.14, 2026-09-26. These are measured results.
This release adds richer single-cell representations, separately validated compact
panels, native single-cell external evaluation, temporal splitting, expanded count
DE and larger public demonstrations. Historical 0.1.1 results remain in
`docs/validation/v0.1.1/VALIDATION.md`; they are not relabeled as 0.2.0 executions.

## Software verification

| Check | Observed result | Evidence beneath `docs/validation/v0.2.0/` |
|---|---|---|
| Full suite against installed wheel, outside source tree | **152 passed**, 197.40 seconds | `clean_pytest.txt` |
| Branch-aware coverage | **87.10%** overall; 90.15% statements, 78.93% branches | `coverage.json` |
| Ruff lint and formatting | Pass, 58 Python files | `ruff.txt` |
| Strict mypy | Pass, 25 source modules | `mypy.txt` |
| Wheel/source build, clean installation, dependency consistency | Pass | Build/install/dependency logs |
| Final installed source identity | All 25 Python modules match delivered source | `provenance_verification.json` |
| Installed synthetic quickstart with compact-panel evaluation | Four planted genes recovered; four pass gates | `validation_summary.json` and archived toy report |
| Public analysis and frozen full/panel evaluation | Three analysis configurations, six external evaluations complete | Tables below and archived reports |
| Report inspection | Coefficient/stability, path and DE figures inspected; readable labels/units | Self-contained HTML and PNG artifacts |
| Remote GitHub Actions / Docker execution | **Not executed here**; configurations supplied | `.github/workflows/ci.yml`, `Dockerfile` in repository |

The one pytest warning comes from constructing an AnnData fixture with initially
duplicated observation names before assigning unique IDs. It is not a warning from
a public analysis. Exact verified dependency versions are in `runtime_constraints.txt`
in the evidence directory above. They include scikit-learn 1.9.1, PyDESeq2 0.5.4 and
AnnData 0.12.19. CI is configured for Python 3.12/3.13; only the locally executed
environment is claimed as verified here.

## Public predictive results

All models tune lambda and L1 mixing inside nested biological-group CV. Panel
discovery, including stability selection, repeats independently in each outer
training fold. Full-model and panel scores evaluate different training procedures.
Final panels are frozen before external prediction.

| Representation | Training biological groups | Input features | Full nested AUC | Panel nested AUC | Final nonzero panel features | Full-model features passing all gates |
|---|---:|---:|---:|---:|---:|---:|
| Yeast bulk cultures | 66 | 7,126 genes | 1.000 | 1.000 | 10 | 146 |
| Lupus classical-monocyte pseudobulk | 136 | 30,172 genes | 0.887 | 0.882 | 3 | 1 |
| Lupus cell-state/program distributions, corrected | 137 | 35 features | 0.945 | 0.943 | 5 | **0** |
| Synthetic quickstart, seed 7 | 80 | 200 genes | 1.000 | 1.000 | 4 | 4 |

All four full and panel permutation p-values are 0.05, the minimum resolvable with
19 permutations. This coarse demonstration setting is not high-precision significance
estimation. Public stability uses 8/9 subsamples and a small grid; follow-up research
needs a larger prespecified budget.

| Frozen evaluation | Held-out units | AUC | Average precision | Log loss | Brier | Sensitivity | Specificity |
|---|---:|---:|---:|---:|---:|---:|---:|
| Yeast full model | 20 cultures | 1.000 | 1.000 | 0.01085 | 0.00012 | 1.000 | 1.000 |
| Yeast 10-gene panel | 20 cultures | 1.000 | 1.000 | 0.01573 | 0.00025 | 1.000 | 1.000 |
| Monocyte full model | 96 donors | 0.833 | 0.884 | 0.65150 | 0.22929 | 0.962 | 0.273 |
| Monocyte 3-gene panel | 96 donors | 0.747 | 0.803 | 0.71947 | 0.25376 | 0.923 | **0.227** |
| Cell-distribution full model | 96 donors | 0.788 | 0.854 | 0.59131 | 0.20423 | 0.788 | 0.523 |
| Cell-distribution 5-feature panel | 96 donors | 0.804 | 0.863 | 0.58584 | 0.20190 | 0.827 | **0.500** |

Threshold-dependent metrics use the predeclared threshold 0.5. Calibration bins,
balanced accuracy and MCC are retained in each report. The compact machine-readable
`docs/validation/v0.2.0/validation_summary.json` preserves exact values.

**Yeast:** Schurch et al.'s 96 independent cultures become 86 after published QC
exclusions; sequencing lanes are combined within culture. Seed 2026 reserves ten
cultures per genotype before fitting. Both models correctly classified all 20.
Nine of the ten panel genes pass their own panel gates; one fails outer-fold
frequency. This is a strong genotype contrast within one experiment, with no
cross-laboratory batch evidence. The 146 full-model candidates are not 146
biologically validated biomarkers.

**Single-cell source:** the pinned Perez et al. deposit contains 1,263,676 cells
from 261 donors. Preparation retains up to 400 cells per specimen using raw counts
and all source genes. Cohort 4 supplies 96 held-out donors. Every other aliquot of
those donors is purged from training. One aliquot per remaining eligible donor gives
137 training donors; the monocyte minimum-cell rule excludes one. Prepared H5ADs
contain 54,800 and 38,400 cells. This tests sparse preparation, not million-cell
model runtime. Both representations and their fixed programs were declared before
external comparison; neither is retrospectively selected as a winner.

The monocyte analysis has one full-model candidate passing internal gates:
ENSG00000135114 (OASL). Its final panel retains three genes; one passes panel gates.
External panel specificity is only 22.7%, and probabilities transfer poorly.
Internal stability does not establish a clinically useful marker.

The richer representation uses five fixed states and two fixed programs, retaining
abundance, mean, variance and upper-tail features. Both full and panel procedures
lose to their baselines in training-batch holdouts. **No feature passes all gates**,
despite favorable random donor CV. The failed batch checks remain prominent.

The initial distribution external results were inspected before a synthetic
capture-depth counterexample exposed an abundance-smoothing defect. The correction
uses a fixed proportion offset instead of a cell-count pseudocount. No reserved
donor entered fitting/tuning, and programs/grids stayed fixed, but the same holdout
was evaluated again. Corrected results are **development verification, not pristine
preregistered confirmation**. Both versions are retained. An independent new cohort
is required for a stronger generalization claim.

## Scientific and adversarial tests

The final passing suite executes these original and expanded controls:

- Strong planted bulk genes recur under group subsampling. Multinomial and continuous
  planted outcomes improve on baselines. Imbalance, depth outliers and correlated
  signal have explicit fixtures.
- Three independent pure-noise seeds cannot bypass global gates even with hypothetical
  perfect feature stability. A separate compact-panel null repeats the entire
  selection/refit policy in permutations; a single-cell null tests richer summaries.
  These finite controls do not prove a zero false-positive rate.
- An unsafe global screen of 2,000 noise features is a positive leakage control.
  Safe screening and real/simulated DE fit inside training boundaries, verified
  through recorded sample IDs. Normalization, scaling and imputation are tested
  against changes in held-out data.
- Donor overlap, duplicated bulk profiles, target copies and identifier predictors
  fail. Batch-only and batch-specific signals cannot claim cross-batch robustness.
  Public preparation purges every aliquot of an external donor.
- Planted abundance-only and heterogeneity-only single-cell signals are recovered.
  Features do not depend on other donors' data or outcome labels. Missing state
  expression remains unavailable, not zero.
- Repeating identical cell distributions at different capture depths leaves
  abundance unchanged. When only recovered cell count encodes phenotype, abundance
  prediction stays at AUC 0.5 with no selected coefficients.
- Configuration rejects six direct outcome/batch/ID/time encodings as cell-state
  annotations. Aliased or externally learned encodings still need provenance review.
- Continuous and multiclass count-DE tests preserve effect units and contrast
  families. Restricting predictor genes cannot shrink the supplied RNA count-DE
  normalization universe.
- A panel cannot inherit confidence from a full-model weight when its own weight
  is zero, reversed or unstable. Panel outer coefficients are exported.
- Group-purged temporal inner/outer splits respect time order; unrestricted temporal
  permutation inference cannot produce positive gates.
- Frozen binary, multinomial and regression predictions match fitted estimators.
  Native single-cell external evaluation reuses fixed features/transforms and
  rejects biological overlap. Sparse import preserves exact IDs and numeric types.

The substitution control alternates identical-expression genes across fits:
individual frequencies are 0.5 but group any-member frequency is 1.0. Empty selected
sets do not acquire perfect Jaccard stability. `docs/REVIEW_0.2.0.md` records findings,
fixes and remaining assumptions.

## Evidence identity and reproduction

All **150 recorded output hashes** across the four current report directories were
verified, along with prepared inputs and external-model artifact hashes. Every
recorded execution source hash resolves to final source or an exact archived module
under `release_evidence/v0.2.0/execution_sources/`.

Gene-level public runs began before final typing/input-contract guards. The corrected
distribution run preceded the supporting-DE-universe safeguard, which is inactive
in that program-only analysis. These later changes do not alter those valid numerical
paths. Original manifests were not rewritten to pretend all runs used final bytes.
The final installed wheel passes all 152 tests and matches all delivered modules.

| Purpose | Directory beneath `release_evidence/v0.2.0/` |
|---|---|
| Bulk training | `replicated_bulk/` |
| Bulk external full/panel | `replicated_bulk_external/`, `replicated_bulk_panel_external/` |
| Monocyte training | `replicated_single_cell/` |
| Monocyte external full/panel | `replicated_single_cell_external/`, `replicated_single_cell_panel_external/` |
| Corrected distribution training | `replicated_cell_distributions_corrected/` |
| Corrected external full/panel | `replicated_cell_distributions_corrected_full_external/`, `replicated_cell_distributions_corrected_panel_external/` |
| Initial distribution development records | `development_distribution_v1/` |
| Installed quickstart | `toy/` |

Exact download/preparation and analysis commands are in
[docs/PUBLIC_DEMOS.md](docs/PUBLIC_DEMOS.md). [RELEASE_GUIDE.md](RELEASE_GUIDE.md) adds
installation and quickstart commands. Source data retain their own terms in
`docs/DATA_SOURCES.md`; MIT does not relicense data or derived reports.

## Publication assessment

Publication as **general-purpose research software for candidate transcriptomic
signature discovery** is justified. Both modalities, all supported outcome families,
nested/grouped validation, compact panels, stability, negative controls, public
preparation, frozen evaluation and reports have executed. No method can guarantee
a biomarker for every cell line or repair an unidentifiable experiment.
Biological validation, causal claims, clinical qualification, independent
cross-laboratory transfer and targeted-assay validation remain separate work.

See `docs/COMPLETION_MATRIX.md` for feature evidence, `LIMITATIONS.md` for the
scientific contract and `docs/PUBLISHING.md` for research priorities.
