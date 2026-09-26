# Completion evidence for 0.3.0

This ledger distinguishes executed behavior from biological claims. Measured public
results and limitations are in [VALIDATION.md](../VALIDATION.md).

| Requirement | Executed implementation/evidence |
|---|---|
| General-purpose organisms/cell lines | Opaque gene IDs, explicit phenotypes/culture groups; nonhuman yeast public run; `CELL_LINE_GUIDE.md` and generic configuration |
| Bulk binary, multiclass, continuous | Nested logistic/multinomial/regression elastic net; all families in final passing scientific tests |
| Single-cell expression baseline | Sparse chunked count pseudobulk, donor splitting; current passing adapter tests; historical public monocyte workflow and 96-donor evaluation |
| Single-cell information beyond means | State abundance, detection, variability, rare tails and within-cell gene/program correlations; coupling/rare-cell controls and public 445-feature all-cell run |
| Cell sampling robustness | Separate within-specimen cell perturbations and whole-group resampling; sample-local/chunk/order-invariance/unequal-cell-count tests |
| Compact panel | Outer-training-only selection/refit; separate permutations, batch checks and frozen external states; panel sizes and held-out scores reported in VALIDATION.md |
| Batch, pairing and temporal structure | Purged batch/group holdouts, identifiable paired DE, forward inner/outer group splits; synthetic leakage/confounding controls |
| Learned transformations | Fold-local normalization, filtering, screening, imputation, scaling and encoding; frozen numeric states and adversarial tests |
| Count DE | Executed PyDESeq2 binary/continuous/multiclass tests; paired-design checks; multiplicity and full count normalization universe |
| Correlated predictors | Correlation groups, co-selection and substitution; alternating identical genes retain group stability |
| Null/confounded results | Bulk nulls, compact-panel null, single-cell null, batch-only and batch-specific controls; public distribution candidates blocked by batch failure |
| Reproducibility and reports | Configurations, hashes, splits, fit records, versions, coefficients, evidence tables, HTML/Markdown and figures; all current output hashes verified against execution manifests |
| Public demonstrations | Pinned yeast and Perez preparations, full/panel external results including poor specificity; historical airway/Kang evidence retained |
| Installation and quality | Clean installed-wheel quickstart; **175 passing tests**, Ruff/strict mypy pass; wheel/source distribution built |
| Repository readiness | Research/design/validation/limitations, MIT/CITATION, examples, Docker and CI supplied; local history bundle; no remote push performed |

Cell-level mixed models, neighborhood inference and neural multiple-instance models
were investigated as distinct approaches. They are not implemented backends and
are not claimed as aliases for the sample-summary engine. Censored survival is
outside the supported outcome contract.

Operational gates do not control gene-level FDR or establish a validated assay.
Independent cultures/donors and a varying phenotype remain necessary. Remote CI and
Docker have not been executed here. The public holdout has been inspected during prior development; that interpretation
is explicit in the validation record. No improvement claim is inferred by comparing
representations with different cell counts.
