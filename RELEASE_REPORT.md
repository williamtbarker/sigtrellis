# SigTrellis 0.1.0 — release assessment

An original, local Python research package for stable candidate transcriptomic signatures. All required core workflows are implemented and executed. No previous employer's code, data, or confidential method is represented as part of this project.

## Repository contents

| Path | Contents |
|---|---|
| `README.md` | Installation, quickstarts, scientific interpretation and output contract |
| `RESEARCH.md` | Primary-method references, design decisions, Stabilomics review and name check |
| `ARCHITECTURE.md` | Modality boundaries, module contracts and extension rules |
| `VALIDATION.md` | Executed checks, numerical results and exact reproduction commands |
| `LIMITATIONS.md` | Supported scope and assumptions that can invalidate conclusions |
| `CHANGELOG.md`, `LICENSE`, `CITATION.cff` | Version history, MIT software terms and citation metadata |
| `pyproject.toml`, `Dockerfile` | Installation, typed/linted/tested package settings and local container recipe |
| `src/sigtrellis/` | Twenty Python modules, CLI, modality adapters, modeling, audits and reporting |
| `tests/` | 89 unit, scientific, integration and adversarial cases |
| `examples/` | Public download/preparation script, scientific benchmark, reproducible study YAMLs |
| `docs/` | Configuration, dataset terms, adversarial review and machine-readable validation evidence |
| `.github/workflows/` | Python 3.12/3.13 checks and wheel-install smoke workflow |
| `dist/` | Built wheel and source distribution |
| `release_evidence/` | Three complete executed reports/audits and public preparation provenance; omitted from Git/source distributions |

## Architecture and statistical methods

Bulk matrices are aligned to metadata; H5AD raw counts are read in chunks and summed within sample/cell type. A shared sample-level engine keeps complete biological groups together through nested CV, elastic-net tuning and retuned stability subsampling. Learned preprocessing and supervised DE screening are fitted separately inside each training boundary.

Binary and multinomial models use logistic elastic net; continuous phenotypes use least-squares elastic net. Both mean-loss lambda and L1 mixing are tuned. Optional PyDESeq2 provides identifiable binary raw-count contrasts with covariates, batch and paired fixed effects. Full-cohort DE is labeled exploratory and cannot feed evaluated prediction.

Feature evidence combines coefficients, signs, selection/rank variation, normalization perturbations, outer-fold stability, batch holdouts and descriptive correlated-feature substitution. Whole-procedure group-aware permutations and baseline comparisons prevent stable-looking noise from automatically becoming a confident panel. Gates remain empirical research criteria, without a gene-level false-discovery guarantee.

## Executed results

- 89 tests pass in both development and fresh installed-wheel environments; 88% branch-aware coverage. Ruff and strict mypy pass.
- The documented 80-sample synthetic quickstart recovers exactly four planted genes as gated candidates.
- Public airway bulk completes with eight samples/four paired donors, count DE and reports; no candidate passes all gates.
- Public Kang single-cell analysis completes from 24,673 input cells, with 5,697 monocytes aggregated into 16 samples/eight donors; no candidate passes all gates.
- Both public examples have AUC 1.0, but the bulk permutation result is 0.25 and the single-cell processing library is perfectly confounded with condition. Reports preserve these reasons to withhold confidence.
- Independent synthetic noise, batch-specific and batch-confounded controls, correlated substitution, imbalance, outliers, paired cells, multiclass and regression checks behave as documented.

Full metrics, figures, fit/split trails, source/input/output hashes, configurations and versions are included. Remote CI and the Docker recipe are supplied but were not executed here.

## Reproduce the public demonstrations

```bash
python -m pip install '.[de,demo]'
python examples/prepare_public.py --dataset all --output data/public
sigtrellis bulk --expression data/public/bulk/counts.tsv \
  --metadata data/public/bulk/metadata.tsv --config examples/public_bulk.yaml \
  --output results/public_bulk
sigtrellis single-cell --input data/public/single_cell/kang.h5ad \
  --config examples/public_single_cell.yaml --output results/public_single_cell
```

Use Python 3.12+ in a virtual environment and fresh output directories. `VALIDATION.md` includes the exact optional dependency constraints used in the fresh verification. Raw public downloads are external; source-specific attribution/reuse terms accompany the derived evidence.

## Limitations, research questions and publication recommendation

Nested procedure performance does not validate the final compact panel. Selection is predictive association, not biological validation, causality, clinical utility or qualification. Public studies are small; perfect batch confounding cannot be resolved computationally. Unseen provenance problems and near-duplicates remain possible. Covariates are jointly penalized; conditional gene inference is deferred. General mixed models, temporal validation, pathway scoring and formal gene-level error control are not implemented.

Future research should address conditional stability with nuisance factors, panel-specific nested validation, larger independently collected cohorts, agreement with R DE reference workflows, hierarchical designs, biological gene-program stability, and prospective assay transfer. Atlas-scale computational claims require separate benchmarking.

**Recommendation:** publish publicly as a tested, clearly scoped research software release with these limitations and the adversarial review attached. It is not ready to be represented as a clinically validated biomarker discovery/qualification platform. The repository and package name have not been published or reserved by this work.
