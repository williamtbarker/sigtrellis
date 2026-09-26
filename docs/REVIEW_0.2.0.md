# Adversarial review of the expanded implementation

Reviewed 2026-09-26. This review evaluates code boundaries and scientific claims;
it is not an independent laboratory replication or external security certification.
Execution counts and public measurements are recorded separately in `VALIDATION.md`.

| Attack or misleading interpretation | Implemented defense and review result |
|---|---|
| Select a stable panel globally, then claim its full-model CV score | `panel.py` runs selection/refit inside each outer training set; panel has separate predictions and permutation/batch results. Boundary spies check training IDs. |
| A gene has a stable full-model weight but a zero/reversed panel weight | Panel-specific nonzero, sign and outer-frequency gates are required; `panel_outer` coefficients are exported. A direct adversarial test covers zero, reversed and unstable panel weights. |
| Cell-rich donors become thousands of independent replicates | Adapters emit specimen rows; splits and subsampling operate on complete biological groups. Within-specimen cell perturbations only alter measurements. |
| A richer cell representation learns from held-out donors | All normalization and moments are cell/specimen-local. State labels, program members and thresholds are fixed inputs. Perturbing another donor or all outcome labels leaves a specimen's features unchanged in tests. |
| Outcome/batch/ID columns masquerade as cell-state annotations | Configuration rejects outcome, sample/group ID, batch, time and permutation stratum as the cell-type column. Six regression cases cover these direct label encodings. Aliased upstream annotations still require provenance review. |
| Missing cell state becomes zero expression | Unavailable expression/program features are NaN; only abundance can represent an absent state. Training-only median imputation and missingness filters are tested independently. |
| Gene/program feature names collide or acquire wrong biological labels | Feature IDs encode kind/state/source as JSON; outputs retain separate source gene, cell type, program and unit. Explicit collision tests and report-schema checks pass. |
| Feature variability is advertised as deconvolved biology | Units and limitations explicitly include measurement noise; supporting DE tests mean counts only. No Memento/NEBULA/Milo/scCODA implementation is claimed. |
| More feature types turn pure noise into a confident panel | Independent null fixtures and full-procedure panel permutations are tested. Noise cannot bypass the operational gates in the executed fixtures; no zero false-positive theorem is claimed. |
| Abundance smoothing converts cell yield into a phenotype marker | The initial count pseudocount depended on capture depth even for absent/pure states. A fixed offset on proportions replaces it. New tests require invariance under duplicated cell distributions and chance prediction when only recovered cell count encodes phenotype. The affected public representation is rerun and its reused holdout is explicitly disclosed. |
| Cohort split includes another aliquot of the same donor | Public preparation exposed reused specimen IDs across processing cohorts. The corrected code defines aliquots first and purges all held-out donors. A dedicated preparation regression test reproduces the trap. |
| Numeric-looking IDs collapse during sparse export import | All IDs/categorical metadata remain exact strings; numeric columns must be explicit. Barcode alignment and counts are checked. |
| H5AD provenance paths break serialization | A test revealed that slash-containing dictionary keys create HDF5 paths. Import provenance now serializes path/hash records as JSON values. Native AnnData round-trip is tested. |
| Continuous DE effect is confused with predictive coefficient scale | DE records per-outcome-unit effects; model coefficients remain per training feature SD. Multiclass plots retain contrast identity and use a joint testing family. |
| A small predictor dictionary silently defines the count-DE normalization reference | Distribution-mode supporting DE preserves all supplied RNA genes while rendering only requested gene predictors. A direct raw-sum test checks both the full count universe and restricted feature dictionary. |
| Forward CV leaks future observations from a returning donor | Group start/end times purge overlapping donors; every inner/outer split is checked. Unrestricted time-series permutation inference is blocked. |
| Batch effects support a compact panel that only predicts one cohort | Full and panel procedures are separately tested across batches; unavailable/failed batch checks block positive gates. |
| External evaluation quietly refits preprocessing or selects features | Native H5AD reuses the stored feature contract; numeric model state is frozen and integrity-checked. Tests compare serialized predictions and reject overlapping groups. |

The original global-screen, target-copy, normalization, duplicate, paired-DE,
multiclass and batch-confounding attacks remain in the suite. Existing source
identifiers never become ordinary predictors. Full-data PCA, correlations, paths
and supporting DE are generated after predictive evaluation and do not feed back
into it. The final report separates association, biological validation and causality.

## Remaining assumptions that matter

- Metadata must correctly identify biological independence. Exact profile checks
  cannot establish the absence of hidden relatives, relabeled near-duplicates,
  or reordered copies of cell collections. Low-dimensional summary ties are not
  sufficient evidence of duplication.
- Fixed imported cell labels may already reflect phenotype-aware integration.
  The package cannot prove upstream annotation independence.
- Cell recovery, depth and dropout can affect abundance/variance features.
  Cell perturbations measure some sampling sensitivity, not all technical bias.
- Covariate-adjusted conditional gene significance and temporal exchangeability
  are not supplied by a global permutation. Such runs remain exploratory.
- Operational feature gates do not control gene-level FDR. Correlated-group
  stability is descriptive; a group may be stable while its representative varies.
- Cross-study and targeted-assay transfer remain empirical questions. This
  release's larger external sets are held-out cultures/cohorts from the same
  source studies, not independent clinical validation.

The code can therefore be published as general-purpose **research software with
an explicit input/replication contract**. A claim that it always produces a valid
biomarker from any arbitrary cell-line dataset would be scientifically unsupported.
