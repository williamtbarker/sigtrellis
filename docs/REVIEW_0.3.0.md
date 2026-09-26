# Adversarial improvement pass — 0.3.0

The goal was to make both bulk and native single-cell analysis defensible and
executable, while retaining cell distribution information that count sums lose.
The review covered adapters, study contracts, preprocessing, grouped nesting,
stability/panel selection, DE, external prediction, reporting and preparation.
It was performed against the delivered 0.2.0 implementation, not a hypothetical
design. Current execution results are in [VALIDATION.md](../VALIDATION.md).

## Verified defects and repairs

Five new counterexamples were executed against the installed 0.2.0 wheel. All
five failed there and pass in 0.3.0. The baseline failure log is preserved in
`docs/validation/v0.3.0/v020_counterexamples.txt`.

| Counterexample | Why it mattered | Repair |
|---|---|---|
| Reorder cells/genes | Perturbations and copied-specimen identity changed with serialization order | Cell-ID-keyed inclusion; canonical gene indices and multiset cell fingerprints |
| Remove an unrelated donor | Sequential random draws changed another donor's cellular perturbations | Inclusion depends only on cell ID, seed and perturbation index |
| Rename and reorder a copied specimen | Exact duplicated cell collections could evade duplicate checks | Multiplicity-preserving raw-cell digest independent of cell/specimen names and row order |
| Change the normalization gene universe but retain selected features | Native external predictions could use a different per-cell denominator | Save and verify the complete RNA gene-universe contract; reject legacy missing contracts |
| Reorder a perturbation matrix without aligning metadata | Training rows could receive another specimen's measurements | Validate exact row/column identity for every perturbed frame |

Additional fixes preserve canonical program-pair names through sorted JSON,
reject ambiguous pairs/thresholds and ignored state-selection arguments, permit
one informative predictor, and mark correlations/variances with insufficient
cellular information as unavailable. Pair features do not inherit a mean-DE
p-value from one of their members.

## New native-cell capabilities and controls

- Fixed high-expression fractions capture rare tails. A controlled fixture keeps
  each specimen's bulk total equal while changing a four-cell subpopulation.
- Declared gene-pair and program-pair correlations preserve pairing within the
  same cells. A fixture changes joint expression while retaining uninformative
  marginal mechanisms; nested native prediction succeeds while marginal and
  pseudobulk baselines do not.
- Direct numerical checks compare streaming statistics with explicit cell-array
  calculations. Constant-member correlations remain missing.
- Chunk-level accumulation replaces repeated per-group sparse products. A fixed
  interleaved-donor microbenchmark decreased median update time from 3.408 to
  0.727 seconds (2,048 cells, 685 possible specimen/state groups, 21 gene pairs;
  three repetitions, seed 27). This is a local microbenchmark, not an end-to-end
  runtime guarantee. The final installed suite checks the same numerical results.
- Means, variances, activation fractions and correlations do not retain cell
  arrays. Exact quantiles explicitly retain scores and enter the memory budget.
- Public all-cell preparation is sparse and chunked. Its AnnData round-trip test
  checks counts, cell identity and gene identity.
- Reports show specimen-level distributions with units and technical cell QC.
  These are post-selection descriptions, not independent evidence or cell-level
  significance tests.
- Filtering audits retain a shared feature universe plus each fit's retained
  list, preserving exclusions without repeatedly serializing every excluded gene.

The existing bulk, paired/grouped, multinomial, continuous, temporal, real-DE,
null, confounding, duplicated-profile, target-copy and fold-boundary tests remain
part of the full installed-wheel suite. The new module adds 23 test cases.

## What this pass does not establish

The representation is a fixed interpretable embedding of each specimen's cells,
followed by elastic net. It is not a neural attention/MIL model, a cell-level
phenotype classifier, deconvolved biological covariance, or differential
abundance inference. Finite distribution features cannot preserve every possible
cell pattern. More cells improve measurement of a specimen; they do not increase
the number of independent donors.

Known risks remain partial cell overlap, unknown donor relationships, annotation
leakage, compositional effects, measurement noise, observational confounding and
overfitting through repeated analyst choices. Exact duplicate checks cannot
resolve all of these. Feature gates are operational filters, not formal gene-level
FDR or biological validation. Previously inspected public holdouts remain
development evidence. Independent cross-laboratory cohorts, technical-noise
models and comparisons with correctly nested learned cell encoders are the next
research steps.
