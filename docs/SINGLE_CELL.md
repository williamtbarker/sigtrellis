# Single-cell models without pseudoreplication

The prediction target belongs to a biological specimen or donor. The package
therefore produces one feature row per specimen and keeps related specimens in
one validation group. This holds in both supported representations.

Pseudobulk is the baseline for sample-level gene expression and negative-binomial
count DE. It is not an inferior substitute merely because it reduces cell count:
the independent replication in a phenotype comparison comes from biological units.
Research comparisons and assumptions are in [RESEARCH.md](../RESEARCH.md).

Distribution mode retains other interpretable information from the cells:

| Kind | Exact definition per specimen/state | Interpretation limit |
|---|---|---|
| `abundance` | `log((p_state + 0.001)/(1 - p_state + 0.001))`, where `p_state = n_state/n_all_QC` | Relative recovery; compositional, not absolute tissue counts |
| `gene_mean` | Mean of cell-local `log1p(count/library * 10000)` | Mean transformed expression, not pseudobulk log fold change |
| `gene_detection` | Fraction of state cells with a positive count | Sensitive to depth and dropout |
| `gene_variance` | Unbiased sample variance of transformed expression | Includes measurement noise; not deconvolved biological variance |
| `program_mean` | Mean of a fixed equal-weight program score | Defined gene set, not an inferred pathway mechanism |
| `program_variance` | Sample variance of that score | Biological and technical variability combined |
| `program_q90` | Empirical 90th percentile | Sampling sensitive, especially with few cells |
| `program_fraction` | Fraction of scores above a supplied fixed threshold | Threshold must be chosen independently of held-out outcomes |

A cell's program score is the mean transformed expression of the explicitly
listed genes. All program genes must exist; the adapter refuses silent identifier
substitution. Gene names, states and programs are user-defined and species agnostic.
The human immune programs in one example are **demonstration hypotheses**, not
hard-coded biology or a general-purpose pathway database.

The fixed abundance offset acts on proportions. Duplicating the same recovered
cell distribution leaves this feature unchanged. Absent/pure states are constant
regardless of total cell number; a count pseudocount would incorrectly encode
capture depth at these boundaries. Finite-cell sampling noise still depends on
cell number and is not removed by this transformation.

`mitochondrial_prefix` matches matrix gene identifiers. A symbol prefix such as
`MT-` does not identify mitochondrial genes in an Ensembl-ID matrix. Use appropriate
identifiers or perform annotated mitochondrial QC upstream; no organism annotation
is downloaded or silently inferred by the modeling core.

Declare `single_cell_mode: distribution`, `input_scale: features`,
`normalization: none`, `imputation: median`, and `feature_blocks` in YAML/TOML.
If using `cell_type`, declare its `cell_states` before fitting. `feature_genes`
optionally limits gene-level features; normalization still uses each cell's full
RNA library. Without gene blocks or supporting DE, gene accumulators are not
allocated. Programs-only runs can therefore read large sparse inputs in chunks.
Accumulator/output estimates respect `max_dense_mb`; sparse chunks and Python
object overhead are additional memory. This is not an atlas-scale benchmark.

States below `min_cells` have missing expression/program measurements. A missing
state is not assigned zero expression. Abundance remains measurable when a state
is absent. Imputation medians, feature availability filters and scaling are fitted
only on training specimens; missingness may nevertheless reflect collection bias.
Inspect cell numbers in the QC report. Detection, variance and proportions can
be driven by technical differences even when implemented without leakage.

`cell_resamples` and `cell_fraction` generate deterministic Bernoulli subsamples
within specimens. They change the feature estimates, not the number of donors.
Stability cycles through baseline and perturbed matrices while subsampling whole
biological groups and retuning. Frequencies are saved separately per perturbation;
positive feature gates require the declared frequency threshold in each. More
resamples improve descriptive resolution but do not create independent trials.

`supporting_de: true` performs separate count DE on eligible specimen/state sums.
It retains the full supplied RNA gene universe for count normalization even when
`feature_genes` restricts the predictive dictionary. This can require substantially
more accumulator memory; the memory guard fails explicitly rather than changing
the count-model reference silently.
It measures gene mean-count effects, **not differential variance or differential
abundance**. All state/contrast tests share a BH family. It remains same-cohort
supporting evidence and never changes predictive inputs. Distribution predictors
can use fold-local association screening; mixed feature types cannot enter DESeq2.

Known alternatives include donor-aware negative-binomial mixed models, Memento,
Milo neighborhood abundance, scCODA composition models and hierarchical
multiple-instance models. These answer partly different questions. SigTrellis does
not label its simple summaries as those methods and does not claim their formal
inference. Adding a learned embedding or neighborhood backend would require
training-only fitting and frozen assignment of held-out cells, then independent
comparative validation. Current typed features are a transparent alternative for
elastic-net phenotype prediction, not a claim to supersede these methods.
