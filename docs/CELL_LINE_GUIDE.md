# Bring your own cell line

SigTrellis asks which RNA measurements reproducibly predict a phenotype **across
independent biological experiments**. It makes no human, cancer, infection or
species assumption in the modeling core. Gene IDs are opaque strings.

## Define the comparison first

Examples include treatment versus control, resistant versus sensitive cultures,
or a measured continuous productivity/yield. The phenotype must vary. RNA from
one unlabelled culture cannot identify a phenotype-associated biomarker by itself.

For bulk, provide a gene count matrix and one metadata row per sequenced specimen.
For single-cell, provide sparse AnnData with integer RNA counts in `X` or a named
layer, and cell metadata that identify each specimen and its biological group.
Cell number does not replace replication. Technical libraries of the same
material should normally be combined at the count level; separately measured
conditions from the same source culture must share a group.

| Field | Meaning | Example |
|---|---|---|
| `sample_id` | Exact sequenced specimen ID | `experiment_01_control` |
| `group` | Unit that must stay together when testing generalization | `experiment_01` |
| `outcome` | Label or measured numeric response | `treated`, or `12.7` |
| `batch` | Observed processing batch, if known | `sequencing_run_02` |
| `covariates` | Explicit measured nuisance/predictive variables | starting density |
| `cell_type` | Fixed cell-state annotation, when available | `cycling` |

Do not set `group` to the cell-line name when all experiments use the same line:
that creates one group. Conversely, renaming aliquots as independent cultures
does not make them independent. Choose the grouping level to match the intended
generalization: new cultures, new experimental days, or new cell lines. A marker
that predicts new cultures from one line has not been shown to transfer to other
lines, laboratories, organisms, or measurement platforms.

## Start with a declared configuration

`examples/cell_line.yaml` is a generic binary template. Edit the phenotype and
metadata column names. It assumes independent culture-constant labels; for paired
conditions declare the shared group, `permutation_scheme: within_group`, and
`de_pair_group: true` if requesting paired count DE. Within-group exchangeability
still requires a defensible experiment. Do not shuffle treatment order blindly.

```bash
sigtrellis bulk --expression counts.tsv --metadata metadata.tsv \
  --config examples/cell_line.yaml --output results/my_cultures
```

For a continuous outcome use `outcome_type: continuous`, a numeric outcome column,
and remove `positive_class`. For multiple categories use `multiclass`; every
training fold must contain every class. Reduce folds only when the design can
still support the intended question; a smaller fold count does not create power.

The single-cell baseline sums counts per specimen and declared state. Distribution
mode additionally models abundance and within-specimen variation. Without state
annotations, omit `cell_type` to use one `all_cells` population, and select gene
or program features; its abundance is not informative. See [SINGLE_CELL.md](SINGLE_CELL.md).

## Read the evidence, not just the ranked list

1. Check biological group counts, cell exclusions, batch confounding and warnings.
2. Check held-out log loss or regression error against the baseline, together with
   discrimination and calibration. AUC alone can conceal poor probabilities.
3. Read selection frequency, sign consistency, perturbation sensitivity and
   correlated-group substitution. Genes can trade places while group signal persists.
4. Inspect `passes_robustness_gates` and its blockers. A completed run with no
   supported candidates is a valid result. Do not lower gates until a panel appears.
5. If using the compact panel, read its own metrics and `panel.csv`. Its weights
   and outer-fold stability can differ from the full model.
6. Freeze the complete preprocessing/model state and evaluate an untouched culture
   set with `sigtrellis external --model-kind panel` before choosing assay targets.

Selection gives a **predictive candidate**. Orthogonal RNA/protein assays,
independent experiments, functional perturbations and intended-use studies provide
different kinds of evidence. They cannot be inferred from a nonzero coefficient.

The included yeast experiment demonstrates a nonhuman cell-culture workflow. Its
same-study held-out cultures do not establish cross-laboratory transfer. The airway
example demonstrates why four paired cell lines remain too few biological groups
for a confident general-purpose panel even with impressive apparent prediction.
