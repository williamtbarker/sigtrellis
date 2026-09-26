# Research and architectural decisions

Reviewed 2026-09-26. SigTrellis is a new general-purpose implementation of public methods. It is not a reconstruction of an employer's software, data, or methodology. This review guided implementation; it is not a claim that one workflow is optimal for every experiment.

## Biological unit and modality boundary

Bulk input is a sample-by-gene matrix. Single-cell input becomes either **summed raw-count pseudobulk per specimen/type** or typed sample-local abundance/expression-distribution features. A sample denotes one specimen/condition/processing aliquot. A separate grouping column identifies the donor or biological unit shared across samples. Every inner/outer split and stability subsample operates on complete groups; batch holdouts purge overlapping groups.

Squair et al. benchmarked false discoveries caused by neglecting biological replication [1]. The 2025 nested-settings benchmark supports conventional pseudobulk methods for individual studies while identifying different tradeoffs for atlas-level analyses [2]. Seurat's current DE vignette likewise documents sample-level summed counts [3]. These findings support a pragmatic default, not a universal prohibition on cell-level mixed models. Mixed models, hierarchical designs beyond one grouping level, and atlas-level DREAM analyses are deferred rather than approximated incorrectly.

Cell QC uses explicit per-cell thresholds; it does not learn phenotype-specific gates. HDF5 chunks avoid densifying cells × genes. Summed sample-level counts are dense with a declared memory limit. Missing sample/cell-type combinations are not fabricated as zero-expression samples. Cell counts and exclusions are preserved. Doublet removal and annotation quality remain upstream responsibilities.

## Differential expression

DESeq2 models counts with negative-binomial GLMs, estimates dispersion using information across genes, and supports nuisance/paired designs [4]. edgeR and limma-voom remain established alternatives; voom explicitly models the mean–variance relationship using precision weights [5,6]. A t-test on log counts is not represented as an equivalent replacement.

**Implementation choice:** optional PyDESeq2 0.5.x integration provides a real count-aware NB/Wald implementation locally in Python [7,8]. It supports binary contrasts, declared covariates, batch, and optional fixed donor effects. Numeric design matrices prevent formula injection. Nonidentifiable or nearly saturated designs fail explicitly. Repeated groups require a declared paired design; donor-constant outcomes cannot be estimated alongside saturated donor fixed effects. Such inputs need a different design or replicate aggregation.

The count integration retains Cook's filtering and reports unshrunk effects and BH-adjusted Wald p-values. Binary contrasts and continuous slopes use backend independent filtering; continuous effects are log2 expression change per original outcome unit. Multiclass uses reference-class contrasts and BH correction across the full gene-by-contrast family, with independent filtering disabled for that joint family. Multiple distribution-mode cell types receive a joint correction across states/contrasts. Missing p-values are not converted into discoveries. Low-count thresholds and actual retained genes are recorded. Exported counts/metadata permit edgeR/DESeq2/limma interoperability; a supplied global DE table cannot enter predictive screening.

There are two intentionally separate DE uses:

1. `candidate_method: deseq2`: every inner-training and outer-training dataset runs its own DE. Its selected genes belong exclusively to that fit.
2. `supporting_de: true`: an explicitly exploratory full-cohort DE analysis occurs after predictive validation. Its evidence is from the same cohort and is **not independent confirmation**.

The faster `association` screen uses a univariate F statistic only as a predictive ranking inside training folds. It emits no DE p-values, does not replace a count model, and does not claim covariate-adjusted inference. The default is outcome-independent training-fold variance ranking.

## Preprocessing and leakage

Varma and Simon demonstrate the optimism introduced by using the same CV results for tuning and performance estimation [9]. Ambroise and McLachlan specifically examine feature-selection bias in high-dimensional expression prediction [10]. These motivate nested validation around the complete feature-selection procedure.

| Operation | Boundary and reason |
|---|---|
| Parse/align IDs; finite/count checks | Input contract; no fitted model |
| Fixed per-cell QC and raw pseudobulk sums | Sample-local, independent of held-out phenotype |
| Gene prevalence and variance filters | Refit on training rows in every inner/outer/subsample fit |
| Library-total log2(1+CPM) | Per-sample arithmetic across the full declared gene universe, before gene selection |
| Frozen median-ratio normalization | Geometric reference learned on training data; reused unchanged on held-out samples |
| Standard scaling and covariate categories | Training-only fitted objects |
| Supervised association/DE screen | Training-only, including during hyperparameter tuning |
| Hyperparameters | Inner CV only |
| Evaluation | Outer held-out biological groups |
| Full-cohort PCA/DE/correlation/path | Descriptive post-evaluation output; no feedback into CV |
| External cohort | Frozen normalization, retained genes, scales, coefficients; no refitting |

Missing raw counts are rejected. Transformed expression/distribution inputs may explicitly request median imputation, fitted inside each training boundary; wholly unavailable training features cannot enter the model. Exact duplicated expression profiles across groups are rejected. Distribution inputs use source-cell fingerprints when available because genuine low-dimensional summaries can tie. Potential target copies trigger provenance review. IDs, grouping, batch, time and outcome columns cannot be ordinary predictors. No audit can identify every concealed target encoding, undocumented preprocessing operation, mislabeled duplicate, or unknown related donor.

## Elastic net and parameter convention

Zou and Hastie motivate the elastic net's combination of sparsity and correlated-feature grouping [11]. SigTrellis uses scikit-learn's tested solvers: SAGA logistic/multinomial elastic net and coordinate-descent least-squares elastic net [12,13]. It does not implement a new optimizer.

Both lambda and `l1_ratio` are tuned. SigTrellis defines lambda on **mean loss**. Regression uses sklearn `alpha=lambda`; classification uses `C=1/(n*lambda)` with weights normalized to sum to n. `l1_ratio` equals glmnet's mixing `alpha`, not sklearn's regression `alpha`. The intercept is unpenalized. Declared covariates are jointly penalized predictors; they are not claimed to provide unpenalized confounder adjustment.

Tuning minimizes proper held-out loss (log loss or squared error), and the optional one-standard-error rule prefers a sparse candidate within a heuristic tolerance of the best loss. Fold overlap means this is a parsimony heuristic, not a formal inferential standard error. Probability calibration, discrimination, class-sensitive metrics and regression errors are reported separately. Classification uses an explicitly declared threshold; it is not optimized on outer outcomes. Calibration plots use held-out probabilities; no test-fitted calibrator is applied [14].

Multinomial models are supported when all training folds retain all classes. Coefficients are reported as logit differences relative to the declared reference class, with contrast-specific stability. Class order and the binary positive class are explicit in every run.

## Stability, correlated predictors, and claims

Meinshausen and Bühlmann's stability selection provides error control under specific assumptions [15]. Retuned elastic-net fits, correlated transcriptomic features, multiple preprocessing choices and grouped subsampling do not automatically satisfy those assumptions. SigTrellis therefore labels frequencies as **empirical/descriptive stability**, not formal PFER control or gene-level FDR.

Every subsample selects complete biological groups and retunes within that subsample. Optional declared normalization perturbations are rotated deterministically; stability is also reported within each normalization. Outputs include nonzero frequency, sign consistency, coefficient distributions, selected-rank variation, outer-fold selection, and nonempty Jaccard summaries. Empty/empty selections do not receive a reassuring stability score.

Marginal expression correlation uses bounded complete-linkage clustering on absolute correlation. It is descriptive, not a learned predictor or pathway inference. A group's frequency of having at least one member selected is reported alongside pairwise exclusive selection and co-selection. This distinguishes stable redundant signal from unstable individual representatives. Phenotype-induced marginal correlation and the candidate cap are stated explicitly.

## Batch, negative controls, and validation targets

Batch effects can invalidate otherwise persuasive genomic signatures [16]. SigTrellis reports batch–outcome association, refuses to identify a perfectly confounded count contrast, and provides purged leave-batch-out validation. It does not globally apply ComBat, integrated embeddings, or outcome-assisted batch correction. Frozen correction is a separate research problem with explicit deployment assumptions [17]. Unknown batches are not made identifiable by a software switch.

Nested label-permutation controls rerun preprocessing, screening, tuning, splitting and evaluation. Group-constant phenotypes are shuffled as group labels; exchangeable paired labels can be permuted within groups. Optional strata must be constant within each group. The finite Monte Carlo p-value uses `(1 + exceedances)/(1 + permutations)` [18]. These are tests of a global predictive statistic under a declared exchangeability assumption, **not** feature-specific tests. With covariates, a global permutation test is not a valid conditional test of genes given the covariates; robustness gates remain exploratory in that case.

Reported predictive improvement is against an intercept/prior baseline or a separately tuned covariate-only baseline. No held-out gain, inadequate permutation resolution, failed controls, small biological cohorts, or unavailable/failed requested batch checks prevent a positive gate. Passing the operational gates does not establish clinical utility. The FDA–NIH BEST resource distinguishes biomarker measurement, validation, and intended use [19]; software selection alone establishes none of them.

Leave-group-out, purged batch holdouts and frozen external evaluation are implemented. Forward temporal splitting uses group start/end times: each training group's end must precede the earliest held-out start minus the declared gap. The same rule applies inside tuning, and temporal CV uses one forward pass. This predicts future samples under the declared ordering; it does not model censoring or survival. Unrestricted temporal permutations are labeled exchangeability-unverified and cannot create positive robustness gates.

## Single-cell information beyond a mean

Pseudobulk with appropriate offsets can match relevant mixed-model properties under specific case-control assumptions [21]; this does not imply that every pseudobulk implementation is equivalent to every cell-level model. The practical reason to retain richer cell information is a different estimand: abundance, detection or distribution shape rather than only mean expression. Milo models neighborhood abundance [22], and scCODA addresses compositional cell populations [23]. Neither treats more recovered cells as more independent donors.

The implemented extension computes fixed state abundance and sample-local means, variances, detection fractions, program quantiles and activation fractions. It is **predictive feature engineering**, not a replacement for these methods' differential-abundance inference. Features retain kinds, units, state and source genes. Programs are predefined equal-weight means of cell-local log1p(CP10K), with no phenotype-selected genes or learned cohort reference. Variability includes technical noise. Abundance is relative recovery. Missing state expression is explicitly unavailable. Cell subsampling perturbs these measurements without multiplying biological sample size. The accompanying scientific fixtures separately plant abundance and variability signals.

We chose this representation because its training boundary is inspectable and frozen deployment is well defined. Hierarchical cell models and learned embeddings remain distinct extensions requiring comparative validation and held-out-cell assignment. Imported annotations may already encode an integrated or phenotype-informed analysis; the package records this upstream limitation rather than claiming to reverse it.

## 0.3.0: distribution shape and within-cell dependence

Memento demonstrates why variance and gene-gene correlation can carry information
beyond expression means, while explicitly modeling single-cell sampling noise
[26]. We adopt the scientific question, not its inference machinery. SigTrellis
computes observed log1p(CP10K) moments/correlations as prediction features; it does
not deconvolve technical variance or produce Memento p-values. Fixed tail fractions
capture rare high-expression populations. Correlations retain pairing within the
same cells, which marginal summaries and count sums can erase. Low-information
correlations remain missing rather than being assigned zero.

Multiple-instance learning treats each sample as a bag of cells. The CELLECTION
preprint explores this phenotype-prediction setting [27]. A neural attention model
would require additional regularization, donor-level nested evaluation, frozen
cell encoders and substantially broader comparative validation. This release uses
a transparent fixed distribution embedding followed by elastic net; it does not
implement that neural method or equate attention with biological causation.

Distribution-preserving sketching is another scale strategy [28]. We use streaming
moments and, in the current public demonstration, all cells in the selected
specimens. We do not claim to implement kernel herding or to preserve the complete
joint distribution. Exact program quantiles retain small per-cell score arrays;
other moments need no retained cell matrix. The preparation writer follows the
AnnData sparse on-disk encoding [29] and is round-trip tested against AnnData.

Independent units remain donors/cultures. Hypothesized gene pairs, states,
programs and thresholds are fixed before validation; data-derived screening,
imputation, scaling and elastic-net tuning stay inside training folds. Cellular
subsampling uses cell-ID-keyed randomness so an unrelated donor or row reordering
cannot alter a cell's inclusion. Exact multiset fingerprints and an explicit full
RNA gene-universe contract make frozen deployment auditable. These engineering
properties prevent specific leakage paths, not hidden upstream confounding.

## Compact-panel validation

Global stability summaries do not inherit the full model's CV performance. With `panel_validation`, each outer training set repeats subsample tuning/selection and applies fixed frequency/sign thresholds plus a maximum panel size. The modal subsample hyperparameters define its refit; there is no second inner-CV search that reuses inner labels after selecting the panel. Only outer predictions estimate this complete policy. The full panel policy is rerun for permutation and batch holdouts. Final panel coefficients and panel-specific outer stability must independently satisfy feature gates. Empty panels become explicit baseline/covariate models. `external --model-kind panel` evaluates the frozen final state without re-selection. Count normalization still requires the original gene universe; a smaller targeted assay needs its own measurement/normalization validation.

## Stabilomics inspection and reuse decision

Reviewed `williamtbarker/stabilomics` at commit `c6676a3dbe09713186faf4980b6f7ca0f1733243`: README, MIT license, `stability.py`, `preprocess.py`, and `model.py` [20].

| Concept/component | Decision |
|---|---|
| Seeded subsampling; selection/sign summaries | Generalized through new group-aware, contrast-aware implementations |
| Jaccard reporting | Retained concept; empty/empty becomes undefined rather than 1.0 |
| Deterministic execution and synthetic validation | Adopted as engineering principles |
| LAD-LASSO linear-program solver | Not reused; different loss, outcome families and scaling requirements |
| Full-table robust scaling | Not reused; violates this project's strict predictive boundary |
| Unpenalized covariate solver | Not ported; no unsupported equivalence claimed |
| Existing report/schema and thresholds | Redesigned for nested evidence and count-specific assumptions |

No Stabilomics source files or data are copied. The repository is MIT-licensed, so appropriately attributed reuse would be permitted; conceptual reuse with independently written code is sufficient here. No proprietary source was accessed.

## Name check

`SigTrellis` combines signature discovery with support for correlated gene groups. On 2026-09-26, exact-name web searches found no established bioinformatics package; GitHub repository search returned zero matches and the PyPI JSON endpoint returned 404. This is a collision check, not trademark clearance or a reserved package name.

## References

1. Squair et al. (2021). Confronting false discoveries in single-cell differential expression. [Nature Communications](https://doi.org/10.1038/s41467-021-25960-2).
2. Hafner et al. (2025). Single-cell differential expression analysis between conditions within nested settings. [Briefings in Bioinformatics](https://doi.org/10.1093/bib/bbaf397).
3. Seurat. [Differential expression vignette](https://satijalab.org/seurat/articles/de_vignette).
4. Love, Huber, Anders (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. [Genome Biology](https://doi.org/10.1186/s13059-014-0550-8). [Current vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).
5. Robinson, McCarthy, Smyth (2010). edgeR. [Bioinformatics](https://doi.org/10.1093/bioinformatics/btp616).
6. Law et al. (2014). voom. [Genome Biology](https://doi.org/10.1186/gb-2014-15-2-r29).
7. Muzellec et al. (2023). PyDESeq2. [Bioinformatics](https://doi.org/10.1093/bioinformatics/btad547).
8. PyDESeq2. [Workflow and design documentation](https://pydeseq2.readthedocs.io/en/stable/auto_examples/plot_minimal_pydeseq2_pipeline.html).
9. Varma, Simon (2006). Bias in error estimation when using cross-validation for model selection. [BMC Bioinformatics](https://doi.org/10.1186/1471-2105-7-91).
10. Ambroise, McLachlan (2002). Selection bias in gene extraction on the basis of microarray gene-expression data. [PNAS](https://doi.org/10.1073/pnas.102102699).
11. Zou, Hastie (2005). Regularization and variable selection via the elastic net. [JRSS B](https://doi.org/10.1111/j.1467-9868.2005.00503.x).
12. scikit-learn. [LogisticRegression](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html).
13. scikit-learn. [ElasticNet](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.ElasticNet.html).
14. scikit-learn. [Probability calibration](https://scikit-learn.org/stable/modules/calibration.html).
15. Meinshausen, Bühlmann (2010). Stability selection. [JRSS B](https://doi.org/10.1111/j.1467-9868.2010.00740.x).
16. Leek et al. (2010). Tackling the widespread and critical impact of batch effects. [Nature Reviews Genetics](https://doi.org/10.1038/nrg2825).
17. Parker, Corrada Bravo, Leek. [Frozen surrogate variable analysis](https://arxiv.org/abs/1301.3947).
18. Phipson, Smyth (2010). Permutation P-values should never be zero. [Statistical Applications in Genetics and Molecular Biology](https://doi.org/10.2202/1544-6115.1585).
19. FDA–NIH. [BEST: Validation](https://www.ncbi.nlm.nih.gov/books/NBK464453/).
20. Barker. [Stabilomics source at inspected commit](https://github.com/williamtbarker/stabilomics/tree/c6676a3dbe09713186faf4980b6f7ca0f1733243).
21. Lee and Han (2024). Pseudobulk with proper offsets has the same statistical properties as generalized linear mixed models in single-cell case-control studies. [Bioinformatics](https://doi.org/10.1093/bioinformatics/btae498).
22. Dann et al. (2022). Differential abundance testing on single-cell data using k-nearest neighbor graphs. [Nature Biotechnology](https://doi.org/10.1038/s41587-021-01033-z).
23. Büttner et al. (2021). scCODA is a Bayesian model for compositional single-cell data analysis. [Nature Communications](https://doi.org/10.1038/s41467-021-27150-6).
24. Schurch et al. (2016). How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? [RNA](https://doi.org/10.1261/rna.053959.115).
25. Perez et al. (2022). Single-cell RNA-seq reveals cell type-specific molecular and genetic associations to lupus. [Science](https://doi.org/10.1126/science.abf1970).
26. Kim et al. (2024). Method of moments framework for differential expression analysis of single-cell RNA sequencing data. [Cell paper](https://doi.org/10.1016/j.cell.2024.09.044), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC11556465/).
27. Hu et al. (2025). Predicting emergent phenotypes from single cell populations using CELLECTION. **bioRxiv preprint**, [doi:10.1101/2025.09.02.673886](https://doi.org/10.1101/2025.09.02.673886), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC12424641/).
28. Distribution-based Sketching of Single-Cell Samples (2022). [arXiv:2207.00584](https://arxiv.org/abs/2207.00584).
29. AnnData. [On-disk format specification](https://anndata.readthedocs.io/en/stable/fileformat-prose.html).
