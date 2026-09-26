# Public demonstration data and reuse terms

The MIT license applies to SigTrellis software. It does **not** relicense these independently published datasets. No raw dataset is included in the source repository or source distribution. `examples/prepare_public.py` downloads public files, verifies SHA-256, validates their schemas, and records every preparation step. A changed source fails verification rather than silently changing the demonstration.

## Replicated bulk: independent yeast cultures

- Schurch et al. (2016), [RNA, doi:10.1261/rna.053959.115](https://doi.org/10.1261/rna.053959.115), accession ERP004763.
- Source: [bartongroup/profDGE48](https://github.com/bartongroup/profDGE48/tree/375dc0d57d9d1fa96a4245a6530e0fda34305891), pinned commit `375dc0d57d9d1fa96a4245a6530e0fda34305891`. The distribution has an MIT license; its original copyright/license notice is retained in `SOURCE_NOTICES/yeast_MIT.txt`.
- `Preprocessed_data/Snf2_countdata.tar.gz`: SHA-256 `8cafcaa7052a760cbec9196605f1e588c58c76f772db064ace599c0bf07d727a`.
- `Preprocessed_data/WT_countdata.tar.gz`: SHA-256 `11b33cf3fb8bd9b7c229c81ec00f3ab502f61a1715d0f59821b9bc7a7eb625ac`.
- `Bad_replicate_identification/exclude.lst`: SHA-256 `57d5e5a7512388b7df2255095b0cb0938438658256262e29b60cd12056214a5e`.
- Ninety-six independently grown cultures, 86 after published QC exclusions, with technical lanes already pooled. Preparation removes HTSeq summary rows, applies the source exclusions, aligns genes, and reserves 20 cultures using seed 2026 before fitting. The inherited published QC cannot be audited as though this tool performed it prospectively.
- Source and processing provenance are written by `examples/prepare_replicated.py`; full commands are in `PUBLIC_DEMOS.md`.

## Replicated single-cell: lupus case-control study

- Perez et al. (2022), [Science, doi:10.1126/science.abf1970](https://doi.org/10.1126/science.abf1970).
- [CELLxGENE collection](https://cellxgene.cziscience.com/collections/436154da-bcf1-4130-9c8b-120ff9a888f2), dataset `218acb0f-9f2f-4f76-b90b-15a4b7c7f629`.
- [Pinned H5AD asset](https://datasets.cellxgene.cziscience.com/c55dc602-d168-4d15-acc1-5de4f2f5d551.h5ad): 12,218,105,530 bytes; SHA-256 `3c0b74d54c03838a49817edce95314e6a4fa048ef35b74d1e697b8b2fd07cc03`.
- **CC BY 4.0**, attributed to Perez et al. and the data contributors, curated/distributed by CZ CELLxGENE Discover. [Public submission terms](https://cellxgene.cziscience.com/docs/032__Contribute%20and%20Publish%20Data) explicitly describe this license. Included derived reports retain this attribution; changes are recorded in preparation provenance.
- Source: 261 donors, 1,263,676 cells, 30,172 genes. `X` is normalized, so preparation explicitly uses **raw/X** integer counts. Source cell labels and filtering are inherited.
- Some specimen IDs occur in multiple processing cohorts. Preparation defines specimen/processing aliquots before assigning cohort 4 to holdout, then purges every holdout donor from training cohorts 2/3. One aliquot per donor gives 137 training and 96 held-out donors. The other control-only processing cohort is excluded by design, not based on model performance.
- Within each chosen aliquot, seed 2026 selects at most 400 cells without expression-based filtering; all source genes are retained. The prepared files contain 54,800 training and 38,400 held-out cells. The monocyte-specific minimum-cell rule may exclude additional specimens and records these exclusions.
- Frozen holdout evaluation concerns another processing cohort of the same observational study. It is not an independent clinical validation, and it cannot separate every effect of disease, treatment, ancestry and other nuisance factors. No source phenotype or donor label is an ordinary gene predictor.

## Bulk: paired airway smooth-muscle RNA-seq

- Study: Himes et al. (2014), *RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells*. [DOI](https://doi.org/10.1371/journal.pone.0099625), [GEO GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778).
- Public count/metadata distribution: [bioconnector/workshops](https://github.com/bioconnector/workshops), commit `68ff3e4868b0e5673a1c4399828a3e0adc72f8e4`. Count source is `data/airway_rawcounts.csv`, not the length-scaled alternative. The [Bioconductor airway documentation](https://bioconductor.org/packages/release/data/experiment/html/airway.html) describes the four cell lines and paired treatment design.
- This distribution's [license](https://github.com/bioconnector/workshops/blob/68ff3e4868b0e5673a1c4399828a3e0adc72f8e4/LICENSE) is **CC BY-NC-SA 4.0**. Use the example as a noncommercial reproducibility demonstration and retain attribution. Included derived bulk demonstration reports carry those source reuse terms separately from MIT code. The underlying study article is open access; do not assume that overrides a downstream distributor's terms for its files.
- 64,102 gene rows, eight samples, four paired cell-line/donor identifiers. `celltype` in the distribution identifies the biological cell line and is renamed to `donor`; it is not an scRNA-seq cell-type annotation.
- Preparation only converts CSV to TSV and renames metadata columns. No phenotype-informed gene selection, normalization, or sample exclusion occurs in preparation.
- Raw counts SHA-256: `504f3148e23e84061f3f28d109a37732ccc9f1f8e4833b00b105307f2f8eb326`.
- Metadata SHA-256: `05bd7e78a0ca5b2a2f60ec715296381ca19f31a8c5a51e70c9f60583fa9fcf97`.

The demonstration holds the two conditions from a cell line together. Four independent biological groups are insufficient to establish a transferable biomarker. Nineteen Monte Carlo permutations cannot create more independent paired labelings than this experiment contains. Count DE uses an identifiable paired design, but its estimates remain fragile in such a small study.

## Single-cell: paired Kang IFN-beta experiment

- Study: Kang et al. (2018), *Multiplexed droplet single-cell RNA-sequencing using natural genetic variation*. [DOI](https://doi.org/10.1038/nbt.4042), [GEO GSE96583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583).
- Open deposit: Daniel Dimitrov (2023), *SingleCellExperiment Object for Kang et al., 2018 PBMCs data*. [Zenodo DOI 10.5281/zenodo.10069528](https://doi.org/10.5281/zenodo.10069528). The deposit declares **CC BY 4.0**, and explicitly identifies its relationship to the original Pertpy AnnData deposit. [Pertpy documentation](https://pertpy.readthedocs.io/en/latest/api/data/pertpy.data.kang_2018.html) provides experimental context.
- Download: `kang_counts_25k.RDS`, 37,150,344 bytes; SHA-256 `eb8a88259d84f757fcdb9a011085f62c73d36a59487249066c550001b29ea6a6`; deposit MD5 `ef6531d1cc2825caa9133c3ef416f76c`.
- 24,673 cells × 15,706 genes; eight biological donors, 16 donor-condition samples. The supplied raw-count dgCMatrix is transposed to sparse cells × genes H5AD. Sample, donor, condition and supplied cell-type annotations are preserved. Selected deterministic chunks are compared after serialization.
- R is not required. The optional `rdata` dependency parses this exact, checksum-verified S4 schema. This example adapter is not a general arbitrary-Seurat/S4 importer.
- The demonstrated cell type is the supplied `CD14+ Monocytes` annotation. No type is selected based on predictive performance.

The two conditions were processed as condition-specific pooled libraries. The preparation script explicitly records that known design as `condition_library`; it is **derived from the published condition/library design**, not an invented independently measured technical batch. Condition and this library factor are perfectly confounded. Donor holdout tests prediction across donors under that same protocol, not separation of treatment biology from pooled-library effects. The package therefore blocks positive robustness claims and reports library-adjusted DE as nonidentifiable.

The deposit already contains cell/gene filtering and cell-type labels. The demonstration cannot retrospectively prove that those upstream choices were independent of phenotype. It validates the executable downstream workflow, not the whole original study or the biological truth of its selected genes.

## Included example evidence

Release bundles place executed reports and audits under `release_evidence/`, outside the source-distribution payload. The raw data stay external. Bulk-derived evidence is attributed to Himes et al. and the bioconnector workshop distribution under CC BY-NC-SA 4.0; single-cell-derived evidence is attributed to Kang et al. and Daniel Dimitrov's deposit under CC BY 4.0. Synthetic fixtures, software, and original prose documentation are MIT-licensed.
