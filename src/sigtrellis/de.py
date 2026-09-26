"""Optional count-aware DE; exploratory results never flow back into validation."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from sigtrellis.config import Config
from sigtrellis.domain import Dataset, FloatArray


def differential_expression(data: Dataset, y: FloatArray, config: Config) -> pd.DataFrame:
    """Binary NB-GLM Wald contrast, with explicit nuisance design and optional pairs.

    Identifiers are never interpolated into formulas. Rank-deficient or saturated
    designs fail rather than silently dropping confounders or pairing terms.
    """
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError as exc:
        raise ValueError("Install sigtrellis[de] for PyDESeq2") from exc
    if config.outcome_type != "binary" or config.input_scale != "counts":
        raise ValueError("DE requires binary phenotype and integer counts")
    if config.group and data.metadata[config.group].duplicated().any() and not config.de_pair_group:
        raise ValueError(
            "DE with repeated biological groups requires de_pair_group; collapse technical replicates first"
        )
    if config.de_pair_group and config.group:
        replication = pd.crosstab(data.metadata[config.group].astype(str).to_numpy(), y)
        if (replication > 1).any().any():
            raise ValueError(
                "Paired DE supports one observation per biological group and condition; "
                "combine technical replicates or use a model for the actual repeated-measures design"
            )
    nuisance = list(config.covariates)
    if config.batch:
        nuisance.append(config.batch)
    if config.de_pair_group and config.group:
        nuisance.append(config.group)
    columns: list[pd.DataFrame] = [pd.DataFrame({"intercept": 1.0}, index=data.metadata.index)]
    for i, name in enumerate(dict.fromkeys(nuisance)):
        v = data.metadata[name]
        categorical = name in {config.batch, config.group} or not pd.api.types.is_numeric_dtype(v)
        if categorical:
            columns.append(
                pd.get_dummies(v.astype(str), prefix=f"n{i}", drop_first=True, dtype=float)
            )
        elif v.nunique() > 1:
            columns.append(pd.DataFrame({f"n{i}": (v - v.mean()) / v.std()}, index=v.index))
    columns.append(pd.DataFrame({"condition": y}, index=data.metadata.index))
    design = pd.concat(columns, axis=1).astype(float)
    if np.linalg.matrix_rank(design.to_numpy()) < design.shape[1]:
        raise ValueError(
            "DE design is rank deficient: phenotype, batch or pairing are not identifiable"
        )
    if len(design) - design.shape[1] < 2:
        raise ValueError("DE design has fewer than two residual degrees of freedom")
    counts = data.expression
    eligible = (counts >= config.min_count).mean(axis=0) >= config.min_prevalence
    counts = counts.loc[:, eligible].astype(np.int64)
    if counts.shape[1] < 2:
        raise ValueError("Too few expressed genes for count-aware DE")
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        dds = DeseqDataSet(
            counts=counts,
            metadata=data.metadata.copy(),
            design=design,
            size_factors_fit_type="poscounts",
            n_cpus=1,
            quiet=True,
            refit_cooks=True,
            low_memory=True,
        )
        dds.deseq2()
        contrast = np.zeros(design.shape[1], dtype=float)
        contrast[-1] = 1.0
        stats = DeseqStats(
            dds, contrast=contrast, n_cpus=1, quiet=True, independent_filter=True, cooks_filter=True
        )
        stats.summary()
    result: pd.DataFrame = stats.results_df.copy()
    result.index.name = "gene_id"
    result = result.rename(
        columns={
            "log2FoldChange": "de_log2_fold_change",
            "padj": "de_adjusted_pvalue",
            "pvalue": "de_pvalue",
        }
    )
    result.attrs["warnings"] = sorted({str(w.message) for w in captured})
    result.attrs["design_columns"] = list(design.columns)
    result.attrs["interpretation"] = "Unshrunk count-model LFC; BH-adjusted Wald p-values"
    return result
