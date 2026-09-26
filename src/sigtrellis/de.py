"""Optional count-aware DE; exploratory results never flow back into validation."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

from sigtrellis.config import Config
from sigtrellis.domain import Dataset, FloatArray


def differential_expression(data: Dataset, y: FloatArray, config: Config) -> pd.DataFrame:
    """NB-GLM Wald contrasts with nuisance terms and optional fixed pairing.

    Identifiers are never interpolated into formulas. Rank-deficient or saturated
    designs fail rather than silently dropping confounders or pairing terms.
    """
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError as exc:
        raise ValueError("Install sigtrellis[de] for PyDESeq2") from exc
    if config.input_scale != "counts":
        raise ValueError("DE requires integer counts")
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
    if config.outcome_type == "multiclass":
        condition_columns = [f"condition_{i}" for i in sorted(np.unique(y).astype(int))[1:]]
        columns.append(
            pd.DataFrame(
                {
                    f"condition_{i}": (y == i).astype(float)
                    for i in sorted(np.unique(y).astype(int))[1:]
                },
                index=data.metadata.index,
            )
        )
    else:
        condition_columns = ["condition"]
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
        results: list[pd.DataFrame] = []
        from sigtrellis.qc import encode_outcome

        _, labels = encode_outcome(data, config)
        contrasts = (
            ["response_per_training_sd"]
            if not labels
            else [f"{v}_vs_{labels[0]}" for v in labels[1:]]
        )
        for name, label in zip(condition_columns, contrasts, strict=True):
            contrast = np.zeros(design.shape[1], dtype=float)
            contrast[list(design.columns).index(name)] = 1.0
            stats = DeseqStats(
                dds,
                contrast=contrast,
                n_cpus=1,
                quiet=True,
                independent_filter=config.outcome_type != "multiclass",
                cooks_filter=True,
            )
            stats.summary()
            frame = stats.results_df.copy()
            frame["contrast"] = label
            frame["de_effect_unit"] = (
                "log2_expression_change_per_outcome_unit" if not labels else "log2_fold_change"
            )
            results.append(frame)
    result = pd.concat(results)
    result.index.name = "gene_id"
    result = result.rename(
        columns={
            "log2FoldChange": "de_log2_fold_change",
            "padj": "de_adjusted_pvalue",
            "pvalue": "de_pvalue",
        }
    )
    if config.outcome_type == "multiclass":
        # Correct the full gene-by-reference-contrast family, rather than
        # choosing a favorable uncorrected contrast for each gene.
        valid = result.de_pvalue.notna().to_numpy()
        result["de_adjusted_pvalue"] = np.nan
        result.loc[valid, "de_adjusted_pvalue"] = false_discovery_control(
            result.loc[valid, "de_pvalue"].to_numpy()
        )
    result.attrs["warnings"] = sorted({str(w.message) for w in captured})
    result.attrs["design_columns"] = list(design.columns)
    result.attrs["interpretation"] = (
        "Unshrunk count-model effects; Wald tests; multiclass BH across genes and reference contrasts"
    )
    result.attrs["outcome_unit"] = (
        "original numeric outcome unit"
        if config.outcome_type == "continuous"
        else "reference-class contrast"
    )
    return result
