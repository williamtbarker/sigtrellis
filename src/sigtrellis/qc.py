"""Fail-closed data checks plus explicitly descriptive confounding diagnostics."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

from sigtrellis.config import Config
from sigtrellis.domain import Dataset, FloatArray, array_hash


def group_values(data: Dataset, config: Config) -> np.ndarray[Any, Any]:
    return data.metadata[config.group or config.sample_id].astype(str).to_numpy()


def encode_outcome(data: Dataset, config: Config) -> tuple[FloatArray, list[str]]:
    series = data.metadata[config.outcome]
    if series.isna().any():
        raise ValueError("Missing outcome values")
    if config.outcome_type == "continuous":
        y = series.to_numpy(dtype=np.float64)
        if not np.isfinite(y).all() or np.std(y) < 1e-12:
            raise ValueError("Continuous outcome must be finite and variable")
        return y, []
    labels = sorted(series.astype(str).unique())
    if config.outcome_type == "binary":
        if len(labels) != 2:
            raise ValueError("Binary outcome requires exactly two classes")
        if config.positive_class:
            if config.positive_class not in labels:
                raise ValueError("positive_class does not occur in outcome")
            labels.remove(config.positive_class)
            labels.append(config.positive_class)
    elif len(labels) < 3:
        raise ValueError("Multiclass outcome requires >= 3 classes")
    mapping = {label: i for i, label in enumerate(labels)}
    return series.astype(str).map(mapping).to_numpy(dtype=np.float64), labels


def validate_dataset(data: Dataset, config: Config) -> dict[str, Any]:
    config.validate()
    x, m = data.expression, data.metadata
    if not x.index.equals(m.index):
        raise ValueError("Expression/metadata row order mismatch")
    if x.index.has_duplicates or x.columns.has_duplicates:
        raise ValueError("Duplicate sample or gene identifiers")
    if x.shape[0] < 4 or x.shape[1] < 2:
        raise ValueError("Need at least four samples and two candidate genes")
    if any(not str(v).strip() or str(v).lower() == "nan" for v in [*x.index, *x.columns]):
        raise ValueError("Empty/invalid sample or gene identifier")
    required = {config.outcome, config.sample_id, *config.covariates}
    required.update(v for v in (config.group, config.batch, config.permutation_strata) if v)
    missing = required - set(m.columns)
    if missing:
        raise ValueError(f"Missing metadata columns: {sorted(missing)}")
    if m[list(required)].isna().any().any():
        raise ValueError("Missing required metadata; no automatic imputation")
    if m[config.sample_id].astype(str).tolist() != x.index.astype(str).tolist():
        raise ValueError("Sample column must equal expression row identifiers")
    reserved = required | {"sample_id", "donor", "batch", "phenotype", "outcome"}
    if set(x.columns) & reserved:
        raise ValueError("Metadata/identifier/outcome column appears in expression feature matrix")
    values = x.to_numpy(dtype=np.float64)
    if not np.isfinite(values).all():
        raise ValueError("Expression must be finite; missing values are rejected")
    if values.nbytes > config.max_dense_mb * 1024**2:
        raise ValueError("Sample-level matrix exceeds max_dense_mb")
    if config.input_scale == "counts":
        if (values < 0).any() or not np.allclose(values, np.rint(values), atol=1e-6, rtol=0):
            raise ValueError(
                "Raw counts must be nonnegative integers; normalized data are not counts"
            )
        if (values.sum(axis=1) <= 0).any():
            raise ValueError("Zero-library samples must be reviewed before modeling")
    y, classes = encode_outcome(data, config)
    groups = group_values(data, config)
    n_groups = len(np.unique(groups))
    if n_groups < 4:
        raise ValueError("Fewer than four biological groups cannot support nested validation")
    hashes: dict[str, str] = {}
    for row, group in zip(values, groups, strict=True):
        key = array_hash(row)
        if key in hashes and hashes[key] != group:
            raise ValueError("Identical expression samples occur in distinct biological groups")
        hashes[key] = group
    messages: list[str] = []
    blockers: list[str] = []
    if n_groups < config.min_groups_gate:
        messages.append(f"SMALL_COHORT: {n_groups} biological groups; uncertainty is substantial")
        blockers.append("too_few_biological_groups")
    if config.input_scale == "log_expression":
        messages.append(
            "UPSTREAM_PREPROCESSING: provenance cannot establish external leakage safety"
        )
    if not config.group:
        messages.append(
            "GROUP_ASSUMPTION: each declared sample is assumed biologically independent"
        )
    if config.outcome_type != "continuous":
        for cls in range(len(classes)):
            if len(np.unique(groups[y == cls])) < 4:
                messages.append(f"SMALL_CLASS: {classes[cls]} has <4 biological groups")
        counts = np.bincount(y.astype(int))
        if counts.min() / counts.sum() < 0.2:
            messages.append(
                "CLASS_IMBALANCE: inspect PR-AUC, sensitivity, specificity and calibration"
            )
    # Perfect association cannot itself establish that a gene is a leaked label.
    # Low-count genes that cannot enter ANY fold need no encoding alarm. Other
    # exact low-cardinality/affine copies are review blockers, never removed using y.
    suspicious: list[str] = []
    for col in range(values.shape[1]):
        feature_values = values[:, col]
        if config.input_scale == "counts" and feature_values.max() < config.min_count:
            continue
        if config.outcome_type != "continuous" and np.unique(feature_values).size > len(classes):
            continue
        if np.std(feature_values) > 0 and abs(np.corrcoef(feature_values, y)[0, 1]) > 1 - 1e-12:
            suspicious.append(str(x.columns[col]))
    if suspicious:
        messages.append(
            "POSSIBLE_OUTCOME_ENCODING: exact low-cardinality/affine associations require provenance review; genuine biology cannot be excluded"
        )
        blockers.append("possible_outcome_encoding_requires_review")
    for name in config.covariates:
        v = m[name]
        if not pd.api.types.is_numeric_dtype(v):
            if v.nunique() > max(10, len(v) // 2):
                raise ValueError(f"Identifier-like high-cardinality covariate: {name}")
            if config.outcome_type != "continuous":
                table = pd.crosstab(v, y)
                if (table.gt(0).sum(axis=1) == 1).all():
                    raise ValueError(f"Covariate deterministically encodes outcome: {name}")
        elif not np.isfinite(v.to_numpy(dtype=float)).all():
            raise ValueError(f"Nonfinite covariate: {name}")
        elif np.std(v) > 0 and abs(np.corrcoef(v, y)[0, 1]) > 1 - 1e-12:
            raise ValueError(f"Covariate is an affine copy of outcome: {name}")
    for sample, label in zip(x.index.astype(str), m[config.outcome].astype(str), strict=True):
        if len(label) >= 3 and label.casefold() in sample.casefold():
            messages.append("ID_CONTAINS_OUTCOME: IDs are excluded from all predictor matrices")
            break
    batch_info: dict[str, Any] = {"status": "not_assessed"}
    if config.batch:
        b = m[config.batch].astype(str)
        if config.outcome_type != "continuous":
            table = pd.crosstab(b, y)
            purity = float(table.max(axis=1).sum() / len(y))
            cramers_v = 0.0
            if min(table.shape) > 1:
                chi = chi2_contingency(table, correction=False)[0]
                cramers_v = float(np.sqrt(chi / (len(y) * (min(table.shape) - 1))))
            batch_info = {
                "status": "observed",
                "cramers_v": cramers_v,
                "purity": purity,
                "table": table.to_dict(),
            }
            if purity == 1:
                messages.append(
                    "PERFECT_BATCH_CONFOUNDING: phenotype cannot be separated from batch"
                )
                blockers.append("perfect_batch_confounding")
            elif cramers_v > 0.5:
                messages.append("BATCH_ASSOCIATION: phenotype and batch strongly associated")
        else:
            means = pd.Series(y, index=m.index).groupby(b).transform("mean").to_numpy()
            eta2 = float(np.var(means) / np.var(y))
            batch_info = {"status": "observed", "eta_squared": eta2}
            if eta2 > 0.5:
                messages.append("BATCH_ASSOCIATION: batch explains >50% of outcome variance")
            if eta2 > 1 - 1e-10:
                blockers.append("perfect_batch_confounding")
    library = values.sum(axis=1) if config.input_scale == "counts" else None
    if library is not None:
        log_library = np.log1p(library)
        mad = float(np.median(abs(log_library - np.median(log_library))))
        outlier = abs(log_library - np.median(log_library)) > max(5 * 1.4826 * mad, 1.5)
        if outlier.any():
            messages.append("LIBRARY_OUTLIERS: flagged for review; no data-driven sample removal")
    else:
        outlier = np.zeros(len(x), dtype=bool)
    return {
        "n_samples": len(x),
        "n_genes": x.shape[1],
        "n_biological_groups": n_groups,
        "classes": classes,
        "warnings": messages,
        "gate_blockers": blockers,
        "batch": batch_info,
        "library_sizes": library,
        "outlier_samples": x.index[outlier].tolist(),
        "possible_outcome_encoded_features": suspicious,
        "upstream": data.upstream_qc,
    }
