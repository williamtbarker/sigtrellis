"""All cohort-learned transformations live in this training-only object."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd
from sklearn.feature_selection import f_classif, f_regression
from sklearn.preprocessing import StandardScaler

from sigtrellis.config import Config
from sigtrellis.de import differential_expression
from sigtrellis.domain import Audit, Dataset, FloatArray, array_hash
from sigtrellis.qc import group_values


@dataclass
class Normalizer:
    method: str
    reference: FloatArray | None = None
    columns: list[str] = field(default_factory=list)

    def fit(self, x: pd.DataFrame) -> Normalizer:
        self.columns = list(x.columns)
        if self.method == "median_ratio":
            values = x.to_numpy(dtype=float)
            usable = (values > 0).all(axis=0)
            if usable.sum() < 2:
                raise ValueError(
                    "Too few positive reference genes for frozen median-ratio normalization"
                )
            self.reference = np.full(values.shape[1], np.nan)
            self.reference[usable] = np.exp(np.log(values[:, usable]).mean(axis=0))
        return self

    def transform(self, x: pd.DataFrame) -> FloatArray:
        if list(x.columns) != self.columns:
            raise ValueError("Gene universe/order differs from the training normalization contract")
        values = x.to_numpy(dtype=np.float64)
        if self.method == "none":
            return values.copy()
        if self.method == "logcpm":
            totals = values.sum(axis=1)
            if (totals <= 0).any():
                raise ValueError("Zero-library sample")
            return np.asarray(np.log2(1 + values / totals[:, None] * 1e6), dtype=np.float64)
        if self.reference is None:
            raise ValueError("Normalizer has not been fitted")
        usable = np.isfinite(self.reference)
        ratios = values[:, usable] / self.reference[usable]
        ratios[ratios <= 0] = np.nan
        factors = np.nanmedian(ratios, axis=1)
        if not np.isfinite(factors).all() or (factors <= 0).any():
            raise ValueError("Sample has no valid positive training-reference ratios")
        return np.asarray(np.log2(1 + values / factors[:, None]), dtype=np.float64)


@dataclass
class CovariateEncoder:
    names: tuple[str, ...]
    categories: dict[str, list[str]] = field(default_factory=dict)
    columns: list[str] = field(default_factory=list)
    scaler: Any = None

    def _encode(self, metadata: pd.DataFrame) -> pd.DataFrame:
        result = pd.DataFrame(index=metadata.index)
        for name in self.names:
            if name in self.categories:
                categories = self.categories[name]
                if not set(metadata[name].astype(str)) <= set(categories):
                    raise ValueError(f"Covariate {name} contains categories absent from training")
                for category in categories[1:]:
                    result[f"covariate:{name}={category}"] = (
                        metadata[name].astype(str) == category
                    ).astype(float)
            else:
                result[f"covariate:{name}"] = metadata[name].astype(float)
        return result

    def fit(self, metadata: pd.DataFrame) -> CovariateEncoder:
        self.categories = {
            name: sorted(metadata[name].astype(str).unique())
            for name in self.names
            if not pd.api.types.is_numeric_dtype(metadata[name])
        }
        raw = self._encode(metadata)
        self.columns = list(raw.columns)
        if self.columns:
            self.scaler = StandardScaler().fit(raw)
        return self

    def transform(self, metadata: pd.DataFrame) -> FloatArray:
        raw = self._encode(metadata)
        if not self.columns:
            return np.empty((len(metadata), 0), dtype=float)
        return np.asarray(self.scaler.transform(raw), dtype=float)


@dataclass
class Prepared:
    config: Config
    normalizer: Normalizer = field(init=False)
    covariates: CovariateEncoder = field(init=False)
    scaler: Any = None
    selected: np.ndarray[Any, Any] = field(default_factory=lambda: np.array([], dtype=int))
    genes: list[str] = field(default_factory=list)
    all_genes: list[str] = field(default_factory=list)

    def fit(self, data: Dataset, y: FloatArray, audit: Audit, context: str) -> Prepared:
        x = data.expression
        self.all_genes = list(x.columns)
        self.normalizer = Normalizer(self.config.normalization).fit(x)
        normalized = self.normalizer.transform(x)
        values = x.to_numpy(dtype=float)
        if self.config.input_scale == "counts":
            prevalence = (values >= self.config.min_count).mean(axis=0)
            eligible = prevalence >= self.config.min_prevalence
        else:
            eligible = np.ones(values.shape[1], dtype=bool)
        variance = normalized.var(axis=0)
        eligible &= variance > 1e-12
        indices = np.flatnonzero(eligible)
        score = variance.copy()
        method = self.config.candidate_method
        de_warnings: list[str] = []
        if method == "association" and len(indices):
            # Ranking statistic only: no DE p-values or replication claims are emitted.
            if self.config.outcome_type == "continuous":
                local_score = f_regression(normalized[:, indices], y)[0]
            else:
                local_score = f_classif(normalized[:, indices], y)[0]
            score[indices] = np.nan_to_num(local_score, nan=-np.inf, posinf=np.finfo(float).max)
        elif method == "deseq2":
            de = differential_expression(data, y, self.config)
            de_warnings = de.attrs["warnings"]
            padj = de["de_adjusted_pvalue"].reindex(x.columns).fillna(1).to_numpy()
            eligible &= padj <= self.config.de_fdr
            indices = np.flatnonzero(eligible)
            score = -padj
        # Stable ties follow the declared gene order, never held-out values.
        if method != "none":
            indices = indices[
                np.argsort(-score[indices], kind="stable")[: self.config.max_features]
            ]
        self.selected = np.sort(indices)
        self.genes = [self.all_genes[i] for i in self.selected]
        if len(indices):
            self.scaler = StandardScaler().fit(normalized[:, self.selected])
        self.covariates = CovariateEncoder(self.config.covariates).fit(data.metadata)
        reference = self.normalizer.reference
        retained_set = set(self.genes)
        audit.fits.append(
            {
                "context": context,
                "sample_ids": x.index.tolist(),
                "group_ids": sorted(set(group_values(data, self.config))),
                "input_expression_hash": array_hash(values),
                "outcome_hash": array_hash(y),
                "normalization": self.config.normalization,
                "reference_hash": array_hash(reference) if reference is not None else None,
                "candidate_method": method,
                "n_input_features": x.shape[1],
                "n_prevalence_variance_eligible": int(eligible.sum()),
                "retained_genes": self.genes,
                "removed_genes": [g for g in self.all_genes if g not in retained_set],
                "gene_scale_hash": array_hash(np.asarray(self.scaler.scale_, dtype=float))
                if self.scaler is not None
                else None,
            }
        )
        audit.warnings.extend(f"DE_WARNING [{context}]: {message}" for message in de_warnings)
        return self

    def transform(self, data: Dataset) -> FloatArray:
        normalized = self.normalizer.transform(data.expression)
        if len(self.selected):
            genes = np.asarray(self.scaler.transform(normalized[:, self.selected]), dtype=float)
        else:
            genes = np.empty((len(data.expression), 0))
        return np.column_stack([genes, self.covariates.transform(data.metadata)])
