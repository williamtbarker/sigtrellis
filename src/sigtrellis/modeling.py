"""Trusted scikit-learn solvers with an explicit mean-loss penalty convention."""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from sklearn.dummy import DummyClassifier, DummyRegressor
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import ElasticNet, LogisticRegression

from sigtrellis.config import Config
from sigtrellis.domain import Audit, Dataset, FloatArray
from sigtrellis.metrics import group_weights, loss
from sigtrellis.preprocessing import Prepared
from sigtrellis.qc import group_values
from sigtrellis.splits import make_splits, record_split


@dataclass(frozen=True)
class Hyperparameters:
    strength: float
    l1_ratio: float


@dataclass
class FittedModel:
    prepared: Prepared
    estimator: Any
    parameters: Hyperparameters
    n_classes: int
    covariates_only: bool = False
    dummy: bool = False

    def predict(self, data: Dataset) -> FloatArray:
        x = self.prepared.transform(data)
        if self.covariates_only:
            x = x[:, len(self.prepared.genes) :]
        if self.dummy:
            x = np.zeros((len(x), 1))
        if self.prepared.config.outcome_type == "continuous":
            return np.asarray(self.estimator.predict(x), dtype=float)
        return np.asarray(self.estimator.predict_proba(x), dtype=float)

    def coefficients(self) -> FloatArray:
        """Per-training-SD gene coefficients; multiclass logits relative to class 0."""
        p = len(self.prepared.all_genes)
        n_contrasts = max(1, self.n_classes - 1)
        result = np.zeros((n_contrasts, p), dtype=float)
        if self.dummy or self.covariates_only:
            return result
        raw = np.atleast_2d(np.asarray(self.estimator.coef_, dtype=float))
        if self.n_classes > 2:
            raw = raw[1:] - raw[0]
        result[:, self.prepared.selected] = raw[:, : len(self.prepared.genes)]
        return result

    def to_state(self) -> dict[str, Any]:
        """Portable numeric state; no unsafe pickle is needed to preserve evidence."""
        prep = self.prepared
        return {
            "genes": prep.all_genes,
            "selected_indices": prep.selected,
            "normalization": prep.normalizer.method,
            "reference": prep.normalizer.reference,
            "imputation": prep.normalizer.imputation,
            "fill_values": prep.normalizer.fill_values,
            "gene_center": prep.scaler.mean_ if prep.scaler is not None else [],
            "gene_scale": prep.scaler.scale_ if prep.scaler is not None else [],
            "covariates": prep.covariates.names,
            "covariate_categories": prep.covariates.categories,
            "covariate_columns": prep.covariates.columns,
            "covariate_center": prep.covariates.scaler.mean_
            if prep.covariates.scaler is not None
            else [],
            "covariate_scale": prep.covariates.scaler.scale_
            if prep.covariates.scaler is not None
            else [],
            "raw_coefficients": getattr(self.estimator, "coef_", None),
            "intercept": getattr(self.estimator, "intercept_", None),
            "class_prior": getattr(self.estimator, "class_prior_", None),
            "constant": getattr(self.estimator, "constant_", None),
            "outcome_type": prep.config.outcome_type,
            "n_classes": self.n_classes,
            "strength": self.parameters.strength,
            "l1_ratio": self.parameters.l1_ratio,
            "dummy": self.dummy,
            "covariates_only": self.covariates_only,
        }


def fit_prepared(
    prepared: Prepared,
    data: Dataset,
    y: FloatArray,
    parameters: Hyperparameters,
    *,
    covariates_only: bool = False,
) -> FittedModel:
    config = prepared.config
    x = prepared.transform(data)
    if covariates_only:
        x = x[:, len(prepared.genes) :]
    weights = group_weights(group_values(data, config))
    n_classes = 0 if config.outcome_type == "continuous" else len(np.unique(y))
    dummy = x.shape[1] == 0
    if dummy:
        x = np.zeros((len(y), 1))
        estimator: Any = (
            DummyRegressor(strategy="mean") if not n_classes else DummyClassifier(strategy="prior")
        )
    elif not n_classes:
        estimator = ElasticNet(
            alpha=parameters.strength,
            l1_ratio=parameters.l1_ratio,
            max_iter=config.max_iter,
            tol=config.tolerance,
            selection="cyclic",
        )
    else:
        # sklearn penalizes summed weighted log loss. sum(weights)=n, so this
        # converts lambda on MEAN log loss to sklearn C exactly.
        estimator = LogisticRegression(
            solver="saga",
            l1_ratio=parameters.l1_ratio,
            C=1 / (len(y) * parameters.strength),
            max_iter=config.max_iter,
            tol=config.tolerance,
            random_state=config.seed,
        )
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        estimator.fit(x, y, sample_weight=weights)
    return FittedModel(prepared, estimator, parameters, n_classes, covariates_only, dummy)


def fit_model(
    data: Dataset,
    y: FloatArray,
    config: Config,
    parameters: Hyperparameters,
    audit: Audit,
    context: str,
    *,
    allowed_features: set[str] | None = None,
) -> FittedModel:
    prepared = Prepared(config)
    if allowed_features is None:
        prepared.fit(data, y, audit, context)
    else:
        prepared.fit(data, y, audit, context, allowed_features)
    return fit_prepared(prepared, data, y, parameters)


def tune(
    data: Dataset,
    y: FloatArray,
    config: Config,
    audit: Audit,
    context: str,
    *,
    covariates_only: bool = False,
) -> tuple[Hyperparameters, pd.DataFrame]:
    parameters = [Hyperparameters(s, r) for s in config.strengths for r in config.l1_ratios]
    records: list[dict[str, Any]] = []
    for split in make_splits(data, y, config, inner=True):
        label = f"{context}/inner:{split.fold}"
        record_split(audit, data, split, label, config)
        train, test = data.subset(split.train), data.subset(split.test)
        prepared = Prepared(config).fit(train, y[split.train], audit, label)
        for param in parameters:
            try:
                model = fit_prepared(
                    prepared, train, y[split.train], param, covariates_only=covariates_only
                )
                prediction = model.predict(test)
                error = loss(
                    y[split.test], prediction, config, group_weights(group_values(test, config))
                )
                n_selected = int(
                    (np.abs(model.coefficients()) > config.coefficient_tolerance).any(axis=0).sum()
                )
                failure = ""
            except ConvergenceWarning as exc:
                error, n_selected, failure = float("inf"), 0, str(exc)
                audit.warnings.append(f"SOLVER_NONCONVERGENCE [{label}]: {param}")
            records.append(
                {
                    "context": context,
                    "inner_fold": split.fold,
                    "strength": param.strength,
                    "l1_ratio": param.l1_ratio,
                    "loss": error,
                    "n_selected": n_selected,
                    "failure": failure,
                }
            )
    table = pd.DataFrame(records)
    aggregate = table.groupby(["strength", "l1_ratio"], sort=True).agg(
        mean_loss=("loss", "mean"),
        sd_loss=("loss", "std"),
        n=("loss", "size"),
        mean_features=("n_selected", "mean"),
    )
    aggregate = aggregate[np.isfinite(aggregate["mean_loss"])].copy()
    if aggregate.empty:
        raise ValueError("Every hyperparameter setting failed to converge")
    best = aggregate.loc[aggregate["mean_loss"].idxmin()]
    cutoff = float(best["mean_loss"])
    if config.tuning_rule == "one_se":
        # A parsimony heuristic; overlapping CV fits do not yield a formal SE.
        cutoff += float(best["sd_loss"] / np.sqrt(best["n"]))
    allowed = aggregate[aggregate["mean_loss"] <= cutoff + 1e-12].reset_index()
    allowed = allowed.sort_values(
        ["mean_features", "strength", "l1_ratio"], ascending=[True, False, False], kind="stable"
    )
    row = allowed.iloc[0]
    chosen = Hyperparameters(float(row["strength"]), float(row["l1_ratio"]))
    table["chosen"] = (table["strength"] == chosen.strength) & (
        table["l1_ratio"] == chosen.l1_ratio
    )
    return chosen, table
