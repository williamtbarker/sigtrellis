"""Nested grouped validation and full-procedure permutation controls."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

import numpy as np
import pandas as pd

from sigtrellis.config import Config
from sigtrellis.domain import Audit, Dataset, FloatArray
from sigtrellis.metrics import evaluate, group_weights, loss
from sigtrellis.modeling import Hyperparameters, fit_model, fit_prepared, tune
from sigtrellis.qc import group_values
from sigtrellis.splits import make_splits, permute_outcome, record_split


@dataclass
class ValidationResult:
    predictions: pd.DataFrame
    tuning: pd.DataFrame
    coefficients: FloatArray
    metrics: dict[str, Any]


def nested_validate(
    data: Dataset, y: FloatArray, config: Config, audit: Audit, context: str = "evaluation"
) -> ValidationResult:
    n_classes = 0 if config.outcome_type == "continuous" else len(np.unique(y))
    prediction_rows: list[dict[str, Any]] = []
    all_tuning: list[pd.DataFrame] = []
    coefficients: list[FloatArray] = []
    fold_metrics: list[dict[str, Any]] = []
    for split in make_splits(data, y, config):
        label = f"{context}/repeat:{split.repeat}/outer:{split.fold}"
        record_split(audit, data, split, label, config)
        train, test = data.subset(split.train), data.subset(split.test)
        parameters, table = tune(train, y[split.train], config, audit, label)
        all_tuning.append(table)
        model = fit_model(train, y[split.train], config, parameters, audit, label + "/refit")
        pred = model.predict(test)
        coefficients.append(model.coefficients())
        if config.covariates:
            baseline_param, baseline_tuning = tune(
                train,
                y[split.train],
                config,
                audit,
                label + "/covariate_baseline",
                covariates_only=True,
            )
            all_tuning.append(baseline_tuning)
        else:
            baseline_param = Hyperparameters(config.strengths[0], config.l1_ratios[0])
        baseline_model = fit_prepared(
            model.prepared, train, y[split.train], baseline_param, covariates_only=True
        )
        baseline = baseline_model.predict(test)
        w = group_weights(group_values(test, config))
        measured = evaluate(y[split.test], pred, config, w)
        fold_metrics.append(
            {
                "context": label,
                "n_test_groups": len(np.unique(group_values(test, config))),
                "metrics": measured,
                "model_loss": loss(y[split.test], pred, config, w),
                "baseline_loss": loss(y[split.test], baseline, config, w),
                "parameters": {"strength": parameters.strength, "l1_ratio": parameters.l1_ratio},
            }
        )
        for local, row in enumerate(split.test):
            record: dict[str, Any] = {
                "sample_id": str(data.expression.index[row]),
                "group_id": str(group_values(data, config)[row]),
                "repeat": split.repeat,
                "fold": split.fold,
                "observed": float(y[row]),
            }
            if not n_classes:
                record.update(prediction=float(pred[local]), baseline=float(baseline[local]))
            else:
                record.update({f"p_{c}": float(pred[local, c]) for c in range(n_classes)})
                record.update(
                    {f"baseline_{c}": float(baseline[local, c]) for c in range(n_classes)}
                )
            prediction_rows.append(record)
    frame = pd.DataFrame(prediction_rows)
    # Pool repeated predictions as repeated OOF rows with equal total weight per
    # biological group. They are NOT independent replications or confidence intervals.
    observed = frame["observed"].to_numpy(dtype=float)
    pred_all = (
        frame[[f"p_{i}" for i in range(n_classes)]].to_numpy()
        if n_classes
        else frame["prediction"].to_numpy()
    )
    base_all = (
        frame[[f"baseline_{i}" for i in range(n_classes)]].to_numpy()
        if n_classes
        else frame["baseline"].to_numpy()
    )
    weights = group_weights(frame["group_id"].to_numpy())
    model_loss, baseline_loss = (
        loss(observed, pred_all, config, weights),
        loss(observed, base_all, config, weights),
    )
    metrics: dict[str, Any] = {
        "out_of_fold": evaluate(observed, pred_all, config, weights),
        "baseline": evaluate(observed, base_all, config, weights),
        "model_loss": model_loss,
        "baseline_loss": baseline_loss,
        "loss_improvement": baseline_loss - model_loss,
        "folds": fold_metrics,
        "estimand": "Equal weight per biological group; prediction of a sample from an unseen group",
        "performance_scope": "Entire training/tuning procedure, not the post-hoc consensus panel",
        "uncertainty": "Fold dispersion is descriptive; overlapping CV folds are not independent",
    }
    return ValidationResult(
        frame,
        pd.concat(all_tuning, ignore_index=True),
        np.asarray(coefficients, dtype=float),
        metrics,
    )


def permutation_control(
    data: Dataset, y: FloatArray, config: Config, observed: float, audit: Audit
) -> dict[str, Any]:
    if not config.permutations:
        return {"status": "not_run", "pvalue": None, "null_improvements": []}
    rng = np.random.default_rng(config.seed + 7919)
    null: list[float] = []
    permutation_audits: list[dict[str, Any]] = []
    for index in range(config.permutations):
        permuted = permute_outcome(data, y, config, rng)
        local_audit = Audit()
        # Label-dependent stratification is rebuilt using permuted labels. Reusing
        # observed-label stratification would condition the null on leaked labels.
        result = nested_validate(
            data,
            permuted,
            replace(config, permutations=0),
            local_audit,
            context=f"permutation:{index}",
        )
        null.append(float(result.metrics["loss_improvement"]))
        permutation_audits.append(
            {
                "index": index,
                "labels": permuted.tolist(),
                "splits": local_audit.splits,
                "fit_count": len(local_audit.fits),
                "fits": [
                    {key: value for key, value in fit.items() if key != "removed_genes"}
                    for fit in local_audit.fits
                ],
                "filtering_note": "Excluded genes are the complement of retained_genes in the declared input universe",
                "tuning_results": result.tuning.to_dict("records"),
                "outer_parameters": [
                    {"context": fold["context"], "parameters": fold["parameters"]}
                    for fold in result.metrics["folds"]
                ],
            }
        )
        audit.warnings.extend(local_audit.warnings)
    pvalue = (1 + sum(value >= observed - 1e-12 for value in null)) / (1 + len(null))
    # A global label-permutation test is not a conditional gene test given covariates.
    status = "global_only_with_covariates" if config.covariates else "tested"
    return {
        "status": status,
        "pvalue": pvalue,
        "null_improvements": null,
        "scheme": config.permutation_scheme,
        "strata": config.permutation_strata,
        "statistic": "baseline minus model held-out proper loss; larger is better",
        "exchangeability": "Must be justified by study design; no automatic observational causal inference",
        "resolution": 1 / (1 + len(null)),
        "audits": permutation_audits,
    }
