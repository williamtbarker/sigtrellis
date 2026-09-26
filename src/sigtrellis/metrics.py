"""Proper-loss tuning and complementary held-out performance measures."""

from __future__ import annotations

from typing import Any

import numpy as np
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    log_loss,
    matthews_corrcoef,
    mean_absolute_error,
    mean_squared_error,
    r2_score,
    roc_auc_score,
)

from sigtrellis.config import Config
from sigtrellis.domain import FloatArray


def group_weights(groups: np.ndarray[Any, Any]) -> FloatArray:
    _, inverse, counts = np.unique(groups, return_inverse=True, return_counts=True)
    weights = 1 / counts[inverse].astype(float)
    return np.asarray(weights / weights.mean(), dtype=np.float64)


def loss(y: FloatArray, prediction: FloatArray, config: Config, weights: FloatArray) -> float:
    if config.outcome_type == "continuous":
        return float(mean_squared_error(y, prediction, sample_weight=weights))
    return float(
        log_loss(y, prediction, labels=np.arange(prediction.shape[1]), sample_weight=weights)
    )


def calibration(
    y: FloatArray, probabilities: FloatArray, weights: FloatArray, bins: int = 5
) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    assignment = np.minimum((probabilities * bins).astype(int), bins - 1)
    for bin_id in range(bins):
        mask = assignment == bin_id
        if mask.any():
            rows.append(
                {
                    "mean_probability": float(
                        np.average(probabilities[mask], weights=weights[mask])
                    ),
                    "observed_fraction": float(np.average(y[mask], weights=weights[mask])),
                    "n_rows": int(mask.sum()),
                    "weight": float(weights[mask].sum()),
                }
            )
    return rows


def evaluate(
    y: FloatArray, prediction: FloatArray, config: Config, weights: FloatArray
) -> dict[str, Any]:
    if config.outcome_type == "continuous":
        return {
            "rmse": float(np.sqrt(mean_squared_error(y, prediction, sample_weight=weights))),
            "mae": float(mean_absolute_error(y, prediction, sample_weight=weights)),
            "r2": float(r2_score(y, prediction, sample_weight=weights)) if len(y) > 1 else None,
        }
    binary = config.outcome_type == "binary"
    hard = (
        (prediction[:, 1] >= config.decision_threshold).astype(int)
        if binary
        else prediction.argmax(axis=1)
    )
    metrics: dict[str, Any] = {
        "log_loss": loss(y, prediction, config, weights),
        "balanced_accuracy": float(balanced_accuracy_score(y, hard, sample_weight=weights)),
        "mcc": float(matthews_corrcoef(y, hard, sample_weight=weights)),
    }
    if binary:
        prob = prediction[:, 1]
        tn, fp, fn, tp = confusion_matrix(y, hard, labels=[0, 1], sample_weight=weights).ravel()
        both = len(np.unique(y)) == 2
        metrics.update(
            {
                "roc_auc": float(roc_auc_score(y, prob, sample_weight=weights)) if both else None,
                "pr_auc_average_precision": float(
                    average_precision_score(y, prob, sample_weight=weights)
                )
                if both
                else None,
                "sensitivity": float(tp / (tp + fn)) if tp + fn else None,
                "specificity": float(tn / (tn + fp)) if tn + fp else None,
                "brier_score": float(brier_score_loss(y, prob, sample_weight=weights)),
                "prevalence": float(np.average(y, weights=weights)),
                "calibration": calibration(y, prob, weights),
            }
        )
    else:
        k = prediction.shape[1]
        onehot = np.eye(k)[y.astype(int)]
        all_classes = len(np.unique(y)) == k
        metrics.update(
            {
                "roc_auc_macro_ovr": float(
                    roc_auc_score(
                        y, prediction, multi_class="ovr", average="macro", sample_weight=weights
                    )
                )
                if all_classes
                else None,
                "pr_auc_macro_ovr": float(
                    average_precision_score(
                        onehot, prediction, average="macro", sample_weight=weights
                    )
                )
                if all_classes
                else None,
                "brier_multiclass_sum": float(
                    np.average(((prediction - onehot) ** 2).sum(axis=1), weights=weights)
                ),
                "calibration_by_class": {
                    str(i): calibration(onehot[:, i], prediction[:, i], weights) for i in range(k)
                },
            }
        )
    return metrics
