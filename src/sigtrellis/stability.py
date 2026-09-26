"""Empirical selection stability, not a claim of formal PFER control."""

from __future__ import annotations

from dataclasses import dataclass, replace
from itertools import combinations
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from sigtrellis.config import Config
from sigtrellis.domain import Audit, Dataset, FloatArray
from sigtrellis.modeling import fit_model, tune
from sigtrellis.splits import subsample_groups


@dataclass
class StabilityResult:
    table: pd.DataFrame
    coefficients: FloatArray
    tuning: pd.DataFrame
    resamples: list[dict[str, Any]]
    summary: dict[str, Any]


def jaccard_scores(selected: np.ndarray[Any, Any]) -> FloatArray:
    scores: list[float] = []
    for a, b in combinations(selected, 2):
        union = int(np.count_nonzero(a | b))
        # Empty/empty is not evidence that a biomarker signature is reproducible.
        scores.append(float(np.count_nonzero(a & b) / union) if union else float("nan"))
    return np.asarray(scores, dtype=float)


def summarize_coefficients(
    coefficients: FloatArray, genes: list[str], contrasts: list[str], tolerance: float
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for contrast_index, contrast in enumerate(contrasts):
        coefs = coefficients[:, contrast_index, :]
        selected = abs(coefs) > tolerance
        ranks = np.asarray([rankdata(-abs(row), method="average") for row in coefs])
        for j, gene in enumerate(genes):
            nonzero = coefs[selected[:, j], j]
            selected_ranks = ranks[selected[:, j], j]
            n = len(nonzero)
            q1, q3 = np.quantile(coefs[:, j], [0.25, 0.75])
            rows.append(
                {
                    "gene_id": gene,
                    "contrast": contrast,
                    "selection_frequency": float(selected[:, j].mean()),
                    "sign_consistency": float(max((nonzero > 0).sum(), (nonzero < 0).sum()) / n)
                    if n
                    else 0.0,
                    "coefficient_mean": float(coefs[:, j].mean()),
                    "coefficient_median": float(np.median(coefs[:, j])),
                    "coefficient_median_selected": float(np.median(nonzero)) if n else 0.0,
                    "coefficient_q25": float(q1),
                    "coefficient_q75": float(q3),
                    "rank_median_selected": float(np.median(selected_ranks)) if n else float("nan"),
                    "rank_iqr_selected": float(
                        np.diff(np.quantile(selected_ranks, [0.25, 0.75]))[0]
                    )
                    if n
                    else float("nan"),
                    "n_selecting_fits": n,
                }
            )
    return pd.DataFrame(rows)


def stability_select(
    data: Dataset, y: FloatArray, config: Config, contrasts: list[str], audit: Audit
) -> StabilityResult:
    rng = np.random.default_rng(config.seed + 101)
    values: list[FloatArray] = []
    tuning: list[pd.DataFrame] = []
    resamples: list[dict[str, Any]] = []
    normalizations = config.stability_normalizations or (config.normalization,)
    for index in range(config.stability_resamples):
        rows = subsample_groups(data, y, config, rng)
        normalization = normalizations[index % len(normalizations)]
        local = replace(config, normalization=normalization)
        subset = data.subset(rows)
        context = f"stability:{index}/{normalization}"
        params, table = tune(subset, y[rows], local, audit, context)
        model = fit_model(subset, y[rows], local, params, audit, context + "/refit")
        values.append(model.coefficients())
        tuning.append(table)
        resamples.append(
            {
                "index": index,
                "sample_ids": subset.expression.index.tolist(),
                "normalization": normalization,
                "strength": params.strength,
                "l1_ratio": params.l1_ratio,
            }
        )
    coefficients = np.asarray(values, dtype=float)
    table = summarize_coefficients(
        coefficients, list(data.expression.columns), contrasts, config.coefficient_tolerance
    )
    for norm in normalizations:
        mask = np.array([r["normalization"] == norm for r in resamples])
        evidence = summarize_coefficients(
            coefficients[mask],
            list(data.expression.columns),
            contrasts,
            config.coefficient_tolerance,
        )
        table[f"frequency_{norm}"] = evidence["selection_frequency"]
    any_selected = np.asarray((abs(coefficients) > config.coefficient_tolerance).any(axis=1))
    jaccard = jaccard_scores(any_selected)
    finite = jaccard[np.isfinite(jaccard)]
    summary = {
        "resamples": len(values),
        "subsample_fraction": config.stability_fraction,
        "normalizations": normalizations,
        "n_selected_per_fit": any_selected.sum(axis=1),
        "median_pairwise_jaccard": float(np.median(finite)) if len(finite) else None,
        "empty_empty_comparisons": int((~np.isfinite(jaccard)).sum()),
        "interpretation": "Descriptive overlapping group subsamples with inner retuning; no PFER/FDR guarantee",
    }
    return StabilityResult(
        table, coefficients, pd.concat(tuning, ignore_index=True), resamples, summary
    )
