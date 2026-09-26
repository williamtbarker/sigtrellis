"""Bounded, descriptive correlation groups and coefficient-substitution evidence."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Any

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from sigtrellis.config import Config
from sigtrellis.domain import Dataset, FloatArray
from sigtrellis.preprocessing import Normalizer


@dataclass
class CorrelationResult:
    groups: pd.DataFrame
    substitutions: pd.DataFrame
    matrix: pd.DataFrame
    membership: dict[str, str]


def correlation_diagnostics(
    data: Dataset,
    coefficients: FloatArray,
    ranking: pd.DataFrame,
    config: Config,
    *,
    contrasts: list[str] | None = None,
) -> CorrelationResult:
    if contrasts is None:
        if coefficients.shape[1] != 1:
            raise ValueError("Explicit contrast order is required for multiclass coefficients")
        contrasts = ranking.contrast.drop_duplicates().tolist()
    if len(contrasts) != coefficients.shape[1] or len(set(contrasts)) != len(contrasts):
        raise ValueError("Contrast labels do not match coefficient axes")
    ranked = ranking.sort_values(
        ["selection_frequency", "sign_consistency", "rank_median_selected", "gene_id"],
        ascending=[False, False, True, True],
        kind="stable",
    )
    genes = ranked.loc[ranked.selection_frequency > 0, "gene_id"].drop_duplicates().tolist()
    genes = genes[: config.correlation_max_features]
    normalized = Normalizer(config.normalization).fit(data.expression).transform(data.expression)
    all_genes = list(data.expression.columns)
    indices = [all_genes.index(g) for g in genes]
    variable = [i for i in indices if normalized[:, i].std() > 1e-12]
    genes = [all_genes[i] for i in variable]
    if not genes:
        return CorrelationResult(
            pd.DataFrame(columns=["group_id", "contrast", "members", "selection_frequency"]),
            pd.DataFrame(columns=["contrast", "gene_a", "gene_b", "exclusive_given_any"]),
            pd.DataFrame(),
            {},
        )
    corr = (
        np.atleast_2d(np.corrcoef(normalized[:, variable], rowvar=False))
        if len(genes) > 1
        else np.ones((1, 1))
    )
    corr = np.clip(corr, -1, 1)
    distance = np.clip(1 - abs(corr), 0, 1)
    np.fill_diagonal(distance, 0)
    labels = (
        fcluster(
            linkage(squareform(distance, checks=False), method="complete"),
            1 - config.correlation_threshold,
            criterion="distance",
        )
        if len(genes) > 1
        else np.array([1])
    )
    # Preserve the contrast axis: evidence for class A must not inflate class B.
    selected = abs(coefficients) > config.coefficient_tolerance
    group_rows: list[dict[str, Any]] = []
    pair_rows: list[dict[str, Any]] = []
    membership: dict[str, str] = {}
    for label in sorted(set(labels)):
        positions = np.flatnonzero(labels == label)
        members = [genes[i] for i in positions]
        source_indices = [all_genes.index(g) for g in members]
        group_id = f"correlation_{int(label):03d}"
        membership.update({g: group_id for g in members})
        for ci, contrast in enumerate(contrasts):
            contrast_selected = selected[:, ci, :]
            group_rows.append(
                {
                    "group_id": group_id,
                    "contrast": contrast,
                    "members": ";".join(members),
                    "n_members": len(members),
                    "selection_frequency": float(
                        contrast_selected[:, source_indices].any(axis=1).mean()
                    ),
                    "minimum_member_frequency": float(
                        contrast_selected[:, source_indices].mean(axis=0).min()
                    ),
                    "interpretation": "Contrast-specific selection in a descriptive marginal coexpression group",
                }
            )
            for a, b in combinations(positions, 2):
                sa = contrast_selected[:, all_genes.index(genes[a])]
                sb = contrast_selected[:, all_genes.index(genes[b])]
                union = sa | sb
                pair_rows.append(
                    {
                        "group_id": group_id,
                        "contrast": contrast,
                        "gene_a": genes[a],
                        "gene_b": genes[b],
                        "expression_correlation": float(corr[a, b]),
                        "frequency_any": float(union.mean()),
                        "frequency_both": float((sa & sb).mean()),
                        "exclusive_given_any": float((sa ^ sb).sum() / union.sum())
                        if union.any()
                        else None,
                        "selection_phi": float(np.corrcoef(sa, sb)[0, 1])
                        if sa.std() and sb.std()
                        else None,
                    }
                )
    return CorrelationResult(
        pd.DataFrame(group_rows),
        pd.DataFrame(
            pair_rows,
            columns=[
                "group_id",
                "contrast",
                "gene_a",
                "gene_b",
                "expression_correlation",
                "frequency_any",
                "frequency_both",
                "exclusive_given_any",
                "selection_phi",
            ],
        ),
        pd.DataFrame(corr, index=genes, columns=genes),
        membership,
    )
