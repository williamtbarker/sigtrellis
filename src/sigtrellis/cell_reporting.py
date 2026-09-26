"""Show biological-replicate distributions without pooling cells into fake replicates."""

from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sigtrellis.config import Config
from sigtrellis.domain import Dataset


def cell_evidence(data: Dataset, table: pd.DataFrame, config: Config, output: Path) -> bool:
    if not data.features:
        return False
    if data.upstream_qc.get("sample_cell_quality"):
        quality = pd.DataFrame.from_dict(data.upstream_qc["sample_cell_quality"], orient="index")
        quality["outcome"] = data.metadata[config.outcome]
        quality["group_id"] = data.metadata[config.group or config.sample_id]
        quality["state_counts"] = [
            str(
                {
                    state: counts[str(sample)]
                    for state, counts in data.upstream_qc["cells_per_state"].items()
                }
            )
            for sample in quality.index
        ]
        quality.to_csv(output / "cell_quality_by_sample.csv", index_label="sample_id")
    selected = table.gene_id.drop_duplicates().head(6).tolist()
    rows = []
    fig, axes = plt.subplots(len(selected), 1, figsize=(10, 3.2 * len(selected)), squeeze=False)
    categories = sorted(data.metadata[config.outcome].astype(str).unique())
    groups = data.metadata[config.group or config.sample_id].astype(str)
    rng = np.random.default_rng(config.seed)
    for axis, feature in zip(axes.ravel(), selected, strict=True):
        spec = data.features[feature]
        values = data.expression[feature].to_numpy(dtype=float)
        if config.outcome_type == "continuous":
            axis.scatter(data.metadata[config.outcome].to_numpy(dtype=float), values, alpha=0.6)
            axis.set_xlabel(config.outcome)
        else:
            labels = data.metadata[config.outcome].astype(str).to_numpy()
            for i, label in enumerate(categories):
                mask = labels == label
                axis.scatter(
                    i + rng.uniform(-0.13, 0.13, mask.sum()), values[mask], alpha=0.6, s=16
                )
            axis.set_xticks(range(len(categories)), [textwrap.fill(v, 25) for v in categories])
        axis.set_title(textwrap.fill(spec.label, 85), fontsize=10)
        axis.set_ylabel(textwrap.fill(spec.unit, 24), fontsize=8)
        for sample, group, value in zip(data.expression.index, groups, values, strict=True):
            rows.append(
                {
                    "sample_id": sample,
                    "group_id": group,
                    "feature_id": feature,
                    "feature_label": spec.label,
                    "observed_value": value,
                    "outcome": data.metadata.loc[sample, config.outcome],
                    "batch": data.metadata.loc[sample, config.batch] if config.batch else None,
                    "state_cells": data.upstream_qc.get("cells_per_state", {})
                    .get(spec.cell_type, {})
                    .get(str(sample)),
                    "feature_unit": spec.unit,
                }
            )
    pd.DataFrame(rows).to_csv(output / "cell_distribution_evidence.csv", index=False)
    fig.suptitle(
        "Cell distributions: each point is a specimen; missing values are omitted", fontsize=11
    )
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(output / "figures/cell_distributions.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    return True
