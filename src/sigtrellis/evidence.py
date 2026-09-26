"""Join scientific identities and supporting tests without altering model inputs."""

from __future__ import annotations

from dataclasses import asdict, replace
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

from sigtrellis.config import Config
from sigtrellis.de import differential_expression
from sigtrellis.domain import Audit, Dataset, FloatArray, write_json


def annotate_features(table: pd.DataFrame, data: Dataset) -> pd.DataFrame:
    result = table.copy()
    result["feature_id"] = result.gene_id
    if not data.features:
        result["source_gene_id"] = result.gene_id
        result["feature_kind"] = "gene_expression"
        result["feature_label"] = result.gene_id
        result["cell_type"] = data.cell_type or "bulk"
        return result
    result["source_gene_id"] = result.gene_id.map(
        {k: v.gene_id or "" for k, v in data.features.items()}
    )
    result["feature_kind"] = result.gene_id.map({k: v.kind for k, v in data.features.items()})
    result["cell_type"] = result.gene_id.map({k: v.cell_type for k, v in data.features.items()})
    result["program"] = result.gene_id.map({k: v.program or "" for k, v in data.features.items()})
    result["feature_unit"] = result.gene_id.map({k: v.unit for k, v in data.features.items()})
    result["gene_partner"] = result.gene_id.map(
        {k: v.gene_partner or "" for k, v in data.features.items()}
    )
    result["program_partner"] = result.gene_id.map(
        {k: v.program_partner or "" for k, v in data.features.items()}
    )
    result["threshold"] = result.gene_id.map({k: v.threshold for k, v in data.features.items()})
    result["feature_label"] = result.gene_id.map({k: v.label for k, v in data.features.items()})
    return result


def public_table(table: pd.DataFrame, data: Dataset) -> pd.DataFrame:
    """Use gene_id for genes and feature_id for typed distribution predictors."""
    if not data.features:
        return table
    return table.drop(columns=["gene_id"]).rename(columns={"source_gene_id": "gene_id"})


def supporting_evidence(
    data: Dataset, y: FloatArray, config: Config, audit: Audit, output: Path
) -> pd.DataFrame | None:
    frames: list[pd.DataFrame] = []
    designs: dict[str, Any] = {}
    if not data.counts_by_cell_type:
        sources = {data.cell_type or "bulk": data}
    else:
        sources = {}
        for state, counts in data.counts_by_cell_type.items():
            numbers = data.upstream_qc["cells_per_state"][state]
            rows = np.flatnonzero(
                np.array([numbers[str(s)] for s in counts.index]) >= config.min_cells
            ).astype(np.int64)
            sources[state] = Dataset(counts.iloc[rows], data.metadata.iloc[rows], state)
    for state, source in sources.items():
        try:
            local = replace(
                config,
                input_scale="counts",
                normalization="logcpm",
                imputation="reject",
                single_cell_mode="pseudobulk",
                stability_normalizations=(),
                cell_resamples=0,
            )
            positions = data.expression.index.get_indexer(source.expression.index)
            result = differential_expression(source, y[positions], local)
            designs[state] = result.attrs
            audit.warnings.extend(
                f"EXPLORATORY_DE [{state}]: {w}" for w in result.attrs["warnings"]
            )
            result = result.copy()
            result["cell_type"] = state
            frames.append(result)
        except ValueError as exc:
            audit.warnings.append(f"EXPLORATORY_DE_UNAVAILABLE [{state}]: {exc}")
            designs[state] = {"status": "unavailable", "reason": str(exc)}
    write_json(output / "de_design.json", designs)
    if not frames:
        return None
    result = pd.concat(frames)
    if len(frames) > 1:
        valid = result.de_pvalue.notna().to_numpy()
        result["de_adjusted_pvalue"] = np.nan
        result.loc[valid, "de_adjusted_pvalue"] = false_discovery_control(
            result.loc[valid, "de_pvalue"].to_numpy()
        )
    result.to_csv(output / "differential_expression.csv")
    return result


def write_feature_contract(data: Dataset, output: Path) -> None:
    write_json(
        output / "feature_schema.json", {name: asdict(spec) for name, spec in data.features.items()}
    )
    if data.features:
        data.expression.to_csv(output / "sample_features.tsv.gz", sep="\t", index_label="sample_id")
