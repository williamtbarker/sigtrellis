"""Evaluate frozen full-data model state on a truly untouched external cohort."""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.special import expit, softmax

from sigtrellis.config import load_config
from sigtrellis.domain import Dataset, FloatArray, array_hash, write_json
from sigtrellis.integrity import verify_run_artifacts
from sigtrellis.io import load_bulk
from sigtrellis.metrics import evaluate, group_weights
from sigtrellis.preprocessing import Normalizer
from sigtrellis.qc import group_values, validate_dataset


def frozen_predict(state: dict[str, Any], data: Dataset) -> FloatArray:
    genes = state["genes"]
    if set(genes) != set(data.expression.columns):
        raise ValueError("External expression must contain the exact training gene universe")
    x = data.expression.loc[:, genes]
    norm = Normalizer(
        state["normalization"],
        np.asarray(state["reference"], dtype=float) if state["reference"] is not None else None,
        genes,
        state.get("imputation", "reject"),
        np.asarray(state["fill_values"], dtype=float)
        if state.get("fill_values") is not None
        else None,
    )
    normalized = norm.transform(x)
    gene_x = normalized[:, state["selected_indices"]]
    gene_x = (gene_x - np.asarray(state["gene_center"])) / np.asarray(state["gene_scale"])
    cov: list[FloatArray] = []
    for name in state["covariates"]:
        if name in state["covariate_categories"]:
            categories = state["covariate_categories"][name]
            if not set(data.metadata[name].astype(str)) <= set(categories):
                raise ValueError(f"External covariate {name} contains unseen categories")
            cov.extend(
                (data.metadata[name].astype(str) == level).to_numpy(dtype=float)
                for level in categories[1:]
            )
        else:
            cov.append(data.metadata[name].to_numpy(dtype=float))
    cov_x = np.column_stack(cov) if cov else np.empty((len(x), 0))
    if len(cov):
        cov_x = (cov_x - np.asarray(state["covariate_center"])) / np.asarray(
            state["covariate_scale"]
        )
    values = np.column_stack([gene_x, cov_x])
    if state["covariates_only"]:
        values = cov_x
    if state["dummy"]:
        if state["outcome_type"] == "continuous":
            return np.full(len(x), float(np.asarray(state["constant"]).ravel()[0]))
        return np.tile(state["class_prior"], (len(x), 1))
    linear = values @ np.atleast_2d(
        np.asarray(state["raw_coefficients"], dtype=float)
    ).T + np.asarray(state["intercept"])
    if state["outcome_type"] == "continuous":
        return np.asarray(linear.ravel(), dtype=float)
    if state["n_classes"] == 2:
        positive = expit(linear.ravel())
        return np.column_stack([1 - positive, positive])
    return np.asarray(softmax(linear, axis=1), dtype=float)


def external_validate(
    run: Path,
    expression: Path | None,
    metadata: Path | None,
    output: Path,
    orientation: str = "auto",
    *,
    single_cell: Path | None = None,
    model_kind: str = "full",
) -> dict[str, Any]:
    if output.exists() and any(output.iterdir()):
        raise ValueError("External output directory must be empty")
    manifest = json.loads((run / "run_manifest.json").read_text())
    if manifest["status"] != "complete":
        raise ValueError("Training run did not complete")
    if model_kind not in {"full", "panel"}:
        raise ValueError("model_kind must be full or panel")
    state_name = "model_state.json" if model_kind == "full" else "panel_state.json"
    verified = verify_run_artifacts(
        run, manifest, ("configuration.json", state_name, "sample_metadata.csv")
    )
    config = load_config(run / "configuration.json")
    if single_cell is not None:
        if expression is not None or metadata is not None:
            raise ValueError("Supply a single-cell input or an expression/metadata pair")
        if config.single_cell_mode == "distribution":
            from sigtrellis.cell_features import distribution_features

            data = distribution_features(single_cell, replace(config, cell_resamples=0))
        else:
            from sigtrellis.singlecell import pseudobulk

            datasets = pseudobulk(
                single_cell, replace(config, cell_type_value=manifest["cell_type"])
            )
            if len(datasets) != 1:
                raise ValueError("Frozen external evaluation requires the single fitted cell type")
            data = datasets[0]
    else:
        if expression is None or metadata is None:
            raise ValueError("Supply both expression and metadata")
        data = load_bulk(expression, metadata, config, orientation)
    validate_dataset(data, config)
    train_metadata = pd.read_csv(
        run / "sample_metadata.csv",
        dtype={config.sample_id: str, config.group or config.sample_id: str},
        keep_default_na=False,
    )
    if set(data.expression.index) & set(train_metadata[config.sample_id].astype(str)):
        raise ValueError("External sample identifiers overlap training")
    if set(group_values(data, config)) & set(
        train_metadata[config.group or config.sample_id].astype(str)
    ):
        raise ValueError("External biological groups overlap training")
    state = json.loads((run / state_name).read_text())
    if set(state["genes"]) != set(data.expression.columns):
        raise ValueError("External expression must contain the exact training gene universe")
    training_hashes = set(manifest["sample_expression_hashes"].values())
    hashes_to_check = [
        array_hash(row) for row in data.expression.loc[:, state["genes"]].to_numpy(dtype=float)
    ]
    if config.input_scale == "features":
        training_hashes = set(manifest.get("sample_cell_hashes", {}).values())
        hashes_to_check = list(data.upstream_qc.get("sample_cell_hashes", {}).values())
    for profile_hash in hashes_to_check:
        if profile_hash in training_hashes:
            raise ValueError(
                "External cohort contains an exact training expression profile, even under a different ID"
            )
    prediction = frozen_predict(state, data)
    if config.outcome_type == "continuous":
        y = data.metadata[config.outcome].to_numpy(dtype=float)
    else:
        mapping = {label: i for i, label in enumerate(manifest["classes"])}
        labels = data.metadata[config.outcome].astype(str)
        if not set(labels) <= set(mapping):
            raise ValueError("Unknown external outcome class")
        y = labels.map(mapping).to_numpy(dtype=float)
    metrics = evaluate(y, prediction, config, group_weights(group_values(data, config)))
    output.mkdir(parents=True, exist_ok=True)
    columns = (
        ["prediction"]
        if config.outcome_type == "continuous"
        else [f"probability:{label}" for label in manifest["classes"]]
    )
    frame = pd.DataFrame(prediction, index=data.expression.index, columns=columns)
    frame.to_csv(output / "external_predictions.csv", index_label="sample_id")
    record = {
        "metrics": metrics,
        "input_hashes": data.input_hashes,
        "training_expression_hash": manifest["expression_hash"],
        "training_model_sha256": verified[state_name],
        "verified_training_artifacts": verified,
        "classes": manifest["classes"],
        "scope": "Frozen final compact panel"
        if model_kind == "panel"
        else "Frozen final full-data model; does not validate a separately refitted consensus panel",
        "model_kind": model_kind,
        "cell_profile_duplicate_check": bool(hashes_to_check)
        if config.input_scale == "features"
        else None,
        "training_transform_refit": False,
    }
    write_json(output / "external_metrics.json", record)
    return record
