"""Strict, serializable run contracts. Unknown keys are errors."""

from __future__ import annotations

import json
import math
import tomllib
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Any, Literal

import yaml

Outcome = Literal["binary", "multiclass", "continuous"]


@dataclass(frozen=True)
class Config:
    outcome: str = "phenotype"
    outcome_type: Outcome = "binary"
    sample_id: str = "sample_id"
    group: str | None = None
    batch: str | None = None
    covariates: tuple[str, ...] = ()
    positive_class: str | None = None
    input_scale: str = "counts"
    normalization: str = "logcpm"
    candidate_method: str = "variance"
    max_features: int = 1000
    min_count: float = 5.0
    min_prevalence: float = 0.2
    de_fdr: float = 0.1
    supporting_de: bool = False
    de_pair_group: bool = False
    outer_folds: int = 5
    inner_folds: int = 3
    repeats: int = 1
    cv_strategy: str = "grouped"
    seed: int = 42
    strengths: tuple[float, ...] = (0.01, 0.05, 0.2)
    l1_ratios: tuple[float, ...] = (0.2, 0.5, 0.9)
    max_iter: int = 10000
    tolerance: float = 0.0001
    tuning_rule: str = "one_se"
    coefficient_tolerance: float = 1e-7
    decision_threshold: float = 0.5
    stability_resamples: int = 50
    stability_fraction: float = 0.75
    stability_normalizations: tuple[str, ...] = ()
    selection_threshold: float = 0.8
    sign_threshold: float = 0.9
    outer_selection_threshold: float = 0.5
    permutations: int = 99
    permutation_scheme: str = "group"
    permutation_strata: str | None = None
    permutation_alpha: float = 0.05
    min_groups_gate: int = 20
    correlation_threshold: float = 0.85
    correlation_max_features: int = 100
    assess_batches: bool = True
    cell_type: str | None = None
    cell_type_value: str | None = None
    layer: str | None = None
    min_cells: int = 20
    cell_min_counts: int = 1
    cell_min_genes: int = 1
    mitochondrial_prefix: str | None = None
    max_mito_fraction: float = 1.0
    chunk_size: int = 2048
    max_dense_mb: int = 2048

    def validate(self) -> None:
        if not isinstance(self.outcome, str) or not isinstance(self.sample_id, str):
            raise ValueError("outcome and sample_id must be nonempty strings")
        for key in (
            "outcome",
            "sample_id",
            "group",
            "batch",
            "positive_class",
            "cell_type",
            "cell_type_value",
            "layer",
            "permutation_strata",
            "mitochondrial_prefix",
        ):
            value = getattr(self, key)
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise ValueError(f"{key} must be a nonempty string")
        for key in ("supporting_de", "de_pair_group", "assess_batches"):
            if not isinstance(getattr(self, key), bool):
                raise ValueError(f"{key} must be a boolean")
        choices = {
            "outcome_type": {"binary", "multiclass", "continuous"},
            "input_scale": {"counts", "log_expression"},
            "normalization": {"logcpm", "median_ratio", "none"},
            "candidate_method": {"variance", "association", "deseq2", "none"},
            "cv_strategy": {"grouped", "leave_group_out", "leave_batch_out"},
            "permutation_scheme": {"group", "within_group"},
            "tuning_rule": {"one_se", "minimum_loss"},
        }
        for key, allowed in choices.items():
            if getattr(self, key) not in allowed:
                raise ValueError(f"{key} must be one of {sorted(allowed)}")
        protected = {self.outcome, self.sample_id, self.group, self.batch}
        if any(c in protected for c in self.covariates):
            raise ValueError("Outcome, sample/group identifiers, and batch cannot be predictors")
        if len(set(self.covariates)) != len(self.covariates):
            raise ValueError("Duplicate covariates")
        if any(not isinstance(c, str) or not c.strip() for c in self.covariates):
            raise ValueError("Covariates must be nonempty column names")
        if self.outcome in {self.sample_id, self.group}:
            raise ValueError("Outcome cannot be a sample or group identifier")
        if self.input_scale == "log_expression" and self.normalization != "none":
            raise ValueError("Declared log_expression requires normalization=none")
        if self.input_scale == "counts" and self.normalization == "none":
            raise ValueError("Counts require logcpm or median_ratio normalization")
        if self.candidate_method == "deseq2" or self.supporting_de:
            if self.input_scale != "counts" or self.outcome_type != "binary":
                raise ValueError("PyDESeq2 integration currently requires binary raw-count data")
        if self.de_pair_group and not self.group:
            raise ValueError("Paired DE requires an explicit group column")
        if self.cv_strategy == "leave_batch_out" and not self.batch:
            raise ValueError("leave_batch_out requires batch metadata")
        for key in ("outer_folds", "inner_folds", "stability_resamples", "min_groups_gate"):
            if not isinstance(getattr(self, key), int) or getattr(self, key) < 2:
                raise ValueError(f"{key} must be an integer >= 2")
        for key in (
            "repeats",
            "max_iter",
            "max_features",
            "min_cells",
            "chunk_size",
            "correlation_max_features",
            "max_dense_mb",
            "cell_min_counts",
            "cell_min_genes",
        ):
            if not isinstance(getattr(self, key), int) or getattr(self, key) < 1:
                raise ValueError(f"{key} must be a positive integer")
        for key in ("seed", "permutations"):
            if not isinstance(getattr(self, key), int) or getattr(self, key) < 0:
                raise ValueError(f"{key} must be a nonnegative integer")
        for key in (
            "min_prevalence",
            "selection_threshold",
            "sign_threshold",
            "outer_selection_threshold",
            "de_fdr",
            "permutation_alpha",
            "correlation_threshold",
            "max_mito_fraction",
        ):
            if (
                not isinstance(getattr(self, key), (float, int))
                or not math.isfinite(getattr(self, key))
                or not 0 < getattr(self, key) <= 1
            ):
                raise ValueError(f"{key} must be in (0,1]")
        for key in ("stability_fraction", "decision_threshold"):
            if (
                not isinstance(getattr(self, key), (float, int))
                or not math.isfinite(getattr(self, key))
                or not 0 < getattr(self, key) < 1
            ):
                raise ValueError(f"{key} must be in (0,1)")
        for key in ("min_count", "tolerance", "coefficient_tolerance"):
            if (
                not isinstance(getattr(self, key), (float, int))
                or not math.isfinite(getattr(self, key))
                or getattr(self, key) <= 0
            ):
                raise ValueError(f"{key} must be finite and positive")
        if not self.strengths or any(
            not isinstance(x, (int, float)) or not math.isfinite(x) or x <= 0
            for x in self.strengths
        ):
            raise ValueError("strengths must be finite and positive")
        if not self.l1_ratios or any(
            not isinstance(x, (int, float)) or not math.isfinite(x) or not 0 < x <= 1
            for x in self.l1_ratios
        ):
            raise ValueError("l1_ratios must be in (0,1]; pure ridge is not sparse selection")
        for norm in self.stability_normalizations:
            if norm not in choices["normalization"]:
                raise ValueError("Invalid stability normalization")
            if (norm == "none") != (self.input_scale == "log_expression"):
                raise ValueError("Stability normalization is incompatible with input_scale")
        if len(set(self.stability_normalizations)) != len(self.stability_normalizations):
            raise ValueError("Duplicate stability normalizations")
        if len(self.stability_normalizations) > self.stability_resamples:
            raise ValueError("Every stability normalization needs at least one resample")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def load_config(path: Path | None, overrides: dict[str, Any] | None = None) -> Config:
    raw: dict[str, Any] = {}
    if path:
        if path.suffix == ".json":
            raw = json.loads(path.read_text())
        elif path.suffix == ".toml":
            raw = tomllib.loads(path.read_text())
        else:
            raw = yaml.safe_load(path.read_text())
        if not isinstance(raw, dict):
            raise ValueError("Configuration must be a mapping")
    raw.update({k: v for k, v in (overrides or {}).items() if v is not None})
    # YAML 1.1 can parse unadorned scientific notation as a string. Accept it
    # only in declared floating-point fields; the validator still rejects NaN.
    for item in fields(Config):
        if isinstance(item.default, float) and isinstance(raw.get(item.name), str):
            raw[item.name] = float(raw[item.name])
    unknown = set(raw) - {f.name for f in fields(Config)}
    if unknown:
        raise ValueError(f"Unknown configuration keys: {sorted(unknown)}")
    for key in ("covariates", "strengths", "l1_ratios", "stability_normalizations"):
        if key in raw:
            if not isinstance(raw[key], (list, tuple)):
                raise ValueError(f"{key} must be a sequence")
            raw[key] = tuple(raw[key])
    config = Config(**raw)
    config.validate()
    # Reject non-JSON values before a costly run starts.
    json.dumps(config.to_dict(), allow_nan=False)
    return config
