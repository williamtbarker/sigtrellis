"""Strict, serializable run contracts. Unknown keys are errors."""

from __future__ import annotations

import json
import math
import tomllib
from dataclasses import asdict, dataclass, field, fields
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
    imputation: str = "reject"
    max_missing_fraction: float = 0.5
    single_cell_mode: str = "pseudobulk"
    cell_states: tuple[str, ...] = ()
    feature_genes: tuple[str, ...] = ()
    gene_thresholds: tuple[float, ...] = (1.0, 2.0, 3.0)
    gene_pairs: tuple[tuple[str, str], ...] = ()
    feature_blocks: tuple[str, ...] = (
        "abundance",
        "gene_mean",
        "gene_detection",
        "program_mean",
        "program_variance",
    )
    programs: dict[str, tuple[str, ...]] = field(default_factory=dict)
    program_thresholds: dict[str, float] = field(default_factory=dict)
    cell_resamples: int = 0
    cell_fraction: float = 0.75
    panel_validation: bool = False
    panel_max_features: int = 20
    time: str | None = None
    temporal_gap: float = 0.0
    temporal_train_fraction: float = 0.5

    def validate(self) -> None:
        for item in fields(self):
            if (
                isinstance(item.default, (int, float))
                and not isinstance(item.default, bool)
                and isinstance(getattr(self, item.name), bool)
            ):
                raise ValueError(f"{item.name} must be numeric, not a boolean")
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
            "time",
        ):
            value = getattr(self, key)
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise ValueError(f"{key} must be a nonempty string")
        for key in ("supporting_de", "de_pair_group", "assess_batches", "panel_validation"):
            if not isinstance(getattr(self, key), bool):
                raise ValueError(f"{key} must be a boolean")
        choices = {
            "outcome_type": {"binary", "multiclass", "continuous"},
            "input_scale": {"counts", "log_expression", "features"},
            "normalization": {"logcpm", "median_ratio", "none"},
            "candidate_method": {"variance", "association", "deseq2", "none"},
            "cv_strategy": {"grouped", "leave_group_out", "leave_batch_out", "temporal"},
            "permutation_scheme": {"group", "within_group"},
            "tuning_rule": {"one_se", "minimum_loss"},
            "single_cell_mode": {"pseudobulk", "distribution"},
            "imputation": {"reject", "median"},
        }
        for key, allowed in choices.items():
            if getattr(self, key) not in allowed:
                raise ValueError(f"{key} must be one of {sorted(allowed)}")
        protected = {self.outcome, self.sample_id, self.group, self.batch, self.time}
        if self.cell_type and self.cell_type in protected | {self.permutation_strata}:
            raise ValueError(
                "Cell-type annotation cannot be the outcome, an identifier, batch, time or permutation stratum"
            )
        if any(c in protected for c in self.covariates):
            raise ValueError("Outcome, sample/group identifiers, and batch cannot be predictors")
        if len(set(self.covariates)) != len(self.covariates):
            raise ValueError("Duplicate covariates")
        if any(not isinstance(c, str) or not c.strip() for c in self.covariates):
            raise ValueError("Covariates must be nonempty column names")
        if self.outcome in {self.sample_id, self.group}:
            raise ValueError("Outcome cannot be a sample or group identifier")
        if self.input_scale in {"log_expression", "features"} and self.normalization != "none":
            raise ValueError("Declared log_expression/features requires normalization=none")
        if self.input_scale == "counts" and self.normalization == "none":
            raise ValueError("Counts require logcpm or median_ratio normalization")
        if self.candidate_method == "deseq2" or self.supporting_de:
            if self.input_scale != "counts" and not (
                self.single_cell_mode == "distribution"
                and self.supporting_de
                and self.candidate_method != "deseq2"
            ):
                raise ValueError(
                    "PyDESeq2 requires raw counts; mixed features use association/variance screening"
                )
        if self.imputation == "median" and self.input_scale == "counts":
            raise ValueError("Missing raw counts cannot be imputed as expression")
        if self.de_pair_group and not self.group:
            raise ValueError("Paired DE requires an explicit group column")
        if self.cv_strategy == "leave_batch_out" and not self.batch:
            raise ValueError("leave_batch_out requires batch metadata")
        if self.cv_strategy == "temporal" and (not self.time or self.repeats != 1):
            raise ValueError("Temporal validation requires time metadata and repeats=1")
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
            "panel_max_features",
        ):
            if not isinstance(getattr(self, key), int) or getattr(self, key) < 1:
                raise ValueError(f"{key} must be a positive integer")
        for key in ("seed", "permutations", "cell_resamples"):
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
            "max_missing_fraction",
        ):
            if (
                not isinstance(getattr(self, key), (float, int))
                or not math.isfinite(getattr(self, key))
                or not 0 < getattr(self, key) <= 1
            ):
                raise ValueError(f"{key} must be in (0,1]")
        for key in (
            "stability_fraction",
            "decision_threshold",
            "cell_fraction",
            "temporal_train_fraction",
        ):
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
            isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) or x <= 0
            for x in self.strengths
        ):
            raise ValueError("strengths must be finite and positive")
        if not self.l1_ratios or any(
            isinstance(x, bool)
            or not isinstance(x, (int, float))
            or not math.isfinite(x)
            or not 0 < x <= 1
            for x in self.l1_ratios
        ):
            raise ValueError("l1_ratios must be in (0,1]; pure ridge is not sparse selection")
        if len(set(self.strengths)) != len(self.strengths) or len(set(self.l1_ratios)) != len(
            self.l1_ratios
        ):
            raise ValueError("Duplicate hyperparameters would duplicate inner-fold evidence")
        for norm in self.stability_normalizations:
            if norm not in choices["normalization"]:
                raise ValueError("Invalid stability normalization")
            if (norm == "none") != (self.input_scale != "counts"):
                raise ValueError("Stability normalization is incompatible with input_scale")
        if len(set(self.stability_normalizations)) != len(self.stability_normalizations):
            raise ValueError("Duplicate stability normalizations")
        if len(self.stability_normalizations) > self.stability_resamples:
            raise ValueError("Every stability normalization needs at least one resample")
        if (
            not isinstance(self.temporal_gap, (int, float))
            or not math.isfinite(self.temporal_gap)
            or self.temporal_gap < 0
        ):
            raise ValueError("temporal_gap must be finite and nonnegative")
        if self.cell_resamples + 1 > self.stability_resamples:
            raise ValueError("Every cell perturbation needs at least one stability resample")
        if self.cell_resamples and self.single_cell_mode != "distribution":
            raise ValueError("Cell perturbations require single_cell_mode=distribution")
        blocks = {
            "abundance",
            "gene_mean",
            "gene_detection",
            "gene_variance",
            "gene_tail",
            "gene_correlation",
            "program_mean",
            "program_variance",
            "program_q90",
            "program_fraction",
            "program_correlation",
        }
        if not self.feature_blocks or set(self.feature_blocks) - blocks:
            raise ValueError("Unknown or empty single-cell feature_blocks")
        for key in ("cell_states", "feature_genes", "feature_blocks"):
            values = getattr(self, key)
            if len(set(values)) != len(values) or any(
                not isinstance(v, str) or not v.strip() for v in values
            ):
                raise ValueError(f"{key} must contain unique nonempty strings")
        if not isinstance(self.programs, dict) or not isinstance(self.program_thresholds, dict):
            raise ValueError("Programs and program_thresholds must be mappings")
        if (
            not self.gene_thresholds
            or any(
                isinstance(t, bool)
                or not isinstance(t, (float, int))
                or not math.isfinite(t)
                or t <= 0
                for t in self.gene_thresholds
            )
            or len(set(self.gene_thresholds)) != len(self.gene_thresholds)
        ):
            raise ValueError("gene_thresholds must be unique finite positive log1p(CP10K) values")
        pairs: set[tuple[str, ...]] = set()
        for pair in self.gene_pairs:
            if (
                not isinstance(pair, (tuple, list))
                or len(pair) != 2
                or any(not isinstance(g, str) or not g.strip() for g in pair)
                or pair[0] == pair[1]
            ):
                raise ValueError("gene_pairs must contain two distinct nonempty gene IDs per pair")
            pair_key = tuple(sorted(pair))
            if pair_key in pairs:
                raise ValueError("Duplicate or reversed gene pair")
            pairs.add(pair_key)
        if "gene_correlation" in self.feature_blocks and not self.gene_pairs:
            raise ValueError("gene_correlation requires prespecified gene_pairs")
        if "program_correlation" in self.feature_blocks and len(self.programs) < 2:
            raise ValueError("program_correlation requires at least two programs")
        if self.cell_states and not self.cell_type and self.cell_states != ("all_cells",):
            raise ValueError("Named cell_states require a cell_type annotation column")
        if self.single_cell_mode == "distribution" and self.cell_type_value:
            raise ValueError("Distribution mode uses cell_states, not cell_type_value")
        for name, genes in self.programs.items():
            if isinstance(genes, str):
                raise ValueError("Program genes must be a sequence, not a string")
            if not isinstance(name, str) or not name.strip() or not genes:
                raise ValueError("Programs require nonempty names and gene lists")
            if len(set(genes)) != len(genes) or any(
                not isinstance(g, str) or not g.strip() for g in genes
            ):
                raise ValueError("Program genes must be unique nonempty strings")
        if set(self.program_thresholds) - set(self.programs):
            raise ValueError("Threshold names must match declared programs")
        if any(
            isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v)
            for v in self.program_thresholds.values()
        ):
            raise ValueError("Program thresholds must be finite numbers")
        if "program_fraction" in self.feature_blocks and set(self.program_thresholds) != set(
            self.programs
        ):
            raise ValueError(
                "Activation fractions require a predefined threshold for every program"
            )

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
    for key in (
        "covariates",
        "strengths",
        "l1_ratios",
        "stability_normalizations",
        "cell_states",
        "feature_genes",
        "feature_blocks",
        "gene_thresholds",
        "gene_pairs",
    ):
        if key in raw:
            if not isinstance(raw[key], (list, tuple)):
                raise ValueError(f"{key} must be a sequence")
            raw[key] = tuple(raw[key])
    if "gene_pairs" in raw:
        if any(not isinstance(p, (tuple, list)) for p in raw["gene_pairs"]):
            raise ValueError("gene_pairs must be a sequence of pairs")
        raw["gene_pairs"] = tuple(tuple(p) for p in raw["gene_pairs"])
    if "programs" in raw:
        if not isinstance(raw["programs"], dict) or any(
            not isinstance(v, (list, tuple)) for v in raw["programs"].values()
        ):
            raise ValueError("programs must map names to gene lists")
        raw["programs"] = {name: tuple(genes) for name, genes in raw["programs"].items()}
    config = Config(**raw)
    config.validate()
    # Reject non-JSON values before a costly run starts.
    json.dumps(config.to_dict(), allow_nan=False)
    return config
