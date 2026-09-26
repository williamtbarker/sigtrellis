"""Small shared domain types and reproducibility helpers."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd

FloatArray = npt.NDArray[np.float64]
IntArray = npt.NDArray[np.int64]


@dataclass(frozen=True)
class Feature:
    """Scientific identity of a predictor, separate from its matrix column name."""

    feature_id: str
    kind: str = "gene_expression"
    cell_type: str | None = None
    gene_id: str | None = None
    program: str | None = None
    unit: str = "expression"


@dataclass
class Dataset:
    expression: pd.DataFrame
    metadata: pd.DataFrame
    cell_type: str | None = None
    input_hashes: dict[str, str] = field(default_factory=dict)
    upstream_qc: dict[str, Any] = field(default_factory=dict)
    symbols: dict[str, str] = field(default_factory=dict)
    features: dict[str, Feature] = field(default_factory=dict)
    cell_resamples: tuple[pd.DataFrame, ...] = ()
    counts_by_cell_type: dict[str, pd.DataFrame] = field(default_factory=dict)

    def subset(self, rows: IntArray) -> Dataset:
        return Dataset(
            self.expression.iloc[rows],
            self.metadata.iloc[rows],
            self.cell_type,
            self.input_hashes,
            self.upstream_qc,
            self.symbols,
            self.features,
            tuple(frame.iloc[rows] for frame in self.cell_resamples),
            {name: frame.iloc[rows] for name, frame in self.counts_by_cell_type.items()},
        )


@dataclass(frozen=True)
class Split:
    train: IntArray
    test: IntArray
    repeat: int
    fold: int


@dataclass
class Audit:
    fits: list[dict[str, Any]] = field(default_factory=list)
    splits: list[dict[str, Any]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def file_hash(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def array_hash(values: FloatArray) -> str:
    # Numerically equal profiles must fingerprint equally: signed zero must not
    # let a duplicate specimen cross the biological validation boundary.
    canonical = np.array(values, dtype="<f8", order="C", copy=True)
    canonical[canonical == 0] = 0.0
    return hashlib.sha256(canonical.tobytes()).hexdigest()


def write_json(path: Path, value: Any) -> None:
    def clean(obj: Any) -> Any:
        if isinstance(obj, dict):
            return {str(k): clean(v) for k, v in obj.items()}
        if isinstance(obj, (tuple, list, np.ndarray)):
            return [clean(v) for v in obj]
        if isinstance(obj, (np.integer, np.bool_)):
            return obj.item()
        if isinstance(obj, (float, np.floating)):
            return float(obj) if np.isfinite(obj) else None
        return obj

    path.write_text(json.dumps(clean(value), indent=2, sort_keys=True, allow_nan=False) + "\n")
