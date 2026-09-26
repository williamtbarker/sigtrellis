"""Expression-only adapters; metadata never implicitly becomes a gene predictor."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pandas as pd

from sigtrellis.config import Config
from sigtrellis.domain import Dataset, file_hash


def read_table(path: Path, text_columns: tuple[str, ...] = ()) -> pd.DataFrame:
    sep = "\t" if ".tsv" in path.name else ","
    # pandas otherwise silently mangles duplicate headers into distinct feature names.
    import gzip

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        header = next(csv.reader(handle, delimiter=sep), [])
    if len(header) < 2:
        raise ValueError(f"Empty or invalid table: {path.name}")
    if len(header) != len(set(header)):
        raise ValueError(f"Duplicate headers in {path.name}")
    # Identifiers are strings, including '001' and literal 'NA'. Blank values
    # are rejected by the scientific input contract, not guessed as categories.
    frame = pd.read_csv(
        path,
        sep=sep,
        keep_default_na=False,
        converters={0: str},
        dtype={name: str for name in text_columns if name != header[0]},
    )
    # Construct the index afterward: the CSV index parser can infer numeric
    # identifiers even when the column converter returns strings.
    return frame.set_index(frame.columns[0])


def load_bulk(
    expression: Path, metadata: Path, config: Config, orientation: str = "auto"
) -> Dataset:
    x = read_table(expression)
    identities = tuple(
        name for name in (config.group, config.batch, config.permutation_strata) if name is not None
    )
    if config.outcome_type != "continuous":
        identities += (config.outcome,)
    m = read_table(metadata, identities)
    # First metadata column is the identifier, explicitly checked by name.
    if m.index.name != config.sample_id:
        raise ValueError(f"First metadata column must be {config.sample_id!r}")
    x.index = x.index.astype(str)
    x.columns = x.columns.astype(str)
    m.index = m.index.astype(str)
    if x.index.has_duplicates or x.columns.has_duplicates or m.index.has_duplicates:
        raise ValueError("Duplicate sample or gene identifiers")
    rows_match = set(x.index) == set(m.index)
    cols_match = set(x.columns) == set(m.index)
    if orientation == "auto":
        if rows_match == cols_match:
            raise ValueError("Ambiguous/mismatched sample axes; specify orientation and exact IDs")
        orientation = "samples_by_genes" if rows_match else "genes_by_samples"
    if orientation not in {"samples_by_genes", "genes_by_samples"}:
        raise ValueError("Invalid orientation")
    if orientation == "genes_by_samples":
        x = x.T
    if set(x.index) != set(m.index):
        raise ValueError("Expression and metadata sample IDs must match exactly")
    x = x.loc[m.index]
    if config.imputation == "median":
        x = x.replace({"": np.nan, "NaN": np.nan, "nan": np.nan})
    x = x.astype(np.float64)
    m[config.sample_id] = m.index
    return Dataset(
        x,
        m,
        input_hashes={str(expression): file_hash(expression), str(metadata): file_hash(metadata)},
    )
