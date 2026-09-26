"""Chunked raw-count pseudobulk, never independent-cell phenotype modeling."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset
from scipy import sparse

from sigtrellis.config import Config
from sigtrellis.domain import Dataset, file_hash


def pseudobulk(path: Path, config: Config) -> list[Dataset]:
    """Sum QC-passing cells per experimental sample and cell type.

    Read only obs/var and count chunks. Layers are accessed directly in HDF5,
    avoiding read_h5ad(backed=...) materializing an entire selected layer.
    """
    if config.input_scale != "counts":
        raise ValueError("Pseudobulk requires raw counts, never log/integrated expression")
    config.validate()
    results: list[Dataset] = []
    with h5py.File(path, "r") as handle:
        obs: pd.DataFrame = read_elem(handle["obs"])
        var: pd.DataFrame = read_elem(handle["var"])
        if obs.index.has_duplicates or var.index.has_duplicates:
            raise ValueError("Duplicate cell or gene identifiers")
        required = [config.sample_id, config.outcome, *config.covariates]
        required += [
            v
            for v in (config.group, config.batch, config.cell_type, config.permutation_strata)
            if v
        ]
        if set(required) - set(obs.columns):
            raise ValueError(f"Missing cell metadata: {sorted(set(required) - set(obs.columns))}")
        if obs[required].isna().any().any():
            raise ValueError("Missing sample, phenotype, grouping, or cell-type metadata")
        obs[config.sample_id] = obs[config.sample_id].astype(str)
        # A sample ID denotes ONE biological specimen/condition, not a pooled donor label.
        invariant = list(
            dict.fromkeys(v for v in required if v not in {config.sample_id, config.cell_type})
        )
        if (obs.groupby(config.sample_id, observed=True)[invariant].nunique() > 1).any().any():
            raise ValueError(
                "Within-sample metadata vary; use unique sample IDs and a separate donor group"
            )
        if config.cell_type:
            types = obs[config.cell_type].astype(str)
        else:
            types = pd.Series("all_cells", index=obs.index)
        type_names = sorted(types.unique())
        if config.cell_type_value:
            if config.cell_type_value not in type_names:
                raise ValueError("Requested cell type is absent")
            type_names = [config.cell_type_value]
        samples = sorted(obs[config.sample_id].unique())
        pairs = [(t, s) for t in type_names for s in samples]
        lookup = {pair: i for i, pair in enumerate(pairs)}
        required_mb = len(pairs) * len(var) * 8 / 1024**2
        if required_mb > config.max_dense_mb:
            raise ValueError("Pseudobulk matrix exceeds max_dense_mb; select one cell type")
        sums = np.zeros((len(pairs), len(var)), dtype=np.float64)
        n_cells = np.zeros(len(pairs), dtype=np.int64)
        removed = 0
        total_selected = 0
        key = f"layers/{config.layer}" if config.layer else "X"
        if key not in handle:
            raise ValueError(f"Count source {key!r} is absent")
        element = handle[key]
        matrix: Any = sparse_dataset(element) if isinstance(element, h5py.Group) else element
        mito = np.array(
            [
                str(g).startswith(config.mitochondrial_prefix)
                if config.mitochondrial_prefix
                else False
                for g in var.index
            ]
        )
        for start in range(0, len(obs), config.chunk_size):
            end = min(start + config.chunk_size, len(obs))
            chunk = matrix[start:end, :]
            vals = chunk.data if sparse.issparse(chunk) else np.asarray(chunk)
            if (
                not np.isfinite(vals).all()
                or (vals < 0).any()
                or not np.allclose(vals, np.rint(vals), atol=1e-6, rtol=0)
            ):
                raise ValueError("Single-cell source is not finite nonnegative integer counts")
            lib = np.asarray(chunk.sum(axis=1)).ravel()
            detected = np.asarray((chunk > 0).sum(axis=1)).ravel()
            mito_counts = np.asarray(chunk[:, mito].sum(axis=1)).ravel()
            fractions = np.divide(
                mito_counts, lib, out=np.zeros_like(lib, dtype=float), where=lib > 0
            )
            keep = (
                (lib >= config.cell_min_counts)
                & (detected >= config.cell_min_genes)
                & (fractions <= config.max_mito_fraction)
            )
            local_types = types.iloc[start:end].to_numpy()
            local_samples = obs[config.sample_id].iloc[start:end].to_numpy()
            selected = np.isin(local_types, type_names)
            total_selected += int(selected.sum())
            removed += int((selected & ~keep).sum())
            positions = np.flatnonzero(keep & selected)
            dest = np.array([lookup[(local_types[i], local_samples[i])] for i in positions])
            if len(positions):
                incidence = sparse.csr_matrix(
                    (np.ones(len(positions)), (dest, positions)), shape=(len(pairs), end - start)
                )
                aggregate = incidence @ chunk
                sums += aggregate.toarray() if sparse.issparse(aggregate) else aggregate
                n_cells += np.bincount(dest, minlength=len(pairs))
        meta = obs.drop_duplicates(config.sample_id).set_index(config.sample_id, drop=False)
        for t in type_names:
            rows = [lookup[(t, s)] for s in samples if n_cells[lookup[(t, s)]] >= config.min_cells]
            retained_samples = [pairs[i][1] for i in rows]
            if len(rows) < 4:
                raise ValueError(f"Cell type {t!r} has fewer than four eligible samples")
            x = pd.DataFrame(sums[rows], index=retained_samples, columns=var.index.astype(str))
            m = meta.loc[retained_samples, [config.sample_id, *invariant]].copy()
            m["n_cells"] = n_cells[rows]
            dropped = [s for s in samples if s not in retained_samples]
            upstream: dict[str, Any] = {
                "modality": "single_cell_pseudobulk",
                "count_source": key,
                "input_cells": len(obs),
                "cells_in_requested_types": total_selected,
                "qc_failed_cells": removed,
                "retained_cells_this_type": int(n_cells[rows].sum()),
                "cells_per_sample": dict(
                    zip(retained_samples, n_cells[rows].tolist(), strict=True)
                ),
                "samples_below_min_cells": dropped,
                "note": "Cell counts do not increase the number of biological replicates",
            }
            symbols = (
                {str(k): str(v) for k, v in var["gene_symbol"].items()}
                if "gene_symbol" in var
                else {}
            )
            results.append(Dataset(x, m, t, {str(path): file_hash(path)}, upstream, symbols))
    return results
