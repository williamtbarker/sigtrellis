"""Sample-local single-cell distributions with explicit biological feature identity.

No phenotype is used to define features. Cell states and gene programs are
declared in advance; cohort-trained embeddings are deliberately not implied.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset
from scipy import sparse

from sigtrellis.config import Config
from sigtrellis.domain import Dataset, Feature, FloatArray, file_hash


def feature_id(kind: str, state: str, source: str = "") -> str:
    """JSON encoding avoids collisions when source identifiers contain separators."""
    return json.dumps([kind, state, source], ensure_ascii=False, separators=(",", ":"))


@dataclass
class CellMoments:
    n_pairs: int
    n_genes: int
    n_programs: int
    counts: FloatArray = field(init=False)
    sums: FloatArray = field(init=False)
    squares: FloatArray = field(init=False)
    detected: FloatArray = field(init=False)
    n: FloatArray = field(init=False)
    programs: list[list[FloatArray]] = field(init=False)

    def __post_init__(self) -> None:
        self.counts = np.zeros((self.n_pairs, self.n_genes))
        self.sums = np.zeros_like(self.counts)
        self.squares = np.zeros_like(self.counts)
        self.detected = np.zeros_like(self.counts)
        self.n = np.zeros(self.n_pairs)
        self.programs = [[] for _ in range(self.n_pairs)]

    def add(self, counts: Any, log_values: Any, scores: FloatArray, dest: Any) -> None:
        for pair in np.unique(dest):
            rows = np.flatnonzero(dest == pair)
            c, v = counts[rows], log_values[rows]
            self.n[pair] += len(rows)
            self.counts[pair] += np.asarray(c.sum(axis=0)).ravel()
            self.sums[pair] += np.asarray(v.sum(axis=0)).ravel()
            self.squares[pair] += np.asarray(v.multiply(v).sum(axis=0)).ravel()
            self.detected[pair] += np.asarray((c > 0).sum(axis=0)).ravel()
            if self.n_programs:
                self.programs[pair].append(scores[rows])


def _render(
    moments: CellMoments,
    totals: FloatArray,
    samples: list[str],
    states: tuple[str, ...],
    genes: list[str],
    program_names: list[str],
    config: Config,
) -> tuple[pd.DataFrame, dict[str, Feature]]:
    columns: dict[str, FloatArray] = {}
    identities: dict[str, Feature] = {}
    n_samples = len(samples)
    rendered_genes = set(config.feature_genes or genes)
    for ti, state in enumerate(states):
        positions = np.arange(ti * n_samples, (ti + 1) * n_samples)
        numbers = moments.n[positions]
        available = numbers >= config.min_cells
        safe_n = np.maximum(numbers, 1)
        means = moments.sums[positions] / safe_n[:, None]
        variances = (
            np.maximum(
                moments.squares[positions] - moments.sums[positions] ** 2 / safe_n[:, None], 0
            )
            / np.maximum(numbers - 1, 1)[:, None]
        )
        for kind in config.feature_blocks:
            if kind == "abundance":
                key = feature_id(kind, state)
                # Regularize proportions, not counts: otherwise an absent/pure
                # population becomes a predictor of the number of captured cells.
                fraction = numbers / np.maximum(totals, 1)
                values = np.log((fraction + 0.001) / (1 - fraction + 0.001))
                values[totals == 0] = np.nan
                columns[key] = values
                identities[key] = Feature(key, kind, state, unit="log_relative_odds")
            elif kind.startswith("gene_"):
                values = {
                    "gene_mean": means,
                    "gene_variance": variances,
                    "gene_detection": moments.detected[positions] / safe_n[:, None],
                }[kind]
                for gi, gene in enumerate(genes):
                    if gene not in rendered_genes:
                        continue
                    key = feature_id(kind, state, gene)
                    columns[key] = np.where(available, values[:, gi], np.nan)
                    unit = {
                        "gene_mean": "mean_log1p_cp10k",
                        "gene_variance": "variance_log1p_cp10k",
                        "gene_detection": "detected_cell_fraction",
                    }[kind]
                    identities[key] = Feature(key, kind, state, gene_id=gene, unit=unit)
            else:
                for pi, program in enumerate(program_names):
                    values = np.full(n_samples, np.nan)
                    for si, pair in enumerate(positions):
                        if not available[si]:
                            continue
                        scores = np.concatenate(moments.programs[pair], axis=0)[:, pi]
                        if kind == "program_mean":
                            values[si] = scores.mean()
                        elif kind == "program_variance":
                            values[si] = scores.var(ddof=1) if len(scores) > 1 else np.nan
                        elif kind == "program_q90":
                            values[si] = np.quantile(scores, 0.9)
                        else:
                            values[si] = (scores > config.program_thresholds[program]).mean()
                    key = feature_id(kind, state, program)
                    columns[key] = values
                    identities[key] = Feature(
                        key,
                        kind,
                        state,
                        program=program,
                        unit="fraction"
                        if kind == "program_fraction"
                        else "squared_log1p_cp10k_program_score"
                        if kind == "program_variance"
                        else "log1p_cp10k_program_score",
                    )
    return pd.DataFrame(columns, index=pd.Index(samples)), identities


def distribution_features(path: Path, config: Config) -> Dataset:
    """Retain state abundance and within-sample moments without cell pseudoreplication.

    Raw counts are read in sparse chunks. Only small program-score arrays are
    retained for quantiles. Cell subsampling is outcome-blind and does not create
    additional sample rows. Each normalization uses that cell's full library.
    """
    config.validate()
    if config.input_scale != "features" or config.normalization != "none":
        raise ValueError("Distribution mode requires input_scale=features and normalization=none")
    if not config.cell_states and config.cell_type:
        raise ValueError(
            "Declare cell_states in advance for a fixed, externally transferable representation"
        )
    states = config.cell_states or ("all_cells",)
    with h5py.File(path, "r") as handle:
        obs: pd.DataFrame = read_elem(handle["obs"])
        var: pd.DataFrame = read_elem(handle["var"])
        if obs.index.has_duplicates or var.index.has_duplicates:
            raise ValueError("Duplicate cell or gene identifiers")
        required = list(
            dict.fromkeys(
                [
                    config.sample_id,
                    config.outcome,
                    *config.covariates,
                    *[
                        v
                        for v in (
                            config.group,
                            config.batch,
                            config.cell_type,
                            config.permutation_strata,
                            config.time,
                        )
                        if v
                    ],
                ]
            )
        )
        if set(required) - set(obs.columns):
            raise ValueError("Missing required single-cell metadata")
        if obs[required].isna().any().any() or any(
            obs[k].astype(str).str.strip().eq("").any() for k in required
        ):
            raise ValueError("Missing or blank required single-cell metadata")
        obs[config.sample_id] = obs[config.sample_id].astype(str)
        invariant = [v for v in required if v not in {config.sample_id, config.cell_type}]
        if (obs.groupby(config.sample_id, observed=True)[invariant].nunique() > 1).any().any():
            raise ValueError(
                "Within-sample metadata vary; use specimen IDs and a separate biological group"
            )
        samples = sorted(obs[config.sample_id].unique())
        sample_indices = {s: i for i, s in enumerate(samples)}
        all_genes = var.index.astype(str).tolist()
        needs_genes = config.supporting_de or any(
            block.startswith("gene_") for block in config.feature_blocks
        )
        # Count-DE normalization must retain the full supplied RNA universe,
        # even when the predictive gene-feature dictionary is deliberately small.
        genes = (
            list(all_genes if config.supporting_de else config.feature_genes or all_genes)
            if needs_genes
            else []
        )
        if set(config.feature_genes) - set(all_genes):
            raise ValueError("A declared feature gene is absent from the count matrix")
        gene_lookup = {g: i for i, g in enumerate(all_genes)}
        gene_indices = [gene_lookup[g] for g in genes]
        program_names = list(config.programs)
        for name, members in config.programs.items():
            if set(members) - set(all_genes):
                raise ValueError(
                    f"Program {name!r} contains absent genes; map identifiers explicitly"
                )
        program_weights = sparse.lil_matrix((len(all_genes), len(program_names)))
        for pi, name in enumerate(program_names):
            for gene in config.programs[name]:
                program_weights[gene_lookup[gene], pi] = 1 / len(config.programs[name])
        program_weights = program_weights.tocsr()
        n_pairs = len(samples) * len(states)
        n_gene_blocks = sum(b.startswith("gene_") for b in config.feature_blocks)
        n_program_blocks = sum(b.startswith("program_") for b in config.feature_blocks)
        output_columns = len(states) * (
            n_gene_blocks * len(config.feature_genes or genes)
            + n_program_blocks * len(program_names)
            + int("abundance" in config.feature_blocks)
        )
        # Budget accumulator and output copies as well as retained program scores.
        # Sparse input chunks and Python object overhead remain additional costs.
        required_bytes = (1 + config.cell_resamples) * (
            n_pairs * len(genes) * 8 * 4
            + len(obs) * len(program_names) * 8
            + len(samples) * output_columns * 8 * 3
        )
        if required_bytes > config.max_dense_mb * 1024**2:
            raise ValueError(
                "Cell feature accumulators exceed max_dense_mb; use feature_genes/programs or fewer resamples"
            )
        moments = [
            CellMoments(n_pairs, len(genes), len(program_names))
            for _ in range(1 + config.cell_resamples)
        ]
        totals = np.zeros((len(moments), len(samples)))
        generators = [
            np.random.default_rng(config.seed + 65537 + i) for i in range(config.cell_resamples)
        ]
        sample_hashes = {s: hashlib.sha256() for s in samples}
        state_indices = {s: i for i, s in enumerate(states)}
        labels = (
            obs[config.cell_type].astype(str).to_numpy()
            if config.cell_type
            else np.repeat("all_cells", len(obs))
        )
        specimen = obs[config.sample_id].to_numpy()
        key = f"layers/{config.layer}" if config.layer else "X"
        if key not in handle:
            raise ValueError("Declared raw count source is absent")
        element = handle[key]
        matrix: Any = sparse_dataset(element) if isinstance(element, h5py.Group) else element
        if tuple(matrix.shape) != (len(obs), len(var)):
            raise ValueError("Count matrix shape does not match cell/gene metadata")
        mito = np.array(
            [
                g.startswith(config.mitochondrial_prefix) if config.mitochondrial_prefix else False
                for g in all_genes
            ]
        )
        removed = 0
        for start in range(0, len(obs), config.chunk_size):
            end = min(start + config.chunk_size, len(obs))
            chunk = sparse.csr_matrix(matrix[start:end, :], dtype=float)
            chunk.sum_duplicates()
            chunk.eliminate_zeros()
            chunk.sort_indices()
            if (
                not np.isfinite(chunk.data).all()
                or (chunk.data < 0).any()
                or not np.allclose(chunk.data, np.rint(chunk.data), atol=1e-6, rtol=0)
            ):
                raise ValueError("Single-cell features require finite nonnegative integer counts")
            libraries = np.asarray(chunk.sum(axis=1)).ravel()
            detected = np.diff(chunk.indptr)
            mito_counts = np.asarray(chunk[:, mito].sum(axis=1)).ravel()
            keep = (libraries >= config.cell_min_counts) & (detected >= config.cell_min_genes)
            keep &= mito_counts / np.maximum(libraries, 1) <= config.max_mito_fraction
            removed += int((~keep).sum())
            normalized = sparse.diags(1e4 / np.maximum(libraries, 1)) @ chunk
            normalized.data = np.log1p(normalized.data)
            scores = np.asarray((normalized @ program_weights).toarray(), dtype=float)
            local_samples, local_labels = specimen[start:end], labels[start:end]
            sample_dest = np.array([sample_indices[s] for s in local_samples])
            for ci in np.flatnonzero(keep):
                lo, hi = chunk.indptr[ci : ci + 2]
                digest = sample_hashes[local_samples[ci]]
                digest.update(np.asarray([hi - lo], dtype="<i8").tobytes())
                digest.update(chunk.indices[lo:hi].astype("<i8").tobytes())
                digest.update(chunk.data[lo:hi].astype("<f8").tobytes())
                digest.update(str(local_labels[ci]).encode() + b"\0")
            selections = [
                keep,
                *[keep & (rng.random(end - start) < config.cell_fraction) for rng in generators],
            ]
            for ri, selected in enumerate(selections):
                totals[ri] += np.bincount(sample_dest[selected], minlength=len(samples))
                positions = np.flatnonzero(selected & np.isin(local_labels, states))
                dest = np.array(
                    [
                        state_indices[local_labels[i]] * len(samples) + sample_dest[i]
                        for i in positions
                    ],
                    dtype=int,
                )
                if len(positions):
                    moments[ri].add(
                        chunk[positions][:, gene_indices],
                        normalized[positions][:, gene_indices],
                        scores[positions],
                        dest,
                    )
        if (totals[0] < config.min_cells).any():
            raise ValueError(
                "A biological sample has too few QC-passing cells; review sample QC before analysis"
            )
        frames = [
            _render(m, total, samples, states, genes, program_names, config)
            for m, total in zip(moments, totals, strict=True)
        ]
        expression, identities = frames[0]
        if expression.shape[1] < 2:
            raise ValueError("Declare at least two identifiable single-cell features")
        meta = obs.drop_duplicates(config.sample_id).set_index(config.sample_id, drop=False)
        meta = meta.loc[samples, [config.sample_id, *invariant]].copy()
        meta["n_cells"] = totals[0].astype(int)
        count_tables = {
            state: pd.DataFrame(
                moments[0].counts[ti * len(samples) : (ti + 1) * len(samples)],
                index=samples,
                columns=genes,
            )
            for ti, state in enumerate(states)
            if genes
        }
        qc: dict[str, Any] = {
            "modality": "single_cell_distributions",
            "input_cells": len(obs),
            "qc_failed_cells": removed,
            "cell_states": states,
            "unmodeled_states": sorted(set(labels) - set(states)),
            "sample_cell_hashes": {s: h.hexdigest() for s, h in sample_hashes.items()},
            "cells_per_sample": dict(zip(samples, totals[0].astype(int).tolist(), strict=True)),
            "cells_per_state": {
                state: dict(
                    zip(
                        samples,
                        moments[0]
                        .n[ti * len(samples) : (ti + 1) * len(samples)]
                        .astype(int)
                        .tolist(),
                        strict=True,
                    )
                )
                for ti, state in enumerate(states)
            },
            "cell_resamples": config.cell_resamples,
            "cell_fraction": config.cell_fraction,
            "feature_definition": "Fixed state annotations/programs; sample-local log1p(CP10K) moments; no phenotype-trained embedding",
            "variance_interpretation": "Observed cell variability includes measurement noise; not a deconvolved biological variance or causal effect",
            "abundance_interpretation": "Relative recovery among QC-passing cells; not absolute tissue cell counts",
            "abundance_definition": "log((state_fraction + 0.001)/(1 - state_fraction + 0.001)); fixed proportion offset, independent of cell count",
        }
        symbols = (
            {str(k): str(v) for k, v in var["gene_symbol"].items()} if "gene_symbol" in var else {}
        )
        return Dataset(
            expression,
            meta,
            "multiple_states",
            {str(path): file_hash(path)},
            qc,
            symbols,
            identities,
            tuple(frame for frame, _ in frames[1:]),
            count_tables,
        )
