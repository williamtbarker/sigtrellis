"""Sample-local single-cell distributions with explicit biological feature identity.

No phenotype is used to define features. Cell states and gene programs are
declared in advance; cohort-trained embeddings are deliberately not implied.
"""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset
from scipy import sparse

from sigtrellis.cell_identity import CellFingerprint, cell_uniforms, gene_universe_hash
from sigtrellis.cell_moments import CellMoments, correlations
from sigtrellis.config import Config
from sigtrellis.domain import Dataset, Feature, FloatArray, file_hash


def feature_id(kind: str, state: str, source: str = "") -> str:
    """JSON encoding avoids collisions when source identifiers contain separators."""
    return json.dumps([kind, state, source], ensure_ascii=False, separators=(",", ":"))


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
            elif kind in {"gene_correlation", "program_correlation"}:
                is_gene = kind == "gene_correlation"
                pairs = moments.gene_pairs if is_gene else moments.program_pairs
                names = genes if is_gene else program_names
                values = correlations(
                    (moments.sums if is_gene else moments.program_sums)[positions],
                    (moments.squares if is_gene else moments.program_squares)[positions],
                    (moments.cross if is_gene else moments.program_cross)[positions],
                    numbers,
                    pairs,
                )
                for pi, (a, b) in enumerate(pairs):
                    first, second = names[a], names[b]
                    key = feature_id(kind, state, json.dumps([first, second]))
                    columns[key] = np.where(available, values[:, pi], np.nan)
                    identities[key] = Feature(
                        key,
                        kind,
                        state,
                        gene_id=first if is_gene else None,
                        gene_partner=second if is_gene else None,
                        program=first if not is_gene else None,
                        program_partner=second if not is_gene else None,
                        unit="within_cell_pearson_correlation",
                    )
            elif kind == "gene_tail":
                for ti, threshold in enumerate(moments.thresholds):
                    for gi, gene in enumerate(genes):
                        if gene not in rendered_genes:
                            continue
                        key = feature_id(kind, state, json.dumps([gene, float(threshold)]))
                        columns[key] = np.where(
                            available, moments.tails[positions, ti, gi] / safe_n, np.nan
                        )
                        identities[key] = Feature(
                            key,
                            kind,
                            state,
                            gene_id=gene,
                            threshold=float(threshold),
                            unit="fraction_cells_above_log1p_cp10k_threshold",
                        )
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
                    columns[key] = np.where(
                        available & ((numbers >= 2) if kind == "gene_variance" else True),
                        values[:, gi],
                        np.nan,
                    )
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
                        n = numbers[si]
                        total = moments.program_sums[pair, pi]
                        if kind == "program_mean":
                            values[si] = total / n
                        elif kind == "program_variance":
                            values[si] = (
                                max(moments.program_squares[pair, pi] - total**2 / n, 0) / (n - 1)
                                if n > 1
                                else np.nan
                            )
                        elif kind == "program_q90":
                            scores = np.concatenate(moments.programs[pair], axis=0)[:, pi]
                            values[si] = np.quantile(scores, 0.9)
                        else:
                            values[si] = moments.program_active[pair, pi] / n
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
        gene_blocks = set(config.feature_blocks) & {
            "gene_mean",
            "gene_detection",
            "gene_variance",
            "gene_tail",
        }
        needs_genes = config.supporting_de or bool(gene_blocks)
        pair_members = (
            {g for pair in config.gene_pairs for g in pair}
            if "gene_correlation" in config.feature_blocks
            else set()
        )
        # Count-DE normalization must retain the full supplied RNA universe,
        # even when the predictive gene-feature dictionary is deliberately small.
        genes = (
            list(all_genes if config.supporting_de else config.feature_genes or all_genes)
            if needs_genes
            else []
        )
        genes = list(dict.fromkeys([*genes, *sorted(pair_members)]))
        if (set(config.feature_genes) | pair_members) - set(all_genes):
            raise ValueError("A declared feature gene is absent from the count matrix")
        gene_lookup = {g: i for i, g in enumerate(all_genes)}
        gene_indices = [gene_lookup[g] for g in genes]
        program_names = sorted(config.programs)
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
        gene_positions = {g: i for i, g in enumerate(genes)}
        gene_pairs = (
            tuple(
                (gene_positions[sorted(pair)[0]], gene_positions[sorted(pair)[1]])
                for pair in config.gene_pairs
            )
            if "gene_correlation" in config.feature_blocks
            else ()
        )
        program_pairs = (
            tuple(combinations(range(len(program_names)), 2))
            if "program_correlation" in config.feature_blocks
            else ()
        )
        tail_thresholds = config.gene_thresholds if "gene_tail" in config.feature_blocks else ()
        retain_scores = "program_q90" in config.feature_blocks
        n_gene_blocks = len(gene_blocks - {"gene_tail"}) + len(tail_thresholds)
        n_program_blocks = sum(
            b.startswith("program_") and b != "program_correlation" for b in config.feature_blocks
        )
        output_columns = len(states) * (
            n_gene_blocks * len(config.feature_genes or genes)
            + n_program_blocks * len(program_names)
            + int("abundance" in config.feature_blocks)
            + len(gene_pairs)
            + len(program_pairs)
        )
        # Budget accumulator and output copies as well as retained program scores.
        # Sparse input chunks and Python object overhead remain additional costs.
        required_bytes = (1 + config.cell_resamples) * (
            n_pairs
            * (
                len(genes) * (4 + len(tail_thresholds))
                + len(program_names) * 3
                + len(gene_pairs)
                + len(program_pairs)
            )
            * 8
            + (len(obs) * len(program_names) * 16 if retain_scores else 0)
            + len(samples) * output_columns * 8 * 3
        )
        if required_bytes > config.max_dense_mb * 1024**2:
            raise ValueError(
                "Cell feature accumulators exceed max_dense_mb; use feature_genes/programs or fewer resamples"
            )
        moments = [
            CellMoments(
                n_pairs,
                len(genes),
                len(program_names),
                tail_thresholds,
                gene_pairs,
                program_pairs,
                retain_scores,
                np.array([config.program_thresholds[name] for name in program_names])
                if "program_fraction" in config.feature_blocks
                else None,
            )
            for _ in range(1 + config.cell_resamples)
        ]
        totals = np.zeros((len(moments), len(samples)))
        technical_totals = np.zeros((len(samples), 3))
        failed_per_sample = np.zeros(len(samples), dtype=np.int64)
        universe_hash = gene_universe_hash(all_genes)
        canonical_gene_indices = np.argsort(np.argsort(np.array(all_genes))).astype(np.int64)
        sample_hashes = {s: CellFingerprint(universe_hash) for s in samples}
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
            failed_per_sample += np.bincount(sample_dest[~keep], minlength=len(samples))
            for qi, measurement in enumerate(
                (libraries, detected, mito_counts / np.maximum(libraries, 1))
            ):
                technical_totals[:, qi] += np.bincount(
                    sample_dest[keep], weights=measurement[keep], minlength=len(samples)
                )
            for ci in np.flatnonzero(keep):
                lo, hi = chunk.indptr[ci : ci + 2]
                digest = sample_hashes[local_samples[ci]]
                digest.add(
                    canonical_gene_indices[chunk.indices[lo:hi]],
                    chunk.data[lo:hi],
                    str(local_labels[ci]),
                )
            selections = [
                keep,
                *[
                    keep
                    & (
                        cell_uniforms(obs.index[start:end].astype(str).tolist(), config.seed, i)
                        < config.cell_fraction
                    )
                    for i in range(config.cell_resamples)
                ],
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
        if expression.shape[1] < 1:
            raise ValueError("Declare at least one identifiable single-cell feature")
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
            if genes and config.supporting_de
        }
        qc: dict[str, Any] = {
            "modality": "single_cell_distributions",
            "input_cells": len(obs),
            "qc_failed_cells": removed,
            "cell_states": states,
            "unmodeled_states": sorted(set(labels) - set(states)),
            "sample_cell_hashes": {s: h.hexdigest() for s, h in sample_hashes.items()},
            "cells_per_sample": dict(zip(samples, totals[0].astype(int).tolist(), strict=True)),
            "sample_cell_quality": {
                sample: {
                    "mean_counts_per_cell": technical_totals[i, 0] / totals[0, i],
                    "mean_detected_genes": technical_totals[i, 1] / totals[0, i],
                    "mean_mito_fraction": technical_totals[i, 2] / totals[0, i]
                    if config.mitochondrial_prefix
                    else None,
                    "failed_cells": int(failed_per_sample[i]),
                }
                for i, sample in enumerate(samples)
            },
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
            "cell_feature_contract": {
                "version": 2,
                "rna_gene_universe_sha256": universe_hash,
                "n_rna_genes": len(all_genes),
                "normalization": "per_cell_log1p_cp10k_full_rna_universe",
                "fingerprint": "sha256_multiset_v2",
            },
            "retained_program_scores": retain_scores,
            "estimated_dense_bytes": required_bytes,
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
