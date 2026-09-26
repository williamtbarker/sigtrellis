"""Sparse streaming distribution statistics; raw cells are never phenotype replicates."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from sigtrellis.domain import FloatArray


@dataclass
class CellMoments:
    n_pairs: int
    n_genes: int
    n_programs: int
    thresholds: tuple[float, ...] = ()
    gene_pairs: tuple[tuple[int, int], ...] = ()
    program_pairs: tuple[tuple[int, int], ...] = ()
    retain_scores: bool = False
    activation_thresholds: FloatArray | None = None
    counts: FloatArray = field(init=False)
    sums: FloatArray = field(init=False)
    squares: FloatArray = field(init=False)
    detected: FloatArray = field(init=False)
    tails: FloatArray = field(init=False)
    cross: FloatArray = field(init=False)
    program_sums: FloatArray = field(init=False)
    program_squares: FloatArray = field(init=False)
    program_cross: FloatArray = field(init=False)
    program_active: FloatArray = field(init=False)
    n: FloatArray = field(init=False)
    programs: list[list[FloatArray]] = field(init=False)

    def __post_init__(self) -> None:
        self.counts = np.zeros((self.n_pairs, self.n_genes))
        self.sums = np.zeros_like(self.counts)
        self.squares = np.zeros_like(self.counts)
        self.detected = np.zeros_like(self.counts)
        self.tails = np.zeros((self.n_pairs, len(self.thresholds), self.n_genes))
        self.cross = np.zeros((self.n_pairs, len(self.gene_pairs)))
        self.program_sums = np.zeros((self.n_pairs, self.n_programs))
        self.program_squares = np.zeros_like(self.program_sums)
        self.program_active = np.zeros_like(self.program_sums)
        self.program_cross = np.zeros((self.n_pairs, len(self.program_pairs)))
        self.n = np.zeros(self.n_pairs)
        self.programs = [[] for _ in range(self.n_pairs)]

    def add(self, counts: Any, log_values: Any, scores: FloatArray, dest: Any) -> None:
        # Compute each product once per input chunk, then reduce by specimen/state.
        # Repeating sparse column slicing for every group is prohibitively slow
        # when pooled sequencing interleaves hundreds of donors in one chunk.
        for pi, (a, b) in enumerate(self.gene_pairs):
            products = np.asarray(log_values[:, a].multiply(log_values[:, b]).toarray()).ravel()
            self.cross[:, pi] += np.bincount(dest, weights=products, minlength=self.n_pairs)
        for pi, (a, b) in enumerate(self.program_pairs):
            self.program_cross[:, pi] += np.bincount(
                dest, weights=scores[:, a] * scores[:, b], minlength=self.n_pairs
            )
        for pair in np.unique(dest):
            rows = np.flatnonzero(dest == pair)
            c, v, p = counts[rows], log_values[rows], scores[rows]
            self.n[pair] += len(rows)
            self.counts[pair] += np.asarray(c.sum(axis=0)).ravel()
            self.sums[pair] += np.asarray(v.sum(axis=0)).ravel()
            self.squares[pair] += np.asarray(v.multiply(v).sum(axis=0)).ravel()
            self.detected[pair] += np.asarray((c > 0).sum(axis=0)).ravel()
            for ti, threshold in enumerate(self.thresholds):
                self.tails[pair, ti] += np.asarray((v > threshold).sum(axis=0)).ravel()
            self.program_sums[pair] += p.sum(axis=0)
            self.program_squares[pair] += (p * p).sum(axis=0)
            if self.activation_thresholds is not None:
                self.program_active[pair] += (p > self.activation_thresholds).sum(axis=0)
            if self.n_programs and self.retain_scores:
                self.programs[pair].append(p)


def correlations(
    sums: FloatArray,
    squares: FloatArray,
    cross: FloatArray,
    numbers: FloatArray,
    pairs: tuple[tuple[int, int], ...],
) -> FloatArray:
    result = np.full(cross.shape, np.nan)
    safe_n = np.maximum(numbers, 1)
    centered = np.maximum(squares - sums**2 / safe_n[:, None], 0)
    for pi, (a, b) in enumerate(pairs):
        denominator = np.sqrt(centered[:, a] * centered[:, b])
        valid = (numbers >= 3) & (centered[:, a] > 1e-10) & (centered[:, b] > 1e-10)
        numerator = cross[:, pi] - sums[:, a] * sums[:, b] / safe_n
        np.divide(numerator, denominator, out=result[:, pi], where=valid)
    return np.asarray(np.clip(result, -1.0, 1.0), dtype=float)
