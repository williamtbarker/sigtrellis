"""Order-invariant cell identity and reproducible cell-local perturbations."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass

import numpy as np

from sigtrellis.domain import FloatArray, IntArray


def gene_universe_hash(genes: list[str]) -> str:
    return hashlib.sha256(json.dumps(sorted(genes), ensure_ascii=False).encode()).hexdigest()


def cell_uniforms(ids: list[str], seed: int, perturbation: int) -> FloatArray:
    """A fixed cell ID's inclusion never depends on other cells or row order."""
    prefix = f"sigtrellis-cell:{seed}:{perturbation}:".encode()
    return np.array(
        [
            (
                int.from_bytes(
                    hashlib.blake2b(prefix + s.encode(), digest_size=8).digest(), "little"
                )
                >> 11
            )
            / 2**53
            for s in ids
        ],
        dtype=float,
    )


@dataclass
class CellFingerprint:
    """Cryptographic multiset digest, independent of cell/sample names and order.

    Modular addition preserves multiplicity (unlike XOR). This detects exact
    reordered cell collections, not partial overlap or biological relatives.
    """

    universe: str
    count: int = 0
    total: int = 0

    def add(self, indices: IntArray, values: FloatArray, state: str) -> None:
        order = np.argsort(indices)
        digest = hashlib.sha256()
        digest.update(self.universe.encode())
        digest.update(json.dumps(state, ensure_ascii=False).encode())
        digest.update(indices[order].astype("<i8").tobytes())
        digest.update(values[order].astype("<f8").tobytes())
        self.total = (self.total + int.from_bytes(digest.digest(), "little")) % 2**256
        self.count += 1

    def hexdigest(self) -> str:
        payload = self.count.to_bytes(8, "little") + self.total.to_bytes(32, "little")
        return hashlib.sha256(payload).hexdigest()
