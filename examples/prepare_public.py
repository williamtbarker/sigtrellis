"""Checksum-pinned public data preparation; no datasets are bundled with the code.

Run from the repository root after installing .[de,demo]. See DATA_SOURCES.md
for study provenance, source licenses, and the limits of both demonstrations.
"""

from __future__ import annotations

import argparse
import urllib.request
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from sigtrellis.domain import file_hash, write_json

BULK_BASE = "https://raw.githubusercontent.com/bioconnector/workshops/68ff3e4868b0e5673a1c4399828a3e0adc72f8e4/data/"
SOURCES = {
    "airway_counts.csv": (
        BULK_BASE + "airway_rawcounts.csv",
        "504f3148e23e84061f3f28d109a37732ccc9f1f8e4833b00b105307f2f8eb326",
    ),
    "airway_metadata.csv": (
        BULK_BASE + "airway_metadata.csv",
        "05bd7e78a0ca5b2a2f60ec715296381ca19f31a8c5a51e70c9f60583fa9fcf97",
    ),
    "kang_counts_25k.rds": (
        "https://zenodo.org/api/records/10069528/files/kang_counts_25k.RDS/content",
        "eb8a88259d84f757fcdb9a011085f62c73d36a59487249066c550001b29ea6a6",
    ),
}


def download(name: str, directory: Path) -> Path:
    url, expected = SOURCES[name]
    path = directory / name
    directory.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        temporary = path.with_suffix(path.suffix + ".part")
        try:
            with (
                urllib.request.urlopen(url, timeout=120) as response,
                temporary.open("wb") as target,
            ):
                while chunk := response.read(1024 * 1024):
                    target.write(chunk)
            if file_hash(temporary) != expected:
                raise ValueError(f"Checksum mismatch for {name}; refusing changed source")
            temporary.replace(path)
        finally:
            temporary.unlink(missing_ok=True)
    if file_hash(path) != expected:
        raise ValueError(f"Cached file checksum mismatch: {path}")
    return path


def bulk(output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    counts_path = download("airway_counts.csv", output / "downloads")
    metadata_path = download("airway_metadata.csv", output / "downloads")
    counts = pd.read_csv(counts_path, index_col=0)
    metadata = pd.read_csv(metadata_path).rename(columns={"id": "sample_id", "celltype": "donor"})
    if counts.shape != (64102, 8) or metadata.donor.nunique() != 4:
        raise ValueError("Unexpected pinned airway schema")
    if not np.equal(counts.to_numpy() % 1, 0).all():
        raise ValueError("Airway source is not raw integer counts")
    counts.to_csv(output / "counts.tsv", sep="\t", index_label="gene_id")
    metadata.to_csv(output / "metadata.tsv", sep="\t", index=False)
    write_json(
        output / "preparation.json",
        {
            "study": "Himes et al. 2014; GSE52778",
            "samples": 8,
            "biological_groups": 4,
            "input_sha256": {p.name: file_hash(p) for p in (counts_path, metadata_path)},
            "source_urls": {
                name: SOURCES[name][0] for name in ("airway_counts.csv", "airway_metadata.csv")
            },
            "license": "Source workshop distribution CC BY-NC-SA 4.0; separate from MIT software",
            "transformations": ["CSV to TSV", "Rename id to sample_id and celltype to donor"],
            "warning": "Only four biological donors: workflow demonstration, not biomarker validation",
        },
    )


def single_cell(output: Path) -> None:
    import rdata

    output.mkdir(parents=True, exist_ok=True)
    source = download("kang_counts_25k.rds", output / "downloads")
    # rdata preserves unknown S4 objects as namespaces. Only this pinned, verified
    # schema is supported here; this is not a general-purpose arbitrary RDS importer.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Missing constructor for R class.*")
        sce = rdata.read_rds(source)
    matrix = sce.assays.data.listData["counts"]
    if list(matrix.Dim) != [15706, 24673]:
        raise ValueError("Unexpected Kang matrix dimensions")
    x = sparse.csc_matrix((matrix.x, matrix.i, matrix.p), shape=tuple(matrix.Dim)).T.tocsr()
    if not np.isfinite(x.data).all() or (x.data < 0).any() or (x.data % 1 != 0).any():
        raise ValueError("Kang count source is not raw integer counts")
    obs = pd.DataFrame(sce.colData.listData, index=pd.Index(sce.colData.rownames.astype(str)))
    obs = obs[["sample", "patient", "condition", "cell_type"]].copy()
    # The original experiment's two condition-specific pooled libraries are a
    # design fact, not independently observed per-donor batch measurements.
    obs["condition_library"] = (
        obs["condition"].astype(str).map({"ctrl": "pooled_control", "stim": "pooled_stimulated"})
    )
    gene_names = np.asarray(sce.rowRanges.partitioning.NAMES, dtype=str)
    var = pd.DataFrame({"gene_symbol": gene_names}, index=gene_names)
    if obs.patient.nunique() != 8 or obs["sample"].nunique() != 16:
        raise ValueError("Unexpected donor/sample metadata")
    prepared = ad.AnnData(x, obs=obs, var=var)
    prepared.uns["preparation_source"] = "https://doi.org/10.5281/zenodo.10069528"
    prepared.write_h5ad(output / "kang.h5ad", compression="gzip")
    check = ad.read_h5ad(output / "kang.h5ad", backed="r")
    try:
        if check.shape != prepared.shape:
            raise ValueError("Prepared AnnData shape mismatch")
        # Cell counts and values in sampled deterministic chunks must survive serialization.
        for start in (0, 10000, 24000):
            block = check.X[start : start + 100, :]
            if (block != x[start : start + 100, :]).nnz:
                raise ValueError("Count corruption in H5AD conversion")
    finally:
        check.file.close()
    write_json(
        output / "preparation.json",
        {
            "study": "Kang et al. 2018; GSE96583",
            "deposit_author": "Daniel Dimitrov",
            "deposit_doi": "10.5281/zenodo.10069528",
            "license": "CC BY 4.0",
            "input_sha256": file_hash(source),
            "output_sha256": file_hash(output / "kang.h5ad"),
            "cells": len(obs),
            "genes": len(var),
            "donors": 8,
            "samples": 16,
            "raw_count_sum": int(x.sum()),
            "raw_count_nonzeros": int(x.nnz),
            "transformations": [
                "Read pinned dgCMatrix raw counts",
                "Transpose genes×cells to cells×genes",
                "Preserve donor, condition, sample and supplied cell-type annotations",
                "Record condition-specific pooled-library design",
                "Write sparse H5AD",
            ],
            "warning": "Condition and pooled-library effects are inseparable; annotations were supplied by the deposit",
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=["bulk", "single-cell", "all"], default="all")
    parser.add_argument("--output", type=Path, default=Path("data/public"))
    args = parser.parse_args()
    if args.dataset in {"bulk", "all"}:
        bulk(args.output / "bulk")
    if args.dataset in {"single-cell", "all"}:
        single_cell(args.output / "single_cell")
    print(f"Prepared public data in {args.output}")


if __name__ == "__main__":
    main()
