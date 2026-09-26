"""Explicitly aligned Matrix Market/Seurat exports to native sparse AnnData."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from threadpoolctl import threadpool_limits

from sigtrellis.domain import file_hash
from sigtrellis.io import read_table


def import_matrix_market(
    matrix: Path,
    genes: Path,
    barcodes: Path,
    metadata: Path,
    output: Path,
    numeric_columns: tuple[str, ...] = (),
) -> None:
    """Matrix rows=genes, columns=cells; both identifier files are headerless TSV."""
    if output.exists():
        raise ValueError("Refusing to overwrite an existing H5AD")
    features = pd.read_csv(genes, sep="\t", header=None, dtype=str, keep_default_na=False)
    cells = pd.read_csv(barcodes, sep="\t", header=None, dtype=str, keep_default_na=False)[0]
    gene_ids = features[0]
    if cells.duplicated().any() or gene_ids.duplicated().any():
        raise ValueError("Duplicate cell/gene identifiers in Matrix Market export")
    if cells.str.strip().eq("").any() or gene_ids.str.strip().eq("").any():
        raise ValueError("Blank cell/gene identifier")
    if features.shape[1] > 2 and not features[2].eq("Gene Expression").all():
        raise ValueError("Supply RNA Gene Expression features only; antibody counts are not RNA")
    header = pd.read_csv(metadata, sep="\t" if ".tsv" in metadata.name else ",", nrows=0)
    obs = read_table(metadata, tuple(header.columns))
    for column in numeric_columns:
        if column not in obs:
            raise ValueError(f"Numeric metadata column {column!r} is absent")
        obs[column] = pd.to_numeric(obs[column], errors="raise")
    if obs.index.has_duplicates or set(obs.index) != set(cells):
        raise ValueError("Cell metadata must match the barcode set exactly")
    with threadpool_limits(limits=1):
        counts = sparse.csr_matrix(mmread(matrix).T)
    if counts.shape != (len(cells), len(gene_ids)):
        raise ValueError("Matrix dimensions disagree with gene/barcode files")
    if (
        not np.isfinite(counts.data).all()
        or (counts.data < 0).any()
        or not np.allclose(counts.data, np.rint(counts.data), atol=1e-6, rtol=0)
    ):
        raise ValueError("Matrix Market input must contain finite nonnegative integer counts")
    obs = obs.loc[cells.tolist()].copy()
    var = pd.DataFrame(index=pd.Index(gene_ids.to_numpy()))
    if features.shape[1] > 1:
        var["gene_symbol"] = features[1].to_numpy()
    result = ad.AnnData(counts, obs=obs, var=var)
    result.uns["sigtrellis_import"] = {
        # HDF5 interprets slashes in dictionary keys as paths. JSON keeps exact
        # source paths as values without changing the AnnData hierarchy.
        "input_hashes_json": json.dumps(
            {str(p): file_hash(p) for p in (matrix, genes, barcodes, metadata)}, sort_keys=True
        ),
        "alignment": "Metadata reordered to explicit matrix barcodes; no normalization",
        "numeric_columns": list(numeric_columns),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    result.write_h5ad(output, compression="gzip")
