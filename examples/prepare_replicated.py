"""Pinned replicated studies, with untouched biological-sample holdouts.

Run from the source repository. Downloads are explicit; no raw data enter Git.
"""

from __future__ import annotations

import argparse
import re
import tarfile
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from sigtrellis.domain import file_hash, write_json

YEAST_BASE = "https://raw.githubusercontent.com/bartongroup/profDGE48/375dc0d57d9d1fa96a4245a6530e0fda34305891/"
YEAST = {
    "Preprocessed_data/Snf2_countdata.tar.gz": "8cafcaa7052a760cbec9196605f1e588c58c76f772db064ace599c0bf07d727a",
    "Preprocessed_data/WT_countdata.tar.gz": "11b33cf3fb8bd9b7c229c81ec00f3ab502f61a1715d0f59821b9bc7a7eb625ac",
    "Bad_replicate_identification/exclude.lst": "57d5e5a7512388b7df2255095b0cb0938438658256262e29b60cd12056214a5e",
}
LUPUS_URL = "https://datasets.cellxgene.cziscience.com/c55dc602-d168-4d15-acc1-5de4f2f5d551.h5ad"
LUPUS_SHA256 = "3c0b74d54c03838a49817edce95314e6a4fa048ef35b74d1e697b8b2fd07cc03"


def download_checked(url: str, expected: str, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        temporary = path.with_suffix(path.suffix + ".part")
        with urllib.request.urlopen(url, timeout=180) as response, temporary.open("wb") as target:
            while block := response.read(4 * 1024**2):
                target.write(block)
        if file_hash(temporary) != expected:
            raise ValueError("Downloaded source checksum mismatch")
        temporary.replace(path)
    if file_hash(path) != expected:
        raise ValueError(f"Source checksum mismatch: {path}")
    return path


def yeast(output: Path) -> None:
    paths = {
        name: download_checked(YEAST_BASE + name, checksum, output / "downloads" / Path(name).name)
        for name, checksum in YEAST.items()
    }
    excluded = {
        name + ".gbgout"
        for name in paths["Bad_replicate_identification/exclude.lst"].read_text().splitlines()
    }
    values: dict[str, pd.Series] = {}
    metadata = []
    for name, path in paths.items():
        if not name.endswith(".tar.gz"):
            continue
        with tarfile.open(path) as archive:
            for member in sorted(archive.getmembers(), key=lambda m: m.name):
                if not member.isfile() or member.name in excluded:
                    continue
                handle = archive.extractfile(member)
                assert handle is not None
                frame = pd.read_csv(
                    handle, sep="\t", header=None, names=["gene", "count"], dtype={"gene": str}
                )
                frame = frame[
                    ~frame.gene.isin(
                        [
                            "no_feature",
                            "ambiguous",
                            "too_low_aQual",
                            "not_aligned",
                            "alignment_not_unique",
                        ]
                    )
                ]
                replicate = int(re.search(r"rep(\d+)_", member.name).group(1))
                condition = "wildtype" if member.name.startswith("WT") else "snf2_knockout"
                sample = f"culture_{len(metadata):03d}"
                values[sample] = frame.set_index("gene")["count"]
                metadata.append(
                    {
                        "sample_id": sample,
                        "culture": sample,
                        "condition": condition,
                        "replicate": replicate,
                        "source_file": member.name,
                    }
                )
    counts, meta = pd.DataFrame(values), pd.DataFrame(metadata).set_index("sample_id", drop=False)
    if counts.shape[1] != 86 or counts.isna().any().any():
        raise ValueError("Unexpected replicated yeast source schema")
    rng = np.random.default_rng(2026)
    holdout: list[str] = []
    for label in sorted(meta.condition.unique()):
        ids = meta.index[meta.condition == label].to_numpy()
        holdout.extend(rng.choice(ids, 10, replace=False).tolist())
    for split, ids in (
        ("train", [s for s in meta.index if s not in holdout]),
        ("external", sorted(holdout)),
    ):
        destination = output / split
        destination.mkdir(parents=True, exist_ok=True)
        counts.loc[:, ids].to_csv(destination / "counts.tsv.gz", sep="\t", index_label="gene_id")
        meta.loc[ids].to_csv(destination / "metadata.tsv", sep="\t", index=False)
    write_json(
        output / "preparation.json",
        {
            "study": "Schurch et al. 2016; doi:10.1261/rna.053959.115; ERP004763",
            "source_urls": {name: YEAST_BASE + name for name in YEAST},
            "sha256": {name: file_hash(path) for name, path in paths.items()},
            "license": "MIT source distribution; retain Barton group attribution",
            "n_cultures": 86,
            "n_training": 66,
            "n_untouched": 20,
            "genes": len(counts),
            "published_qc_exclusions": sorted(excluded),
            "holdout_policy": "10 independently grown cultures per genotype chosen by seed 2026 before fitting",
            "limitations": "Same-study held-out cultures, not a separate laboratory/cohort. Genotype response demonstration, not clinical biomarker qualification. Published upstream QC is inherited.",
        },
    )


def lupus_partitions(obs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    obs = obs.copy()
    obs["specimen_id"] = obs.sample_uuid + ":cohort:" + obs.Processing_Cohort
    if (obs.groupby("specimen_id")[["donor_id", "disease"]].nunique() != 1).any().any():
        raise ValueError("Source specimen metadata are inconsistent")
    meta = obs.drop_duplicates("specimen_id")
    external_donors = set(meta.loc[meta.Processing_Cohort == "4.0", "donor_id"])
    train = meta[meta.Processing_Cohort.isin(["2.0", "3.0"]) & ~meta.donor_id.isin(external_donors)]
    external = meta[meta.Processing_Cohort == "4.0"]
    train = train.sort_values(["Processing_Cohort", "specimen_id"]).drop_duplicates("donor_id")
    external = external.sort_values("specimen_id").drop_duplicates("donor_id")
    if set(train.donor_id) & set(external.donor_id):
        raise ValueError("External donors leaked into preparation")
    return obs, train, external


def write_count_subset(matrix, rows, obs, var, target: Path, policy: str) -> None:
    """Write selected raw cells in bounded sparse chunks using the AnnData schema."""
    import h5py
    from anndata.io import write_elem

    temporary = target.with_suffix(".partial.h5ad")
    with h5py.File(temporary, "w") as destination:
        destination.attrs.update({"encoding-type": "anndata", "encoding-version": "0.1.0"})
        write_elem(destination, "obs", obs)
        write_elem(destination, "var", var)
        write_elem(destination, "uns", {"source": LUPUS_URL, "subset_policy": policy})
        group = destination.create_group("X")
        group.attrs.update(
            {
                "encoding-type": "csr_matrix",
                "encoding-version": "0.1.0",
                "shape": (len(rows), len(var)),
            }
        )
        options = {"shape": (0,), "maxshape": (None,), "compression": "gzip", "chunks": True}
        stored_data = group.create_dataset("data", dtype="int64", **options)
        stored_indices = group.create_dataset("indices", dtype="int32", **options)
        stored_indptr = group.create_dataset("indptr", shape=(len(rows) + 1,), dtype="int64")
        nnz = 0
        stored_indptr[0] = 0
        for start in range(0, len(rows), 2048):
            end = min(start + 2048, len(rows))
            values = matrix[rows[start:end], :].tocsr()
            if (
                not np.isfinite(values.data).all()
                or (values.data < 0).any()
                or not np.allclose(values.data, np.rint(values.data), atol=1e-6, rtol=0)
            ):
                raise ValueError("The pinned raw matrix is not integer counts")
            stored_data.resize((nnz + values.nnz,))
            stored_indices.resize((nnz + values.nnz,))
            stored_data[nnz:] = values.data.astype(np.int64)
            stored_indices[nnz:] = values.indices
            stored_indptr[start + 1 : end + 1] = nnz + values.indptr[1:]
            nnz += values.nnz
    temporary.replace(target)


def lupus(source: Path, output: Path, cells_per_sample: int = 400) -> None:
    import h5py
    from anndata.io import read_elem, sparse_dataset

    if file_hash(source) != LUPUS_SHA256:
        raise ValueError("Lupus source does not match the pinned release checksum")
    output.mkdir(parents=True, exist_ok=True)
    with h5py.File(source, "r") as handle:
        columns = [
            "donor_id",
            "sample_uuid",
            "library_uuid",
            "Processing_Cohort",
            "disease",
            "sex",
            "author_cell_type",
        ]
        obs = pd.DataFrame({c: read_elem(handle["obs"][c]).astype(str) for c in columns})
        # A specimen was sometimes processed in multiple cohorts. Preserve the
        # aliquot boundary before reserving cohort 4 and purging its donors.
        obs, train, external = lupus_partitions(obs)
        var = read_elem(handle["raw/var"])
        var["gene_symbol"] = var["feature_name"].astype(str)
        matrix = sparse_dataset(handle["raw/X"])
        records = {}
        for label, selected in (("train", train), ("external", external)):
            rng = np.random.default_rng(2026)
            rows = []
            for sample in sorted(selected.specimen_id):
                candidates = np.flatnonzero(obs.specimen_id.to_numpy() == sample)
                rows.extend(
                    candidates.tolist()
                    if cells_per_sample == 0
                    else rng.choice(
                        candidates, min(cells_per_sample, len(candidates)), replace=False
                    ).tolist()
                )
            rows = np.sort(rows)
            cells = obs.iloc[rows].copy()
            cells.index = pd.Index([f"source_cell_{i}" for i in rows])
            policy = (
                "All source cells in the declared specimens; no cell downsampling"
                if cells_per_sample == 0
                else (
                    "Seed 2026 sample-local cell subsampling; all genes retained; no expression/phenotype-based cell selection"
                )
            )
            write_count_subset(matrix, rows, cells, var.copy(), output / f"{label}.h5ad", policy)
            selected.to_csv(output / f"{label}_samples.tsv", sep="\t", index=False)
            records[label] = {
                "samples": len(selected),
                "donors": selected.donor_id.nunique(),
                "cells": len(cells),
                "classes": selected.disease.value_counts().to_dict(),
                "sha256": file_hash(output / f"{label}.h5ad"),
            }
            print(label, records[label], flush=True)
    write_json(
        output / "preparation.json",
        {
            "study": "Perez et al. 2022; doi:10.1126/science.abf1970",
            "collection": "https://cellxgene.cziscience.com/collections/436154da-bcf1-4130-9c8b-120ff9a888f2",
            "source": LUPUS_URL,
            "source_sha256": LUPUS_SHA256,
            "license": "CC BY 4.0; data contributed by Perez et al., curated/distributed by CZ CELLxGENE Discover",
            "source_donors": 261,
            "source_cells": 1263676,
            "count_source": "raw/X",
            "seed": 2026,
            "holdout_policy": "Processing cohort 4 external; cohorts 2/3 training after purging every external donor; one specimen per donor",
            "cohort_1_exclusion": "Only healthy controls; omitted by predeclared case/control design",
            "cells_per_sample_cap": cells_per_sample,
            "records": records,
            "limitations": "Processing-cohort holdout from one study; observational disease/treatment associations. Supplied cell annotations may reflect upstream integration. Not an independent clinical validation study.",
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=["yeast", "lupus"], required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source", type=Path)
    parser.add_argument(
        "--cells-per-sample",
        type=int,
        default=400,
        help="0 retains every source cell in each declared specimen",
    )
    parser.add_argument(
        "--download", action="store_true", help="Explicitly allow the 12.2 GB lupus source download"
    )
    args = parser.parse_args()
    if args.cells_per_sample < 0:
        parser.error("cells-per-sample must be nonnegative")
    if args.dataset == "yeast":
        yeast(args.output)
    else:
        source = args.source
        if source is None:
            if not args.download:
                parser.error("Lupus requires --source or --download (12.2 GB)")
            source = download_checked(
                LUPUS_URL, LUPUS_SHA256, args.output / "downloads/source.h5ad"
            )
        lupus(source, args.output, args.cells_per_sample)


if __name__ == "__main__":
    main()
