"""Command-line interface. No network or cloud dependency in analysis commands."""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any

from sigtrellis import __version__
from sigtrellis.cell_features import distribution_features
from sigtrellis.config import load_config
from sigtrellis.domain import file_hash
from sigtrellis.io import load_bulk
from sigtrellis.prediction import external_validate
from sigtrellis.simulate import write_example
from sigtrellis.singlecell import pseudobulk
from sigtrellis.workflow import run_analysis


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(
        description="SigTrellis: candidate transcriptomic signatures with auditable validation"
    )
    root.add_argument("--version", action="version", version=__version__)
    commands = root.add_subparsers(dest="command", required=True)
    for name in ("bulk", "single-cell"):
        sub = commands.add_parser(name)
        sub.add_argument("--config", type=Path)
        sub.add_argument("--output", type=Path, required=True)
        for flag in ("outcome", "outcome-type", "sample-id", "group", "batch", "positive-class"):
            sub.add_argument("--" + flag)
        sub.add_argument("--covariate", action="append", dest="covariates")
        sub.add_argument("--seed", type=int)
        sub.add_argument("--permutations", type=int)
        sub.add_argument("--stability-resamples", type=int)
        sub.add_argument("--outer-folds", type=int)
        sub.add_argument("--inner-folds", type=int)
        sub.add_argument("--supporting-de", action="store_true", default=None)
        sub.add_argument("--panel-validation", action="store_true", default=None)
        if name == "bulk":
            sub.add_argument("--expression", type=Path, required=True)
            sub.add_argument("--metadata", type=Path, required=True)
            sub.add_argument(
                "--orientation",
                choices=["auto", "samples_by_genes", "genes_by_samples"],
                default="auto",
            )
            sub.add_argument(
                "--annotations",
                type=Path,
                help="CSV with gene_id,gene_symbol; joins only after modeling",
            )
        else:
            sub.add_argument("--input", type=Path, required=True)
            sub.add_argument("--cell-type")
            sub.add_argument("--cell-type-value")
            sub.add_argument("--layer")
            sub.add_argument("--single-cell-mode", choices=["pseudobulk", "distribution"])
            sub.add_argument("--cell-state", action="append", dest="cell_states")
            sub.add_argument("--feature-block", action="append", dest="feature_blocks")
            sub.add_argument("--feature-gene", action="append", dest="feature_genes")
            sub.add_argument("--cell-resamples", type=int)
    simulate = commands.add_parser("simulate")
    simulate.add_argument("--output", type=Path, required=True)
    simulate.add_argument("--modality", choices=["bulk", "single-cell"], default="bulk")
    simulate.add_argument(
        "--samples", type=int, default=80, help="Biological donors for single-cell simulation"
    )
    simulate.add_argument("--features", type=int, default=200)
    simulate.add_argument("--seed", type=int, default=7)
    simulate.add_argument(
        "--scenario",
        choices=[
            "signal",
            "noise",
            "correlated",
            "batch_confounded",
            "batch_specific",
            "imbalance",
            "outliers",
            "multiclass",
            "continuous",
        ],
        default="signal",
    )
    external = commands.add_parser("external")
    for flag in ("run", "output"):
        external.add_argument("--" + flag, type=Path, required=True)
    for flag in ("expression", "metadata", "input"):
        external.add_argument("--" + flag, type=Path)
    external.add_argument("--model-kind", choices=["full", "panel"], default="full")
    external.add_argument(
        "--orientation", default="auto", choices=["auto", "samples_by_genes", "genes_by_samples"]
    )
    imported = commands.add_parser(
        "import-mtx", help="Convert raw Matrix Market/Seurat exports to H5AD"
    )
    imported.add_argument("--numeric-column", action="append", default=[], dest="numeric_columns")
    for flag in ("matrix", "genes", "barcodes", "metadata", "output"):
        imported.add_argument("--" + flag, type=Path, required=True)
    return root


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        if args.command == "import-mtx":
            from sigtrellis.matrix_import import import_matrix_market

            import_matrix_market(
                args.matrix,
                args.genes,
                args.barcodes,
                args.metadata,
                args.output,
                tuple(args.numeric_columns),
            )
            print(f"Aligned raw-count H5AD written to {args.output}")
            return 0
        if args.command == "simulate":
            write_example(
                args.output, args.modality, args.samples, args.features, args.seed, args.scenario
            )
            print(f"Synthetic fixture written to {args.output}")
            return 0
        if args.command == "external":
            external_validate(
                args.run,
                args.expression,
                args.metadata,
                args.output,
                args.orientation,
                single_cell=args.input,
                model_kind=args.model_kind,
            )
            print(f"Frozen external validation written to {args.output}")
            return 0
        excluded = {
            "command",
            "config",
            "output",
            "expression",
            "metadata",
            "orientation",
            "input",
            "annotations",
        }
        overrides: dict[str, Any] = {k: v for k, v in vars(args).items() if k not in excluded}
        config = load_config(args.config, overrides)
        if args.command == "bulk":
            datasets = [load_bulk(args.expression, args.metadata, config, args.orientation)]
            if args.annotations:
                import pandas as pd

                annotations = pd.read_csv(args.annotations, dtype=str)
                if (
                    not {"gene_id", "gene_symbol"} <= set(annotations)
                    or annotations.gene_id.duplicated().any()
                ):
                    raise ValueError("Annotations require unique gene_id and gene_symbol columns")
                datasets[0].symbols = {
                    str(k): str(v) for k, v in annotations.set_index("gene_id").gene_symbol.items()
                }
                datasets[0].input_hashes[str(args.annotations)] = file_hash(args.annotations)
        else:
            if config.single_cell_mode == "distribution":
                config = replace(
                    config, input_scale="features", normalization="none", imputation="median"
                )
                datasets = [distribution_features(args.input, config)]
            else:
                datasets = pseudobulk(args.input, config)
        if args.output.exists() and any(args.output.iterdir()):
            raise ValueError("Output directory must be empty")
        for index, data in enumerate(datasets):
            local_config = config
            if len(datasets) > 1:
                local_config = replace(
                    config, permutation_alpha=config.permutation_alpha / len(datasets)
                )
                data.upstream_qc["cell_type_family_tests"] = len(datasets)
                data.upstream_qc["family_permutation_alpha"] = config.permutation_alpha
            output = args.output
            if len(datasets) > 1:
                safe = re.sub(r"[^A-Za-z0-9_.-]", "_", data.cell_type or "all")[:60]
                output = output / f"{index:03d}_{safe}"
            print(
                f"Analyzing {len(data.expression)} samples; cell type: {data.cell_type or 'bulk'}",
                flush=True,
            )
            manifest = run_analysis(data, local_config, output)
            n = len({r["gene_id"] for r in manifest["candidate_features"]})
            print(
                f"COMPLETE: {n} features pass candidate-association gates; report: {output / 'report.html'}",
                flush=True,
            )
        return 0
    except (ValueError, OSError, ImportError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
