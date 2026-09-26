"""Reproduce the controlled scientific benchmark table, independently of public data."""

from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path
from typing import Any

from threadpoolctl import threadpool_limits

from sigtrellis.config import Config
from sigtrellis.domain import Audit, write_json
from sigtrellis.qc import encode_outcome, validate_dataset
from sigtrellis.simulate import synthetic_bulk
from sigtrellis.stability import stability_select
from sigtrellis.validation import nested_validate, permutation_control
from sigtrellis.workflow import batch_diagnostics


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=Path("docs/validation/scientific_benchmarks.json")
    )
    args = parser.parse_args()
    config = Config(
        group="donor",
        outer_folds=3,
        inner_folds=2,
        strengths=(0.02, 0.15),
        l1_ratios=(0.5,),
        stability_resamples=12,
        permutations=19,
        max_features=150,
    )
    records: list[dict[str, Any]] = []
    with threadpool_limits(limits=1):
        for scenario in [
            "signal",
            "noise",
            "correlated",
            "batch_confounded",
            "batch_specific",
            "imbalance",
            "outliers",
            "multiclass",
            "continuous",
        ]:
            data, truth = synthetic_bulk(80, 200, 17, scenario)
            local = replace(
                config,
                outcome_type=scenario if scenario in ("multiclass", "continuous") else "binary",
                batch="batch" if scenario.startswith("batch_") else None,
            )
            y, _ = encode_outcome(data, local)
            qc, audit = validate_dataset(data, local), Audit()
            validated = nested_validate(data, y, local, audit)
            metrics = {
                k: v
                for k, v in validated.metrics["out_of_fold"].items()
                if not isinstance(v, (list, dict))
            }
            row: dict[str, Any] = {
                "scenario": scenario,
                "seed": 17,
                "n_samples": 80,
                "n_features": 200,
                "metrics": metrics,
                "loss_improvement": validated.metrics["loss_improvement"],
                "qc_warnings": qc["warnings"],
                "qc_gate_blockers": qc["gate_blockers"],
            }
            if scenario in ("signal", "noise"):
                row["permutation_pvalue"] = permutation_control(
                    data, y, local, validated.metrics["loss_improvement"], audit
                )["pvalue"]
                stability = stability_select(data, y, local, ["1_vs_0"], audit)
                stable = stability.table.query(
                    "selection_frequency>=.8 and sign_consistency>=.9"
                ).gene_id.tolist()
                row.update(
                    stable_genes=stable,
                    planted_genes=truth,
                    planted_recovered=len(set(stable) & set(truth)),
                    stable_nonplanted=len(set(stable) - set(truth)),
                )
            if scenario.startswith("batch_"):
                row["batch_validation"] = batch_diagnostics(data, y, local, audit)[0]
            records.append(row)
            print(scenario, metrics, flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_json(
        args.output,
        {
            "scope": "Controlled synthetic validation; no biological claims",
            "config": config.to_dict(),
            "records": records,
        },
    )


if __name__ == "__main__":
    main()
