"""Record exact metrics for the controlled coupling fixture used by pytest.

Run with the package installed: python tests/record_cell_benchmarks.py --output PATH
This is validation evidence, not a biological analysis or a production adapter.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

from test_cell_adversarial import joint_config, joint_fixture
from threadpoolctl import threadpool_limits

from sigtrellis.cell_features import distribution_features
from sigtrellis.domain import Audit, write_json
from sigtrellis.qc import encode_outcome
from sigtrellis.singlecell import pseudobulk
from sigtrellis.validation import nested_validate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    records = {}
    with TemporaryDirectory() as temporary, threadpool_limits(limits=1):
        path = joint_fixture(Path(temporary) / "coupling.h5ad")
        joint = joint_config()
        configs = {
            "within_cell_coupling": joint,
            "marginal_summaries": joint_config(
                feature_blocks=("gene_mean", "gene_variance", "gene_detection")
            ),
            "count_pseudobulk": replace(
                joint,
                single_cell_mode="pseudobulk",
                input_scale="counts",
                normalization="logcpm",
                imputation="reject",
                cell_type_value="A",
            ),
        }
        for name, config in configs.items():
            data = (
                pseudobulk(path, config)[0]
                if config.single_cell_mode == "pseudobulk"
                else distribution_features(path, config)
            )
            y, _ = encode_outcome(data, config)
            audit = Audit()
            result = nested_validate(data, y, config, audit)
            assert all(not set(s["train_groups"]) & set(s["test_groups"]) for s in audit.splits)
            records[name] = {"metrics": result.metrics, "configuration": config.to_dict()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_json(
        args.output,
        {
            "scope": "Controlled synthetic capability test, not general performance superiority",
            "fixture": "tests/test_cell_adversarial.py:joint_fixture",
            "seed": 314,
            "biological_samples": 48,
            "cells_per_sample": 80,
            "records": records,
        },
    )
    for name, record in records.items():
        print(name, record["metrics"]["out_of_fold"]["roc_auc"])


if __name__ == "__main__":
    main()
