import json
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from sigtrellis.cli import main
from sigtrellis.de import differential_expression
from sigtrellis.prediction import external_validate
from sigtrellis.qc import encode_outcome
from sigtrellis.simulate import synthetic_bulk, synthetic_single_cell
from sigtrellis.singlecell import pseudobulk
from sigtrellis.workflow import run_analysis


@pytest.mark.integration
def test_complete_report_and_external_frozen_validation(data, config, tmp_path):
    output = tmp_path / "results"
    result = run_analysis(data, config, output)
    assert result["status"] == "complete"
    expected = [
        "run_manifest.json",
        "qc_report.json",
        "qc_report.html",
        "model_metrics.json",
        "biomarkers.csv",
        "biomarker_stability.csv",
        "coefficients.csv",
        "cv_results.csv",
        "report.md",
        "report.html",
        "model_state.json",
        "audit.json",
        "resample_coefficients.npz",
    ]
    assert all((output / name).is_file() for name in expected)
    assert len(list((output / "figures").glob("*.png"))) >= 8
    html = (output / "report.html").read_text()
    assert "data:image/png;base64," in html
    assert "post-hoc consensus panel" in html
    assert not pd.read_csv(output / "biomarkers.csv").passes_robustness_gates.any()
    with pytest.raises(ValueError, match="not empty"):
        run_analysis(data, config, output)
    external, _ = synthetic_bulk(24, 60, seed=90)
    ids = [f"external_{i}" for i in range(24)]
    external.expression.index = ids
    external.metadata.index = ids
    external.metadata["sample_id"] = ids
    external.metadata["donor"] = ids
    external.expression.to_csv(tmp_path / "external.tsv", sep="\t", index_label="sample_id")
    external.metadata.to_csv(tmp_path / "external_meta.tsv", sep="\t", index=False)
    eval_result = external_validate(
        output,
        tmp_path / "external.tsv",
        tmp_path / "external_meta.tsv",
        tmp_path / "external_results",
    )
    assert not eval_result["training_transform_refit"]
    assert eval_result["metrics"]["roc_auc"] > 0.8
    # Renaming both sample and donor identifiers must not disguise copied training data.
    external.expression.iloc[0] = data.expression.iloc[0]
    external.expression.to_csv(tmp_path / "external.tsv", sep="\t", index_label="sample_id")
    with pytest.raises(ValueError, match="exact training expression profile"):
        external_validate(
            output,
            tmp_path / "external.tsv",
            tmp_path / "external_meta.tsv",
            tmp_path / "copied_results",
        )


@pytest.mark.integration
def test_singlecell_cli_complete(tmp_path):
    path = tmp_path / "cells.h5ad"
    synthetic_single_cell(path, donors=8, genes=25, cells=5)
    config = tmp_path / "single.yaml"
    config.write_text(
        "outcome: phenotype\ngroup: donor\ncell_type: cell_type\ncell_type_value: type_A\nmin_cells: 3\nouter_folds: 2\ninner_folds: 2\nstrengths: [0.05]\nl1_ratios: [0.5]\nstability_resamples: 2\npermutations: 0\n"
    )
    output = tmp_path / "out"
    assert (
        main(
            ["single-cell", "--input", str(path), "--config", str(config), "--output", str(output)]
        )
        == 0
    )
    assert (output / "pseudobulk_counts.tsv.gz").is_file()
    manifest = json.loads((output / "run_manifest.json").read_text())
    assert manifest["n_biological_groups"] == 8
    assert manifest["n_samples"] == 16


@pytest.mark.integration
def test_pydeseq2_actual_fit_and_design_errors(data, config):
    pytest.importorskip("pydeseq2")
    y, _ = encode_outcome(data, config)
    result = differential_expression(data, y, replace(config, supporting_de=True))
    planted = result.loc[["gene_0000", "gene_0001"]]
    assert (planted.de_log2_fold_change > 0).all()
    assert (planted.de_adjusted_pvalue < 0.05).all()
    assert result.attrs["design_columns"] == ["intercept", "condition"]
    data.metadata["batch"] = y
    with pytest.raises(ValueError, match="rank deficient"):
        differential_expression(data, y, replace(config, batch="batch"))


@pytest.mark.integration
def test_paired_pseudobulk_de(tmp_path, config):
    pytest.importorskip("pydeseq2")
    path = tmp_path / "cells.h5ad"
    synthetic_single_cell(path, donors=6, genes=20, cells=6)
    local = replace(
        config, cell_type="cell_type", cell_type_value="type_A", min_cells=4, de_pair_group=True
    )
    data = pseudobulk(path, local)[0]
    y, _ = encode_outcome(data, local)
    with pytest.raises(ValueError, match="repeated biological groups"):
        differential_expression(data, y, replace(local, de_pair_group=False))
    result = differential_expression(data, y, local)
    assert len(result.attrs["design_columns"]) == 7
    assert (result.loc[["gene_0000", "gene_0001"], "de_log2_fold_change"] > 0).all()


def test_failed_run_preserves_manifest(data, config, tmp_path):
    data.expression.iloc[0, 0] = np.nan
    with pytest.raises(ValueError):
        run_analysis(data, config, tmp_path / "bad")
    manifest = json.loads((tmp_path / "bad" / "run_manifest.json").read_text())
    assert manifest["status"] == "failed"
    assert "finite" in manifest["error"]


def test_cli_error_and_simulation(tmp_path):
    output = tmp_path / "fixture"
    assert main(["simulate", "--output", str(output), "--samples", "24", "--features", "20"]) == 0
    assert main(["simulate", "--output", str(output)]) == 1
    assert (output / "truth.json").exists()


@pytest.mark.integration
@pytest.mark.parametrize("kind", ["continuous", "multiclass"])
def test_nonbinary_report_paths(config, tmp_path, kind):
    data, _ = synthetic_bulk(48, 30, seed=34, scenario=kind)
    local = replace(config, outcome_type=kind, strengths=(0.05,), stability_resamples=2)
    run_analysis(data, local, tmp_path / kind)
    metrics = json.loads((tmp_path / kind / "model_metrics.json").read_text())
    assert metrics["loss_improvement"] > 0
    assert (tmp_path / kind / "report.html").stat().st_size > 10000
