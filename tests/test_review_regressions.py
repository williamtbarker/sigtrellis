"""Reproductions of defects found in the second adversarial source review."""

import json
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from sigtrellis.correlation import correlation_diagnostics
from sigtrellis.de import differential_expression
from sigtrellis.domain import Audit
from sigtrellis.io import load_bulk
from sigtrellis.metrics import evaluate
from sigtrellis.prediction import external_validate
from sigtrellis.qc import encode_outcome, validate_dataset
from sigtrellis.simulate import synthetic_bulk
from sigtrellis.stability import summarize_coefficients
from sigtrellis.validation import nested_validate, permutation_control
from sigtrellis.workflow import evidence_gates, run_analysis


@pytest.mark.parametrize(
    "change",
    [
        {"seed": True},
        {"max_features": True},
        {"min_count": True},
        {"strengths": (0.05, 0.05)},
        {"l1_ratios": (0.5, 0.5)},
        {"strengths": (True,)},
        {"l1_ratios": (True,)},
    ],
)
def test_invalid_numeric_configuration_is_not_silently_coerced(config, change):
    with pytest.raises(ValueError):
        replace(config, **change).validate()


def test_inconsistent_h5ad_dimensions_fail_before_aggregation(config, tmp_path):
    import h5py

    from sigtrellis.simulate import synthetic_single_cell
    from sigtrellis.singlecell import pseudobulk

    path = tmp_path / "inconsistent.h5ad"
    synthetic_single_cell(path, donors=6, genes=20, cells=5)
    with h5py.File(path, "r+") as handle:
        n, p = handle["X"].attrs["shape"]
        handle["X"].attrs["shape"] = (n + 1, p)
    with pytest.raises(ValueError, match="shape does not match"):
        pseudobulk(path, replace(config, min_cells=3))


@pytest.mark.parametrize("orientation", ["samples_by_genes", "genes_by_samples"])
def test_identifiers_survive_csv_parsing_verbatim(data, config, tmp_path, orientation):
    ids = [f"{i:04d}" for i in range(len(data.expression))]
    data.expression.index = ids
    data.metadata.index = ids
    data.metadata["sample_id"] = ids
    data.metadata["donor"] = ["01", "1", "NA", *ids[3:]]
    data.metadata["phenotype"] = np.where(data.metadata.phenotype == 1, "01", "1")
    x = data.expression if orientation == "samples_by_genes" else data.expression.T
    x.to_csv(tmp_path / "x.tsv", sep="\t")
    data.metadata.to_csv(tmp_path / "m.tsv", sep="\t", index=False)
    loaded = load_bulk(tmp_path / "x.tsv", tmp_path / "m.tsv", config, orientation)
    assert loaded.expression.index.tolist() == ids
    assert loaded.metadata.donor.tolist() == data.metadata.donor.tolist()
    assert loaded.metadata.phenotype.tolist() == data.metadata.phenotype.tolist()


def test_blank_biological_identity_is_rejected(data, config):
    data.metadata.loc[data.metadata.index[0], "donor"] = "  "
    with pytest.raises(ValueError, match="[Bb]lank"):
        validate_dataset(data, config)


def test_signed_zero_cannot_disguise_an_exact_duplicate(data, config):
    data.expression = data.expression.astype(float)
    data.expression.iloc[0, 0] = 0.0
    data.expression.iloc[1] = data.expression.iloc[0]
    data.expression.iloc[1, 0] = -0.0
    np.testing.assert_array_equal(data.expression.iloc[0], data.expression.iloc[1])
    with pytest.raises(ValueError, match="Identical expression"):
        validate_dataset(data, config)


def test_paired_de_rejects_repeated_donor_condition_before_backend(data, config, monkeypatch):
    dds = pytest.importorskip("pydeseq2.dds")
    data.metadata["donor"] = [f"d{i // 2}" for i in range(len(data.expression))]
    data.metadata["phenotype"] = np.tile([0, 1], len(data.expression) // 2)
    data.metadata.loc[data.metadata.index[1], "phenotype"] = 0
    local = replace(config, de_pair_group=True)

    def forbidden_backend(**kwargs):
        raise AssertionError("Count backend reached an unsupported repeated-condition design")

    monkeypatch.setattr(dds, "DeseqDataSet", forbidden_backend)
    with pytest.raises(ValueError, match="one observation per biological group and condition"):
        differential_expression(data, encode_outcome(data, local)[0], local)


def test_multiclass_correlated_group_evidence_keeps_contrasts_separate(data, config):
    coefs = np.zeros((6, 2, data.expression.shape[1]))
    coefs[:, 0, 0] = 1
    coefs[0, 1, 0] = 1
    ranking = summarize_coefficients(
        coefs, list(data.expression.columns), ["a_vs_c", "b_vs_c"], 1e-7
    )
    result = correlation_diagnostics(
        data,
        coefs,
        ranking,
        replace(config, outcome_type="multiclass"),
        contrasts=["a_vs_c", "b_vs_c"],
    )
    assert "contrast" in result.groups.columns
    freq = result.groups.set_index("contrast").selection_frequency
    assert freq["a_vs_c"] == 1
    assert freq["b_vs_c"] == pytest.approx(1 / 6)


@pytest.mark.parametrize(
    "bad_field,bad_value",
    [
        ("loss_improvement", np.nan),
        ("loss_improvement", np.inf),
        ("pvalue", np.nan),
        ("pvalue", -0.01),
        ("pvalue", 0.0),
    ],
)
def test_invalid_statistics_never_pass_evidence_gates(config, bad_field, bad_value):
    table = pd.DataFrame(
        {
            "selection_frequency": [1.0],
            "sign_consistency": [1.0],
            "outer_selection_frequency": [1.0],
            "outer_sign_consistency": [1.0],
            "frequency_logcpm": [1.0],
            "coefficient": [1.0],
            "coefficient_median_selected": [1.0],
            "outer_dominant_sign": [1.0],
        }
    )
    metrics = {"loss_improvement": 1.0, "permutation": {"status": "tested", "pvalue": 0.01}}
    if bad_field == "pvalue":
        metrics["permutation"][bad_field] = bad_value
    else:
        metrics[bad_field] = bad_value
    gated, blockers = evidence_gates(
        table,
        replace(config, permutations=99),
        {"gate_blockers": []},
        metrics,
        {"status": "not_assessed"},
    )
    assert blockers
    assert not gated.passes_robustness_gates.any()


def test_nonfinite_observed_permutation_statistic_fails(data, config, monkeypatch):
    import sigtrellis.validation as validation

    def forbidden_validation(*args, **kwargs):
        raise AssertionError("Nonfinite statistic reached expensive permutation fitting")

    monkeypatch.setattr(validation, "nested_validate", forbidden_validation)
    with pytest.raises(ValueError, match="finite"):
        permutation_control(
            data, encode_outcome(data, config)[0], replace(config, permutations=19), np.nan, Audit()
        )


def test_single_class_fold_does_not_claim_balanced_discrimination(config):
    y = np.zeros(4)
    metrics = evaluate(y, np.tile([0.8, 0.2], (4, 1)), config, np.ones(4))
    assert metrics["balanced_accuracy"] is None
    assert metrics["mcc"] is None
    assert metrics["brier_score"] == pytest.approx(0.04)


def test_constant_regression_fold_has_undefined_r_squared(config):
    metrics = evaluate(
        np.ones(4),
        np.arange(4, dtype=float),
        replace(config, outcome_type="continuous"),
        np.ones(4),
    )
    assert metrics["r2"] is None
    assert metrics["rmse"] > 0


@pytest.mark.integration
@pytest.mark.parametrize(
    "artifact", ["model_state.json", "configuration.json", "sample_metadata.csv"]
)
def test_changed_training_artifacts_cannot_be_externally_validated(
    data, config, tmp_path, artifact
):
    run = tmp_path / "training"
    run_analysis(data, replace(config, strengths=(0.05,), stability_resamples=2), run)
    external, _ = synthetic_bulk(24, data.expression.shape[1], seed=73)
    ids = [f"external_{i}" for i in range(len(external.expression))]
    external.expression.index = ids
    external.metadata.index = ids
    external.metadata["sample_id"] = ids
    external.metadata["donor"] = ids
    xp, mp = tmp_path / "x.tsv", tmp_path / "m.tsv"
    external.expression.to_csv(xp, sep="\t")
    external.metadata.to_csv(mp, sep="\t", index=False)
    target = run / artifact
    if artifact == "model_state.json":
        state = json.loads(target.read_text())
        state["intercept"][0] += 100
        target.write_text(json.dumps(state))
    elif artifact == "configuration.json":
        state = json.loads(target.read_text())
        state["strengths"] = [1.0]
        target.write_text(json.dumps(state))
    else:
        target.write_text(target.read_text() + "\n")
    with pytest.raises(ValueError, match="[Ii]ntegrity"):
        external_validate(run, xp, mp, tmp_path / "external")


@pytest.mark.integration
def test_actual_count_de_screen_runs_inside_nested_cv(data, config):
    pytest.importorskip("pydeseq2")
    local = replace(config, candidate_method="deseq2", strengths=(0.05,))
    audit = Audit()
    result = nested_validate(data, encode_outcome(data, local)[0], local, audit)
    assert result.metrics["out_of_fold"]["roc_auc"] > 0.8
    assert len(audit.fits) == local.outer_folds * (local.inner_folds + 1)
    assert all(len(fit["sample_ids"]) < len(data.expression) for fit in audit.fits)
    assert all(fit["candidate_method"] == "deseq2" for fit in audit.fits)


@pytest.mark.integration
def test_batch_holdout_hyperparameters_and_coefficients_are_preserved(data, config, tmp_path):
    output = tmp_path / "batch_run"
    local = replace(config, batch="batch", strengths=(0.05,), stability_resamples=2)
    run_analysis(data, local, output)
    metrics = json.loads((output / "model_metrics.json").read_text())
    batch = metrics["batch_validation"]
    assert batch["status"] == "tested"
    for fold in batch["folds"]:
        assert fold["parameters"]["strength"] in local.strengths
        assert fold["tuning_results"]
        assert "baseline_parameters" in fold
    with np.load(output / "resample_coefficients.npz", allow_pickle=False) as stored:
        assert stored["batch"].shape == (len(batch["folds"]), 1, data.expression.shape[1])
        assert stored["batch_labels"].tolist() == [row["batch"] for row in batch["folds"]]
