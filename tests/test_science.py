from dataclasses import replace

import numpy as np
import pandas as pd
import pytest
from sklearn.linear_model import LogisticRegression

from sigtrellis.correlation import correlation_diagnostics
from sigtrellis.domain import Audit, write_json
from sigtrellis.metrics import group_weights
from sigtrellis.modeling import Hyperparameters, fit_model
from sigtrellis.prediction import frozen_predict
from sigtrellis.qc import encode_outcome, validate_dataset
from sigtrellis.simulate import synthetic_bulk
from sigtrellis.stability import jaccard_scores, stability_select, summarize_coefficients
from sigtrellis.validation import nested_validate, permutation_control
from sigtrellis.workflow import batch_diagnostics, evidence_gates


@pytest.mark.scientific
def test_signal_recovery_and_permutation_control(config):
    data, truth = synthetic_bulk(80, 200, seed=17)
    local = replace(config, stability_resamples=12, permutations=19, max_features=150)
    y, _ = encode_outcome(data, local)
    audit = Audit()
    val = nested_validate(data, y, local, audit)
    assert val.metrics["out_of_fold"]["roc_auc"] > 0.9
    perm = permutation_control(data, y, local, val.metrics["loss_improvement"], audit)
    assert perm["pvalue"] <= 0.05
    stability = stability_select(data, y, local, ["1_vs_0"], audit)
    assert set(truth) <= set(stability.table.query("selection_frequency >= 0.8").gene_id)


@pytest.mark.scientific
@pytest.mark.parametrize("seed", [3, 29, 101])
def test_pure_noise_has_no_credible_panel(config, seed):
    data, _ = synthetic_bulk(64, 120, seed=seed, scenario="noise")
    local = replace(config, permutations=19, stability_resamples=4, max_features=80)
    y, _ = encode_outcome(data, local)
    audit = Audit()
    val = nested_validate(data, y, local, audit)
    perm = permutation_control(data, y, local, val.metrics["loss_improvement"], audit)
    # Even a hypothetical perfectly stable selected feature cannot bypass the null gate.
    fabricated = pd.DataFrame(
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
    val.metrics["permutation"] = perm
    gated, blockers = evidence_gates(
        fabricated, local, validate_dataset(data, local), val.metrics, {"status": "not_assessed"}
    )
    assert not gated.passes_robustness_gates.any()
    assert blockers
    assert val.metrics["out_of_fold"]["roc_auc"] < 0.75


@pytest.mark.scientific
def test_batch_only_signal_cannot_pass(config):
    data, _ = synthetic_bulk(60, 100, seed=17, scenario="batch_confounded")
    local = replace(config, batch="batch")
    qc = validate_dataset(data, local)
    assert "perfect_batch_confounding" in qc["gate_blockers"]
    y, _ = encode_outcome(data, local)
    diagnostic, coefficients = batch_diagnostics(data, y, local, Audit())
    assert diagnostic["status"] == "incomplete"
    assert coefficients is None


@pytest.mark.scientific
def test_one_batch_biomarkers_fail_batch_robustness(config):
    data, _ = synthetic_bulk(90, 100, seed=3, scenario="batch_specific")
    local = replace(config, batch="batch")
    y, _ = encode_outcome(data, local)
    diagnostic, coefficients = batch_diagnostics(data, y, local, Audit())
    assert diagnostic["status"] == "tested"
    assert any(row["loss_improvement"] <= 0 for row in diagnostic["folds"])
    assert coefficients is not None


@pytest.mark.scientific
@pytest.mark.parametrize(
    "scenario", ["continuous", "multiclass", "imbalance", "outliers", "correlated"]
)
def test_supported_scientific_scenarios(config, scenario):
    data, truth = synthetic_bulk(72, 100, seed=23, scenario=scenario)
    local = replace(
        config, outcome_type=scenario if scenario in {"continuous", "multiclass"} else "binary"
    )
    qc = validate_dataset(data, local)
    y, labels = encode_outcome(data, local)
    result = nested_validate(data, y, local, Audit())
    assert result.metrics["loss_improvement"] > 0
    if scenario == "outliers":
        assert qc["outlier_samples"]
    if scenario == "continuous":
        assert result.metrics["out_of_fold"]["r2"] > 0.4
    if scenario == "multiclass":
        assert result.coefficients.shape[1] == 2
    if scenario == "correlated":
        assert truth


def test_coefficient_sign_rank_and_empty_jaccard():
    values = np.array([[[1.0, 0.0, -2.0]], [[0.0, 2.0, -3.0]], [[1.0, -1.0, 0.0]]])
    t = summarize_coefficients(values, ["a", "b", "c"], ["positive_vs_negative"], 1e-7)
    np.testing.assert_allclose(t.selection_frequency, [2 / 3] * 3)
    np.testing.assert_allclose(t.sign_consistency, [1, 0.5, 1])
    assert np.isnan(jaccard_scores(np.zeros((2, 3), dtype=bool))[0])


def test_correlated_substitution_preserves_group_evidence(config):
    data, _ = synthetic_bulk(48, 20, seed=11)
    data.expression.iloc[:, 1] = data.expression.iloc[:, 0]
    coefs = np.zeros((10, 1, 20))
    coefs[::2, 0, 0] = 1
    coefs[1::2, 0, 1] = 1
    ranked = summarize_coefficients(coefs, list(data.expression.columns), ["x"], 1e-7)
    result = correlation_diagnostics(data, coefs, ranked, config)
    assert result.membership["gene_0000"] == result.membership["gene_0001"]
    assert result.groups.selection_frequency.max() == 1
    assert result.substitutions.exclusive_given_any.max() == 1


@pytest.mark.parametrize("kind", ["binary", "multiclass", "continuous"])
def test_frozen_serialization_matches_sklearn(config, tmp_path, kind):
    data, _ = synthetic_bulk(48, 40, seed=12, scenario=kind if kind != "binary" else "signal")
    local = replace(config, outcome_type=kind)
    y, _ = encode_outcome(data, local)
    model = fit_model(data, y, local, Hyperparameters(0.05, 0.5), Audit(), "all")
    import json

    path = tmp_path / "model.json"
    write_json(path, model.to_state())
    np.testing.assert_allclose(
        frozen_predict(json.loads(path.read_text()), data), model.predict(data), atol=1e-10
    )


def test_solver_matches_declared_mean_loss_convention(data, config):
    y, _ = encode_outcome(data, config)
    model = fit_model(data, y, config, Hyperparameters(0.05, 0.5), Audit(), "all")
    x = model.prepared.transform(data)
    reference = LogisticRegression(
        solver="saga",
        l1_ratio=0.5,
        C=1 / (len(y) * 0.05),
        random_state=config.seed,
        max_iter=config.max_iter,
        tol=config.tolerance,
    ).fit(x, y)
    np.testing.assert_allclose(reference.coef_, model.estimator.coef_)


def test_reproducible_nested_predictions(data, config):
    y, _ = encode_outcome(data, config)
    a = nested_validate(data, y, config, Audit())
    b = nested_validate(data, y, config, Audit())
    pd.testing.assert_frame_equal(a.predictions, b.predictions)
    np.testing.assert_array_equal(a.coefficients, b.coefficients)


def test_group_equal_weights():
    weights = group_weights(np.array(["a", "a", "a", "b"]))
    assert weights[:3].sum() == weights[3]


def test_sign_disagreement_and_zero_final_fit_cannot_pass(config):
    table = pd.DataFrame(
        {
            "selection_frequency": [1.0] * 4,
            "sign_consistency": [1.0] * 4,
            "outer_selection_frequency": [1.0] * 4,
            "outer_sign_consistency": [1.0] * 4,
            "frequency_logcpm": [1.0] * 4,
            "coefficient": [1.0, -1.0, 0.0, 1.0],
            "coefficient_median_selected": [1.0] * 4,
            "outer_dominant_sign": [1.0] * 4,
            "batch_selection_frequency": [1.0] * 4,
            "batch_sign_consistency": [1.0] * 4,
            "batch_dominant_sign": [1.0, 1.0, 1.0, -1.0],
        }
    )
    result, blockers = evidence_gates(
        table,
        replace(config, permutations=99, batch="batch"),
        {"gate_blockers": []},
        {"permutation": {"status": "tested", "pvalue": 0.01}, "loss_improvement": 0.1},
        {"status": "tested", "folds": [{"loss_improvement": 0.1}]},
    )
    assert not blockers
    assert result.passes_robustness_gates.tolist() == [True, False, False, False]


def test_covariate_baseline_and_global_null_limit(data, config):
    rng = np.random.default_rng(91)
    data.metadata["age"] = rng.uniform(20, 80, len(data.expression))
    data.metadata["site"] = np.tile(["a", "b"], len(data.expression) // 2)
    local = replace(config, covariates=("age", "site"), permutations=2)
    y, _ = encode_outcome(data, local)
    audit = Audit()
    result = nested_validate(data, y, local, audit)
    assert any("baseline" in split["context"] for split in audit.splits)
    assert result.metrics["loss_improvement"] > 0
    perm = permutation_control(data, y, local, result.metrics["loss_improvement"], audit)
    assert perm["status"] == "global_only_with_covariates"
    model = fit_model(data, y, local, Hyperparameters(0.05, 0.5), Audit(), "all")
    np.testing.assert_allclose(frozen_predict(model.to_state(), data), model.predict(data))
    data.metadata.loc[data.metadata.index[0], "site"] = "never_seen"
    with pytest.raises(ValueError, match="absent from training"):
        model.predict(data)
