import json
from dataclasses import replace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.io import mmwrite
from test_expanded_workflows import cell_config, cell_fixture

from sigtrellis.cell_features import distribution_features
from sigtrellis.config import Config
from sigtrellis.domain import Audit
from sigtrellis.matrix_import import import_matrix_market
from sigtrellis.prediction import external_validate
from sigtrellis.qc import encode_outcome
from sigtrellis.simulate import synthetic_bulk
from sigtrellis.validation import nested_validate
from sigtrellis.workflow import panel_feature_gates, run_analysis


def test_matrix_import_preserves_ids_and_explicit_numeric_types(tmp_path):
    mmwrite(tmp_path / "matrix.mtx", np.array([[2, 4], [5, 8]]))
    (tmp_path / "genes.tsv").write_text("g0\tG0\ng1\tG1\n")
    (tmp_path / "barcodes.tsv").write_text("001\n002\n")
    (tmp_path / "obs.tsv").write_text("cell\tdonor\tresponse\n002\t020\t1.5\n001\t010\t2.5\n")
    args = [
        tmp_path / name
        for name in ("matrix.mtx", "genes.tsv", "barcodes.tsv", "obs.tsv", "out.h5ad")
    ]
    import_matrix_market(*args, numeric_columns=("response",))
    data = ad.read_h5ad(args[-1])
    assert data.obs.index.tolist() == ["001", "002"]
    assert data.obs.donor.tolist() == ["010", "020"]
    assert data.obs.response.tolist() == [2.5, 1.5]
    np.testing.assert_array_equal(data.X.toarray(), [[2, 5], [4, 8]])
    with pytest.raises(ValueError, match="overwrite"):
        import_matrix_market(*args)
    args[-1] = tmp_path / "bad.h5ad"
    (tmp_path / "barcodes.tsv").write_text("001\n003\n")
    with pytest.raises(ValueError, match="match"):
        import_matrix_market(*args)


def test_frozen_native_single_cell_panel_on_untouched_donors(tmp_path):
    source = ad.read_h5ad(cell_fixture(tmp_path / "all.h5ad"))
    training = source.obs.sample_id.astype(str) < "s032"
    source[training].copy().write_h5ad(tmp_path / "train.h5ad")
    source[~training].copy().write_h5ad(tmp_path / "external.h5ad")
    config = cell_config(
        feature_blocks=("abundance", "program_mean"),
        panel_validation=True,
        panel_max_features=2,
        outer_folds=2,
        stability_resamples=2,
        strengths=(0.03,),
    )
    data = distribution_features(tmp_path / "train.h5ad", config)
    run = tmp_path / "run"
    run_analysis(data, config, run)
    before = (run / "panel_state.json").read_bytes()
    result = external_validate(
        run,
        None,
        None,
        tmp_path / "evaluation",
        single_cell=tmp_path / "external.h5ad",
        model_kind="panel",
    )
    assert result["training_transform_refit"] is False
    assert result["model_kind"] == "panel"
    assert result["cell_profile_duplicate_check"]
    assert result["metrics"]["roc_auc"] > 0.9
    assert (run / "panel_state.json").read_bytes() == before
    with pytest.raises(ValueError, match="overlap"):
        external_validate(
            run,
            None,
            None,
            tmp_path / "bad",
            single_cell=tmp_path / "train.h5ad",
            model_kind="panel",
        )


def test_compact_panel_noise_control_cannot_pass(config, tmp_path):
    data, _ = synthetic_bulk(samples=48, features=50, seed=27, scenario="noise")
    local = replace(
        config,
        outer_folds=2,
        inner_folds=2,
        stability_resamples=2,
        panel_validation=True,
        panel_max_features=3,
        permutations=19,
        strengths=(0.1,),
        l1_ratios=(0.8,),
    )
    manifest = run_analysis(data, local, tmp_path / "noise")
    metrics = json.loads((tmp_path / "noise/model_metrics.json").read_text())
    assert not manifest["candidate_features"]
    assert metrics["permutation"]["panel_pvalue"] > 0.05
    assert len(metrics["permutation"]["panel_null_improvements"]) == 19


def test_temporal_nested_predictions_are_future_samples(data, config):
    data.metadata["time"] = np.repeat(np.arange(12), 4)
    local = replace(
        config, cv_strategy="temporal", time="time", outer_folds=2, inner_folds=2, strengths=(0.05,)
    )
    y, _ = encode_outcome(data, local)
    audit = Audit()
    result = nested_validate(data, y, local, audit)
    assert len(result.predictions) == 24
    assert all(s["train_time_max"] < s["test_time_min"] for s in audit.splits)


def test_program_only_features_do_not_allocate_gene_accumulators(tmp_path):
    data = distribution_features(
        cell_fixture(tmp_path / "cells.h5ad", samples=12),
        cell_config(feature_blocks=("abundance", "program_variance"), max_dense_mb=1),
    )
    assert not data.counts_by_cell_type
    assert data.expression.shape == (12, 4)


def test_single_cell_null_does_not_become_signal(tmp_path):
    data = distribution_features(
        cell_fixture(tmp_path / "cells.h5ad", "noise", 60),
        cell_config(feature_blocks=("gene_mean", "gene_variance")),
    )
    config = cell_config(feature_blocks=("gene_mean", "gene_variance"))
    y, _ = encode_outcome(data, config)
    result = nested_validate(data, y, config, Audit())
    assert 0.2 < result.metrics["out_of_fold"]["roc_auc"] < 0.8


def test_public_preparation_purges_aliquots_of_external_donors():
    import importlib.util
    from pathlib import Path

    path = Path(__file__).resolve().parents[1] / "examples/prepare_replicated.py"
    spec = importlib.util.spec_from_file_location("preparation", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    source = pd.DataFrame(
        {
            "sample_uuid": ["same", "same", "unseen", "train", "train"],
            "donor_id": ["d0", "d0", "d1", "d2", "d2"],
            "Processing_Cohort": ["2.0", "4.0", "4.0", "2.0", "3.0"],
            "disease": ["case", "case", "control", "control", "control"],
        }
    )
    obs, train, external = module.lupus_partitions(source)
    assert set(train.donor_id) == {"d2"}
    assert set(external.donor_id) == {"d0", "d1"}
    assert len(train) == 1
    assert obs.specimen_id.nunique() == 5
    assert set(external.Processing_Cohort) == {"4.0"}


def test_panel_cannot_inherit_confidence_for_zero_or_reversed_weights(config):
    table = pd.DataFrame(
        {
            "passes_robustness_gates": [True] * 4,
            "panel_coefficient": [1.0, 0.0, -1.0, 1.0],
            "coefficient_median_selected": [1.0] * 4,
            "panel_outer_selection_frequency": [1.0, 1.0, 1.0, 0.0],
            "panel_outer_sign_consistency": [1.0] * 4,
            "panel_outer_dominant_sign": [1.0] * 4,
        }
    )
    assert panel_feature_gates(table, config).tolist() == [True, False, False, False]


def test_unsupported_cell_perturbation_mode_fails_explicitly():
    with pytest.raises(ValueError, match="distribution"):
        Config(cell_resamples=1).validate()


def test_feature_csv_cannot_silently_omit_declared_cell_perturbations(data, tmp_path):
    with pytest.raises(ValueError, match="single-cell adapter"):
        run_analysis(data, cell_config(cell_resamples=1), tmp_path / "invalid")
    assert json.loads((tmp_path / "invalid/run_manifest.json").read_text())["status"] == "failed"


@pytest.mark.parametrize("column", ["phenotype", "sample_id", "donor", "batch", "time", "stratum"])
def test_cell_annotation_cannot_reencode_protected_metadata(column):
    with pytest.raises(ValueError, match="Cell-type annotation"):
        cell_config(
            cell_type=column, batch="batch", time="time", permutation_strata="stratum"
        ).validate()


def test_abundance_does_not_change_when_the_same_cell_distribution_is_repeated(tmp_path):
    source = ad.read_h5ad(cell_fixture(tmp_path / "cells.h5ad", samples=12))
    repeated = source[np.repeat(np.arange(source.n_obs), 3)].copy()
    repeated.obs_names = [f"replicated_cell_{i}" for i in range(repeated.n_obs)]
    repeated.write_h5ad(tmp_path / "repeated.h5ad")
    config = cell_config(feature_blocks=("abundance",))
    original = distribution_features(tmp_path / "cells.h5ad", config)
    changed = distribution_features(tmp_path / "repeated.h5ad", config)
    np.testing.assert_allclose(original.expression, changed.expression, atol=1e-12)


def test_capture_depth_alone_cannot_manufacture_an_abundance_biomarker(tmp_path):
    source = ad.read_h5ad(cell_fixture(tmp_path / "cells.h5ad", "noise", 24))
    # Both outcomes consist exclusively of state A, with no state B. Only the
    # number of recovered cells differs systematically (20 versus 80).
    keep = (source.obs.phenotype.astype(str) == "1").to_numpy() | (
        np.arange(source.n_obs) % 80 < 20
    )
    source[keep].copy().write_h5ad(tmp_path / "depth.h5ad")
    config = cell_config(feature_blocks=("abundance",))
    data = distribution_features(tmp_path / "depth.h5ad", config)
    assert (data.expression.nunique() == 1).all()
    y, _ = encode_outcome(data, config)
    result = nested_validate(data, y, config, Audit())
    assert result.metrics["out_of_fold"]["roc_auc"] == pytest.approx(0.5)
    assert not np.any(result.coefficients)


def test_predictor_dictionary_cannot_shrink_count_de_normalization_universe(tmp_path):
    path = cell_fixture(tmp_path / "cells.h5ad", samples=12)
    config = cell_config(
        feature_blocks=("gene_mean",), feature_genes=("g0", "g1"), supporting_de=True
    )
    data = distribution_features(path, config)
    assert data.expression.shape[1] == 4  # two predictive genes in two states
    assert all(frame.shape[1] == 12 for frame in data.counts_by_cell_type.values())
    source = ad.read_h5ad(path)
    for state, counts in data.counts_by_cell_type.items():
        keep = ((source.obs.sample_id == "s000") & (source.obs.state == state)).to_numpy()
        np.testing.assert_array_equal(
            counts.loc["s000"], np.asarray(source.X[keep].sum(axis=0)).ravel()
        )
