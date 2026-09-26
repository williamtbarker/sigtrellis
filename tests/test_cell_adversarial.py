from __future__ import annotations

import json
from dataclasses import replace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from test_expanded_workflows import cell_config, cell_fixture

from sigtrellis.cell_features import distribution_features
from sigtrellis.config import load_config
from sigtrellis.domain import Audit
from sigtrellis.prediction import external_validate
from sigtrellis.qc import encode_outcome, validate_dataset
from sigtrellis.singlecell import pseudobulk
from sigtrellis.validation import nested_validate
from sigtrellis.workflow import run_analysis


def joint_fixture(path, samples=48, seed=314):
    """Identical per-gene marginal mechanisms; only within-cell coupling has signal."""
    rng = np.random.default_rng(seed)
    values, rows = [], []
    for i in range(samples):
        n = 80
        y = i % 2
        a = np.repeat([10, 60], n // 2) + rng.integers(0, 10, n)
        b = np.sort(rng.integers(10, 70, n))
        if not y:
            b = b[::-1]
        x = rng.poisson(8, (n, 10))
        x[:, 0], x[:, 1] = a, b
        x[:, 2], x[:, 3] = 100 - a, 100 - b
        x[:, -1] = 1000 - x[:, :-1].sum(axis=1)
        values.append(x)
        rows.extend(
            {"sample_id": f"s{i:03d}", "donor": f"d{i:03d}", "phenotype": str(y), "state": "A"}
            for _ in range(n)
        )
    obs = pd.DataFrame(rows, index=[f"cell_{i}" for i in range(samples * 80)])
    ad.AnnData(
        sparse.csr_matrix(np.vstack(values)),
        obs=obs,
        var=pd.DataFrame(index=[f"g{i}" for i in range(10)]),
    ).write_h5ad(path)
    return path


def joint_config(**kw):
    return cell_config(
        cell_states=("A",),
        feature_genes=("g0", "g1"),
        gene_pairs=(("g0", "g1"), ("g2", "g3")),
        programs=kw.pop("programs", {"first": ("g0",), "second": ("g1",)}),
        feature_blocks=kw.pop("feature_blocks", ("gene_correlation", "program_correlation")),
        **kw,
    )


def test_reordering_cells_preserves_perturbations_and_fingerprints(tmp_path):
    path = cell_fixture(tmp_path / "cells.h5ad", samples=12)
    config = cell_config(cell_resamples=2)
    original = distribution_features(path, config)
    cells = ad.read_h5ad(path)
    perm = np.random.default_rng(15).permutation(cells.n_obs)
    cells[perm, :].copy()[:, np.arange(cells.n_vars)[::-1]].copy().write_h5ad(
        tmp_path / "reordered.h5ad"
    )
    reordered = distribution_features(tmp_path / "reordered.h5ad", config)
    assert original.upstream_qc["sample_cell_hashes"] == reordered.upstream_qc["sample_cell_hashes"]
    for a, b in zip(
        (original.expression, *original.cell_resamples),
        (reordered.expression, *reordered.cell_resamples),
        strict=True,
    ):
        np.testing.assert_allclose(a, b, equal_nan=True, atol=1e-10)


def test_unrelated_donor_removal_does_not_change_cell_resampling(tmp_path):
    path = cell_fixture(tmp_path / "cells.h5ad", samples=12)
    config = cell_config(cell_resamples=2)
    original = distribution_features(path, config)
    cells = ad.read_h5ad(path)
    cells[cells.obs.sample_id != "s000"].copy().write_h5ad(tmp_path / "subset.h5ad")
    subset = distribution_features(tmp_path / "subset.h5ad", config)
    for a, b in zip(original.cell_resamples, subset.cell_resamples, strict=True):
        np.testing.assert_allclose(a.loc[b.index], b, equal_nan=True, atol=1e-10)


def test_renamed_reordered_cell_copy_is_rejected(tmp_path):
    path = cell_fixture(tmp_path / "cells.h5ad", samples=12)
    cells = ad.read_h5ad(path)
    original = cells[cells.obs.sample_id == "s000"].copy()
    original.obs["sample_id"] = "copy"
    original.obs["donor"] = "copy_donor"
    original.obs_names = [f"new_{i}" for i in range(len(original))]
    combined = ad.concat([cells, original[::-1]], index_unique="-")
    combined.write_h5ad(tmp_path / "duplicate.h5ad")
    data = distribution_features(tmp_path / "duplicate.h5ad", cell_config())
    with pytest.raises(ValueError, match="Identical expression"):
        validate_dataset(data, cell_config())


def test_tail_and_joint_features_match_direct_cell_calculations(tmp_path):
    path = joint_fixture(tmp_path / "cells.h5ad", samples=12)
    config = joint_config(feature_blocks=("gene_tail", "gene_correlation", "program_correlation"))
    data = distribution_features(path, config)
    cells = ad.read_h5ad(path)
    raw = cells.X[cells.obs.sample_id.to_numpy() == "s000"].toarray()
    log = np.log1p(raw / raw.sum(axis=1)[:, None] * 1e4)
    for feature, spec in data.features.items():
        if spec.kind == "gene_tail":
            expected = (log[:, int(spec.gene_id[1:])] > spec.threshold).mean()
        elif spec.kind == "gene_correlation":
            expected = np.corrcoef(
                log[:, int(spec.gene_id[1:])], log[:, int(spec.gene_partner[1:])]
            )[0, 1]
        else:
            expected = np.corrcoef(log[:, 0], log[:, 1])[0, 1]
        assert data.expression.loc["s000", feature] == pytest.approx(expected)


def test_coupling_signal_survives_when_marginal_summaries_cannot_predict(tmp_path):
    path = joint_fixture(tmp_path / "cells.h5ad")
    config = joint_config()
    joint = distribution_features(path, config)
    y, _ = encode_outcome(joint, config)
    validate_dataset(joint, config)
    audit = Audit()
    full = nested_validate(joint, y, config, audit)
    marginal_config = joint_config(feature_blocks=("gene_mean", "gene_variance", "gene_detection"))
    marginal = distribution_features(path, marginal_config)
    baseline = nested_validate(marginal, y, marginal_config, Audit())
    assert full.metrics["out_of_fold"]["roc_auc"] > 0.95
    assert baseline.metrics["out_of_fold"]["roc_auc"] < 0.75
    bulk_config = replace(
        config,
        single_cell_mode="pseudobulk",
        input_scale="counts",
        normalization="logcpm",
        imputation="reject",
        cell_type_value="A",
    )
    bulk = pseudobulk(path, bulk_config)[0]
    bulk_validation = nested_validate(bulk, y, bulk_config, Audit())
    assert bulk_validation.metrics["out_of_fold"]["roc_auc"] < 0.75
    assert all(not set(s["train_groups"]) & set(s["test_groups"]) for s in audit.splits)


def test_single_correlation_feature_can_run_complete_report(tmp_path):
    data = distribution_features(
        joint_fixture(tmp_path / "cells.h5ad", 24),
        joint_config(feature_blocks=("program_correlation",)),
    )
    config = joint_config(
        feature_blocks=("program_correlation",), stability_resamples=2, outer_folds=2
    )
    manifest = run_analysis(data, config, tmp_path / "result")
    assert manifest["status"] == "complete"
    assert manifest["n_features"] == 1
    assert (tmp_path / "result/figures/cell_distributions.png").is_file()
    table = pd.read_csv(tmp_path / "result/cell_distribution_evidence.csv")
    assert len(table) == 24
    assert table.group_id.nunique() == 24
    assert "within_cell_pearson_correlation" in table.feature_unit.to_list()


def test_low_information_correlations_are_unavailable(tmp_path):
    path = joint_fixture(tmp_path / "cells.h5ad", 12)
    cells = ad.read_h5ad(path)
    x = cells.X.toarray()
    # Constant genes within a sample, with equal cell totals.
    x[:, 0] = 25
    x[:, -1] = 1000 - x[:, :-1].sum(axis=1)
    cells.X = sparse.csr_matrix(x)
    cells.write_h5ad(path)
    config = joint_config(feature_blocks=("gene_correlation",))
    data = distribution_features(path, config)
    first = [k for k, v in data.features.items() if v.gene_id == "g0"][0]
    assert data.expression[first].isna().all()


def test_external_requires_original_cell_normalization_universe(tmp_path):
    path = cell_fixture(tmp_path / "train.h5ad", samples=12)
    config = cell_config(stability_resamples=2, outer_folds=2)
    data = distribution_features(path, config)
    run_analysis(data, config, tmp_path / "run")
    cells = ad.read_h5ad(path)
    cells.obs["sample_id"] = "external_" + cells.obs.sample_id.astype(str)
    cells.obs["donor"] = "external_" + cells.obs.donor.astype(str)
    # Removing an unmodeled background gene changes every cell's library size.
    cells[:, :-1].copy().write_h5ad(tmp_path / "external.h5ad")
    with pytest.raises(ValueError, match="normalization/identity contract"):
        external_validate(
            tmp_path / "run",
            None,
            None,
            tmp_path / "external",
            single_cell=tmp_path / "external.h5ad",
        )


def test_permuted_perturbation_metadata_cannot_cross_boundary(tmp_path):
    config = cell_config(cell_resamples=1)
    data = distribution_features(cell_fixture(tmp_path / "cells.h5ad", samples=12), config)
    data.cell_resamples = (data.cell_resamples[0].iloc[::-1],)
    with pytest.raises(ValueError, match="perturbation rows/columns"):
        validate_dataset(data, config)


def test_program_moments_do_not_retain_cell_arrays(tmp_path):
    config = cell_config(feature_blocks=("program_mean", "program_variance"))
    data = distribution_features(cell_fixture(tmp_path / "cells.h5ad", samples=12), config)
    assert data.upstream_qc["retained_program_scores"] is False
    assert not data.counts_by_cell_type


@pytest.mark.parametrize(
    "changes",
    [
        {"gene_pairs": (("g0", "g0"),)},
        {"gene_pairs": (("g0", "g1"), ("g1", "g0"))},
        {"gene_thresholds": (float("nan"),)},
        {"gene_thresholds": (True,)},
        {"gene_thresholds": (1.0, 1.0)},
        {"cell_type_value": "A"},
        {"feature_blocks": ("gene_correlation",), "gene_pairs": ()},
        {"feature_blocks": ("program_correlation",), "programs": {}},
    ],
)
def test_joint_configuration_rejects_ambiguous_contract(changes):
    with pytest.raises(ValueError):
        cell_config(**changes).validate()


def test_joint_configuration_round_trip(tmp_path):
    config = joint_config(feature_blocks=("gene_tail", "gene_correlation", "program_correlation"))
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config.to_dict()))
    assert load_config(path) == config


def test_streamed_public_preparation_preserves_each_raw_cell(tmp_path):
    import importlib.util
    from pathlib import Path

    path = Path(__file__).resolve().parents[1] / "examples/prepare_replicated.py"
    spec = importlib.util.spec_from_file_location("prep_streamed", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    cells = ad.read_h5ad(joint_fixture(tmp_path / "source.h5ad", 12))
    rows = np.arange(0, cells.n_obs, 3)
    module.write_count_subset(
        cells.X, rows, cells.obs.iloc[rows], cells.var, tmp_path / "prepared.h5ad", "test subset"
    )
    prepared = ad.read_h5ad(tmp_path / "prepared.h5ad")
    np.testing.assert_array_equal(prepared.X.toarray(), cells.X[rows].toarray())
    pd.testing.assert_frame_equal(prepared.obs, cells.obs.iloc[rows])
    assert prepared.var_names.equals(cells.var_names)


def test_rare_tail_signal_is_visible_without_a_bulk_mean_difference(tmp_path):
    path = joint_fixture(tmp_path / "rare.h5ad")
    cells = ad.read_h5ad(path)
    values = cells.X.toarray()
    for i, sample in enumerate(cells.obs.sample_id.unique()):
        indices = np.flatnonzero(cells.obs.sample_id.to_numpy() == sample)
        values[indices, 0] = 5
        if i % 2:
            # Four of 80 cells express 100 counts: same count sum as 5 per cell.
            values[indices, 0] = 0
            values[indices[:4], 0] = 100
    values[:, -1] = 1000 - values[:, :-1].sum(axis=1)
    cells.X = sparse.csr_matrix(values)
    cells.write_h5ad(path)
    config = joint_config(feature_blocks=("gene_tail",), gene_thresholds=(5.0,))
    data = distribution_features(path, config)
    y, _ = encode_outcome(data, config)
    result = nested_validate(data, y, config, Audit())
    assert result.metrics["out_of_fold"]["roc_auc"] > 0.95
    for sample in cells.obs.sample_id.unique():
        assert values[cells.obs.sample_id.to_numpy() == sample, 0].sum() == 400


def test_filter_audit_is_reconstructable_without_repeating_excluded_features(
    data, config, tmp_path
):
    from sigtrellis.domain import identifier_hash
    from sigtrellis.preprocessing import Prepared

    audit = Audit()
    prepared = Prepared(replace(config, max_features=3)).fit(
        data, encode_outcome(data, config)[0], audit, "test"
    )
    record = audit.fits[0]
    assert record["feature_universe_sha256"] == identifier_hash(data.expression.columns.tolist())
    assert record["removed_gene_count"] == len(set(data.expression) - set(record["retained_genes"]))
    assert record["retained_genes"] == prepared.genes
    assert "removed_genes" not in record


def test_program_pair_identity_survives_sorted_json_configuration(tmp_path):
    path = joint_fixture(tmp_path / "cells.h5ad", 12)
    config = joint_config(
        programs={"z_program": ("g0",), "a_program": ("g1",)},
        feature_blocks=("program_correlation",),
    )
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config.to_dict(), sort_keys=True))
    first = distribution_features(path, config)
    restored = distribution_features(path, load_config(config_path))
    pd.testing.assert_frame_equal(first.expression, restored.expression)
    assert first.features == restored.features
