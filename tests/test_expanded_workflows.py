from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from sigtrellis.cell_features import distribution_features, feature_id
from sigtrellis.config import Config, load_config
from sigtrellis.de import differential_expression
from sigtrellis.domain import Audit, Dataset
from sigtrellis.panel import fit_panel
from sigtrellis.prediction import frozen_predict
from sigtrellis.preprocessing import Prepared
from sigtrellis.qc import encode_outcome, validate_dataset
from sigtrellis.simulate import synthetic_bulk
from sigtrellis.splits import group_values, make_splits
from sigtrellis.stability import stability_select
from sigtrellis.validation import nested_validate
from sigtrellis.workflow import run_analysis


def cell_fixture(path: Path, scenario: str = "abundance", samples: int = 48) -> Path:
    rng = np.random.default_rng(912)
    xs, metadata = [], []
    for si in range(samples):
        y = si % 2
        n = 80
        values = rng.poisson(8, (n, 12)).astype(float)
        if scenario == "abundance":
            proportion = np.clip((0.8 if y else 0.2) + rng.normal(0, 0.03), 0.05, 0.95)
            labels = np.where(np.arange(n) < int(n * proportion), "A", "B")
            values[:, 5] += (labels == "A") * 15
        else:
            labels = np.repeat("A", n)
            baseline = int(rng.integers(10, 20))
            values[:, :3] = baseline
            if y and scenario == "heterogeneity":
                values[: n // 2, :3] = 0
                values[n // 2 :, :3] = 2 * baseline
        # A depth-balancing background makes planted raw-count means exactly
        # identical for the two heterogeneity distributions, conditional on baseline.
        values[:, -1] = 1000 - values[:, :-1].sum(axis=1)
        xs.append(values)
        metadata.extend(
            {"sample_id": f"s{si:03d}", "donor": f"d{si:03d}", "phenotype": str(y), "state": state}
            for state in labels
        )
    obs = pd.DataFrame(metadata, index=[f"c{i:06d}" for i in range(samples * 80)])
    ad.AnnData(
        sparse.csr_matrix(np.vstack(xs)),
        obs=obs,
        var=pd.DataFrame(index=[f"g{i}" for i in range(12)]),
    ).write_h5ad(path)
    return path


def cell_config(**kwargs):
    return replace(
        Config(
            group="donor",
            cell_type="state",
            cell_states=("A", "B"),
            single_cell_mode="distribution",
            input_scale="features",
            normalization="none",
            imputation="median",
            min_cells=5,
            feature_genes=("g0", "g1", "g2"),
            max_features=30,
            programs={"response": ("g0", "g1", "g2")},
            feature_blocks=(
                "abundance",
                "gene_mean",
                "gene_detection",
                "gene_variance",
                "program_mean",
                "program_variance",
                "program_q90",
            ),
            outer_folds=3,
            inner_folds=2,
            strengths=(0.03, 0.12),
            l1_ratios=(0.5,),
            stability_resamples=4,
            permutations=0,
        ),
        **kwargs,
    )


def test_cell_features_are_exact_sample_local_summaries(tmp_path):
    path = cell_fixture(tmp_path / "cells.h5ad", samples=12)
    config = cell_config(cell_resamples=2)
    data = distribution_features(path, config)
    cells = ad.read_h5ad(path)
    selected = (cells.obs.sample_id == "s000") & (cells.obs.state == "A")
    raw = cells.X[selected.to_numpy()].toarray()
    log = np.log1p(raw / raw.sum(axis=1)[:, None] * 1e4)
    assert data.expression.loc["s000", feature_id("gene_mean", "A", "g0")] == pytest.approx(
        log[:, 0].mean()
    )
    assert data.expression.loc[
        "s000", feature_id("program_variance", "A", "response")
    ] == pytest.approx(log[:, :3].mean(axis=1).var(ddof=1))
    assert data.expression.loc["s000", feature_id("gene_detection", "A", "g0")] == pytest.approx(
        (raw[:, 0] > 0).mean()
    )
    assert len(data.cell_resamples) == 2
    assert all(frame.index.equals(data.expression.index) for frame in data.cell_resamples)
    assert validate_dataset(data, config)["n_biological_groups"] == 12


def test_cell_features_do_not_learn_from_other_donors_or_outcomes(tmp_path):
    path = cell_fixture(tmp_path / "cells.h5ad", samples=12)
    config = cell_config()
    original = distribution_features(path, config)
    cells = ad.read_h5ad(path)
    cells.obs["phenotype"] = np.where(cells.obs["phenotype"].astype(str) == "1", "0", "1")
    values = cells.X.toarray()
    values[cells.obs.sample_id == "s011"] *= 100
    cells.X = sparse.csr_matrix(values)
    cells.write_h5ad(tmp_path / "changed.h5ad")
    changed = distribution_features(tmp_path / "changed.h5ad", config)
    pd.testing.assert_frame_equal(original.expression.iloc[:-1], changed.expression.iloc[:-1])


def test_missing_state_is_unavailable_not_zero_expression(tmp_path):
    data = distribution_features(
        cell_fixture(tmp_path / "cells.h5ad", "heterogeneity", 12), cell_config()
    )
    assert data.expression[feature_id("gene_mean", "B", "g0")].isna().all()
    assert np.isfinite(data.expression[feature_id("abundance", "B")]).all()


@pytest.mark.parametrize(
    "scenario,blocks,expected",
    [
        ("abundance", ("abundance", "gene_mean"), "abundance"),
        ("heterogeneity", ("gene_variance", "program_variance"), "variance"),
    ],
)
def test_single_cell_specific_planted_signal_is_recovered(tmp_path, scenario, blocks, expected):
    config = cell_config(feature_blocks=blocks, cell_resamples=2)
    data = distribution_features(cell_fixture(tmp_path / "cells.h5ad", scenario), config)
    y, _ = encode_outcome(data, config)
    validate_dataset(data, config)
    audit = Audit()
    result = nested_validate(data, y, config, audit)
    assert result.metrics["out_of_fold"]["roc_auc"] > 0.9
    stable = stability_select(data, y, config, ["1_vs_0"], audit)
    assert any(
        expected in gene
        for gene in stable.table.loc[stable.table.selection_frequency >= 0.75, "gene_id"]
    )
    assert "frequency_cells_1" in stable.table
    for split in audit.splits:
        assert not set(split["train_groups"]) & set(split["test_groups"])


def test_cell_chunk_size_does_not_change_features_or_subsamples(tmp_path):
    path = cell_fixture(tmp_path / "cells.h5ad", samples=12)
    a = distribution_features(path, cell_config(cell_resamples=2, chunk_size=37))
    b = distribution_features(path, cell_config(cell_resamples=2, chunk_size=512))
    np.testing.assert_allclose(a.expression, b.expression, atol=1e-11, equal_nan=True)
    assert a.upstream_qc["sample_cell_hashes"] == b.upstream_qc["sample_cell_hashes"]
    for first, second in zip(a.cell_resamples, b.cell_resamples, strict=True):
        np.testing.assert_allclose(first, second, atol=1e-11, equal_nan=True)


def test_imputation_is_fitted_on_training_rows_only(config):
    metadata = pd.DataFrame(
        {
            "sample_id": ["a", "b", "c", "d"],
            "donor": ["a", "b", "c", "d"],
            "phenotype": [0, 1, 0, 1],
        }
    ).set_index("sample_id", drop=False)
    x = pd.DataFrame(
        {"a_feature": [1.0, 3.0, np.nan, 1000.0], "absent_train": [np.nan, np.nan, 12.0, 100.0]},
        index=metadata.index,
    )
    data = Dataset(x, metadata)
    local = replace(config, input_scale="features", normalization="none", imputation="median")
    prepared = Prepared(local).fit(
        data.subset(np.array([0, 1])), np.array([0.0, 1.0]), Audit(), "train"
    )
    assert prepared.normalizer.fill_values.tolist() == [2.0, 0.0]
    assert prepared.genes == ["a_feature"]
    transformed = prepared.transform(data.subset(np.array([2, 3])))
    assert transformed[0, 0] == pytest.approx(0)
    assert prepared.normalizer.fill_values.tolist() == [2.0, 0.0]


def test_compact_panel_selection_stays_within_outer_training(data, config):
    local = replace(config, panel_validation=True, panel_max_features=3, stability_resamples=3)
    y, _ = encode_outcome(data, local)
    audit = Audit()
    result = nested_validate(data, y, local, audit)
    assert len(result.metrics["panel"]["folds"]) == 3
    for fold in result.metrics["panel"]["folds"]:
        split = next(s for s in audit.splits if s["context"] == fold["context"])
        assert fold["n_features"] <= 3
        for fit in audit.fits:
            if fit["context"].startswith(fold["context"] + "/panel"):
                assert set(fit["sample_ids"]) <= set(split["train_samples"])
                assert not set(fit["sample_ids"]) & set(split["test_samples"])
    panel = fit_panel(data, y, local, Audit(), "full")
    import json

    # JSON serialization is part of the frozen numeric contract.
    state = json.loads(json.dumps(panel.model.to_state(), default=lambda v: v.tolist()))
    np.testing.assert_allclose(panel.model.predict(data), frozen_predict(state, data), atol=1e-12)


def test_temporal_splits_purge_donors_and_preserve_order(data, config):
    data.metadata["day"] = np.repeat(np.arange(12), 4)
    local = replace(config, cv_strategy="temporal", time="day", outer_folds=2, inner_folds=2)
    y, _ = encode_outcome(data, local)
    for split in make_splits(data, y, local):
        train, test = data.subset(split.train), data.subset(split.test)
        assert train.metadata.day.max() < test.metadata.day.min()
        assert not set(group_values(train, local)) & set(group_values(test, local))
        for inner in make_splits(train, y[split.train], local, inner=True):
            assert (
                train.metadata.day.iloc[inner.train].max()
                < train.metadata.day.iloc[inner.test].min()
            )


@pytest.mark.parametrize("scenario", ["continuous", "multiclass"])
def test_count_de_supports_continuous_and_multiclass(scenario, config):
    data, _ = synthetic_bulk(samples=60, features=30, seed=7, scenario=scenario)
    local = replace(config, outcome_type=scenario, candidate_method="deseq2")
    local.validate()
    y, _ = encode_outcome(data, local)
    result = differential_expression(data, y, local)
    assert result.contrast.nunique() == (2 if scenario == "multiclass" else 1)
    assert result.de_adjusted_pvalue.notna().sum() > 10
    assert (
        result.de_adjusted_pvalue.dropna()
        >= result.de_pvalue[result.de_adjusted_pvalue.notna()] - 1e-12
    ).all()
    if scenario == "continuous":
        assert (result.loc[["gene_0000", "gene_0001"], "de_log2_fold_change"] > 0).all()


def test_rich_single_cell_report_preserves_feature_identity(tmp_path):
    config = cell_config(
        panel_validation=True,
        panel_max_features=3,
        cell_resamples=1,
        feature_blocks=("abundance", "program_mean"),
    )
    data = distribution_features(cell_fixture(tmp_path / "cells.h5ad"), config)
    manifest = run_analysis(data, config, tmp_path / "result")
    assert manifest["status"] == "complete"
    table = pd.read_csv(tmp_path / "result/biomarkers.csv")
    assert {"gene_id", "feature_id", "feature_kind", "feature_unit"} <= set(table.columns)
    assert (tmp_path / "result/panel_state.json").exists()
    assert (tmp_path / "result/feature_schema.json").exists()
    assert not (tmp_path / "result/pseudobulk_counts.tsv.gz").exists()


@pytest.mark.parametrize(
    "changes",
    [
        {"program_thresholds": {"response": float("nan")}},
        {"feature_blocks": ("unknown",)},
        {"programs": {"bad": ("g0", "g0")}},
    ],
)
def test_malformed_feature_contract_fails(changes):
    with pytest.raises(ValueError):
        cell_config(**changes).validate()


def test_program_configuration_roundtrips(tmp_path):
    path = tmp_path / "run.yaml"
    path.write_text("programs:\n  generic_response: [organism_gene_A, organism_gene_B]\n")
    assert load_config(path).programs == {
        "generic_response": ("organism_gene_A", "organism_gene_B")
    }


def test_feature_names_cannot_collide():
    assert feature_id("gene_mean", "A|B", "C") != feature_id("gene_mean", "A", "B|C")
