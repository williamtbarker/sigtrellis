from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from sigtrellis.domain import Audit
from sigtrellis.preprocessing import Prepared
from sigtrellis.qc import encode_outcome, validate_dataset
from sigtrellis.simulate import synthetic_bulk
from sigtrellis.splits import assert_disjoint, make_splits, permute_outcome
from sigtrellis.validation import nested_validate


def test_duplicates_cannot_cross_folds(data, config):
    data.expression.iloc[1] = data.expression.iloc[0]
    with pytest.raises(ValueError, match="Identical expression"):
        validate_dataset(data, config)


def test_direct_outcome_feature_blocks_claims(data, config):
    data.expression["covert_label"] = 10 * data.metadata.phenotype.to_numpy(dtype=float) + 5
    qc = validate_dataset(data, config)
    assert "covert_label" in qc["possible_outcome_encoded_features"]
    assert "possible_outcome_encoding_requires_review" in qc["gate_blockers"]


def test_low_count_coincidences_are_not_false_encoding_errors(data, config):
    data.expression["low_count_gene"] = data.metadata.phenotype.to_numpy(dtype=float)
    qc = validate_dataset(data, config)
    assert "low_count_gene" not in qc["possible_outcome_encoded_features"]


def test_batch_encodes_phenotype_blocks_confidence(data, config):
    data.metadata["batch"] = data.metadata.phenotype.astype(str)
    qc = validate_dataset(data, replace(config, batch="batch"))
    assert "perfect_batch_confounding" in qc["gate_blockers"]
    with pytest.raises(ValueError, match="cannot be predictors"):
        replace(config, batch="batch", covariates=("batch",)).validate()


def test_covariate_label_copy_rejected(data, config):
    data.metadata["secret"] = data.metadata.phenotype * 10 + 3
    with pytest.raises(ValueError, match="affine copy"):
        validate_dataset(data, replace(config, covariates=("secret",)))


def test_identifier_labels_never_predict(data, config):
    y, _ = encode_outcome(data, config)
    before = Prepared(config).fit(data, y, Audit(), "before").transform(data)
    ids = [f"condition_{v}_sample_{i}" for i, v in enumerate(y)]
    data.expression.index = ids
    data.metadata.index = ids
    data.metadata["sample_id"] = ids
    data.metadata["donor"] = ids
    after = Prepared(config).fit(data, y, Audit(), "after").transform(data)
    np.testing.assert_array_equal(before, after)


def test_every_supervised_fit_is_inside_its_training_boundary(data, config, monkeypatch):
    audit = Audit()
    calls = []
    original = Prepared.fit

    def guarded(self, subset, y, active_audit, context):
        parent = context.removesuffix("/refit")
        matches = [s for s in active_audit.splits if s["context"] == parent]
        assert matches, f"Unaudited or global preprocessing at {context}"
        split = matches[-1]
        assert set(subset.expression.index) == set(split["train_samples"])
        assert not set(subset.expression.index) & set(split["test_samples"])
        calls.append(context)
        return original(self, subset, y, active_audit, context)

    monkeypatch.setattr(Prepared, "fit", guarded)
    nested_validate(
        data,
        encode_outcome(data, config)[0],
        replace(config, candidate_method="association"),
        audit,
    )
    assert len(calls) == config.outer_folds * (config.inner_folds + 1)


def test_de_is_recomputed_for_each_training_split(data, config, monkeypatch):
    import sigtrellis.preprocessing as preprocessing

    observed = []

    def fake_de(subset, y, local):
        observed.append(set(subset.expression.index))
        frame = pd.DataFrame(
            {"de_adjusted_pvalue": np.linspace(0.001, 0.09, subset.expression.shape[1])},
            index=subset.expression.columns,
        )
        frame.attrs["warnings"] = []
        return frame

    monkeypatch.setattr(preprocessing, "differential_expression", fake_de)
    audit = Audit()
    local = replace(config, candidate_method="deseq2")
    nested_validate(data, encode_outcome(data, local)[0], local, audit)
    assert len(observed) == local.outer_folds * (local.inner_folds + 1)
    assert all(len(ids) < len(data.expression) for ids in observed)
    valid_training_sets = [set(s["train_samples"]) for s in audit.splits]
    assert all(ids in valid_training_sets for ids in observed)


def test_donor_overlap_is_rejected(data, config):
    data.metadata.iloc[1, data.metadata.columns.get_loc("donor")] = data.metadata.donor.iloc[0]
    with pytest.raises(ValueError, match="Biological groups"):
        assert_disjoint(data, np.array([0]), np.array([1]), config)


def test_grouped_and_leave_group_splits(data, config):
    y, _ = encode_outcome(data, config)
    for strategy in ("grouped", "leave_group_out"):
        for split in make_splits(data, y, replace(config, cv_strategy=strategy)):
            assert_disjoint(data, split.train, split.test, config)


def test_batch_holdout_purges_shared_donors(data, config):
    data.metadata["batch"] = np.repeat(["a", "b", "c"], 16)
    data.metadata.loc[data.metadata.index[16], "donor"] = data.metadata.donor.iloc[0]
    local = replace(config, batch="batch", cv_strategy="leave_batch_out")
    for split in make_splits(data, encode_outcome(data, local)[0], local):
        assert_disjoint(data, split.train, split.test, local)
        assert len(set(data.metadata.batch.iloc[split.test])) == 1
        assert not set(data.metadata.batch.iloc[split.test]) & set(
            data.metadata.batch.iloc[split.train]
        )
        if 0 in split.test:
            assert 16 not in split.train


def test_permutation_respects_declared_strata(data, config):
    data.metadata["site"] = np.repeat(["first", "second"], 24)
    local = replace(config, permutation_strata="site")
    y, _ = encode_outcome(data, local)
    perm = permute_outcome(data, y, local, np.random.default_rng(91))
    for site in data.metadata.site.unique():
        rows = data.metadata.site == site
        assert sorted(y[rows]) == sorted(perm[rows])
    assert not np.array_equal(y, perm)


def test_permutation_preserves_biological_unit(data, config):
    # Two rows per donor with a constant phenotype.
    y, _ = encode_outcome(data, config)
    data.metadata["donor"] = [f"d{i // 2}" for i in range(48)]
    y = np.repeat(np.tile([0.0, 1.0], 12), 2)
    perm = permute_outcome(data, y, config, np.random.default_rng(3))
    assert np.all(perm[::2] == perm[1::2])
    assert sorted(perm) == sorted(y)
    paired = np.tile([0.0, 1.0], 24)
    with pytest.raises(ValueError, match="constant outcome"):
        permute_outcome(data, paired, config, np.random.default_rng(1))
    perm = permute_outcome(
        data, paired, replace(config, permutation_scheme="within_group"), np.random.default_rng(1)
    )
    assert np.all(perm.reshape(-1, 2).sum(axis=1) == 1)


@pytest.mark.scientific
def test_global_feature_screening_positive_leak_control(config):
    data, _ = synthetic_bulk(60, 2000, seed=31, scenario="noise")
    # Independent continuous noise removes count-composition issues from this attack.
    rng = np.random.default_rng(22)
    data.expression[:] = rng.normal(size=data.expression.shape)
    local = replace(
        config,
        input_scale="log_expression",
        normalization="none",
        candidate_method="association",
        max_features=5,
    )
    y, _ = encode_outcome(data, local)
    safe = nested_validate(data, y, local, Audit()).metrics["out_of_fold"]["roc_auc"]
    global_screen = Prepared(local).fit(data, y, Audit(), "deliberately_unsafe")
    data.expression = data.expression.loc[:, global_screen.genes]
    unsafe = nested_validate(data, y, replace(local, candidate_method="none"), Audit()).metrics[
        "out_of_fold"
    ]["roc_auc"]
    assert safe < 0.75
    assert unsafe > safe + 0.12
