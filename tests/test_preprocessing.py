from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from sigtrellis.domain import Audit
from sigtrellis.preprocessing import Normalizer, Prepared
from sigtrellis.qc import encode_outcome


def test_logcpm_library_before_feature_subset():
    x = pd.DataFrame([[10, 90], [20, 180]], columns=["a", "b"])
    norm = Normalizer("logcpm").fit(x)
    out = norm.transform(x)
    np.testing.assert_allclose(out[0], out[1])
    np.testing.assert_allclose(out[0], np.log2(1 + np.array([100000, 900000])))


def test_frozen_reference_cannot_learn_test_distribution(data):
    train = data.expression.iloc[:30]
    test = data.expression.iloc[30:].copy()
    normalizer = Normalizer("median_ratio").fit(train)
    reference = normalizer.reference.copy()
    normalizer.transform(test * 10000)
    np.testing.assert_array_equal(reference, normalizer.reference)
    np.testing.assert_allclose(normalizer.transform(test.iloc[:1]), normalizer.transform(test)[:1])


def test_no_valid_ratio_reference():
    with pytest.raises(ValueError, match="reference"):
        Normalizer("median_ratio").fit(pd.DataFrame([[0, 1], [1, 0]]))


def test_normalization_does_not_use_outcome(data, config):
    y, _ = encode_outcome(data, config)
    p = Prepared(config).fit(data, y, Audit(), "one")
    q = Prepared(config).fit(data, y[::-1], Audit(), "two")
    np.testing.assert_array_equal(p.transform(data), q.transform(data))


def test_filter_scaler_train_only(data, config):
    y, _ = encode_outcome(data, config)
    train = data.subset(np.arange(24))
    test = data.subset(np.arange(24, 48))
    train.expression = train.expression.copy()
    test.expression = test.expression.copy()
    gene = train.expression.columns[0]
    train.expression[gene] = 0
    test.expression[gene] = 10000000
    prep = Prepared(config).fit(train, y[:24], Audit(), "training")
    assert gene not in prep.genes
    old = prep.scaler.mean_.copy()
    prep.transform(test)
    np.testing.assert_array_equal(prep.scaler.mean_, old)


def test_unsupported_gene_universe_rejected(data):
    n = Normalizer("logcpm").fit(data.expression)
    with pytest.raises(ValueError, match="universe"):
        n.transform(data.expression.iloc[:, :-1])


def test_covariate_encoding_no_metadata_predictors(data, config):
    data.metadata["age"] = np.arange(len(data.expression), dtype=float)
    data.metadata["site"] = np.tile(["A", "B"], 24)
    local = replace(config, covariates=("age", "site"))
    y, _ = encode_outcome(data, local)
    prep = Prepared(local).fit(data, y, Audit(), "test")
    assert prep.transform(data).shape[1] == len(prep.genes) + 2
    assert prep.covariates.columns == ["covariate:age", "covariate:site=B"]
