from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from sigtrellis.config import load_config
from sigtrellis.io import load_bulk, read_table
from sigtrellis.qc import encode_outcome, validate_dataset


@pytest.mark.parametrize(
    "change",
    [
        {"outcome_type": "invalid"},
        {"normalization": "phenotype_scaled"},
        {"outer_folds": 1},
        {"l1_ratios": (0,)},
        {"strengths": (float("nan"),)},
        {"covariates": ("phenotype",)},
        {"input_scale": "log_expression"},
        {"stability_fraction": 1.0},
        {"permutations": -1},
        {"candidate_method": "global_de"},
        {"cv_strategy": "leave_batch_out"},
        {"supporting_de": True, "input_scale": "features", "normalization": "none"},
        {"de_pair_group": True, "group": None},
        {"min_count": 0},
        {"stability_normalizations": ("none",)},
        {"stability_normalizations": ("logcpm", "logcpm")},
        {"sample_id": None},
        {"outcome": None},
        {"covariates": ("",)},
    ],
)
def test_bad_config(config, change):
    with pytest.raises(ValueError):
        replace(config, **change).validate()


def test_config_unknown_and_formats(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text("unknown: 1\n")
    with pytest.raises(ValueError, match="Unknown"):
        load_config(p)
    p.write_text("outcome: status\nstrengths: [0.1]\n")
    assert load_config(p).strengths == (0.1,)
    t = tmp_path / "c.toml"
    t.write_text('outcome="status"\nstrengths=[0.1]\n')
    assert load_config(t) == load_config(p)


@pytest.mark.parametrize("orientation", ["samples_by_genes", "genes_by_samples", "auto"])
def test_bulk_axis_alignment(data, config, tmp_path, orientation):
    xp = tmp_path / "x.tsv"
    mp = tmp_path / "m.tsv"
    x = data.expression.iloc[::-1]
    if orientation == "genes_by_samples":
        x = x.T
    x.to_csv(xp, sep="\t", index_label="sample_id")
    data.metadata.to_csv(mp, sep="\t", index=False)
    loaded = load_bulk(xp, mp, config, orientation)
    pd.testing.assert_frame_equal(loaded.expression, data.expression, check_names=False)
    assert loaded.input_hashes


def test_duplicate_headers_rejected(tmp_path):
    p = tmp_path / "bad.csv"
    p.write_text("sample_id,g,g\ns,1,2\n")
    with pytest.raises(ValueError, match="Duplicate"):
        read_table(p)


@pytest.mark.parametrize(
    "kind", ["nan", "negative", "fractional", "zero_library", "metadata_feature", "duplicate_gene"]
)
def test_bad_expression(data, config, kind):
    if kind == "nan":
        data.expression.iloc[0, 0] = np.nan
    if kind == "negative":
        data.expression.iloc[0, 0] = -1
    if kind == "fractional":
        data.expression.iloc[0, 0] = 0.5
    if kind == "zero_library":
        data.expression.iloc[0] = 0
    if kind == "metadata_feature":
        data.expression["phenotype"] = data.metadata.phenotype
    if kind == "duplicate_gene":
        data.expression.columns = ["duplicate"] * data.expression.shape[1]
    with pytest.raises(ValueError):
        validate_dataset(data, config)


def test_missing_metadata(data, config):
    data.metadata.loc[data.metadata.index[0], "donor"] = None
    with pytest.raises(ValueError, match="Missing"):
        validate_dataset(data, config)


def test_positive_class_explicit(data, config):
    data.metadata["phenotype"] = np.where(data.metadata.phenotype == 1, "control", "case")
    y, labels = encode_outcome(data, replace(config, positive_class="case"))
    assert labels == ["control", "case"]
    assert np.array_equal(y, 1 - synthetic_labels(data))


def synthetic_labels(data):
    return (data.metadata.phenotype == "control").astype(int).to_numpy()


def test_no_human_specific_assumptions(data, config):
    data.expression.columns = [f"fungal_locus_{i}" for i in range(data.expression.shape[1])]
    assert validate_dataset(data, config)["n_genes"] == 60
