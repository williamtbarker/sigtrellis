from __future__ import annotations

import pytest

from sigtrellis.config import Config
from sigtrellis.domain import Dataset
from sigtrellis.simulate import synthetic_bulk


@pytest.fixture
def data() -> Dataset:
    return synthetic_bulk(samples=48, features=60, seed=11)[0]


@pytest.fixture
def config() -> Config:
    return Config(
        group="donor",
        outer_folds=3,
        inner_folds=2,
        strengths=(0.02, 0.15),
        l1_ratios=(0.5,),
        stability_resamples=4,
        permutations=0,
        max_features=50,
        min_groups_gate=20,
    )
