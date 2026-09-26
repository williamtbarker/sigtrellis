"""Biological-group split contracts, including purged batch holdouts."""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np
from sklearn.model_selection import GroupKFold, KFold, StratifiedKFold

from sigtrellis.config import Config
from sigtrellis.domain import Audit, Dataset, FloatArray, IntArray, Split
from sigtrellis.qc import group_values


def assert_disjoint(data: Dataset, train: IntArray, test: IntArray, config: Config) -> None:
    groups = group_values(data, config)
    if np.intersect1d(train, test).size or set(groups[train]) & set(groups[test]):
        raise ValueError("Biological groups cross the train/test boundary")
    if not len(train) or not len(test):
        raise ValueError("Empty train/test split after group purging")


def make_splits(
    data: Dataset, y: FloatArray, config: Config, *, inner: bool = False, seed: int | None = None
) -> list[Split]:
    groups = group_values(data, config)
    unique = np.unique(groups)
    n_folds = config.inner_folds if inner else config.outer_folds
    strategy = "grouped" if inner else config.cv_strategy
    repeats = 1 if inner else config.repeats
    result: list[Split] = []
    random_seed = config.seed if seed is None else seed
    if strategy == "grouped" and len(unique) < n_folds:
        raise ValueError(f"{len(unique)} biological groups cannot support {n_folds} folds")
    homogeneous = all(len(np.unique(y[groups == group])) == 1 for group in unique)
    for repeat in range(repeats):
        pairs: list[tuple[Any, Any]] = []
        if strategy == "leave_batch_out":
            assert config.batch
            batches = data.metadata[config.batch].astype(str).to_numpy()
            if len(np.unique(batches)) < 2:
                raise ValueError("Batch holdout needs >=2 batches")
            for batch in np.unique(batches):
                test = np.flatnonzero(batches == batch)
                train = np.flatnonzero(~np.isin(groups, groups[test]))
                pairs.append((train, test))
        elif strategy == "leave_group_out":
            pairs = [(np.flatnonzero(groups != g), np.flatnonzero(groups == g)) for g in unique]
        elif homogeneous:
            group_y = np.array([y[groups == g][0] for g in unique])
            if config.outcome_type != "continuous":
                if np.bincount(group_y.astype(int)).min() < n_folds:
                    raise ValueError(
                        "Too few biological groups per class for requested stratified folds"
                    )
                splitter: Any = StratifiedKFold(
                    n_folds, shuffle=True, random_state=random_seed + repeat
                )
            else:
                splitter = KFold(n_folds, shuffle=True, random_state=random_seed + repeat)
            for train_g, test_g in splitter.split(unique, group_y):
                pairs.append(
                    (
                        np.flatnonzero(np.isin(groups, unique[train_g])),
                        np.flatnonzero(np.isin(groups, unique[test_g])),
                    )
                )
        else:
            splitter = GroupKFold(n_folds, shuffle=True, random_state=random_seed + repeat)
            pairs = list(splitter.split(np.zeros(len(y)), y, groups))
        for fold, (tr, te) in enumerate(pairs):
            train, test = np.asarray(tr, dtype=np.int64), np.asarray(te, dtype=np.int64)
            assert_disjoint(data, train, test, config)
            if config.outcome_type != "continuous" and set(y[train]) != set(y):
                raise ValueError(
                    "Training fold lacks an outcome class; redesign CV or collect replicates"
                )
            result.append(Split(train, test, repeat, fold))
    return result


def record_split(audit: Audit, data: Dataset, split: Split, context: str, config: Config) -> None:
    assert_disjoint(data, split.train, split.test, config)
    audit.splits.append(
        {
            "context": context,
            "repeat": split.repeat,
            "fold": split.fold,
            "train_samples": data.expression.index[split.train].tolist(),
            "test_samples": data.expression.index[split.test].tolist(),
            "train_groups": sorted(set(group_values(data, config)[split.train])),
            "test_groups": sorted(set(group_values(data, config)[split.test])),
        }
    )


def subsample_groups(
    data: Dataset, y: FloatArray, config: Config, rng: np.random.Generator
) -> IntArray:
    groups = group_values(data, config)
    unique = np.unique(groups)
    n = min(
        len(unique) - 1,
        max(config.inner_folds, int(np.ceil(len(unique) * config.stability_fraction))),
    )
    # Rejection sampling preserves complete biological groups and avoids impossible inner CV.
    for _ in range(200):
        chosen = rng.choice(unique, n, replace=False)
        rows = np.flatnonzero(np.isin(groups, chosen)).astype(np.int64)
        try:
            make_splits(data.subset(rows), y[rows], replace(config, repeats=1), inner=True)
            if config.outcome_type == "continuous" or set(y[rows]) == set(y):
                return rows
        except ValueError:
            continue
    raise ValueError("Stability subsamples cannot support inner CV; increase fraction/groups")


def permute_outcome(
    data: Dataset, y: FloatArray, config: Config, rng: np.random.Generator
) -> FloatArray:
    groups = group_values(data, config)
    result = y.copy()
    unique = np.unique(groups)
    if config.permutation_scheme == "within_group":
        if config.permutation_strata:
            raise ValueError("within_group and permutation_strata cannot be combined")
        if all(len(np.unique(y[groups == g])) == 1 for g in unique):
            raise ValueError("Within-group permutation cannot alter a group-constant outcome")
        for g in unique:
            rows = np.flatnonzero(groups == g)
            result[rows] = rng.permutation(y[rows])
    else:
        if any(len(np.unique(y[groups == g])) != 1 for g in unique):
            raise ValueError(
                "Group permutation requires constant outcome per group; declare within_group if exchangeable"
            )
        strata: list[str] = []
        for g in unique:
            if config.permutation_strata:
                s = data.metadata.loc[groups == g, config.permutation_strata].astype(str).unique()
                if len(s) != 1:
                    raise ValueError(
                        "Permutation stratum must be constant within each biological group"
                    )
                strata.append(str(s[0]))
            else:
                strata.append("all")
        group_y = np.array([y[groups == g][0] for g in unique])
        shuffled = group_y.copy()
        for s in sorted(set(strata)):
            rows = np.flatnonzero(np.asarray(strata) == s)
            shuffled[rows] = rng.permutation(group_y[rows])
        for g, value in zip(unique, shuffled, strict=True):
            result[groups == g] = value
    return result
