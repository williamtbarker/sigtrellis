"""Known-truth test fixtures; these are not public biological demonstrations."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse

from sigtrellis.domain import Dataset, write_json


def synthetic_bulk(
    samples: int = 80, features: int = 200, seed: int = 7, scenario: str = "signal"
) -> tuple[Dataset, list[str]]:
    if samples < 12 or features < 10:
        raise ValueError("Synthetic examples require >=12 samples and >=10 genes")
    rng = np.random.default_rng(seed)
    y = np.tile([0, 1], int(np.ceil(samples / 2)))[:samples]
    rng.shuffle(y)
    batch = np.array([f"batch_{i % 3}" for i in range(samples)])
    outcome: np.ndarray = y.copy()
    if scenario == "imbalance":
        outcome = np.zeros(samples, dtype=int)
        outcome[rng.choice(samples, max(6, samples // 5), replace=False)] = 1
        y = outcome.copy()
    if scenario == "multiclass":
        outcome = np.arange(samples) % 3
        rng.shuffle(outcome)
    if scenario == "continuous":
        outcome = rng.normal(size=samples)
    base = rng.uniform(20, 100, size=features)
    log_mu = np.tile(np.log(base), (samples, 1))
    truth = [f"gene_{i:04d}" for i in range(4)]
    signal = np.asarray(outcome, dtype=float)
    if scenario == "noise":
        truth = []
    elif scenario == "batch_confounded":
        batch = np.where(y == 1, "batch_case", "batch_control")
        log_mu[:, :20] += y[:, None] * 1.6
        truth = []
    elif scenario == "batch_specific":
        log_mu[:, :4] += (signal * (batch == "batch_0"))[:, None] * 2
    elif scenario == "multiclass":
        for cls in range(3):
            log_mu[:, cls * 2 : cls * 2 + 2] += (outcome == cls)[:, None] * 1.6
        truth = [f"gene_{i:04d}" for i in range(6)]
    else:
        log_mu[:, :2] += signal[:, None] * 1.7
        log_mu[:, 2:4] -= signal[:, None] * 1.7
    if scenario == "correlated":
        latent = signal * 1.7 + rng.normal(0, 0.25, size=samples)
        log_mu[:, :6] = np.log(50) + latent[:, None] + rng.normal(0, 0.04, size=(samples, 6))
        truth = [f"gene_{i:04d}" for i in range(6)]
    log_mu += rng.normal(0, 0.2, size=samples)[:, None]
    mu = np.exp(log_mu)
    x = rng.negative_binomial(20, 20 / (20 + mu)).astype(float)
    if scenario == "outliers":
        x[0, :] *= 1000
    ids = [f"sample_{i:04d}" for i in range(samples)]
    genes = [f"gene_{i:04d}" for i in range(features)]
    expression = pd.DataFrame(x, index=ids, columns=genes)
    metadata = pd.DataFrame(
        {"sample_id": ids, "donor": ids, "phenotype": outcome, "batch": batch}, index=ids
    )
    return Dataset(expression, metadata), truth


def synthetic_single_cell(
    path: Path, donors: int = 24, genes: int = 100, cells: int = 30, seed: int = 7
) -> None:
    rng = np.random.default_rng(seed)
    counts: list[np.ndarray] = []
    records: list[dict[str, str]] = []
    for donor in range(donors):
        donor_effect = rng.normal(0, 0.35, size=genes)
        for condition in (0, 1):
            for cell_type in ("type_A", "type_B"):
                log_mu = np.log(2) + donor_effect
                log_mu = log_mu.copy()
                log_mu[:4] += condition * 1.7
                mu = np.exp(log_mu)
                values = rng.negative_binomial(5, 5 / (5 + mu), size=(cells, genes))
                counts.append(values)
                for _ in range(cells):
                    records.append(
                        {
                            "sample_id": f"d{donor:03d}_s{condition}",
                            "donor": f"d{donor:03d}",
                            "phenotype": str(condition),
                            "cell_type": cell_type,
                            "batch": f"batch_{donor % 3}",
                        }
                    )
    obs = pd.DataFrame(records, index=[f"cell_{i:06d}" for i in range(len(records))])
    var = pd.DataFrame(index=[f"gene_{i:04d}" for i in range(genes)])
    ad.AnnData(sparse.csr_matrix(np.vstack(counts)), obs=obs, var=var).write_h5ad(
        path, compression="gzip"
    )


def write_example(
    output: Path, modality: str, samples: int, features: int, seed: int, scenario: str
) -> None:
    if output.exists() and any(output.iterdir()):
        raise ValueError("Synthetic output directory must be empty")
    output.mkdir(parents=True, exist_ok=True)
    config = {
        "outcome": "phenotype",
        "sample_id": "sample_id",
        "group": "donor",
        "outer_folds": 3,
        "inner_folds": 2,
        "stability_resamples": 20,
        "strengths": [0.01, 0.05, 0.2],
        "l1_ratios": [0.5, 0.9],
        "permutations": 19,
        "seed": seed,
        "max_features": min(features, 1000),
    }
    if modality == "bulk":
        data, truth = synthetic_bulk(samples, features, seed, scenario)
        data.expression.to_csv(output / "counts.tsv", sep="\t", index_label="sample_id")
        data.metadata.to_csv(output / "metadata.tsv", sep="\t", index=False)
        config["outcome_type"] = scenario if scenario in {"continuous", "multiclass"} else "binary"
        if scenario in {"batch_confounded", "batch_specific"}:
            config["batch"] = "batch"
    else:
        synthetic_single_cell(output / "experiment.h5ad", samples, features, seed=seed)
        truth = [f"gene_{i:04d}" for i in range(4)]
        config.update(
            cell_type="cell_type", cell_type_value="type_A", permutation_scheme="within_group"
        )
    (output / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=True))
    write_json(
        output / "truth.json", {"synthetic": True, "scenario": scenario, "planted_genes": truth}
    )
