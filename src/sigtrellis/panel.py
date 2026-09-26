"""Select and fit an explicitly compact panel entirely inside a training boundary."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, replace

import pandas as pd

from sigtrellis.config import Config
from sigtrellis.domain import Audit, Dataset, FloatArray
from sigtrellis.modeling import FittedModel, Hyperparameters, fit_model
from sigtrellis.qc import encode_outcome
from sigtrellis.stability import StabilityResult, stability_select


@dataclass
class PanelResult:
    model: FittedModel
    features: list[str]
    evidence: pd.DataFrame
    stability: StabilityResult


def fit_panel(
    data: Dataset,
    y: FloatArray,
    config: Config,
    audit: Audit,
    context: str,
    stability: StabilityResult | None = None,
) -> PanelResult:
    """Fixed panel policy: stable features, predeclared cap, modal tuning choice.

    Hyperparameters come from the training subsample searches. We do not tune
    again on an inner split after using that split's labels to select the panel.
    An empty selected panel becomes an explicit covariate/intercept model.
    """
    if stability is None:
        _, classes = encode_outcome(data, config)
        contrasts = (
            [f"{c}_vs_{classes[0]}" for c in classes[1:]]
            if classes
            else ["response_per_training_sd"]
        )
        stability = stability_select(data, y, config, contrasts, audit, context + "/selection")
    table = stability.table
    eligible = table[
        (table.selection_frequency >= config.selection_threshold)
        & (table.sign_consistency >= config.sign_threshold)
    ]
    ranked = eligible.sort_values(
        ["selection_frequency", "sign_consistency", "rank_median_selected", "gene_id"],
        ascending=[False, False, True, True],
        kind="stable",
    )
    selected = ranked.gene_id.drop_duplicates().head(config.panel_max_features).tolist()
    counts = Counter((float(r["strength"]), float(r["l1_ratio"])) for r in stability.resamples)
    strength, ratio = sorted(counts, key=lambda key: (-counts[key], -key[0], -key[1]))[0]
    local = replace(config, candidate_method="none")
    model = fit_model(
        data,
        y,
        local,
        Hyperparameters(strength, ratio),
        audit,
        context + "/refit",
        allowed_features=set(selected),
    )
    actual = model.prepared.genes
    evidence = table[table.gene_id.isin(actual)].copy()
    coefficients = model.coefficients()
    contrast_ids = {str(value): i for i, value in enumerate(table.contrast.drop_duplicates())}
    gene_ids = {value: i for i, value in enumerate(model.prepared.all_genes)}
    evidence["panel_coefficient"] = [
        float(coefficients[contrast_ids[str(contrast)], gene_ids[str(gene)]])
        for gene, contrast in zip(evidence.gene_id, evidence.contrast, strict=True)
    ]
    return PanelResult(model, actual, evidence, stability)
