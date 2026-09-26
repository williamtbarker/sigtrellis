"""Orchestrate separate predictive, descriptive, and evidence-gating paths."""

from __future__ import annotations

import hashlib
import importlib.metadata
import platform
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from sigtrellis import __version__
from sigtrellis.config import Config
from sigtrellis.correlation import correlation_diagnostics
from sigtrellis.de import differential_expression
from sigtrellis.domain import Audit, Dataset, FloatArray, Split, array_hash, file_hash, write_json
from sigtrellis.metrics import group_weights, loss
from sigtrellis.modeling import Hyperparameters, fit_model, fit_prepared, tune
from sigtrellis.qc import encode_outcome, group_values, validate_dataset
from sigtrellis.splits import record_split
from sigtrellis.stability import stability_select
from sigtrellis.validation import nested_validate, permutation_control


def batch_diagnostics(
    data: Dataset, y: FloatArray, config: Config, audit: Audit
) -> tuple[dict[str, Any], FloatArray | None]:
    if not config.batch or not config.assess_batches:
        return {"status": "not_assessed", "folds": []}, None
    batches = data.metadata[config.batch].astype(str).to_numpy()
    groups = group_values(data, config)
    if len(set(batches)) < 2:
        return {"status": "not_assessed_single_batch", "folds": []}, None
    coefficients: list[FloatArray] = []
    rows: list[dict[str, Any]] = []
    for index, batch in enumerate(sorted(set(batches))):
        test = np.flatnonzero(batches == batch).astype(np.int64)
        train = np.flatnonzero(~np.isin(groups, groups[test])).astype(np.int64)
        try:
            split = Split(train, test, 0, index)
            record_split(audit, data, split, f"batch:{batch}", config)
            if config.outcome_type != "continuous" and set(y[train]) != set(y):
                raise ValueError("Training data lack an outcome class")
            tr, te = data.subset(train), data.subset(test)
            params, tuning = tune(tr, y[train], config, audit, f"batch:{batch}")
            model = fit_model(tr, y[train], config, params, audit, f"batch:{batch}/refit")
            coefficients.append(model.coefficients())
            baseline_params = Hyperparameters(config.strengths[0], config.l1_ratios[0])
            baseline_tuning = pd.DataFrame()
            if config.covariates:
                baseline_params, baseline_tuning = tune(
                    tr, y[train], config, audit, f"batch:{batch}/baseline", covariates_only=True
                )
            baseline = fit_prepared(
                model.prepared, tr, y[train], baseline_params, covariates_only=True
            )
            w = group_weights(group_values(te, config))
            improvement = loss(y[test], baseline.predict(te), config, w) - loss(
                y[test], model.predict(te), config, w
            )
            if not np.isfinite(improvement):
                raise ValueError("Batch holdout produced nonfinite loss improvement")
            rows.append(
                {
                    "batch": batch,
                    "status": "tested",
                    "loss_improvement": improvement,
                    "n_train_groups": len(set(groups[train])),
                    "n_test_groups": len(set(groups[test])),
                    "parameters": {"strength": params.strength, "l1_ratio": params.l1_ratio},
                    "tuning_results": tuning.to_dict("records"),
                    "baseline_parameters": {
                        "strength": baseline_params.strength,
                        "l1_ratio": baseline_params.l1_ratio,
                    },
                    "baseline_tuning_results": baseline_tuning.to_dict("records"),
                }
            )
        except ValueError as exc:
            rows.append({"batch": batch, "status": "unidentifiable", "reason": str(exc)})
    complete = all(r["status"] == "tested" for r in rows)
    return {"status": "tested" if complete else "incomplete", "folds": rows}, (
        np.asarray(coefficients, dtype=float) if complete else None
    )


def evidence_gates(
    table: pd.DataFrame,
    config: Config,
    qc: dict[str, Any],
    metrics: dict[str, Any],
    batch: dict[str, Any],
) -> tuple[pd.DataFrame, list[str]]:
    blockers = list(qc["gate_blockers"])
    permutation = metrics["permutation"]
    if permutation["status"] != "tested":
        blockers.append("no_eligible_permutation_test")
    if config.permutations < 19 or 1 / (1 + config.permutations) > config.permutation_alpha + 1e-12:
        blockers.append("insufficient_permutation_resolution")
    pvalue = permutation.get("pvalue")
    if pvalue is None or not np.isfinite(pvalue) or not 0 < pvalue <= 1:
        blockers.append("invalid_or_unavailable_permutation_pvalue")
    elif pvalue > config.permutation_alpha:
        blockers.append("permutation_control_not_passed")
    if not np.isfinite(metrics["loss_improvement"]):
        blockers.append("nonfinite_held_out_improvement")
    elif metrics["loss_improvement"] <= 0:
        blockers.append("no_held_out_improvement_over_baseline")
    if config.batch and batch["status"] != "tested":
        blockers.append("cross_batch_robustness_not_established")
    if batch["status"] == "tested" and any(
        not np.isfinite(f["loss_improvement"]) or f["loss_improvement"] <= 0 for f in batch["folds"]
    ):
        blockers.append("cross_batch_performance_failure")
    eligible = (
        (table.selection_frequency >= config.selection_threshold)
        & (table.sign_consistency >= config.sign_threshold)
        & (table.outer_selection_frequency >= config.outer_selection_threshold)
        & (table.outer_sign_consistency >= config.sign_threshold)
    )
    eligible &= table.coefficient.abs() > config.coefficient_tolerance
    eligible &= np.sign(table.coefficient) == np.sign(table.coefficient_median_selected)
    eligible &= np.sign(table.coefficient) == table.outer_dominant_sign
    for norm in config.stability_normalizations or (config.normalization,):
        eligible &= table[f"frequency_{norm}"] >= config.selection_threshold
    if batch["status"] == "tested":
        eligible &= table.batch_selection_frequency >= config.selection_threshold
        eligible &= table.batch_sign_consistency >= config.sign_threshold
        eligible &= np.sign(table.coefficient) == table.batch_dominant_sign
    result = table.copy()
    result["passes_feature_stability_gates"] = eligible
    result["passes_robustness_gates"] = eligible & (not blockers)
    result["evidence_status"] = np.where(
        result.passes_robustness_gates, "candidate_predictive_association", "exploratory_only"
    )
    result["biological_validation"] = "not_established"
    result["causal_claim"] = "not_established"
    return result, blockers


def _frequency_sign(coefficients: FloatArray, tolerance: float) -> tuple[FloatArray, FloatArray]:
    selected = abs(coefficients) > tolerance
    count = selected.sum(axis=0)
    sign = np.maximum(
        (coefficients > tolerance).sum(axis=0), (coefficients < -tolerance).sum(axis=0)
    )
    return selected.mean(axis=0), np.divide(
        sign, count, out=np.zeros_like(sign, dtype=float), where=count > 0
    )


def run_analysis(data: Dataset, config: Config, output: Path) -> dict[str, Any]:
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output directory is not empty; choose a fresh directory")
    output.mkdir(parents=True, exist_ok=True)
    audit = Audit()
    manifest: dict[str, Any] = {
        "status": "running",
        "package": "sigtrellis",
        "version": __version__,
        "config": config.to_dict(),
        "input_hashes": data.input_hashes,
        "expression_hash": array_hash(data.expression.to_numpy(dtype=float)),
        "sample_expression_hashes": {
            str(sample): array_hash(row)
            for sample, row in zip(
                data.expression.index, data.expression.to_numpy(dtype=float), strict=True
            )
        },
        "metadata_hash": hashlib.sha256(data.metadata.to_csv().encode()).hexdigest(),
        "python": sys.version,
        "platform": platform.platform(),
        "seed": config.seed,
        "thread_limit": 1,
        "cell_type": data.cell_type,
        "software_versions": {},
    }
    for name in (
        "numpy",
        "scipy",
        "pandas",
        "scikit-learn",
        "anndata",
        "h5py",
        "pydeseq2",
        "matplotlib",
    ):
        try:
            manifest["software_versions"][name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            manifest["software_versions"][name] = "not_installed"
    manifest["source_hashes"] = {
        p.name: file_hash(p) for p in sorted(Path(__file__).parent.glob("*.py"))
    }
    write_json(output / "run_manifest.json", manifest)
    write_json(output / "configuration.json", config.to_dict())
    try:
        with threadpool_limits(limits=1):
            result = _run(data, config, output, audit, manifest)
    except Exception as exc:
        manifest.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        write_json(output / "run_manifest.json", manifest)
        write_json(
            output / "audit.json",
            {"fits": audit.fits, "splits": audit.splits, "warnings": audit.warnings},
        )
        raise
    return result


def _run(
    data: Dataset, config: Config, output: Path, audit: Audit, manifest: dict[str, Any]
) -> dict[str, Any]:
    from sigtrellis.reporting import generate_report

    qc = validate_dataset(data, config)
    write_json(output / "qc_report.json", qc)
    y, classes = encode_outcome(data, config)
    contrasts = (
        ["response_per_training_sd"]
        if config.outcome_type == "continuous"
        else [f"{c}_vs_{classes[0]}" for c in classes[1:]]
    )
    validation = nested_validate(data, y, config, audit)
    metrics = validation.metrics
    metrics["classes"] = classes
    metrics["permutation"] = permutation_control(
        data, y, config, metrics["loss_improvement"], audit
    )
    stability = stability_select(data, y, config, contrasts, audit)
    params, final_tuning = tune(data, y, config, audit, "final_refit")
    final = fit_model(data, y, config, params, audit, "final_refit/all_samples")
    final_coefficients = final.coefficients()
    batch, batch_coefficients = batch_diagnostics(data, y, config, audit)
    if any("SOLVER_NONCONVERGENCE" in message for message in audit.warnings):
        qc["gate_blockers"].append("solver_nonconvergence_in_search")
        qc["warnings"].append("SOLVER_REVIEW: one or more tuning candidates did not converge")
        write_json(output / "qc_report.json", qc)
    metrics["batch_validation"] = batch
    outer_freq, outer_sign = _frequency_sign(validation.coefficients, config.coefficient_tolerance)
    table = stability.table.copy()
    table["coefficient"] = final_coefficients.ravel()
    table["coefficient_direction"] = np.where(
        table.coefficient > config.coefficient_tolerance,
        "positive",
        np.where(table.coefficient < -config.coefficient_tolerance, "negative", "zero"),
    )
    table["outer_selection_frequency"] = outer_freq.ravel()
    table["outer_sign_consistency"] = outer_sign.ravel()
    table["outer_dominant_sign"] = np.sign(
        np.sign(
            np.where(
                abs(validation.coefficients) > config.coefficient_tolerance,
                validation.coefficients,
                0,
            )
        ).sum(axis=0)
    ).ravel()
    table["batch_selection_frequency"] = np.nan
    table["batch_sign_consistency"] = np.nan
    table["batch_dominant_sign"] = np.nan
    if batch_coefficients is not None:
        batch_freq, batch_sign = _frequency_sign(batch_coefficients, config.coefficient_tolerance)
        table["batch_selection_frequency"] = batch_freq.ravel()
        table["batch_sign_consistency"] = batch_sign.ravel()
        table["batch_dominant_sign"] = np.sign(
            np.sign(
                np.where(
                    abs(batch_coefficients) > config.coefficient_tolerance, batch_coefficients, 0
                )
            ).sum(axis=0)
        ).ravel()
    table["gene_symbol"] = table.gene_id.map(data.symbols).fillna("")
    table["expression_prevalence"] = (
        table.gene_id.map((data.expression > 0).mean(axis=0))
        if config.input_scale == "counts"
        else np.nan
    )
    table["cell_type"] = data.cell_type or "bulk"
    table["phenotype_association"] = config.outcome
    de: pd.DataFrame | None = None
    table["de_log2_fold_change"] = np.nan
    table["de_adjusted_pvalue"] = np.nan
    if config.supporting_de:
        try:
            de = differential_expression(data, y, config)
            de.to_csv(output / "differential_expression.csv")
            for column in ("de_log2_fold_change", "de_adjusted_pvalue"):
                table[column] = table.gene_id.map(de[column])
            audit.warnings.extend(f"EXPLORATORY_DE: {w}" for w in de.attrs["warnings"])
            write_json(output / "de_design.json", de.attrs)
        except ValueError as exc:
            audit.warnings.append(f"EXPLORATORY_DE_UNAVAILABLE: {exc}")
    correlation = correlation_diagnostics(
        data, stability.coefficients, table, config, contrasts=contrasts
    )
    table["correlated_feature_group"] = table.gene_id.map(correlation.membership).fillna(
        "not_assessed"
    )
    frequencies = correlation.groups.set_index(["group_id", "contrast"])["selection_frequency"]
    table["correlated_group_frequency"] = [
        frequencies.get((group, contrast), np.nan)
        for group, contrast in zip(table.correlated_feature_group, table.contrast, strict=True)
    ]
    table, blockers = evidence_gates(table, config, qc, metrics, batch)
    table = table.sort_values(
        [
            "passes_robustness_gates",
            "selection_frequency",
            "sign_consistency",
            "outer_selection_frequency",
            "rank_median_selected",
            "gene_id",
            "contrast",
        ],
        ascending=[False, False, False, False, True, True, True],
        kind="stable",
    )
    table.insert(0, "rank", np.arange(1, len(table) + 1))
    table.to_csv(output / "biomarkers.csv", index=False)
    stability.table.to_csv(output / "biomarker_stability.csv", index=False)
    table[["gene_id", "contrast", "coefficient", "coefficient_direction"]].to_csv(
        output / "coefficients.csv", index=False
    )
    pd.concat([validation.tuning, final_tuning, stability.tuning], ignore_index=True).to_csv(
        output / "cv_results.csv", index=False
    )
    validation.predictions.to_csv(output / "out_of_fold_predictions.csv", index=False)
    correlation.groups.to_csv(output / "correlation_groups.csv", index=False)
    correlation.substitutions.to_csv(output / "feature_substitution.csv", index=False)
    correlation.matrix.to_csv(output / "correlations.csv")
    write_json(output / "stability_summary.json", stability.summary)
    write_json(output / "stability_resamples.json", stability.resamples)
    np.savez_compressed(
        output / "resample_coefficients.npz",
        stability=stability.coefficients,
        outer=validation.coefficients,
        batch=batch_coefficients
        if batch_coefficients is not None
        else np.empty((0, len(contrasts), data.expression.shape[1])),
        batch_labels=np.asarray(
            [row["batch"] for row in batch["folds"]] if batch_coefficients is not None else [],
            dtype=str,
        ),
        genes=np.asarray(data.expression.columns, dtype=str),
        contrasts=np.asarray(contrasts, dtype=str),
    )
    # Coefficient path: exploratory full-data fits after all predictive evaluation.
    path_rows: list[dict[str, Any]] = []
    for strength in config.strengths:
        path_model = fit_prepared(
            final.prepared, data, y, Hyperparameters(strength, params.l1_ratio)
        )
        for ci, contrast in enumerate(contrasts):
            for gene, coefficient in zip(
                data.expression.columns, path_model.coefficients()[ci], strict=True
            ):
                path_rows.append(
                    {
                        "strength": strength,
                        "l1_ratio": params.l1_ratio,
                        "contrast": contrast,
                        "gene_id": gene,
                        "coefficient": coefficient,
                    }
                )
    path_table = pd.DataFrame(path_rows)
    path_table.to_csv(output / "coefficient_path.csv", index=False)
    write_json(output / "model_state.json", final.to_state())
    data.metadata.to_csv(output / "sample_metadata.csv", index=False)
    if data.cell_type:
        data.expression.to_csv(
            output / "pseudobulk_counts.tsv.gz",
            sep="\t",
            compression="gzip",
            index_label=config.sample_id,
        )
    # The permutation audit is large and separated from the reader-facing metrics.
    perm_audit = metrics["permutation"].pop("audits", [])
    write_json(output / "permutation_audit.json", perm_audit)
    write_json(output / "model_metrics.json", metrics)
    write_json(
        output / "audit.json",
        {"fits": audit.fits, "splits": audit.splits, "warnings": audit.warnings},
    )
    manifest.update(
        status="complete",
        n_samples=len(data.expression),
        n_features=data.expression.shape[1],
        n_biological_groups=qc["n_biological_groups"],
        classes=classes,
        final_hyperparameters={"strength": params.strength, "l1_ratio": params.l1_ratio},
        selected_features=table.loc[
            table.coefficient.abs() > config.coefficient_tolerance, ["gene_id", "contrast"]
        ].to_dict("records"),
        candidate_features=table.loc[
            table.passes_robustness_gates, ["gene_id", "contrast"]
        ].to_dict("records"),
        gate_blockers=blockers,
        warnings=sorted(set(qc["warnings"] + audit.warnings)),
    )
    generate_report(
        data,
        y,
        config,
        output,
        qc,
        metrics,
        table,
        correlation,
        validation,
        path_table,
        stability.summary,
        manifest,
    )
    manifest["output_hashes"] = {
        str(p.relative_to(output)): file_hash(p)
        for p in sorted(output.rglob("*"))
        if p.is_file() and p.name != "run_manifest.json"
    }
    write_json(output / "run_manifest.json", manifest)
    return manifest
