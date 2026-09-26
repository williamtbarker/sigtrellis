"""Self-contained scientific reports and purpose-labeled diagnostic figures."""

from __future__ import annotations

import base64
import html
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.metrics import precision_recall_curve, roc_curve

from sigtrellis.cell_reporting import cell_evidence
from sigtrellis.config import Config
from sigtrellis.correlation import CorrelationResult
from sigtrellis.domain import Dataset, FloatArray
from sigtrellis.metrics import group_weights
from sigtrellis.preprocessing import descriptive_normalized
from sigtrellis.validation import ValidationResult


def _save(output: Path, name: str, caption: str, figures: list[tuple[str, str]]) -> None:
    plt.tight_layout()
    plt.savefig(
        output / "figures" / name, dpi=140, bbox_inches="tight", metadata={"Software": "SigTrellis"}
    )
    plt.close()
    figures.append((name, caption))


def make_figures(
    data: Dataset,
    y: FloatArray,
    config: Config,
    output: Path,
    metrics: dict[str, Any],
    table: pd.DataFrame,
    correlation: CorrelationResult,
    validation: ValidationResult,
    path: pd.DataFrame,
) -> list[tuple[str, str]]:
    (output / "figures").mkdir()
    figures: list[tuple[str, str]] = []
    normalized = descriptive_normalized(data, config)
    variable = normalized.var(axis=0) > 1e-12
    if variable.sum() >= 2:
        pc = PCA(n_components=2, svd_solver="full").fit(normalized[:, variable])
        scores = pc.transform(normalized[:, variable])
        plt.figure(figsize=(7, 5))
        if config.outcome_type == "continuous":
            plt.scatter(scores[:, 0], scores[:, 1], c=y, cmap="viridis")
            plt.colorbar(label=config.outcome)
        else:
            for i, label in enumerate(metrics["classes"]):
                mask = y == i
                plt.scatter(scores[mask, 0], scores[mask, 1], label=label)
            plt.legend()
        plt.xlabel(f"PC1 ({pc.explained_variance_ratio_[0]:.1%})")
        plt.ylabel(f"PC2 ({pc.explained_variance_ratio_[1]:.1%})")
        plt.title("Sample structure • phenotype")
        _save(
            output,
            "pca_phenotype.png",
            "Full-cohort normalized PCA for exploration only; never supplied to predictive fits.",
            figures,
        )
        if config.batch:
            plt.figure(figsize=(7, 5))
            for label in sorted(data.metadata[config.batch].astype(str).unique()):
                mask = data.metadata[config.batch].astype(str).to_numpy() == label
                plt.scatter(scores[mask, 0], scores[mask, 1], label=label)
            plt.legend()
            plt.xlabel("PC1")
            plt.ylabel("PC2")
            plt.title("Sample structure • batch")
            _save(
                output,
                "pca_batch.png",
                "The same descriptive PCA colored by batch; visual separation is not proof of confounding.",
                figures,
            )
    top = table.head(20).iloc[::-1]
    labels = (top.feature_label if "feature_label" in top else top.gene_id).astype(str) + (
        " | " + top.contrast.astype(str) if config.outcome_type == "multiclass" else ""
    )
    plt.figure(figsize=(8, max(4, len(top) * 0.25)))
    plt.barh(labels, top.selection_frequency, color="#327c8c")
    plt.axvline(config.selection_threshold, color="#b04728", linestyle="--")
    plt.xlabel("Selection frequency across biological-group subsamples")
    plt.xlim(0, 1)
    plt.title("Empirical feature stability")
    _save(
        output,
        "selection_frequency.png",
        "Overlapping subsamples; frequencies are not independent-trial probabilities or adjusted p-values.",
        figures,
    )
    plt.figure(figsize=(8, max(4, len(top) * 0.25)))
    for i, row in enumerate(
        top[["coefficient_q25", "coefficient_q75", "coefficient_median"]].to_numpy(dtype=float)
    ):
        plt.plot(row[:2], [i, i], color="#327c8c", linewidth=3)
        plt.scatter(row[2], i, color="#b04728", s=18)
    plt.yticks(range(len(top)), labels.tolist())
    plt.axvline(0, color="grey", linewidth=0.8)
    plt.xlabel("Coefficient per training-fold SD of the feature")
    plt.title("Median and interquartile range across subsamples")
    _save(
        output,
        "coefficient_stability.png",
        "IQRs describe coefficient variation including zero fits; they are not confidence intervals.",
        figures,
    )
    candidate_genes = table.gene_id.drop_duplicates().head(8)
    first_contrast = table.contrast.iloc[0]
    plt.figure(figsize=(8, 5))
    for gene in candidate_genes:
        sub = path[(path.gene_id == gene) & (path.contrast == first_contrast)].sort_values(
            "strength"
        )
        label = data.features[gene].label if gene in data.features else gene
        plt.plot(sub.strength, sub.coefficient, marker="o", label=label)
    plt.xscale("log")
    plt.xlabel("Lambda on mean loss")
    plt.ylabel("Standardized coefficient")
    plt.title(f"Exploratory regularization path • {first_contrast}")
    plt.legend(fontsize=7, ncol=2)
    _save(
        output,
        "coefficient_path.png",
        "Full-data diagnostic after validation; mixing ratio fixed at the final chosen value.",
        figures,
    )
    tuning = validation.tuning
    tuning = tuning[~tuning.context.str.contains("covariate_baseline")]
    surface = tuning.groupby(["strength", "l1_ratio"]).loss.mean().unstack()
    plt.figure(figsize=(6, 4))
    plt.imshow(surface, aspect="auto", cmap="viridis_r")
    plt.xticks(range(len(surface.columns)), [str(v) for v in surface.columns])
    plt.yticks(range(len(surface.index)), [str(v) for v in surface.index])
    plt.xlabel("L1 mixing ratio")
    plt.ylabel("Lambda")
    plt.colorbar(label="Mean inner held-out loss")
    plt.title("Tuning landscape across outer training folds")
    _save(
        output,
        "tuning_surface.png",
        "Inner-fold tuning losses, not an independent performance estimate.",
        figures,
    )
    if not correlation.matrix.empty:
        display = correlation.matrix.iloc[:30, :30]
        plt.figure(figsize=(9, 8))
        plt.imshow(display, vmin=-1, vmax=1, cmap="RdBu_r")
        plt.xticks(range(len(display)), display.columns.tolist(), rotation=90, fontsize=6)
        plt.yticks(range(len(display)), display.index.tolist(), fontsize=6)
        plt.colorbar(label="Pearson r")
        plt.title("Descriptive marginal feature correlations")
        _save(
            output,
            "feature_correlation.png",
            "Up to 30 candidate features. Phenotype itself can induce marginal correlation; clusters are not established pathways.",
            figures,
        )
    de_path = output / "differential_expression.csv"
    de = pd.read_csv(de_path) if de_path.exists() else pd.DataFrame()
    de_groups: list[tuple[Any, pd.DataFrame]] = (
        list(de.groupby(["cell_type", "contrast"], sort=True)) if len(de) else []
    )
    for di, ((state, contrast), contrast_table) in enumerate(de_groups):
        shown = contrast_table.dropna(subset=["de_log2_fold_change", "de_adjusted_pvalue"])
        if shown.empty:
            continue
        plt.figure(figsize=(7, 5))
        plt.scatter(
            shown.de_log2_fold_change,
            -np.log10(np.maximum(shown.de_adjusted_pvalue, 1e-300)),
            s=6,
            alpha=0.5,
        )
        plt.xlabel(
            "Log2 expression change per outcome unit"
            if config.outcome_type == "continuous"
            else "Count-model log2 fold change"
        )
        plt.ylabel("−log10 adjusted p-value")
        plt.title(f"Exploratory differential expression • {state} • {contrast}")
        _save(
            output,
            "de_volcano.png" if di == 0 else f"de_volcano_{di}.png",
            "Same-cohort supporting evidence, not independent biological validation; unshrunk effects.",
            figures,
        )
    frame = validation.predictions
    observed = frame.observed.to_numpy(dtype=float)
    weights = group_weights(frame.group_id.to_numpy())
    if config.outcome_type != "continuous":
        class_ids = [1] if config.outcome_type == "binary" else list(range(len(metrics["classes"])))
        plt.figure(figsize=(6, 5))
        for i in class_ids:
            fpr, tpr, _ = roc_curve(
                (observed == i).astype(int), frame[f"p_{i}"], sample_weight=weights
            )
            plt.plot(fpr, tpr, label=metrics["classes"][i])
        plt.plot([0, 1], [0, 1], "--", color="grey")
        plt.xlabel("False positive rate")
        plt.ylabel("True positive rate")
        plt.legend()
        plt.title("Nested out-of-fold ROC")
        _save(
            output,
            "roc.png",
            "Biological-group-weighted held-out predictions; repeated folds do not create independent donors.",
            figures,
        )
        plt.figure(figsize=(6, 5))
        for i in class_ids:
            precision, recall, _ = precision_recall_curve(
                (observed == i).astype(int), frame[f"p_{i}"], sample_weight=weights
            )
            plt.plot(recall, precision, label=metrics["classes"][i])
            plt.axhline(
                np.average(observed == i, weights=weights), color="grey", linestyle=":", alpha=0.6
            )
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.legend()
        plt.title("Nested out-of-fold precision–recall")
        _save(
            output,
            "precision_recall.png",
            "Dotted lines show class prevalence; numerical PR-AUC uses average precision.",
            figures,
        )
        plt.figure(figsize=(6, 5))
        for i in class_ids:
            bins = (
                metrics["out_of_fold"]["calibration"]
                if config.outcome_type == "binary"
                else metrics["out_of_fold"]["calibration_by_class"][str(i)]
            )
            plt.plot(
                [b["mean_probability"] for b in bins],
                [b["observed_fraction"] for b in bins],
                marker="o",
                label=metrics["classes"][i],
            )
        plt.plot([0, 1], [0, 1], "--", color="grey")
        plt.xlabel("Predicted probability")
        plt.ylabel("Observed fraction")
        plt.legend()
        plt.title("Held-out calibration • 5 fixed bins")
        _save(
            output,
            "calibration.png",
            "Descriptive reliability only; no calibrator is fitted using outer test outcomes.",
            figures,
        )
    else:
        plt.figure(figsize=(6, 5))
        plt.scatter(observed, frame.prediction, s=15)
        low, high = (
            min(observed.min(), frame.prediction.min()),
            max(observed.max(), frame.prediction.max()),
        )
        plt.plot([low, high], [low, high], "--", color="grey")
        plt.xlabel("Observed outcome")
        plt.ylabel("Held-out prediction")
        _save(
            output,
            "regression.png",
            "Nested held-out regression predictions in original outcome units.",
            figures,
        )
    null = metrics["permutation"]["null_improvements"]
    if null:
        plt.figure(figsize=(7, 4))
        plt.hist(null, bins=min(15, len(null)), color="#327c8c")
        plt.axvline(metrics["loss_improvement"], color="#b04728", label="Observed")
        plt.xlabel("Held-out loss improvement over baseline")
        plt.ylabel("Permutation count")
        plt.legend()
        _save(
            output,
            "permutation_control.png",
            "The full nested procedure is repeated after exchangeability-aware label permutation.",
            figures,
        )
    if cell_evidence(data, table, config, output):
        figures.append(
            (
                "cell_distributions.png",
                "Post-selection descriptive evidence in the original feature units. Each point is a specimen, not a cell; related specimens remain dependent. No post-selection p-values are inferred. See cell_distribution_evidence.csv for group IDs and state cell counts.",
            )
        )
    return figures


def generate_report(
    data: Dataset,
    y: FloatArray,
    config: Config,
    output: Path,
    qc: dict[str, Any],
    metrics: dict[str, Any],
    table: pd.DataFrame,
    correlation: CorrelationResult,
    validation: ValidationResult,
    path: pd.DataFrame,
    stability: dict[str, Any],
    manifest: dict[str, Any],
) -> None:
    figures = make_figures(data, y, config, output, metrics, table, correlation, validation, path)
    count = int(table.loc[table.passes_robustness_gates, "gene_id"].nunique())
    noun = "features" if data.features else "genes"
    title = f"SigTrellis • {data.cell_type or 'bulk RNA-seq'}"
    summary = (
        f"{len(data.expression)} samples from {qc['n_biological_groups']} biological groups; "
        f"{data.expression.shape[1]} input {noun}. {count} {noun} pass the predefined candidate-association gates."
    )
    limits = [
        "Statistical selection is predictive association. Biological validation, mechanism, causality, clinical utility, and biomarker qualification are not established.",
        "Nested metrics evaluate the training procedure. The final fitted model and post-hoc consensus panel have not been independently validated.",
        "Full-cohort DE, PCA, correlations, and coefficient paths are exploratory and never supplied to CV fitting. Same-cohort DE is not independent confirmation.",
        "Selection frequencies and coefficient/rank IQRs describe overlapping fits. They do not provide formal false-discovery or PFER control.",
        "Batch omission or unknown sample relationships can invalidate validation. More cells do not replace more donors. A perfect phenotype/batch confound is unidentifiable.",
        "Log-relative expression is compositional; apparent decreases may reflect other genes increasing. Correlation groups need biological annotation and independent replication.",
        "The permutation scheme requires exchangeability justified by study design. Covariate-adjusted conditional significance is not implemented; with covariates, gates remain exploratory.",
        "External preprocessing, annotation, target-derived features, near-duplicate specimens, and unrecorded confounding cannot be ruled out by automated checks.",
        "No automatic imputation or batch correction is applied. Declared covariates are jointly regularized predictors, not unpenalized adjustment terms.",
    ]
    if config.imputation == "median":
        limits[-1] = (
            "Missing feature values are imputed from training-fold medians only. Missing cell populations can reflect recovery bias. No global batch correction is applied."
        )
    if config.panel_validation:
        limits[1] = (
            "Compact-panel selection and refitting are repeated within outer training folds. Panel metrics evaluate that discovery policy; the final frozen panel still needs untouched external evaluation."
        )
    if data.features:
        limits.append(
            "State proportions reflect relative sampled-cell recovery. Cell-level variances include measurement noise. Program/state definitions are declared in advance; supplied annotations can carry upstream bias."
        )
    simple_metrics = {
        k: v for k, v in metrics["out_of_fold"].items() if not isinstance(v, (dict, list))
    }
    preview = table[
        [
            "feature_label" if data.features else "gene_id",
            "contrast",
            "selection_frequency",
            "sign_consistency",
            "outer_selection_frequency",
            "coefficient",
            "passes_robustness_gates",
        ]
    ].head(25)
    warning_lines = manifest["warnings"] + manifest["gate_blockers"]
    md = [
        f"# {title}",
        "",
        summary,
        "",
        "## Held-out performance",
        "",
        "```json",
        json.dumps(simple_metrics, indent=2),
        "```",
        "",
        f"Proper-loss improvement over baseline: {metrics['loss_improvement']:.6g}.",
        f"Permutation p-value: {metrics['permutation'].get('pvalue')}; status: {metrics['permutation']['status']}.",
        "",
        "## Evidence gates and warnings",
        "",
        *(f"- {w}" for w in warning_lines),
        "",
        "## Interpretation boundaries",
        "",
        *(f"- {s}" for s in limits),
        "",
        "## Leading candidate evidence",
        "",
        "```csv",
        preview.to_csv(index=False).strip(),
        "```",
        "",
        "## Diagnostics",
        "",
    ]
    if config.panel_validation:
        panel_metrics = {
            k: v
            for k, v in metrics["panel"]["out_of_fold"].items()
            if not isinstance(v, (dict, list))
        }
        md.extend(
            [
                "## Compact-panel discovery",
                "",
                json.dumps(panel_metrics, indent=2),
                "",
                f"Final panel: {metrics['panel']['final_panel_size']} features. Nested panel permutation p: {metrics['panel']['permutation_pvalue']}.",
                "",
                "Selection occurs inside each outer training fold. This evaluates the panel-discovery policy; external validation uses the frozen panel_state.json.",
                "",
            ]
        )
    md.extend(f"![{caption}](figures/{name})\n\n{caption}\n" for name, caption in figures)
    (output / "report.md").write_text("\n".join(md) + "\n")
    style = "body{font:16px/1.55 system-ui;max-width:1100px;margin:40px auto;padding:0 24px;color:#1f2937}h1,h2{color:#164e63}table{border-collapse:collapse;font-size:13px}td,th{padding:6px;border-bottom:1px solid #ddd;text-align:left}figure{margin:24px 0}img{max-width:100%;height:auto}figcaption{color:#475569}pre{white-space:pre-wrap;background:#f1f5f9;padding:16px}.notice{border-left:5px solid #b45309;padding:10px 18px;background:#fffbeb}"
    content = [
        f"<!doctype html><html lang='en'><meta charset='utf-8'><title>{html.escape(title)}</title><style>{style}</style><body>",
        f"<h1>{html.escape(title)}</h1><p>{html.escape(summary)}</p>",
        "<div class='notice'><strong>Candidate association only.</strong> No clinical, biological-validation, or causal claims.</div>",
        "<h2>Held-out performance</h2>",
        pd.DataFrame([simple_metrics]).to_html(index=False, escape=True),
        f"<p>Loss improvement: {metrics['loss_improvement']:.6g}. Permutation p: {metrics['permutation'].get('pvalue')} ({html.escape(metrics['permutation']['status'])}).</p>",
        "<h2>Evidence gates and warnings</h2><ul>",
        *(f"<li>{html.escape(w)}</li>" for w in warning_lines),
        "</ul>",
        "<h2>Interpretation boundaries</h2><ul>",
        *(f"<li>{html.escape(s)}</li>" for s in limits),
        "</ul>",
        "<h2>Leading evidence</h2>",
        preview.to_html(index=False, escape=True, float_format=lambda x: f"{x:.3g}"),
        "<h2>Diagnostics</h2>",
    ]
    if config.panel_validation:
        content.extend(
            [
                "<h2>Compact-panel discovery</h2>",
                pd.DataFrame([panel_metrics]).to_html(index=False, escape=True),
                f"<p>Final panel: {metrics['panel']['final_panel_size']} features; nested permutation p: {metrics['panel']['permutation_pvalue']}. Selection was repeated inside the outer training folds.</p>",
            ]
        )
    for name, caption in figures:
        encoded = base64.b64encode((output / "figures" / name).read_bytes()).decode()
        content.append(
            f"<figure><img alt='{html.escape(caption, quote=True)}' src='data:image/png;base64,{encoded}'><figcaption>{html.escape(caption)}</figcaption></figure>"
        )
    content.append("</body></html>")
    (output / "report.html").write_text("\n".join(content))
    qc_text = (output / "qc_report.json").read_text()
    (output / "qc_report.html").write_text(
        f"<!doctype html><html lang='en'><meta charset='utf-8'><title>SigTrellis QC</title><style>{style}</style><body><h1>QC and input contract</h1><pre>{html.escape(qc_text)}</pre></body></html>"
    )
