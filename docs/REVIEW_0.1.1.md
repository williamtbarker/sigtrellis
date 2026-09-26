# Adversarial source review and corrections — 0.1.1

Reviewed 2026-09-26, following the first delivered release. The original archive was recovered and all 186 recorded file hashes matched before inspection. This review examined scientific assumptions, identity handling, validation boundaries, numerical failure modes, persisted artifacts and reporting. It is a second review by the implementing assistant, **not an independent human or external laboratory review**.

## Findings and disposition

| Finding | Severity / consequence | Correction and targeted evidence |
|---|---|---|
| CSV inference changed identifiers such as `001`, `01`, `1` and literal `NA` | High: incorrect biological grouping, mismatched axes or changed outcome categories | Parse identity columns as exact strings and construct the index afterward. Test both matrix orientations. Reject blank required metadata. |
| Exact duplicate detection distinguished `0.0` from `-0.0` | High edge case: numerically identical profiles could evade duplicate checks | Canonicalize signed zero before hashing. A cross-group duplicate with a changed zero sign must fail QC. |
| Paired fixed-effect DE accepted multiple observations in one donor/condition | High: an unsupported repeated-measures design could overstate independent residual evidence | Require at most one observation per biological group/condition. Reject before reaching PyDESeq2; technical replicates must be combined appropriately. |
| External validation trusted changed model/configuration/metadata files | Medium: a modified model could be reported as the original frozen fit | Verify all three artifacts against the saved manifest before cohort evaluation. Tests modify each artifact separately. Record the verified model fingerprint and label probability columns by class. |
| Correlation summaries pooled multinomial contrast selection | Medium: selection for one contrast inflated another contrast's group evidence | Keep the contrast axis in group and substitution frequencies; explicit coefficient-axis labels are required for multiclass diagnostics. A 100%-versus-16.7% contrast fixture must remain distinct. |
| NaN/infinite statistics bypassed comparison-based gates | High latent failure: nonfinite comparisons could produce a minimum permutation p-value or permit a positive evidence gate | Reject nonfinite observed/null statistics and invalid p-values; fail confidence gates explicitly. These were adversarial injected-statistic cases, not failures observed in the public demonstrations. |
| Missing-class folds and constant regression folds emitted misleading metric values | Medium: apparent balanced discrimination or finite R² without the required outcome variation | Report balanced accuracy/MCC as undefined when classes are missing and R² as undefined for constant outcomes. Proper loss and Brier remain available. |
| H5AD matrix dimensions were not compared with annotation dimensions | Medium: extra matrix rows could be silently ignored | Verify cell/gene dimensions before aggregation. A deliberately inconsistent sparse shape must fail. |
| Boolean numeric settings and duplicate grid entries were accepted | Medium: accidental YAML booleans or duplicated tuning evidence | Reject both. Duplicate entries must not inflate the fold count used by the one-SE heuristic. |
| Batch model tuning/coefficients were incompletely preserved | Medium: reported robustness was harder to reproduce independently | Save chosen and baseline hyperparameters, complete batch tuning rows, coefficient arrays and batch-axis labels. A complete batch-run test checks artifacts. |
| Filtering-stage counts and final nonzero manifest definitions were inconsistent | Low: misleading audit summaries | Record pre-screen and post-screen counts separately; use the declared coefficient tolerance for manifest selection. |

Every finding above is corrected in the delivered source. None is resolved by loosening scientific evidence thresholds or selecting a more favorable dataset.

## Reproduction and verification

`tests/test_review_regressions.py` adds **27 cases**. Against the original 0.1.0 code, **26 adversarial cases fail**, while a real PyDESeq2 nested-screen integration test passes as a positive control. The baseline evidence is preserved in `docs/validation/baseline_regressions.log` (25 failures, one pass) and `baseline_signed_zero.log` (one further failure). These are test-case counts across the finding classes above, not a claim of 26 unrelated defects.

The revised full suite has **116 passing tests**, including all original tests and every added review case, in both development and a clean installed-wheel environment. Branch-aware coverage is **89%**. Ruff and strict mypy pass. The count-DE integration now runs real negative-binomial fits inside nested training folds, complementing the existing boundary spies.

The final wheel is installed outside the source tree. Both public workflows and the synthetic quickstart run again with delivered source hashes. Current metrics, manifests, complete fit/split evidence and outputs are in `release_evidence/` and summarized in `VALIDATION.md`. Original execution records are retained under `docs/validation/v0.1.0/`; they are historical, not presented as 0.1.1 results.

## Scientific conclusions that survived review

- The main nested path refits normalization references, filters, supervised screens, scales and covariate categories within each training boundary.
- Donor/group splitting and subsampling preserve biological units. Single-cell counts are aggregated before modeling.
- The synthetic quickstart still recovers exactly the four planted genes.
- Public bulk and single-cell workflows still complete, while their small cohorts and the single-cell library confound prevent positive confidence gates.
- Same-cohort DE remains exploratory. Procedure CV is not assigned to the post-hoc compact panel.

## Remaining limitations

Passing this battery is not proof of universal correctness. Hidden upstream preprocessing, mislabeled or related donors, near-duplicates, partial confounding, invalid permutation exchangeability and analyst-level selection across many runs remain material risks. Artifact hashes assume the manifest itself is trusted; a person able to change both the manifest and files can defeat unsigned integrity checks. The software does not claim adversarial cryptographic authentication.

Covariates remain jointly penalized predictors, and conditional gene inference is not implemented. General repeated-measures count models, temporal validation, formal feature-level error control, panel-specific nested discovery and larger independent cohorts remain research work. Remote GitHub Actions and Docker execution were not performed in this review.

**Recommendation:** release publicly as a scoped research package with this review, negative-control results and limitations. Independent statistical review and external panel validation remain necessary before stronger scientific or clinical claims.
