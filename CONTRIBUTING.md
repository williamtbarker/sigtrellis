# Contributing

Use Python 3.12+ and install `python -m pip install -e '.[dev,de,demo]'`.
Run Ruff lint/format, strict mypy and pytest before proposing a change. CI also
builds and installs a wheel outside the source tree. Public downloads stay in
ignored `data/`; generated reports stay in ignored `results/` or release assets.

For a statistical change, state the estimand, independent biological unit,
training boundary and assumptions. Add a meaningful planted-signal or negative
control when behavior changes. A faster method must retain leakage-safe fitting;
a higher AUC is not by itself evidence of an improvement. Do not adjust example
settings against the held-out cohorts or hide failed/empty panels.

For a new single-cell representation, use explicit feature identities and units.
Any learned reference, embedding, graph or program must be trained within folds
and have a frozen rule for assigning held-out specimens/cells. More cells must
not be counted as additional phenotype replicates. Include sparse-memory behavior
and appropriate controls for cell recovery, batch and donor effects.

Bug reports should include the version, configuration, traceback, expected
behavior and a small synthetic reproducer where possible. Reports and manifests
can contain specimen IDs and metadata; review those before making them public.
Never contribute proprietary source or data without appropriate rights.
