"""Verify persisted scientific inputs against a trusted run manifest."""

from pathlib import Path
from typing import Any

from sigtrellis.domain import file_hash


def verify_run_artifacts(
    run: Path, manifest: dict[str, Any], names: tuple[str, ...]
) -> dict[str, str]:
    """Detect missing/changed files; this is integrity checking, not authentication."""
    recorded = manifest.get("output_hashes", {})
    verified: dict[str, str] = {}
    for name in names:
        expected = recorded.get(name)
        path = run / name
        if not isinstance(expected, str) or not path.is_file() or file_hash(path) != expected:
            raise ValueError(f"Integrity verification failed for training artifact: {name}")
        verified[name] = expected
    return verified
