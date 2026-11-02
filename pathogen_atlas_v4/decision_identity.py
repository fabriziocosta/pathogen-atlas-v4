"""Deterministic input and cache identities with scope-aware host keys."""

from __future__ import annotations

import hashlib
import json
from typing import Any

from pathogen_atlas_v4.schemas import validate_scope


def stable_hash(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":"), default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def decision_cache_key(*, scope: str, pathogen_id: str, input_hash: str,
                       feature_id: str, question_version: str, profile_version: str,
                       model_version: str, host_id: str | None = None,
                       relationship_context: str | None = None, query_version: str = "",
                       prompt_version: str = "", feature_version: str = "",
                       run_id: str = "", run_config_version: str = "") -> str:
    scope = validate_scope(scope)
    # Run ID is provenance, not decision identity; matching inputs can be reused
    # across smoke and full runs when the versioned decision context is equal.
    identity = {"scope": scope, "pathogen_id": str(pathogen_id), "input_hash": str(input_hash),
                "feature_id": str(feature_id), "question_version": str(question_version),
                "feature_version": str(feature_version), "profile_version": str(profile_version),
                "query_version": str(query_version), "prompt_version": str(prompt_version),
                "model_version": str(model_version),
                "run_config_version": str(run_config_version)}
    if scope in {"host", "pathogen_host_pair"} and not host_id:
        raise ValueError(f"host_id is required for {scope} decisions")
    if scope == "host":
        identity.pop("pathogen_id")
        identity["host_id"] = str(host_id)
    elif scope == "pathogen_host_pair":
        identity["host_id"] = str(host_id)
        identity["relationship_context"] = relationship_context or "disease"
    return stable_hash(identity)
