"""Resumable, provider-agnostic decision ledger runner."""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable

import pandas as pd

from pathogen_atlas_v4.decision_identity import decision_cache_key
from pathogen_atlas_v4.evidence import normalize_pair_gate

REQUIRED_TASK_FIELDS = ("scope", "pathogen_id", "input_hash", "feature_id", "query_version",
                        "prompt_version", "feature_version", "question_version", "profile_version",
                        "model_version", "run_id", "run_config_version", "pmid")


def _identity(task: dict[str, Any]) -> str:
    missing = [key for key in REQUIRED_TASK_FIELDS if key not in task]
    if missing:
        raise ValueError(f"Decision task missing required fields: {missing}")
    return decision_cache_key(scope=task["scope"], pathogen_id=task["pathogen_id"],
                              host_id=task.get("host_id"), relationship_context=task.get("relationship_context"),
                              input_hash=task["input_hash"], feature_id=task["feature_id"],
                              query_version=task["query_version"], prompt_version=task["prompt_version"],
                              feature_version=task["feature_version"], question_version=task["question_version"],
                              profile_version=task["profile_version"], model_version=task["model_version"],
                              run_id=task["run_id"], run_config_version=task["run_config_version"])


def _atomic_parquet(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    frame.to_parquet(temporary, index=False)
    temporary.replace(path)


def run_decisions(tasks: Iterable[dict[str, Any]], handler: Callable[[dict[str, Any]], dict[str, Any]], *,
                  ledger_path: str | Path, raw_response_path: str | Path,
                  max_attempts: int = 2) -> pd.DataFrame:
    """Run missing tasks, retry failures, append raw responses, and atomically checkpoint.

    A handler returns a normalized result mapping. For pair-gate tasks it includes
    `outcome`, `confidence`, and `evidence_span`; all other scopes can return a
    feature-specific mapping. Failed rows remain retryable on the next invocation.
    """
    if max_attempts < 1:
        raise ValueError("max_attempts must be positive")
    path, raw_path = Path(ledger_path), Path(raw_response_path)
    existing = pd.read_parquet(path) if path.exists() else pd.DataFrame()
    records = existing.to_dict("records") if not existing.empty else []
    completed_keys = {str(row["cache_key"]) for row in records
                      if row.get("status") == "completed" and row.get("cache_key")}
    seen_tasks: set[str] = set()
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    with raw_path.open("a", encoding="utf-8") as raw_stream:
        for task in tasks:
            task = dict(task)
            cache_key = _identity(task)
            if cache_key in seen_tasks or cache_key in completed_keys:
                continue
            seen_tasks.add(cache_key)
            result: dict[str, Any] = {}
            error = None
            attempts = 0
            for attempts in range(1, max_attempts + 1):
                try:
                    result = dict(handler(task))
                    error = None
                    break
                except Exception as exc:  # provider errors remain data, not silent omissions
                    error = f"{type(exc).__name__}: {exc}"
            status = "completed" if error is None else "failed"
            if status == "completed" and task["scope"] == "pathogen_host_pair":
                try:
                    result = normalize_pair_gate(
                        pathogen_id=str(task["pathogen_id"]), host_id=str(task["host_id"]),
                        pmid=str(task["pmid"]), outcome=result["outcome"],
                        confidence=result.get("confidence"), evidence_text=result.get("evidence_span"),
                        title_hash=str(task.get("title_hash", "")), abstract_hash=str(task.get("abstract_hash", "")),
                        run_id=str(task["run_id"]), relationship_context=str(task.get("relationship_context", "disease")),
                        provenance_status="v4_native")
                except Exception as exc:
                    status = "failed"
                    error = f"{type(exc).__name__}: {exc}"
            raw_record = {"cache_key": cache_key, "task": task, "result": result,
                          "error": error, "attempts": attempts,
                          "recorded_at": datetime.now(timezone.utc).isoformat()}
            raw_stream.write(json.dumps(raw_record, ensure_ascii=False, default=str) + "\n")
            raw_stream.flush()
            records.append({**task, **result, "cache_key": cache_key, "status": status,
                            "attempts": attempts, "error": error,
                            "raw_response_path": str(raw_path),
                            "completed_at": datetime.now(timezone.utc).isoformat() if status == "completed" else None})
            _atomic_parquet(pd.DataFrame(records), path)
            if status == "completed":
                completed_keys.add(cache_key)
    return pd.DataFrame(records)
