"""Resumable, provider-agnostic decision ledger runner."""

from __future__ import annotations

import json
import os
import time
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Callable, Iterable

import pandas as pd

from pathogen_atlas_v4.decision_identity import decision_cache_key
from pathogen_atlas_v4.evidence import normalize_pair_gate, normalize_pathogen_screen

REQUIRED_TASK_FIELDS = ("scope", "pathogen_id", "input_hash", "feature_id", "query_version",
                        "prompt_version", "feature_version", "question_version", "profile_version",
                        "model_version", "run_id", "run_config_version", "pmid")


def decision_cache_key_for_task(task: dict[str, Any]) -> str:
    """Return the runner's stable identity for one fully versioned task."""
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


def _identity(task: dict[str, Any]) -> str:
    return decision_cache_key_for_task(task)


def _atomic_parquet(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    frame.to_parquet(temporary, index=False)
    temporary.replace(path)


def run_decisions(tasks: Iterable[dict[str, Any]], handler: Callable[[dict[str, Any]], dict[str, Any]], *,
                  ledger_path: str | Path, raw_response_path: str | Path,
                  max_attempts: int = 2, max_workers: int = 1,
                  checkpoint_every: int = 1) -> pd.DataFrame:
    """Run missing tasks and retry transient failures until completion.

    A handler returns a normalized result mapping. For pair-gate tasks it includes
    `outcome`, `confidence`, and `evidence_span`; all other scopes can return a
    feature-specific mapping. Retryable failures are checkpointed and retried with
    capped backoff; a non-retryable failure stops the run with an explicit error.
    """
    if max_attempts < 1:
        raise ValueError("max_attempts must be positive")
    if max_workers < 1:
        raise ValueError("max_workers must be positive")
    if checkpoint_every < 1:
        raise ValueError("checkpoint_every must be positive")
    path, raw_path = Path(ledger_path), Path(raw_response_path)
    existing = pd.read_parquet(path) if path.exists() else pd.DataFrame()
    records = existing.to_dict("records") if not existing.empty else []
    completed_keys = {str(row["cache_key"]) for row in records
                      if row.get("status") == "completed" and row.get("cache_key")}
    seen_tasks: set[str] = set()
    pending = []
    for original in tasks:
        task = dict(original)
        cache_key = _identity(task)
        if cache_key in seen_tasks or cache_key in completed_keys:
            continue
        seen_tasks.add(cache_key)
        pending.append((task, cache_key))

    def execute(task: dict[str, Any]) -> tuple[dict[str, Any], str | None, int, str, bool]:
        result: dict[str, Any] = {}
        error = None
        attempts = 0
        retryable = True
        for attempts in range(1, max_attempts + 1):
            try:
                result = dict(handler(task))
                error = None
                break
            except Exception as exc:  # provider errors remain data, not silent omissions
                error = f"{type(exc).__name__}: {exc}"
                retryable = bool(getattr(exc, "retryable", not isinstance(exc, (ValueError, TypeError, KeyError))))
                if not retryable:
                    break
        status = "completed" if error is None else "failed"
        if status == "completed" and task["scope"] == "pathogen_host_pair":
            try:
                if task.get("decision_kind") == "pathogen_screen":
                    result = normalize_pathogen_screen(
                        pathogen_id=str(task["pathogen_id"]), host_id=str(task["host_id"]),
                        pmid=str(task["pmid"]), answer=result["answer"], confidence=result.get("confidence"),
                        evidence_text=result.get("evidence_span"), title=str(task.get("title_text", "")),
                        abstract=str(task.get("abstract_text", "")), run_id=str(task["run_id"]))
                else:
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
                retryable = False
        return result, error, attempts, status, retryable

    attempts_by_key = {str(row.get("cache_key")): int(row.get("attempts") or 0)
                       for row in records if row.get("cache_key")}
    raw_path.parent.mkdir(parents=True, exist_ok=True)

    def checkpoint() -> None:
        if not records:
            return
        latest = pd.DataFrame(records)
        if "cache_key" in latest:
            latest = latest.drop_duplicates("cache_key", keep="last")
        _atomic_parquet(latest, path)

    retry_round = 0
    with raw_path.open("a", encoding="utf-8") as raw_stream:
        while pending:
            retry_round += 1
            if retry_round > 1:
                delay = min(300, 10 * (2 ** min(retry_round - 2, 5)))
                print(f"Decision retry round {retry_round}: {len(pending):,} tasks remain; waiting {delay}s.", flush=True)
                time.sleep(delay)
            else:
                print(f"Decision screening: {len(pending):,} tasks pending; workers={max_workers}.", flush=True)

            retry_pending: list[tuple[dict[str, Any], str]] = []
            terminal_errors: list[tuple[str, str]] = []
            round_completed = 0
            batch_size = max_workers * 2
            with ThreadPoolExecutor(max_workers=max_workers) as pool:
                for batch_start in range(0, len(pending), batch_size):
                    batch = pending[batch_start:batch_start + batch_size]
                    futures = {pool.submit(execute, task): (task, cache_key)
                               for task, cache_key in batch}
                    for future in as_completed(futures):
                        task, cache_key = futures[future]
                        result, error, attempts, status, retryable = future.result()
                        total_attempts = attempts_by_key.get(cache_key, 0) + attempts
                        attempts_by_key[cache_key] = total_attempts
                        raw_record = {"cache_key": cache_key, "task": task, "result": result,
                                      "error": error, "attempts": attempts, "retryable": retryable,
                                      "retry_round": retry_round,
                                      "recorded_at": datetime.now(timezone.utc).isoformat()}
                        raw_stream.write(json.dumps(raw_record, ensure_ascii=False, default=str) + "\n")
                        raw_stream.flush()
                        records.append({**task, **result, "cache_key": cache_key, "status": status,
                                        "attempts": total_attempts, "error": error, "retryable": retryable,
                                        "retry_round": retry_round, "raw_response_path": str(raw_path),
                                        "completed_at": datetime.now(timezone.utc).isoformat()
                                        if status == "completed" else None})
                        if status == "completed":
                            completed_keys.add(cache_key)
                        elif retryable:
                            retry_pending.append((task, cache_key))
                        else:
                            terminal_errors.append((str(task.get("pmid", "unknown")), str(error)))
                        round_completed += 1
                        if round_completed % checkpoint_every == 0:
                            checkpoint()
                    checkpoint()
            if terminal_errors:
                examples = "; ".join(f"PMID {pmid}: {error}" for pmid, error in terminal_errors[:3])
                raise RuntimeError(f"{len(terminal_errors):,} screening tasks failed with non-retryable errors; "
                                   f"results remain checkpointed. Examples: {examples}")
            pending = retry_pending
    if records:
        checkpoint()
        return pd.DataFrame(records).drop_duplicates("cache_key", keep="last")
    return pd.DataFrame()
