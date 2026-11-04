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
                  checkpoint_every: int = 1, progress_every: int = 5_000,
                  starting_credit_balance_usd: float | None = None,
                  progress_callback: Callable[[str], None] | None = None,
                  stop_requested: Callable[[], bool] | None = None) -> pd.DataFrame:
    """Run missing tasks and retry transient failures until completion.

    A pair-gate handler returns `outcome`, `confidence`, and `evidence_span`.
    System 1 pathogen screens return a binary `answer` and `p_yes`; text-generating
    screens may also return an exact `evidence_span`. Retryable failures are
    checkpointed and retried with capped backoff; a non-retryable failure stops
    the run with an explicit error.
    """
    if max_attempts < 1:
        raise ValueError("max_attempts must be positive")
    if max_workers < 1:
        raise ValueError("max_workers must be positive")
    if checkpoint_every < 1:
        raise ValueError("checkpoint_every must be positive")
    if progress_every < 1:
        raise ValueError("progress_every must be positive")

    def report(message: str) -> None:
        if progress_callback is not None:
            progress_callback(message)
        else:
            print(message, flush=True)
    path, raw_path = Path(ledger_path), Path(raw_response_path)
    existing = pd.read_parquet(path) if path.exists() else pd.DataFrame()
    records = existing.to_dict("records") if not existing.empty else []
    completed_keys = {str(row["cache_key"]) for row in records
                      if row.get("status") == "completed" and row.get("cache_key")}
    seen_tasks: set[str] = set()
    cached_completed = 0
    pending = []
    for original in tasks:
        task = dict(original)
        cache_key = _identity(task)
        if cache_key in seen_tasks:
            continue
        seen_tasks.add(cache_key)
        if cache_key in completed_keys:
            cached_completed += 1
            continue
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
        if status == "completed" and task["scope"] == "pathogen_host_pair" and task.get("decision_kind") != "feature_screen":
            try:
                if task.get("decision_kind") == "pathogen_screen":
                    result = normalize_pathogen_screen(
                        pathogen_id=str(task["pathogen_id"]), host_id=str(task["host_id"]),
                        pmid=str(task["pmid"]), answer=result["answer"], confidence=result.get("confidence"),
                        evidence_text=result.get("evidence_span"), title=str(task.get("title_text", "")),
                        abstract=str(task.get("abstract_text", "")), run_id=str(task["run_id"]),
                        p_yes=result.get("p_yes"))
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
    total_pending = len(pending)
    total_tasks = len(seen_tasks)
    completed_this_run = 0
    cost_known_tasks = 0
    reported_cost_usd = 0.0
    next_progress_at = progress_every
    run_started = time.monotonic()

    def format_duration(seconds: float) -> str:
        total_seconds = max(0, int(seconds))
        hours, remainder = divmod(total_seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        if hours:
            return f"{hours}h {minutes}m"
        if minutes:
            return f"{minutes}m {seconds}s"
        return f"{seconds}s"

    with raw_path.open("a", encoding="utf-8") as raw_stream:
        stopped = False
        while pending:
            if stop_requested is not None and stop_requested():
                checkpoint()
                report("Pause requested before the next decision batch; saved decisions will resume later.")
                return pd.DataFrame(records)
            retry_round += 1
            if retry_round > 1:
                delay = min(300, 10 * (2 ** min(retry_round - 2, 5)))
                overall_completed = cached_completed + completed_this_run
                report(
                    f"Decision retry round {retry_round}: {len(pending):,} pending; "
                    f"{completed_this_run:,} completed this run; {cached_completed:,} cached; "
                    f"overall {overall_completed:,}/{total_tasks:,}; waiting {delay}s.")
                time.sleep(delay)
            else:
                start_message = (
                    f"Decision screening: {len(pending):,} pending; {cached_completed:,} cached; "
                    f"{total_tasks:,} total; workers={max_workers}.")
                if starting_credit_balance_usd is not None:
                    start_message += f" OpenRouter credits at start: ${starting_credit_balance_usd:,.4f}."
                else:
                    start_message += " OpenRouter balance unavailable (management key required for lookup)."
                report(start_message)

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
                            task_cost = result.get("cost_usd")
                            if task_cost is not None:
                                try:
                                    reported_cost_usd += max(0.0, float(task_cost))
                                    cost_known_tasks += 1
                                except (TypeError, ValueError):
                                    pass
                            completed_keys.add(cache_key)
                            completed_this_run += 1
                            if completed_this_run >= next_progress_at or completed_this_run == total_pending:
                                elapsed = max(time.monotonic() - run_started, 0.001)
                                rate = completed_this_run / elapsed
                                remaining = max(0, total_pending - completed_this_run)
                                eta = remaining / rate if rate > 0 else 0
                                percent = completed_this_run / total_pending * 100 if total_pending else 100
                                overall_completed = cached_completed + completed_this_run
                                overall_percent = overall_completed / total_tasks * 100 if total_tasks else 100
                                projected_run_cost = (reported_cost_usd / cost_known_tasks * total_pending
                                                      if cost_known_tasks else None)
                                projected_left = (starting_credit_balance_usd - projected_run_cost
                                                  if starting_credit_balance_usd is not None
                                                  and projected_run_cost is not None else None)
                                cost_estimate = (
                                    f" | reported cost ${reported_cost_usd:.4f}; "
                                    f"projected run cost ${projected_run_cost:.4f}; "
                                    f"predicted credits left ${max(0.0, projected_left):.4f}"
                                    if projected_run_cost is not None and projected_left is not None
                                    else (f" | reported cost ${reported_cost_usd:.4f}; "
                                          "credit forecast unavailable")
                                )
                                report(
                                    f"Decision progress: {completed_this_run:,}/{total_pending:,} this run "
                                    f"({percent:.1f}%) | cached {cached_completed:,} | "
                                    f"overall {overall_completed:,}/{total_tasks:,} ({overall_percent:.1f}%) | "
                                    f"{rate:,.1f} tasks/s | "
                                    f"elapsed {format_duration(elapsed)} | ETA ~{format_duration(eta)}"
                                    f"{cost_estimate}.")
                                next_progress_at = completed_this_run + progress_every
                        elif retryable:
                            retry_pending.append((task, cache_key))
                        else:
                            terminal_errors.append((str(task.get("pmid", "unknown")), str(error)))
                        round_completed += 1
                        if round_completed % checkpoint_every == 0:
                            checkpoint()
                    checkpoint()
                    if stop_requested is not None and stop_requested():
                        stopped = True
                        report("Pause requested; checkpointed completed decisions and will resume remaining tasks later.")
                        break
                if stopped:
                    break
            if terminal_errors:
                examples = "; ".join(f"PMID {pmid}: {error}" for pmid, error in terminal_errors[:3])
                raise RuntimeError(f"{len(terminal_errors):,} screening tasks failed with non-retryable errors; "
                                   f"results remain checkpointed. Examples: {examples}")
            pending = retry_pending
        if stopped:
            return pd.DataFrame(records)
    projected_end_balance = (starting_credit_balance_usd - reported_cost_usd
                             if starting_credit_balance_usd is not None else None)
    if projected_end_balance is not None:
        report(
            f"Decision screening complete: {completed_this_run:,} completed this run; "
            f"{cached_completed:,} cached; reported cost ${reported_cost_usd:.4f}; "
            f"predicted OpenRouter credits left ${max(0.0, projected_end_balance):.4f}.")
    else:
        report(
            f"Decision screening complete: {completed_this_run:,} completed this run; "
            f"{cached_completed:,} cached; reported cost ${reported_cost_usd:.4f}; "
            "credit balance unavailable for end-of-run estimate.")
    if records:
        checkpoint()
        return pd.DataFrame(records).drop_duplicates("cache_key", keep="last")
    return pd.DataFrame()
