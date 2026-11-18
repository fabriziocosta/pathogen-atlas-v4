"""Baseline discovery and plan creation for incremental dashboard runs."""

from __future__ import annotations

import json
from datetime import datetime, timedelta, timezone
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from pathogen_atlas_v4.pubmed_client import literature_cache_key


COMPATIBLE_CONFIG_PATHS = (
    ("taxonomy", "staging_snapshot"),
    ("features", "profile_version"),
    ("acquisition", "query_strategy_version"),
    ("acquisition", "pathogen_literature_query_version"),
    ("decisions", "prompt_set_version"),
    ("decisions", "pathogen_screen_question_version"),
    ("decisions", "pathogen_screen_prompt_version"),
    ("decisions", "pair_feature_prompt_version"),
    ("decisions", "model_version"),
)


def _nested(config: dict[str, Any], path: tuple[str, ...]) -> Any:
    value: Any = config
    for part in path:
        value = value.get(part) if isinstance(value, dict) else None
    return value


def _load_run_config(run_root: Path) -> dict[str, Any]:
    dashboard_config = run_root / "dashboard-run.yaml"
    if dashboard_config.is_file():
        value = yaml.safe_load(dashboard_config.read_text(encoding="utf-8")) or {}
        if isinstance(value, dict):
            return value
    resolved_path = run_root / "resolved_config.json"
    if resolved_path.is_file():
        value = json.loads(resolved_path.read_text(encoding="utf-8"))
        return value.get("source_config", {})
    return {}


def _compatible(candidate: dict[str, Any], active: dict[str, Any]) -> bool:
    return all(_nested(candidate, path) == _nested(active, path) for path in COMPATIBLE_CONFIG_PATHS)


def _cache_key(config: dict[str, Any]) -> str:
    acquisition = config.get("acquisition", {})
    return literature_cache_key(
        query_version=str(acquisition.get("query_strategy_version", "")),
        pathogen_literature_query_version=str(acquisition.get("pathogen_literature_query_version", "")),
        taxonomy_snapshot=str(config.get("taxonomy", {}).get("staging_snapshot", "")),
        host_expansion_hash=str(acquisition.get("host_expansion_hash", "")),
    )


def _pair_rows(path: Path) -> list[tuple[str, str]]:
    if not path.is_file():
        return []
    try:
        frame = pd.read_parquet(path, columns=["pathogen_id", "host_id"])
    except (OSError, ValueError, KeyError):
        return []
    return list(dict.fromkeys((str(p), str(h)) for p, h in
                              frame[["pathogen_id", "host_id"]].dropna().itertuples(index=False, name=None)))


def _pair_search_index(path: Path) -> dict[tuple[str, str], str]:
    if not path.is_file():
        return {}
    try:
        frame = pd.read_parquet(path, columns=["pathogen_id", "host_id", "query_hash", "status"])
    except (OSError, ValueError, KeyError):
        return {}
    if not frame.empty:
        frame = frame.loc[frame.status.astype(str).eq("completed")]
    return {(str(row.pathogen_id), str(row.host_id)): str(row.query_hash)
            for row in frame.itertuples(index=False)}


@lru_cache(maxsize=12)
def _pathogen_ids_for_groups(project_root: str, taxonomy_snapshot: str,
                             groups: tuple[str, ...]) -> tuple[str, ...]:
    """Resolve the same pinned staged plus supplemental IDs used by notebook 00."""
    project = Path(project_root)
    taxonomy_root = project / "data/v4/staging" / taxonomy_snapshot / "taxonomy"
    base_pathogens = pd.read_parquet(taxonomy_root / "pathogens.parquet")
    base_aliases = pd.read_parquet(taxonomy_root / "pathogen_aliases.parquet")
    manifest = json.loads((project / "data/source_taxonomy/ncbi/source_manifest.json").read_text(encoding="utf-8"))
    release = f"NCBI Taxonomy {manifest['retrieved_utc']} sha256:{manifest['archive_sha256']}"
    from pathogen_atlas_v4.pathogen_taxonomy import load_ncbi_pathogen_additions, merge_pathogen_catalogues
    additions, addition_aliases, scoped_taxa = load_ncbi_pathogen_additions(
        project / "data/source_taxonomy/ncbi", groups, release=release,
        source_version=f"ncbi-taxdump:{manifest['archive_sha256']}")
    pathogens, _ = merge_pathogen_catalogues(
        base_pathogens, base_aliases, additions, addition_aliases, scoped_taxa)
    selected = pathogens.loc[pathogens.pathogen_group.astype(str).isin(groups), "pathogen_id"]
    return tuple(dict.fromkeys(selected.dropna().astype(str).tolist()))


def build_update_plan(*, root: str | Path, state_root: str | Path,
                      hosts: list[str], groups: list[str], taxonomy_snapshot: str,
                      overlap_days: int = 7, now: datetime | None = None) -> dict[str, Any]:
    """Choose latest compatible completed coverage independently for each pair."""
    project = Path(root)
    states = Path(state_root)
    current_time = now or datetime.now(timezone.utc)
    if current_time.tzinfo is None:
        current_time = current_time.replace(tzinfo=timezone.utc)
    end = current_time.astimezone(timezone.utc).date()
    active = yaml.safe_load((project / "configs/pilot_v1.yaml").read_text(encoding="utf-8")) or {}
    active = json.loads(json.dumps(active))
    active.setdefault("taxonomy", {})["staging_snapshot"] = taxonomy_snapshot

    selected_pathogens = _pathogen_ids_for_groups(str(project.resolve()), taxonomy_snapshot, tuple(sorted(groups)))
    requested_pairs = [(str(pathogen), str(host)) for pathogen in selected_pathogens for host in hosts]
    requested = set(requested_pairs)

    candidates = []
    for manifest_path in states.glob("*/manifest.json"):
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        config_meta = manifest.get("config", {})
        if (manifest.get("status") != "completed" or config_meta.get("mode") != "full"
                or config_meta.get("run_type", "full") not in {"full", "incremental_update"}):
            continue
        if config_meta.get("taxonomy_snapshot") != taxonomy_snapshot:
            continue
        run_id = str(manifest.get("run_id") or manifest_path.parent.name)
        run_root = project / "data/v4/runs" / run_id
        run_config = _load_run_config(run_root)
        if not run_config or not _compatible(run_config, active):
            continue
        cache_key = _cache_key(run_config)
        cache_root = project / "data/v4/cache/literature" / cache_key
        pair_table = cache_root / "tables/pair_searches.parquet"
        query_hashes = _pair_search_index(pair_table)
        if not query_hashes:
            continue
        search_pairs = set(query_hashes)
        prior_pairs = _pair_rows(run_root / "pathogen_host_selection.parquet")
        covered = requested & set(prior_pairs) & search_pairs
        if not covered:
            continue
        expansions = run_config.get("acquisition", {}).get("host_search_expansions", {})
        window_end = config_meta.get("update_window_end") or manifest.get("completed_at")
        if not window_end:
            continue
        candidates.append({
            "run_id": run_id,
            "completed_at": str(manifest.get("completed_at", manifest.get("created_at", ""))),
            "window_end": str(window_end),
            "cache_key": cache_key,
            "cache_root": str(cache_root),
            "config": run_config,
            "host_expansions": expansions,
            "query_hashes": query_hashes,
            "covered_pairs": covered,
        })
    candidates.sort(key=lambda item: (item["completed_at"], item["run_id"]), reverse=True)

    host_expansions: dict[str, dict[str, Any]] = {}
    for candidate in candidates:
        for host_id, expansion in candidate["host_expansions"].items():
            host_expansions.setdefault(str(host_id), expansion)

    baseline_by_pair = {}
    for candidate in candidates:
        for pair in candidate["covered_pairs"]:
            if pair not in requested or pair in baseline_by_pair:
                continue
            expansion = candidate["host_expansions"].get(pair[1])
            if expansion is None or expansion != host_expansions.get(pair[1]):
                continue
            baseline_by_pair[pair] = candidate

    plan_pairs = []
    for pair in requested_pairs:
        candidate = baseline_by_pair.get(pair)
        if candidate is None:
            continue
        prior_config = candidate["config"]
        end_value = candidate["window_end"]
        try:
            prior_end = datetime.fromisoformat(end_value.replace("Z", "+00:00"))
        except ValueError:
            prior_end = datetime.strptime(end_value[:10], "%Y-%m-%d").replace(tzinfo=timezone.utc)
        if prior_end.tzinfo is None:
            prior_end = prior_end.replace(tzinfo=timezone.utc)
        start = (prior_end.astimezone(timezone.utc).date() - timedelta(days=overlap_days))
        plan_pairs.append({
            "pathogen_id": pair[0], "host_id": pair[1],
            "baseline_run_id": candidate["run_id"],
            "baseline_cache_key": candidate["cache_key"],
            "baseline_query_hash": candidate["query_hashes"][pair],
            "window_start": start.isoformat(), "window_end": end.isoformat(),
            "previous_window_end": candidate["window_end"],
            "baseline_profile_version": _nested(prior_config, ("features", "profile_version")),
        })
    covered_keys = {(item["pathogen_id"], item["host_id"]) for item in plan_pairs}
    missing = requested - covered_keys
    starts = [item["window_start"] for item in plan_pairs]
    return {
        "created_at": current_time.astimezone(timezone.utc).isoformat(),
        "window_end": end.isoformat(), "overlap_days": overlap_days,
        "taxonomy_snapshot": taxonomy_snapshot,
        "hosts": hosts, "groups": groups,
        "pairs": plan_pairs,
        "host_search_expansions": {key: value for key, value in host_expansions.items() if key in hosts},
        "baseline_run_ids": sorted({item["baseline_run_id"] for item in plan_pairs}),
        "requested_pair_count": len(requested_pairs),
        "covered_pair_count": len(plan_pairs),
        "uncovered_pair_count": len(missing),
        "uncovered_pair_examples": [list(pair) for pair in sorted(missing)[:20]],
        "earliest_window_start": min(starts) if starts else None,
    }


def write_update_plan(path: str | Path, plan: dict[str, Any]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    meta = {key: value for key, value in plan.items() if key != "pairs"}
    with target.open("w", encoding="utf-8") as stream:
        stream.write(json.dumps({"_meta": meta}, ensure_ascii=False) + "\n")
        for pair in plan.get("pairs", []):
            stream.write(json.dumps(pair, ensure_ascii=False, sort_keys=True) + "\n")


def read_update_plan(path: str | Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    meta: dict[str, Any] = {}
    pairs = []
    with Path(path).open(encoding="utf-8") as stream:
        for line in stream:
            if not line.strip():
                continue
            item = json.loads(line)
            if "_meta" in item:
                meta = item["_meta"]
            else:
                pairs.append(item)
    return meta, pairs
