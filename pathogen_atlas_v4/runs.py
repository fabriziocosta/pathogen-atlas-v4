"""Run configuration validation and explicit pair expansion."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

from pathogen_atlas_v4.decision_identity import stable_hash


@dataclass(frozen=True)
class ResolvedRun:
    run_id: str
    config_version: str
    host_ids: tuple[str, ...]
    pathogen_ids: tuple[str, ...]
    pairs: tuple[tuple[str, str], ...]
    pair_generation: str
    resolved_config_hash: str


def resolve_run(config: dict[str, Any], *, known_host_ids: set[str], known_pathogen_ids: set[str],
                host_group_ids: Mapping[str, Iterable[str]] | None = None,
                pathogen_group_ids: Mapping[str, Iterable[str]] | None = None,
                pathogen_taxon_ids: Mapping[str, Iterable[str]] | None = None) -> ResolvedRun:
    """Validate selectors and resolve only explicitly requested pair rules."""
    run = config.get("run", {})
    version = str(run.get("config_version", "")).strip()
    run_id = str(run.get("run_id", "")).strip()
    if not version or not run_id:
        raise ValueError("run.config_version and run.run_id are required")
    selected_hosts = list(map(str, run.get("host_ids", [])))
    for group in run.get("host_groups", []):
        if host_group_ids is None or str(group) not in host_group_ids:
            raise ValueError(f"Unknown host_group selector {group!r}")
        selected_hosts.extend(map(str, host_group_ids[str(group)]))
    selected_pathogens = list(map(str, run.get("pathogen_ids", [])))
    for group in run.get("pathogen_groups", []):
        if pathogen_group_ids is None or str(group) not in pathogen_group_ids:
            raise ValueError(f"Unknown pathogen_group selector {group!r}")
        selected_pathogens.extend(map(str, pathogen_group_ids[str(group)]))
    for taxon_id in run.get("pathogen_taxon_ids", []):
        if pathogen_taxon_ids is None or str(taxon_id) not in pathogen_taxon_ids:
            raise ValueError(f"Unknown pathogen taxon selector {taxon_id!r}")
        selected_pathogens.extend(map(str, pathogen_taxon_ids[str(taxon_id)]))
    hosts = tuple(dict.fromkeys(selected_hosts))
    pathogens = tuple(dict.fromkeys(selected_pathogens))
    missing_hosts = sorted(set(hosts) - known_host_ids)
    missing_pathogens = sorted(set(pathogens) - known_pathogen_ids)
    if missing_hosts or missing_pathogens:
        raise ValueError(f"Unknown selected IDs: hosts={missing_hosts}; pathogens={missing_pathogens}")
    pair_generation = str(run.get("pair_generation", "explicit_pairs"))
    explicit_pairs = run.get("pairs")
    if explicit_pairs is not None:
        pairs = tuple(dict.fromkeys((str(p["pathogen_id"]), str(p["host_id"])) for p in explicit_pairs))
        if any(pathogen not in known_pathogen_ids or host not in known_host_ids for pathogen, host in pairs):
            raise ValueError("Explicit pairs must reference known pathogen and host IDs")
        if hosts and set(host for _, host in pairs) - set(hosts):
            raise ValueError("Explicit pairs include a host outside run.host_ids")
        if pathogens and set(pathogen for pathogen, _ in pairs) - set(pathogens):
            raise ValueError("Explicit pairs include a pathogen outside run.pathogen_ids")
    elif pair_generation == "cartesian_product":
        if not hosts or not pathogens:
            raise ValueError("cartesian_product requires nonempty host_ids and pathogen_ids")
        pairs = tuple((p, h) for p, h in product(pathogens, hosts))
    else:
        raise ValueError("Provide explicit run.pairs or set pair_generation=cartesian_product")
    return ResolvedRun(run_id, version, hosts, pathogens, pairs, pair_generation,
                       stable_hash({"source_config": config, "host_ids": hosts,
                                    "pathogen_ids": pathogens, "pairs": pairs,
                                    "pair_generation": pair_generation}))


def save_resolved_run(config: dict[str, Any], resolved: ResolvedRun, path: str | Path) -> None:
    """Persist selected entities and the exact pair list before acquisition begins."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    record = {"source_config": config, "run_id": resolved.run_id,
              "config_version": resolved.config_version, "resolved_config_hash": resolved.resolved_config_hash,
              "host_ids": list(resolved.host_ids), "pathogen_ids": list(resolved.pathogen_ids),
              "pair_generation": resolved.pair_generation,
              "pairs": [{"pathogen_id": pathogen_id, "host_id": host_id, "relationship_context": "disease"}
                        for pathogen_id, host_id in resolved.pairs]}
    destination.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
