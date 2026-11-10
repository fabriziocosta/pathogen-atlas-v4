#!/usr/bin/env python3
"""Restore pathogenicity fractions from retained screened sample counts."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SITE = ROOT / "sites/pathogen-atlas-explorer"
INDEX_VERSION = "pathogenic-interest-index-v4.6-pair-query-fraction-scores"


def _read_gzip_json(path: Path) -> Any:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        return json.load(stream)


def _write_gzip_json(path: Path, payload: Any) -> None:
    encoded = json.dumps(payload, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with temporary.open("wb") as output:
            with gzip.GzipFile(fileobj=output, mode="wb", mtime=0) as compressed:
                compressed.write(encoded)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _sample_fraction(profile: dict[str, Any]) -> float | None:
    if str(profile.get("pair_selection_status", "")).lower() not in {"accepted", "supported"}:
        return None
    positive = int(profile.get("n_positive_records") or profile.get("n_jev_accepted_pmids")
                   or profile.get("n_pathogen_pmids") or 0)
    sampled = int(profile.get("n_sampled_records") or profile.get("n_abstracts_fetched") or 0)
    if sampled <= 0:
        return None
    if positive < 0 or positive > sampled:
        raise ValueError(f"Invalid sample counts: {positive} positive of {sampled} sampled")
    return positive / sampled


def restore(site_root: Path = DEFAULT_SITE) -> dict[str, Any]:
    site_root = Path(site_root)
    source_path = site_root / "source/human_pair_data.json.gz"
    source = _read_gzip_json(source_path)
    source_restored = 0
    human_host_id = str(source.get("host_id") or "")
    human_profiles = {str(profile.get("pathogen_id") or profile.get("organism_id") or ""): profile
                      for profile in source.get("profiles", [])}
    for profile in source.get("profiles", []):
        sampled = int(profile.get("n_abstracts_fetched") or profile.get("n_sampled_records") or 0)
        positive = int(profile.get("n_pathogen_pmids") or profile.get("n_jev_accepted_pmids") or 0)
        if sampled <= 0:
            continue
        if positive < 0 or positive > sampled:
            raise ValueError(f"Invalid Human sample counts for {profile.get('pathogen_id')}: {positive}/{sampled}")
        profile.setdefault("legacy_pathogenicity_fraction", profile.get("pathogenicity_fraction"))
        profile["pathogenicity_fraction"] = positive / sampled
        pair_total = profile.get("n_pair_search_records")
        if pair_total is None:
            pair_total = profile.get("n_pubmed_search_records")
        if pair_total is not None:
            profile["n_pair_search_records"] = int(pair_total)
            profile["estimated_pathogen_publications"] = positive / sampled * int(pair_total)
        source_restored += 1
    _write_gzip_json(source_path, source)

    data_path = site_root / "explorer-data.json.gz"
    data = _read_gzip_json(data_path)
    columns = data["species_columns"]
    by_group: dict[str, list[list[Any]]] = {}
    fraction_logs: dict[str, list[float]] = {}
    volume_logs: dict[str, list[float]] = {}
    restored_by_group: dict[str, int] = {}
    for group in data["groups"]:
        shard_path = site_root / "taxa" / f"{group}.json.gz"
        rows = _read_gzip_json(shard_path)
        for values in rows:
            species = dict(zip(columns, values))
            for host_id, pair in (species.get("pair_profiles") or {}).items():
                if str(host_id) == human_host_id and species["organism_id"] in human_profiles:
                    human_profile = human_profiles[species["organism_id"]]
                    pair_total = human_profile.get("n_pair_search_records")
                    if pair_total is not None:
                        pair["n_pair_search_records"] = int(pair_total)
                fraction = _sample_fraction(pair)
                if fraction is None:
                    pair["pathogenicity_fraction"] = None
                    pair["pathogenicity_fraction_index"] = None
                    pair["pathogenic_interest_index"] = None
                    pair["estimated_pathogen_publications"] = None
                    pair["estimated_pathogen_publications_index"] = None
                    continue
                pair["pathogenicity_fraction"] = fraction
                fraction_logs.setdefault(group, []).append(math.log1p(fraction))
                restored_by_group[group] = restored_by_group.get(group, 0) + 1
                pair_total = pair.get("n_pair_search_records")
                if pair_total is not None and int(pair_total) > 0:
                    estimate = fraction * int(pair_total)
                    pair["estimated_pathogen_publications"] = estimate
                    volume_logs.setdefault(group, []).append(math.log1p(estimate))
                else:
                    pair["estimated_pathogen_publications"] = None
                    pair["estimated_pathogen_publications_index"] = None
                pair["index_version"] = INDEX_VERSION
        by_group[group] = rows

    fraction_max = {group: max(values) for group, values in fraction_logs.items()}
    volume_max = {group: max(values) for group, values in volume_logs.items()}
    versions = {}
    for group, rows in by_group.items():
        for values in rows:
            species = dict(zip(columns, values))
            for pair in (species.get("pair_profiles") or {}).values():
                fraction = pair.get("pathogenicity_fraction")
                maximum = fraction_max.get(group, 0.0)
                if fraction is not None and maximum > 0:
                    score = math.log1p(float(fraction)) / maximum
                    pair["pathogenicity_fraction_index"] = score
                    pair["pathogenic_interest_index"] = score
                estimate = pair.get("estimated_pathogen_publications")
                maximum = volume_max.get(group, 0.0)
                if estimate is not None and maximum > 0:
                    pair["estimated_pathogen_publications_index"] = math.log1p(float(estimate)) / maximum
        encoded = json.dumps(rows, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode()
        versions[group] = hashlib.sha256(encoded).hexdigest()[:16]
        _write_gzip_json(site_root / "taxa" / f"{group}.json.gz", rows)
    data["taxa_shard_versions"] = versions
    _write_gzip_json(data_path, data)
    return {"source_profiles_restored": source_restored, "restored_by_group": restored_by_group,
            "taxa_shard_versions": versions}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site-root", type=Path, default=DEFAULT_SITE)
    args = parser.parse_args()
    print(json.dumps(restore(args.site_root), indent=2))


if __name__ == "__main__":
    main()
