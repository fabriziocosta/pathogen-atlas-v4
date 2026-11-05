"""Integrity checks for the static Explorer release payload."""

from __future__ import annotations

import gzip
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any


REQUIRED_FILES = ("index.html", "explorer-data.json.gz", "evidence.json.gz", "feature-catalog.json")


def _read_json(path: Path, *, compressed: bool = False) -> Any:
    raw = path.read_bytes()
    if compressed:
        raw = gzip.decompress(raw)
    return json.loads(raw)


def _finite_number(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError, OverflowError):
        return False


def validate_explorer_release(site: Path, *, expected_version: str | None = None) -> dict[str, int]:
    """Validate all data files that the Explorer page loads before/after deployment."""
    site = Path(site)
    errors: list[str] = []
    missing = [name for name in REQUIRED_FILES if not (site / name).is_file()]
    if missing:
        raise ValueError(f"Explorer release is missing required files: {', '.join(missing)}")

    html = (site / "index.html").read_text(encoding="utf-8")
    if expected_version:
        if f'<span class="title-version">v{expected_version}</span>' not in html:
            errors.append(f"index.html does not display expected release version {expected_version}")
    for asset in ("explorer-data.json.gz", "feature-catalog.json"):
        if asset not in html:
            errors.append(f"index.html does not reference {asset}")
    if ("const separator=path.includes('?')?'&':'?'" not in html
            or "DATA.evidence_version" not in html
            or "fetch(url,{cache:'no-store'})" not in html
            or "fetch('explorer-data.json.gz',{cache:'no-store'})" not in html):
        errors.append("index.html does not fetch release metadata and evidence with cache-busting keys and a no-store policy")

    try:
        metadata = _read_json(site / "explorer-data.json.gz", compressed=True)
        catalog = _read_json(site / "feature-catalog.json")
        evidence = _read_json(site / "evidence.json.gz", compressed=True)
    except (OSError, EOFError, gzip.BadGzipFile, json.JSONDecodeError) as exc:
        raise ValueError(f"Explorer core payload could not be read: {exc}") from exc

    if not isinstance(metadata, dict) or not isinstance(evidence, list) or not isinstance(catalog, dict):
        raise ValueError("Explorer core payload has an unexpected JSON shape")
    groups = metadata.get("groups")
    columns = metadata.get("species_columns")
    hosts = metadata.get("hosts")
    if not isinstance(groups, list) or not groups or not isinstance(columns, list) or "pair_profiles" not in columns:
        errors.append("explorer-data.json.gz is missing groups or pair_profiles columns")
        groups, columns = [], []
    if not isinstance(hosts, list) or not hosts:
        errors.append("explorer-data.json.gz contains no hosts")

    feature_ids = set((catalog.get("features") or {}).keys()) if isinstance(catalog.get("features"), dict) else set()
    if not feature_ids:
        errors.append("feature-catalog.json contains no feature definitions")
    for dimension in catalog.get("feature_hierarchy", []):
        for category in dimension.get("categories", []):
            unknown = set(category.get("feature_ids", [])) - feature_ids
            if unknown:
                errors.append(f"feature hierarchy references undefined feature IDs: {sorted(unknown)[:5]}")
                break

    accepted_docs: dict[tuple[str, str], list[dict[str, Any]]] = {}
    doc_keys: set[tuple[str, str, str]] = set()
    for i, row in enumerate(evidence):
        if not isinstance(row, dict):
            errors.append(f"evidence record {i} is not an object")
            continue
        pmid = str(row.get("pmid") or "").strip()
        organism_id = str(row.get("organism_id") or row.get("pathogen_id") or "").strip()
        host_id = str(row.get("host_id") or "").strip()
        if not pmid.isdigit():
            errors.append(f"evidence record {i} has no numeric PMID")
        if not organism_id or not host_id:
            errors.append(f"evidence record {i} has no pathogen or host ID")
        key = (organism_id, host_id, pmid)
        if key in doc_keys:
            errors.append(f"duplicate evidence document for {organism_id}, {host_id}, PMID {pmid}")
        doc_keys.add(key)
        if row.get("accepted_pathogen_evidence") and row.get("evidence_kind") != "literature_feature":
            accepted_docs.setdefault((organism_id, host_id), []).append(row)
        unknown_yes = set(map(str, row.get("feature_ids_yes") or [])) - feature_ids
        if unknown_yes:
            errors.append(f"evidence record {i} references unknown feature IDs: {sorted(unknown_yes)[:5]}")

    evidence_digest = hashlib.sha256(
        json.dumps(evidence, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()[:16]
    if metadata.get("evidence_version") and metadata["evidence_version"] != evidence_digest:
        errors.append("evidence.json.gz does not match explorer-data.json.gz evidence_version")

    if not groups or not columns:
        groups = []
    else:
        pair_index = columns.index("pair_profiles")
        id_index = columns.index("organism_id") if "organism_id" in columns else 0
        count_by_group = (metadata.get("stats") or {}).get("pathogen_species_by_group") or {}
        coordinate_hosts = set(map(str, metadata.get("coordinate_host_ids") or []))
        coordinates_by_host: Counter[str] = Counter()
        seen_species: set[str] = set()
        total_species = total_profiles = 0
        for group in groups:
            shard_path = site / "taxa" / f"{group}.json.gz"
            if not shard_path.is_file():
                errors.append(f"missing taxonomy shard taxa/{group}.json.gz")
                continue
            try:
                rows = _read_json(shard_path, compressed=True)
            except (OSError, EOFError, gzip.BadGzipFile, json.JSONDecodeError) as exc:
                errors.append(f"taxonomy shard {group} could not be read: {exc}")
                continue
            if not isinstance(rows, list):
                errors.append(f"taxonomy shard {group} is not a list")
                continue
            expected_rows = count_by_group.get(group)
            if expected_rows is not None and len(rows) != int(expected_rows):
                errors.append(f"taxonomy shard {group} has {len(rows)} species; metadata reports {expected_rows}")
            shard_digest = hashlib.sha256(
                json.dumps(rows, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode()
            ).hexdigest()[:16]
            declared_digest = (metadata.get("taxa_shard_versions") or {}).get(group)
            if declared_digest and declared_digest != shard_digest:
                errors.append(f"taxonomy shard {group} does not match its declared version")
            total_species += len(rows)
            for row_number, row in enumerate(rows):
                if not isinstance(row, list) or len(row) != len(columns):
                    errors.append(f"taxonomy shard {group} row {row_number} does not match species_columns")
                    continue
                organism_id = str(row[id_index] or "")
                if not organism_id or organism_id in seen_species:
                    errors.append(f"taxonomy shard {group} has an empty or duplicate organism ID: {organism_id}")
                seen_species.add(organism_id)
                profiles = row[pair_index] or {}
                if not isinstance(profiles, dict):
                    errors.append(f"{organism_id} pair_profiles is not an object")
                    continue
                for host_id, profile in profiles.items():
                    total_profiles += 1
                    host_id = str(host_id)
                    if not isinstance(profile, dict):
                        errors.append(f"{organism_id}/{host_id} profile is not an object")
                        continue
                    evidence_file = profile.get("evidence_file")
                    evidence_path = str(evidence_file or "").split("?", 1)[0]
                    if evidence_path and not (site / evidence_path).is_file():
                        errors.append(f"{organism_id}/{host_id} references missing document file {evidence_file}")
                    x, y = profile.get("x"), profile.get("y")
                    if (x is None) != (y is None):
                        errors.append(f"{organism_id}/{host_id} has only one 2D coordinate")
                    elif x is not None:
                        if not (_finite_number(x) and _finite_number(y)):
                            errors.append(f"{organism_id}/{host_id} has non-finite 2D coordinates")
                        else:
                            coordinates_by_host[host_id] += 1
                    for embedding, point in (profile.get("embedding_coordinates") or {}).items():
                        if not isinstance(point, dict) or not (_finite_number(point.get("x")) and _finite_number(point.get("y"))):
                            errors.append(f"{organism_id}/{host_id} has invalid {embedding} coordinates")
                    expected_pmids = int(profile.get("n_jev_accepted_pmids") or profile.get("n_positive_records") or 0)
                    docs = accepted_docs.get((organism_id, host_id), [])
                    if expected_pmids > 0 and not docs:
                        errors.append(f"{organism_id}/{host_id} reports {expected_pmids} accepted PMIDs but has no documents")

        if metadata.get("coordinates_ready") is not True:
            errors.append("Explorer metadata says 2D coordinates are not ready")
        for host_id in sorted(coordinate_hosts):
            if not coordinates_by_host[host_id]:
                errors.append(f"coordinate host {host_id} has no finite 2D points in taxonomy shards")
        if total_species and not seen_species:
            errors.append("taxonomy shards contain no species IDs")

    if errors:
        preview = "\n - ".join(errors[:40])
        remaining = f"\n - ... and {len(errors) - 40} more" if len(errors) > 40 else ""
        raise ValueError(f"Explorer release validation failed ({len(errors)} issue(s)):\n - {preview}{remaining}")

    return {
        "species": total_species,
        "pair_profiles": total_profiles,
        "evidence_documents": len(evidence),
        "coordinate_hosts": len(coordinates_by_host),
        "features": len(feature_ids),
    }
