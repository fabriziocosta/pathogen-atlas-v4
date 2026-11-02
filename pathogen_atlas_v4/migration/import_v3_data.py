"""Checksum, snapshot, and deterministically convert a V3 data tree to V4 staging."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow.parquet as pq

from pathogen_atlas_v4.decision_identity import stable_hash
from pathogen_atlas_v4.host_taxonomy.ncbi import load_ncbi_hosts
from pathogen_atlas_v4.schemas import PAIR_RELATIONSHIP

SOURCE_COMMIT = "59409a9136bb1ce8477bc80179c853a9132188a8"
MIGRATION_CONFIG_VERSION = "v3-to-v4-import-v5-host-cohorts-ncbi-f278697a8ba6"
HOST_TAXON_IDS = ["9606", "9913", "9940", "9031", "9823", "4577", "4530", "4565", "3847"]
HUMAN_HOST = {
    "host_id": "ncbi:9606", "accepted_name": "Homo sapiens", "aliases": ["human"],
    "host_group": "human", "taxon_id": "9606", "parent_taxon_id": "9605",
    "rank": "species", "taxonomy_source": "ncbi_taxonomy",
}


def sha256_file(path: Path) -> tuple[int, str]:
    digest, size = hashlib.sha256(), 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            size += len(block)
            digest.update(block)
    return size, digest.hexdigest()


def build_manifest(data_root: Path, source_commit: str = SOURCE_COMMIT) -> dict[str, Any]:
    root = data_root.resolve()
    symlinks = sorted(p for p in root.rglob("*") if p.is_symlink())
    if symlinks:
        raise ValueError("Review source symlinks before snapshotting: " + ", ".join(map(str, symlinks[:20])))
    files = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        size, digest = sha256_file(path)
        files.append({"path": path.relative_to(root).as_posix(), "size_bytes": size, "sha256": digest})
    return {"migration_id": "v3-to-v4", "created_utc": datetime.now(timezone.utc).isoformat(),
            "source_repository": "pathogen-atlas-v3",
            "source_commit": source_commit, "source_data_root": "data/", "file_count": len(files),
            "total_bytes": sum(row["size_bytes"] for row in files), "files": files}


def verify_manifest(manifest: dict[str, Any], data_root: Path) -> dict[str, Any]:
    root = data_root.resolve()
    expected = {row["path"]: row for row in manifest["files"]}
    actual = {p.relative_to(root).as_posix(): p for p in root.rglob("*") if p.is_file()}
    errors = []
    for rel, row in expected.items():
        path = actual.get(rel)
        if path is None:
            errors.append(f"missing: {rel}")
            continue
        size, digest = sha256_file(path)
        if size != row["size_bytes"]:
            errors.append(f"size mismatch: {rel}")
        elif digest != row["sha256"]:
            errors.append(f"SHA-256 mismatch: {rel}")
    errors.extend(f"unexpected: {rel}" for rel in sorted(set(actual) - set(expected)))
    return {"ok": not errors, "errors": errors, "verified_files": len(expected)}


def create_snapshot(source: Path, snapshot_dir: Path) -> dict[str, Any]:
    """Copy a frozen V3 data root to an immutable-input location and verify it."""
    manifest = build_manifest(source)
    target = snapshot_dir / "data"
    if target.exists() and any(target.iterdir()):
        raise FileExistsError(f"Snapshot destination is not empty: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(source, target, dirs_exist_ok=True, copy_function=shutil.copy2)
    (snapshot_dir / "transfer_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    result = verify_manifest(manifest, target)
    if not result["ok"]:
        raise RuntimeError("Snapshot verification failed: " + "\n".join(result["errors"][:30]))
    return manifest


def _read_table(path: Path, *, required: bool = False) -> pd.DataFrame:
    if not path.exists():
        if required:
            raise FileNotFoundError(path)
        return pd.DataFrame()
    return pd.read_parquet(path)


def _write(frame: pd.DataFrame, root: Path, relative: str) -> int:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(path, index=False)
    return len(frame)


def _renamed(frame: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    frame = frame.copy()
    return frame.rename(columns={old: new for old, new in mapping.items() if old in frame.columns})


def _articles(frame: pd.DataFrame, conflicts: list[dict[str, Any]]) -> pd.DataFrame:
    if frame.empty:
        return frame
    if "pmid" not in frame:
        raise ValueError("V3 abstracts table has no pmid column")
    frame = frame.copy()
    frame["pmid"] = frame["pmid"].astype("string").str.strip()
    content_cols = [col for col in ("title", "abstract", "publication_date", "publication_year", "doi") if col in frame]
    chosen = []
    for pmid, group in frame.groupby("pmid", sort=True, dropna=False):
        signatures = {stable_hash({col: row.get(col) for col in content_cols}) for row in group.to_dict("records")}
        if len(signatures) > 1:
            conflicts.append({"pmid": str(pmid), "duplicate_rows": len(group), "metadata_signatures": len(signatures)})
        # Select deterministically independent of source row order.
        rows = group.to_dict("records")
        rows.sort(key=lambda row: stable_hash(row))
        selected = rows[0]
        selected["source_row_hash"] = stable_hash(selected)
        selected["metadata_conflict"] = len(signatures) > 1
        selected["abstract_fetch_status"] = "fetched"
        selected["source_version"] = "v3"
        chosen.append(selected)
    return pd.DataFrame(chosen).sort_values("pmid", kind="stable").reset_index(drop=True)


def _legacy_decisions(frame: pd.DataFrame, source_file: str) -> pd.DataFrame:
    if frame.empty:
        return frame.assign(source_file=pd.Series(dtype="string"), provenance_status=pd.Series(dtype="string"))
    frame = frame.copy()
    frame["source_file"] = source_file
    frame["source_project"] = "pathogen-atlas-v3"
    frame["provenance_status"] = "legacy_v3"
    frame["native_v4_decision"] = False
    if "organism_id" in frame:
        frame = frame.rename(columns={"organism_id": "pathogen_id", "organism_group": "pathogen_group"})
    return frame


def convert_snapshot(snapshot_data: Path, output: Path, *, migration_config_version: str = MIGRATION_CONFIG_VERSION,
                     dry_run: bool = False, host_taxdump: Path | None = None,
                     host_taxonomy_release: str | None = None,
                     source_manifest: dict[str, Any] | None = None) -> dict[str, Any]:
    """Convert canonical V3 tables while keeping original source data untouched."""
    if not snapshot_data.is_dir():
        raise FileNotFoundError(snapshot_data)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Output staging directory is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    report: dict[str, Any] = {"migration_id": "v3-to-v4", "migration_config_version": migration_config_version,
                              "source_commit": SOURCE_COMMIT, "source_data": str(snapshot_data),
                              "conversion_id": stable_hash({"source_commit": SOURCE_COMMIT, "config": migration_config_version})[:20],
                              "tables": {}, "article_metadata_conflicts": [], "warnings": [], "rejected_rows": [],
                              "duplicate_counts": {}, "unmapped_counts": {}}
    base = snapshot_data

    species = _read_table(base / "01_taxonomy/release/species.parquet", required=True)
    species = _renamed(species, {"organism_id": "pathogen_id", "organism_group": "pathogen_group",
                                  "source_taxon_id": "taxon_id"})
    species["source_version"] = "v3"
    report["tables"]["taxonomy/pathogens.parquet"] = {"source_rows": len(species), "migrated_rows": len(species)}
    report["taxonomy_mapping_status_counts"] = (
        species["mapping_status"].fillna("missing").astype(str).value_counts().sort_index().to_dict()
        if "mapping_status" in species else {"missing_status_column": len(species)})
    if not dry_run:
        _write(species, output, "taxonomy/pathogens.parquet")
    pathogen_ids = set(species["pathogen_id"].dropna().astype(str))

    aliases = _read_table(base / "01_taxonomy/release/aliases.parquet", required=True)
    aliases = _renamed(aliases, {"organism_id": "pathogen_id", "organism_group": "pathogen_group"})
    bad_alias = aliases["pathogen_id"].notna() & ~aliases["pathogen_id"].astype(str).isin(pathogen_ids)
    report["unmapped_counts"]["pathogen_aliases"] = int(bad_alias.sum())
    report["alias_mapping_status_counts"] = (
        aliases["mapping_status"].fillna("missing").astype(str).value_counts().sort_index().to_dict()
        if "mapping_status" in aliases else {"missing_status_column": len(aliases)})
    report["tables"]["taxonomy/pathogen_aliases.parquet"] = {"source_rows": len(aliases), "migrated_rows": len(aliases) - int(bad_alias.sum()), "unmapped_rows": int(bad_alias.sum())}
    if not dry_run:
        _write(aliases, output, "taxonomy/pathogen_aliases.parquet")

    variants = _read_table(base / "01_taxonomy/release/virus_variants.parquet")
    if not variants.empty:
        variants = _renamed(variants, {"organism_id": "pathogen_id", "organism_group": "pathogen_group"})
    report["tables"]["taxonomy/pathogen_variants.parquet"] = {"source_rows": len(variants), "migrated_rows": len(variants)}
    if not dry_run and not variants.empty:
        _write(variants, output, "taxonomy/pathogen_variants.parquet")

    articles_source = _read_table(base / "02_pubmed/release/abstracts.parquet", required=True)
    articles = _articles(articles_source, report["article_metadata_conflicts"])
    report["duplicate_counts"]["articles"] = len(articles_source) - len(articles)
    report["tables"]["literature/articles.parquet"] = {"source_rows": len(articles_source), "migrated_rows": len(articles), "duplicate_rows": len(articles_source) - len(articles), "conflicting_pmids": len(report["article_metadata_conflicts"])}

    matches = _read_table(base / "02_pubmed/release/species_abstracts.parquet", required=True)
    matches = _renamed(matches, {"organism_id": "pathogen_id", "organism_group": "pathogen_group"})
    article_ids_before_placeholders = set(articles["pmid"].dropna().astype(str))
    missing_pmids = sorted(set(matches["pmid"].dropna().astype(str)) - article_ids_before_placeholders)
    if missing_pmids:
        placeholders = pd.DataFrame({"pmid": missing_pmids, "title": None, "abstract": None,
                                    "source_row_hash": None, "metadata_conflict": False,
                                    "abstract_fetch_status": "not_fetched_in_v3", "source_version": "v3_search_match"})
        articles = pd.concat([articles, placeholders], ignore_index=True, sort=False).sort_values("pmid", kind="stable").reset_index(drop=True)
    report["tables"]["literature/articles.parquet"].update({
        "migrated_rows": len(articles), "source_abstract_rows": len(articles_source),
        "search_match_pmids_without_fetched_abstract": len(missing_pmids),
        "placeholder_article_rows": len(missing_pmids), "output_rows_including_placeholders": len(articles),
    })
    article_ids = set(articles["pmid"].dropna().astype(str)) if "pmid" in articles else set()
    if not dry_run:
        _write(articles, output, "literature/articles.parquet")
    invalid_pathogen = ~matches["pathogen_id"].astype(str).isin(pathogen_ids)
    invalid_article = ~matches["pmid"].astype(str).isin(article_ids)
    duplicate_match_rows = int(matches.duplicated(["pathogen_id", "pmid"], keep="first").sum())
    report["duplicate_counts"]["pathogen_article_match_keys"] = duplicate_match_rows
    matches["source_version"] = "v3"
    report["unmapped_counts"]["literature_match_pathogen"] = int(invalid_pathogen.sum())
    report["unmapped_counts"]["literature_match_article"] = int(invalid_article.sum())
    report["tables"]["literature/pathogen_article_matches.parquet"] = {"source_rows": len(matches), "migrated_rows": int((~invalid_pathogen & ~invalid_article).sum()), "unmapped_pathogen_rows": int(invalid_pathogen.sum()), "unmapped_article_rows": int(invalid_article.sum())}
    valid_matches = matches.loc[~invalid_pathogen & ~invalid_article].copy()
    rejected_matches = matches.loc[invalid_pathogen | invalid_article].copy()
    if not rejected_matches.empty:
        rejected_matches["migration_rejection_reason"] = [
            "unknown_pathogen_and_article" if p and a else "unknown_pathogen" if p else "unknown_article"
            for p, a in zip(invalid_pathogen.loc[rejected_matches.index], invalid_article.loc[rejected_matches.index])
        ]
    report["tables"]["literature/pathogen_article_matches.parquet"]["migrated_rows"] = len(valid_matches)
    report["tables"]["literature/pathogen_article_matches.parquet"]["rejected_rows"] = len(rejected_matches)
    if not dry_run:
        _write(valid_matches, output, "literature/pathogen_article_matches.parquet")
        if not rejected_matches.empty:
            _write(rejected_matches, output, "literature/unresolved_pathogen_article_matches.parquet")

    if host_taxdump:
        host = load_ncbi_hosts(host_taxdump, HOST_TAXON_IDS,
                               release=host_taxonomy_release or "NCBI Taxonomy taxdump (release unspecified)")
    else:
        host = pd.DataFrame([{**HUMAN_HOST, "taxonomy_release": "configured NCBI Taxonomy taxid 9606", "mapping_status": "explicit_legacy_mapping"}])
    report["tables"]["taxonomy/hosts.parquet"] = {"source_rows": len(host), "migrated_rows": len(host)}
    if not dry_run:
        _write(host, output, "taxonomy/hosts.parquet")

    ledger = _read_table(base / "03_decisions/parsed/decision_ledger.parquet", required=True)
    gate = ledger.loc[ledger.get("decision_kind", pd.Series(index=ledger.index, dtype="string")).astype(str).eq("pathogen_evidence")].copy()
    features = ledger.loc[ledger.get("decision_kind", pd.Series(index=ledger.index, dtype="string")).astype(str).ne("pathogen_evidence")].copy()
    gate = _legacy_decisions(gate, "03_decisions/parsed/decision_ledger.parquet")
    features = _legacy_decisions(features, "03_decisions/parsed/decision_ledger.parquet")
    v2 = _read_table(base / "03_decisions/legacy_v2/decisions.parquet")
    if not v2.empty:
        v2 = _legacy_decisions(v2, "03_decisions/legacy_v2/decisions.parquet")
    all_legacy = pd.concat([gate, v2], ignore_index=True, sort=False)
    report["legacy_gate_status_counts"] = (
        gate["status"].fillna("missing").astype(str).value_counts().sort_index().to_dict()
        if "status" in gate else {"missing_status_column": len(gate)})
    report["legacy_feature_status_counts"] = (
        features["status"].fillna("missing").astype(str).value_counts().sort_index().to_dict()
        if "status" in features else {"missing_status_column": len(features)})
    report["tables"]["decisions/legacy_v3_decisions.parquet"] = {"source_rows": len(gate) + len(v2), "migrated_rows": len(all_legacy), "v3_ledger_rows": len(gate), "v2_legacy_rows": len(v2)}
    report["tables"]["decisions/legacy_v3_feature_assessments.parquet"] = {"source_rows": len(features), "migrated_rows": len(features)}
    if not dry_run:
        _write(all_legacy, output, "decisions/legacy_v3_decisions.parquet")
        _write(features, output, "decisions/legacy_v3_feature_assessments.parquet")

    accepted_path = base / "04_release/accepted_pathogen_evidence.parquet"
    accepted = _read_table(accepted_path)
    if not accepted.empty:
        flag = accepted.get("accepted_pathogen_evidence", pd.Series(False, index=accepted.index)).fillna(False).astype(bool)
        accepted = accepted.loc[flag].copy()
        if "organism_id" in accepted:
            accepted = accepted.rename(columns={"organism_id": "pathogen_id", "organism_group": "pathogen_group"})
        accepted["host_id"] = HUMAN_HOST["host_id"]
        accepted["relationship_context"] = PAIR_RELATIONSHIP
        accepted["decision"] = "supported"
        accepted["provenance_status"] = "legacy_v3_accepted"
        accepted["source_project"] = "pathogen-atlas-v3"
        accepted["source_question_version"] = accepted.get("question_version")
        accepted["source_threshold"] = 0.75
        article_lookup = articles.set_index(articles["pmid"].astype(str), drop=False).to_dict("index") if "pmid" in articles else {}
        accepted["title_hash"] = accepted["pmid"].astype(str).map(
            lambda pmid: stable_hash(article_lookup.get(pmid, {}).get("title")))
        accepted["abstract_hash"] = accepted["pmid"].astype(str).map(
            lambda pmid: stable_hash(article_lookup.get(pmid, {}).get("abstract")))
        # V3 stored model probabilities but no extracted offsets. Preserve the complete
        # source abstract as the legacy evidence passage and label its coarse span type.
        accepted["evidence_span"] = accepted["pmid"].astype(str).map(
            lambda pmid: article_lookup.get(pmid, {}).get("abstract"))
        accepted["evidence_span_type"] = "full_abstract_legacy"
        accepted["source_record_id"] = accepted.get("input_hash")
        evidence = accepted[[c for c in ("pathogen_id", "host_id", "relationship_context", "pmid", "title_hash", "abstract_hash", "evidence_span", "evidence_span_type", "decision", "p_yes", "provenance_status", "source_project", "source_question_version", "source_threshold", "source_record_id") if c in accepted]].rename(columns={"p_yes": "confidence"})
        missing_span = evidence["evidence_span"].fillna("").astype(str).str.strip().eq("")
        bad_evidence = (~evidence["pathogen_id"].astype(str).isin(pathogen_ids)
                        | ~evidence["pmid"].astype(str).isin(article_ids) | missing_span)
        report["unmapped_counts"]["accepted_legacy_evidence"] = int(bad_evidence.sum())
        report["unmapped_counts"]["accepted_legacy_evidence_missing_abstract"] = int(missing_span.sum())
        rejected_evidence = evidence.loc[bad_evidence].copy()
        if not rejected_evidence.empty:
            rejected_evidence["migration_rejection_reason"] = [
                "missing_source_abstract" if missing else "unknown_pathogen_or_article"
                for missing in missing_span.loc[rejected_evidence.index]
            ]
            if not dry_run:
                _write(rejected_evidence, output, "evidence/unresolved_legacy_pathogen_host_evidence.parquet")
        evidence = evidence.loc[~bad_evidence]
    else:
        evidence = pd.DataFrame(columns=["pathogen_id", "host_id", "relationship_context", "pmid", "title_hash", "abstract_hash", "evidence_span", "evidence_span_type", "decision", "confidence", "provenance_status", "source_project", "source_question_version", "source_threshold", "source_record_id"])
    report["tables"]["evidence/pathogen_host_evidence.parquet"] = {"source_rows": len(accepted), "migrated_rows": len(evidence), "unmapped_rows": len(accepted) - len(evidence), "evidence_span_available": bool(evidence.empty or evidence["evidence_span"].fillna("").astype(str).str.strip().ne("").all()), "evidence_span_type": "full_abstract_legacy"}
    if not dry_run:
        _write(evidence, output, "evidence/pathogen_host_evidence.parquet")

    # The immutable snapshot is the archive for cache trees, raw provider responses,
    # old releases, embeddings, and V2 snapshots; index files without rewriting them.
    if source_manifest:
        file_index = [{"source_path": row["path"], "size_bytes": row["size_bytes"],
                       "sha256": row["sha256"], "archive_status": "preserved_in_immutable_snapshot"}
                      for row in source_manifest["files"]]
    else:
        file_index = []
        for path in sorted(p for p in base.rglob("*") if p.is_file()):
            size, digest = sha256_file(path)
            file_index.append({"source_path": path.relative_to(base).as_posix(), "size_bytes": size,
                               "sha256": digest, "archive_status": "preserved_in_immutable_snapshot"})
    record_index = [{**entry, "record_type": "source_file", "line_number": None,
                     "query_hash": None, "query_id": None, "pmid": None, "input_hash": None,
                     "source_record_id": entry["source_path"]} for entry in file_index]
    jsonl_records = jsonl_pmids = malformed_jsonl = 0
    for entry in file_index:
        rel = entry["source_path"]
        if not rel.endswith(".jsonl"):
            continue
        with (base / rel).open(encoding="utf-8") as stream:
            for line_number, line in enumerate(stream, 1):
                if not line.strip():
                    continue
                try:
                    item = json.loads(line)
                except json.JSONDecodeError:
                    malformed_jsonl += 1
                    record_index.append({"source_path": rel, "record_type": "malformed_jsonl_record",
                                         "line_number": line_number, "query_hash": None, "query_id": None,
                                         "pmid": None, "input_hash": None, "source_record_id": f"{rel}:{line_number}"})
                    continue
                jsonl_records += 1
                pmids = item.get("pmids") if isinstance(item, dict) else None
                pmids = [str(value) for value in pmids] if isinstance(pmids, list) else []
                if isinstance(item, dict) and item.get("pmid") is not None:
                    pmids = [str(item["pmid"])]
                jsonl_pmids += len(pmids)
                indexed_pmids: list[str | None] = pmids or [None]
                for pmid in indexed_pmids:
                    record_index.append({"source_path": rel, "record_type": "jsonl_record",
                                         "line_number": line_number, "query_hash": item.get("query_hash") if isinstance(item, dict) else None,
                                         "query_id": item.get("query_id") if isinstance(item, dict) else None,
                                         "pmid": pmid, "input_hash": item.get("input_hash") if isinstance(item, dict) else None,
                                         "source_record_id": str((item.get("pmid") or item.get("query_id") or item.get("segment_id") or f"{rel}:{line_number}") if isinstance(item, dict) else f"{rel}:{line_number}")})
    provenance = pd.DataFrame(record_index)
    report["tables"]["provenance/source_records.parquet"] = {
        "source_files": len(file_index), "indexed_files": len(file_index),
        "jsonl_records_indexed": jsonl_records, "jsonl_pmids_indexed": jsonl_pmids,
        "malformed_jsonl_records": malformed_jsonl,
    }
    archived_rows = []
    converted_sources = {
        "01_taxonomy/release/species.parquet", "01_taxonomy/release/aliases.parquet",
        "01_taxonomy/release/virus_variants.parquet", "02_pubmed/release/abstracts.parquet",
        "02_pubmed/release/species_abstracts.parquet", "03_decisions/parsed/decision_ledger.parquet",
        "03_decisions/legacy_v2/decisions.parquet", "04_release/accepted_pathogen_evidence.parquet",
    }
    for entry in file_index:
        rel = entry["source_path"]
        if not rel.endswith(".parquet") or rel in converted_sources:
            continue
        rows = int(pq.ParquetFile(base / rel).metadata.num_rows)
        archived_rows.append({"source_path": rel, "source_rows": rows,
                              "disposition": "archived_comparison_only"})
    report["archived_parquet_tables"] = archived_rows
    report["archived_parquet_row_count"] = sum(row["source_rows"] for row in archived_rows)
    if not dry_run:
        _write(provenance, output, "provenance/source_records.parquet")
        (output / "migration_report.json").write_text(json.dumps(report, indent=2, default=str) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-data", type=Path, required=True, help="V3 repository data/ directory")
    parser.add_argument("--snapshot-dir", type=Path, required=True, help="V4 data/migration/v3_snapshot directory")
    parser.add_argument("--output", type=Path, required=True, help="Empty V4 staging output directory")
    parser.add_argument("--snapshot", action="store_true", help="Create and verify an immutable source snapshot first")
    parser.add_argument("--dry-run", action="store_true", help="Read and report conversion without writing converted tables")
    parser.add_argument("--config-version", default=MIGRATION_CONFIG_VERSION)
    parser.add_argument("--host-taxdump", type=Path, help="Verified directory containing NCBI names.dmp and nodes.dmp")
    parser.add_argument("--host-taxonomy-release", help="Pinned NCBI taxonomy release/checksum label")
    args = parser.parse_args()
    if args.snapshot:
        manifest = create_snapshot(args.source_data, args.snapshot_dir)
        print(f"Snapshot verified: {manifest['file_count']} files, {manifest['total_bytes']} bytes")
    manifest_path = args.snapshot_dir / "transfer_manifest.json"
    snapshot_data = args.snapshot_dir / "data"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    verification = verify_manifest(manifest, snapshot_data)
    if not verification["ok"]:
        raise SystemExit("Snapshot verification failed: " + "; ".join(verification["errors"][:20]))
    report = convert_snapshot(snapshot_data, args.output, migration_config_version=args.config_version, dry_run=args.dry_run,
                              host_taxdump=args.host_taxdump, host_taxonomy_release=args.host_taxonomy_release,
                              source_manifest=manifest)
    if args.dry_run:
        print(json.dumps(report, indent=2, default=str))
    else:
        print(f"Converted staging {report['conversion_id']} to {args.output}")


if __name__ == "__main__":
    main()
