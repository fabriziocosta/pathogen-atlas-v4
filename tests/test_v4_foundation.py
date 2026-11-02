from __future__ import annotations

import json
import gzip
import importlib.util
from pathlib import Path

import pandas as pd
import pytest

from pathogen_atlas_v4.decision_identity import decision_cache_key
from pathogen_atlas_v4.decision_runner import run_decisions
from pathogen_atlas_v4.evidence import normalize_pair_gate
from pathogen_atlas_v4.features import resolve_profile
from pathogen_atlas_v4.host_taxonomy.ncbi import load_ncbi_hosts
from pathogen_atlas_v4.migration.import_v3_data import build_manifest, convert_snapshot, verify_manifest
from pathogen_atlas_v4.pubmed import build_pair_query
from pathogen_atlas_v4.release import validate_pair_release
from pathogen_atlas_v4.runs import resolve_run, save_resolved_run

_SITE_BUILDER_PATH = Path(__file__).parents[1] / "sites/pathogen-atlas-explorer/scripts/build_v4_site_data.py"
_SITE_BUILDER_SPEC = importlib.util.spec_from_file_location("v4_site_data_builder", _SITE_BUILDER_PATH)
assert _SITE_BUILDER_SPEC and _SITE_BUILDER_SPEC.loader
_SITE_BUILDER = importlib.util.module_from_spec(_SITE_BUILDER_SPEC)
_SITE_BUILDER_SPEC.loader.exec_module(_SITE_BUILDER)


def test_cache_keys_follow_feature_scope():
    common = dict(pathogen_id="p1", input_hash="input", feature_id="f1",
                  question_version="q1", profile_version="pr1", model_version="m1")
    pathogen = decision_cache_key(scope="pathogen", **common)
    assert pathogen == decision_cache_key(scope="pathogen", host_id="h2", **common)
    assert decision_cache_key(scope="pathogen_host_pair", host_id="h1", **common) != decision_cache_key(
        scope="pathogen_host_pair", host_id="h2", **common)
    with pytest.raises(ValueError):
        decision_cache_key(scope="pathogen_host_pair", **common)


def test_run_pair_generation_is_explicit_and_validated():
    run = resolve_run({"run": {"run_id": "r1", "config_version": "v1", "host_ids": ["h1", "h2"],
                                "pathogen_ids": ["p1"], "pair_generation": "cartesian_product"}},
                      known_host_ids={"h1", "h2"}, known_pathogen_ids={"p1"})
    assert run.pairs == (("p1", "h1"), ("p1", "h2"))
    with pytest.raises(ValueError):
        resolve_run({"run": {"run_id": "r1", "config_version": "v1", "host_ids": ["bad"],
                              "pathogen_ids": ["p1"], "pair_generation": "cartesian_product"}},
                    known_host_ids={"h1"}, known_pathogen_ids={"p1"})


def test_run_groups_resolve_and_persist_exact_config(tmp_path: Path):
    config = {"run": {"run_id": "r2", "config_version": "v1", "host_groups": ["crop"],
                      "pathogen_groups": ["fungi"], "pair_generation": "cartesian_product"}}
    resolved = resolve_run(config, known_host_ids={"h1", "h2"}, known_pathogen_ids={"p1", "p2"},
                           host_group_ids={"crop": ["h1", "h2"]},
                           pathogen_group_ids={"fungi": ["p1", "p2"]})
    assert len(resolved.pairs) == 4
    path = tmp_path / "resolved.json"
    save_resolved_run(config, resolved, path)
    saved = json.loads(path.read_text())
    assert saved["resolved_config_hash"] == resolved.resolved_config_hash
    assert len(saved["pairs"]) == 4


def test_pubmed_pair_query_is_host_aware_and_versioned():
    query = build_pair_query(["Pathogen alpha", "P. alpha"], ["Oryza sativa", "rice"],
                             strategy_version="pair-query-v1")
    assert "Pathogen alpha" in query["query"]
    assert "Oryza sativa" in query["query"]
    assert query["query_version"] == "pair-query-v1"
    assert query["query_hash"] == build_pair_query(["Pathogen alpha", "P. alpha"], ["Oryza sativa", "rice"],
                                                   strategy_version="pair-query-v1")["query_hash"]


def test_decision_runner_resumes_only_matching_pair_cache(tmp_path: Path):
    base = {"scope": "pathogen_host_pair", "pathogen_id": "p1", "host_id": "h1",
            "relationship_context": "disease", "pmid": "1", "input_hash": "i1", "feature_id": "gate",
            "query_version": "q1", "prompt_version": "pr1", "feature_version": "f1",
            "question_version": "qv1", "profile_version": "pv1", "model_version": "m1",
            "run_id": "r1", "run_config_version": "rc1", "title_hash": "th", "abstract_hash": "ah"}
    ledger, raw = tmp_path / "ledger.parquet", tmp_path / "raw.jsonl"
    calls = []
    def handler(task):
        calls.append(task["host_id"])
        return {"outcome": "supported", "confidence": .93, "evidence_span": "Disease in this host."}
    first = run_decisions([base], handler, ledger_path=ledger, raw_response_path=raw)
    assert first.iloc[0]["status"] == "completed"
    resumed = run_decisions([base], handler, ledger_path=ledger, raw_response_path=raw)
    assert calls == ["h1"]
    second_host = {**base, "host_id": "h2"}
    run_decisions([second_host], handler, ledger_path=ledger, raw_response_path=raw)
    assert calls == ["h1", "h2"]
    assert pd.read_parquet(ledger).cache_key.nunique() == 2


def test_supported_gate_requires_passage_and_release_validates_links():
    with pytest.raises(ValueError, match="evidence span"):
        normalize_pair_gate(pathogen_id="p", host_id="h", pmid="1", outcome="supported",
                            confidence=.9, evidence_text=None, title_hash="t", abstract_hash="a", run_id="r")
    record = normalize_pair_gate(pathogen_id="p", host_id="h", pmid="1", outcome="supported",
                                 confidence=.9, evidence_text="The pathogen caused disease in host.",
                                 title_hash="t", abstract_hash="a", run_id="r")
    validate_pair_release(pd.DataFrame([{"pathogen_id": "p"}]), pd.DataFrame([{"host_id": "h"}]),
                          pd.DataFrame([{"pmid": "1"}]), pd.DataFrame([record]))


def test_profile_resolves_ordered_modules():
    catalogue = {"profile_id": "test", "modules": {
        "core": [{"feature_id": "a", "scope": "pathogen"}],
        "pathogen_group:fungi": [{"feature_id": "b", "scope": "host"}],
        "host_group:crop": [{"feature_id": "c", "scope": "pathogen_host_pair"}],
        "host_taxon:4530": [{"feature_id": "d", "scope": "host"}],
        "pair_context:disease": [{"feature_id": "e", "scope": "pathogen_host_pair"}],
    }}
    result = resolve_profile(catalogue, pathogen_group="fungi", host_group="crop", host_taxon_id="4530",
                             relationship_context="disease", profile_version="v1")
    assert [x["feature_id"] for x in result["ordered_features"]] == ["a", "b", "c", "d", "e"]


def test_ncbi_adapter_resolves_names_and_nodes(tmp_path: Path):
    (tmp_path / "names.dmp").write_text("9606\t| Homo sapiens\t|\t| scientific name\t|\n9606\t| human\t|\t| genbank common name\t|\n")
    (tmp_path / "nodes.dmp").write_text("9606\t| 9605\t| species\t|\n")
    hosts = load_ncbi_hosts(tmp_path, ["9606"], release="fixture")
    assert hosts.iloc[0]["host_id"] == "ncbi:9606"
    assert hosts.iloc[0]["aliases"] == ["human"]


def _fixture_snapshot(root: Path) -> Path:
    def write(relative: str, rows: list[dict]):
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_parquet(path, index=False)
    write("01_taxonomy/release/species.parquet", [{"organism_id": "bacteria:id1", "organism_group": "bacteria",
          "accepted_name": "Bugus testus", "source_taxon_id": "10", "mapping_status": "accepted"}])
    write("01_taxonomy/release/aliases.parquet", [{"organism_id": "bacteria:id1", "organism_group": "bacteria", "name": "Bugus testus"}])
    write("02_pubmed/release/abstracts.parquet", [
        {"pmid": "1", "title": "Title", "abstract": "Disease reported."},
        {"pmid": "1", "title": "Title", "abstract": "Disease reported."},
        {"pmid": "2", "title": "Other", "abstract": "Second source."},
        {"pmid": "2", "title": "Conflict", "abstract": "Different source."},
    ])
    write("02_pubmed/release/species_abstracts.parquet", [{"organism_id": "bacteria:id1", "organism_group": "bacteria",
          "pmid": "1", "selected_for_fetch": True, "query_hash": "q"},
          {"organism_id": "bacteria:id1", "organism_group": "bacteria", "pmid": "3",
           "selected_for_fetch": False, "query_hash": "q2"}])
    write("03_decisions/parsed/decision_ledger.parquet", [{"organism_id": "bacteria:id1", "organism_group": "bacteria",
          "pmid": "1", "decision_kind": "pathogen_evidence", "p_yes": .9, "status": "completed"},
          {"organism_id": "bacteria:id1", "organism_group": "bacteria", "pmid": "1",
           "decision_kind": "feature_assessment", "feature_id": "f1", "p_yes": .8, "status": "completed"}])
    write("04_release/accepted_pathogen_evidence.parquet", [{"organism_id": "bacteria:id1", "organism_group": "bacteria",
          "pmid": "1", "accepted_pathogen_evidence": True, "p_yes": .9, "question_version": "q1", "input_hash": "i1"}])
    cache = root / "02_pubmed/cache/search.jsonl"
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps({"query_hash": "hash-q", "query_id": "query-1", "pmids": ["1", "3"]}) + "\n")
    return root


def test_snapshot_manifest_and_converter_account_legacy_evidence(tmp_path: Path):
    source = _fixture_snapshot(tmp_path / "source")
    manifest = build_manifest(source, "commit")
    copied = tmp_path / "copy"
    import shutil
    shutil.copytree(source, copied)
    assert verify_manifest(manifest, copied)["ok"]
    output = tmp_path / "staging"
    report = convert_snapshot(copied, output)
    assert report["duplicate_counts"]["articles"] == 2
    assert len(report["article_metadata_conflicts"]) == 1
    articles = pd.read_parquet(output / "literature/articles.parquet")
    assert len(articles) == 3 and articles.loc[articles.pmid.eq("2"), "metadata_conflict"].item()
    # Every search match gets an article identity; unmatched abstracts are explicit stubs.
    assert articles.loc[articles.pmid.eq("3"), "abstract_fetch_status"].item() == "not_fetched_in_v3"
    matches = pd.read_parquet(output / "literature/pathogen_article_matches.parquet")
    assert set(matches.pmid) == {"1", "3"}
    assert set(matches.pmid) <= set(articles.pmid)
    evidence = pd.read_parquet(output / "evidence/pathogen_host_evidence.parquet")
    assert evidence.iloc[0]["host_id"] == "ncbi:9606"
    assert evidence.iloc[0]["provenance_status"] == "legacy_v3_accepted"
    assert evidence.iloc[0]["evidence_span_type"] == "full_abstract_legacy"
    legacy = pd.read_parquet(output / "decisions/legacy_v3_decisions.parquet")
    assert not legacy["native_v4_decision"].any()
    assert report["unmapped_counts"]["accepted_legacy_evidence"] == 0
    source_records = pd.read_parquet(output / "provenance/source_records.parquet")
    indexed = source_records.loc[source_records.query_hash.eq("hash-q")]
    assert set(indexed.pmid.dropna()) == {"1", "3"}
    second_output = tmp_path / "staging_repeat"
    second_report = convert_snapshot(copied, second_output)
    assert report == second_report
    pd.testing.assert_frame_equal(
        pd.read_parquet(output / "evidence/pathogen_host_evidence.parquet"),
        pd.read_parquet(second_output / "evidence/pathogen_host_evidence.parquet"),
    )


def test_manifest_detects_changes(tmp_path: Path):
    source = tmp_path / "source"
    source.mkdir()
    (source / "one.txt").write_text("original")
    manifest = build_manifest(source)
    (source / "one.txt").write_text("changed")
    assert not verify_manifest(manifest, source)["ok"]


def test_human_explorer_import_keeps_legacy_parity_without_copying_abstracts(tmp_path: Path):
    staging, export, output = tmp_path / "staging", tmp_path / "legacy", tmp_path / "site"
    for relative, rows in {
        "taxonomy/pathogens.parquet": [{"pathogen_id": "bacteria:p1"}],
        "taxonomy/hosts.parquet": [
            {"host_id": "ncbi:9606", "accepted_name": "Homo sapiens"},
            {"host_id": "ncbi:9913", "accepted_name": "Bos taurus"},
            {"host_id": "ncbi:4577", "accepted_name": "Zea mays"},
        ],
    }.items():
        path = staging / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_parquet(path, index=False)

    (export / "evidence").mkdir(parents=True)
    profile = {
        "organism_id": "bacteria:p1", "organism_group": "bacteria", "accepted_name": "Pathogen alpha",
        "n_pathogen_pmids": 1, "n_jev_accepted_pmids": 1, "n_pubmed_search_records": 3,
        "n_abstracts_fetched": 1, "pathogenic_interest_index": 0.42, "x": 1.0, "y": 2.0,
        "feature_values": {"f1": 0.5}, "feature_counts": {"f1": {"n_yes": 1, "n_scored": 2}},
        "evidence_file": "evidence/p1.json", "taxonomy": {"genus": "Pathogen"},
    }
    (export / "evidence/p1.json").write_text(json.dumps([{
        "organism_id": "bacteria:p1", "pmid": "123", "title": "A human disease report",
        "accepted_pathogen_evidence": True, "feature_ids_yes": ["f1"],
        "abstract": "Do not copy this abstract.",
    }]))
    legacy = {
        "release": "Pathogen Atlas V3", "species": [profile], "groups": ["bacteria"],
        "features": {"f1": {"label": "Feature 1"}},
        "feature_hierarchy": [{"id": "d1", "label": "Dimension", "categories": [{"id": "c1", "label": "Category", "feature_ids": ["f1"]}]}],
        "feature_profiles": {"bacteria": ["f1"]}, "coordinates_ready": True,
        "stats": {"accepted_species_pmid_records": 1, "feature_question_set_version": "legacy-profile-v1"},
    }
    with gzip.open(export / "explorer-data.json.gz", "wt", encoding="utf-8") as stream:
        json.dump(legacy, stream)

    result = _SITE_BUILDER.build_site_data(staging, output, export)
    payload = json.loads(gzip.decompress((output / "explorer-data.json.gz").read_bytes()))
    evidence = json.loads(gzip.decompress((output / "evidence.json.gz").read_bytes()))
    assert result["profile_species"] == 1 and result["accepted_species_pmid_records"] == 1
    assert payload["species"][0]["host_id"] == "ncbi:9606"
    assert payload["species"][0]["profile_provenance_status"] == "legacy_import"
    assert payload["species"][0]["feature_values"] == {"f1": 0.5}
    assert evidence[0]["title"] == "A human disease report"
    assert evidence[0]["host_id"] == "ncbi:9606"
    assert "abstract" not in evidence[0] and "evidence_span" not in evidence[0]
    assert payload["hosts"][0]["host_group"] == "Animals"
    assert {host["host_group"] for host in payload["hosts"]} == {"Animals", "Plants"}
