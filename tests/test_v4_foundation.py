from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from pathogen_atlas_v4.decision_identity import decision_cache_key
from pathogen_atlas_v4.decision_runner import run_decisions
from pathogen_atlas_v4.evidence import normalize_pair_gate
from pathogen_atlas_v4.features import resolve_profile
from pathogen_atlas_v4.host_taxonomy.ncbi import load_ncbi_hosts
from pathogen_atlas_v4.pubmed import build_pair_query
from pathogen_atlas_v4.release import validate_pair_release
from pathogen_atlas_v4.runs import resolve_run, save_resolved_run

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
