from pathlib import Path

import pandas as pd
import yaml

from pathogen_atlas_v4.feature_screening import build_feature_screen_tasks


def test_feature_tasks_score_every_resolved_question_for_each_accepted_pair_abstract():
    profiles = [{"resolved_profile": {
        "profile_version": "profile-1",
        "selectors": {"pathogen_group": "viruses", "host_group": "crop", "host_taxon_id": "4530"},
        "ordered_features": [
            {"feature_id": "viral_x", "question": "Does the abstract report X?", "question_version": "q1",
             "scope": "pathogen", "applies_to_label": "All viruses"},
            {"feature_id": "plant_x", "question": "Does the abstract report Y?", "question_version": "q1",
             "scope": "host", "applies_to_label": "All plants"},
            {"feature_id": "rice_x", "question": "Does the abstract report Z?", "question_version": "q1",
             "scope": "pathogen_host_pair", "applies_to_label": "Rice viruses"},
        ],
    }}]
    matches = pd.DataFrame([{"pathogen_id": "v1", "host_id": "ncbi:4530", "pmid": "100",
                             "query_version": "query-1"}])
    articles = pd.DataFrame([{"pmid": "100", "title": "Title", "abstract": "An abstract."}])
    selection = pd.DataFrame([{"pathogen_id": "v1", "host_id": "ncbi:4530", "selection_status": "accepted",
                               "positive_pmids": ["100"]}])
    pathogens = pd.DataFrame([{"pathogen_id": "v1", "pathogen_group": "viruses", "accepted_name": "Virus one"}])
    hosts = pd.DataFrame([{"host_id": "ncbi:4530", "host_group": "crop", "taxon_id": 4530,
                           "accepted_name": "Rice"}])
    prompts = {"pair_feature": {"prompt_version": "prompt-1", "type": "noul", "instructions": "Assess it.",
                                "criteria": {"true": "yes", "false": "no"}}}

    tasks = build_feature_screen_tasks(matches, articles, selection, pathogens, hosts, profiles, prompts,
        run_id="r1", run_config_version="config-1", model_version="model-1", profile_version="profile-1")

    assert len(tasks) == 3
    assert [task["feature_id"] for task in tasks] == ["viral_x", "plant_x", "rice_x"]
    assert [task["scope"] for task in tasks] == ["pathogen", "host", "pathogen_host_pair"]
    assert tasks[0]["host_id"] is None
    assert tasks[1]["host_id"] == "ncbi:4530"
    assert tasks[2]["host_id"] == "ncbi:4530"
    assert all(task["prompt"]["questions"] for task in tasks)


def test_feature_tasks_exclude_unaccepted_pairs_and_abstractless_records():
    profiles = [{"resolved_profile": {"profile_version": "p1",
        "selectors": {"pathogen_group": "fungi", "host_group": "livestock", "host_taxon_id": "9913"},
        "ordered_features": [{"feature_id": "f", "question": "Question?", "question_version": "q1",
                              "scope": "pathogen"}]}}]
    matches = pd.DataFrame([
        {"pathogen_id": "f1", "host_id": "ncbi:9913", "pmid": "100", "query_version": "q"},
        {"pathogen_id": "f1", "host_id": "ncbi:9913", "pmid": "101", "query_version": "q"},
    ])
    articles = pd.DataFrame([{"pmid": "100", "title": "No abstract", "abstract": ""},
                             {"pmid": "101", "title": "Rejected pair", "abstract": "Text."}])
    selection = pd.DataFrame([{"pathogen_id": "f1", "host_id": "ncbi:9913", "selection_status": "rejected",
                               "positive_pmids": []}])
    pathogens = pd.DataFrame([{"pathogen_id": "f1", "pathogen_group": "fungi", "accepted_name": "Fungus"}])
    hosts = pd.DataFrame([{"host_id": "ncbi:9913", "host_group": "livestock", "taxon_id": 9913,
                           "accepted_name": "Cattle"}])
    prompts = {"pair_feature": {"prompt_version": "prompt-1", "type": "noul", "instructions": "Assess it.",
                                "criteria": {"true": "yes", "false": "no"}}}
    tasks = build_feature_screen_tasks(matches, articles, selection, pathogens, hosts, profiles, prompts,
        run_id="r1", run_config_version="config-1", model_version="model-1", profile_version="p1")
    assert tasks == []


def test_catalogue_questions_are_interrogative_and_explorer_displays_the_question():
    root = Path(__file__).resolve().parents[1]
    catalogue = yaml.safe_load((root / "configs/features.yaml").read_text(encoding="utf-8"))
    assert all(str(question.get("question", "")).strip().endswith("?")
               for question in catalogue["questions"].values())
    assert "Does this abstract report that the specified virus is transmitted to rice by the soil-borne fungus Polymyxa graminis?" == catalogue["questions"]["rice_soilborne_vector"]["question"]
    source = (root / "sites/pathogen-atlas-explorer/index.html").read_text(encoding="utf-8")
    assert "textContent=f.question||f.definition" in source
