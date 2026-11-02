import pandas as pd

from pathogen_atlas_v4.pathogen_screening import build_pathogen_screen_tasks


def test_pathogen_screen_prompt_names_the_taxon_checked_literature_alias():
    matches = pd.DataFrame([{
        "pathogen_id": "viruses:rice-virus", "host_id": "ncbi:4530", "pmid": "12345678",
        "query_hash": "query-hash", "query_version": "ictv-vmr-v3", "sample_selection": "oldest",
    }])
    articles = pd.DataFrame([{"pmid": "12345678", "title": "Rice stripe disease",
                              "abstract": "Rice stripe virus causes disease in rice."}])
    pathogens = pd.DataFrame([{"pathogen_id": "viruses:rice-virus", "accepted_name": "Tenuivirus oryzaclavatae"}])
    hosts = pd.DataFrame([{"host_id": "ncbi:4530", "accepted_name": "Oryza sativa"}])
    prompts = {"pathogen_screen": {"type": "noul", "instructions": "Assess the named pathogen and host.",
                                    "question_version": "screen-v1", "prompt_version": "noul-v1"}}

    tasks = build_pathogen_screen_tasks(
        matches, articles, pathogens, hosts, prompts,
        run_id="test", run_config_version="test-v1", model_version="test-model",
        pathogen_names_by_id={"viruses:rice-virus": ["rice stripe virus", "RSV"]},
    )

    pathogen_label = tasks[0]["prompt"]["state"]["pathogen"]["name"]
    assert pathogen_label.startswith("Tenuivirus oryzaclavatae (also known as ")
    assert "rice stripe virus" in pathogen_label
    assert "RSV" in pathogen_label
