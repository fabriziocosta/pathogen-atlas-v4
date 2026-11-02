import json

import pandas as pd

from pathogen_atlas_v4.pubmed_client import pubmed_aliases_by_pathogen


def test_pubmed_aliases_include_verified_vmr_names_but_not_isolates_or_homonyms():
    pathogens = pd.DataFrame([
        {"pathogen_id": "virus:1", "accepted_name": "Tenuivirus oryzabrevis", "taxon_id": "ICTV1"},
        {"pathogen_id": "virus:2", "accepted_name": "Tenuivirus oryzae", "taxon_id": "ICTV2"},
    ])
    vmr_one = {
        "Species": "Tenuivirus oryzabrevis", "ICTV_ID": "ICTV1",
        "Virus name(s)": "rice grassy stunt virus",
        "Virus name abbreviation(s)": "RGSV; RVS",
        "Virus isolate designation": "IRRI isolate",
    }
    vmr_two = {
        "Species": "Tenuivirus oryzae", "ICTV_ID": "ICTV2",
        "Virus name(s)": "rice stripe virus",
        "Virus name abbreviation(s)": "RSV; RVS",
        "Virus isolate designation": "China: isolate 2",
    }
    aliases = pd.DataFrame([
        {"pathogen_id": "virus:1", "name": "rice grassy stunt virus", "name_type": "virus_or_strain_alias", "mapping_status": "source_alias", "source": "ictv_vmr", "source_record_json": json.dumps(vmr_one)},
        {"pathogen_id": "virus:1", "name": "IRRI isolate", "name_type": "virus_or_strain_alias", "mapping_status": "source_alias", "source": "ictv_vmr", "source_record_json": json.dumps(vmr_one)},
        {"pathogen_id": "virus:2", "name": "rice stripe virus", "name_type": "virus_or_strain_alias", "mapping_status": "source_alias", "source": "ictv_vmr", "source_record_json": json.dumps(vmr_two)},
        {"pathogen_id": "virus:1", "name": "ambiguous old name", "name_type": "synonym", "mapping_status": "synonym", "source": "curated", "source_record_json": None},
        {"pathogen_id": "virus:2", "name": "ambiguous old name", "name_type": "synonym", "mapping_status": "synonym", "source": "curated", "source_record_json": None},
        {"pathogen_id": "virus:1", "name": "misassigned alias", "name_type": "virus_or_strain_alias", "mapping_status": "source_alias", "source": "ictv_vmr", "source_record_json": json.dumps({**vmr_one, "Species": "Some other species", "ICTV_ID": "ICTV-other", "Virus name(s)": "misassigned alias"})},
    ])

    searchable, counts = pubmed_aliases_by_pathogen(pathogens, aliases)

    assert set(searchable["virus:1"]) == {"rice grassy stunt virus", "RGSV"}
    assert set(searchable["virus:2"]) == {"rice stripe virus", "RSV"}
    assert "IRRI isolate" not in searchable["virus:1"]
    assert "ambiguous old name" not in searchable.get("virus:1", [])
    assert "ambiguous old name" not in searchable.get("virus:2", [])
    assert "misassigned alias" not in searchable.get("virus:1", [])
    assert counts["aliases_vmr_records_checked"] == 3
    assert counts["aliases_shared_across_pathogens"] >= 2
