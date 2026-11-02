"""Create a deterministic 30-pathogen, four-host pilot with explicit pair IDs."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd
import yaml

HOST_IDS = ["ncbi:9913", "ncbi:9940", "ncbi:4577", "ncbi:4530"]
GROUPS = ["bacteria", "fungi", "viruses"]
SEED = "v4-pilot-v1"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", type=Path, default=Path("../version3/data/01_taxonomy/release/species.parquet"))
    parser.add_argument("--output", type=Path, default=Path("configs/pilot_v1.yaml"))
    args = parser.parse_args()
    species = pd.read_parquet(args.species)
    selected = []
    for group in GROUPS:
        candidates = species.loc[species["organism_group"].astype(str).eq(group), "organism_id"].astype(str).unique()
        chosen = sorted(candidates, key=lambda value: hashlib.sha256(f"{SEED}|{value}".encode()).hexdigest())[:10]
        if len(chosen) != 10:
            raise ValueError(f"Expected at least 10 pathogen IDs for {group}, got {len(chosen)}")
        selected.extend(chosen)
    pairs = [{"pathogen_id": pathogen_id, "host_id": host_id, "relationship_context": "disease"}
             for pathogen_id in selected for host_id in HOST_IDS]
    config = {
        "project": {"id": "pathogen-atlas-v4", "schema_version": 1},
        "run": {"run_id": "pilot-v1", "config_version": "pilot-config-v1",
                "selection_seed": SEED, "pathogen_ids": selected, "host_ids": HOST_IDS,
                "pathogen_groups": GROUPS, "pair_generation": "explicit_pairs", "pairs": pairs},
        "taxonomy": {"host_source": "ncbi_taxonomy", "host_release_manifest": "data/source_taxonomy/ncbi/source_manifest.json"},
        "features": {"catalogue": "configs/features.yaml", "profile_version": "host-aware-profile-v1"},
        "acquisition": {"query_strategy_version": "pathogen-host-title-abstract-v1", "relationship_context": "disease"},
        "decisions": {"gate_question_version": "pathogen-host-disease-gate-v1", "outcomes": ["supported", "unsupported", "uncertain", "insufficient_evidence"]},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(yaml.safe_dump(config, sort_keys=False))
    print(f"Wrote {len(selected)} pathogens and {len(pairs)} explicit host–pathogen pairs to {args.output}")


if __name__ == "__main__":
    main()
