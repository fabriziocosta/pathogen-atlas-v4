"""NCBI Taxonomy dump adapter for normalized host records."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd

from pathogen_atlas_v4.schemas import HOST_COLUMNS

PILOT_HOSTS = {
    "9606": "human",
    "9913": "livestock",  # cattle
    "9940": "livestock",  # sheep
    "9031": "livestock",  # chicken
    "9823": "livestock",  # pig
    "4577": "crop",       # maize
    "4530": "crop",       # rice
    "4565": "crop",       # wheat
    "3847": "crop",       # soybean
    "8049": "fish",       # Atlantic cod
    "8030": "fish",       # Atlantic salmon
}


def load_ncbi_hosts(taxdump_dir: str | Path, taxon_ids: Iterable[str], *, release: str) -> pd.DataFrame:
    """Read selected hosts and their scientific/synonym names from names.dmp and nodes.dmp."""
    root = Path(taxdump_dir)
    names_path, nodes_path = root / "names.dmp", root / "nodes.dmp"
    if not names_path.is_file() or not nodes_path.is_file():
        raise FileNotFoundError("NCBI Taxonomy taxdump must contain names.dmp and nodes.dmp")
    wanted = {str(taxon_id) for taxon_id in taxon_ids}
    names: dict[str, dict[str, list[str]]] = {}
    with names_path.open(encoding="utf-8") as stream:
        for line in stream:
            fields = [field.strip() for field in line.split("|")]
            if len(fields) < 4 or fields[0] not in wanted:
                continue
            names.setdefault(fields[0], {}).setdefault(fields[3], []).append(fields[1])
    nodes: dict[str, tuple[str, str]] = {}
    with nodes_path.open(encoding="utf-8") as stream:
        for line in stream:
            fields = [field.strip() for field in line.split("|")]
            if len(fields) >= 3 and fields[0] in wanted:
                nodes[fields[0]] = (fields[1], fields[2])
    rows = []
    for taxon_id in sorted(wanted, key=int):
        record_names = names.get(taxon_id, {})
        accepted = sorted(record_names.get("scientific name", []), key=str.casefold)
        if len(accepted) != 1 or taxon_id not in nodes:
            raise ValueError(f"Taxon {taxon_id} missing a unique scientific name or node record")
        host_group = PILOT_HOSTS.get(taxon_id, "other")
        rows.append({"host_id": f"ncbi:{taxon_id}", "accepted_name": accepted[0],
                     "aliases": sorted({name for label, values in record_names.items() if label != "scientific name" for name in values}, key=str.casefold),
                     "host_group": host_group, "taxon_id": taxon_id, "parent_taxon_id": nodes[taxon_id][0],
                     "rank": nodes[taxon_id][1], "taxonomy_source": "ncbi_taxonomy",
                     "taxonomy_release": release, "mapping_status": "accepted"})
    return pd.DataFrame(rows, columns=HOST_COLUMNS)
