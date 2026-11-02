"""Supplemental NCBI Taxonomy scopes for V4 pathogen groups."""

from __future__ import annotations

import json
from array import array
from pathlib import Path
from typing import Iterable

import pandas as pd


# Parasites are polyphyletic, so this selection explicitly covers major
# protozoan and helminth lineages that infect animals. Oomycetes are a
# separate taxonomy group despite their historical treatment as fungi.
PATHOGEN_TAXONOMY_ROOTS = {
    "parasites": ("5794", "5878", "6157", "6231", "554915", "33682", "2611341", "10232", "35581", "10215", "33209"),
    "oomycetes": ("4762",),
}

PATHOGEN_ADDITION_COLUMNS = [
    "pathogen_id", "pathogen_group", "accepted_name", "taxon_id", "taxonomy_source",
    "taxonomy_release", "rank", "parent_taxon_id", "taxonomy_path",
    "source_record_json", "mapping_status", "source_version",
]
ALIAS_COLUMNS = [
    "pathogen_id", "pathogen_group", "name", "name_type", "source",
    "source_record_id", "source_record_json", "mapping_status",
]


def load_ncbi_pathogen_additions(taxdump_dir: str | Path, groups: Iterable[str], *,
                                 release: str, source_version: str) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str]]:
    """Load species under selected supplemental roots from a pinned NCBI taxdump.

    Returns pathogen rows, alias rows, and all scoped taxon IDs mapped to group (including
    existing V3 taxa) so callers can move legacy misclassified records into
    the correct group while preserving their stable IDs.
    """
    selected = set(map(str, groups)) & set(PATHOGEN_TAXONOMY_ROOTS)
    if not selected:
        return (pd.DataFrame(columns=PATHOGEN_ADDITION_COLUMNS),
                pd.DataFrame(columns=ALIAS_COLUMNS), {})
    root_groups: dict[int, int] = {}
    code_groups: dict[int, str] = {}
    for code, group in enumerate(sorted(selected), start=1):
        code_groups[code] = group
        for tax_id in PATHOGEN_TAXONOMY_ROOTS[group]:
            root_groups[int(tax_id)] = code

    root = Path(taxdump_dir)
    nodes_path, names_path = root / "nodes.dmp", root / "names.dmp"
    if not nodes_path.is_file() or not names_path.is_file():
        raise FileNotFoundError("Pinned NCBI taxdump must contain names.dmp and nodes.dmp")

    parents = array("I", [0])
    species = bytearray(1)
    for line in nodes_path.open(encoding="utf-8"):
        fields = [item.strip() for item in line.split("|")]
        if len(fields) < 3:
            continue
        tax_id, parent_id = int(fields[0]), int(fields[1])
        if tax_id >= len(parents):
            extension = tax_id + 1 - len(parents)
            parents.extend(array("I", [0]) * extension)
            species.extend(b"\0" * extension)
        parents[tax_id] = parent_id
        species[tax_id] = 1 if fields[2] == "species" else 0

    scope = bytearray(len(parents))  # 0 unknown, 255 outside scope, group code inside
    for tax_id in range(1, len(species)):
        if not species[tax_id]:
            continue
        current, lineage = tax_id, []
        while current and current != parents[current]:
            if current in root_groups:
                code = root_groups[current]
                break
            cached = scope[current]
            if cached:
                code = 0 if cached == 255 else cached
                break
            lineage.append(current)
            current = parents[current]
        else:
            code = 0
        mark = code or 255
        for member in lineage:
            scope[member] = mark

    grouped_taxa = {tax_id: code_groups[scope[tax_id]] for tax_id in range(1, len(scope))
                    if species[tax_id] and scope[tax_id] not in (0, 255)}
    scoped_taxa = set(grouped_taxa)
    if not scoped_taxa:
        raise ValueError(f"Selected NCBI Taxonomy roots have no species records: {sorted(selected)}")
    required_names = set(scoped_taxa)
    for tax_id in scoped_taxa:
        current = parents[tax_id]
        while current and current != parents[current]:
            required_names.add(current)
            current = parents[current]

    names: dict[int, list[tuple[str, str]]] = {}
    ranks: dict[int, str] = {}
    for line in nodes_path.open(encoding="utf-8"):
        fields = [item.strip() for item in line.split("|")]
        if len(fields) >= 3 and int(fields[0]) in required_names:
            ranks[int(fields[0])] = fields[2]
    for line in names_path.open(encoding="utf-8"):
        fields = [item.strip() for item in line.split("|")]
        if len(fields) >= 4:
            tax_id = int(fields[0])
            if tax_id in required_names:
                names.setdefault(tax_id, []).append((fields[1], fields[3]))
    scientific: dict[int, str] = {}
    for tax_id, rows in names.items():
        accepted = sorted({name for name, label in rows if label == "scientific name"}, key=str.casefold)
        if accepted:
            scientific[tax_id] = accepted[0]

    pathogens, aliases = [], []
    root_ids_by_group = {group: list(PATHOGEN_TAXONOMY_ROOTS[group]) for group in selected}
    for tax_id, group in sorted(grouped_taxa.items()):
        accepted_name = scientific.get(tax_id)
        if not accepted_name:
            continue
        pathogen_id = f"ncbi:{tax_id}"
        lineage, current = [], tax_id
        while current and current != parents[current]:
            name = scientific.get(current)
            if name:
                lineage.append(f"{ranks.get(current, 'no rank')}={name}")
            current = parents[current]
        record_json = json.dumps({"group_roots": root_ids_by_group[group]}, sort_keys=True)
        pathogens.append({
            "pathogen_id": pathogen_id, "pathogen_group": group, "accepted_name": accepted_name,
            "taxon_id": str(tax_id), "taxonomy_source": "ncbi_taxonomy", "taxonomy_release": release,
            "rank": "species", "parent_taxon_id": str(parents[tax_id]),
            "taxonomy_path": ";".join(reversed(lineage)), "source_record_json": record_json,
            "mapping_status": "accepted", "source_version": source_version,
        })
        names_for_taxon = sorted(set(names.get(tax_id, [])), key=lambda row: (row[1] != "scientific name", row[0].casefold()))
        for name, name_type in names_for_taxon:
            aliases.append({"pathogen_id": pathogen_id, "pathogen_group": group, "name": name,
                            "name_type": "accepted" if name_type == "scientific name" else name_type,
                            "source": "ncbi_taxonomy", "source_record_id": str(tax_id),
                            "source_record_json": record_json,
                            "mapping_status": "accepted" if name_type == "scientific name" else "synonym"})
    return (pd.DataFrame(pathogens, columns=PATHOGEN_ADDITION_COLUMNS),
            pd.DataFrame(aliases, columns=ALIAS_COLUMNS),
            {str(tax_id): group for tax_id, group in grouped_taxa.items()})


def merge_pathogen_catalogues(base_pathogens: pd.DataFrame, base_aliases: pd.DataFrame,
                              additions: pd.DataFrame, addition_aliases: pd.DataFrame,
                              scoped_taxa: dict[str, str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Append supplemental taxa and reclassify overlapping legacy taxa per run."""
    base = base_pathogens.copy()
    aliases = base_aliases.copy()
    base["taxon_id"] = base["taxon_id"].astype("string")
    base["pathogen_group"] = base.apply(
        lambda row: scoped_taxa.get(str(row["taxon_id"]), row["pathogen_group"]), axis=1)
    existing_taxa = set(base["taxon_id"].dropna().astype(str))
    additions = additions.loc[~additions["taxon_id"].astype(str).isin(existing_taxa)].copy()
    existing_ids = set(base["pathogen_id"].astype(str))
    additions = additions.loc[~additions["pathogen_id"].astype(str).isin(existing_ids)]
    pathogens = pd.concat([base, additions], ignore_index=True, sort=False)
    group_by_id = pathogens.set_index(pathogens.pathogen_id.astype(str)).pathogen_group.to_dict()
    if not aliases.empty:
        aliases["pathogen_group"] = aliases["pathogen_id"].astype(str).map(group_by_id).fillna(aliases["pathogen_group"])
    addition_aliases = addition_aliases.loc[addition_aliases["pathogen_id"].astype(str).isin(set(additions.pathogen_id.astype(str)))]
    aliases = pd.concat([aliases, addition_aliases], ignore_index=True, sort=False)
    pathogens = pathogens.drop_duplicates("pathogen_id", keep="first").sort_values(
        ["pathogen_group", "accepted_name"], key=lambda col: col.astype(str).str.casefold()).reset_index(drop=True)
    aliases = aliases.drop_duplicates(["pathogen_id", "name", "name_type"], keep="first").reset_index(drop=True)
    return pathogens, aliases
