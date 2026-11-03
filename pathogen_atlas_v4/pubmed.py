"""Versioned host-aware PubMed query construction."""

from __future__ import annotations

import hashlib
import re
from collections.abc import Iterable


def _quote(term: str) -> str:
    return '"' + " ".join(str(term).split()).replace('"', "") + '"'


def build_pair_query(pathogen_names: Iterable[str], host_names: Iterable[str], *, strategy_version: str,
                     host_mention_terms: Iterable[str] = (), host_mesh_terms: Iterable[str] = ()) -> dict[str, str]:
    """Build a pathogen AND host query with explicit alias and mention branches.

    Taxonomic names and host-mention cues are quoted Title/Abstract phrases and
    are ORed within the host clause. Optional controlled-vocabulary headings
    are ORed into that same clause using PubMed's MeSH field.
    """
    pathogen_terms = sorted({" ".join(str(x).split()) for x in pathogen_names if len(" ".join(str(x).split())) >= 3}, key=str.casefold)
    host_terms = sorted({" ".join(str(x).split()) for x in [*host_names, *host_mention_terms]
                         if len(" ".join(str(x).split())) >= 3}, key=str.casefold)
    mesh_terms = sorted({" ".join(str(x).split()) for x in host_mesh_terms
                         if len(" ".join(str(x).split())) >= 3}, key=str.casefold)
    if any(any(char in term for char in '"[]{}') for term in [*pathogen_terms, *host_terms, *mesh_terms]):
        raise ValueError("Search terms cannot contain quotes or PubMed field-tag brackets")
    if not pathogen_terms or not (host_terms or mesh_terms):
        raise ValueError("At least one searchable pathogen name and host name, mention, or MeSH term are required")
    pathogens = " OR ".join(f"{_quote(name)}[Title/Abstract]" for name in pathogen_terms)
    host_clauses = [f"{_quote(name)}[Title/Abstract]" for name in host_terms]
    host_clauses.extend(f"{_quote(name)}[MeSH Terms]" for name in mesh_terms)
    hosts = " OR ".join(host_clauses)
    query = f"({pathogens}) AND ({hosts})"
    digest = hashlib.sha256(f"{strategy_version}\n{query}".encode()).hexdigest()
    return {"query": query, "query_hash": digest, "query_version": strategy_version}


def build_pathogen_query(pathogen_names: Iterable[str], *, strategy_version: str) -> dict[str, str]:
    """Build a host-independent PubMed query for papers mentioning a pathogen."""
    names = sorted({" ".join(str(name).split()) for name in pathogen_names
                    if len(" ".join(str(name).split())) >= 3}, key=str.casefold)
    if not names:
        raise ValueError("At least one searchable pathogen name is required")
    terms = " OR ".join(f"{_quote(name)}[Title/Abstract]" for name in names)
    query = f"({terms})"
    digest = hashlib.sha256(f"{strategy_version}\n{query}".encode()).hexdigest()
    return {"query": query, "query_hash": digest, "query_version": strategy_version}
