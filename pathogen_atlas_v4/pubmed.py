"""Versioned host-aware PubMed query construction."""

from __future__ import annotations

import hashlib
import re
from collections.abc import Iterable


def _quote(term: str) -> str:
    return '"' + " ".join(str(term).split()).replace('"', "") + '"'


def build_pair_query(pathogen_names: Iterable[str], host_names: Iterable[str], *, strategy_version: str) -> dict[str, str]:
    """Build a reproducible pathogen AND host query; aliases remain explicit inputs."""
    pathogen_terms = sorted({" ".join(str(x).split()) for x in pathogen_names if len(" ".join(str(x).split())) >= 3}, key=str.casefold)
    host_terms = sorted({" ".join(str(x).split()) for x in host_names if len(" ".join(str(x).split())) >= 3}, key=str.casefold)
    if not pathogen_terms or not host_terms:
        raise ValueError("At least one searchable pathogen name and host name are required")
    pathogens = " OR ".join(f"{_quote(name)}[Title/Abstract]" for name in pathogen_terms)
    hosts = " OR ".join(f"{_quote(name)}[Title/Abstract]" for name in host_terms)
    query = f"({pathogens}) AND ({hosts})"
    digest = hashlib.sha256(f"{strategy_version}\n{query}".encode()).hexdigest()
    return {"query": query, "query_hash": digest, "query_version": strategy_version}
