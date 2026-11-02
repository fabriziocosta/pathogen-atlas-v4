"""Versioned V4 table contracts and evidence scope definitions."""

from __future__ import annotations

from typing import Final

PATHOGEN_GROUPS: Final[tuple[str, ...]] = ("bacteria", "fungi", "oomycetes", "parasites", "viruses")
FEATURE_SCOPES: Final[tuple[str, ...]] = ("pathogen", "host", "pathogen_host_pair")
GATE_OUTCOMES: Final[tuple[str, ...]] = ("supported", "unsupported", "uncertain", "insufficient_evidence")
PAIR_RELATIONSHIP: Final[str] = "disease"

PATHOGEN_COLUMNS: Final[tuple[str, ...]] = (
    "pathogen_id", "accepted_name", "pathogen_group", "taxon_id", "parent_taxon_id",
    "rank", "taxonomy_source", "taxonomy_release", "source_version", "mapping_status",
)
HOST_COLUMNS: Final[tuple[str, ...]] = (
    "host_id", "accepted_name", "aliases", "host_group", "taxon_id", "parent_taxon_id",
    "rank", "taxonomy_source", "taxonomy_release", "mapping_status",
)
PAIR_COLUMNS: Final[tuple[str, ...]] = (
    "pathogen_id", "host_id", "relationship_context", "pair_id", "run_id",
)
PAIR_EVIDENCE_COLUMNS: Final[tuple[str, ...]] = (
    "pathogen_id", "host_id", "relationship_context", "pmid", "title_hash", "abstract_hash",
    "evidence_span", "evidence_span_type", "decision", "confidence", "p_yes", "provenance_status", "source_project",
    "source_question_version", "source_threshold", "source_record_id",
)
DECISION_VERSION_COLUMNS: Final[tuple[str, ...]] = (
    "query_version", "prompt_version", "feature_version", "profile_version", "model_version",
    "input_hash", "run_id", "provenance_status",
)


def validate_group(value: str) -> str:
    group = str(value).strip().lower()
    if group not in PATHOGEN_GROUPS:
        raise ValueError(f"Unsupported pathogen group {value!r}; expected {PATHOGEN_GROUPS}.")
    return group


def validate_scope(value: str) -> str:
    scope = str(value).strip().lower()
    if scope not in FEATURE_SCOPES:
        raise ValueError(f"Unsupported feature scope {value!r}; expected {FEATURE_SCOPES}.")
    return scope
