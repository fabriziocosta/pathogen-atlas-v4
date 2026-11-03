"""LLM-assisted, versioned discovery of host mention cues for PubMed searches."""
from __future__ import annotations

import hashlib
import json
import re
from typing import Any, Literal

from openai import OpenAI
from pydantic import BaseModel, ConfigDict, Field


DEFAULT_HOST_EXPANSION_MODEL = "gpt-6.1-terra"
HOST_EXPANSION_PROMPT_VERSION = "host-mention-expansion-v1"


class HostMentionTerm(BaseModel):
    """A candidate phrase that may signal the target host in an abstract."""

    model_config = ConfigDict(extra="forbid")

    term: str = Field(min_length=3, max_length=120)
    category: Literal[
        "common_name",
        "person_or_patient_descriptor",
        "age_or_sex_descriptor",
        "host_context_phrase",
    ]
    precision_risk: Literal["low", "medium", "high"]
    rationale: str = Field(min_length=1, max_length=400)


class HostMentionExpansion(BaseModel):
    """Structured candidate phrases returned by the host expansion prompt."""

    model_config = ConfigDict(extra="forbid")

    terms: list[HostMentionTerm] = Field(max_length=40)


SYSTEM_PROMPT = """You are preparing a high-recall biomedical literature search.
Suggest short phrases that authors may use in a scientific title, abstract, or
author keyword to indicate that a study concerns the specified host. The phrases
will be used as retrieval cues, not treated as biological synonyms or evidence.

The search will put each phrase in a quoted PubMed Title/Abstract field clause.
That host clause is ORed with the host's accepted scientific name and vetted
taxonomy aliases, then ANDed with an OR clause of pathogen accepted names and
vetted synonyms. A separate configured MeSH clause may be ORed into the host
clause. Retrieved records are candidates; a later abstract-level screen decides
whether the named pathogen is reported to cause disease in the named host.

Return useful host names and real-world ways of indicating the host. For a human
host, consider cues such as person/people/patient, child/children, girl/boy,
woman/man, infant, and age/sex descriptions; use such forms only when they
genuinely indicate the specified host. Prefer phrases
likely to identify the host in biomedical writing. Include a broad cue such as
'patient' only when it can add recall in combination with the pathogen clause,
and mark its precision risk. Return singular/plural or spelling variants as
separate phrases where authors commonly use both, and include concise age/sex
and case-description forms when they identify this host. Do not return generic study words, anatomy terms,
unrelated diseases, PubMed syntax, MeSH headings, or claims that a phrase proves
the host relationship. Do not repeat the accepted name or supplied aliases.
Return only the structured fields requested."""


def build_host_mention_expansion_prompt(*, host_id: str, accepted_name: str,
                                        aliases: list[str] | tuple[str, ...],
                                        common_label: str = "") -> dict[str, Any]:
    """Return the generic expansion prompt and a deterministic provenance hash."""
    state = {
        "host_id": str(host_id),
        "host_label": str(common_label or accepted_name),
        "accepted_scientific_name": str(accepted_name),
        "known_taxonomy_aliases": sorted({str(value).strip() for value in aliases
                                           if str(value).strip()}, key=str.casefold),
        "target_database": "PubMed",
        "purpose": "High-recall discovery of abstracts about disease in this host.",
        "search_use": {
            "text_field": "Each returned phrase is quoted and searched in [Title/Abstract].",
            "boolean_logic": "Host phrases are ORed together, then ANDed with pathogen-name alternatives.",
            "validation": "Retrieved records are candidates and undergo abstract-level pathogen–host disease screening.",
        },
    }
    encoded = json.dumps(state, sort_keys=True, ensure_ascii=False)
    return {
        "prompt_version": HOST_EXPANSION_PROMPT_VERSION,
        "system": SYSTEM_PROMPT,
        "state": state,
        "prompt_sha256": hashlib.sha256((SYSTEM_PROMPT + "\n" + encoded).encode("utf-8")).hexdigest(),
    }


def _validated_terms(terms: list[HostMentionTerm], *, known_aliases: set[str]) -> list[dict[str, str]]:
    """Normalize terms and prevent model text from injecting PubMed syntax."""
    seen = {" ".join(alias.split()).casefold() for alias in known_aliases}
    result = []
    for item in terms:
        term = " ".join(item.term.split()).strip(" ,;.")
        if len(term) < 3 or len(term) > 120 or any(char in term for char in '"[]{}'):
            continue
        if re.search(r"[\r\n\t]", term) or term.casefold() in seen:
            continue
        seen.add(term.casefold())
        result.append({
            "term": term,
            "category": item.category,
            "precision_risk": item.precision_risk,
            "rationale": item.rationale.strip(),
        })
    return result


def expand_host_mentions(*, host_id: str, accepted_name: str, aliases: list[str] | tuple[str, ...],
                          common_label: str = "", model: str = DEFAULT_HOST_EXPANSION_MODEL,
                          client: OpenAI | None = None) -> dict[str, Any]:
    """Generate and validate host mention candidates with a configured OpenAI model."""
    if not model.strip():
        raise ValueError("A configured host-expansion model name is required")
    prompt = build_host_mention_expansion_prompt(
        host_id=host_id, accepted_name=accepted_name, aliases=aliases, common_label=common_label)
    api = client or OpenAI(max_retries=5, timeout=120.0)
    response = api.responses.parse(
        model=model,
        input=[
            {"role": "system", "content": prompt["system"]},
            {"role": "user", "content": json.dumps(prompt["state"], ensure_ascii=False)},
        ],
        text_format=HostMentionExpansion,
        reasoning={"effort": "low"},
        store=False,
    )
    parsed = response.output_parsed
    if parsed is None:
        raise RuntimeError("OpenAI returned no structured host mention expansion")
    known = {str(accepted_name), *map(str, aliases)}
    terms = _validated_terms(parsed.terms, known_aliases=known)
    return {
        "host_id": str(host_id),
        "host_label": str(common_label or accepted_name),
        "accepted_name": str(accepted_name),
        "model": model,
        "prompt_version": prompt["prompt_version"],
        "prompt_sha256": prompt["prompt_sha256"],
        "terms": terms,
    }


def host_expansion_hash(expansions: dict[str, Any]) -> str:
    """Hash only the PubMed-effective terms, independent of model prose/order."""
    effective = {}
    for host_id, record in sorted(expansions.items()):
        effective[str(host_id)] = {
            "title_abstract_terms": sorted({str(item["term"]).casefold()
                                             for item in record.get("terms", [])}),
            "mesh_terms": sorted({str(term).casefold()
                                  for term in record.get("mesh_terms", [])}),
        }
    return hashlib.sha256(json.dumps(effective, sort_keys=True, ensure_ascii=False).encode("utf-8")).hexdigest()
