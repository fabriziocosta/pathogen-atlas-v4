"""Pair-gate normalization and evidence-span validation."""

from __future__ import annotations

from typing import Any

from pathogen_atlas_v4.schemas import GATE_OUTCOMES, PAIR_RELATIONSHIP

PATHOGEN_SCREEN_YES_THRESHOLD = 0.7


def normalize_pair_gate(*, pathogen_id: str, host_id: str, pmid: str, outcome: str,
                        confidence: float | None, evidence_text: str | None,
                        title_hash: str, abstract_hash: str, run_id: str,
                        relationship_context: str = PAIR_RELATIONSHIP,
                        provenance_status: str = "v4_native",
                        p_yes: float | None = None) -> dict[str, Any]:
    """Create a disease-pair evidence row from text-generating or System 1 decisions."""
    outcome = str(outcome).strip().lower()
    if outcome not in GATE_OUTCOMES:
        raise ValueError(f"Invalid gate outcome {outcome!r}")
    if relationship_context != PAIR_RELATIONSHIP:
        raise ValueError("V4 pair evidence currently supports disease relationships only")
    if not str(pmid).strip() or not pathogen_id or not host_id:
        raise ValueError("Pair evidence requires pathogen_id, host_id, and PMID")
    span = str(evidence_text or "").strip()
    if outcome == "supported" and not span and provenance_status != "system1_jev":
        raise ValueError("Supported pair evidence requires a nonempty evidence span")
    if confidence is not None and not 0 <= float(confidence) <= 1:
        raise ValueError("confidence must be between zero and one")
    if p_yes is not None and not 0 <= float(p_yes) <= 1:
        raise ValueError("p_yes must be between zero and one")
    if provenance_status == "system1_jev" and p_yes is None:
        raise ValueError("System 1 pair evidence must retain Jev's yes probability")
    if p_yes is not None and (outcome == "supported") != (float(p_yes) >= PATHOGEN_SCREEN_YES_THRESHOLD):
        raise ValueError("Pair outcome does not match the 0.7 Jev yes-probability threshold")
    return {"pathogen_id": str(pathogen_id), "host_id": str(host_id),
            "relationship_context": relationship_context, "pmid": str(pmid),
            "title_hash": title_hash, "abstract_hash": abstract_hash,
            "evidence_span": span if outcome == "supported" and span else None,
            "evidence_span_type": "exact_quote" if span and outcome == "supported" else (
                "no_text_generated" if outcome == "supported" and provenance_status == "system1_jev" else None
            ),
            "decision": outcome, "confidence": confidence,
            "p_yes": p_yes, "provenance_status": provenance_status, "run_id": run_id}


def normalize_pathogen_screen(*, pathogen_id: str, host_id: str, pmid: str, answer: str,
                              confidence: float | None, evidence_text: str | None,
                              title: str, abstract: str, run_id: str,
                              p_yes: float | None = None) -> dict[str, Any]:
    """Validate a binary host-specific screen, optionally backed by a Jev probability."""
    answer = str(answer).strip().casefold()
    if answer not in {"yes", "no"}:
        raise ValueError(f"Invalid pathogen-screen answer {answer!r}; expected yes or no")
    if not pathogen_id or not host_id or not str(pmid).strip():
        raise ValueError("Pathogen screening requires pathogen_id, host_id, and PMID")
    if confidence is not None and not 0 <= float(confidence) <= 1:
        raise ValueError("confidence must be between zero and one")
    span = str(evidence_text or "").strip()
    if answer == "yes" and not span and p_yes is None:
        raise ValueError("A text-generating yes pathogen-screen answer requires an evidence span")
    if answer == "yes" and not str(abstract or "").strip():
        raise ValueError("A yes pathogen-screen answer requires a supplied abstract")
    if span and span not in str(abstract or ""):
        raise ValueError("Pathogen-screen evidence span must quote the supplied abstract exactly")
    if p_yes is not None and not 0 <= float(p_yes) <= 1:
        raise ValueError("p_yes must be between zero and one")
    if p_yes is not None and answer != ("yes" if float(p_yes) >= PATHOGEN_SCREEN_YES_THRESHOLD else "no"):
        raise ValueError("Pathogen-screen answer does not match the 0.7 Jev yes-probability threshold")
    provenance = "system1_jev" if p_yes is not None and not span else "v4_native"
    return {"pathogen_id": str(pathogen_id), "host_id": str(host_id), "pmid": str(pmid),
            "relationship_context": "disease", "answer": answer,
            "decision": "supported" if answer == "yes" else "unsupported",
            "confidence": confidence, "p_yes": p_yes,
            "evidence_span": span if answer == "yes" and span else None,
            "evidence_span_type": "exact_quote" if answer == "yes" and span else (
                "no_text_generated" if answer == "yes" and provenance == "system1_jev" else None
            ), "run_id": str(run_id), "provenance_status": provenance}
