"""Pair-gate normalization and evidence-span validation."""

from __future__ import annotations

from typing import Any

from pathogen_atlas_v4.schemas import GATE_OUTCOMES, PAIR_RELATIONSHIP


def normalize_pair_gate(*, pathogen_id: str, host_id: str, pmid: str, outcome: str,
                        confidence: float | None, evidence_text: str | None,
                        title_hash: str, abstract_hash: str, run_id: str,
                        relationship_context: str = PAIR_RELATIONSHIP,
                        provenance_status: str = "v4_native") -> dict[str, Any]:
    """Create a disease-pair evidence row, requiring a quoted span for positives."""
    outcome = str(outcome).strip().lower()
    if outcome not in GATE_OUTCOMES:
        raise ValueError(f"Invalid gate outcome {outcome!r}")
    if relationship_context != PAIR_RELATIONSHIP:
        raise ValueError("V4 pair evidence currently supports disease relationships only")
    if not str(pmid).strip() or not pathogen_id or not host_id:
        raise ValueError("Pair evidence requires pathogen_id, host_id, and PMID")
    if outcome == "supported" and not str(evidence_text or "").strip():
        raise ValueError("Supported pair evidence requires a nonempty evidence span")
    if confidence is not None and not 0 <= float(confidence) <= 1:
        raise ValueError("confidence must be between zero and one")
    return {"pathogen_id": str(pathogen_id), "host_id": str(host_id),
            "relationship_context": relationship_context, "pmid": str(pmid),
            "title_hash": title_hash, "abstract_hash": abstract_hash,
            "evidence_span": evidence_text if outcome == "supported" else None,
            "decision": outcome, "confidence": confidence,
            "provenance_status": provenance_status, "run_id": run_id}
