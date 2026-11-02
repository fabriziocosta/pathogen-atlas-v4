"""Deterministic, versioned decision prompt previews for V4."""
from __future__ import annotations

from typing import Any


def build_pathogen_screen_question(*, pathogen_name: str, pathogen_id: str, host_name: str, host_id: str,
                                  article_title: str = "", article_abstract: str = "",
                                  prompt_config: dict[str, Any]) -> dict[str, Any]:
    """Build Jev state and a typed yes/no question for pathogen screening."""
    spec = prompt_config["pathogen_screen"]
    state = {
        "pathogen": {"name": pathogen_name, "id": pathogen_id},
        "host": {"name": host_name, "id": host_id},
        "article": {"title": str(article_title or ""), "abstract": str(article_abstract or "")},
    }
    question = {"type": spec["type"], "instructions": spec["instructions"],
                "criteria": spec.get("criteria", {})}
    return {"question_version": spec["question_version"], "prompt_version": spec["prompt_version"],
            "state": state, "questions": {"pathogen_screen": question}}


def build_pair_feature_question(*, pathogen_name: str, pathogen_id: str, host_name: str, host_id: str,
                              feature: dict[str, Any], article_title: str = "", article_abstract: str = "",
                              prompt_config: dict[str, Any]) -> dict[str, Any]:
    """Build Jev state and a typed yes/no question for a pair-specific feature."""
    if feature.get("scope") != "pathogen_host_pair":
        raise ValueError("Pair feature preview requires a pathogen_host_pair feature")
    spec = prompt_config["pair_feature"]
    state = {
        "feature_id": feature["feature_id"],
        "feature_question_version": feature.get("question_version", ""),
        "feature_question": feature["question"],
        "pathogen": {"name": pathogen_name, "id": pathogen_id},
        "host": {"name": host_name, "id": host_id},
        "relationship_context": "disease",
        "article": {"title": str(article_title or ""), "abstract": str(article_abstract or "")},
    }
    question = {"type": spec["type"], "instructions": spec["instructions"],
                "criteria": spec.get("criteria", {})}
    return {"feature_id": feature["feature_id"], "question_version": feature.get("question_version", ""),
            "prompt_version": spec["prompt_version"], "state": state,
            "questions": {"pair_feature": question}}
