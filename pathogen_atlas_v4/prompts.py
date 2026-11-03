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
    result = build_feature_question(pathogen_name=pathogen_name, pathogen_id=pathogen_id,
        host_name=host_name, host_id=host_id, feature=feature, article_title=article_title,
        article_abstract=article_abstract, prompt_config=prompt_config)
    result["questions"] = {"pair_feature": result["questions"]["feature"]}
    return result


def build_feature_question(*, pathogen_name: str, pathogen_id: str, host_name: str, host_id: str,
                           feature: dict[str, Any], article_title: str = "", article_abstract: str = "",
                           prompt_config: dict[str, Any]) -> dict[str, Any]:
    """Build a scope-aware Jev question for pathogen, host, or pair features.

    The question asks whether this supplied abstract explicitly reports the
    claim. It does not ask the model to infer biological truth from absence.
    """
    scope = feature.get("scope")
    if scope not in {"pathogen", "host", "pathogen_host_pair"}:
        raise ValueError(f"Unsupported feature scope: {scope!r}")
    spec = prompt_config["pair_feature"]
    state: dict[str, Any] = {
        "feature_id": feature["feature_id"],
        "feature_question_version": feature.get("question_version", ""),
        "feature_question": feature["question"],
        "feature_scope": scope,
        "article": {"title": str(article_title or ""), "abstract": str(article_abstract or "")},
    }
    if scope in {"pathogen", "pathogen_host_pair"}:
        state["pathogen"] = {"name": pathogen_name, "id": pathogen_id}
    if scope in {"host", "pathogen_host_pair"}:
        state["host"] = {"name": host_name, "id": host_id}
    if scope == "pathogen_host_pair":
        state["relationship_context"] = "disease"
    question = {"type": spec["type"], "instructions": spec["instructions"],
                "criteria": spec.get("criteria", {})}
    return {"feature_id": feature["feature_id"], "question_version": feature.get("question_version", ""),
            "prompt_version": spec["prompt_version"], "state": state,
            "questions": {"feature": question}}
