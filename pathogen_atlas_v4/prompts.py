"""Deterministic, versioned decision prompt previews for V4."""
from __future__ import annotations

from typing import Any


def _article_block(title: str, abstract: str) -> str:
    title = str(title or "").strip() or "[PubMed title will be inserted here]"
    abstract = str(abstract or "").strip() or "[PubMed abstract will be inserted here]"
    return f"Article title:\n{title}\n\nArticle abstract:\n{abstract}"


def build_pathogen_screen_prompt(*, pathogen_name: str, pathogen_id: str, host_name: str, host_id: str,
                                  article_title: str = "", article_abstract: str = "",
                                  prompt_config: dict[str, Any]) -> dict[str, Any]:
    """Build a host-specific yes/no pathogen evidence request."""
    spec = prompt_config["pathogen_screen"]
    user = (
        f"Question: {spec['question']}\n\n"
        f"Named organism: {pathogen_name}\nPathogen ID: {pathogen_id}\n"
        f"Selected host: {host_name}\nHost ID: {host_id}\n\n"
        f"{_article_block(article_title, article_abstract)}\n\n"
        "Return one JSON object with keys answer, confidence, evidence_span, and rationale. "
        "The answer must be exactly yes or no and must apply to this pathogen–host pair. Answer yes only when the abstract directly supports the disease relationship. "
        "The evidence span must quote the supplied abstract exactly; a title alone cannot support yes. "
        "use an empty string when the answer is no."
    )
    return {"question_version": spec["question_version"], "prompt_version": spec["prompt_version"],
            "system": spec["system_instructions"], "user": user, "output_schema": spec["output_schema"]}


def build_pair_feature_prompt(*, pathogen_name: str, pathogen_id: str, host_name: str, host_id: str,
                              feature: dict[str, Any], article_title: str = "", article_abstract: str = "",
                              prompt_config: dict[str, Any]) -> dict[str, Any]:
    """Build a pair-specific feature assessment request from a resolved feature."""
    if feature.get("scope") != "pathogen_host_pair":
        raise ValueError("Pair feature preview requires a pathogen_host_pair feature")
    spec = prompt_config["pair_feature"]
    user = (
        f"Feature ID: {feature['feature_id']}\nFeature question version: {feature.get('question_version', '')}\n"
        f"Pathogen: {pathogen_name} ({pathogen_id})\nHost: {host_name} ({host_id})\n"
        f"Relationship context: disease\n\nQuestion: {feature['question']}\n\n"
        f"{_article_block(article_title, article_abstract)}\n\n"
        "Return one JSON object with keys outcome, confidence, evidence_span, and rationale. "
        "For supported, evidence_span must quote the supplied text exactly and establish the named pathogen–host disease relationship."
    )
    return {"feature_id": feature["feature_id"], "question_version": feature.get("question_version", ""),
            "prompt_version": spec["prompt_version"], "system": spec["system_instructions"],
            "user": user, "output_schema": spec["output_schema"]}
