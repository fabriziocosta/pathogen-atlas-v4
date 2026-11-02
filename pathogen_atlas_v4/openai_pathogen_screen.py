"""OpenAI Responses API adapter for pathogen–host abstract screening."""
from __future__ import annotations

import json
from typing import Literal

from openai import OpenAI
from pydantic import BaseModel, ConfigDict, Field


class PathogenScreenResult(BaseModel):
    """Structured decision expected by the pair-screening ledger."""

    model_config = ConfigDict(extra="forbid")

    answer: Literal["yes", "no"]
    confidence: float = Field(ge=0, le=1)
    evidence_span: str
    rationale: str


def make_openai_pathogen_screen_handler(*, model: str, client: OpenAI | None = None):
    """Build a resumable runner handler using the requested OpenAI model."""
    if not model.strip():
        raise ValueError("A configured OpenAI model name is required")
    api = client or OpenAI(max_retries=5, timeout=120.0)

    def handler(task: dict) -> dict:
        prompt = task.get("prompt") or {}
        state = prompt.get("state")
        questions = prompt.get("questions") or {}
        question = questions.get("pathogen_screen")
        if not isinstance(state, dict) or not isinstance(question, dict):
            raise ValueError("Pathogen screen task is missing its state or typed question")
        criteria = question.get("criteria", {})
        system = (
            str(question.get("instructions", "")) + "\n\n"
            f"Yes criterion: {criteria.get('true', '')}\n"
            f"No criterion: {criteria.get('false', '')}\n"
            "Return JSON with answer, confidence, evidence_span, and rationale. "
            "Use yes or no. For yes, evidence_span must exactly quote the abstract in the state."
        )
        response = api.responses.parse(
            model=model,
            input=[
                {"role": "system", "content": system},
                {"role": "user", "content": json.dumps(state, ensure_ascii=False)},
            ],
            text_format=PathogenScreenResult,
            reasoning={"effort": "low"},
            store=False,
        )
        parsed = response.output_parsed
        if parsed is None:
            raise RuntimeError("OpenAI returned no structured pathogen-screen decision")
        return parsed.model_dump()

    return handler
