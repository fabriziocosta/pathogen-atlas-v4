"""OpenAI Responses API adapter for pathogen–host abstract screening."""
from __future__ import annotations

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
        if not prompt.get("system") or not prompt.get("user"):
            raise ValueError("Pathogen screen task is missing its system or user prompt")
        response = api.responses.parse(
            model=model,
            input=[
                {"role": "system", "content": prompt["system"]},
                {"role": "user", "content": prompt["user"]},
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
