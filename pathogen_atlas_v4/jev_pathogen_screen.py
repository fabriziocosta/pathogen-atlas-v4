"""OpenRouter Jev adapter for pair-scoped pathogen screening."""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import time
from threading import local
from typing import Any

import requests

DEFAULT_JEV_MODEL = "typesafe/jev-1.14"
DEFAULT_DECISIONS_URL = "https://openrouter.ai/api/alpha/decisions"


def _openrouter_key_from_keychain() -> str | None:
    if sys.platform != "darwin" or shutil.which("security") is None:
        return None
    service = os.environ.get("OPENROUTER_KEYCHAIN_SERVICE", "OpenRouter").strip()
    account = os.environ.get("OPENROUTER_KEYCHAIN_ACCOUNT", "OPENROUTER_API_KEY").strip()
    server = os.environ.get("OPENROUTER_KEYCHAIN_INTERNET_SERVER", "openrouter.ai").strip()
    for kind, name in (("find-generic-password", service), ("find-generic-password", server),
                       ("find-internet-password", server)):
        base = ["security", kind, "-s", name]
        commands = [base + ["-a", account, "-w"], base + ["-w"]] if account else [base + ["-w"]]
        for command in commands:
            try:
                result = subprocess.run(command, capture_output=True, text=True, check=False, timeout=15)
            except (OSError, subprocess.TimeoutExpired):
                continue
            if result.returncode == 0 and result.stdout.strip():
                return result.stdout.strip()
    return None


def openrouter_key_available() -> bool:
    """Check for an OpenRouter key without printing or returning its value."""
    return bool(os.environ.get("OPENROUTER_API_KEY", "").strip() or _openrouter_key_from_keychain())


class JevRequestError(RuntimeError):
    def __init__(self, message: str, *, attempts: int, retryable: bool):
        super().__init__(message)
        self.attempts = attempts
        self.retryable = retryable


class JevResponseError(RuntimeError):
    retryable = True


class OpenRouterJevClient:
    def __init__(self, *, api_key: str | None = None, model: str = DEFAULT_JEV_MODEL,
                 endpoint: str = DEFAULT_DECISIONS_URL, timeout: float = 180,
                 retries: int = 5, backoff_seconds: float = 2.0):
        self.api_key = (api_key or os.environ.get("OPENROUTER_API_KEY", "")).strip() or _openrouter_key_from_keychain() or ""
        if not self.api_key:
            raise ValueError("Set OPENROUTER_API_KEY or add it to the OpenRouter Keychain entry before screening.")
        if not model.strip() or retries < 0 or timeout <= 0 or backoff_seconds < 0:
            raise ValueError("Invalid Jev model, retry, timeout, or backoff setting.")
        self.model = str(model)
        self.endpoint = str(endpoint)
        self.timeout = float(timeout)
        self.retries = int(retries)
        self.backoff_seconds = float(backoff_seconds)
        self._thread_state = local()

    def decide_many(self, state: Any, *, questions: dict[str, dict[str, Any]]) -> dict[str, Any]:
        if not questions:
            raise ValueError("At least one typed Jev question is required.")
        payload = {
            "model": self.model,
            "state": state,
            "questions": {
                key: {"type": question["type"], "instructions": question["instructions"],
                      "criteria": question.get("criteria", {})}
                for key, question in questions.items()
            },
        }
        headers = {"Authorization": f"Bearer {self.api_key}", "Content-Type": "application/json"}
        last_error: Exception | None = None
        retryable = False
        for attempt in range(self.retries + 1):
            session = getattr(self._thread_state, "session", None)
            if session is None:
                session = requests.Session()
                self._thread_state.session = session
            try:
                response = session.post(self.endpoint, json=payload, headers=headers, timeout=self.timeout)
                response.raise_for_status()
                value = response.json()
                if not isinstance(value, dict):
                    raise ValueError("Jev response must be a JSON object.")
                return value
            except Exception as exc:
                last_error = exc
                transient = isinstance(exc, (requests.Timeout, requests.ConnectionError,
                                               requests.exceptions.JSONDecodeError))
                if isinstance(exc, requests.HTTPError):
                    status = exc.response.status_code if exc.response is not None else 0
                    transient = status == 429 or 500 <= status < 600
                retryable = transient
                if not transient or attempt == self.retries:
                    break
                time.sleep(self.backoff_seconds * (2 ** attempt))
        raise JevRequestError(f"Jev request failed: {last_error}", attempts=attempt + 1,
                              retryable=retryable) from last_error


def _abstract_sentences(abstract: str) -> list[str]:
    text = str(abstract or "").strip()
    if not text:
        return []
    sentences = [part.strip() for part in re.split(r"(?<=[.!?])\s+", text) if part.strip()]
    return sentences or [text]


def make_jev_pathogen_screen_handler(*, model: str = DEFAULT_JEV_MODEL,
                                     client: OpenRouterJevClient | None = None):
    """Return exact-span yes/no decisions using Jev's typed choice output.

    Each positive choice maps to one exact sentence from the supplied abstract,
    allowing the regular V4 evidence validator to keep enforcing quote provenance.
    """
    api = client or OpenRouterJevClient(model=model)

    def handler(task: dict[str, Any]) -> dict[str, Any]:
        prompt = task.get("prompt") or {}
        abstract = str(task.get("abstract_text") or "")
        if not prompt.get("system") or not prompt.get("user"):
            raise ValueError("Pathogen screen task is missing its system or user prompt")
        if not abstract.strip():
            raise ValueError("Jev pathogen screening requires a nonempty abstract")

        sentences = _abstract_sentences(abstract)
        sentence_keys = {f"yes_{index:03d}": sentence for index, sentence in enumerate(sentences)}
        criteria = {
            "no": "The abstract does not directly establish that the named organism causes disease in the selected host. Choose no for title-only evidence, background mentions, another organism or host, or unclear claims."
        }
        for key, sentence in sentence_keys.items():
            criteria[key] = f"Yes, and this exact abstract sentence directly supports the named pathogen–host disease relationship: {sentence}"

        response = api.decide_many(
            state={"screening_rules": prompt["system"], "record": prompt["user"]},
            questions={"pathogen_screen": {
                "type": "choice",
                "instructions": (
                    "Determine whether this record provides direct evidence that the named organism causes disease "
                    "in the selected host. Choose no unless one listed abstract sentence directly supports that "
                    "specific relationship. If yes, choose the single strongest supporting sentence option. "
                    "A title or an abstract sentence about a different organism or host is insufficient."
                ),
                "criteria": criteria,
            }},
        )
        answers = response.get("answers")
        decision = answers.get("pathogen_screen") if isinstance(answers, dict) else None
        if not isinstance(decision, dict) or decision.get("type") != "choice":
            raise JevResponseError("Jev response is missing the pathogen_screen choice")
        choice = str(decision.get("choice") or "")
        confidence = decision.get("confidence")
        if confidence is None and isinstance(decision.get("probabilities"), dict):
            confidence = decision["probabilities"].get(choice)
        if choice == "no":
            return {"answer": "no", "confidence": confidence, "evidence_span": "",
                    "rationale": "Jev found no directly supporting sentence in the abstract."}
        evidence_span = sentence_keys.get(choice)
        if evidence_span is None:
            raise JevResponseError(f"Jev returned an unknown pathogen evidence choice: {choice!r}")
        if evidence_span not in abstract:
            raise JevResponseError("Selected Jev evidence sentence is not an exact substring of the source abstract")
        return {"answer": "yes", "confidence": confidence, "evidence_span": evidence_span,
                "rationale": "Jev selected this exact abstract sentence as direct pathogen–host evidence."}

    return handler
