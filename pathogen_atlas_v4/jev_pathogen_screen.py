"""OpenRouter Jev adapter for pair-scoped pathogen screening."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
from threading import local
from typing import Any

import requests

from pathogen_atlas_v4.evidence import PATHOGEN_SCREEN_YES_THRESHOLD

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


def make_jev_pathogen_screen_handler(*, model: str = DEFAULT_JEV_MODEL,
                                     client: OpenRouterJevClient | None = None):
    """Return binary decisions using Jev's typed yes/no (noul) output."""
    api = client or OpenRouterJevClient(model=model)

    def handler(task: dict[str, Any]) -> dict[str, Any]:
        prompt = task.get("prompt") or {}
        abstract = str(task.get("abstract_text") or "")
        if not prompt.get("state") or not prompt.get("questions"):
            raise ValueError("Pathogen screen task is missing its Jev state or typed question")
        if not abstract.strip():
            raise ValueError("Jev pathogen screening requires a nonempty abstract")
        response = api.decide_many(
            state=prompt["state"], questions=prompt["questions"],
        )
        answers = response.get("answers")
        decision = answers.get("pathogen_screen") if isinstance(answers, dict) else None
        if not isinstance(decision, dict) or decision.get("type") != "noul":
            raise JevResponseError("Jev response is missing the pathogen_screen yes/no answer")
        try:
            p_yes = float(decision["noul"])
        except (KeyError, TypeError, ValueError) as exc:
            raise JevResponseError("Jev response is missing the pathogen_screen yes probability") from exc
        if not 0 <= p_yes <= 1:
            raise JevResponseError(f"Jev returned an invalid yes probability: {p_yes!r}")
        return {"answer": "yes" if p_yes >= PATHOGEN_SCREEN_YES_THRESHOLD else "no", "p_yes": p_yes}

    return handler
