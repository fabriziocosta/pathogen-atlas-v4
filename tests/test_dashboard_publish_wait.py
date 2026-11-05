import json
import subprocess

from pathogen_atlas_v4 import dashboard


def test_publish_waits_for_an_active_pages_workflow(monkeypatch):
    responses = iter([
        [{"databaseId": 41, "status": "in_progress", "url": "https://github.com/run/41"}],
        [{"databaseId": 41, "status": "completed", "url": "https://github.com/run/41"}],
    ])
    calls = []

    def run(command, **kwargs):
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, json.dumps(next(responses)), "")

    waited = []
    monkeypatch.setattr(dashboard.shutil, "which", lambda name: "/usr/bin/gh")
    monkeypatch.setattr(dashboard.subprocess, "run", run)
    monkeypatch.setattr(dashboard, "_publish_state", lambda *args, **kwargs: waited.append((args, kwargs)))
    monkeypatch.setattr(dashboard.time, "sleep", lambda seconds: None)

    handler = object.__new__(dashboard.Handler)
    assert handler._wait_for_previous_pages("run-1", timeout_seconds=10) == {}
    assert len(calls) == 2
    assert waited[0][0][1] == "waiting"


def test_publish_fails_closed_when_github_cli_is_unavailable(monkeypatch):
    monkeypatch.setattr(dashboard.shutil, "which", lambda name: None)
    handler = object.__new__(dashboard.Handler)
    result = handler._wait_for_previous_pages("run-1")
    assert "nothing was pushed" in result["error"]
