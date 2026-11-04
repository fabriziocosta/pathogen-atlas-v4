"""Local-only V4 pipeline dashboard server."""

from __future__ import annotations

import json
import gzip
import os
import shutil
import subprocess
import sys
import threading
import time
import uuid
from functools import lru_cache
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
STATE_ROOT = ROOT / "data/v4/dashboard/runs"
UI = Path(__file__).with_name("dashboard.html")
DASHBOARD_MODE = os.environ.get("PATHOGEN_ATLAS_DASHBOARD_MODE", "").strip()
SMOKE_DASHBOARD_PORT = os.environ.get("PATHOGEN_ATLAS_SMOKE_PORT", "8765")
DASHBOARD_VERSION = "1.0.1"
PHASES = [
    ("setup", "Taxonomy, setup & host query expansion"), ("taxonomy", "Taxonomy snapshot review"),
    ("prompts", "Prompt preview"), ("search", "PubMed search & sampling"),
    ("screen", "Pathogen–host screening"), ("estimates", "Publication estimates"),
    ("handoff", "Evidence handoff"), ("profiles", "Feature profiles"),
    ("scoring", "Feature scoring"), ("validation", "Release validation"),
    ("wikipedia", "Wikipedia resolution"), ("gbif", "GBIF resolution"),
    ("build", "Explorer build & embeddings"),
]


def _json(path: Path, default):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return default


def _save(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2), encoding="utf-8")
    temp.replace(path)


@lru_cache(maxsize=1)
def _explorer_metadata():
    path = ROOT / "sites/pathogen-atlas-explorer/explorer-data.json.gz"
    try:
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            return json.load(stream)
    except (OSError, ValueError):
        return {}


@lru_cache(maxsize=5)
def _explorer_pathogens(group: str):
    path = ROOT / "sites/pathogen-atlas-explorer/taxa" / f"{group}.json.gz"
    try:
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            rows = json.load(stream)
        columns = _explorer_metadata().get("species_columns", [])
        id_col, name_col, rank_col = (columns.index(key) for key in ("organism_id", "accepted_name", "rank"))
        return tuple((str(row[id_col]), str(row[name_col])) for row in rows
                     if len(row) > max(id_col, name_col, rank_col) and str(row[rank_col]).casefold() == "species")
    except (OSError, ValueError, KeyError, StopIteration):
        return ()


def _taxonomy_options():
    snapshots = sorted((ROOT / "data/v4/staging").glob("*/taxonomy"))
    hosts, groups, counts = [], [], {}
    explorer_hosts = _explorer_metadata().get("hosts", [])
    explorer_labels = {str(row.get("host_id")): str(row.get("label")) for row in explorer_hosts}
    explorer_order = {str(row.get("host_id")): index for index, row in enumerate(explorer_hosts)}
    host_config = yaml.safe_load((ROOT / "configs/host_choices.yaml").read_text()) or {}
    host_mapping = host_config.get("host_name_to_id", host_config) if isinstance(host_config, dict) else {}
    seen_host_ids = set()
    for label, host_id in host_mapping.items():
        host_id = str(host_id)
        if host_id not in seen_host_ids:
            hosts.append({"id": host_id, "label": explorer_labels.get(host_id, str(label)),
                          "config_name": str(label)})
            seen_host_ids.add(host_id)
    hosts.sort(key=lambda host: explorer_order.get(host["id"], len(explorer_order)))

    group_config = yaml.safe_load((ROOT / "configs/pathogen_taxonomies.yaml").read_text()) or {}
    group_registry = group_config.get("pathogen_groups", {})
    if snapshots:
        tax = snapshots[-1]
        path_frame = pd.read_parquet(tax / "pathogens.parquet")
        counts = path_frame.groupby("pathogen_group").pathogen_id.nunique().to_dict()
    explorer_counts = _explorer_metadata().get("stats", {}).get("pathogen_species_by_group", {})
    if explorer_counts:
        counts = explorer_counts
    groups = [
        {"id": str(group), "label": (f"{config.get('label', str(group).title())} · {int(counts[group]):,} taxa"
                                      if group in counts else config.get("label", str(group).title())),
         "count": int(counts.get(group, 0))}
        for group, config in group_registry.items()
    ]
    return {"hosts": hosts, "groups": groups, "snapshot": snapshots[-1].parent.name if snapshots else None,
            "pathogen_counts": {str(k): int(v) for k, v in counts.items()}}


def _event(run_dir: Path, kind: str, **data):
    row = {"time": time.strftime("%Y-%m-%d %H:%M:%S"), "kind": kind, **data}
    with (run_dir / "events.jsonl").open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(row, ensure_ascii=False) + "\n")
    manifest = _json(run_dir / "manifest.json", {})
    manifest["updated_at"] = row["time"]
    _save(run_dir / "manifest.json", manifest)


def _worker_alive(run_dir: Path):
    pid = _json(run_dir / "worker.json", {}).get("pid")
    if not pid:
        return False


def _other_worker_active(except_run_id: str | None = None):
    for folder in STATE_ROOT.iterdir() if STATE_ROOT.exists() else []:
        if folder.is_dir() and folder.name != except_run_id and _worker_alive(folder):
            return True
    return False
    try:
        os.kill(int(pid), 0)
        return True
    except (OSError, ValueError):
        return False


class Handler(BaseHTTPRequestHandler):
    server_version = "PathogenAtlasDashboard/1"

    def log_message(self, fmt, *args):
        return

    def _body(self):
        length = int(self.headers.get("content-length", "0"))
        return json.loads(self.rfile.read(length) or b"{}")

    def _reply(self, data, status=200, content_type="application/json"):
        body = data if isinstance(data, bytes) else json.dumps(data).encode()
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        path = urlparse(self.path).path
        if path == "/":
            html = UI.read_text(encoding="utf-8")
            mode_labels = {"smoke_test": "Smoke test", "full": "Full cohort"}
            html = html.replace("__DASHBOARD_VERSION__", DASHBOARD_VERSION)
            html = html.replace("__DASHBOARD_MODE__", DASHBOARD_MODE)
            html = html.replace("__DASHBOARD_MODE_LABEL__", mode_labels.get(DASHBOARD_MODE, "Choose run mode"))
            html = html.replace('href="http://127.0.0.1:8765/"',
                                f'href="http://127.0.0.1:{SMOKE_DASHBOARD_PORT}/"')
            return self._reply(html.encode("utf-8"), content_type="text/html; charset=utf-8")
        if path == "/logo.png":
            logo = ROOT / "sites/pathogen-atlas-explorer/assets/pathogen-atlas-v4-logo.png"
            return self._reply(logo.read_bytes(), content_type="image/png")
        if path == "/api/options":
            return self._reply(_taxonomy_options())
        if path == "/api/pathogens":
            query = parse_qs(urlparse(self.path).query)
            group = str((query.get("group") or [""])[0])
            search = str((query.get("q") or [""])[0]).strip()[:100]
            options = _taxonomy_options()
            if group not in options["pathogen_counts"]:
                return self._reply({"pathogens": []})
            pathogens = [(pathogen_id, name) for pathogen_id, name in _explorer_pathogens(group)
                         if not search or search.casefold() in name.casefold()]
            pathogens = sorted(pathogens, key=lambda item: item[1].casefold())[:40]
            return self._reply({"pathogens": [
                {"id": pathogen_id, "label": f"{name} · {pathogen_id}"}
                for pathogen_id, name in pathogens
            ]})
        if path == "/api/runs":
            items = []
            for folder in sorted(STATE_ROOT.glob("*"), reverse=True):
                if folder.is_dir():
                    items.append(_json(folder / "manifest.json", {}))
            return self._reply(items)
        if path.startswith("/api/run/"):
            run_id = path.rsplit("/", 1)[-1]
            run_dir = STATE_ROOT / run_id
            manifest = _json(run_dir / "manifest.json", None)
            if manifest is None:
                return self._reply({"error": "Run not found"}, 404)
            events = []
            event_path = run_dir / "events.jsonl"
            if event_path.exists():
                for line in event_path.read_text(encoding="utf-8").splitlines()[-600:]:
                    try: events.append(json.loads(line))
                    except ValueError: pass
            manifest["events"] = events
            manifest["worker_alive"] = _worker_alive(run_dir)
            return self._reply(manifest)
        if path.startswith("/api/stream/"):
            run_id = path.rsplit("/", 1)[-1]
            run_dir = STATE_ROOT / run_id
            self.send_response(200)
            self.send_header("Content-Type", "text/event-stream")
            self.send_header("Cache-Control", "no-cache")
            self.send_header("Connection", "keep-alive")
            self.end_headers()
            offset = 0
            try:
                while True:
                    file = run_dir / "events.jsonl"
                    if file.exists():
                        with file.open("r", encoding="utf-8") as stream:
                            stream.seek(offset)
                            for line in stream:
                                self.wfile.write(f"data: {line.strip()}\n\n".encode())
                                self.wfile.flush()
                            offset = stream.tell()
                    time.sleep(1)
            except (BrokenPipeError, ConnectionResetError):
                return
        return self._reply({"error": "Not found"}, 404)

    def do_POST(self):
        path = urlparse(self.path).path
        if path == "/api/runs":
            if _other_worker_active():
                return self._reply({"error": "Another V4 pipeline run is active. Pause or finish it before starting a new run."}, 409)
            data = self._body()
            hosts = list(dict.fromkeys(map(str, data.get("hosts", []))))
            groups = list(dict.fromkeys(map(str, data.get("groups", []))))
            mode = str(data.get("mode", "smoke_test"))
            smoke_pathogen_id = str(data.get("smoke_pathogen_id", "")).strip()
            sample = data.get("sample", 10 if mode == "smoke_test" else 200)
            if isinstance(sample, str) and sample.isdecimal():
                sample = int(sample)
            if mode not in {"smoke_test", "full"}:
                return self._reply({"error": "Mode must be 'smoke_test' or 'full'."}, 400)
            if DASHBOARD_MODE and mode != DASHBOARD_MODE:
                return self._reply({"error": f"This address is reserved for {DASHBOARD_MODE.replace('_', ' ')} runs."}, 400)
            if not hosts or not groups or isinstance(sample, bool) or not isinstance(sample, int) or sample < 1:
                return self._reply({"error": "Choose at least one host and pathogen group; sample must be a positive integer."}, 400)
            maximum_sample = 10 if mode == "smoke_test" else 200
            if sample > maximum_sample:
                return self._reply({"error": f"{mode} allows at most {maximum_sample} PMIDs per pair."}, 400)
            if mode == "smoke_test" and (len(hosts) != 1 or len(groups) != 1):
                return self._reply({"error": "Smoke test requires exactly one host and one pathogen group."}, 400)
            options = _taxonomy_options()
            valid_hosts = {h["id"] for h in options["hosts"]}
            valid_groups = {group["id"] for group in options["groups"]}
            if not set(hosts) <= valid_hosts or not set(groups) <= valid_groups:
                return self._reply({"error": "Selection contains an unknown host or pathogen group."}, 400)
            if mode == "smoke_test":
                if not smoke_pathogen_id:
                    return self._reply({"error": "Choose one pathogen species for the smoke test."}, 400)
                if smoke_pathogen_id not in {pathogen_id for pathogen_id, _ in _explorer_pathogens(groups[0])}:
                    return self._reply({"error": "The smoke pathogen must be a species in the selected pathogen group."}, 400)
            run_id = f"dashboard-{time.strftime('%Y%m%d-%H%M%S')}-{uuid.uuid4().hex[:6]}"
            run_dir = STATE_ROOT / run_id
            run_dir.mkdir(parents=True, exist_ok=False)
            manifest = {"run_id": run_id, "status": "queued", "config": {"hosts": hosts, "groups": groups,
                "host_names": [h.get("config_name", h["label"]) for h in options["hosts"] if h["id"] in hosts],
                "mode": mode, "smoke_pathogen_id": smoke_pathogen_id,
                "sample_per_pair": sample, "wikipedia": bool(data.get("wikipedia", True)),
                "gbif": bool(data.get("gbif", True)), "taxonomy_snapshot": options["snapshot"]},
                "phases": {key: {"label": label, "status": "pending", "progress": 0, "completed": 0,
                                  "total": None, "eta_seconds": None, "cache_hits": 0, "cost_usd": None,
                                  "estimate_usd": None} for key, label in PHASES}, "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
                "publish_ready": False}
            _save(run_dir / "manifest.json", manifest)
            _event(run_dir, "run", message="Run configuration saved.")
            self._launch(run_id)
            return self._reply({"run_id": run_id}, 201)
        pieces = path.strip("/").split("/")
        if len(pieces) >= 4 and pieces[:2] == ["api", "run"]:
            run_id, action = pieces[2], pieces[3]
            run_dir = STATE_ROOT / run_id
            manifest = _json(run_dir / "manifest.json", None)
            if not manifest:
                return self._reply({"error": "Run not found"}, 404)
            if action == "pause":
                (run_dir / "pause.request").touch()
                manifest["status"] = "pausing"
                _save(run_dir / "manifest.json", manifest)
                _event(run_dir, "control", message="Pause requested; current provider batch will checkpoint.")
                return self._reply({"ok": True})
            if action == "resume":
                if _worker_alive(run_dir):
                    return self._reply({"error": "Worker is already active."}, 409)
                if _other_worker_active(run_id):
                    return self._reply({"error": "Another V4 pipeline run is active. Pause or finish it before resuming this run."}, 409)
                (run_dir / "pause.request").unlink(missing_ok=True)
                self._launch(run_id)
                return self._reply({"ok": True})
            if action == "publish-preview":
                return self._reply(self._publish_preview(run_dir, manifest))
            if action == "publish":
                data = self._body()
                if data.get("confirm") is not True:
                    return self._reply({"error": "Explicit confirmation required."}, 400)
                return self._reply(self._publish(run_dir, manifest))
        return self._reply({"error": "Not found"}, 404)

    def _launch(self, run_id):
        run_dir = STATE_ROOT / run_id
        command = [sys.executable, "-m", "pathogen_atlas_v4.dashboard_worker", "--run-id", run_id]
        env = dict(os.environ, PYTHONPATH=str(ROOT) + os.pathsep + os.environ.get("PYTHONPATH", ""))
        process = subprocess.Popen(command, cwd=ROOT, env=env, stdin=subprocess.DEVNULL,
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                   start_new_session=True)
        _save(run_dir / "worker.json", {"pid": process.pid})

    def _publish_preview(self, run_dir, manifest):
        if manifest.get("status") != "completed" or not manifest.get("publish_ready"):
            return {"error": "Publishing is enabled after Explorer validation and packaging complete."}
        paths = _release_paths()
        return {"run_id": manifest["run_id"], "version": _next_version(), "files": paths,
                "summary": "Additive V4 Explorer update. Only listed Explorer payloads and release metadata will be staged."}

    def _publish(self, run_dir, manifest):
        preview = self._publish_preview(run_dir, manifest)
        if "error" in preview:
            return preview
        branch = subprocess.run(["git", "branch", "--show-current"], cwd=ROOT, text=True, capture_output=True, check=True).stdout.strip()
        if branch != "main":
            return {"error": f"Publishing requires the Pages deployment branch main; current branch is {branch}."}
        upstream = subprocess.run(["git", "rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{upstream}"],
                                  cwd=ROOT, text=True, capture_output=True, check=False)
        remote_branch = upstream.stdout.strip()
        if upstream.returncode != 0 or "/" not in remote_branch:
            return {"error": "The current branch has no configured upstream remote."}
        remote, remote_main = remote_branch.split("/", 1)
        if remote_main != "main":
            return {"error": f"Configured upstream {remote_branch} is not the Pages branch main."}
        fetched = subprocess.run(["git", "fetch", remote], cwd=ROOT, text=True, capture_output=True)
        if fetched.returncode:
            return {"error": fetched.stderr[-3000:] or fetched.stdout[-3000:]}
        counts = subprocess.run(["git", "rev-list", "--left-right", "--count", f"HEAD...{remote}/main"],
                                cwd=ROOT, text=True, capture_output=True, check=True).stdout.split()
        if len(counts) == 2 and int(counts[1]) > 0:
            return {"error": f"{remote}/main has {counts[1]} commit(s) missing locally. Sync the checkout before publishing."}
        # Only after the release target and remote are verified do we bump release metadata.
        result = subprocess.run([sys.executable, "-m", "pathogen_atlas_v4.dashboard_worker", "--bump-release",
                                 manifest["run_id"]], cwd=ROOT, text=True, capture_output=True)
        if result.returncode:
            return {"error": result.stderr[-3000:] or result.stdout[-3000:]}
        paths = _release_paths() + ["CHANGELOG.md", "pyproject.toml", "pathogen_atlas_v4/__init__.py",
                                    "sites/pathogen-atlas-explorer/dist/index.html"]
        # A path-limited commit excludes any existing staged or working-tree changes outside the allowlist.
        version = _version()
        commit = subprocess.run(["git", "commit", "--only", "-m", f"Release Pathogen Atlas V4 {version}", "--", *paths], cwd=ROOT,
                                text=True, capture_output=True)
        if commit.returncode:
            return {"error": commit.stderr[-3000:] or commit.stdout[-3000:]}
        push = subprocess.run(["git", "push", remote, "main"], cwd=ROOT, text=True, capture_output=True)
        if push.returncode:
            return {"error": push.stderr[-3000:] or push.stdout[-3000:], "commit": commit.stdout.strip()}
        manifest["publish_status"] = "pushed"
        _save(run_dir / "manifest.json", manifest)
        sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=True).stdout.strip()
        return {"status": "pushed", "version": version, "commit": sha,
                "workflow": self._wait_for_pages(sha)}

    def _wait_for_pages(self, sha):
        if not shutil.which("gh"):
            return "Push succeeded; GitHub CLI is unavailable, so Pages workflow status was not polled."
        for _ in range(30):
            result = subprocess.run(["gh", "run", "list", "--workflow", "deploy-pages.yml", "--branch", "main", "--commit", sha,
                                     "--limit", "1", "--json", "status,conclusion,url"],
                                    cwd=ROOT, text=True, capture_output=True)
            if result.returncode == 0:
                try:
                    rows = json.loads(result.stdout)
                    if rows and rows[0].get("status") == "completed":
                        return {"status": rows[0].get("conclusion"), "url": rows[0].get("url")}
                except (ValueError, TypeError):
                    pass
            time.sleep(8)
        return "Push succeeded; Pages workflow is still running or its status is unavailable."


def _release_paths():
    site = ROOT / "sites/pathogen-atlas-explorer"
    candidates = [site / "index.html", site / "explorer-data.json.gz", site / "evidence.json.gz", site / "feature-catalog.json"]
    candidates += sorted((site / "dist/taxa").glob("*.json.gz")) if (site / "dist/taxa").exists() else []
    return [str(path.relative_to(ROOT)) for path in candidates if path.is_file()]


def _version():
    import re
    text = (ROOT / "pyproject.toml").read_text()
    return re.search(r'^version = "([^"]+)"', text, re.M).group(1)


def _next_version():
    import re
    match = re.fullmatch(r"(\d+)\.(\d+)\.(\d+)", _version())
    if not match: return _version()
    return f"{match[1]}.{match[2]}.{int(match[3]) + 1}"


def main():
    os.chdir(ROOT)
    STATE_ROOT.mkdir(parents=True, exist_ok=True)
    server = ThreadingHTTPServer(("127.0.0.1", int(os.environ.get("PATHOGEN_ATLAS_DASHBOARD_PORT", "8765"))), Handler)
    mode_label = {"smoke_test": "Smoke", "full": "Full"}.get(DASHBOARD_MODE, "Pipeline")
    print(f"Pathogen Atlas V4 {mode_label} dashboard: http://{server.server_address[0]}:{server.server_address[1]}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
