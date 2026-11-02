"""Small helpers shared by the V4 workflow notebooks."""

from __future__ import annotations

import logging
import os
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone


def find_project_root(start: Path | None = None) -> Path:
    root = (start or Path.cwd()).resolve()
    while root != root.parent and not (root / "pathogen_atlas_v4").is_dir():
        root = root.parent
    if not (root / "pathogen_atlas_v4").is_dir():
        raise RuntimeError("Could not locate the V4 repository root")
    return root


def configured_run_id(root: Path, config_path: Path | None = None) -> str:
    """Read the run ID written by phase 0 for phase-specific logs."""
    import yaml

    path = config_path or (root / "configs" / "pilot_v1.yaml")
    if not path.exists():
        raise FileNotFoundError(f"Run configuration is missing; run 00_workflow.ipynb first: {path}")
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    return str(config.get("run", {}).get("run_id", "interactive"))


def ncbi_contact_email(root: Path, *, prompt: bool = True) -> str:
    """Resolve the required NCBI contact email from env, Git, or a prompt."""
    email = os.environ.get("NCBI_EMAIL", "").strip()
    if email and "@" in email:
        return email
    result = subprocess.run(
        ["git", "config", "--get", "user.email"], cwd=root,
        capture_output=True, text=True, check=False,
    )
    git_email = result.stdout.strip()
    if git_email and "@" in git_email:
        print("Using the configured Git email as the NCBI contact address.")
        return git_email
    if prompt:
        email = input("NCBI contact email (saved only in your local environment): ").strip()
        if email and "@" in email:
            return email
    raise ValueError(
        "Set NCBI_EMAIL or configure a valid `git config user.email` before PubMed acquisition."
    )


def configure_notebook(name: str, run_id: str = "interactive") -> tuple[Path, logging.Logger, Path]:
    root = find_project_root()
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    log_dir = root / "data" / "v4" / "logs" / run_id
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{name}.log"
    logger = logging.getLogger(f"pathogen_atlas_v4.notebooks.{run_id}.{name}")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    for handler in list(logger.handlers):
        handler.close()
        logger.removeHandler(handler)
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s", datefmt="%Y-%m-%d %H:%M")
    console = logging.StreamHandler()
    console.setFormatter(formatter)
    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(console)
    logger.addHandler(file_handler)
    logger.info("Notebook session started at %s", datetime.now(timezone.utc).isoformat())
    logger.info("Project root: %s", root)
    logger.info("Run log: %s", log_path)
    return root, logger, log_path


def log_step(logger: logging.Logger, name: str, **values: object) -> None:
    details = " | ".join(f"{key}={value}" for key, value in values.items())
    logger.info("%s%s", name, f" | {details}" if details else "")
