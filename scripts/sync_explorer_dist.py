"""Rebuild the exact static bundle consumed by the GitHub Pages workflow."""

from __future__ import annotations

import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SITE_ROOT = ROOT / "sites/pathogen-atlas-explorer"
DIST_ROOT = SITE_ROOT / "dist"
ROOT_FILES = (
    "index.html",
    "explorer-data.json.gz",
    "evidence.json.gz",
    "feature-catalog.json",
)
DIRECTORIES = ("assets", "taxa")


def sync_explorer_dist(site_root: Path = SITE_ROOT, dist_root: Path = DIST_ROOT) -> None:
    """Replace dist with a clean byte-for-byte copy of the canonical site files."""
    site_root = Path(site_root)
    dist_root = Path(dist_root)
    missing = [name for name in ROOT_FILES if not (site_root / name).is_file()]
    missing.extend(name for name in DIRECTORIES if not (site_root / name).is_dir())
    if missing:
        raise FileNotFoundError(f"Cannot build Explorer dist; canonical site is missing: {', '.join(missing)}")

    if dist_root.exists():
        shutil.rmtree(dist_root)
    dist_root.mkdir(parents=True)
    for name in ROOT_FILES:
        shutil.copy2(site_root / name, dist_root / name)
    for name in DIRECTORIES:
        shutil.copytree(site_root / name, dist_root / name)


if __name__ == "__main__":
    sync_explorer_dist()
    print(f"Explorer dist synchronized from {SITE_ROOT}.")
