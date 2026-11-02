import re
import tomllib
from pathlib import Path

from pathogen_atlas_v4 import __version__


def test_package_version_uses_semver_and_matches_project_metadata():
    root = Path(__file__).resolve().parents[1]
    project = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
    project_version = project["project"]["version"]

    assert re.fullmatch(r"\d+\.\d+\.\d+", project_version)
    assert __version__ == project_version

    explorer = (root / "sites/pathogen-atlas-explorer/index.html").read_text(encoding="utf-8")
    assert f'<span class="title-version">v{project_version}</span>' in explorer
