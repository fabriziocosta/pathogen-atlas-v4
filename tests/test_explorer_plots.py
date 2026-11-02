from __future__ import annotations

import gzip
import json
import math
import shutil
import threading
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest
from playwright.sync_api import sync_playwright


REPOSITORY = Path(__file__).resolve().parents[1]
SITE_SOURCE = REPOSITORY / "sites/pathogen-atlas-explorer"
SITE_DIST = SITE_SOURCE / "dist"


def _metadata() -> dict:
    return json.loads(gzip.decompress((SITE_SOURCE / "explorer-data.json.gz").read_bytes()))


def _shard_rows(group: str, *, deployed: bool = True) -> list[list]:
    root = SITE_DIST if deployed else SITE_SOURCE
    path = root / "taxa" / f"{group}.json.gz"
    return json.loads(gzip.decompress(path.read_bytes()))


def _human_coordinate_count(group: str, *, deployed: bool = True) -> int:
    columns = _metadata()["species_columns"]
    pair_profiles_index = columns.index("pair_profiles")
    return sum(
        1
        for row in _shard_rows(group, deployed=deployed)
        if (profile := (row[pair_profiles_index] or {}).get("ncbi:9606"))
        and profile.get("x") is not None
        and profile.get("y") is not None
        and math.isfinite(float(profile["x"]))
        and math.isfinite(float(profile["y"]))
    )


def _feature_profile_counts_from_shards() -> dict[str, dict[str, int]]:
    metadata = _metadata()
    columns = metadata["species_columns"]
    pair_profiles_index = columns.index("pair_profiles")
    counts = {
        host["host_id"]: {group: 0 for group in metadata["groups"]}
        for host in metadata["hosts"]
    }
    for group in metadata["groups"]:
        for row in _shard_rows(group):
            for host_id, profile in (row[pair_profiles_index] or {}).items():
                values = profile.get("feature_values") or {}
                if any(value is not None for value in values.values()):
                    counts.setdefault(host_id, {}).setdefault(group, 0)
                    counts[host_id][group] += 1
    return counts


class QuietRequestHandler(SimpleHTTPRequestHandler):
    def log_message(self, format: str, *args: object) -> None:
        pass


@pytest.fixture(scope="module")
def explorer_url(tmp_path_factory: pytest.TempPathFactory):
    """Serve the same mix of source and dist files used by the Pages workflow."""
    site = tmp_path_factory.mktemp("explorer-site")
    for filename in ("index.html", "explorer-data.json.gz", "feature-catalog.json", "evidence.json.gz"):
        shutil.copy2(SITE_SOURCE / filename, site / filename)
    shutil.copytree(SITE_DIST / "taxa", site / "taxa")
    if (SITE_SOURCE / "assets").exists():
        shutil.copytree(SITE_SOURCE / "assets", site / "assets")

    handler = partial(QuietRequestHandler, directory=str(site))
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{server.server_port}/"
    finally:
        server.shutdown()
        thread.join(timeout=5)
        server.server_close()


@pytest.mark.parametrize("group", ["viruses", "fungi"])
def test_human_2d_plot_renders_after_switching_pathogen_group(explorer_url: str, group: str):
    expected = _human_coordinate_count(group)
    expected_bacteria = _human_coordinate_count("bacteria")
    assert expected > 0
    assert expected_bacteria > 0

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": 1440, "height": 1000})
        page_errors: list[str] = []
        page.on("pageerror", lambda error: page_errors.append(str(error)))
        try:
            page.goto(explorer_url, wait_until="domcontentloaded", timeout=30_000)
            page.wait_for_function(
                """({expected}) => document.querySelector('#group')?.value === 'bacteria'
                    && !!document.querySelector('#map .plot-container')
                    && !document.querySelector('#map')?.classList.contains('catalog-mode')
                    && Number(document.querySelector('#catalog-visible')?.textContent.replace(/[^0-9]/g, '') || 0) === expected
                    && Number(document.querySelector('#catalog-total')?.textContent.replace(/[^0-9]/g, '') || 0) === expected""",
                arg={"expected": expected_bacteria},
                timeout=30_000,
            )

            page.locator("#group").select_option(group)
            page.wait_for_function(
                """({group, expected}) => {
                    const map = document.querySelector('#map');
                    const data = map?.data || [];
                    const plotted = data.reduce((count, trace) => {
                        const sample = trace.customdata?.[0];
                        return count + (Array.isArray(sample) && sample.length > 1 ? trace.x.length : 0);
                    }, 0);
                    return document.querySelector('#host')?.value === 'ncbi:9606'
                        && document.querySelector('#group')?.value === group
                        && Number(document.querySelector('#catalog-visible')?.textContent.replace(/[^0-9]/g, '') || 0) === expected
                        && Number(document.querySelector('#catalog-total')?.textContent.replace(/[^0-9]/g, '') || 0) === expected
                        && !map?.classList.contains('catalog-mode')
                        && !map?.innerText.includes('Loading pathogen taxonomy')
                        && !!map?.querySelector('.plot-container')
                        && !!map?.querySelector('.gl-container canvas')
                        && plotted === expected;
                }""",
                arg={"group": group, "expected": expected},
                timeout=30_000,
            )

            assert page_errors == []
        finally:
            browser.close()


def test_catalog_denominator_remains_full_taxonomy_without_coordinates(explorer_url: str):
    metadata = _metadata()
    group = next((group for group in metadata["groups"] if _human_coordinate_count(group) == 0), None)
    if group is None:
        pytest.skip("Every current pathogen group has human coordinates; no catalog-only group is available")
    expected_catalog_size = len(_shard_rows(group))
    assert expected_catalog_size > 0

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": 1440, "height": 1000})
        try:
            page.goto(explorer_url, wait_until="domcontentloaded", timeout=30_000)
            page.wait_for_function("!!document.querySelector('#map .plot-container')", timeout=30_000)
            page.locator("#group").select_option(group)
            page.wait_for_function(
                """({group, expected}) => document.querySelector('#group')?.value === group
                    && document.querySelector('#map')?.classList.contains('catalog-mode')
                    && Number(document.querySelector('#catalog-visible')?.textContent.replace(/[^0-9]/g, '') || 0) === expected
                    && Number(document.querySelector('#catalog-total')?.textContent.replace(/[^0-9]/g, '') || 0) === expected""",
                arg={"group": group, "expected": expected_catalog_size},
                timeout=30_000,
            )
            assert "matching pathogens" in page.locator("#map").inner_text()
        finally:
            browser.close()


def test_host_and_pathogen_options_are_muted_without_feature_profiles(explorer_url: str):
    metadata = _metadata()
    expected = metadata["profile_availability_counts"]
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": 1440, "height": 1000})
        try:
            page.goto(explorer_url, wait_until="domcontentloaded", timeout=30_000)
            page.wait_for_function("!!document.querySelector('#map .plot-container')", timeout=30_000)

            page.locator("#host-trigger").click()
            host_items = page.locator("#host-menu [role=option]").evaluate_all(
                "items => items.map(item => ({id:item.dataset.value, muted:item.classList.contains('feature-unavailable'), disabled:item.getAttribute('aria-disabled')==='true', color:getComputedStyle(item).color}))"
            )
            assert host_items
            for item in host_items:
                has_profiles = any(value > 0 for value in expected.get(item["id"], {}).values())
                assert item["muted"] is (not has_profiles)
                assert item["disabled"] is False
            host_muted_colors = {item["color"] for item in host_items if item["muted"]}
            host_available_colors = {item["color"] for item in host_items if not item["muted"]}
            if host_muted_colors and host_available_colors:
                assert host_muted_colors.isdisjoint(host_available_colors)

            page.locator("#host-trigger").click()
            page.locator("#group-trigger").click()
            group_items = page.locator("#group-menu [role=option]").evaluate_all(
                "items => items.map(item => ({id:item.dataset.value, muted:item.classList.contains('feature-unavailable'), disabled:item.getAttribute('aria-disabled')==='true', color:getComputedStyle(item).color}))"
            )
            assert group_items
            for item in group_items:
                host_group_counts = expected.get("ncbi:9606", {})
                has_profiles = (any(value > 0 for value in host_group_counts.values())
                                if item["id"] == "all" else host_group_counts.get(item["id"], 0) > 0)
                assert item["muted"] is (not has_profiles)
                assert item["disabled"] is False
            group_muted_colors = {item["color"] for item in group_items if item["muted"]}
            group_available_colors = {item["color"] for item in group_items if not item["muted"]}
            if group_muted_colors and group_available_colors:
                assert group_muted_colors.isdisjoint(group_available_colors)
        finally:
            browser.close()


def test_available_human_estimate_replaces_host_pair_search_zero(explorer_url: str):
    source = json.loads(gzip.decompress((SITE_SOURCE / "source/human_pair_data.json.gz").read_bytes()))
    profile = next(
        profile for profile in source["profiles"]
        if profile.get("estimated_pathogen_publications") is not None
        and profile.get("estimated_pathogen_publications_index") is not None
        and int(profile.get("n_abstracts_fetched") or 0) > 0
        and profile.get("n_pubmed_search_records") is not None
    )
    group = profile["organism_group"]
    pathogen_id = profile.get("pathogen_id") or profile["organism_id"]

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": 1440, "height": 1000})
        try:
            page.goto(explorer_url, wait_until="domcontentloaded", timeout=30_000)
            page.wait_for_function("!!document.querySelector('#map .plot-container')", timeout=30_000)
            if group != "bacteria":
                page.locator("#group").select_option(group)
                page.wait_for_function(
                    "group => document.querySelector('#group')?.value === group && !!document.querySelector('#map .plot-container')",
                    arg=group,
                    timeout=30_000,
                )
            metrics = page.evaluate(
                """async pathogenId => {
                    const species = DATA.species.find(row => row.organism_id === pathogenId);
                    await showProfile(profileForHost(species));
                    return {
                        profile: [...document.querySelectorAll('#profile .metric-box')].map(box => ({
                            label: box.querySelector('span')?.textContent,
                            value: box.querySelector('strong')?.textContent
                        })),
                        groupLabel: document.querySelector('#nsearch-label').textContent,
                        groupEstimate: Number(document.querySelector('#nsearch').textContent.replaceAll(',', '')),
                        expectedGroupEstimate: DATA.species.reduce((total, row) => {
                            const value = profileForHost(row).estimated_pathogen_publications;
                            return total + (value == null ? 0 : Number(value));
                        }, 0)
                    };
                }""",
                pathogen_id,
            )
            labels = {item["label"]: item["value"] for item in metrics["profile"]}
            expected_publications = round(float(profile["estimated_pathogen_publications"]))
            expected_sample_size = int(profile["n_abstracts_fetched"])
            expected_accepted_pmids = int(profile.get("n_pathogen_pmids") or profile.get("n_jev_accepted_pmids") or 0)
            expected_pubmed_records = int(profile["n_pubmed_search_records"])
            assert int("".join(filter(str.isdigit, labels["Estimated pathogen-supporting publications"]))) == expected_publications
            assert labels["Positive sampled records"] == f"{expected_accepted_pmids}/{expected_sample_size}"
            assert int("".join(filter(str.isdigit, labels["Organism PubMed records"]))) == expected_pubmed_records
            assert float(labels["Publication volume score · 0–1"]) == round(float(profile["estimated_pathogen_publications_index"]), 3)
            expected_fraction_score = profile.get("pathogenicity_fraction_index", profile.get("pathogenic_interest_index"))
            assert float(labels["Pathogenicity fraction score · 0–1"]) == round(float(expected_fraction_score), 3)
            expected_fraction = float(profile["pathogenicity_fraction"])
            assert labels["Estimated pathogenicity fraction"] == f"{expected_fraction * 100:.1f}%"
            assert set(labels) == {
                "Publication volume score · 0–1",
                "Pathogenicity fraction score · 0–1",
                "Positive sampled records",
                "Estimated pathogenicity fraction",
                "Estimated pathogen-supporting publications",
                "Organism PubMed records",
            }
            assert metrics["groupLabel"] == "estimated pathogen-supporting publications"
            assert math.isclose(metrics["groupEstimate"], metrics["expectedGroupEstimate"], abs_tol=0.001)
        finally:
            browser.close()


def test_human_coordinate_profiles_are_present_in_deployed_shards():
    metadata = _metadata()
    for group in metadata["groups"]:
        source_shard = SITE_SOURCE / "taxa" / f"{group}.json.gz"
        deployed_shard = SITE_DIST / "taxa" / f"{group}.json.gz"
        assert deployed_shard.read_bytes() == source_shard.read_bytes()
        assert _human_coordinate_count(group, deployed=False) == _human_coordinate_count(group)

    for group in ("bacteria", "viruses", "fungi"):
        assert _human_coordinate_count(group) > 0

    assert metadata["feature_profile_counts"] == _feature_profile_counts_from_shards()


def test_human_estimates_extrapolate_the_sampled_abstract_fraction():
    source = json.loads(gzip.decompress((SITE_SOURCE / "source/human_pair_data.json.gz").read_bytes()))
    checked_profiles = 0
    for profile in source["profiles"]:
        sampled = int(profile.get("n_abstracts_fetched") or 0)
        accepted = int(profile.get("n_pathogen_pmids") or profile.get("n_jev_accepted_pmids") or 0)
        total = profile.get("n_pubmed_search_records")
        if sampled == 0:
            continue

        checked_profiles += 1
        fraction = accepted / sampled
        assert math.isclose(float(profile["pathogenicity_fraction"]), fraction, abs_tol=1e-12)
        if total is not None:
            estimate = fraction * float(total)
            assert math.isclose(float(profile["estimated_pathogen_publications"]), estimate, abs_tol=1e-8)

    assert checked_profiles > 0

    explorer_html = (SITE_SOURCE / "index.html").read_text(encoding="utf-8")
    assert "positive pathogen PMIDs divided by abstracts fetched" in explorer_html
    assert "extrapolate that sample fraction to pathogen-wide PubMed records" in explorer_html
    assert "Publication volume score · 0–1" in explorer_html
    assert "Pathogenicity fraction score · 0–1" in explorer_html

