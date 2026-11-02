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
EXPECTED_COORDINATE_COUNTS = {"bacteria": 1703, "viruses": 257, "fungi": 440}


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
    expected = EXPECTED_COORDINATE_COUNTS[group]

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": 1440, "height": 1000})
        page_errors: list[str] = []
        page.on("pageerror", lambda error: page_errors.append(str(error)))
        try:
            page.goto(explorer_url, wait_until="domcontentloaded", timeout=30_000)
            page.wait_for_function(
                """() => document.querySelector('#group')?.value === 'bacteria'
                    && !!document.querySelector('#map .plot-container')
                    && document.querySelector('#catalog-visible')?.textContent.replaceAll(',', '') === '1703'
                    && document.querySelector('#catalog-total')?.textContent.replaceAll(',', '') === '1703'""",
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
                        && document.querySelector('#catalog-visible')?.textContent.replaceAll(',', '') === String(expected)
                        && document.querySelector('#catalog-total')?.textContent.replaceAll(',', '') === String(expected)
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
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": 1440, "height": 1000})
        try:
            page.goto(explorer_url, wait_until="domcontentloaded", timeout=30_000)
            page.wait_for_function("!!document.querySelector('#map .plot-container')", timeout=30_000)
            page.locator("#group").select_option("oomycetes")
            page.wait_for_function(
                """() => document.querySelector('#group')?.value === 'oomycetes'
                    && document.querySelector('#map')?.classList.contains('catalog-mode')
                    && document.querySelector('#catalog-visible')?.textContent.replaceAll(',', '') === '3548'
                    && document.querySelector('#catalog-total')?.textContent.replaceAll(',', '') === '3548'""",
                timeout=30_000,
            )
            assert "Showing 100 of 3,548 matching pathogens" in page.locator("#map").inner_text()
        finally:
            browser.close()


def test_human_coordinate_profiles_are_present_in_deployed_shards():
    metadata = json.loads(gzip.decompress((SITE_SOURCE / "explorer-data.json.gz").read_bytes()))
    columns = metadata["species_columns"]
    pair_profiles_index = columns.index("pair_profiles")
    for group, expected in EXPECTED_COORDINATE_COUNTS.items():
        source_shard = SITE_SOURCE / "taxa" / f"{group}.json.gz"
        deployed_shard = SITE_DIST / "taxa" / f"{group}.json.gz"
        assert deployed_shard.read_bytes() == source_shard.read_bytes()
        rows = json.loads(gzip.decompress(deployed_shard.read_bytes()))
        count = sum(
            1
            for row in rows
            if (profile := (row[pair_profiles_index] or {}).get("ncbi:9606"))
            and profile.get("x") is not None
            and profile.get("y") is not None
        )
        assert count == expected


def test_human_estimates_extrapolate_the_sampled_abstract_fraction():
    source = json.loads(gzip.decompress((SITE_SOURCE / "source/human_pair_data.json.gz").read_bytes()))
    for profile in source["profiles"]:
        sampled = int(profile.get("n_abstracts_fetched") or 0)
        accepted = int(profile.get("n_pathogen_pmids") or profile.get("n_jev_accepted_pmids") or 0)
        total = profile.get("n_pubmed_search_records")
        if sampled == 0:
            continue

        fraction = accepted / sampled
        assert math.isclose(float(profile["pathogenicity_fraction"]), fraction, abs_tol=1e-12)
        if total is not None:
            estimate = fraction * float(total)
            assert math.isclose(float(profile["estimated_pathogen_publications"]), estimate, abs_tol=1e-8)

    e_coli = next(profile for profile in source["profiles"] if profile["accepted_name"] == "Escherichia coli")
    assert e_coli["n_pathogen_pmids"] == 58
    assert e_coli["n_abstracts_fetched"] == 200
    assert e_coli["n_pubmed_search_records"] == 20_529
    assert e_coli["pathogenicity_fraction"] == 0.29
    assert e_coli["estimated_pathogen_publications"] == 5_953.41

    explorer_html = (SITE_SOURCE / "index.html").read_text(encoding="utf-8")
    assert "accepted pathogen PMIDs divided by abstracts fetched" in explorer_html
    assert "extrapolate that fraction to pathogen-wide PubMed records" in explorer_html
    assert "a zero means no pair-specific query records were returned" in explorer_html

