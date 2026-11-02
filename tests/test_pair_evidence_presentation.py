from pathlib import Path


def test_pair_evidence_uses_separate_readable_paper_cards():
    root = Path(__file__).resolve().parents[1]
    html = (root / "sites/pathogen-atlas-explorer/index.html").read_text(encoding="utf-8")

    assert 'class="paper pair-evidence-record"' in html
    assert 'class="paper-head"><a class="paper-pmid"' in html
    assert 'class="paper-confidence">p(yes)' in html
    assert 'class="paper-title">${escapeHtml(record.title' in html
    assert "sorted.map(renderPairEvidencePaper).join('')" in html
    assert "class=\"pmid\"><span><a" not in html
    assert 'class="feature-evidence-toolbar"><label class="small" for="sort-pmids">' in html
    assert 'class="paper-sort" id="sort-pmids"' in html
