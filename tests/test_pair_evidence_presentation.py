from pathlib import Path


def test_pair_evidence_uses_separate_readable_paper_cards():
    root = Path(__file__).resolve().parents[1]
    html = (root / "sites/pathogen-atlas-explorer/index.html").read_text(encoding="utf-8")

    assert 'class="paper pair-evidence-record"' in html
    assert 'class="paper-head"><a class="paper-pmid" href="https://pubmed.ncbi.nlm.nih.gov/${escapeHtml(record.pmid)}/"' in html
    assert 'class="paper-confidence">Conf ${Number(confidence||0).toFixed(2)}' in html
    assert "confidence=record.confidence??record.p_yes" in html
    assert 'class="paper-title">${escapeHtml(record.title' in html
    assert 'class="paper-meta">${escapeHtml(meta)}' in html
    assert "sorted.map(renderPairEvidencePaper).join('')" in html
    assert "class=\"pmid\"><span><a" not in html
    assert ".paper-pmid{color:#176c68" in html
    assert "text-decoration:underline" in html
    assert ".paper{padding:14px 14px;border:1px solid #e0e9ed;border-radius:12px" in html
    assert 'class="feature-evidence-toolbar"><label class="small" for="sort-pmids">' in html
    assert 'class="paper-sort" id="sort-pmids"' in html


def test_deployed_evidence_page_matches_the_source_renderer():
    root = Path(__file__).resolve().parents[1]
    source = (root / "sites/pathogen-atlas-explorer/index.html").read_bytes()
    deployed = (root / "sites/pathogen-atlas-explorer/dist/index.html").read_bytes()

    assert deployed == source


def test_feature_evidence_uses_the_same_confidence_label_and_card_layout():
    root = Path(__file__).resolve().parents[1]
    html = (root / "sites/pathogen-atlas-explorer/index.html").read_text(encoding="utf-8")

    assert 'class="paper-confidence">Conf ${Number(x.feature_confidences?.[id]||0).toFixed(2)}' in html
    assert '<article class="paper"><div class="paper-head"><a class="paper-pmid"' in html
    assert 'class="paper-title">${x.title||\'Title unavailable\'}</div>' in html
    assert 'class="paper-meta">${meta}' in html
