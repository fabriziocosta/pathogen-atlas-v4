from pathlib import Path


def test_feature_dropdown_exposes_only_applicable_features_with_positive_support():
    explorer = Path(__file__).resolve().parents[1] / "sites/pathogen-atlas-explorer/index.html"
    source = explorer.read_text(encoding="utf-8")

    assert "option.disabled=!applicable.has(id)||!positive.has(id)" in source
    assert "if(selectId==='color')options=options.filter(option=>!option.disabled)" in source
    assert "(s.feature_sources||{})[id]==='literature'" in source
    assert "(species.feature_sources||{})[key]" in source
