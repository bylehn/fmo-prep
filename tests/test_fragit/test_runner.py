"""Tests for FragIt config rendering."""

import pytest
from pathlib import Path


@pytest.mark.skip(reason="render_config not yet implemented")
def test_render_config_basic(tmp_path):
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.runner import render_config

    cfg = FragitConfig(
        central_fragment_resname="LIG",
        boundaries=2.0,
        basis="6-31G*",
        fmo_level=2,
        mp2_level=True,
    )
    out = render_config(cfg, central_fragment_id=42, output_path=tmp_path / "fragit.ini")
    text = out.read_text()

    assert "writer = GAMESS-FMO" in text
    assert "boundaries = 2.0" in text
    assert "basis = 6-31G*" in text
    assert "mp2level = True" in text
    assert "fmolevel = 2" in text
    assert "centralfragment = 42" in text
