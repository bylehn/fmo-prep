"""Tests for FragIt config rendering."""

import pytest
from pathlib import Path


def test_render_config_basic(tmp_path):
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.runner import render_config

    cfg = FragitConfig(
        central_fragment_resname="LIG",
        boundaries=2.0,
        basis="6-31G*",
        calc_mode="2layer",
    )
    out = render_config(cfg, central_fragment_id=42, output_path=tmp_path / "fragit.ini")
    text = out.read_text()

    assert "writer = GAMESS-FMO" in text
    assert "boundaries = 2.0" in text
    assert "basis = 6-31G*" in text
    assert "mp2level = True" in text
    assert "fmolevel = 2" in text
    assert "centralfragment = 42" in text


def test_render_config_mp2_full(tmp_path):
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.runner import render_config

    cfg = FragitConfig(central_fragment_resname="LIG", calc_mode="mp2")
    out = render_config(cfg, central_fragment_id=0, output_path=tmp_path / "fragit.ini")
    text = out.read_text()

    assert "mp2level = True" in text
    assert "fmolevel = 1" in text


def test_render_config_hf(tmp_path):
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.runner import render_config

    cfg = FragitConfig(central_fragment_resname="LIG", calc_mode="hf")
    out = render_config(cfg, central_fragment_id=0, output_path=tmp_path / "fragit.ini")
    text = out.read_text()

    assert "mp2level" not in text


def test_render_config_central_fragment_zero(tmp_path):
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.runner import render_config

    cfg = FragitConfig(central_fragment_resname="LIG")
    out = render_config(cfg, central_fragment_id=0, output_path=tmp_path / "fragit.ini")
    text = out.read_text()
    assert "centralfragment = 0" in text


def test_render_config_returns_path(tmp_path):
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.runner import render_config

    cfg = FragitConfig(central_fragment_resname="LIG")
    out_path = tmp_path / "fragit.ini"
    result = render_config(cfg, central_fragment_id=1, output_path=out_path)
    assert result == out_path
    assert result.exists()


def test_find_central_fragment_id(tmp_path):
    from fmo_prep.fragit.runner import find_central_fragment_id

    inp = tmp_path / "test.inp"
    inp.write_text("""\
 $FMO
      NFRAG=3
      FRGNAM(1)= ALA001,  GLY002,  LIG276
      INDAT(1)=0
            1    -5      0
            6    -10     0
           11    -15     0
 $END
""")
    fid = find_central_fragment_id(inp, "LIG")
    assert fid == 3


def test_find_central_fragment_id_case_insensitive(tmp_path):
    from fmo_prep.fragit.runner import find_central_fragment_id

    inp = tmp_path / "test.inp"
    inp.write_text("""\
 $FMO
      NFRAG=2
      FRGNAM(1)= ALA001,  LZ1276
      INDAT(1)=0
            1    -5      0
            6    -10     0
 $END
""")
    fid = find_central_fragment_id(inp, "lz1")
    assert fid == 2


def test_find_central_fragment_id_not_found(tmp_path):
    from fmo_prep.fragit.runner import find_central_fragment_id

    inp = tmp_path / "test.inp"
    inp.write_text("""\
 $FMO
      NFRAG=2
      FRGNAM(1)= ALA001,  GLY002
      INDAT(1)=0
            1    -5      0
            6    -10     0
 $END
""")
    with pytest.raises(ValueError, match="not found"):
        find_central_fragment_id(inp, "LIG")
