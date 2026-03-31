"""Tests for GAMESS input postprocessing."""

import pytest
from pathlib import Path


@pytest.fixture
def sample_inp(tmp_path):
    """Minimal FragIt-style .inp without header blocks."""
    content = """\
 $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END
 $FMO
      NFRAG=3
      NBODY=2
      ICHARG(1)=  0,  0,  0
      FRGNAM(1)= ALA001,  GLY002,  LIG001
      INDAT(1)=0
            1    -5      0
            6    -10     0
           11    -15     0
 $END
 $DATA
...
 $END
"""
    p = tmp_path / "raw.inp"
    p.write_text(content)
    return p


@pytest.mark.skip(reason="patch_inp not yet implemented")
def test_patch_inp_adds_header(sample_inp, tmp_path):
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.postprocess import patch_inp

    cfg = FragitConfig(central_fragment_resname="LIG", mwords=100, ngroup=1,
                       resdim=2.0, rcorsd=2.0, mp2_level=True)
    out = tmp_path / "patched.inp"
    patch_inp(sample_inp, cfg, output_path=out)

    text = out.read_text()
    assert "$SYSTEM MWORDS=100" in text
    assert "$GDDI NGROUP=1" in text
    assert "$FMOPRP" in text
    assert "RESDIM=2.0" in text
    assert "RCORSD=2.0" in text
    assert "NLAYER=2" in text
    assert "MPLEVL(1)=2,0" in text
