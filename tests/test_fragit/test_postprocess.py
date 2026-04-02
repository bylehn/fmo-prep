"""Tests for GAMESS input postprocessing."""

import pytest
from pathlib import Path

FIXTURES = Path(__file__).parent.parent / "fixtures"


@pytest.fixture
def raw_inp(tmp_path):
    """Minimal FragIt-style .inp that already contains header blocks (as FragIt writes them)."""
    content = """\
 $SYSTEM MWORDS=50 $END
 $GDDI NGROUP=1 $END
 $SCF CONV=1E-8 $END
 $CONTRL NPRINT=-5 ISPHER=1
         RUNTYP=ENERGY
 $END
 $BASIS GBASIS=STO NGAUSS=3 $END
 $FMOPRP NPRINT=9 $END
 $FMO
      NFRAG=3
      NBODY=2
      RESDIM=2.0
      RCORSD=2.0
      NLAYER=2
      MPLEVL(1)=0,2
      ICHARG(1)=  0,  0,  1
      FRGNAM(1)= ALA001,  GLY002,  LIG001
      INDAT(1)=0
            1    -5      0
            6    -10     0
           11    -15     0
 $END
 $DATA
 title
 C1
...
 $END
"""
    p = tmp_path / "raw.inp"
    p.write_text(content)
    return p


# --- 2layer mode ---

def test_patch_inp_2layer_header(raw_inp, tmp_path):
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.postprocess import patch_inp

    cfg = FragitConfig(central_fragment_resname="LIG", mwords=100, ngroup=2,
                       basis="6-31G*", calc_mode="2layer")
    out = tmp_path / "patched.inp"
    patch_inp(raw_inp, cfg, output_path=out)
    text = out.read_text()

    assert " $SYSTEM MWORDS=100 $END" in text
    assert " $GDDI NGROUP=2 $END" in text
    assert " $SCF CONV=1E-6 DIRSCF=.T. NPUNCH=0 DIIS=.T. SOSCF=.F. $END" in text
    assert " $CONTRL NPRINT=-5 ISPHER=1 MAXIT=100" in text
    assert "SCFTYP=RHF" in text
    assert " $FMOPRP NPRINT=9 NGUESS=2 MAXIT=100 $END" in text
    assert "IPIEDA" not in text
    assert "$PCM" not in text
    assert "MPLEVL(1)=0,2" in text


# --- mp2 mode ---

def test_patch_inp_mp2_header(raw_inp, tmp_path):
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.postprocess import patch_inp

    cfg = FragitConfig(central_fragment_resname="LIG", mwords=125, ngroup=1,
                       basis="6-31G*", calc_mode="mp2")
    out = tmp_path / "patched.inp"
    patch_inp(raw_inp, cfg, output_path=out)
    text = out.read_text()

    assert " $SCF CONV=1E-7 DIRSCF=.T. NPUNCH=0 DIIS=.F. SOSCF=.T. $END" in text
    assert "SCFTYP" not in text
    assert "MAXIT" not in text.split("$FMOPRP")[0]  # MAXIT only in FMOPRP, not CONTRL
    assert "PRTDST(1)=100.0,0.5,0.6,0.0" in text
    assert "IPIEDA=2" in text
    assert "$PCM" not in text
    assert "MPLEVL(1)=2" in text
    assert "MPLEVL(1)=0,2" not in text


# --- hf mode ---

def test_patch_inp_hf_header(raw_inp, tmp_path):
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.postprocess import patch_inp

    cfg = FragitConfig(central_fragment_resname="LIG", mwords=125, ngroup=1,
                       basis="6-31G*", calc_mode="hf")
    out = tmp_path / "patched.inp"
    patch_inp(raw_inp, cfg, output_path=out)
    text = out.read_text()

    assert " $SCF CONV=1E-6 DIRSCF=.T. NPUNCH=0 DIIS=.F. SOSCF=.T. $END" in text
    assert "SCFTYP" not in text
    assert "IPIEDA" not in text
    assert "PRTDST" not in text
    assert "$PCM" not in text
    assert "MPLEVL(1)=0" in text
    assert "MPLEVL(1)=0,2" not in text


# --- implicit solvent ---

def test_patch_inp_implicit_solvent(raw_inp, tmp_path):
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.postprocess import patch_inp

    cfg = FragitConfig(central_fragment_resname="LIG", mwords=125, ngroup=1,
                       basis="6-31G*", calc_mode="mp2", implicit_solvent=True)
    out = tmp_path / "patched.inp"
    patch_inp(raw_inp, cfg, output_path=out)
    text = out.read_text()

    assert " $PCM SOLVNT=WATER IFMO=1 ICOMP=0 $END" in text
    assert "IPIEDA=1" in text
    assert " $SCF CONV=1E-6 DIRSCF=.T. NPUNCH=0 DIIS=.F. SOSCF=.T. $END" in text
    assert "IPIEDA=2" not in text
    assert "PRTDST" not in text


def test_patch_inp_pcm_appears_after_basis(raw_inp, tmp_path):
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.postprocess import patch_inp

    cfg = FragitConfig(central_fragment_resname="LIG", mwords=125, ngroup=1,
                       calc_mode="mp2", implicit_solvent=True)
    out = tmp_path / "patched.inp"
    patch_inp(raw_inp, cfg, output_path=out)
    text = out.read_text()

    basis_pos = text.index("$BASIS")
    pcm_pos   = text.index("$PCM")
    fmoprp_pos = text.index("$FMOPRP")
    assert basis_pos < pcm_pos < fmoprp_pos


# --- shared behaviour ---

def test_patch_inp_no_duplicate_blocks(raw_inp, tmp_path):
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.postprocess import patch_inp

    cfg = FragitConfig(central_fragment_resname="LIG", mwords=100, ngroup=1, calc_mode="2layer")
    out = tmp_path / "patched.inp"
    patch_inp(raw_inp, cfg, output_path=out)
    text = out.read_text()

    for block in ["$SYSTEM", "$GDDI", "$SCF", "$CONTRL", "$BASIS", "$FMOPRP"]:
        assert text.count(block) == 1, f"{block} appears more than once"


def test_patch_inp_fmo_block_preserved(raw_inp, tmp_path):
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.postprocess import patch_inp

    cfg = FragitConfig(central_fragment_resname="LIG", mwords=100, ngroup=1, calc_mode="2layer")
    out = tmp_path / "patched.inp"
    patch_inp(raw_inp, cfg, output_path=out)
    text = out.read_text()

    assert "NFRAG=3" in text
    assert "NLAYER=2" in text
    assert "MPLEVL(1)=0,2" in text
    assert "LIG001" in text


def test_patch_inp_replaces_resdim(raw_inp, tmp_path):
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.postprocess import patch_inp

    cfg = FragitConfig(central_fragment_resname="LIG", mwords=100, ngroup=1,
                       calc_mode="2layer", resdim=4.5)
    out = tmp_path / "patched.inp"
    patch_inp(raw_inp, cfg, output_path=out)
    assert "RESDIM=4.5" in out.read_text()


def test_patch_inp_basis_sto3g(raw_inp, tmp_path):
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.postprocess import patch_inp

    cfg = FragitConfig(central_fragment_resname="LIG", mwords=50, ngroup=1,
                       calc_mode="hf", basis="STO-3G")
    out = tmp_path / "patched.inp"
    patch_inp(raw_inp, cfg, output_path=out)
    assert " $BASIS GBASIS=STO NGAUSS=3 $END" in out.read_text()


def test_patch_inp_raises_on_missing_fmo(tmp_path):
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.postprocess import patch_inp

    no_fmo = tmp_path / "no_fmo.inp"
    no_fmo.write_text(" $SYSTEM MWORDS=50 $END\n $BASIS GBASIS=N31 $END\n")
    cfg = FragitConfig(central_fragment_resname="LIG", mwords=50, ngroup=1)
    with pytest.raises(ValueError, match="No .FMO block"):
        patch_inp(no_fmo, cfg)


def test_patch_inp_overwrites_in_place(raw_inp):
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.postprocess import patch_inp

    cfg = FragitConfig(central_fragment_resname="LIG", mwords=100, ngroup=1, calc_mode="2layer")
    result = patch_inp(raw_inp, cfg)
    assert result == raw_inp
    assert "$SYSTEM MWORDS=100" in raw_inp.read_text()
