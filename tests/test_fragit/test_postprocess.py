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
      ICHARG(1)=  0,  0,  0
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
    assert "IPIEDA=2" in text
    assert "PRTDST" not in text
    assert "$PCM" not in text
    assert "NLAYER=1" in text
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
    assert "IPIEDA=2" in text
    assert "PRTDST" not in text
    assert "$PCM" not in text
    assert "NLAYER=1" in text
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


# --- ICHARG auto-correction ---

def _make_inp(tmp_path, icharg_line, frgnam_line, fmobnd_lines=""):
    """Helper: write a minimal .inp with custom ICHARG, FRGNAM, and optional FMOBND.

    INDAT has 3 fragments: atoms 1-5, 6-10, 11-15.
    Default FMOBND (empty) means no backbone cuts → no artifact correction.
    Pass fmobnd_lines to add cuts, e.g. '      -4        6 6-31G* 6-31G*'
    which means BDA=4 (in frag 0, atoms 1-5) → BAA=6 (first atom of frag 1).
    """
    fmobnd_block = f" $FMOBND\n{fmobnd_lines}\n $END\n" if fmobnd_lines else ""
    content = f"""\
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
      {icharg_line}
      {frgnam_line}
      INDAT(1)=0
            1    -5      0
            6    -10     0
           11    -15     0
 $END
{fmobnd_block} $DATA
 title
 C1
...
 $END
"""
    p = tmp_path / "icharg.inp"
    p.write_text(content)
    return p


def test_fix_fragment_charges_nonstandard_next_with_cut_corrected(tmp_path):
    """Fragment with non-standard next residue AND FMOBND cut gets charge zeroed."""
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.postprocess import patch_inp

    # PHE014 (frag 1, atoms 6-10) has charge -1; next frag PEX015 is non-standard.
    # FMOBND cut: BDA=9 (in frag 1, atoms 6-10) → BAA=11 (first atom of frag 2).
    inp = _make_inp(
        tmp_path,
        "ICHARG(1)=  0, -1,  0",
        "FRGNAM(1)= ALA013,  PHE014,  PEX015",
        fmobnd_lines="      -9       11 6-31G* 6-31G*",
    )
    cfg = FragitConfig(central_fragment_resname="LIG", mwords=50, ngroup=1, calc_mode="hf")
    out = tmp_path / "fixed.inp"
    patch_inp(inp, cfg, output_path=out)

    assert "ICHARG(1)=  0,  0,  0" in out.read_text()


def test_fix_fragment_charges_nonstandard_next_no_cut_unchanged(tmp_path):
    """Fragment before a non-standard ligand (no FMOBND cut) is left untouched."""
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.postprocess import patch_inp

    # ARG (frag 1, atoms 6-10) has charge -1; next frag LZ1 is non-standard ligand.
    # No FMOBND cut between frag 1 and frag 2 → ligand is a separate complete fragment.
    inp = _make_inp(
        tmp_path,
        "ICHARG(1)=  0, -1,  0",
        "FRGNAM(1)= ALA013,  ARG014,  LZ1015",
        fmobnd_lines="",  # no cut
    )
    cfg = FragitConfig(central_fragment_resname="LIG", mwords=50, ngroup=1, calc_mode="hf")
    out = tmp_path / "fixed.inp"
    patch_inp(inp, cfg, output_path=out)

    assert "ICHARG(1)=  0, -1,  0" in out.read_text()


def test_fix_fragment_charges_standard_next_unchanged(tmp_path):
    """Fragment before a standard charged residue (LYS) is left untouched."""
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.postprocess import patch_inp

    # GLN013 (frag 0) has charge +1 because LYS014 is the inner residue — legitimate.
    inp = _make_inp(
        tmp_path,
        "ICHARG(1)=  1,  0,  0",
        "FRGNAM(1)= GLN013,  LYS014,  ALA015",
        fmobnd_lines="      -4        6 6-31G* 6-31G*",
    )
    cfg = FragitConfig(central_fragment_resname="LIG", mwords=50, ngroup=1, calc_mode="hf")
    out = tmp_path / "fixed.inp"
    patch_inp(inp, cfg, output_path=out)

    assert "ICHARG(1)=  1,  0,  0" in out.read_text()


def test_fix_fragment_charges_no_artifacts_unchanged(tmp_path):
    """No correction when all non-zero charges have standard next residues."""
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.postprocess import patch_inp

    inp = _make_inp(
        tmp_path,
        "ICHARG(1)=  0,  1, -1",
        "FRGNAM(1)= ALA013,  LYS014,  GLU015",
        fmobnd_lines="      -4        6 6-31G* 6-31G*\n      -9       11 6-31G* 6-31G*",
    )
    cfg = FragitConfig(central_fragment_resname="LIG", mwords=50, ngroup=1, calc_mode="hf")
    out = tmp_path / "fixed.inp"
    patch_inp(inp, cfg, output_path=out)

    assert "ICHARG(1)=  0,  1, -1" in out.read_text()


# --- Pass 2: impossible |charge| > 1 clamping ---

def test_pass2_clamps_interior_standard_residue(tmp_path):
    """Interior ALA with charge +2 is clamped to +1 (MMFF94 rounding artifact)."""
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.postprocess import patch_inp

    # ALA014 (frag 1, atoms 6-10): charge +2, has both incoming (BAA=6) and outgoing (BDA=9)
    inp = _make_inp(
        tmp_path,
        "ICHARG(1)=  0,  2,  0",
        "FRGNAM(1)= GLN013,  ALA014,  VAL015",
        fmobnd_lines="      -4        6 6-31G* 6-31G*\n      -9       11 6-31G* 6-31G*",
    )
    cfg = FragitConfig(mwords=50, ngroup=1, calc_mode="hf")
    out = tmp_path / "clamped.inp"
    patch_inp(inp, cfg, output_path=out)
    assert "ICHARG(1)=  0,  1,  0" in out.read_text()


def test_pass2_skips_chain_terminus(tmp_path):
    """Chain-terminal fragment with |charge| > 1 is NOT clamped."""
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.postprocess import patch_inp

    # GLN013 (frag 0, atoms 1-5): charge +2, only outgoing cut (no incoming) → terminus
    inp = _make_inp(
        tmp_path,
        "ICHARG(1)=  2,  0,  0",
        "FRGNAM(1)= GLN013,  ALA014,  VAL015",
        fmobnd_lines="      -4        6 6-31G* 6-31G*",
    )
    cfg = FragitConfig(mwords=50, ngroup=1, calc_mode="hf")
    out = tmp_path / "terminus.inp"
    patch_inp(inp, cfg, output_path=out)
    assert "ICHARG(1)=  2,  0,  0" in out.read_text()


def test_pass2_skips_phosphorylated_residue(tmp_path):
    """SEP (phosphorylated) with charge -2 is NOT clamped."""
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.postprocess import patch_inp

    # SEP014 (frag 1, atoms 6-10): charge -2, interior — but SEP is in _PHOSPHORYLATED_RESIDUES
    inp = _make_inp(
        tmp_path,
        "ICHARG(1)=  0, -2,  0",
        "FRGNAM(1)= GLN013,  SEP014,  VAL015",
        fmobnd_lines="      -4        6 6-31G* 6-31G*\n      -9       11 6-31G* 6-31G*",
    )
    cfg = FragitConfig(mwords=50, ngroup=1, calc_mode="hf")
    out = tmp_path / "phos.inp"
    patch_inp(inp, cfg, output_path=out)
    assert "ICHARG(1)=  0, -2,  0" in out.read_text()


# --- build_fragment_residue_map / find_fragment_by_chain_resname ---

def _make_fmoxyz_fixtures(tmp_path, frags):
    """Build a minimal patched .inp (with $FMOXYZ) and a matching .pdb.

    frags: list of (resname, chain, resid, [(x,y,z), ...])
    """
    all_atoms = []
    for resname, chain, resid, coords in frags:
        for x, y, z in coords:
            all_atoms.append((resname, chain, resid, x, y, z))

    pdb_lines = []
    for i, (resname, chain, resid, x, y, z) in enumerate(all_atoms, start=1):
        record = "HETATM" if resname not in ("ALA", "GLY", "LYS", "ARG") else "ATOM  "
        pdb_lines.append(
            f"{record}{i:5d}  CA  {resname:<3s} {chain}{resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
        )
    pdb_lines.append("END")
    pdb = tmp_path / "test.pdb"
    pdb.write_text("\n".join(pdb_lines) + "\n")

    indat_lines = []
    fmoxyz_lines = []
    atom_idx = 1
    for resname, chain, resid, coords in frags:
        n = len(coords)
        indat_lines.append(f"      {atom_idx:5d} {-(atom_idx + n - 1):6d}     0")
        for x, y, z in coords:
            fmoxyz_lines.append(f"C    6.0  {x:10.3f}{y:10.3f}{z:10.3f}")
        atom_idx += n

    nfrag = len(frags)
    frgnam = ",  ".join(f"{r[0]}{r[2]:03d}" for r in frags)
    inp_text = (
        f" $SYSTEM MWORDS=125 $END\n"
        f" $GDDI NGROUP=1 $END\n"
        f" $SCF CONV=1E-6 DIRSCF=.T. NPUNCH=0 DIIS=.F. SOSCF=.T. $END\n"
        f" $CONTRL NPRINT=-5 ISPHER=1\n"
        f"         RUNTYP=ENERGY\n"
        f" $END\n"
        f" $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END\n"
        f" $FMOPRP NPRINT=9 NGUESS=2 IPIEDA=2 $END\n"
        f" $FMO\n"
        f"      NFRAG={nfrag}\n"
        f"      NBODY=2\n"
        f"      NLAYER=1\n"
        f"      MPLEVL(1)=2\n"
        f"      ICHARG(1)=" + ",".join("  0" for _ in frags) + "\n"
        f"      FRGNAM(1)= {frgnam}\n"
        f"      INDAT(1)=0\n"
        + "\n".join(indat_lines) + "\n"
        f" $END\n"
        f" $FMOXYZ\n"
        + "\n".join(fmoxyz_lines) + "\n"
        f" $END\n"
    )
    inp = tmp_path / "patched.inp"
    inp.write_text(inp_text)
    return inp, pdb


def test_build_fragment_residue_map_basic(tmp_path):
    from fmo_prep.fragit.postprocess import build_fragment_residue_map

    inp, pdb = _make_fmoxyz_fixtures(tmp_path, [
        ("ALA", "A", 1, [(1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]),
        ("GLY", "A", 2, [(3.0, 0.0, 0.0), (4.0, 0.0, 0.0)]),
        ("LIG", "A", 3, [(5.0, 0.0, 0.0), (6.0, 0.0, 0.0)]),
    ])
    fmap = build_fragment_residue_map(inp, pdb)

    assert len(fmap) == 3
    assert fmap[0]["fragment_index"] == 1
    assert fmap[0]["majority_residue"] == "ALA"
    assert fmap[2]["fragment_index"] == 3
    assert fmap[2]["majority_residue"] == "LIG"
    assert fmap[2]["chain"] == "A"
    # fragment_map.json written alongside inp
    json_path = inp.parent / "fragment_map.json"
    assert json_path.exists()


def test_find_fragment_by_chain_resname_no_chain(tmp_path):
    from fmo_prep.fragit.postprocess import build_fragment_residue_map, find_fragment_by_chain_resname

    inp, pdb = _make_fmoxyz_fixtures(tmp_path, [
        ("ALA", "A", 1, [(1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]),
        ("LIG", "A", 3, [(3.0, 0.0, 0.0), (4.0, 0.0, 0.0)]),
    ])
    fmap = build_fragment_residue_map(inp, pdb)
    assert find_fragment_by_chain_resname(fmap, None, "LIG") == 2


def test_find_fragment_by_chain_resname_with_chain(tmp_path):
    from fmo_prep.fragit.postprocess import build_fragment_residue_map, find_fragment_by_chain_resname

    inp, pdb = _make_fmoxyz_fixtures(tmp_path, [
        ("LIG", "A", 1, [(1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]),
        ("LIG", "B", 1, [(3.0, 0.0, 0.0), (4.0, 0.0, 0.0)]),
    ])
    fmap = build_fragment_residue_map(inp, pdb)
    assert find_fragment_by_chain_resname(fmap, "B", "LIG") == 2


def test_find_fragment_by_chain_resname_not_found(tmp_path):
    from fmo_prep.fragit.postprocess import build_fragment_residue_map, find_fragment_by_chain_resname

    inp, pdb = _make_fmoxyz_fixtures(tmp_path, [
        ("ALA", "A", 1, [(1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]),
    ])
    fmap = build_fragment_residue_map(inp, pdb)
    with pytest.raises(ValueError, match="No fragment"):
        find_fragment_by_chain_resname(fmap, None, "LIG")
