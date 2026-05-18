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
    assert "centralfragment = 42" in text
    assert "peptide_methylated" not in text


def test_render_config_mp2_full(tmp_path):
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.runner import render_config

    cfg = FragitConfig(central_fragment_resname="LIG", calc_mode="mp2")
    out = render_config(cfg, central_fragment_id=0, output_path=tmp_path / "fragit.ini")
    text = out.read_text()

    # mp2/2layer settings are handled by postprocess.py, not the FragIt ini
    assert "dohop = True" in text
    assert "efpwaters = 0" in text


def test_render_config_hf(tmp_path):
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.runner import render_config

    cfg = FragitConfig(central_fragment_resname="LIG", calc_mode="hf")
    out = render_config(cfg, central_fragment_id=0, output_path=tmp_path / "fragit.ini")
    text = out.read_text()

    assert "dohop = True" in text


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


def test_render_config_protein_peptide(tmp_path):
    from fmo_prep.config import FragitConfig
    from fmo_prep.fragit.runner import render_config

    cfg = FragitConfig(central_fragment_resname="B", calc_mode="hf")
    out = render_config(cfg, central_fragment_id=0, output_path=tmp_path / "fragit.ini",
                        system_type="protein_peptide")
    text = out.read_text()

    assert "peptide_methylated" in text
    assert "nterminal" not in text


def _make_coord_fixtures(tmp_path, frags):
    """Build a minimal .inp (with $FMOXYZ + $FMO) and a matching .pdb.

    frags: list of (resname, chain, resid, coords)
    coords: list of (x, y, z) tuples — one per atom in the fragment
    """
    all_atoms = []
    for resname, chain, resid, coords in frags:
        for x, y, z in coords:
            all_atoms.append((resname, chain, resid, x, y, z))

    # Build PDB
    pdb_lines = []
    for i, (resname, chain, resid, x, y, z) in enumerate(all_atoms, start=1):
        record = "HETATM" if resname not in ("ALA", "GLY", "LYS") else "ATOM  "
        pdb_lines.append(
            f"{record}{i:5d}  CA  {resname:<3s} {chain}{resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
        )
    pdb_lines.append("END")
    pdb = tmp_path / "test.pdb"
    pdb.write_text("\n".join(pdb_lines) + "\n")

    # Build INDAT ranges and FMOXYZ
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
        f" $FMO\n"
        f"      NFRAG={nfrag}\n"
        f"      NBODY=2\n"
        f"      ICHARG(1)=" + ",".join("  0" for _ in frags) + "\n"
        f"      FRGNAM(1)= {frgnam}\n"
        f"      INDAT(1)=0\n"
        + "\n".join(indat_lines) + "\n"
        f" $END\n"
        f" $FMOXYZ\n"
        + "\n".join(fmoxyz_lines) + "\n"
        f" $END\n"
    )
    inp = tmp_path / "test.inp"
    inp.write_text(inp_text)
    return inp, pdb


def test_find_central_fragment_id(tmp_path):
    from fmo_prep.fragit.runner import find_central_fragment_id

    inp, pdb = _make_coord_fixtures(tmp_path, [
        ("ALA", "A", 1, [(1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]),
        ("GLY", "A", 2, [(3.0, 0.0, 0.0), (4.0, 0.0, 0.0)]),
        ("LIG", "A", 3, [(5.0, 0.0, 0.0), (6.0, 0.0, 0.0)]),
    ])
    fid = find_central_fragment_id(inp, pdb, "LIG")
    assert fid == 3


def test_find_central_fragment_id_with_chain(tmp_path):
    from fmo_prep.fragit.runner import find_central_fragment_id

    inp, pdb = _make_coord_fixtures(tmp_path, [
        ("ALA", "A", 1, [(1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]),
        ("LIG", "A", 2, [(3.0, 0.0, 0.0), (4.0, 0.0, 0.0)]),
        ("LIG", "B", 3, [(5.0, 0.0, 0.0), (6.0, 0.0, 0.0)]),
    ])
    fid = find_central_fragment_id(inp, pdb, "LIG", chain="B")
    assert fid == 3


def test_find_central_fragment_id_not_found(tmp_path):
    from fmo_prep.fragit.runner import find_central_fragment_id

    inp, pdb = _make_coord_fixtures(tmp_path, [
        ("ALA", "A", 1, [(1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]),
        ("GLY", "A", 2, [(3.0, 0.0, 0.0), (4.0, 0.0, 0.0)]),
    ])
    with pytest.raises(ValueError, match="No fragment"):
        find_central_fragment_id(inp, pdb, "LIG")
