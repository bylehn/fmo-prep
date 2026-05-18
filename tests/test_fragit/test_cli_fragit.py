"""Tests for the `fmo-prep fragit` CLI command and FragitOnlyConfig."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest
import yaml


# ---------------------------------------------------------------------------
# FragitOnlyConfig
# ---------------------------------------------------------------------------

def test_fragit_only_config_minimal(tmp_path):
    from fmo_prep.config import FragitOnlyConfig

    cfg_yaml = tmp_path / "config.yaml"
    cfg_yaml.write_text("pdb: ./complex.pdb\noutput_dir: ./out\n")
    cfg = FragitOnlyConfig.from_yaml(cfg_yaml)

    assert cfg.pdb == Path("./complex.pdb")
    assert cfg.output_dir == Path("./out")
    assert cfg.fragit.basis == "6-31G*"
    assert cfg.fragit.calc_mode == "mp2"
    assert cfg.fragit.fmo_level == 2
    assert cfg.fragit.charge_model == "MMFF94"
    assert cfg.fragit.maxfragsize == 100


def test_fragit_only_config_full(tmp_path):
    from fmo_prep.config import FragitOnlyConfig

    raw = {
        "pdb": "./complex.pdb",
        "output_dir": "./out",
        "fragit": {
            "basis": "STO-3G",
            "calc_mode": "2layer",
            "fmo_level": "low",
            "charge_model": "formal",
            "central_fragment_resname": "LIG",
            "ligand_chain": "B",
            "maxfragsize": 50,
        },
    }
    cfg_yaml = tmp_path / "config.yaml"
    cfg_yaml.write_text(yaml.dump(raw))
    cfg = FragitOnlyConfig.from_yaml(cfg_yaml)

    assert cfg.fragit.basis == "STO-3G"
    assert cfg.fragit.fmo_level == 1  # "low" → 1
    assert cfg.fragit.charge_model == "formal"
    assert cfg.fragit.ligand_chain == "B"
    assert cfg.fragit.maxfragsize == 50


def test_fragit_config_fmo_level_high():
    from fmo_prep.config import FragitConfig

    cfg = FragitConfig(fmo_level="high")
    assert cfg.fmo_level == 2

    cfg2 = FragitConfig(fmo_level="low")
    assert cfg2.fmo_level == 1


# ---------------------------------------------------------------------------
# CLI `fragit` command (subprocess mocked)
# ---------------------------------------------------------------------------

def _make_mock_inp(path: Path, fmoxyz: bool = False) -> None:
    """Write a minimal patched .inp (with optional $FMOXYZ for all 4 atoms) for testing.

    Matches the PDB used in test_cli_fragit_mp2: ALA (atoms 1-2) + LIG (atoms 3-4).
    Coords: (1,0,0), (2,0,0) for ALA; (3,0,0), (4,0,0) for LIG.
    """
    fmoxyz_block = (
        " $FMOXYZ\n"
        "C    6.0       1.000     0.000     0.000\n"
        "C    6.0       2.000     0.000     0.000\n"
        "C    6.0       3.000     0.000     0.000\n"
        "O    8.0       4.000     0.000     0.000\n"
        " $END\n"
    ) if fmoxyz else ""
    path.write_text(
        " $SYSTEM MWORDS=125 $END\n"
        " $GDDI NGROUP=1 $END\n"
        " $SCF CONV=1E-7 DIRSCF=.T. NPUNCH=0 DIIS=.F. SOSCF=.T. $END\n"
        " $CONTRL NPRINT=-5 ISPHER=1\n"
        "         RUNTYP=ENERGY\n"
        " $END\n"
        " $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END\n"
        " $FMOPRP NPRINT=9 NGUESS=2 IPIEDA=2 $END\n"
        " $FMO\n"
        "      NFRAG=2\n"
        "      NBODY=2\n"
        "      NLAYER=1\n"
        "      MPLEVL(1)=2\n"
        "      ICHARG(1)=  0,  0\n"
        "      FRGNAM(1)= ALA001,  LIG002\n"
        "      INDAT(1)=0\n"
        "            1    -2      0\n"
        "            3    -4      0\n"
        " $END\n"
        + fmoxyz_block
    )


def test_cli_fragit_mp2(tmp_path):
    """fmo-prep fragit runs for mp2 mode (single FragIt pass)."""
    from click.testing import CliRunner
    from fmo_prep.cli import cli

    pdb = tmp_path / "complex.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  CB  ALA A   1       2.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    3  C1  LIG A   2       3.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    4  O1  LIG A   2       4.000   0.000   0.000  1.00  0.00           O\n"
        "END\n"
    )
    cfg_yaml = tmp_path / "config.yaml"
    cfg_yaml.write_text(
        f"pdb: {pdb}\noutput_dir: {tmp_path / 'out'}\n"
        "fragit:\n  calc_mode: mp2\n"
    )
    out_dir = tmp_path / "out"
    raw_inp = out_dir / "fragit" / "complex.inp"
    final_inp = out_dir / "fmo_run.inp"

    def fake_run(cmd, **kw):
        raw_inp.parent.mkdir(parents=True, exist_ok=True)
        _make_mock_inp(raw_inp, fmoxyz=True)
        result = MagicMock()
        result.returncode = 0
        result.stdout = ""
        return result

    with patch("subprocess.run", side_effect=fake_run):
        runner = CliRunner()
        result = runner.invoke(cli, ["fragit", "-c", str(cfg_yaml)])

    assert result.exit_code == 0, result.output
    assert final_inp.exists()
    assert (out_dir / "fragment_map.json").exists()


def test_cli_fragit_2layer_requires_resname(tmp_path):
    """fmo-prep fragit exits with UsageError when 2layer but no resname."""
    from click.testing import CliRunner
    from fmo_prep.cli import cli

    pdb = tmp_path / "complex.pdb"
    pdb.write_text("END\n")
    cfg_yaml = tmp_path / "config.yaml"
    cfg_yaml.write_text(
        f"pdb: {pdb}\noutput_dir: {tmp_path / 'out'}\n"
        "fragit:\n  calc_mode: 2layer\n"
    )

    def fake_run(cmd, **kw):
        inp_path = tmp_path / "out" / "fragit" / "complex.inp"
        inp_path.parent.mkdir(parents=True, exist_ok=True)
        _make_mock_inp(inp_path, fmoxyz=True)
        r = MagicMock()
        r.returncode = 0
        r.stdout = ""
        return r

    with patch("subprocess.run", side_effect=fake_run):
        runner = CliRunner()
        result = runner.invoke(cli, ["fragit", "-c", str(cfg_yaml)])

    assert result.exit_code != 0
    assert "central_fragment_resname" in result.output
