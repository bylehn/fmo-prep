"""GAMESS input file postprocessing.

FragIt writes a complete GAMESS input including all header blocks:
  $SYSTEM, $GDDI, $SCF, $CONTRL, $BASIS, $FMOPRP, $FMO, $FMOBND, $DATA,
  $FMOHYB, $FMOXYZ

This module replaces those header blocks with mode-appropriate versions
and optionally inserts a $PCM block for implicit solvent.

Calculation modes (cfg.calc_mode):

  hf      - HF only, no MP2.
             SCF: CONV=1E-6, DIIS=.F., SOSCF=.T.
             CONTRL: no MAXIT, no SCFTYP

  mp2     - MP2 on entire system (for PIEDA analysis).
             SCF: CONV=1E-7, DIIS=.F., SOSCF=.T.
             CONTRL: no MAXIT, no SCFTYP
             FMOPRP: PRTDST(1)=100.0,0.5,0.6,0.0 IPIEDA=2

  2layer  - MP2 at active site (layer 2), HF elsewhere (layer 1).
             SCF: CONV=1E-6, DIIS=.T., SOSCF=.F.
             CONTRL: MAXIT=100 SCFTYP=RHF
             FMOPRP: MAXIT=100

Implicit solvent (cfg.implicit_solvent=True):
  - Inserts $PCM SOLVNT=WATER IFMO=1 ICOMP=0 $END after $BASIS
  - Changes FMOPRP to IPIEDA=1 instead of mode default
  - SCF: CONV=1E-6, DIIS=.F., SOSCF=.T. (overrides mp2 mode)
"""

from __future__ import annotations

import re
from pathlib import Path

from fmo_prep.config import FragitConfig

# Blocks written by FragIt that we strip and replace
_STRIP_PATTERNS = [
    r"^ \$SYSTEM\b.*?\$END\n",
    r"^ \$GDDI\b.*?\$END\n",
    r"^ \$SCF\b.*?\$END\n",
    r"^ \$CONTRL\b.*?\$END\n",   # may span 2 lines — handled with DOTALL
    r"^ \$BASIS\b.*?\$END\n",
    r"^ \$FMOPRP\b.*?\$END\n",
    r"^ \$PCM\b.*?\$END\n",
]


def _build_scf(cfg: FragitConfig) -> str:
    if cfg.calc_mode == "2layer":
        return " $SCF CONV=1E-6 DIRSCF=.T. NPUNCH=0 DIIS=.T. SOSCF=.F. $END\n"
    elif cfg.calc_mode == "mp2" and not cfg.implicit_solvent:
        return " $SCF CONV=1E-7 DIRSCF=.T. NPUNCH=0 DIIS=.F. SOSCF=.T. $END\n"
    else:  # hf, or mp2/2layer with implicit_solvent
        return " $SCF CONV=1E-6 DIRSCF=.T. NPUNCH=0 DIIS=.F. SOSCF=.T. $END\n"


def _build_contrl(cfg: FragitConfig) -> str:
    if cfg.calc_mode == "2layer":
        return (
            " $CONTRL NPRINT=-5 ISPHER=1 MAXIT=100\n"
            "         RUNTYP=ENERGY SCFTYP=RHF\n"
            " $END\n"
        )
    else:
        return (
            " $CONTRL NPRINT=-5 ISPHER=1\n"
            "         RUNTYP=ENERGY\n"
            " $END\n"
        )


def _build_fmoprp(cfg: FragitConfig) -> str:
    if cfg.implicit_solvent:
        return " $FMOPRP NPRINT=9 NGUESS=2 IPIEDA=1 $END\n"
    elif cfg.calc_mode == "2layer":
        return " $FMOPRP NPRINT=9 NGUESS=2 MAXIT=100 $END\n"
    elif cfg.calc_mode == "mp2":
        return " $FMOPRP NPRINT=9 NGUESS=2 PRTDST(1)=100.0,0.5,0.6,0.0 IPIEDA=2 $END\n"
    else:  # hf
        return " $FMOPRP NPRINT=9 NGUESS=2 $END\n"


def _build_basis(cfg: FragitConfig) -> str:
    basis_map = {
        "6-31G*":   "GBASIS=N31 NGAUSS=6 NDFUNC=1",
        "6-31G(d)": "GBASIS=N31 NGAUSS=6 NDFUNC=1",
        "6-31G":    "GBASIS=N31 NGAUSS=6",
        "STO-3G":   "GBASIS=STO NGAUSS=3",
        "3-21G":    "GBASIS=N21 NGAUSS=3",
    }
    gbasis = basis_map.get(cfg.basis, "GBASIS=N31 NGAUSS=6 NDFUNC=1")
    return f" $BASIS {gbasis} $END\n"


def _build_pcm() -> str:
    return " $PCM SOLVNT=WATER IFMO=1 ICOMP=0 $END\n"


def patch_inp(inp_path: Path, cfg: FragitConfig, output_path: Path | None = None) -> Path:
    """Patch a FragIt-generated GAMESS .inp file with mode-appropriate header blocks.

    Steps applied:
    1. Strip FragIt's $SYSTEM, $GDDI, $SCF, $CONTRL, $BASIS, $FMOPRP (and $PCM if present).
    2. Prepend our versions of those blocks, chosen based on cfg.calc_mode.
    3. Insert $PCM block after $BASIS when cfg.implicit_solvent=True.
    4. Replace RESDIM/RCORSD in $FMO if non-default values are configured.

    NLAYER and MPLEVL are written correctly by FragIt from the .ini config.

    Args:
        inp_path: Path to the FragIt-generated .inp file.
        cfg: FragitConfig supplying all GAMESS header parameters.
        output_path: Destination path. Defaults to overwriting inp_path.

    Returns:
        Path to the patched file.

    Raises:
        ValueError: If the file contains no $FMO block.
    """
    inp_path = Path(inp_path)
    output_path = Path(output_path) if output_path else inp_path

    text = inp_path.read_text()

    if "$FMO" not in text:
        raise ValueError(f"No $FMO block found in {inp_path} — is this a valid FragIt .inp?")

    # --- Step 1: strip existing header blocks ---
    # $CONTRL may span two lines — handle with DOTALL first
    text = re.sub(r"^ \$CONTRL\b.*?\$END\n", "", text, flags=re.MULTILINE | re.DOTALL)
    # Strip all other single-line blocks (skip index 3 = $CONTRL, already handled)
    for i, pattern in enumerate(_STRIP_PATTERNS):
        if i == 3:
            continue
        text = re.sub(pattern, "", text, flags=re.MULTILINE)

    # --- Step 2: prepend standardised header ---
    header = (
        f" $SYSTEM MWORDS={cfg.mwords} $END\n"
        f" $GDDI NGROUP={cfg.ngroup} $END\n"
        + _build_scf(cfg)
        + _build_contrl(cfg)
        + _build_basis(cfg)
        + (_build_pcm() if cfg.implicit_solvent else "")
        + _build_fmoprp(cfg)
    )
    text = header + text.lstrip("\n")

    # --- Steps 3–4: patch RESDIM / RCORSD if non-default ---
    if cfg.resdim != 2.0:
        text = re.sub(
            r"(^\s*RESDIM\s*=\s*)\S+",
            lambda m: m.group(1) + str(cfg.resdim),
            text, flags=re.MULTILINE,
        )
    if cfg.rcorsd != 2.0:
        text = re.sub(
            r"(^\s*RCORSD\s*=\s*)\S+",
            lambda m: m.group(1) + str(cfg.rcorsd),
            text, flags=re.MULTILINE,
        )

    output_path.write_text(text)
    return output_path
