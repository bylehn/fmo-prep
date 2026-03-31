"""GAMESS input file postprocessing.

FragIt writes a complete GAMESS input including all header blocks:
  $SYSTEM, $GDDI, $SCF, $CONTRL, $BASIS, $FMOPRP, $FMO, $FMOBND, $DATA,
  $FMOHYB, $FMOXYZ

This module replaces those header blocks with our standardised versions
(different SCF convergence, CONTRL settings, FMOPRP settings) and updates
RESDIM/RCORSD in the $FMO block when non-default values are requested.

The $FMO block body (NFRAG, NBODY, RESDIM, RCORSD, NLAYER, MPLEVL, ICHARG,
FRGNAM, INDAT, LAYER) is preserved unchanged — FragIt writes all of these
correctly based on the .ini config.

Reference for target header format:
    cdk2_benchmark/data1/step2_fragit/minimised_complex.inp
"""

from __future__ import annotations

import re
from pathlib import Path

from fmo_prep.config import FragitConfig

# Regex patterns for blocks we strip and replace
_STRIP_PATTERNS = [
    r"^ \$SYSTEM\b.*?\$END\n",   # $SYSTEM ... $END  (single line)
    r"^ \$GDDI\b.*?\$END\n",     # $GDDI   ... $END  (single line)
    r"^ \$SCF\b.*?\$END\n",      # $SCF    ... $END  (single line)
    r"^ \$CONTRL\b.*?\$END\n",   # $CONTRL ... $END  (may span 2 lines)
    r"^ \$BASIS\b.*?\$END\n",    # $BASIS  ... $END  (single line)
    r"^ \$FMOPRP\b.*?\$END\n",   # $FMOPRP ... $END  (single line)
]


def _build_header(cfg: FragitConfig) -> str:
    return (
        f" $SYSTEM MWORDS={cfg.mwords} $END\n"
        f" $GDDI NGROUP={cfg.ngroup} $END\n"
        f" $SCF CONV=1E-6 DIRSCF=.T. NPUNCH=0 DIIS=.T. SOSCF=.F. $END\n"
        f" $CONTRL NPRINT=-5 ISPHER=1 MAXIT=100\n"
        f"         RUNTYP=ENERGY\n"
        f" $END\n"
    )


def _build_basis(cfg: FragitConfig) -> str:
    # Translate config basis string to GAMESS GBASIS notation
    basis_map = {
        "6-31G*":    "GBASIS=N31 NGAUSS=6 NDFUNC=1",
        "6-31G(d)":  "GBASIS=N31 NGAUSS=6 NDFUNC=1",
        "6-31G":     "GBASIS=N31 NGAUSS=6",
        "STO-3G":    "GBASIS=STO NGAUSS=3",
        "3-21G":     "GBASIS=N21 NGAUSS=3",
    }
    gbasis = basis_map.get(cfg.basis, f"GBASIS=N31 NGAUSS=6 NDFUNC=1")
    return f" $BASIS {gbasis} $END\n"


def _build_fmoprp() -> str:
    return " $FMOPRP NPRINT=9 NGUESS=2 MAXIT=100 $END\n"


def patch_inp(inp_path: Path, cfg: FragitConfig, output_path: Path | None = None) -> Path:
    """Patch a FragIt-generated GAMESS .inp file with standardised header blocks.

    Steps applied:
    1. Strip FragIt's $SYSTEM, $GDDI, $SCF, $CONTRL, $BASIS, $FMOPRP lines.
    2. Prepend our versions of those blocks.
    3. If cfg.resdim != 2.0: replace RESDIM value in the $FMO block.
    4. If cfg.rcorsd != 2.0: replace RCORSD value in the $FMO block.

    NLAYER and MPLEVL are already written correctly by FragIt based on the
    .ini config, so no injection is needed.

    Args:
        inp_path: Path to the FragIt-generated .inp file.
        cfg: FragitConfig supplying mwords, ngroup, basis, resdim, rcorsd.
        output_path: Destination path. Defaults to overwriting inp_path.

    Returns:
        Path to the patched file.

    Raises:
        ValueError: If the file contains no $FMO line (not a valid FragIt output).
    """
    inp_path = Path(inp_path)
    output_path = Path(output_path) if output_path else inp_path

    text = inp_path.read_text()

    if "$FMO" not in text:
        raise ValueError(f"No $FMO block found in {inp_path} — is this a valid FragIt .inp?")

    # --- Step 1: strip old header blocks ---
    # $CONTRL may span two lines; handle multiline first
    text = re.sub(
        r"^ \$CONTRL\b.*?\$END\n",
        "",
        text,
        flags=re.MULTILINE | re.DOTALL,
    )
    # Remove remaining single-line blocks (skip index 3 = $CONTRL, already handled above)
    for i, pattern in enumerate(_STRIP_PATTERNS):
        if i == 3:
            continue  # $CONTRL already handled with DOTALL above
        text = re.sub(pattern, "", text, flags=re.MULTILINE)

    # --- Step 2: prepend standardised header ---
    header = _build_header(cfg) + _build_basis(cfg) + _build_fmoprp()
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
