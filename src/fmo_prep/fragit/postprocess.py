"""GAMESS input file postprocessing.

Automates the manual edits previously applied to raw FragIt .inp output:

1. Prepend standard header blocks ($SYSTEM, $GDDI, $SCF, $CONTRL, $FMOPRP)
   before the $BASIS line.
2. Insert RESDIM and RCORSD into the $FMO block after the NBODY line.
3. If mp2_level=True and NLAYER/MPLEVL are absent: insert them after RCORSD.

Reference: cdk2_benchmark/data1/step2_fragit/minimised_complex.inp
"""

from __future__ import annotations

from pathlib import Path

from fmo_prep.config import FragitConfig


_HEADER_TEMPLATE = """\
 $SYSTEM MWORDS={mwords} $END
 $GDDI NGROUP={ngroup} $END
 $SCF CONV=1E-6 DIRSCF=.T. NPUNCH=0 DIIS=.T. SOSCF=.F. $END
 $CONTRL NPRINT=-5 ISPHER=1 MAXIT=100
         RUNTYP=ENERGY
 $END
 $FMOPRP NPRINT=9 NGUESS=2 MAXIT=100 $END
"""


def patch_inp(inp_path: Path, cfg: FragitConfig, output_path: Path | None = None) -> Path:
    """Patch a raw FragIt-generated GAMESS .inp file.

    Steps applied:
    - Prepend header blocks before $BASIS
    - Add RESDIM / RCORSD to $FMO block (after NBODY)
    - Add NLAYER / MPLEVL to $FMO block if mp2_level=True and absent

    Args:
        inp_path: Path to the raw FragIt .inp file.
        cfg: FragitConfig supplying RESDIM, RCORSD, mp2_level, mwords, ngroup.
        output_path: Destination path. Defaults to overwriting inp_path.

    Returns:
        Path to the patched file.
    """
    raise NotImplementedError("fragit/postprocess.patch_inp not yet implemented")
