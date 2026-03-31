"""GAMESS FMO log file parser.

Unified port of:
- fmo-poc/scripts/gamout.py      (PIEDA mode: ES, EX, CT, DI, SOL components)
- fmo-poc/scripts/gamout_nopieda.py (non-PIEDA mode: total pair energies only)

Main entry point::

    df = parse_gamout("fmo_run.log", pieda=True)
    df.to_csv("interactions.csv", index=False)
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def parse_gamout(log_path: str | Path, pieda: bool = True) -> pd.DataFrame:
    """Parse a GAMESS FMO output file and return a tidy DataFrame.

    Args:
        log_path: Path to the GAMESS .log file.
        pieda: If True, parse PIEDA decomposition (ES, EX, CT, DI, SOL).
               If False, parse total pair energies only.

    Returns:
        DataFrame with columns:
        - pieda=True:  I, IFRG, J, JFRG, R, Q, COMPONENT, ENERGY, TOTAL
        - pieda=False: I, IFRG, J, JFRG, R, Q, ENERGY
    """
    raise NotImplementedError("analysis/parser.parse_gamout not yet implemented")
