"""Summary report generation.

Writes per-residue interaction energy tables and macroscopic totals to text files.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def write_summary(
    net_totals: pd.Series,
    output_dir: Path,
    mode: str,
    partner_label: str,
    is_delta: bool = False,
    cmd_str: str = "",
) -> Path:
    """Write a plain-text summary of interaction energies.

    Args:
        net_totals: Series of (residue_label → net energy kcal/mol).
        output_dir: Where to write the summary file.
        mode: 'chain' or 'ligand'.
        partner_label: Display label for the partner (ligand name or 'Chain B').
        is_delta: If True, use delta-H labelling.
        cmd_str: Command string to record at the top of the file.

    Returns:
        Path to the written summary file.
    """
    raise NotImplementedError("analysis/reports.write_summary not yet implemented")
