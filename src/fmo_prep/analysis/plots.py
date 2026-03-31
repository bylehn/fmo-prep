"""Interaction energy visualisation.

Unified port of:
- fmo-poc/scripts/plot_pieda.py       (single-system plots)
- fmo-poc/scripts/plot_pieda_delta.py (delta / comparative plots)

Main entry point for the CLI::

    run_analysis(csv_path, compare_csv, frag_info, ...)
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd


def run_analysis(
    csv_path: Path,
    compare_csv: Optional[Path],
    frag_info: dict,
    ligand_hint: Optional[str],
    interaction_mode: str,
    significant_threshold: float,
    cov_threshold: float,
    output_dir: Path,
) -> None:
    """Top-level analysis runner called by the CLI.

    Generates stacked bar plot, heatmaps, and summary text file.

    Args:
        csv_path: Primary interaction CSV (from parse_gamout).
        compare_csv: Optional second CSV for delta analysis.
        frag_info: Fragment info dict from parse_frag_map_file.
        ligand_hint: Optional ligand label hint for auto-detection.
        interaction_mode: 'auto', 'chain', or 'ligand'.
        significant_threshold: Min |net energy| for bar plot inclusion.
        cov_threshold: Covalent outlier suppression threshold.
        output_dir: Where to write PNG and TXT outputs.
    """
    raise NotImplementedError("analysis/plots.run_analysis not yet implemented")


def process_csv(
    csv_file: str | Path,
    frag_info: dict,
    cov_threshold: float = 150.0,
    ligand_frag_ids=None,
    ligand_hint: Optional[str] = None,
    interaction_mode: str = "auto",
):
    """Process interaction CSV and return pivot DataFrames.

    Port of plot_pieda.py:process_csv().

    Returns:
        (pivot_df, heatmap_pivot, heatmap_all_pivot, name_map,
         ligand_frag_ids, mode_used, partner_label_str)
    """
    raise NotImplementedError("analysis/plots.process_csv not yet implemented")


def detect_ligand_fragments(df: pd.DataFrame, id_to_name: dict, ligand_hint: Optional[str] = None):
    """Auto-detect ligand fragment IDs from a PIEDA DataFrame.

    Port of plot_pieda.py:detect_ligand_fragments().
    """
    raise NotImplementedError("analysis/plots.detect_ligand_fragments not yet implemented")


def build_fragment_maps(df: pd.DataFrame, frag_info: dict):
    """Build fragment ID → name/plot-label maps.

    Port of plot_pieda.py:build_fragment_maps().
    """
    raise NotImplementedError("analysis/plots.build_fragment_maps not yet implemented")
