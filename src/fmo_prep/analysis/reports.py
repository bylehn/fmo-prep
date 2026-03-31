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
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = "delta_H_summary.txt" if is_delta else "interactions_summary.txt"
    out_path = output_dir / filename

    overall_total = net_totals.sum()

    if is_delta:
        energy_header = "[+] Total Delta Interaction Energy (System 2 - System 1) per Protein Residue"
        if mode == "chain":
            energy_header += " vs Peptide:"
            overall_label = f"OVERALL MACROSCOPIC \u0394H (Protein-Peptide, System 2 - System 1): {overall_total:>10.3f} kcal/mol"
        else:
            energy_header += " vs Ligand:"
            overall_label = f"OVERALL MACROSCOPIC \u0394H (Protein-Ligand, System 2 - System 1): {overall_total:>10.3f} kcal/mol"
    else:
        if mode == "chain":
            energy_header = "[+] Total Interaction Energy per Protein Residue vs Peptide:"
            overall_label = f"OVERALL TOTAL PROTEIN-PEPTIDE INTERACTION ENERGY: {overall_total:>10.3f} kcal/mol"
        else:
            energy_header = "[+] Total Interaction Energy per Protein Residue vs Ligand:"
            overall_label = f"OVERALL TOTAL PROTEIN-LIGAND INTERACTION ENERGY: {overall_total:>10.3f} kcal/mol"

    with open(out_path, "w") as f:
        if cmd_str:
            f.write(f"[+] Command run: {cmd_str}\n")
        f.write("=" * 50 + "\n")
        if mode == "chain":
            f.write(f"[+] Partner definition: {partner_label}\n")
        else:
            f.write(f"[+] Ligand fragment(s): {partner_label}\n")
        f.write(energy_header + "\n")
        f.write("-" * 50 + "\n")

        for label, val in net_totals.items():
            if abs(val) > 0.01:
                f.write(f"{label:<15}: {val:>10.3f} kcal/mol\n")

        f.write("=" * 50 + "\n")
        f.write(overall_label + "\n")
        f.write("=" * 50 + "\n")

        if is_delta:
            if overall_total < 0:
                f.write("--> System 2 binds MORE strongly/favorably than System 1.\n")
            elif overall_total > 0:
                f.write("--> System 2 binds LESS strongly/favorably than System 1.\n")

    return out_path
