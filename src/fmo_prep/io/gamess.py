"""GAMESS .inp file I/O utilities.

Ported from:
- fmo-poc/scripts/map_fragments.py    (parse_inp_file, write_fragmapping)
- fmo-poc/inputs/fn001/debug_fragit.py (validate_fragments)
"""

from __future__ import annotations

from pathlib import Path


def parse_inp_file(inp_path: str | Path) -> tuple[list[str], list[list[int]]]:
    """Parse FRGNAM and INDAT sections from a GAMESS FMO .inp file.

    Args:
        inp_path: Path to the GAMESS .inp file.

    Returns:
        (frag_names, fragments) where:
        - frag_names: list of fragment name strings (e.g. ['MET001', 'GLU002', ...])
        - fragments: list of atom index lists (1-based), one per fragment
    """
    raise NotImplementedError("io/gamess.parse_inp_file not yet implemented")


def get_pdb_atoms(pdb_path: str | Path) -> list[dict]:
    """Read ATOM/HETATM records from a PDB file.

    Returns a list of dicts with keys: res_id, atom_name.
    Indices match 1-based GAMESS INDAT atom numbering (TER records ignored).
    """
    raise NotImplementedError("io/gamess.get_pdb_atoms not yet implemented")


def write_fragmapping(
    pdb_path: str | Path,
    inp_path: str | Path,
    output_path: str | Path,
) -> Path:
    """Write a human-readable fragment→residue mapping table.

    Columns: Frag | Name | Residues and Included Atoms

    Args:
        pdb_path: Path to the PDB file (atom order must match .inp).
        inp_path: Path to the GAMESS .inp file.
        output_path: Where to write the mapping text file.

    Returns:
        Path to the written file.
    """
    raise NotImplementedError("io/gamess.write_fragmapping not yet implemented")


def parse_frag_map_file(map_file: str | Path) -> dict:
    """Parse a complex_fragmapping.txt file into a frag_info dict.

    Args:
        map_file: Path to fragmapping file (pipe-delimited, produced by write_fragmapping).

    Returns:
        Dict mapping fragment_id (int) → {'chain': str, 'plot_label': str, 'name': str}
    """
    raise NotImplementedError("io/gamess.parse_frag_map_file not yet implemented")


def validate_fragments(inp_path: str | Path) -> None:
    """Print a fragment charge distribution summary and flag anomalies.

    Port of fmo-poc/inputs/fn001/debug_fragit.py validation logic.
    """
    raise NotImplementedError("io/gamess.validate_fragments not yet implemented")
