"""SDF file I/O helpers."""

from __future__ import annotations

from pathlib import Path


def load_mol(path: str | Path):
    """Load a molecule from SDF using RDKit."""
    from rdkit import Chem
    mol = Chem.SDMolSupplier(str(path), removeHs=False)[0]
    if mol is None:
        raise ValueError(f"Could not parse SDF: {path}")
    return mol
