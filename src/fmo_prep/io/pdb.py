"""PDB file I/O helpers (thin wrappers / utilities)."""

from __future__ import annotations

from pathlib import Path


def load_structure(path: str | Path):
    """Load a PDB file as a parmed Structure."""
    import parmed as pmd
    return pmd.load_file(str(path))


def save_structure(structure, path: str | Path, overwrite: bool = True) -> Path:
    """Save a parmed Structure to a PDB file."""
    structure.save(str(path), overwrite=overwrite)
    return Path(path)
