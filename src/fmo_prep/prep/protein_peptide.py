"""Protein-peptide system preparation.

Workflow
--------
1. Load complex PDB with parmed
2. Remove water molecules
3. cut_and_cap: truncate to residues within cut_radius of the peptide,
   capping backbone cuts with H link atoms
4. Save capped PDB ready for FragIt

Source: PDBProgress tutorial notebooks + fmo-poc/structures/processed/cap_termini.py
Note: ACE/NME termini capping is deferred — H-link atoms only for now.
"""

from __future__ import annotations

import logging
from pathlib import Path

import parmed as pmd

from fmo_prep.config import RunConfig
from fmo_prep.prep.common import cut_and_cap, remove_water

logger = logging.getLogger(__name__)


def run(input_files: dict[str, str], output_dir: Path, cfg: RunConfig) -> Path:
    """Prepare a protein-peptide complex for FMO fragmentation.

    Args:
        input_files: Must contain 'complex' key with path to the full-system PDB.
        output_dir: Directory for intermediate and final outputs.
        cfg: Full run configuration (uses cfg.prep.cut_radius and
             cfg.fragit.central_fragment_resname).

    Returns:
        Path to the capped PDB file ready for FragIt.
    """
    complex_pdb = input_files.get("complex")
    if not complex_pdb:
        raise ValueError("Missing required input 'complex'")

    logger.info(f"Loading structure: {complex_pdb}")
    structure = pmd.load_file(complex_pdb)
    logger.info(f"Loaded {len(structure.atoms)} atoms, {len(structure.residues)} residues")

    logger.info("Removing water molecules...")
    structure = remove_water(structure)
    logger.info(f"After cleaning: {len(structure.atoms)} atoms, {len(structure.residues)} residues")

    selection = cfg.fragit.central_fragment_resname
    radius = cfg.prep.cut_radius
    logger.info(f"Truncating to {radius} Å around '{selection}'...")

    truncated = cut_and_cap(structure, cut_selection=selection, cut_radius=radius)
    logger.info(
        f"Truncated region: {len(truncated.atoms)} atoms, {len(truncated.residues)} residues"
    )

    out_path = output_dir / "capped.pdb"
    truncated.save(str(out_path), overwrite=True)
    logger.info(f"Saved: {out_path}")

    return out_path
