"""Shared structure preparation utilities.

Contains:
- Structure cleaning helpers (remove_water, remove_by_residue, remove_chain)
  Ported from PDBProgress/pdbprogress/cleaning/cleaning.py
- H-link capping / cut_and_cap for QM region selection
  Ported from PDBProgress/pdbprogress/capping/cut_and_cap.py
  and supporting utils (parmed_utils.py, parmed_hacks.py, utils.py)

ACE/NME termini capping is excluded — H-link atoms only for now.
"""

from __future__ import annotations

import itertools
import logging
from typing import List, Optional, Union

import numpy as np
import parmed as pmd

logger = logging.getLogger(__name__)

# Equilibrium H bond distances (Å) at peptide backbone cuts
_H_EQ_DIST = {"C": 1.09, "N": 0.99}


# ---------------------------------------------------------------------------
# Internal geometry / parmed helpers (inlined from PDBProgress utils)
# ---------------------------------------------------------------------------

def _as_range(iterable) -> str:
    in_list = list(iterable)
    if len(in_list) > 1:
        return f"{in_list[0]}-{in_list[-1]}"
    return f"{in_list[0]}"


def _get_dashed_ranges(number_list: Union[list, np.ndarray]) -> str:
    """Convert a list of integers to a compact Amber-mask range string.

    E.g. [1, 2, 3, 10] → '1-3,10'
    """
    number_list = sorted(list(number_list))
    return ",".join(
        _as_range(g)
        for _, g in itertools.groupby(
            number_list, key=lambda n, c=itertools.count(): n - next(c)
        )
    )


def _atom_coords(atom: pmd.Atom) -> np.ndarray:
    return np.array([atom.xx, atom.xy, atom.xz])


def _direction_norm(coords1: np.ndarray, coords2: np.ndarray) -> np.ndarray:
    v = coords1 - coords2
    return v / np.linalg.norm(v)


def _add_atom_to_structure(
    structure: pmd.Structure, residue: pmd.Residue, atom: pmd.Atom
) -> None:
    """Append *atom* to *residue* and insert it at the correct position in *structure*.

    Workaround for a parmed edge case where ``structure.add_atom_to_residue``
    raises unexpectedly (see PDBProgress parmed_hacks.py for details).
    """
    last_atom = residue.atoms[-1]
    residue.add_atom(atom)
    if not structure.atoms or last_atom is structure.atoms[-1]:
        structure.atoms.append(atom)
    else:
        structure.atoms.insert(last_atom.idx + 1, atom)


# ---------------------------------------------------------------------------
# Cleaning
# ---------------------------------------------------------------------------

def remove_by_residue(
    structure: pmd.Structure,
    resids: List = [],
    resnames: List = [],
    remove_non_selected: bool = False,
) -> pmd.Structure:
    """Remove residues from *structure* by residue number or name.

    Args:
        resids: Residue numbers to remove (or keep if remove_non_selected=True).
        resnames: Residue names to remove (or keep if remove_non_selected=True).
        remove_non_selected: If True, keep only matched residues (invert logic).
    """
    cleaned = pmd.Structure()
    for residue in structure.residues:
        matched = residue.number in resids or residue.name in resnames
        include = matched if remove_non_selected else not matched
        if include:
            cleaned.residues.append(residue)
            for atom in residue.atoms:
                cleaned.atoms.append(atom)
    for bond in structure.bonds:
        if bond.atom1 in cleaned.atoms and bond.atom2 in cleaned.atoms:
            cleaned.bonds.append(bond)
    return cleaned


def remove_water(structure: pmd.Structure) -> pmd.Structure:
    """Remove all HOH water molecules from *structure*."""
    return remove_by_residue(structure, resnames=["HOH"])


def remove_chain(structure: pmd.Structure, chain_id: str) -> pmd.Structure:
    """Remove an entire chain by chain ID."""
    cleaned = pmd.Structure()
    for residue in structure.residues:
        if residue.chain != chain_id:
            cleaned.residues.append(residue)
            for atom in residue.atoms:
                cleaned.atoms.append(atom)
    for bond in structure.bonds:
        if bond.atom1 in cleaned.atoms and bond.atom2 in cleaned.atoms:
            cleaned.bonds.append(bond)
    return cleaned


def remove_hydrogens(structure: pmd.Structure) -> pmd.Structure:
    """Remove all hydrogen atoms from *structure*."""
    cleaned = pmd.Structure()
    for residue in structure.residues:
        new_residue = pmd.Residue(name=residue.name, number=residue.number)
        cleaned.residues.append(new_residue)
        for atom in residue.atoms:
            if atom.atomic_number != 1:
                cleaned.atoms.append(atom)
                new_residue.atoms.append(atom)
    for bond in structure.bonds:
        if bond.atom1 in cleaned.atoms and bond.atom2 in cleaned.atoms:
            cleaned.bonds.append(bond)
    return cleaned


# ---------------------------------------------------------------------------
# H-link capping / cut_and_cap
# ---------------------------------------------------------------------------

def cut_and_cap(
    struct: pmd.Structure,
    cut_selection: str,
    cut_radius: Optional[float] = None,
) -> pmd.Structure:
    """Extract a region of *struct* and cap backbone cuts with H link atoms.

    Boundary C–N peptide bonds are capped with a hydrogen placed at the
    equilibrium C–H or N–H distance along the cut bond direction.
    ACE/NME capping is not applied here.

    Args:
        struct: Input ParmEd structure (full system).
        cut_selection: Amber mask or 3-letter residue name defining the QM
            region centre. If *cut_radius* is given this is used as the
            centre of a radial expansion.
        cut_radius: If provided, expand the selection to all residues with
            at least one atom within this distance (Å) of *cut_selection*.

    Returns:
        New ParmEd structure containing the selected region with H caps at
        all backbone boundary bonds.
    """
    if cut_radius is not None:
        cut_selection = _residue_radius_selection(struct, cut_selection, cut_radius)
    elif cut_selection[0] not in [":", "!", "@", "("]:
        if len(cut_selection) > 3:
            raise ValueError(
                f"Expected an Amber mask or a 3-letter residue name, got: {cut_selection!r}"
            )
        cut_selection = ":" + cut_selection

    cut_view = struct.view[cut_selection]
    output_struct = struct[cut_selection]

    qm_atoms = [a.idx for a in cut_view.atoms]

    for qm_res_idx, residue in enumerate(cut_view.residues):
        for atom in residue.atoms:
            for bond in atom.bonds:
                if _is_boundary_bond(bond, qm_atoms):
                    h_cap = _get_h_link_atom(atom, bond)
                    _add_atom_to_structure(
                        output_struct, output_struct.residues[qm_res_idx], h_cap
                    )

    return output_struct


def _is_boundary_bond(bond: pmd.Bond, qm_atoms: list) -> bool:
    """Return True if *bond* crosses the QM region boundary."""
    return (bond.atom1.idx in qm_atoms) != (bond.atom2.idx in qm_atoms)


def _residue_radius_selection(
    structure: pmd.Structure,
    residue_selection: Union[str, int],
    cut_off: float,
) -> str:
    """Return an Amber mask of all residues with an atom within *cut_off* Å of *residue_selection*.

    *residue_selection* can be:
    - A single uppercase letter → treated as a chain ID (e.g. 'B')
    - An Amber mask string starting with ':' or '@' → used as-is
    - A residue name or number string → prefixed with ':' automatically
    """
    sel = str(residue_selection).strip()

    # Single uppercase letter with no digits → chain ID, not an Amber mask token.
    # parmed's ':B' matches residue *name* B, not chain B, so we must resolve
    # chain membership via the parmed API and build an atom index mask instead.
    if sel.isalpha() and len(sel) == 1 and sel.isupper():
        atom_indices = [a.idx + 1 for a in structure.atoms if a.residue.chain == sel]
        if not atom_indices:
            raise ValueError(
                f"Chain '{sel}' not found in structure. "
                f"Available chains: {sorted({r.chain for r in structure.residues})}"
            )
        # Build an @<idx1,idx2,...> mask from the chain atom indices
        sel = "@" + ",".join(str(i) for i in atom_indices)
    elif not sel.startswith(":") and not sel.startswith("@") and not sel.startswith("("):
        sel = ":" + sel

    selected = structure.view[f"{sel}<@{cut_off}"]
    residue_range = [r.idx for r in selected.residues]
    if not residue_range:
        raise ValueError(
            f"Residue selection '{residue_selection}' returns no residues in the structure."
        )
    return ":" + _get_dashed_ranges([r + 1 for r in residue_range])


def _get_h_link_atom(atom: pmd.Atom, bond: pmd.Bond) -> pmd.Atom:
    """Create an H link atom along the cut backbone bond.

    Works for backbone C (cap towards N outside QM) and N (cap towards C
    outside QM). Raises ValueError for non-backbone cuts.
    """
    if atom.name not in _H_EQ_DIST:
        raise ValueError(
            f"H-link capping only supports backbone C and N atoms; got atom '{atom.name}'. "
            "Non-peptide-backbone cuts are not yet supported."
        )

    # Identify the atom in QM region (atom) and the one being replaced by H.
    # bond.atom1/atom2 ordering is not guaranteed, so resolve dynamically.
    other = bond.atom2 if bond.atom1 is atom else bond.atom1

    atom_coords = _atom_coords(atom)
    other_coords = _atom_coords(other)

    dist = _H_EQ_DIST[atom.name]
    # Place H at equilibrium distance from atom, along atom→other direction
    link_coords = atom_coords + dist * _direction_norm(other_coords, atom_coords)

    h_cap = pmd.Atom(name=f"HL{atom.name}", atomic_number=1, mass=1.0079)
    h_cap.xx, h_cap.xy, h_cap.xz = link_coords
    return h_cap


# ---------------------------------------------------------------------------
# Term flag parsing (for protein_peptide CLI config)
# ---------------------------------------------------------------------------

def parse_term_flags(flags: list[str]) -> set[tuple[str, int]]:
    """Convert ['A:1', 'B:50'] → {('A', 1), ('B', 50)}.

    Used to pass true_nterm / true_cterm from config to cap_termini.
    """
    result = set()
    for flag in flags:
        chain, resi = flag.split(":")
        result.add((chain, int(resi)))
    return result
