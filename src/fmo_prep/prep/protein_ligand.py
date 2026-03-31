"""Protein-ligand system preparation.

Workflow
--------
1. Load complex PDB with parmed
2. Identify non-standard residues (ligands) via identify_ligands
3. Prepare structure for Amber (strip Hs, fix His naming, fix halogen atom names)
4. Parameterise each ligand with Antechamber (GAFF2) + parmchk2
5. Align atom names in the complex to the parameterised mol2 files
6. Build Amber topology with TLeap (ff14SB + GAFF2 + TIP3P)
7. Run sander energy minimisation
8. Save the minimised complex PDB ready for FragIt

Source: bison/services/simulation-setup-service/app/task.py
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import parmed as pmd

from fmo_prep.config import RunConfig

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Standard residue sets (used by identify_ligands)
# ---------------------------------------------------------------------------

STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HIE", "HID", "HIP",
}
STANDARD_NUCLEIC_ACIDS = {"DA", "DT", "DG", "DC", "A", "U", "G", "C"}
WATER_RESIDUES = {"HOH", "WAT", "H2O", "TIP3", "SPC", "SOL"}
COMMON_IONS = {
    "NA", "CL", "MG", "CA", "ZN", "FE", "K", "MN", "CU", "CO", "NI",
    "SO4", "PO4", "NO3", "NH4",
}

_STANDARD_RESIDUES = (
    STANDARD_AMINO_ACIDS | STANDARD_NUCLEIC_ACIDS | WATER_RESIDUES | COMMON_IONS
)


# ---------------------------------------------------------------------------
# Structure helpers (also used by other prep modules)
# ---------------------------------------------------------------------------

def identify_ligands(structure: pmd.Structure) -> List[pmd.Residue]:
    """Return all non-standard residues (ligands) in *structure*."""
    return [r for r in structure.residues if r.name not in _STANDARD_RESIDUES]


def prepare_structure_for_amber(
    structure: pmd.Structure, strip_hydrogens: bool = True
) -> pmd.Structure:
    """Strip hydrogens and fix naming conventions for Amber/LEaP compatibility.

    - Strips all H atoms so LEaP/Antechamber regenerate them consistently
    - Fixes histidine protonation names (HIS → HID/HIE/HIP based on H atoms)
    - Fixes halogen atom name conventions (CL→Cl, BR→Br)
    """
    prepared = structure

    if strip_hydrogens:
        logger.info("Stripping hydrogens for LEaP/Antechamber regeneration...")
        prepared = structure["!@H="]
        logger.info(f"After strip: {len(prepared.atoms)} atoms")

    for residue in prepared.residues:
        if residue.name == "HIS":
            h_atoms = [a.name for a in residue.atoms if a.atomic_number == 1]
            if "HD1" in h_atoms and "HE2" in h_atoms:
                residue.name = "HIP"
            elif "HD1" in h_atoms:
                residue.name = "HID"
            elif "HE2" in h_atoms:
                residue.name = "HIE"

    for atom in prepared.atoms:
        if atom.name.startswith("CL"):
            atom.name = atom.name.replace("CL", "Cl", 1)
        elif atom.name.startswith("BR"):
            atom.name = atom.name.replace("BR", "Br", 1)

    return prepared


# ---------------------------------------------------------------------------
# Antechamber / parmchk2 wrapper
# ---------------------------------------------------------------------------

class AntechamberWrapper:
    """Run Antechamber + parmchk2 to parameterise a small-molecule ligand."""

    def __init__(self, param_type: str = "gaff2", charge_model: str = "bcc"):
        if param_type not in ("gaff", "gaff2"):
            raise ValueError("param_type must be 'gaff' or 'gaff2'")
        self.param_type = param_type
        self.charge_model = charge_model
        self.atom_type = param_type  # gaff or gaff2

    def run(
        self,
        input_structure: str,
        net_charge: Optional[int] = None,
        residue_name: str = "UNL",
        target_dir: str = ".",
    ) -> Tuple[str, List[str], int]:
        """Parameterise *input_structure* (PDB) and write mol2 + frcmod.

        Returns:
            (mol2_path, [frcmod_path], net_charge)
        """
        if len(residue_name) > 3:
            raise ValueError("residue_name must be ≤ 3 characters")

        target_dir = Path(target_dir)
        target_dir.mkdir(parents=True, exist_ok=True)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            mol2_file = tmpdir_path / f"{residue_name}.mol2"

            cmd = [
                "antechamber",
                "-i", str(input_structure),
                "-fi", "pdb",
                "-o", str(mol2_file),
                "-fo", "mol2",
                "-at", self.atom_type,
                "-c", self.charge_model,
                "-rn", residue_name,
            ]
            if net_charge is not None:
                cmd.extend(["-nc", str(net_charge)])

            logger.info(f"Running Antechamber for {residue_name}")
            result = subprocess.run(
                cmd, cwd=tmpdir,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            )
            if result.returncode != 0:
                raise RuntimeError(f"Antechamber failed:\n{result.stderr}")

            frcmod_file = tmpdir_path / f"{residue_name}.frcmod"
            parmchk_cmd = [
                "parmchk2",
                "-i", str(mol2_file),
                "-f", "mol2",
                "-o", str(frcmod_file),
                "-s", self.param_type,
            ]
            logger.info(f"Running parmchk2 for {residue_name}")
            parmchk_result = subprocess.run(
                parmchk_cmd, cwd=tmpdir,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            )
            if parmchk_result.returncode != 0:
                raise RuntimeError(f"parmchk2 failed:\n{parmchk_result.stderr}")

            final_mol2 = target_dir / f"{residue_name}.mol2"
            final_frcmod = target_dir / f"{residue_name}.frcmod"
            shutil.copy2(mol2_file, final_mol2)
            shutil.copy2(frcmod_file, final_frcmod)

            if net_charge is None:
                net_charge = self._extract_charge_from_mol2(final_mol2)

            return str(final_mol2), [str(final_frcmod)], net_charge

    @staticmethod
    def _extract_charge_from_mol2(mol2_file: Path) -> int:
        try:
            with open(mol2_file) as f:
                for line in f:
                    if line.startswith("@<TRIPOS>MOLECULE"):
                        next(f)  # molecule name
                        parts = next(f).strip().split()
                        if len(parts) >= 2:
                            return int(float(parts[1]))
        except Exception as e:
            logger.warning(f"Could not extract charge from mol2: {e}")
        return 0


# ---------------------------------------------------------------------------
# TLeap wrapper
# ---------------------------------------------------------------------------

class TLeapWrapper:
    """Build an Amber topology using tleap."""

    def __init__(self, work_dir: str = "."):
        self.script_lines: List[str] = []
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)

    def load_parameters(self, *param_files: str) -> None:
        for param_file in param_files:
            if os.path.isfile(param_file):
                abs_path = os.path.abspath(param_file)
                if param_file.endswith(".frcmod"):
                    self.script_lines.append(f"loadAmberParams {abs_path}")
                elif param_file.endswith((".lib", ".off")):
                    self.script_lines.append(f"loadOff {abs_path}")
                else:
                    self.script_lines.append(f"source {abs_path}")
            else:
                self.script_lines.append(f"source {param_file}")

    def load_unit(self, unit_name: str, file_path: str) -> None:
        abs_path = os.path.abspath(file_path)
        if file_path.endswith(".mol2"):
            self.script_lines.append(f"{unit_name} = loadMol2 {abs_path}")
        elif file_path.endswith(".pdb"):
            self.script_lines.append(f"{unit_name} = loadPdb {abs_path}")
        else:
            raise ValueError(f"Unsupported file format for tleap: {file_path}")

    def save_unit(self, unit_name: str, output_path: str) -> None:
        abs_path = os.path.abspath(output_path)
        if output_path.endswith((".prmtop", ".parm7")):
            coord_ext = ".inpcrd" if output_path.endswith(".prmtop") else ".rst7"
            coord_path = os.path.abspath(str(Path(output_path).with_suffix(coord_ext)))
            self.script_lines.append(f"saveAmberParm {unit_name} {abs_path} {coord_path}")
        elif output_path.endswith(".pdb"):
            self.script_lines.append(f"savePdb {unit_name} {abs_path}")

    def set_pb_radii(self, param_set: str = "mbondi3") -> None:
        if param_set not in ("bondi", "mbondi", "mbondi2", "mbondi3"):
            raise ValueError(f"Invalid PB radii set: {param_set}")
        self.script_lines.append(f"set default PBRadii {param_set}")

    def run(self, script_path: Optional[str] = None) -> List[str]:
        """Execute the accumulated tleap script; return any WARNING lines."""
        script_content = "\n".join(self.script_lines) + "\nquit\n"

        if script_path:
            script_file = Path(script_path)
            script_file.parent.mkdir(parents=True, exist_ok=True)
            script_file.write_text(script_content)
            cmd = ["tleap", "-f", str(script_file)]
            input_data = None
        else:
            cmd = ["tleap"]
            input_data = script_content

        result = subprocess.run(
            cmd,
            input=input_data,
            cwd=str(self.work_dir),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"TLeap failed:\n{result.stdout}")

        return [line.strip() for line in result.stdout.splitlines() if "WARNING:" in line]


# ---------------------------------------------------------------------------
# Minimisation
# ---------------------------------------------------------------------------

def _create_minimisation_input(output_path: Path, level: str = "low") -> Path:
    """Write an Amber sander minimisation input file."""
    if level == "high":
        maxcyc, drms, ncyc = 5000, 0.0001, 500
    else:
        maxcyc, drms, ncyc = 1000, 0.001, 100

    min_in = output_path / "min.in"
    min_in.write_text(
        f"Minimize molecule\n"
        f"&cntrl\n"
        f"    imin=1,\n"
        f"    maxcyc={maxcyc},\n"
        f"    drms={drms},\n"
        f"    ncyc={ncyc},\n"
        f"    ntpr=50,\n"
        f"    cut=999.0,\n"
        f"    igb=0,\n"
        f"    ntb=0,\n"
        f"/\n"
        f"&end\n"
    )
    return min_in


def _run_minimisation(
    topology: Path, coords: Path, min_input: Path, output_dir: Path, prefix: str
) -> Path:
    """Run sander minimisation; return path to output coordinate file."""
    out_coord = output_dir / f"{prefix}_min.rst7"
    mdout = output_dir / f"{prefix}_min.out"

    cmd = [
        "sander",
        "-O",
        "-i", str(min_input),
        "-p", str(topology),
        "-c", str(coords),
        "-o", str(mdout),
        "-r", str(out_coord),
    ]
    logger.info("Running sander minimisation...")
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"sander minimisation failed:\n{result.stdout}")

    logger.info(f"Minimisation complete → {out_coord}")
    return out_coord


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run(input_files: dict[str, str], output_dir: Path, cfg: RunConfig) -> Path:
    """Prepare a protein-ligand complex for FMO fragmentation.

    Args:
        input_files: Must contain 'complex' (PDB). 'ligand' (SDF/PDB) is
            optional — if absent the ligand is identified in the complex PDB.
        output_dir: Directory for all intermediate and final outputs.
        cfg: Full run configuration (uses cfg.prep.amber_level).

    Returns:
        Path to the minimised complex PDB ready for FragIt.
    """
    complex_pdb = input_files.get("complex")
    if not complex_pdb or not os.path.exists(complex_pdb):
        raise FileNotFoundError(f"Complex PDB not found: {complex_pdb!r}")

    amber_level = cfg.prep.amber_level
    prefix = cfg.fragit.central_fragment_resname.lower()

    logger.info(f"Loading structure: {complex_pdb}")
    structure = pmd.load_file(complex_pdb)

    # Identify ligands
    ligands = identify_ligands(structure)
    logger.info(f"Found {len(ligands)} ligand(s): {[r.name for r in ligands]}")

    # Prepare structure (strip Hs, fix naming)
    prepared = prepare_structure_for_amber(structure)

    ligand_mol2s: List[str] = []
    ligand_frcmods: List[str] = []
    ligand_resnames: List[str] = []

    if ligands:
        ligand_dir = output_dir / "ligands"
        ligand_dir.mkdir(exist_ok=True)
        antechamber = AntechamberWrapper(param_type="gaff2")

        for i, lig_residue in enumerate(ligands):
            resname = lig_residue.name[:3].strip() or f"LIG{i+1:02d}"
            logger.info(f"Parameterising ligand {i+1}: {resname}")

            # Extract ligand atoms to a temporary PDB
            lig_pdb = ligand_dir / f"{resname}.pdb"
            start = lig_residue.atoms[0].idx
            end = lig_residue.atoms[-1].idx
            structure[start: end + 1].save(str(lig_pdb), overwrite=True)

            try:
                mol2_file, frcmods, charge = antechamber.run(
                    input_structure=str(lig_pdb),
                    residue_name=resname,
                    target_dir=str(ligand_dir),
                )
            except RuntimeError:
                logger.error(f"Parameterisation failed for {resname}")
                raise

            ligand_mol2s.append(mol2_file)
            ligand_frcmods.extend(frcmods)
            ligand_resnames.append(resname)
            logger.info(f"Parameterised {resname} (charge: {charge})")

            # Align atom names in complex to mol2 so LEaP can map parameters
            mol2_struct = pmd.load_file(mol2_file)
            target_res = next(
                (r for r in prepared.residues
                 if r.name == lig_residue.name and r.idx == lig_residue.idx),
                None,
            )
            if target_res:
                p_heavy = [a for a in target_res.atoms if a.atomic_number != 1]
                m_heavy = [a for a in mol2_struct.atoms if a.atomic_number != 1]
                if len(p_heavy) != len(m_heavy):
                    logger.warning(
                        f"Heavy atom count mismatch for {resname}: "
                        f"PDB={len(p_heavy)}, mol2={len(m_heavy)}"
                    )
                for p_atom, m_atom in zip(p_heavy, m_heavy):
                    p_atom.name = m_atom.name
                    p_atom.xx, p_atom.xy, p_atom.xz = m_atom.xx, m_atom.xy, m_atom.xz
                logger.info(f"Aligned {len(p_heavy)} heavy atoms for {resname}")

    # Save prepared PDB
    prepared_pdb = output_dir / f"{prefix}_prepared.pdb"
    prepared.save(str(prepared_pdb), overwrite=True)

    # Build Amber topology with TLeap
    logger.info("Building Amber topology with TLeap...")
    tleap = TLeapWrapper(work_dir=str(output_dir))
    tleap.load_parameters("leaprc.protein.ff14SB")
    if ligands:
        tleap.load_parameters("leaprc.gaff2")
    tleap.load_parameters("leaprc.water.tip3p")
    for frcmod in ligand_frcmods:
        tleap.load_parameters(frcmod)
    for mol2_file, resname in zip(ligand_mol2s, ligand_resnames):
        tleap.load_unit(resname, mol2_file)
    tleap.set_pb_radii("mbondi3")
    tleap.load_unit("system", str(prepared_pdb))

    topology_file = output_dir / f"{prefix}.prmtop"
    tleap.save_unit("system", str(topology_file))

    leap_script = output_dir / "leap.in"
    warnings = tleap.run(str(leap_script))
    if warnings:
        logger.warning(f"TLeap warnings: {warnings}")

    coord_file = output_dir / f"{prefix}.inpcrd"
    if not coord_file.exists():
        raise RuntimeError(f"TLeap did not create coordinate file: {coord_file}")

    # Energy minimisation
    logger.info("Running energy minimisation...")
    min_input = _create_minimisation_input(output_dir, level=amber_level)
    min_coords = _run_minimisation(topology_file, coord_file, min_input, output_dir, prefix)

    # Convert minimised coordinates back to PDB
    logger.info("Converting minimised coordinates to PDB...")
    min_struct = pmd.load_file(str(topology_file), xyz=str(min_coords))
    out_pdb = output_dir / f"{prefix}_minimised.pdb"
    min_struct.save(str(out_pdb), overwrite=True)
    logger.info(f"Minimised complex: {out_pdb}")

    return out_pdb
