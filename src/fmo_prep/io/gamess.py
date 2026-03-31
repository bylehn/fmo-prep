"""GAMESS .inp file I/O utilities.

Ported from:
- fmo-poc/scripts/map_fragments.py    (parse_inp_file, write_fragmapping)
- fmo-poc/inputs/fn001/debug_fragit.py (validate_fragments)
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# parse_inp_file
# ---------------------------------------------------------------------------

def parse_inp_file(inp_path: str | Path) -> tuple[list[str], list[list[int]]]:
    """Parse FRGNAM and INDAT sections from a GAMESS FMO .inp file.

    Uses the line-by-line state-machine approach from map_fragments.py,
    which is robust against multi-line continuation blocks.

    Args:
        inp_path: Path to the GAMESS .inp file.

    Returns:
        (frag_names, fragments) where:
        - frag_names: list of fragment name strings (e.g. ['MET001', 'GLU002', ...])
        - fragments: list of atom index lists (1-based), one per fragment
    """
    frag_names: list[str] = []
    fragments:  list[list[int]] = []

    parsing_frgnam = False
    parsing_indat  = False
    current_frag:  list[int] = []

    with open(inp_path) as f:
        for line in f:
            # --- FRGNAM block ---
            if "FRGNAM(1)=" in line:
                parsing_frgnam = True
                raw = line.split("=", 1)[1]
                frag_names.extend(
                    x.strip().rstrip(",") for x in re.findall(r"[A-Za-z][A-Za-z0-9]+", raw)
                )
                continue

            if parsing_frgnam:
                if "INDAT(1)=0" in line:
                    parsing_frgnam = False
                    parsing_indat  = True
                    continue
                frag_names.extend(
                    x.strip().rstrip(",") for x in re.findall(r"[A-Za-z][A-Za-z0-9]+", line)
                )
                continue

            # --- INDAT block ---
            if parsing_indat:
                if "$END" in line:
                    parsing_indat = False
                    break

                tokens = line.split()
                if not tokens:
                    continue

                # Stop if we hit a new keyword line (e.g. LAYER(1)=...)
                if tokens[0] and not tokens[0].lstrip("-").isdigit():
                    break

                ends_with_zero = tokens[-1] == "0"
                if ends_with_zero:
                    tokens = tokens[:-1]

                i = 0
                while i < len(tokens):
                    v = int(tokens[i])
                    if v > 0:
                        if i + 1 < len(tokens) and int(tokens[i + 1]) < 0:
                            end_val = abs(int(tokens[i + 1]))
                            current_frag.extend(range(v, end_val + 1))
                            i += 2
                        else:
                            current_frag.append(v)
                            i += 1
                    else:
                        i += 1

                if ends_with_zero and current_frag:
                    fragments.append(current_frag)
                    current_frag = []

    if current_frag:
        fragments.append(current_frag)

    return frag_names, fragments


# ---------------------------------------------------------------------------
# get_pdb_atoms
# ---------------------------------------------------------------------------

def get_pdb_atoms(pdb_path: str | Path) -> list[dict]:
    """Read ATOM/HETATM records from a PDB file (TER records excluded).

    Returns a list of dicts with keys: res_id, atom_name.
    The list index (0-based) corresponds to the 1-based GAMESS INDAT atom index.
    """
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                res_name  = line[17:20].strip()
                chain     = line[21].strip()
                res_num   = line[22:26].strip()
                atoms.append({
                    "res_id": f"{res_name} {chain} {res_num}",
                    "atom_name": atom_name,
                })
    return atoms


# ---------------------------------------------------------------------------
# write_fragmapping
# ---------------------------------------------------------------------------

def write_fragmapping(
    pdb_path: str | Path,
    inp_path: str | Path,
    output_path: str | Path,
) -> Path:
    """Write a human-readable fragment → residue mapping table.

    Format:  Frag | Name       | Residues and Included Atoms

    Args:
        pdb_path: Path to the PDB file (atom order must match the .inp).
        inp_path: Path to the GAMESS .inp file.
        output_path: Where to write the mapping text file.

    Returns:
        Path to the written file.
    """
    pdb_atoms = get_pdb_atoms(pdb_path)
    frag_names, fragments = parse_inp_file(inp_path)
    output_path = Path(output_path)

    header = f"{'Frag':<4} | {'Name':<10} | Residues and Included Atoms\n"
    separator = "-" * 100 + "\n"

    with open(output_path, "w") as f:
        f.write(header)
        f.write(separator)
        for i, (name, atom_indices) in enumerate(zip(frag_names, fragments)):
            residue_map: dict[str, list[str]] = defaultdict(list)
            for idx in atom_indices:
                if idx - 1 < len(pdb_atoms):
                    info = pdb_atoms[idx - 1]
                    residue_map[info["res_id"]].append(info["atom_name"])
            res_strings = [
                f"{res_id} ({', '.join(atom_names)})"
                for res_id, atom_names in residue_map.items()
            ]
            f.write(f"{i+1:<4} | {name:<10} | {'; '.join(res_strings)}\n")

    return output_path


# ---------------------------------------------------------------------------
# parse_frag_map_file
# ---------------------------------------------------------------------------

def parse_frag_map_file(map_file: str | Path) -> dict:
    """Parse a fragmapping.txt into a frag_info dict for the analysis pipeline.

    Expects the pipe-delimited format written by write_fragmapping.

    Returns:
        Dict mapping fragment_id (int) → {'chain': str, 'plot_label': str, 'name': str}
    """
    frag_info: dict[int, dict] = {}
    with open(map_file) as f:
        lines = f.readlines()

    for line in lines[2:]:  # skip header and separator
        if not line.strip():
            continue
        parts = line.split("|")
        if len(parts) < 3:
            continue

        frag_id   = int(parts[0].strip())
        frag_name = parts[1].strip()
        res_data  = parts[2].strip()

        best_label = frag_name
        best_chain = "U"
        max_atoms  = -1

        for res_str in res_data.split(";"):
            res_str = res_str.strip()
            if not res_str:
                continue
            res_parts = res_str.split("(")
            res_info  = res_parts[0].strip().split()
            atom_count = len(res_parts[1].split(",")) if len(res_parts) > 1 else 0
            if atom_count > max_atoms:
                max_atoms = atom_count
                if len(res_info) >= 3:
                    best_label = f"{res_info[0]}{res_info[2]}"
                    best_chain = res_info[1]
                else:
                    best_label = res_parts[0].strip()

        frag_info[frag_id] = {
            "chain": best_chain,
            "plot_label": best_label,
            "name": frag_name,
        }

    return frag_info


# ---------------------------------------------------------------------------
# validate_fragments
# ---------------------------------------------------------------------------

def validate_fragments(inp_path: str | Path, pdb_path: str | Path | None = None) -> None:
    """Print a fragment charge / residue summary and flag anomalies.

    Ported from fmo-poc/inputs/fn001/debug_fragit.py.

    Args:
        inp_path: Path to the GAMESS .inp file.
        pdb_path: Optional PDB file for residue-level labelling. When
            provided, each fragment line shows which residues it contains.
    """
    inp_path = Path(inp_path)
    content  = inp_path.read_text()

    # Parse ICHARG
    icharg_match = re.search(
        r"ICHARG\(1\)\s*=(.*?)(?=\n\s*[A-Z]|\$END)", content, re.DOTALL
    )
    charges: list[int] = []
    if icharg_match:
        charges = [int(x) for x in re.findall(r"-?\d+", icharg_match.group(1))]

    frag_names, frag_atoms = parse_inp_file(inp_path)

    # Optionally map atom serials → residue info from PDB
    serial_to_res: dict[int, tuple[str, str, int]] = {}
    if pdb_path:
        with open(pdb_path) as f:
            for line in f:
                if line[:4] in ("ATOM", "HETA"):
                    try:
                        serial  = int(line[6:11])
                        resname = line[17:20].strip()
                        chain   = line[21]
                        resi    = int(line[22:26])
                        serial_to_res[serial] = (resname, chain, resi)
                    except ValueError:
                        pass

    _NEUTRAL_CAPS   = {"ACE", "NME"}
    _CHARGED_POS    = {"LYS", "ARG", "HIP"}
    _CHARGED_NEG    = {"GLU", "ASP"}

    print(f"  {'Frag':<10} {'Name':<10} {'Charge':>6}  {'Atoms':>5}  Residues")
    print(f"  {'-'*10} {'-'*10} {'-'*6}  {'-'*5}  {'-'*40}")

    for i, atoms in enumerate(frag_atoms):
        name   = frag_names[i] if i < len(frag_names) else f"FRAG{i+1}"
        charge = charges[i]    if i < len(charges)    else "?"

        res_set: set[tuple[str, int, str]] = set()
        for serial in atoms:
            if serial in serial_to_res:
                resname, chain, resi = serial_to_res[serial]
                res_set.add((chain, resi, resname))

        res_str = ", ".join(f"{ch}:{ri} {rn}" for ch, ri, rn in sorted(res_set))

        flag = ""
        if isinstance(charge, int):
            resnames_in_frag = {rn for _, _, rn in res_set}
            if charge != 0 and resnames_in_frag.issubset(_NEUTRAL_CAPS):
                flag = "  *** CAP SHOULD BE NEUTRAL"
            elif charge > 0 and not resnames_in_frag.intersection(_CHARGED_POS | _NEUTRAL_CAPS):
                flag = "  *** UNEXPECTED POSITIVE CHARGE"
            elif charge < 0 and not resnames_in_frag.intersection(_CHARGED_NEG | _NEUTRAL_CAPS):
                flag = "  *** UNEXPECTED NEGATIVE CHARGE"

        print(f"  {i+1:<10} {name:<10} {charge!s:>6}  {len(atoms):>5}  {res_str}{flag}")

    total = sum(c for c in charges[:len(frag_atoms)] if isinstance(c, int))
    print(f"\n  Total system charge: {total}")
