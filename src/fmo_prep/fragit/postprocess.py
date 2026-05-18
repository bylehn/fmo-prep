"""GAMESS input file postprocessing.

FragIt writes a complete GAMESS input including all header blocks:
  $SYSTEM, $GDDI, $SCF, $CONTRL, $BASIS, $FMOPRP, $FMO, $FMOBND, $DATA,
  $FMOHYB, $FMOXYZ

This module replaces those header blocks with mode-appropriate versions
and optionally inserts a $PCM block for implicit solvent.

Calculation modes (cfg.calc_mode):

  hf      - HF only, no MP2.
             SCF: CONV=1E-6, DIIS=.F., SOSCF=.T.
             CONTRL: no MAXIT, no SCFTYP

  mp2     - MP2 on entire system (for PIEDA analysis).
             SCF: CONV=1E-7, DIIS=.F., SOSCF=.T.
             CONTRL: no MAXIT, no SCFTYP
             FMOPRP: PRTDST(1)=100.0,0.5,0.6,0.0 IPIEDA=2

  2layer  - MP2 at active site (layer 2), HF elsewhere (layer 1).
             SCF: CONV=1E-6, DIIS=.T., SOSCF=.F.
             CONTRL: MAXIT=100 SCFTYP=RHF
             FMOPRP: MAXIT=100

Implicit solvent (cfg.implicit_solvent=True):
  - Inserts $PCM SOLVNT=WATER IFMO=1 ICOMP=0 $END after $BASIS
  - Changes FMOPRP to IPIEDA=1 instead of mode default
  - SCF: CONV=1E-6, DIIS=.F., SOSCF=.T. (overrides mp2 mode)
"""

from __future__ import annotations

import json
import logging
import re
from collections import Counter
from pathlib import Path
from typing import Optional

from fmo_prep.config import FragitConfig

logger = logging.getLogger(__name__)

# Standard amino acid residue names (3-letter codes) whose MMFF94 charges are
# reliable. Non-standard residue fragments not in this set are candidates for
# automatic charge correction when the total ICHARG is non-zero.
_STANDARD_RESIDUES = {
    # 20 canonical amino acids
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL",
    # protonation / tautomer variants
    "HIE", "HID", "HIP", "CYX", "CYM", "ASH", "GLH", "LYN",
    # common terminal caps
    "ACE", "NME", "NHE", "NMA", "FOR",
    # other frequently seen modifications
    "MSE",  # selenomethionine
    "SEP", "TPO", "PTR",  # phosphorylated residues
}

# Residues with a phosphate group that can legitimately carry ICHARG = -2 for
# interior fragments.  Phosphate dianion at physiological pH contributes -2.
_PHOSPHORYLATED_RESIDUES = {"SEP", "TPO", "PTR"}

# Blocks written by FragIt that we strip and replace
_STRIP_PATTERNS = [
    r"^ \$SYSTEM\b.*?\$END\n",
    r"^ \$GDDI\b.*?\$END\n",
    r"^ \$SCF\b.*?\$END\n",
    r"^ \$CONTRL\b.*?\$END\n",   # may span 2 lines — handled with DOTALL
    r"^ \$BASIS\b.*?\$END\n",
    r"^ \$FMOPRP\b.*?\$END\n",
    r"^ \$PCM\b.*?\$END\n",
]


def _build_scf(cfg: FragitConfig) -> str:
    if cfg.calc_mode == "2layer":
        return " $SCF CONV=1E-6 DIRSCF=.T. NPUNCH=0 DIIS=.T. SOSCF=.F. $END\n"
    elif cfg.calc_mode == "mp2" and not cfg.implicit_solvent:
        return " $SCF CONV=1E-7 DIRSCF=.T. NPUNCH=0 DIIS=.F. SOSCF=.T. $END\n"
    else:  # hf, or mp2/2layer with implicit_solvent
        return " $SCF CONV=1E-6 DIRSCF=.T. NPUNCH=0 DIIS=.F. SOSCF=.T. $END\n"


def _build_contrl(cfg: FragitConfig, afo: bool = False) -> str:
    local = " LOCAL=BOYS" if afo else ""
    if cfg.calc_mode == "2layer":
        return (
            f" $CONTRL NPRINT=-5 ISPHER=1 MAXIT=100{local}\n"
            "         RUNTYP=ENERGY SCFTYP=RHF\n"
            " $END\n"
        )
    else:
        return (
            f" $CONTRL NPRINT=-5 ISPHER=1{local}\n"
            "         RUNTYP=ENERGY\n"
            " $END\n"
        )


def _build_fmoprp(cfg: FragitConfig) -> str:
    if cfg.implicit_solvent:
        return " $FMOPRP NPRINT=9 NGUESS=2 IPIEDA=1 $END\n"
    elif cfg.calc_mode == "2layer":
        return " $FMOPRP NPRINT=9 NGUESS=2 MAXIT=100 $END\n"
    else:  # hf or mp2
        return " $FMOPRP NPRINT=9 NGUESS=2 IPIEDA=2 $END\n"


def _build_basis(cfg: FragitConfig) -> str:
    basis_map = {
        "6-31G*":   "GBASIS=N31 NGAUSS=6 NDFUNC=1",
        "6-31G(d)": "GBASIS=N31 NGAUSS=6 NDFUNC=1",
        "6-31G":    "GBASIS=N31 NGAUSS=6",
        "STO-3G":   "GBASIS=STO NGAUSS=3",
        "3-21G":    "GBASIS=N21 NGAUSS=3",
    }
    gbasis = basis_map.get(cfg.basis, "GBASIS=N31 NGAUSS=6 NDFUNC=1")
    return f" $BASIS {gbasis} $END\n"


def _build_pcm() -> str:
    return " $PCM SOLVNT=WATER IFMO=1 ICOMP=0 $END\n"


def _format_icharg(charges: list[int]) -> str:
    """Render a list of integer charges as a GAMESS ICHARG(1)= block.

    Returns a string beginning with ``ICHARG(1)=`` (no leading whitespace).
    The caller is responsible for preserving whatever indentation precedes
    the original ``ICHARG`` keyword in the file.
    """
    rows = []
    for i in range(0, len(charges), 10):
        chunk = charges[i : i + 10]
        values = ",".join(f"{c:3d}" for c in chunk)
        if i + 10 < len(charges):
            values += ","
        prefix = "ICHARG(1)=" if i == 0 else "          "
        rows.append(prefix + values)
    return "\n".join(rows)


def _parse_indat(text: str) -> list[list[tuple[int, int]]]:
    """Parse INDAT(1) block → list of range lists, one per fragment.

    Each fragment is terminated by a ``0`` token. A positive integer starts a
    range and the immediately following negative integer ends it. A bare positive
    with no following negative (e.g. a cap hydrogen at a backbone cut) is treated
    as a single-atom range [p, p]. Fragments with many atom ranges may span
    multiple lines — 0 is the delimiter, not the line boundary.

    Returns a list whose i-th element is a list of (start, end) tuples covering
    all atom ranges for fragment i.
    """
    m = re.search(r"INDAT\(1\)\s*=\s*0\n(.*?)(?:\n\s*\n|\s*\$END)", text, re.DOTALL)
    if not m:
        return []

    all_tokens = []
    for line in m.group(1).splitlines():
        all_tokens.extend(int(t) for t in line.split())

    fragments = []
    current_ranges: list[tuple[int, int]] = []
    i = 0
    while i < len(all_tokens):
        v = all_tokens[i]
        if v == 0:
            if current_ranges:
                fragments.append(current_ranges)
                current_ranges = []
            i += 1
        elif v > 0:
            if i + 1 < len(all_tokens) and all_tokens[i + 1] < 0:
                current_ranges.append((v, abs(all_tokens[i + 1])))
                i += 2
            else:
                current_ranges.append((v, v))
                i += 1
        else:
            i += 1
    if current_ranges:
        fragments.append(current_ranges)
    return fragments


def _parse_fmobnd(text: str) -> list[tuple[int, int]]:
    """Parse $FMOBND block → list of (bda, baa) cut-bond atom index pairs.

    Each line has the form ``-bda  baa  basis1  basis2``.  The BDA (Bond
    Detached Atom) is stored negative; BAA (Bond Attached Atom) is positive.
    """
    m = re.search(r"^\s*\$FMOBND\b(.*?)^\s*\$END\b", text, re.MULTILINE | re.DOTALL)
    if not m:
        return []
    pairs = []
    for line in m.group(1).splitlines():
        tokens = line.split()
        if len(tokens) >= 2:
            try:
                bda_tok, baa_tok = int(tokens[0]), int(tokens[1])
                if bda_tok < 0 and baa_tok > 0:
                    pairs.append((abs(bda_tok), baa_tok))
            except ValueError:
                pass
    return pairs


def _fix_fragment_charges(text: str) -> str:
    """Auto-correct ICHARG for fragments whose charge is a MMFF94 artifact.

    Two passes are applied:

    Pass 1 — non-standard residue cut artifacts:
      FragIt cuts each peptide chain at Cα–C backbone bonds. Fragment i
      therefore contains the C=O of the "donor" residue (FRGNAM[i]) and
      the N, Cα, and side chain of the "inner" residue (FRGNAM[i+1]).
      A non-zero ICHARG on fragment i is a MMFF94 artifact when:
        1. FRGNAM[i+1] is a non-standard residue (not in _STANDARD_RESIDUES), AND
        2. $FMOBND contains a cut between fragment i and i+1.
      Those fragments are corrected to ICHARG=0.

    Pass 2 — impossible magnitude for standard residues:
      MMFF94 partial-charge rounding can produce |ICHARG| > 1 for standard
      residue fragments near charged residues (e.g. fragment following ARG may
      accumulate +2). No standard amino acid can legitimately contribute
      |ICHARG| > 1 to its fragment. Those values are clamped to sign(charge).
    """
    icharg_m = re.search(r"(ICHARG\(1\)\s*=)(.*?)(\n\s+FRGNAM)", text, re.DOTALL)
    frgnam_m = re.search(r"FRGNAM\(1\)\s*=(.*?)(\n\s+INDAT)", text, re.DOTALL)
    if not icharg_m or not frgnam_m:
        return text

    charges = [int(x) for x in re.findall(r"-?\d+", icharg_m.group(2))]
    names = re.findall(r"[A-Z][A-Z0-9]+\d+", frgnam_m.group(1))

    def resname(frag: str) -> str:
        return re.match(r"([A-Z]+)", frag).group(1)

    frag_ranges = _parse_indat(text)
    fmobnd_pairs = _parse_fmobnd(text)
    bda_to_baa = {bda: baa for bda, baa in fmobnd_pairs}

    def in_fragment(atom: int, fi: int) -> bool:
        return any(s <= atom <= e for s, e in frag_ranges[fi])

    def has_fmobnd_cut(fi: int, fi1: int) -> bool:
        if fi >= len(frag_ranges) or fi1 >= len(frag_ranges):
            return False
        for bda, baa in bda_to_baa.items():
            if in_fragment(bda, fi) and in_fragment(baa, fi1):
                return True
        return False

    corrected = list(charges)
    changed = False

    # Pass 1: correct artifacts caused by non-standard residues at cut sites
    for i, charge in enumerate(charges):
        if (charge != 0
                and i + 1 < len(names)
                and resname(names[i + 1]) not in _STANDARD_RESIDUES
                and has_fmobnd_cut(i, i + 1)):
            logger.warning(
                "Auto-corrected ICHARG for fragment %s: %+d → 0 "
                "(FMOBND-confirmed cut into non-standard residue %s — likely MMFF94 artifact).",
                names[i], charge, resname(names[i + 1]),
            )
            corrected[i] = 0
            changed = True

    # Pass 2: clamp |charge| > 1 for interior standard residue fragments.
    # Exemptions: chain termini and phosphorylated residues (legitimate -2).
    for i, charge in enumerate(corrected):
        if abs(charge) <= 1 or resname(names[i]) not in _STANDARD_RESIDUES:
            continue
        if resname(names[i]) in _PHOSPHORYLATED_RESIDUES:
            continue
        if i >= len(frag_ranges):
            continue
        has_outgoing = any(in_fragment(bda, i) for bda in bda_to_baa)
        has_incoming = any(in_fragment(baa, i) for baa in bda_to_baa.values())
        if not has_outgoing or not has_incoming:
            continue  # chain terminus — charge may legitimately exceed ±1
        clamped = 1 if charge > 0 else -1
        logger.warning(
            "Auto-corrected ICHARG for fragment %s: %+d → %+d "
            "(|charge| > 1 impossible for interior standard residue — MMFF94 partial-charge rounding artifact).",
            names[i], charge, clamped,
        )
        corrected[i] = clamped
        changed = True

    if not changed:
        return text

    new_icharg = _format_icharg(corrected)
    text = text[: icharg_m.start()] + new_icharg + text[icharg_m.end() - len(icharg_m.group(3)):]
    return text


def patch_inp(inp_path: Path, cfg: FragitConfig, output_path: Path | None = None) -> Path:
    """Patch a FragIt-generated GAMESS .inp file with mode-appropriate header blocks.

    Steps applied:
    1. Strip FragIt's $SYSTEM, $GDDI, $SCF, $CONTRL, $BASIS, $FMOPRP (and $PCM if present).
    2. Prepend our versions of those blocks, chosen based on cfg.calc_mode.
    3. Set NLAYER in $FMO: 2layer→2, hf/mp2→1 (override whatever FragIt wrote).
       Set MPLEVL(1) in $FMO: hf→0, mp2→2, 2layer→0,2.
    4. Replace RESDIM/RCORSD in $FMO if non-default values are configured.
    5. Auto-correct ICHARG if total fragment charge ≠ 0 (see _fix_fragment_charges).

    Args:
        inp_path: Path to the FragIt-generated .inp file.
        cfg: FragitConfig supplying all GAMESS header parameters.
        output_path: Destination path. Defaults to overwriting inp_path.

    Returns:
        Path to the patched file.

    Raises:
        ValueError: If the file contains no $FMO block.
    """
    inp_path = Path(inp_path)
    output_path = Path(output_path) if output_path else inp_path

    text = inp_path.read_text()

    if "$FMO" not in text:
        raise ValueError(f"No $FMO block found in {inp_path} — is this a valid FragIt .inp?")

    # Detect AFO mode before stripping — fragit writes RAFO(1)= in $FMO when dohop=False.
    afo_mode = bool(re.search(r"^\s*RAFO\(1\)\s*=", text, flags=re.MULTILINE))

    # --- Step 1: strip existing header blocks ---
    # $CONTRL may span two lines — handle with DOTALL first
    text = re.sub(r"^ \$CONTRL\b.*?\$END\n", "", text, flags=re.MULTILINE | re.DOTALL)
    # Strip all other single-line blocks (skip index 3 = $CONTRL, already handled)
    for i, pattern in enumerate(_STRIP_PATTERNS):
        if i == 3:
            continue
        text = re.sub(pattern, "", text, flags=re.MULTILINE)

    # Remove RAFO only for HOP mode (dohop=True); AFO mode needs RAFO(1)= in $FMO.
    if not afo_mode:
        text = re.sub(r"^\s*RAFO\(1\)\s*=.*\n", "", text, flags=re.MULTILINE)

    # --- Step 2: prepend standardised header ---
    header = (
        f" $SYSTEM MWORDS={cfg.mwords} $END\n"
        f" $GDDI NGROUP={cfg.ngroup} $END\n"
        + _build_scf(cfg)
        + _build_contrl(cfg, afo=afo_mode)
        + _build_basis(cfg)
        + (_build_pcm() if cfg.implicit_solvent else "")
        + _build_fmoprp(cfg)
    )
    text = header + text.lstrip("\n")

    # --- Step 3: set NLAYER and MPLEVL in $FMO block ---
    # NLAYER: 2layer mode needs 2 (active-site layer + environment), all others need 1.
    # FragIt may write NLAYER=2 for any run when mp2level is set; we always override it.
    nlayer_val = "2" if cfg.calc_mode == "2layer" else "1"
    if re.search(r"^\s*NLAYER\s*=", text, flags=re.MULTILINE):
        text = re.sub(
            r"(^\s*NLAYER\s*=\s*)\S+",
            lambda m: m.group(1) + nlayer_val,
            text, flags=re.MULTILINE,
        )
    else:
        # Insert before ICHARG
        nlayer_anchor = re.search(r"^\s*ICHARG\(1\)\s*=", text, flags=re.MULTILINE)
        if nlayer_anchor:
            indent = "      "
            text = text[:nlayer_anchor.start()] + f"{indent}NLAYER={nlayer_val}\n" + text[nlayer_anchor.start():]

    mplevl_map = {"hf": "0", "mp2": "2", "2layer": "0,2"}
    mplevl_val = mplevl_map[cfg.calc_mode]
    if re.search(r"^\s*MPLEVL\(1\)\s*=", text, flags=re.MULTILINE):
        text = re.sub(
            r"(^\s*MPLEVL\(1\)\s*=\s*)\S+",
            lambda m: m.group(1) + mplevl_val,
            text, flags=re.MULTILINE,
        )
    else:
        # Insert after NLAYER line
        nlayer_line = re.search(r"^\s*NLAYER\s*=\s*\S+\n", text, flags=re.MULTILINE)
        if nlayer_line:
            pos = nlayer_line.end()
        else:
            mplevl_anchor = re.search(r"^\s*ICHARG\(1\)\s*=", text, flags=re.MULTILINE)
            pos = mplevl_anchor.start() if mplevl_anchor else None
        if pos is not None:
            indent = "      "
            text = text[:pos] + f"{indent}MPLEVL(1)={mplevl_val}\n" + text[pos:]

    # --- Steps 4–5: patch RESDIM / RCORSD if non-default ---
    if cfg.resdim != 2.0:
        text = re.sub(
            r"(^\s*RESDIM\s*=\s*)\S+",
            lambda m: m.group(1) + str(cfg.resdim),
            text, flags=re.MULTILINE,
        )
    if cfg.rcorsd != 2.0:
        text = re.sub(
            r"(^\s*RCORSD\s*=\s*)\S+",
            lambda m: m.group(1) + str(cfg.rcorsd),
            text, flags=re.MULTILINE,
        )

    # --- Step 6: auto-correct ICHARG if total charge ≠ 0 ---
    text = _fix_fragment_charges(text)

    output_path.write_text(text)
    return output_path


def _parse_pdb_atoms(pdb_path: Path) -> list[dict]:
    """Parse every ATOM and HETATM record from a PDB file.

    Both record types are included so that ligands (HETATM) are present
    alongside standard residues (ATOM) — FragIt writes all of them into
    $FMOXYZ, so the PDB atom list must contain them for coordinate matching.
    """
    atoms = []
    for line in Path(pdb_path).read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            atoms.append(
                {
                    "index": int(line[6:11].strip()),
                    "name": line[12:16].strip(),
                    "resname": line[17:20].strip(),
                    "chain": line[21].strip(),
                    "resid": int(line[22:26].strip()),
                    "x": float(line[30:38].strip()),
                    "y": float(line[38:46].strip()),
                    "z": float(line[46:54].strip()),
                }
            )
    return atoms


def _parse_fmoxyz_atoms(text: str) -> list[tuple[float, float, float]]:
    """Extract atom coordinates from the $FMOXYZ block.

    Each line: LABEL  NUCLEAR_CHARGE  X  Y  Z (Å).
    Returns a list parallel to INDAT order (element 0 → INDAT atom 1).
    """
    m = re.search(r"^\s*\$FMOXYZ\b(.*?)^\s*\$END\b", text, re.MULTILINE | re.DOTALL)
    if not m:
        return []
    coords = []
    for line in m.group(1).splitlines():
        tokens = line.split()
        if len(tokens) >= 5:
            try:
                coords.append((float(tokens[2]), float(tokens[3]), float(tokens[4])))
            except ValueError:
                pass
    return coords


def _match_fmoxyz_to_pdb(
    fmoxyz_coords: list[tuple[float, float, float]],
    pdb_atoms: list[dict],
) -> list[dict]:
    """Map each FMOXYZ atom to its PDB counterpart by coordinate proximity (1e-3 Å).

    Returns a list parallel to fmoxyz_coords: element i → PDB atom for INDAT atom i+1.

    Raises:
        ValueError: If any FMOXYZ coordinate cannot be matched to a PDB atom.
    """
    tol = 1e-3
    matched = []
    for i, (x, y, z) in enumerate(fmoxyz_coords):
        pdb_atom = next(
            (
                p for p in pdb_atoms
                if abs(p["x"] - x) < tol and abs(p["y"] - y) < tol and abs(p["z"] - z) < tol
            ),
            None,
        )
        if pdb_atom is None:
            raise ValueError(
                f"No PDB atom matches INDAT atom {i + 1} at ({x:.3f}, {y:.3f}, {z:.3f}). "
                f"Ensure the PDB passed to build_fragment_residue_map is the same file "
                f"that was given to FragIt."
            )
        matched.append(pdb_atom)
    return matched


def _map_fragment2residues(
    indat_ranges: list[list[tuple[int, int]]],
    atom_residues: list[dict],
) -> list[dict]:
    """Assign each FMO fragment to its dominant PDB residue using INDAT atom ranges.

    Returns a list of fragment dicts:
    {
        "fragment_index": int,       # 1-based
        "chain": str,
        "majority_residue": str,
        "majority_resnum": int,
        "all_residues": [{"resname": str, "resnum": int, "chain": str}, ...]
    }
    """
    residue_atom_totals = Counter(
        (a["chain"], a["resname"], a["resid"]) for a in atom_residues
    )

    result = []
    for i, ranges in enumerate(indat_ranges):
        fragment_atoms = []
        for start, end in ranges:
            fragment_atoms += atom_residues[start - 1 : end]

        fragment_counts = Counter(
            (a["chain"], a["resname"], a["resid"]) for a in fragment_atoms
        )

        maj_chain, maj_resname, maj_resnum = fragment_counts.most_common(1)[0][0]

        all_residues = []
        for (chain, resname, resnum), count in fragment_counts.items():
            if count / residue_atom_totals[(chain, resname, resnum)] > 0.5:
                all_residues.append({"resname": resname, "resnum": resnum, "chain": chain})

        result.append(
            {
                "fragment_index": i + 1,
                "chain": maj_chain,
                "majority_residue": maj_resname,
                "majority_resnum": maj_resnum,
                "all_residues": all_residues,
            }
        )
    return result


def build_fragment_residue_map(inp_path: Path, pdb_path: Path) -> list[dict]:
    """Build a mapping from FMO fragment indices to PDB residues.

    Matches $FMOXYZ atom coordinates back to PDB atoms (tolerance 1e-3 Å) to
    produce a reliable fragment→residue mapping even in multi-chain systems.

    Writes ``fragment_map.json`` alongside the .inp file and returns the list.

    Args:
        inp_path: Path to the FragIt-generated (and patched) GAMESS .inp file.
        pdb_path: Path to the PDB file given to FragIt — coordinates must match.

    Returns:
        List of fragment dicts as described in _map_fragment2residues.
    """
    text = inp_path.read_text()
    pdb_atoms = _parse_pdb_atoms(pdb_path)
    fmoxyz_coords = _parse_fmoxyz_atoms(text)
    atom_residues = _match_fmoxyz_to_pdb(fmoxyz_coords, pdb_atoms)
    indat_ranges = _parse_indat(text)
    fragment_map = _map_fragment2residues(indat_ranges, atom_residues)

    output_path = inp_path.parent / "fragment_map.json"
    output_path.write_text(json.dumps(fragment_map, indent=2))
    logger.info("Wrote fragment residue map to %s (%d fragments)", output_path, len(fragment_map))

    return fragment_map


def find_fragment_by_chain_resname(
    fragment_map: list[dict],
    chain: Optional[str],
    resname: str,
) -> int:
    """Return the 1-based fragment index whose majority residue matches chain+resname.

    Args:
        fragment_map: Output of build_fragment_residue_map.
        chain: Chain ID to match, or None to ignore chain.
        resname: Residue name to match.

    Raises:
        ValueError: If no matching fragment is found.
    """
    for entry in fragment_map:
        if chain is None:
            if entry["majority_residue"] == resname:
                return entry["fragment_index"]
        else:
            if entry["chain"] == chain and entry["majority_residue"] == resname:
                return entry["fragment_index"]

    if chain is None:
        raise ValueError(f"No fragment with residue {resname}")
    raise ValueError(f"No fragment in chain {chain} with residue {resname}")
