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

import logging
import re
from pathlib import Path

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


def _build_contrl(cfg: FragitConfig) -> str:
    if cfg.calc_mode == "2layer":
        return (
            " $CONTRL NPRINT=-5 ISPHER=1 MAXIT=100\n"
            "         RUNTYP=ENERGY SCFTYP=RHF\n"
            " $END\n"
        )
    else:
        return (
            " $CONTRL NPRINT=-5 ISPHER=1\n"
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


def _parse_indat(text: str) -> list[tuple[int, int]]:
    """Parse INDAT(1) block → list of (start, end) atom index ranges, one per fragment.

    Each line in INDAT has the form ``start  -end  0`` where start and end are
    1-based atom indices (end is stored negative) and 0 terminates the fragment.
    Returns a list whose i-th element is (start, end) for fragment i.
    """
    m = re.search(r"INDAT\(1\)\s*=\s*0\n(.*?)(?:\n\s*\n|\s*\$END)", text, re.DOTALL)
    if not m:
        return []
    fragments = []
    for line in m.group(1).splitlines():
        tokens = line.split()
        if not tokens:
            continue
        # Each line: start -end 0  (may occasionally have multiple range pairs
        # before the trailing 0, but FragIt always writes one pair per line)
        ints = [int(t) for t in tokens]
        positives = [v for v in ints if v > 0]
        negatives = [v for v in ints if v < 0]
        if positives and negatives:
            fragments.append((positives[0], abs(negatives[-1])))
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

    In GAMESS FMO, FragIt cuts each peptide chain at Cα–C backbone bonds.
    Fragment i therefore contains:
      - the C=O of the "donor" residue (named FRGNAM[i])
      - the N, Cα, and side chain of the "inner" residue (named FRGNAM[i+1])

    This split is recorded in $FMOBND as a (BDA, BAA) atom-index pair where
    the BDA belongs to fragment i and the BAA belongs to fragment i+1.

    A non-zero ICHARG on fragment i is a MMFF94 artifact when:
      1. FRGNAM[i+1] is a non-standard residue (not in _STANDARD_RESIDUES), AND
      2. $FMOBND contains a cut between fragment i and fragment i+1.
         (This confirms that the non-standard residue's atoms actually appear
         inside fragment i — as opposed to being a separate ligand/cofactor
         fragment that merely happens to follow fragment i.)

    Charges on fragments NOT meeting both criteria are left untouched,
    preserving legitimate charges from LYS, ARG, GLU, ASP, ligands, etc.
    """
    icharg_m = re.search(r"(ICHARG\(1\)\s*=)(.*?)(\n\s+FRGNAM)", text, re.DOTALL)
    frgnam_m = re.search(r"FRGNAM\(1\)\s*=(.*?)(\n\s+INDAT)", text, re.DOTALL)
    if not icharg_m or not frgnam_m:
        return text

    charges = [int(x) for x in re.findall(r"-?\d+", icharg_m.group(2))]
    names = re.findall(r"[A-Z][A-Z0-9]+\d+", frgnam_m.group(1))

    def resname(frag: str) -> str:
        return re.match(r"([A-Z]+)", frag).group(1)

    # Build atom-index sets per fragment and the set of BDA atom indices.
    frag_ranges = _parse_indat(text)   # list of (start, end)
    fmobnd_pairs = _parse_fmobnd(text) # list of (bda, baa)

    # Map: BDA atom index → BAA atom index (for quick lookup)
    bda_to_baa = {bda: baa for bda, baa in fmobnd_pairs}

    def has_fmobnd_cut(fi: int, fi1: int) -> bool:
        """Return True if FMOBND contains a cut between fragment fi and fi+1."""
        if fi >= len(frag_ranges) or fi1 >= len(frag_ranges):
            return False
        start_i, end_i = frag_ranges[fi]
        start_i1, end_i1 = frag_ranges[fi1]
        for bda, baa in bda_to_baa.items():
            if start_i <= bda <= end_i and start_i1 <= baa <= end_i1:
                return True
        return False

    artifacts = [
        i for i, charge in enumerate(charges)
        if charge != 0
        and i + 1 < len(names)
        and resname(names[i + 1]) not in _STANDARD_RESIDUES
        and has_fmobnd_cut(i, i + 1)
    ]

    if not artifacts:
        return text

    corrected = list(charges)
    for i in artifacts:
        logger.warning(
            "Auto-corrected ICHARG for fragment %s: %+d → 0 "
            "(FMOBND-confirmed cut into non-standard residue %s — likely MMFF94 artifact).",
            names[i], corrected[i], resname(names[i + 1]),
        )
        corrected[i] = 0

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

    # --- Step 1: strip existing header blocks ---
    # $CONTRL may span two lines — handle with DOTALL first
    text = re.sub(r"^ \$CONTRL\b.*?\$END\n", "", text, flags=re.MULTILINE | re.DOTALL)
    # Strip all other single-line blocks (skip index 3 = $CONTRL, already handled)
    for i, pattern in enumerate(_STRIP_PATTERNS):
        if i == 3:
            continue
        text = re.sub(pattern, "", text, flags=re.MULTILINE)

    # --- Step 2: prepend standardised header ---
    header = (
        f" $SYSTEM MWORDS={cfg.mwords} $END\n"
        f" $GDDI NGROUP={cfg.ngroup} $END\n"
        + _build_scf(cfg)
        + _build_contrl(cfg)
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
