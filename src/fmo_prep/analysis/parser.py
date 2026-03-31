"""GAMESS FMO log file parser.

Unified port of:
- fmo-poc/scripts/gamout.py      (PIEDA mode: ES, EX, CT, DI, SOL components)
- fmo-poc/scripts/gamout_nopieda.py (non-PIEDA mode: total pair energies only)

Main entry point::

    df = parse_gamout("fmo_run.log", pieda=True)
    df.to_csv("interactions.csv", index=False)
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


# Fragment name table: "CONV\n ===...\n<rows>\n\n Close fragment pairs"
_FRAGMENT_RE = re.compile(
    r"CONV\n ={70,90}\n(.*?)\n\n Close fragment pairs", re.DOTALL
)

# PIEDA table: " ---...(105 dashes)...\n<rows>\n\n Total energy"
_PIEDA_RE = re.compile(r" -{105}\n(.*?)\n\n Total energy", re.DOTALL)

# Non-PIEDA table: " ---...(120 dashes)...\n<rows>\n\n Total energy"
_NOPIEDA_RE = re.compile(r" -{120}\n(.*?)\n\n Total energy", re.DOTALL)


def _parse_fragments(gamout_str: str, convert: bool = False) -> dict[int, str]:
    """Extract fragment id → name mapping from a GAMESS FMO log string."""
    frgs: dict[int, str] = {}
    for m in _FRAGMENT_RE.finditer(gamout_str):
        for line in m.group(1).split("\n"):
            if not line.strip():
                continue
            els = line.split()
            name = els[1]
            if convert:
                name = _convert_frgname(name)
            frgs[int(els[0])] = name
    return frgs


def _convert_frgname(frgname: str) -> str:
    """Convert FragIt FRGNAM to 'chain:seqresname' format."""
    asym_id = frgname[-1:]
    seq_id  = frgname[-5:-1]
    comp_id = frgname[:-5]
    return f"{asym_id}:{seq_id}{comp_id}"


def parse_gamout(log_path: str | Path, pieda: bool = True) -> pd.DataFrame:
    """Parse a GAMESS FMO output file and return a tidy DataFrame.

    Args:
        log_path: Path to the GAMESS .log file.
        pieda: If True, parse PIEDA decomposition (ES, EX, CT, DI, SOL).
               If False, parse total pair energies only (no PIEDA run).

    Returns:
        DataFrame with columns:
        - pieda=True:  I, IFRG, J, JFRG, R, Q, COMPONENT, ENERGY, TOTAL
        - pieda=False: I, IFRG, J, JFRG, R, Q, ENERGY
    """
    gamout_str = Path(log_path).read_text()
    frgs = _parse_fragments(gamout_str)

    rows: list[dict] = []

    if pieda:
        for m in _PIEDA_RE.finditer(gamout_str):
            for line in m.group(1).split("\n"):
                if not line.strip():
                    continue
                # Fixed-width columns from gamout.py
                f = {
                    "I":       line[:5],
                    "J":       line[5:10],
                    "DL":      line[10:13],
                    "Z":       line[13:16],
                    "R":       line[16:23],
                    "QIJ":     line[23:31],
                    "EIJ":     line[31:41],
                    "dDIJVIJ": line[41:50],
                    "total":   line[50:60],
                    "Ees":     line[60:70],
                    "Eex":     line[70:79],
                    "Ectmix":  line[79:88],
                    "Edisp":   line[88:97],
                    "Gsol":    line[97:106],
                }
                I, J = int(f["I"]), int(f["J"])
                qij   = float(f["QIJ"])
                total = f["total"].strip()
                R     = float(f["R"])
                for tag, raw in [
                    ("ES",  f["Ees"]),
                    ("EX",  f["Eex"]),
                    ("CT",  f["Ectmix"]),
                    ("DI",  f["Edisp"]),
                    ("SOL", f["Gsol"]),
                ]:
                    rows.append({
                        "I":         I,
                        "IFRG":      frgs.get(I, "UNK"),
                        "J":         J,
                        "JFRG":      frgs.get(J, "UNK"),
                        "R":         R,
                        "Q":         qij,
                        "COMPONENT": tag,
                        "ENERGY":    raw.strip(),
                        "TOTAL":     total,
                    })
        df = pd.DataFrame(rows, columns=["I", "IFRG", "J", "JFRG", "R", "Q", "COMPONENT", "ENERGY", "TOTAL"])
        df["ENERGY"] = pd.to_numeric(df["ENERGY"], errors="coerce")
        df["TOTAL"]  = pd.to_numeric(df["TOTAL"],  errors="coerce")
    else:
        for m in _NOPIEDA_RE.finditer(gamout_str):
            for line in m.group(1).split("\n"):
                if not line.strip():
                    continue
                f = {
                    "I":              line[:5],
                    "J":              line[5:10],
                    "DL":             line[10:13],
                    "Z":              line[13:16],
                    "R":              line[16:23],
                    "QIJ":            line[23:31],
                    "Ecorr":          line[31:48],
                    "Euncorr":        line[48:65],
                    "EIJ-EI-EJ,corr": line[65:78],
                    "EIJ-EI-EJ,unc":  line[65:91],
                    "dDIJVIJ,unc":    line[91:103],
                    "Gsol":           line[103:112],
                    "tot,corr":       line[112:121],
                }
                I, J = int(f["I"]), int(f["J"])
                rows.append({
                    "I":     I,
                    "IFRG":  frgs.get(I, "UNK"),
                    "J":     J,
                    "JFRG":  frgs.get(J, "UNK"),
                    "R":     float(f["R"]),
                    "Q":     float(f["QIJ"]),
                    "ENERGY": f["tot,corr"].strip(),
                })
        df = pd.DataFrame(rows, columns=["I", "IFRG", "J", "JFRG", "R", "Q", "ENERGY"])
        df["ENERGY"] = pd.to_numeric(df["ENERGY"], errors="coerce")

    return df
