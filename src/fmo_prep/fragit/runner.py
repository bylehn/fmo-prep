"""FragIt invocation and Jinja2 config rendering.

Functions
---------
render_config   Render a FragIt .ini file from a Jinja2 template + FragitConfig.
run_fragit      Run the fragit CLI and return the generated .inp path.
find_central_fragment_id
                Look up the fragment number for a given residue name in a .inp file.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from jinja2 import Environment, PackageLoader

from fmo_prep.config import FragitConfig


def render_config(
    cfg: FragitConfig,
    central_fragment_id: int,
    output_path: Path,
) -> Path:
    """Render a FragIt .ini config file from the Jinja2 base template.

    Args:
        cfg: FragIt configuration values.
        central_fragment_id: Integer fragment index for the central fragment
            (0 = let FragIt auto-detect / not set).
        output_path: Where to write the rendered .ini file.

    Returns:
        Path to the written .ini file.
    """
    raise NotImplementedError("fragit/runner.render_config not yet implemented")


def run_fragit(pdb_path: Path, config_path: Path, output_dir: Path) -> Path:
    """Run FragIt on a prepared PDB file.

    Args:
        pdb_path: Path to the prepared complex PDB.
        config_path: Path to the rendered .ini config.
        output_dir: Directory where FragIt writes its output.

    Returns:
        Path to the generated GAMESS .inp file.

    Raises:
        RuntimeError: If FragIt exits with a non-zero return code.
    """
    raise NotImplementedError("fragit/runner.run_fragit not yet implemented")


def find_central_fragment_id(inp_path: Path, resname: str) -> int:
    """Return the 1-based fragment index whose FRGNAM matches *resname*.

    Args:
        inp_path: Path to a FragIt-generated GAMESS .inp file.
        resname: Residue name to search for (case-insensitive prefix match).

    Returns:
        Fragment index (1-based).

    Raises:
        ValueError: If no matching fragment is found.
    """
    from fmo_prep.io.gamess import parse_inp_file

    frag_names, _ = parse_inp_file(inp_path)
    target = resname.strip().upper()
    for i, name in enumerate(frag_names, start=1):
        if name.upper().startswith(target):
            return i
    raise ValueError(
        f"Fragment '{resname}' not found in {inp_path}. "
        f"Available: {frag_names[:10]}{'...' if len(frag_names) > 10 else ''}"
    )
