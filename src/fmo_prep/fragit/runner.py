"""FragIt invocation and Jinja2 config rendering.

FragIt is called as:
    fragit --use-config <ini_file> -o <output.inp> <pdb_file>

The raw .inp output from FragIt contains:
  - $BASIS ... $END
  - $FMO ... $END  (with NFRAG, NBODY, RESDIM=2.0, RCORSD=2.0, ICHARG, FRGNAM, INDAT)
                   (plus NLAYER / MPLEVL when mp2level=True in the .ini)
  - $DATA ... $END
  - (optionally) $FMOXYZ ... $END

It does NOT contain $SYSTEM, $GDDI, $SCF, $CONTRL, or $FMOPRP — those are
added by postprocess.patch_inp.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape

from fmo_prep.config import FragitConfig


def render_config(
    cfg: FragitConfig,
    central_fragment_id: int,
    output_path: Path,
) -> Path:
    """Render a FragIt .ini config file from the Jinja2 base template.

    Args:
        cfg: FragIt configuration values.
        central_fragment_id: Fragment index for centralfragment in the .ini.
            Pass 0 to leave it unset (fragit will not apply boundary layers).
        output_path: Where to write the rendered .ini file.

    Returns:
        Path to the written .ini file.
    """
    templates_dir = Path(__file__).parent / "templates"
    env = Environment(
        loader=FileSystemLoader(str(templates_dir)),
        autoescape=select_autoescape([]),
        trim_blocks=True,
        lstrip_blocks=True,
    )
    template = env.get_template("base.ini.j2")
    rendered = template.render(
        boundaries=cfg.boundaries,
        central_fragment_id=central_fragment_id,
        use_atom_names=cfg.use_atom_names,
        basis=cfg.basis,
        calc_mode=cfg.calc_mode,
    )
    output_path = Path(output_path)
    output_path.write_text(rendered)
    return output_path


def run_fragit(pdb_path: Path, config_path: Path, output_dir: Path) -> Path:
    """Run FragIt on a prepared PDB file.

    Args:
        pdb_path: Path to the prepared complex PDB (absolute).
        config_path: Path to the rendered .ini config (absolute).
        output_dir: Directory where the output .inp file will be written.

    Returns:
        Path to the generated GAMESS .inp file.

    Raises:
        RuntimeError: If FragIt exits with a non-zero return code.
        FileNotFoundError: If the expected output .inp is not found.
    """
    pdb_path = Path(pdb_path).resolve()
    config_path = Path(config_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    out_inp = output_dir / (pdb_path.stem + ".inp")

    cmd = [
        "fragit",
        "--use-config", str(config_path),
        "-o", str(out_inp),
        str(pdb_path),
    ]

    result = subprocess.run(
        cmd,
        cwd=str(output_dir),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"FragIt failed (exit {result.returncode}):\n{result.stdout}"
        )

    if not out_inp.exists():
        raise FileNotFoundError(
            f"FragIt ran successfully but expected output not found: {out_inp}\n"
            f"FragIt output:\n{result.stdout}"
        )

    return out_inp


def find_central_fragment_id(inp_path: Path, resname: str) -> int:
    """Return the 1-based fragment index whose FRGNAM starts with *resname*.

    Args:
        inp_path: Path to a FragIt-generated GAMESS .inp file.
        resname: Residue name prefix to search for (case-insensitive).

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
        f"Available names (first 10): {frag_names[:10]}"
    )
