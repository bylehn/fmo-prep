"""CLI entry points for fmo-prep.

Commands
--------
fmo-prep prep   --config config.yaml
    Runs the full prep-to-GAMESS-input pipeline for a given system.

fmo-prep analyze --log fmo_run.log [--compare fmo_run2.log]
    Parses a GAMESS FMO log file and produces interaction energy outputs.
"""

from __future__ import annotations

import sys
from pathlib import Path

import click

from fmo_prep import __version__


@click.group()
@click.version_option(__version__)
def cli() -> None:
    """fmo-prep: FMO preparation and analysis workflow."""


# ---------------------------------------------------------------------------
# prep
# ---------------------------------------------------------------------------

@cli.command()
@click.option(
    "--config", "-c",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to run config YAML.",
)
def prep(config: Path) -> None:
    """Run the full prep → FragIt → GAMESS input pipeline."""
    from fmo_prep.config import RunConfig

    cfg = RunConfig.from_yaml(config)
    output_dir = cfg.output_path
    output_dir.mkdir(parents=True, exist_ok=True)

    click.echo(f"System type : {cfg.system_type}")
    click.echo(f"Output dir  : {output_dir}")

    # Step 1: structure preparation
    if cfg.system_type == "protein_peptide":
        from fmo_prep.prep.protein_peptide import run as prep_run
    else:
        from fmo_prep.prep.protein_ligand import run as prep_run

    prepared_pdb = prep_run(cfg.inputs, output_dir, cfg)
    click.echo(f"Prepared PDB: {prepared_pdb}")

    # Step 2: render FragIt config and run FragIt
    from fmo_prep.fragit.runner import render_config, run_fragit, find_central_fragment_id

    fragit_dir = output_dir / "fragit"
    fragit_dir.mkdir(exist_ok=True)

    # First pass: run FragIt with centralfragment=0 to generate the .inp and
    # learn the fragment names. Central fragment only matters for mp2_level runs
    # (it sets which fragment gets the NLAYER=2 active region).
    ini_path = render_config(cfg.fragit, central_fragment_id=0, output_path=fragit_dir / "fragit.ini", system_type=cfg.system_type)
    raw_inp = run_fragit(prepared_pdb, ini_path, fragit_dir)
    click.echo(f"Raw FragIt input: {raw_inp}")

    # Second pass: re-render with the correct central fragment ID and re-run.
    # Only needed for mp2_level=True (NLAYER/MPLEVL layers require a central fragment).
    if cfg.fragit.calc_mode == "2layer":
        central_id = find_central_fragment_id(raw_inp, cfg.fragit.central_fragment_resname)
        click.echo(f"Central fragment: {cfg.fragit.central_fragment_resname} → fragment {central_id}")
        ini_path = render_config(cfg.fragit, central_fragment_id=central_id, output_path=fragit_dir / "fragit.ini", system_type=cfg.system_type)
        raw_inp = run_fragit(prepared_pdb, ini_path, fragit_dir)

    # Step 3: patch GAMESS input
    from fmo_prep.fragit.postprocess import patch_inp

    final_inp = output_dir / "fmo_run.inp"
    patch_inp(raw_inp, cfg.fragit, output_path=final_inp)
    click.echo(f"Final GAMESS input: {final_inp}")


# ---------------------------------------------------------------------------
# analyze
# ---------------------------------------------------------------------------

@cli.command()
@click.option(
    "--log", "-l",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="GAMESS FMO log file to analyse.",
)
@click.option(
    "--compare",
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Second GAMESS log for delta (compare - log) analysis.",
)
@click.option(
    "--frag-map",
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Optional complex_fragmapping.txt for residue label mapping.",
)
@click.option(
    "--pdb",
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Optional PDB file for fragment-to-residue mapping (requires --inp).",
)
@click.option(
    "--inp",
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Optional GAMESS .inp file for fragment mapping (requires --pdb).",
)
@click.option(
    "--ligand",
    default=None,
    help="Ligand fragment label hint for auto-detection, e.g. 'LZ1276'.",
)
@click.option(
    "--interaction-mode",
    type=click.Choice(["auto", "chain", "ligand"]),
    default="auto",
    show_default=True,
    help="How to define the interacting partner.",
)
@click.option(
    "--no-pieda",
    is_flag=True,
    default=False,
    help="Parse as non-PIEDA output (total pair energies only).",
)
@click.option(
    "--significant-threshold",
    type=float,
    default=1.0,
    show_default=True,
    help="Min |net energy| (kcal/mol) to include in bar plot.",
)
@click.option(
    "--cov-threshold",
    type=float,
    default=150.0,
    show_default=True,
    help="Suppress covalent energy outliers above this threshold.",
)
@click.option(
    "--output-dir", "-o",
    default=".",
    type=click.Path(file_okay=False, path_type=Path),
    help="Directory for output files.",
)
def analyze(
    log: Path,
    compare: Path | None,
    frag_map: Path | None,
    pdb: Path | None,
    inp: Path | None,
    ligand: str | None,
    interaction_mode: str,
    no_pieda: bool,
    significant_threshold: float,
    cov_threshold: float,
    output_dir: Path,
) -> None:
    """Parse a GAMESS FMO log and produce interaction energy outputs."""
    output_dir.mkdir(parents=True, exist_ok=True)

    from fmo_prep.analysis.parser import parse_gamout
    from fmo_prep.analysis.plots import run_analysis
    from fmo_prep.io.gamess import parse_frag_map_file, write_fragmapping

    # Optionally generate fragment mapping from PDB + inp
    if pdb and inp:
        mapping_path = output_dir / (inp.stem + "_fragmapping.txt")
        write_fragmapping(pdb, inp, mapping_path)
        frag_map = frag_map or mapping_path
        click.echo(f"Fragment mapping: {mapping_path}")

    frag_info = parse_frag_map_file(frag_map) if frag_map else {}

    # Parse log → CSV
    pieda = not no_pieda
    csv_path = output_dir / (log.stem + ("_pieda.csv" if pieda else "_nopieda.csv"))
    df = parse_gamout(log, pieda=pieda)
    df.to_csv(csv_path, index=False)
    click.echo(f"Parsed CSV: {csv_path}")

    compare_csv = None
    if compare:
        compare_csv = output_dir / (compare.stem + ("_pieda.csv" if pieda else "_nopieda.csv"))
        df2 = parse_gamout(compare, pieda=pieda)
        df2.to_csv(compare_csv, index=False)
        click.echo(f"Compare CSV: {compare_csv}")

    # Generate plots and summary
    run_analysis(
        csv_path=csv_path,
        compare_csv=compare_csv,
        frag_info=frag_info,
        ligand_hint=ligand,
        interaction_mode=interaction_mode,
        significant_threshold=significant_threshold,
        cov_threshold=cov_threshold,
        output_dir=output_dir,
    )
    click.echo(f"Outputs written to: {output_dir}")
