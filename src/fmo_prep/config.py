"""Pydantic models for fmo-prep run configuration.

Loaded from a YAML file via RunConfig.from_yaml(path).

Example usage::

    config = RunConfig.from_yaml("config.yaml")
    print(config.system_type)  # "protein_ligand"
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, Field, model_validator


class PrepConfig(BaseModel):
    """System-preparation parameters."""

    # protein_peptide
    cut_radius: float = Field(5.0, description="Radius (Å) around selection for truncation")

    # protein_ligand
    amber_level: Literal["low", "high"] = Field(
        "low", description="Amber minimisation level (low=1000 cycles, high=5000)"
    )

    # capping
    true_nterm: list[str] = Field(
        default_factory=list,
        description="Chain:resnum pairs that are genuine N-termini (skip ACE cap), e.g. ['A:1']",
    )
    true_cterm: list[str] = Field(
        default_factory=list,
        description="Chain:resnum pairs that are genuine C-termini (skip NME cap), e.g. ['A:300']",
    )


class FragitConfig(BaseModel):
    """Parameters controlling FragIt fragmentation and GAMESS input generation."""

    central_fragment_resname: str = Field(
        description="Residue name of the ligand/peptide chain to use as central fragment, e.g. 'LZ1276'"
    )

    # FragIt .ini output options
    boundaries: float = Field(2.0, description="Boundary distance (Å) for FragIt output layers")
    basis: str = Field("6-31G*", description="Basis set string passed to FragIt (e.g. '6-31G*', '6-31G(d)')")
    use_atom_names: bool = Field(False, description="FragIt useatomnames setting")

    # FMO calculation options
    calc_mode: Literal["hf", "mp2", "2layer"] = Field(
        "mp2",
        description=(
            "Calculation mode:\n"
            "  hf     - HF only, single FragIt pass\n"
            "  mp2    - MP2 on entire system (PIEDA), single FragIt pass\n"
            "  2layer - MP2 at active site (layer 2), HF elsewhere (layer 1); "
            "requires central_fragment_resname, two FragIt passes"
        ),
    )
    implicit_solvent: bool = Field(
        False,
        description="Add PCM implicit solvent ($PCM SOLVNT=WATER IFMO=1 ICOMP=0 $END)",
    )
    nbody: int = Field(2, description="GAMESS NBODY setting")

    # Postprocessing: GAMESS $FMO block parameters
    resdim: float = Field(2.0, description="GAMESS RESDIM cutoff (Å)")
    rcorsd: float = Field(2.0, description="GAMESS RCORSD cutoff (Å)")

    # Postprocessing: GAMESS $SYSTEM / $GDDI
    mwords: int = Field(125, description="$SYSTEM MWORDS")
    ngroup: int = Field(1, description="$GDDI NGROUP")


class RunConfig(BaseModel):
    """Top-level run configuration loaded from a YAML file."""

    system_type: Literal["protein_ligand", "protein_peptide"]
    inputs: dict[str, str] = Field(
        description="Logical name → file path. Required key: 'complex'. "
        "protein_ligand also requires 'ligand' (SDF or PDB)."
    )
    output_dir: str = Field("outputs", description="Directory for all generated outputs")
    prep: PrepConfig = Field(default_factory=PrepConfig)
    fragit: FragitConfig

    @model_validator(mode="after")
    def _check_inputs(self) -> RunConfig:
        required = {"complex"}
        if self.system_type == "protein_ligand":
            required.add("ligand")
        missing = required - self.inputs.keys()
        if missing:
            raise ValueError(f"Missing required input(s): {missing}")
        return self

    @classmethod
    def from_yaml(cls, path: str | Path) -> RunConfig:
        """Load and validate configuration from a YAML file."""
        raw: dict[str, Any] = yaml.safe_load(Path(path).read_text())
        return cls.model_validate(raw)

    @property
    def output_path(self) -> Path:
        return Path(self.output_dir)
