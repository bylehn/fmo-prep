# fmo-prep — Project Status

## Goal

A single installable Python package (`fmo-prep`) that consolidates two FMO (Fragment Molecular Orbital)
preparation workflows and all downstream analysis scripts into a unified CLI:

```
fmo-prep prep     --config config.yaml   # structure prep → FragIt → GAMESS .inp
fmo-prep analyze  --log fmo_run.log      # parse GAMESS log → CSVs + plots
```

Production runs live **outside** the repo, e.g.:
```
~/fmo-runs/2026-03-31_kinase_imatinib/
├── config.yaml
├── inputs/complex.pdb
└── outputs/fmo_run.inp
```

---

## Source material extracted / ported

| Source | Module | Status |
|--------|--------|--------|
| `bison/services/simulation-setup-service/app/task.py` | `prep/protein_ligand.py` | ✅ Done |
| `fmo-poc/structures/processed/cap_termini.py` | `prep/common.py` (geometry helpers, parse_term_flags) | ✅ Done (ACE/NME deferred) |
| `PDBProgress/pdbprogress/capping/cut_and_cap.py` | `prep/common.py` (cut_and_cap, H-link caps) | ✅ Done |
| `PDBProgress/pdbprogress/cleaning/cleaning.py` | `prep/common.py` (remove_water, remove_by_residue, remove_chain) | ✅ Done |
| `fmo-poc/scripts/gamout.py` + `gamout_nopieda.py` | `analysis/parser.py` | ⬜ Stub only |
| `fmo-poc/scripts/plot_pieda.py` | `analysis/plots.py` | ⬜ Stub only |
| `fmo-poc/scripts/plot_pieda_delta.py` | `analysis/plots.py` (delta mode) | ⬜ Stub only |
| `fmo-poc/scripts/map_fragments.py` | `io/gamess.py` | ⬜ Stub only |
| `fmo-poc/inputs/fn001/debug_fragit.py` | `io/gamess.py` (validate_fragments) | ⬜ Stub only |
| `fmo-poc/inputs/benchmarks/data/myconfig.ini` | `fragit/templates/base.ini.j2` | ✅ Done |

---

## Repo structure

```
fmo-prep/
├── pixi.toml                          # conda-forge env (python≥3.11, rdkit, ambertools, fragit, parmed, …)
├── pyproject.toml                     # hatchling build; entry point: fmo-prep = fmo_prep.cli:cli
├── src/fmo_prep/
│   ├── __init__.py
│   ├── cli.py                         # ✅ Click entry points: prep + analyze (all options wired)
│   ├── config.py                      # ✅ Pydantic v2 models: RunConfig, FragitConfig, PrepConfig
│   ├── prep/
│   │   ├── common.py                  # ✅ cut_and_cap (H-link caps), cleaning helpers
│   │   ├── protein_peptide.py         # ✅ load → clean → cut_and_cap → save capped.pdb
│   │   └── protein_ligand.py          # ✅ AntechamberWrapper, TLeapWrapper, minimisation pipeline
│   ├── fragit/
│   │   ├── runner.py                  # ⬜ render_config, run_fragit stubs
│   │   ├── postprocess.py             # ⬜ patch_inp stub
│   │   └── templates/base.ini.j2      # ✅ Jinja2 FragIt config template
│   ├── analysis/
│   │   ├── parser.py                  # ⬜ parse_gamout stub
│   │   ├── plots.py                   # ⬜ process_csv, run_analysis stubs
│   │   └── reports.py                 # ⬜ write_summary stub
│   └── io/
│       ├── gamess.py                  # ⬜ parse_inp_file, write_fragmapping stubs
│       ├── pdb.py                     # ✅ thin parmed wrappers
│       └── sdf.py                     # ✅ RDKit SDF loader
├── tests/
│   ├── test_fragit/test_postprocess.py  # skipped, documents expected interface
│   ├── test_fragit/test_runner.py       # skipped
│   └── test_analysis/test_parser.py    # skipped
└── examples/
    ├── protein_ligand/config.yaml
    └── protein_peptide/config.yaml
```

---

## What works now

### CLI

```bash
fmo-prep --help
fmo-prep prep --help
fmo-prep analyze --help
```

### Config loading

```python
from fmo_prep.config import RunConfig
cfg = RunConfig.from_yaml("examples/protein_ligand/config.yaml")
```

- Validates `system_type`, required inputs, and all nested fields
- Missing `ligand` for `protein_ligand` raises a clear `ValueError`

### Prep — `protein_peptide`

```python
from fmo_prep.prep.protein_peptide import run
run({"complex": "complex.pdb"}, Path("outputs"), cfg)
```

1. Loads PDB with parmed
2. Strips water (HOH)
3. `cut_and_cap`: selects all residues within `cut_radius` Å of `central_fragment_resname`,
   caps each backbone C–N cut with an H link atom (HLN / HLC) at equilibrium bond distance
   (C–H = 1.09 Å, N–H = 0.99 Å)
4. Saves `outputs/capped.pdb`

**Functional test** (CDK2 benchmark, `:LZ1`, 5 Å radius):
→ 4569 atoms → 246 atoms, 14 residues, 16 H-link cap atoms

### Prep — `protein_ligand`

```python
from fmo_prep.prep.protein_ligand import run
run({"complex": "complex.pdb"}, Path("outputs"), cfg)
```

1. Loads PDB, identifies non-standard residues
2. Strips H atoms for Amber consistency; fixes HIS naming; fixes halogen atom names
3. Antechamber (GAFF2 + BCC charges) → mol2; parmchk2 → frcmod; aligns atom names in complex
4. TLeap: ff14SB + GAFF2 + TIP3P → `.prmtop` / `.inpcrd`
5. Sander minimisation (`low` = 1000 cycles, `high` = 5000 cycles)
6. Converts minimised coords back to PDB → `outputs/<prefix>_minimised.pdb`

---

## Design decisions

- **No PDBProgress dependency** — `cut_and_cap` and cleaning functions copied into `prep/common.py`
  (PDBProgress is Python 3.10, fmo-prep targets ≥3.11)
- **ACE/NME capping deferred** — H-link atoms only for now; `cap_termini.py` logic will be
  added to `prep/common.py` later
- **postprocess.py uses detect-and-patch** for NLAYER/MPLEVL — handles both FragIt versions
  (some write these, some don't)
- **Analysis decoupled** — `fmo-prep analyze` works on any GAMESS FMO log, no config required

---

## Next steps (in priority order)

1. **`fragit/runner.py`** — implement `render_config` (Jinja2 → .ini) and `run_fragit` (subprocess)
2. **`fragit/postprocess.py`** — implement `patch_inp`: prepend `$SYSTEM/$GDDI/$SCF/$CONTRL/$FMOPRP`
   blocks and inject `RESDIM`/`RCORSD`/`NLAYER`/`MPLEVL` into `$FMO` block
3. **`io/gamess.py`** — port `parse_inp_file` and `write_fragmapping` from `map_fragments.py`
4. **`analysis/parser.py`** — port unified `parse_gamout` from `gamout.py` + `gamout_nopieda.py`
5. **`analysis/plots.py`** — port `process_csv` and all plotting from `plot_pieda.py` + `plot_pieda_delta.py`
6. **`analysis/reports.py`** — port `write_summary`
7. **End-to-end test** — run full pipeline on CDK2 benchmark, diff `.inp` against reference

### Config notes for production runs

- `protein_ligand`: set `central_fragment_resname` to the 3-letter residue name of your ligand
  (e.g. `LIG`, `LZ1`) — this is used both for FragIt `centralfragment` and for the output file prefix
- `protein_peptide`: set `central_fragment_resname` to the chain ID or residue selection of the peptide
  (e.g. `:B` as an Amber mask, or just `B` — the code prepends `:` automatically if needed)
