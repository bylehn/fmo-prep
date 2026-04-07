# fmo-prep

Command-line tool for preparing GAMESS FMO input files and analysing FMO output.
Handles two system types:

- **protein_ligand** — small-molecule ligand bound to a protein (full complex, Amber minimisation)
- **protein_peptide** — peptide bound to a protein (truncated by radius, H-link capped)

## Installation

Requires [pixi](https://prefix.dev/docs/pixi/overview).

```bash
git clone <repo>
cd fmo-prep
pixi install
```

This installs all dependencies (Python ≥ 3.11, RDKit, AmberTools, FragIt, parmed, …) into
an isolated conda environment. The `fmo-prep` CLI is available inside that environment.

## Quick start

```bash
cd examples/protein_ligand
pixi run fmo-prep prep --config config.yaml

cd examples/protein_peptide
pixi run fmo-prep prep --config config.yaml
```

The final GAMESS input is written to `outputs/fmo_run.inp`.

## Running your own system

Create a working directory with a `config.yaml` and your input files, then run from that directory:

```bash
mkdir ~/fmo-runs/my_system && cd ~/fmo-runs/my_system
# copy or symlink your PDB (and SDF for protein_ligand) here
pixi run --manifest-path /path/to/fmo-prep/pixi.toml fmo-prep prep --config config.yaml
```

Or add fmo-prep's pixi environment to your `PATH`:

```bash
export PATH="/path/to/fmo-prep/.pixi/envs/default/bin:$PATH"
fmo-prep prep --config config.yaml
```

## Config file

### protein_ligand

```yaml
system_type: protein_ligand

inputs:
  complex: complex.pdb   # protein + ligand, any standard PDB
  ligand: ligand.sdf     # ligand structure for Antechamber (optional; used for charge)

output_dir: outputs

prep:
  amber_level: low       # low = 1000 minimisation cycles; high = 5000

fragit:
  central_fragment_resname: LIG   # 3-letter residue name of the ligand
  boundaries: 2.0                 # FragIt boundary distance (Å)
  basis: "6-31G*"                 # basis set for $BASIS and $FMOBND
  calc_mode: "mp2"                # hf | mp2 | 2layer
  implicit_solvent: false         # true adds $PCM SOLVNT=WATER block
  nbody: 2
  resdim: 2.0
  rcorsd: 2.0
  mwords: 125                     # $SYSTEM MWORDS
  ngroup: 1                       # $GDDI NGROUP
```

**Pipeline:** load PDB → strip H/water → parameterise ligand with Antechamber (GAFF2 + BCC) →
build Amber topology (ff14SB) with TLeap → energy minimise with sander → run FragIt →
patch GAMESS header blocks → `outputs/fmo_run.inp`

### protein_peptide

```yaml
system_type: protein_peptide

inputs:
  complex: complex.pdb   # full protein–peptide complex

output_dir: outputs

prep:
  cut_radius: 5.0        # keep all residues within this distance (Å) of the peptide
  true_nterm: []         # list genuine N-termini to skip capping, e.g. ["A:1"]
  true_cterm: []         # list genuine C-termini to skip capping, e.g. ["A:300"]

fragit:
  central_fragment_resname: B   # chain ID of the peptide (for 2layer mode)
  boundaries: 5.0
  basis: "6-31G(d)"
  calc_mode: "mp2"
  implicit_solvent: false
  nbody: 2
  resdim: 2.0
  rcorsd: 2.0
  mwords: 125
  ngroup: 1
```

**Pipeline:** load PDB → strip water → select all residues within `cut_radius` Å of the peptide →
cap backbone cuts with H-link atoms → run FragIt → patch GAMESS header blocks → `outputs/fmo_run.inp`

### calc_mode

| Mode | Description | NLAYER | MPLEVL |
|------|-------------|--------|--------|
| `hf` | HF throughout | 1 | 0 |
| `mp2` | MP2 throughout (PIEDA) | 1 | 2 |
| `2layer` | MP2 at ligand/peptide, HF elsewhere | 2 | 0,2 |

For `2layer` mode, `central_fragment_resname` must match the fragment that should be the active region.

## Analysing GAMESS output

```bash
fmo-prep analyze \
  --log fmo_run.log \
  --output-dir results/
```

With fragment mapping for residue labels:

```bash
fmo-prep analyze \
  --log fmo_run.log \
  --pdb outputs/lig_minimised.pdb \
  --inp outputs/fmo_run.inp \
  --ligand LZ1276 \
  --output-dir results/
```

Delta analysis (compare two runs):

```bash
fmo-prep analyze \
  --log run_A/fmo_run.log \
  --compare run_B/fmo_run.log \
  --output-dir results/delta/
```

Outputs written to `--output-dir`:

- `*_pieda.csv` — per-pair interaction energies (ES, EX, CT+mix, DI, total)
- `pieda_bar.png` — bar chart of significant interactions with the ligand/peptide
- `pieda_heatmap.png` — full interaction energy heatmap
- `*_summary.txt` — per-residue net energies

### analyze options

| Option | Default | Description |
|--------|---------|-------------|
| `--log` | required | GAMESS FMO log file |
| `--compare` | — | Second log for delta mode |
| `--pdb` | — | PDB for fragment→residue mapping |
| `--inp` | — | GAMESS .inp for fragment mapping (requires `--pdb`) |
| `--frag-map` | — | Pre-computed `*_fragmapping.txt` file |
| `--ligand` | auto | Ligand fragment label hint (e.g. `LZ1276`) |
| `--interaction-mode` | `auto` | `auto` \| `chain` \| `ligand` |
| `--no-pieda` | false | Parse as non-PIEDA output (pair energies only) |
| `--significant-threshold` | 1.0 | Min \|energy\| (kcal/mol) for bar plot |
| `--cov-threshold` | 150.0 | Suppress covalent outliers above this value |
| `--output-dir` | `.` | Directory for output files |

## Running tests

```bash
pixi run test
```

## Output directory layout

After `fmo-prep prep`:

```
outputs/
├── capped.pdb                  # protein_peptide: truncated + H-capped structure
├── lig_minimised.pdb           # protein_ligand: Amber-minimised complex
├── lig_prepared.pdb            # protein_ligand: pre-minimisation prepared structure
├── lig.prmtop / lig.inpcrd     # Amber topology + coordinates
├── ligands/
│   ├── LIG.mol2                # ligand GAFF2 parameters
│   └── LIG.frcmod
├── fragit/
│   ├── fragit.ini              # rendered FragIt config
│   └── *.inp                   # raw FragIt output
└── fmo_run.inp                 # final patched GAMESS input
```
