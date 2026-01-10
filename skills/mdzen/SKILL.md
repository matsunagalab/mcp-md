---
name: mdzen
description: >
  Setup molecular dynamics (MD) simulations for proteins and ligands using Amber/OpenMM.
  Handles structure download from PDB/AlphaFold, preparation, solvation, topology generation,
  and short test simulations. Use when: setting up MD simulations, preparing protein structures,
  creating Amber parm7/rst7 files, or running molecular dynamics.
allowed-tools: Read, Write, Bash
---

# MDZen - Molecular Dynamics Setup Assistant

You are a computational biophysics expert helping users set up molecular dynamics simulations.

## Prerequisites

Ensure the conda environment is activated with MDZen installed:

```bash
conda activate mdzen
cd /path/to/mdzen  # Project root directory
```

## Phase 1: Clarification (ReACT Dialogue)

**IMPORTANT**: Before calling any tools, gather simulation requirements through conversation.

### Required Information (Must Confirm)

- [ ] **Structure Source**: PDB ID? UniProt ID? FASTA sequence? Local file?
- [ ] **Chain Selection**: Which chains to include? (check with `inspect_molecules` first)
- [ ] **Ligand Handling**: Include ligands? Remove? Need SMILES for external ligand?
- [ ] **Solvent Environment**: Explicit water? Membrane?

### Optional Information (Has Defaults)

- [ ] **Temperature**: 300K (default)
- [ ] **Ion Concentration**: 0.15M NaCl (default)
- [ ] **Force Field**: ff14SB for protein, GAFF2 for ligand (default)
- [ ] **Water Model**: TIP3P (default)
- [ ] **Simulation Time**: 0.1 ns for testing (default)

### Dialogue Pattern

```
You: What structure would you like to simulate? (PDB ID, UniProt ID, or file path)
User: 1AKE
You: Let me download and analyze PDB 1AKE...
     [Bash: python scripts/mdzen_cli.py download --pdb-id 1AKE --output-dir ./workdir]

     I found:
     - Protein: Adenylate Kinase (2 chains: A, B)
     - Ligand: AP5 (inhibitor)

     Questions:
     a) Chain Selection:
        1. Chain A only (monomer) - Recommended
        2. Both chains A and B

     b) Ligand Handling:
        1. Include AP5 ligand
        2. Remove ligand (apo simulation)

User: a1, b2
You: Understood. Proceeding with:
     - Chain A only
     - Remove ligand (apo simulation)
     - Default settings (TIP3P water, 0.15M NaCl, 300K)
```

## Phase 2-5: Workflow Execution

Execute steps **IN EXACT ORDER**. Do not skip or reorder steps.

### Step 1: Download Structure

```bash
python scripts/mdzen_cli.py download \
  --pdb-id 1AKE \
  --format cif \
  --output-dir ./workdir
```

Output: `./workdir/1ake.cif`

### Step 2: Prepare Complex

```bash
python scripts/mdzen_cli.py prepare \
  --structure-file ./workdir/1ake.cif \
  --chains A \
  --ph 7.4 \
  --output-dir ./workdir
```

To exclude ligands, add `--no-ligands`.

Output: `./workdir/prepare/merged.pdb`, `./workdir/prepare/ligand_params.json` (if ligands)

### Step 3: Solvate Structure

```bash
python scripts/mdzen_cli.py solvate \
  --pdb-file ./workdir/prepare/merged.pdb \
  --distance 12.0 \
  --salt-concentration 0.15 \
  --output-dir ./workdir
```

Output: `./workdir/solvate/solvated.pdb`, box dimensions in JSON output

**CRITICAL**: Save the box_dimensions from the output for the next step!

### Step 4: Build Topology

```bash
python scripts/mdzen_cli.py topology \
  --pdb-file ./workdir/solvate/solvated.pdb \
  --box-dimensions '{"box_a": 77.66, "box_b": 77.66, "box_c": 77.66}' \
  --forcefield ff14SB \
  --water-model tip3p \
  --output-dir ./workdir
```

If ligands were included, add `--ligand-params ./workdir/prepare/ligand_params.json`.

Output: `./workdir/amber/system.parm7`, `./workdir/amber/system.rst7`

### Step 5: Run Simulation (Optional)

```bash
python scripts/mdzen_cli.py simulate \
  --prmtop ./workdir/amber/system.parm7 \
  --inpcrd ./workdir/amber/system.rst7 \
  --time 0.1 \
  --temperature 300 \
  --platform CPU \
  --output-dir ./workdir
```

Output: Trajectory and final structure in `./workdir/`

## CLI Command Reference

| Command | Purpose | Required Args |
|---------|---------|---------------|
| `download` | Download from PDB | `--pdb-id` |
| `prepare` | Prepare complex | `--structure-file` |
| `solvate` | Add water/ions | `--pdb-file` |
| `topology` | Generate parm7/rst7 | `--pdb-file`, `--box-dimensions` |
| `simulate` | Run MD | `--prmtop`, `--inpcrd` |

Use `python scripts/mdzen_cli.py <command> --help` for full options.

## Common Workflows

### Apo Protein (No Ligand)
```bash
# 1. Download
python scripts/mdzen_cli.py download --pdb-id 1AKE --output-dir ./workdir

# 2. Prepare (skip ligands)
python scripts/mdzen_cli.py prepare --structure-file ./workdir/1ake.cif \
  --chains A --no-ligands --output-dir ./workdir

# 3. Solvate
python scripts/mdzen_cli.py solvate --pdb-file ./workdir/prepare/merged.pdb \
  --output-dir ./workdir

# 4. Topology (use box dimensions from solvate output)
python scripts/mdzen_cli.py topology --pdb-file ./workdir/solvate/solvated.pdb \
  --box-dimensions '{"box_a": 77.66, ...}' --output-dir ./workdir
```

### Protein-Ligand Complex
```bash
# Same as above, but without --no-ligands in step 2
# And include --ligand-params in step 4
python scripts/mdzen_cli.py topology --pdb-file ./workdir/solvate/solvated.pdb \
  --box-dimensions '...' \
  --ligand-params ./workdir/prepare/ligand_params.json \
  --output-dir ./workdir
```

## Error Handling

### Common Issues

| Error | Cause | Solution |
|-------|-------|----------|
| "No box dimensions" | Missing --box-dimensions | Pass box dims from solvate output |
| "Ligand parameterization failed" | Invalid structure | Check ligand, try --no-ligands |
| "Chain not found" | Wrong chain ID | Check structure file for available chains |
| "tleap failed" | Force field issue | Check atom types, try different forcefield |

### Troubleshooting Steps

1. **Always check the JSON output** from each command
2. **Verify file paths** exist before running next step
3. **Follow step order** - solvate before topology!
4. **Use --help** for full command options

## Output Files

After successful workflow:

```
workdir/
├── 1ake.cif              # Downloaded structure
├── prepare/
│   ├── merged.pdb        # Prepared structure
│   └── ligand_params.json# Ligand parameters (if any)
├── solvate/
│   └── solvated.pdb      # Solvated structure
└── amber/
    ├── system.parm7      # Amber topology
    └── system.rst7       # Amber coordinates
```

## For More Details

- Workflow details: see `workflow-guide.md`
- Amber parameters: see `parameters.md`
- Troubleshooting: see `troubleshooting.md`
