# Example Workflows

## 1. Simple Protein MD System

```bash
mcp-md plan "Create MD system for PDB 1ABC with water and ions"
```

This will generate a workflow plan for:
- Fetching PDB structure
- Cleaning the structure
- Building system with tleap
- Generating OpenMM scripts
- Packaging the system

## 2. Protein-Ligand Complex from SMILES

```bash
mcp-md plan "Build protein-ligand complex MD system. Protein: PDB 1ABC. Ligand: aspirin (SMILES: CC(=O)Oc1ccccc1C(=O)O)"
```

Workflow includes:
- PDB structure retrieval
- SMILES to 3D conversion
- GAFF2 parameterization with AM1-BCC
- smina docking
- System building
- OpenMM MD scripts

## 3. Boltz-2 Structure Prediction + MD

```bash
mcp-md plan "Predict structure from FASTA sequence MKTAYIAKQR... and create MD system"
```

Uses Boltz-2 for:
- De novo structure prediction
- Missing residue completion
- Followed by standard MD system preparation

## 4. Boltz-2 Complex with Affinity Prediction

```bash
mcp-md plan "Predict protein-ligand complex from FASTA and SMILES with binding affinity"
```

Advanced workflow:
- Boltz-2 complex prediction
- Affinity calculation (IC50)
- Ligand parameterization
- MD system setup

## Running Individual Servers

Each MCP server can be run independently:

```bash
# Structure operations
uv run python servers/structure_server.py

# Ligand parameterization
uv run python servers/ligand_server.py

# Docking
uv run python servers/docking_server.py

# System assembly
uv run python servers/assembly_server.py

# MD protocol generation
uv run python servers/protocol_server.py

# Export & packaging
uv run python servers/export_server.py
```

## LM Studio Setup

1. Download and install LM Studio from https://lmstudio.ai/
2. Download a model (recommended: Llama-3.1-8B-Instruct)
3. Start local server (default: http://localhost:1234)
4. Set environment variables:

```bash
export LM_STUDIO_BASE_URL="http://localhost:1234/v1"
export LM_STUDIO_MODEL="llama-3.1-8b-instruct"
```

5. Test connection:

```python
from core.llm_client import LMStudioClient

client = LMStudioClient()
if client.is_available():
    print("LM Studio is ready!")
```

