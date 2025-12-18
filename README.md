# MCP-MD: Molecular Dynamics Input File Generation Agent

An AI agent system specialized for Amber-based MD input file generation. Built with LangGraph + FastMCP using a 3-phase workflow (Clarification → Setup → Validation).

## Features

- **LangGraph Integration**: Stateful workflows, persistence, human feedback
  - StateGraph-based implementation compliant with LangChain 1.0
  - Official MCP integration via langchain-mcp-adapters
  - Checkpoint functionality for pause/resume
- **ReAct Pattern**: Phase 1 pre-inspects PDB structures before generating appropriate questions
  - Analyze structures with `fetch_molecules`/`inspect_molecules` tools
  - Auto-detect multi-chain structures and ligand presence
  - Simple single-chain proteins proceed automatically
- **Boltz-2 Integration**: High-accuracy structure prediction and binding affinity prediction from FASTA/SMILES
- **AmberTools Complete**: No external QM software required for ligand parameterization (AM1-BCC charge calculation)
- **FastMCP Integration**: Modular 5 independent servers, type-safe automatic schema generation
- **OpenMM Dedicated**: Python-programmable production-ready script generation

## 📚 Documentation

- **[ARCHITECTURE.md](ARCHITECTURE.md)** - Project architecture, implementation plan, technical specifications
- **[CLAUDE.md](CLAUDE.md)** - Claude Code guidance and development patterns
- **[AGENTS.md](AGENTS.md)** - Cursor AI Agent settings and guidelines
- **[.cursor/rules/](.cursor/rules/)** - Project rules and development workflow

## Installation

### Prerequisites

- Python 3.11 or higher
- [conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/)
- GPU recommended (for Boltz-2, OpenMM acceleration)

### Steps

#### 1. Set up conda environment

```bash
# Create conda environment
conda create -n mcp-md python=3.11
conda activate mcp-md

# Install scientific computing packages
conda install -c conda-forge openmm rdkit mdanalysis biopython pandas numpy scipy openblas pdbfixer

# MD preparation tools
conda install -c conda-forge ambertools packmol smina
```

#### 2. Install Python packages

```bash
# Clone the project
git clone https://github.com/matsunagalab/mcp-md.git
cd mcp-md

# Install package (editable mode)
pip install -e .
```

#### 3. Install Boltz-2 (Optional)

Boltz-2 is used in Phase 2-3 (Setup/Validation). Install when needed:

```bash
# If you have a CUDA-compatible GPU
pip install 'boltz[cuda]' --no-deps

# Then install missing dependencies individually
pip install torch hydra-core pytorch-lightning einops einx mashumaro modelcif wandb

# Or downgrade scipy first then do normal install
conda install -c conda-forge scipy=1.13.1
pip install 'boltz[cuda]'
```

> **Note**: One of Boltz-2's dependencies (fairscale) strictly requires scipy==1.13.1, which may conflict with scipy already installed via conda. Using the `--no-deps` option preserves existing packages while adding only missing ones.

#### 4. Install Ollama (Optional)
Ollama is a local LLM execution environment. By default, the system uses Ollama's `gpt-oss:20b` model.

```bash
# For Mac
brew install ollama
brew pull gpt-oss:20b
brew services start ollama
```

## Usage

### CLI (main.py)

```bash
# Interactive mode - setup while chatting with agent (recommended)
python main.py interactive
python main.py interactive "Setup MD for PDB 1AKE"

# Batch mode - fully automated workflow execution
python main.py batch "Setup MD for PDB 1AKE in explicit water, 1 ns at 300K"

# Batch processing with JSON output
python main.py batch "Setup MD for 1AKE" --output-json results.json

# Resume interrupted session
python main.py resume --thread-id md_session_xxxxx

# Phase 1 only (SimulationBrief generation)
python main.py clarify "Setup MD for PDB 1AKE"

# List MCP servers
python main.py list-servers

# Help
python main.py --help
python main.py info
```

### Notebook Development

```bash
jupyter notebook notebooks/md_agent_v2.ipynb
```

### MCP Server Testing

Each FastMCP server can be tested independently:

```bash
# Launch MCP Inspector (Structure Server example)
mcp dev servers/structure_server.py

# Test other servers
mcp dev servers/genesis_server.py
mcp dev servers/solvation_server.py
mcp dev servers/amber_server.py
mcp dev servers/md_simulation_server.py
```

### MCP Server List

| Server | Description |
|--------|-------------|
| `structure_server` | Structure retrieval from PDB/AlphaFold/PDB-REDO, chain separation, structure repair, ligand GAFF2 parameterization |
| `genesis_server` | Structure prediction from FASTA sequences via Boltz-2 (monomer/multimer support) |
| `solvation_server` | Solvation (water box) and lipid membrane embedding via packmol-memgen |
| `amber_server` | Amber topology (parm7) and coordinate (rst7) file generation via tleap |
| `md_simulation_server` | MD execution with OpenMM, trajectory analysis with MDTraj |

## Directory Structure

```
mcp-md/
├── main.py               # CLI entry point
│
├── src/mcp_md/           # Source code (edit directly)
│   ├── config.py                   # Configuration management (env var support)
│   ├── prompts.py                  # Prompt templates
│   ├── utils.py                    # Utilities
│   ├── state_scope.py              # Phase 1 state definitions
│   ├── state_setup.py              # Phase 2 state definitions
│   ├── state_validation.py         # Phase 3 state definitions
│   ├── state_full.py               # Integrated state definitions
│   ├── clarification_agent.py      # Phase 1: ReAct Agent (structure inspection → questions)
│   ├── setup_agent.py              # Phase 2: ReAct Setup Agent
│   ├── validation_agent.py         # Phase 3: Validation & Report
│   ├── mcp_integration.py          # MCP integration
│   └── full_agent.py               # 3-phase integration
│
├── notebooks/            # For testing and demos
│   ├── 1_clarification.ipynb       # Phase 1 test
│   ├── md_agent_v2.ipynb           # Integration test
│   └── test_*.ipynb                # MCP server tests
│
├── servers/              # FastMCP servers (5 servers)
│   ├── structure_server.py         # PDB retrieval, structure repair, ligand GAFF2 parameterization
│   ├── genesis_server.py           # Boltz-2 structure generation (FASTA → PDB)
│   ├── solvation_server.py         # Solvation and membrane embedding (packmol-memgen)
│   ├── amber_server.py             # Amber topology/coordinate generation (tleap)
│   └── md_simulation_server.py     # MD execution and analysis (OpenMM/MDTraj)
│
├── common/               # Shared libraries
│   ├── base.py                     # BaseToolWrapper
│   ├── errors.py                   # Unified error handling
│   └── utils.py                    # Common utilities
│
├── checkpoints/          # LangGraph checkpoints
├── ARCHITECTURE.md       # Detailed architecture
├── CLAUDE.md             # Claude Code guidance
├── AGENTS.md             # Cursor AI Agent settings
└── README.md             # This file
```

## Development Workflow

### Direct Python Files

This project adopts the **Direct Python Files** pattern:

```
✅ Edit src/mcp_md/ directly
✅ Test and demo in notebooks/
✅ Format check with ruff check src/mcp_md/

🚫 Code generation via %%writefile is not recommended
```

### Code Formatting

```bash
# Format check
ruff check src/mcp_md/

# Auto-fix
ruff check src/mcp_md/ --fix
```

### Test Execution

```bash
# Run unit tests
pytest tests/ -v

# Run specific test file
pytest tests/test_structure_server.py -v

# Run with coverage
pytest tests/ --cov=src/mcp_md --cov-report=html

# Import test (verify new modules)
python -c "from mcp_md.config import settings; from mcp_md.utils import parse_tool_result; print('OK')"
```

## Configuration (Environment Variables)

Settings can be customized via `MCPMD_` prefixed environment variables:

```bash
# Set via .env file or environment variables
export MCPMD_OUTPUT_DIR="./custom_output"
export MCPMD_CLARIFICATION_MODEL="anthropic:claude-haiku-4-5-20251001"
export MCPMD_SETUP_MODEL="anthropic:claude-sonnet-4-20250514"
export MCPMD_DEFAULT_TIMEOUT=300
export MCPMD_MAX_MESSAGE_HISTORY=6
```

Available settings:

| Environment Variable | Default | Description |
|---------------------|---------|-------------|
| `MCPMD_OUTPUT_DIR` | `output` | Output directory |
| `MCPMD_CLARIFICATION_MODEL` | `anthropic:claude-haiku-4-5-20251001` | Phase 1 model |
| `MCPMD_SETUP_MODEL` | `anthropic:claude-sonnet-4-20250514` | Phase 2 model |
| `MCPMD_COMPRESS_MODEL` | `anthropic:claude-haiku-4-5-20251001` | Compression model |
| `MCPMD_DEFAULT_TIMEOUT` | `300` | Default timeout (seconds) |
| `MCPMD_MD_SIMULATION_TIMEOUT` | `3600` | MD execution timeout (seconds) |
| `MCPMD_MAX_MESSAGE_HISTORY` | `6` | Number of message history to retain |

## License

MIT License

## Citations

When using this tool, please cite the following:

### Boltz-2

```
S. Passaro et al., Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction.
bioRxiv (2025). doi:10.1101/2025.06.14.659707
```

### AmberTools

```
D. A. Case et al., AmberTools, J. Chem. Inf. Model. 63, 6183 (2023).
```

### OpenMM

```
P. Eastman et al., OpenMM 8: Molecular Dynamics Simulation with Machine Learning Potentials,
J. Phys. Chem. B 128, 109 (2024).
```
