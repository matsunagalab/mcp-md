# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**MCP-MD** is an AI-powered system for generating molecular dynamics (MD) input files optimized for the Amber/OpenMM ecosystem. It combines:
- **LangGraph 1.0+** for multi-phase workflow orchestration
- **FastMCP** server integration for specialized MD tools
- **Boltz-2** for AI-driven structure prediction
- **AmberTools** for topology generation and parameterization
- **OpenMM** for production-ready MD simulations

## Development Commands

### Environment Setup

```bash
# Create conda environment with scientific packages
conda create -n mcp-md python=3.11
conda activate mcp-md
conda install -c conda-forge openmm rdkit mdanalysis biopython pandas numpy scipy openblas pdbfixer
conda install -c conda-forge ambertools packmol smina

# Install project in editable mode
git clone https://github.com/matsunagalab/mcp-md.git
cd mcp-md
pip install -e .

# Optional: Boltz-2 (for Phase 2-3)
pip install 'boltz[cuda]' --no-deps
pip install torch hydra-core pytorch-lightning einops einx mashumaro modelcif wandb
```

### CLI Usage (main.py)

```bash
# Interactive mode - chat with agent (recommended)
python main.py interactive
python main.py interactive "Setup MD for PDB 1AKE"

# Batch mode - fully automated workflow
python main.py batch "Setup MD for PDB 1AKE in explicit water, 1 ns at 300K"

# Batch with JSON output
python main.py batch "Setup MD for 1AKE" --output-json results.json

# Resume interrupted session
python main.py resume --thread-id md_session_xxxxx

# Phase 1 only (generate SimulationBrief)
python main.py clarify "Setup MD for PDB 1AKE"

# Show available MCP servers
python main.py list-servers

# Show help
python main.py --help
python main.py info
```

### Development Workflow

```bash
# Test individual MCP servers with MCP Inspector
mcp dev servers/structure_server.py
mcp dev servers/genesis_server.py
mcp dev servers/solvation_server.py
mcp dev servers/amber_server.py
mcp dev servers/md_simulation_server.py

# Code quality checks
ruff check src/mcp_md/           # Format checking
ruff check src/mcp_md/ --fix     # Auto-fix format issues
pytest tests/                     # Run tests

# Test in Jupyter notebook
jupyter notebook notebooks/md_agent_v2.ipynb
```

## Architecture

### 3-Phase Workflow

The system follows a **3-phase workflow pattern** adapted from `deep_research_from_scratch`:

1. **Phase 1: Clarification** (`notebooks/1_clarification.ipynb`)
   - Pattern: **ReAct Agent** with MCP tool-calling loop
   - Nodes: `llm_call` → `tool_node` → `combined_router` → `generate_simulation_brief`
   - MCP Tools: `fetch_molecules`, `inspect_molecules` (structure inspection before asking questions)
   - Outputs: Structured `SimulationBrief` (Pydantic model)
   - Implementation: `src/mcp_md/clarification_agent.py`, `src/mcp_md/state_scope.py`

2. **Phase 2: Setup** (`notebooks/2_setup_agent.ipynb`, `notebooks/3_setup_coordinator.ipynb`)
   - Pattern: **Coordinator-Tools** (supervisor pattern)
   - Coordinator node selects tools, Tools node executes them
   - Fixed skeleton: structure_fetch → repair → ligand_param → complex_generation → qc_check
   - Implementation: `src/mcp_md/setup_coordinator.py`, `src/mcp_md/state_setup.py`

3. **Phase 3: Validation** (`notebooks/4_validation.ipynb`)
   - QC checks, format conversion, report generation
   - Implementation: `src/mcp_md/validation_agent.py` (planned)

### Directory Structure

```
mcp-md/
├── notebooks/              # 🎯 PRIMARY development location (source of truth)
│   ├── 1_clarification.ipynb       # Phase 1
│   ├── 2_setup_agent.ipynb         # Phase 2 basic
│   ├── 3_setup_coordinator.ipynb   # Phase 2 advanced
│   ├── 4_validation.ipynb          # Phase 3
│   ├── 5_full_agent.ipynb          # End-to-end integration
│   └── utils.py                    # Rich formatting for notebooks
│
├── src/mcp_md/            # 🚫 Auto-generated source (DO NOT EDIT DIRECTLY)
│   ├── clarification_agent.py      # Phase 1 implementation
│   ├── setup_coordinator.py        # Phase 2 implementation
│   ├── validation_agent.py         # Phase 3 implementation
│   ├── state_scope.py              # Phase 1 state definitions
│   ├── state_setup.py              # Phase 2 state definitions
│   ├── prompts.py                  # LLM prompt templates
│   └── mcp_integration.py          # MCP client setup
│
├── servers/               # FastMCP servers (5 independent servers)
│   ├── structure_server.py         # PDB fetch, clean, parameterize
│   ├── genesis_server.py           # Boltz-2 structure prediction
│   ├── solvation_server.py         # Solvent/membrane embedding
│   ├── amber_server.py             # Amber topology generation
│   └── md_simulation_server.py     # OpenMM MD execution & analysis
│
├── common/                # Shared utilities
│   ├── base.py                     # BaseToolWrapper for external tools
│   └── utils.py                    # Common utilities (logging, etc.)
│
└── checkpoints/           # LangGraph state persistence
```

### FastMCP Servers

The system uses **5 independent FastMCP servers**, each providing specialized MD preparation tools:

1. **structure_server.py** - Structure acquisition and preparation
   - `fetch_molecules()`: Download from PDB/AlphaFold/PDB-REDO
   - `inspect_molecules()`: Analyze structure, classify chains
   - `clean_protein()`: PDBFixer + protonation + pdb4amber
   - `clean_ligand()`: Template matching + geometry optimization
   - `run_antechamber_robust()`: GAFF2 + AM1-BCC parameterization
   - `prepare_complex()`: All-in-one preparation pipeline

2. **genesis_server.py** - AI structure prediction
   - Boltz-2 integration for FASTA → PDB conversion

3. **solvation_server.py** - Solvation and membrane setup
   - packmol-memgen for solvent/membrane embedding

4. **amber_server.py** - Topology generation
   - tleap for prmtop and inpcrd file generation

5. **md_simulation_server.py** - MD execution
   - OpenMM simulation and trajectory analysis

Each server is **independently testable** using `mcp dev servers/<server_name>.py`.

## Development Pattern: Direct Python Files

**Source of truth**: `src/mcp_md/` Python files (NOT notebooks)

### Development Workflow

1. Edit Python files directly in `src/mcp_md/`
2. Test with notebooks (e.g., `notebooks/md_agent_v2.ipynb`)
3. Run `ruff check src/mcp_md/` for linting

### Notebooks Are For:
- Testing and demos only
- NOT for code generation via `%%writefile`

### File Structure

```
src/mcp_md/
├── config.py               # Configuration (env vars with MCPMD_ prefix)
├── prompts.py              # All prompt templates
├── utils.py                # Common utilities (parse_tool_result, extract_output_paths,
│                           #   validate_step_prerequisites)
├── state_scope.py          # Phase 1 state definitions
├── state_setup.py          # Phase 2 state definitions (SETUP_STEPS, reducers)
├── state_validation.py     # Phase 3 state definitions
├── state_full.py           # Full agent state
├── clarification_agent.py  # Phase 1 implementation
├── setup_agent.py          # Phase 2 ReAct agent (STEP_TO_TOOL, workflow tracking,
│                           #   error recovery suggestions)
├── validation_agent.py     # Phase 3 implementation
├── mcp_integration.py      # MCP client factory
└── full_agent.py           # 3-phase integration

common/
├── base.py                 # BaseToolWrapper + timeout configuration functions
│                           #   (get_default_timeout, get_solvation_timeout, etc.)
├── errors.py               # Unified error handling (create_error_result, etc.)
└── utils.py                # Shared utilities (create_unique_subdir, etc.)

tests/
├── test_utils.py           # Unit tests for mcp_md.utils module
└── test_integration.py     # Integration tests
```

### Testing

```bash
# Unit tests
pytest tests/ -v

# Test specific file
pytest tests/test_structure_server.py -v

# Test with coverage
pytest tests/ --cov=src/mcp_md --cov-report=html

# Import verification (quick sanity check)
python -c "from mcp_md.config import settings; from mcp_md.utils import parse_tool_result; from common.errors import create_error_result; print('All imports OK')"

# Test full workflow in notebook
jupyter notebook notebooks/md_agent_v2.ipynb

# Or run main.py for CLI testing
python main.py
```

## Key Technical Patterns

### LangGraph 1.0+ Patterns

This project uses **LangGraph 1.0+** advanced features:

- **ReAct Pattern** (Phase 1): Tool-calling loop with conditional routing
  ```python
  # Graph: START → llm_call → combined_router → tool_node (loop) → END
  def combined_router(state: AgentState) -> Literal["tool_node", "generate_simulation_brief", "__end__"]:
      last_message = state["messages"][-1]
      if hasattr(last_message, "tool_calls") and last_message.tool_calls:
          return "tool_node"  # Execute MCP tools
      return route_after_llm(state)  # Check if ready to proceed
  ```

- **MessagesState**: Proper message accumulation with `add_messages` reducer
  ```python
  class AgentState(MessagesState):
      messages: Annotated[list, add_messages]
      simulation_brief: Optional[SimulationBrief]
  ```

- **Subgraphs**: Phase 2 Setup as independent subgraph with typed input/output
  ```python
  setup_subgraph = StateGraph(SetupState, input=AgentState, output=SetupOutputState)
  ```

- **Persistence**: SqliteSaver for checkpoint storage and state replay
  ```python
  checkpointer = SqliteSaver.from_conn_string("checkpoints/clarification.sqlite")
  ```

### Structured Output Schemas

All LLM decisions use **Pydantic schemas** for type safety:

```python
class ClarifyWithUser(BaseModel):
    """Clarification decision schema."""
    need_clarification: bool
    question: str
    verification: str

class SimulationBrief(BaseModel):
    """Structured MD setup parameters."""
    pdb_id: Optional[str]
    fasta_sequence: Optional[str]
    ligand_smiles: str
    simulation_type: Literal["equilibration", "production"]
    # ... other parameters
```

Key schemas:
- `ClarifyWithUser`: Clarification routing decision
- `SimulationBrief`: Structured MD parameters
- `ExecuteSetupStep`: Task delegation to tools (Coordinator pattern)
- `SetupComplete`: Completion signal

### Coordinator-Tools Pattern (Phase 2)

Phase 2 uses the **Coordinator-Tools pattern** (supervisor pattern from deep_research):

**Structure**: 2 nodes
1. **Coordinator node** (decision maker)
   - Selects which tool to execute next
   - Uses `ExecuteSetupStep` or `SetupComplete` structured outputs
   - Records decisions in decision log

2. **Tools node** (executor)
   - Executes the selected MCP tool
   - Returns results to coordinator
   - Loops back to coordinator until `SetupComplete`

**Key difference from deep_research**: Sequential execution (not parallel)

```python
async def setup_coordinator(state: SetupState) -> Command[Literal["setup_tools"]]:
    """Coordinator decides next action."""
    setup_tools = [ExecuteSetupStep, SetupComplete, think_tool]
    model_with_tools = model.bind_tools(setup_tools)
    response = await model_with_tools.ainvoke([...])
    return Command(goto="setup_tools", update={"setup_messages": [response]})

async def setup_tools(state: SetupState) -> Command[Literal["setup_coordinator", "__end__"]]:
    """Tools execute and return to coordinator."""
    # Execute tool, log decision
    if isinstance(last_message, SetupComplete):
        return Command(goto=END, update={...})
    return Command(goto="setup_coordinator", update={...})
```

### State Management

**State hierarchy** (follows deep_research pattern):

- `AgentInputState`: User input interface (minimal fields)
- `AgentState`: Main multi-phase state (messages, brief, decision logs)
- `SetupState`: Phase 2 subgraph state (extends AgentState)
- `SetupOutputState`: Phase 2 output interface

**Reducers** for proper state merging:
- `add_messages`: For message lists (LangGraph standard)
- `operator.add`: For list accumulators (decision_log, completed_steps, raw_notes)

### Workflow Step Tracking (Phase 2)

Phase 2 setup_agent uses explicit step tracking for reliable workflow execution:

```python
# state_setup.py
SETUP_STEPS = ["prepare_complex", "solvate", "build_topology", "run_simulation"]

# setup_agent.py
STEP_TO_TOOL = {
    "prepare_complex": "prepare_complex",
    "solvate": "solvate_structure",
    "build_topology": "build_amber_system",
    "run_simulation": "run_md_simulation",
}
TOOL_TO_STEP = {v: k for k, v in STEP_TO_TOOL.items()}  # Reverse mapping
```

The prompt includes `<Workflow_Progress>` section:
```
CURRENT STEP: 2 of 4 - solvate
NEXT TOOL TO CALL: solvate_structure
INPUT REQUIREMENTS: Requires: merged_pdb from outputs['merged_pdb']
```

### Unified Error Handling

`common/errors.py` provides standardized error formats for LLM understanding:

```python
from common.errors import create_error_result, create_file_not_found_error

# In MCP tool implementation
try:
    # operation
except FileNotFoundError as e:
    return create_file_not_found_error(file_path, "PDB file")

# Generic errors with recovery hints
return create_error_result(
    error,
    hints=["Check file path", "Verify format"],
    context={"file": path}
)
```

### Tool Result Parsing

`utils.py` provides safe result parsing for MCP tools:

```python
from mcp_md.utils import parse_tool_result, extract_output_paths, validate_step_prerequisites

# Parse any result format (dict, str JSON, other)
result = parse_tool_result(raw_result)

# Extract file paths including box_dimensions and ligand_params
outputs = extract_output_paths(result)
# Returns: {"merged_pdb": "...", "box_dimensions": {...}, "ligand_params": [...]}

# Validate prerequisites before executing a step
is_valid, errors = validate_step_prerequisites("solvate", outputs)
# Checks: merged_pdb exists for solvate, solvated_pdb + box_dimensions for build_topology, etc.
```

### Error Recovery (Phase 2)

`setup_agent.py` automatically adds error recovery suggestions to failed tool results:

```python
# In tool_node, failed results get suggested_action and action_message:
if not result.get("success", True):
    # Categorize error and add recovery hints
    if "not found" in error_text:
        result["suggested_action"] = "check_previous_step"
        result["action_message"] = "Required file missing. Check if previous step completed."
    elif "timeout" in error_text:
        result["suggested_action"] = "retry_with_longer_timeout"
        result["action_message"] = "Operation timed out. May need longer timeout."
```

### Timeout Configuration

`common/base.py` provides configurable timeouts via environment variables:

```python
from common.base import get_default_timeout, get_solvation_timeout, get_membrane_timeout

# Each operation has appropriate timeout
tleap_timeout = get_default_timeout()        # MCPMD_DEFAULT_TIMEOUT (300s)
solvation_timeout = get_solvation_timeout()  # MCPMD_SOLVATION_TIMEOUT (600s)
membrane_timeout = get_membrane_timeout()    # MCPMD_MEMBRANE_TIMEOUT (1800s)
```

### MCP Integration

**FastMCP** provides type-safe tool integration:

- Type-safe schema generation from Python function signatures
- Each server is independent and testable with `mcp dev`
- Standard error format: `{"success": bool, "errors": [], "warnings": [], "output_file": ...}`
- LLM-friendly error messages with recovery hints

**MCP client setup** in `src/mcp_md/mcp_integration.py`:
```python
from langchain_mcp_adapters import create_mcp_client

async def setup_mcp_clients():
    clients = []
    for server in ["structure_server", "genesis_server", ...]:
        client = await create_mcp_client(f"servers/{server}.py")
        clients.append(client)
    return clients
```

## Code Quality Standards

### Format Configuration

- **Line length**: 100 characters (ruff/black)
- **Python target**: 3.11+
- **Pydantic version**: 2.12+
- **Docstring format**: D212 (summary on first line)
- **Import order**: stdlib → third-party → local

### Quality Checklist

After editing Python files in `src/mcp_md/`:
1. Run `ruff check src/mcp_md/`
2. Fix any linting errors directly in the Python file
3. Test in notebooks (e.g., `notebooks/md_agent_v2.ipynb`)
4. Verify imports work correctly

## Important References

- **ARCHITECTURE.md**: Detailed technical architecture (26k+ tokens, very comprehensive)
- **AGENTS.md**: Cursor AI Agent guidelines and development workflow
- **.cursor/rules/project-rules.md**: Critical development rules
- **.cursor/rules/notebook-development.md**: Notebook development patterns
- **deep_research_from_scratch**: Reference implementation pattern (in project for reference)

## Configuration

### Environment Variables

All settings can be configured via `MCPMD_` prefixed environment variables:

```bash
# .env or shell exports
export MCPMD_OUTPUT_DIR="./output"
export MCPMD_CLARIFICATION_MODEL="anthropic:claude-haiku-4-5-20251001"
export MCPMD_SETUP_MODEL="anthropic:claude-sonnet-4-20250514"
export MCPMD_COMPRESS_MODEL="anthropic:claude-haiku-4-5-20251001"
export MCPMD_DEFAULT_TIMEOUT=300
export MCPMD_SOLVATION_TIMEOUT=600
export MCPMD_MEMBRANE_TIMEOUT=1800
export MCPMD_MD_SIMULATION_TIMEOUT=3600
export MCPMD_MAX_MESSAGE_HISTORY=6
```

Access in code:
```python
from mcp_md.config import settings

model_name = settings.setup_model
output_dir = settings.output_dir
```

### API Keys

```bash
OPENAI_API_KEY=...
ANTHROPIC_API_KEY=...
OLLAMA_BASE_URL=http://localhost:11434
```

### Default Models

- **Clarification (Phase 1)**: `anthropic:claude-haiku-4-5-20251001` (fast, cheap)
- **Setup (Phase 2)**: `anthropic:claude-sonnet-4-20250514` (balanced)
- **Compression**: `anthropic:claude-haiku-4-5-20251001` (fast)
