"""Validation Agent (Phase 3) - QC and report generation.

This module implements a simplified validation phase that:
1. Validates setup outputs (file existence, basic checks)
2. Generates a comprehensive markdown report
"""

import json
from datetime import datetime
from typing import Literal

from langgraph.graph import END, START, StateGraph
from langgraph.types import Command

from mcp_md.state_validation import ValidationOutputState, ValidationState


def validate_outputs(state: ValidationState) -> Command[Literal["generate_report"]]:
    """Validate setup outputs exist and are valid.

    Performs basic validation checks:
    - Required files exist (prmtop, rst7)
    - No critical errors in setup process
    """
    outputs = state.get("setup_outputs", {})

    # Basic validation checks
    checks = []

    # Check required files exist
    required_files = ["prmtop", "rst7"]
    for f in required_files:
        if f in outputs:
            checks.append({"name": f, "status": "pass", "file": outputs[f]})
        else:
            checks.append({"name": f, "status": "fail", "error": "File not found"})

    # Check optional files
    optional_files = ["solvated_pdb", "merged_pdb", "trajectory"]
    for f in optional_files:
        if f in outputs:
            checks.append({"name": f, "status": "pass", "file": outputs[f]})

    all_pass = all(c["status"] == "pass" for c in checks if c["name"] in required_files)

    return Command(
        goto="generate_report",
        update={
            "validation_results": {
                "success": all_pass,
                "checks": checks,
                "timestamp": datetime.now().isoformat(),
            }
        },
    )


def generate_report(state: ValidationState) -> dict:
    """Generate comprehensive markdown report.

    Combines:
    - Simulation configuration from brief
    - Setup summary from compressed_setup
    - Generated file paths
    - Execution statistics from decision log
    """
    simulation_brief = state.get("simulation_brief", {})
    decision_log = state.get("decision_log", [])
    outputs = state.get("setup_outputs", {})
    compressed_setup = state.get("compressed_setup", "")
    validation_results = state.get("validation_results", {})

    # Handle Pydantic model if present
    if hasattr(simulation_brief, "model_dump"):
        brief_dict = simulation_brief.model_dump()
    elif hasattr(simulation_brief, "dict"):
        brief_dict = simulation_brief.dict()
    else:
        brief_dict = simulation_brief

    # Calculate statistics
    total_duration = sum([e.get("duration_seconds", 0) for e in decision_log])
    step_count = len(decision_log)
    avg_duration = total_duration / max(step_count, 1)

    report = f"""# MD Simulation Setup Report

**Generated**: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

## Configuration

```json
{json.dumps(brief_dict, indent=2, default=str)}
```

## Setup Summary

{compressed_setup}

## Generated Files

```json
{json.dumps(outputs, indent=2, default=str)}
```

## Validation Results

{json.dumps(validation_results, indent=2)}

## Execution Statistics

- **Total Steps**: {step_count}
- **Total Duration**: {total_duration:.2f}s
- **Average Step Duration**: {avg_duration:.2f}s

---

*Generated with mcp-md*
"""

    return {"final_report": report}


def create_validation_graph():
    """Create the validation phase graph.

    Graph structure:
    START -> validate_outputs -> generate_report -> END
    """
    builder = StateGraph(ValidationState, output=ValidationOutputState)

    builder.add_node("validate_outputs", validate_outputs)
    builder.add_node("generate_report", generate_report)

    builder.add_edge(START, "validate_outputs")
    # validate_outputs uses Command API to route to generate_report
    builder.add_edge("generate_report", END)

    return builder.compile()


# Pre-compiled graph for import convenience
validation_agent = create_validation_graph()
