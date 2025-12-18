"""Validation Agent (Phase 3) - QC and report generation.

This module implements a simplified validation phase that:
1. Validates setup outputs (file existence, basic checks)
2. Generates a comprehensive markdown report
"""

import json
from datetime import datetime
from pathlib import Path
from typing import Literal

from langgraph.graph import END, START, StateGraph
from langgraph.types import Command

from mcp_md.state_validation import ValidationOutputState, ValidationState


def _check_file_exists(file_path: str) -> bool:
    """Check if a file path exists on disk.

    Args:
        file_path: Path string to check

    Returns:
        True if file exists, False otherwise
    """
    if not file_path:
        return False
    return Path(file_path).exists()


def validate_outputs(state: ValidationState) -> Command[Literal["generate_report"]]:
    """Validate setup outputs exist and are valid.

    Performs validation checks:
    - Required files: key exists AND file exists on disk
    - Optional files: warn if key exists but file doesn't exist

    This ensures we don't falsely report success when files are missing.
    """
    outputs = state.get("setup_outputs", {})

    # Basic validation checks
    checks = []

    # Check required files (key must exist AND file must exist on disk)
    # parm7 is the raw output from amber_server, prmtop is the mapped key
    required_files = {
        "topology": ["prmtop", "parm7"],  # Either key is valid
        "coordinates": ["rst7"],
    }

    for name, valid_keys in required_files.items():
        found = False
        found_key = None
        found_path = None

        for key in valid_keys:
            if key in outputs and outputs[key]:
                found_key = key
                found_path = outputs[key]
                # Check if file actually exists on disk
                if _check_file_exists(found_path):
                    checks.append({
                        "name": name,
                        "key": key,
                        "status": "pass",
                        "file": found_path
                    })
                    found = True
                    break
                else:
                    # Key exists but file doesn't - record as fail
                    checks.append({
                        "name": name,
                        "key": key,
                        "status": "fail",
                        "file": found_path,
                        "error": f"File path exists in outputs but file not found on disk: {found_path}"
                    })
                    found = True  # Don't check other keys, we found the entry
                    break

        if not found:
            checks.append({
                "name": name,
                "keys_checked": valid_keys,
                "status": "fail",
                "error": f"No valid file path found (checked keys: {valid_keys})"
            })

    # Check optional files (warn if key exists but file doesn't)
    optional_files = ["solvated_pdb", "merged_pdb", "trajectory", "trajectory_file"]
    for f in optional_files:
        if f in outputs and outputs[f]:
            file_path = outputs[f]
            if _check_file_exists(file_path):
                checks.append({"name": f, "status": "pass", "file": file_path})
            else:
                # Optional file key exists but file missing - warn (not fail)
                checks.append({
                    "name": f,
                    "status": "warn",
                    "file": file_path,
                    "warning": f"Optional file path in outputs but file not found: {file_path}"
                })

    # Only required files determine overall pass/fail
    required_names = set(required_files.keys())
    all_pass = all(
        c["status"] == "pass"
        for c in checks
        if c["name"] in required_names
    )

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
    - Generated file paths (organized by session directory)
    - Execution statistics from decision log
    """
    simulation_brief = state.get("simulation_brief", {})
    decision_log = state.get("decision_log", [])
    outputs = state.get("setup_outputs", {})
    compressed_setup = state.get("compressed_setup", "")
    validation_results = state.get("validation_results", {})
    session_dir = state.get("session_dir", outputs.get("session_dir", ""))

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

    # Format validation checks
    validation_summary = []
    checks = validation_results.get("checks", [])
    for check in checks:
        status = check.get("status", "unknown")
        name = check.get("name", "unknown")
        if status == "pass":
            validation_summary.append(f"- ✓ **{name}**: {check.get('file', 'OK')}")
        elif status == "warn":
            validation_summary.append(f"- ⚠ **{name}**: {check.get('warning', 'Warning')}")
        else:
            validation_summary.append(f"- ✗ **{name}**: {check.get('error', 'Failed')}")
    validation_text = "\n".join(validation_summary) if validation_summary else "No validation checks performed"

    # Format output files with session directory emphasis
    output_files_text = _format_output_files(outputs, session_dir)

    report = f"""# MD Simulation Setup Report

**Generated**: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

## Session Directory

All output files are organized under this session directory:

```
{session_dir}/
```

## Configuration

```json
{json.dumps(brief_dict, indent=2, default=str)}
```

## Setup Summary

{compressed_setup}

## Generated Files

{output_files_text}

## Validation Results

{validation_text}

**Overall Status**: {"✓ All required files present" if validation_results.get("success") else "✗ Some required files missing"}

## Execution Statistics

- **Total Steps**: {step_count}
- **Total Duration**: {total_duration:.2f}s
- **Average Step Duration**: {avg_duration:.2f}s

---

*Generated with mcp-md*
"""

    return {"final_report": report}


def _format_output_files(outputs: dict, session_dir: str) -> str:
    """Format output files showing their location relative to session directory.

    Args:
        outputs: Dictionary of output file paths
        session_dir: Session directory path

    Returns:
        Formatted markdown string showing file structure
    """
    if not outputs:
        return "No output files generated."

    # Categorize outputs
    categories = {
        "Structure Preparation": ["merged_pdb", "output_dir"],
        "Solvation": ["solvated_pdb", "box_dimensions"],
        "Topology": ["prmtop", "parm7", "rst7"],
        "Simulation": ["trajectory", "trajectory_file", "energy_file", "final_structure"],
        "Ligand Parameters": ["ligand_params"],
    }

    lines = []
    session_dir_str = str(session_dir) if session_dir else ""

    # Show session directory structure
    if session_dir_str:
        lines.append(f"**Session Directory**: `{session_dir_str}`")
        lines.append("")
        lines.append("```")
        lines.append(f"{session_dir_str.split('/')[-1]}/")

    # Group files by category
    categorized = {cat: [] for cat in categories}
    other_files = []

    for key, value in outputs.items():
        if key == "session_dir" or value is None:
            continue

        found_category = False
        for cat, keys in categories.items():
            if key in keys:
                categorized[cat].append((key, value))
                found_category = True
                break
        if not found_category:
            other_files.append((key, value))

    # Build tree structure
    for cat, files in categorized.items():
        if files:
            for key, value in files:
                # Show relative path if possible
                if session_dir_str and isinstance(value, str) and session_dir_str in value:
                    rel_path = value.replace(session_dir_str + "/", "")
                    lines.append(f"├── {rel_path}")
                elif isinstance(value, dict):
                    lines.append(f"├── {key}: (dict)")
                elif isinstance(value, list):
                    lines.append(f"├── {key}: ({len(value)} items)")
                else:
                    lines.append(f"├── {key}: {value}")

    if other_files:
        for key, value in other_files:
            if isinstance(value, str) and session_dir_str and session_dir_str in value:
                rel_path = value.replace(session_dir_str + "/", "")
                lines.append(f"├── {rel_path}")
            else:
                lines.append(f"├── {key}: {value}")

    if session_dir_str:
        lines.append("```")

    # Also show full paths for easy access
    lines.append("")
    lines.append("### Full Paths")
    lines.append("")
    for key, value in outputs.items():
        if key == "session_dir" or value is None:
            continue
        if isinstance(value, str):
            lines.append(f"- **{key}**: `{value}`")
        elif isinstance(value, dict):
            lines.append(f"- **{key}**: `{json.dumps(value)}`")
        elif isinstance(value, list):
            lines.append(f"- **{key}**: {len(value)} item(s)")

    return "\n".join(lines)


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
