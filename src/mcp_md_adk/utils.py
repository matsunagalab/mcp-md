"""Common utilities for MCP-MD ADK.

Re-exports utilities from mcp_md and adds ADK-specific helpers.
"""

# Re-export all utilities from mcp_md
from mcp_md.utils import (
    canonical_tool_name,
    get_today_str,
    compress_tool_result,
    parse_tool_result,
    extract_output_paths,
    format_duration,
    validate_step_prerequisites,
)

# Workflow step definitions (same as mcp_md.state_setup)
SETUP_STEPS = ["prepare_complex", "solvate", "build_topology", "run_simulation"]

STEP_TO_TOOL = {
    "prepare_complex": "prepare_complex",
    "solvate": "solvate_structure",
    "build_topology": "build_amber_system",
    "run_simulation": "run_md_simulation",
}

TOOL_TO_STEP = {v: k for k, v in STEP_TO_TOOL.items()}

STEP_INPUTS = {
    "prepare_complex": "Requires: PDB ID or structure file",
    "solvate": "Requires: merged_pdb from outputs['merged_pdb']",
    "build_topology": "Requires: solvated_pdb, box_dimensions",
    "run_simulation": "Requires: prmtop, rst7",
}


def get_current_step_info(completed_steps: list) -> dict:
    """Get information about the current workflow step.

    Args:
        completed_steps: List of completed step names (may have duplicates)

    Returns:
        Dictionary with current_step, next_tool, step_index, and input_requirements
    """
    # Deduplicate completed steps
    completed_set = set(completed_steps)

    # Find first incomplete step in order
    for i, step in enumerate(SETUP_STEPS):
        if step not in completed_set:
            return {
                "current_step": step,
                "next_tool": STEP_TO_TOOL[step],
                "step_index": i + 1,
                "total_steps": len(SETUP_STEPS),
                "input_requirements": STEP_INPUTS[step],
                "is_complete": False,
            }

    # All steps completed
    return {
        "current_step": None,
        "next_tool": None,
        "step_index": len(SETUP_STEPS),
        "total_steps": len(SETUP_STEPS),
        "input_requirements": "",
        "is_complete": True,
    }


def add_error_recovery_hints(result: dict) -> dict:
    """Add recovery suggestions to failed tool results.

    Args:
        result: Tool result dictionary

    Returns:
        Result dictionary with suggested_action and action_message if failed
    """
    if result.get("success", True):
        return result

    errors = result.get("errors", [])
    error_text = " ".join(str(e).lower() for e in errors)

    if "not found" in error_text or "does not exist" in error_text:
        result["suggested_action"] = "check_previous_step"
        result["action_message"] = "Required file missing. Check if previous step completed."
    elif "timeout" in error_text:
        result["suggested_action"] = "retry_with_longer_timeout"
        result["action_message"] = "Operation timed out. May need longer timeout."
    elif "permission" in error_text:
        result["suggested_action"] = "check_permissions"
        result["action_message"] = "Permission denied. Check file/directory permissions."
    elif "memory" in error_text or "oom" in error_text:
        result["suggested_action"] = "reduce_system_size"
        result["action_message"] = "Out of memory. Consider reducing system size."
    else:
        result["suggested_action"] = "report_and_stop"
        result["action_message"] = "Unrecoverable error. Please check the error details."

    return result


__all__ = [
    # Re-exported from mcp_md
    "canonical_tool_name",
    "get_today_str",
    "compress_tool_result",
    "parse_tool_result",
    "extract_output_paths",
    "format_duration",
    "validate_step_prerequisites",
    # ADK-specific
    "SETUP_STEPS",
    "STEP_TO_TOOL",
    "TOOL_TO_STEP",
    "STEP_INPUTS",
    "get_current_step_info",
    "add_error_recovery_hints",
]
