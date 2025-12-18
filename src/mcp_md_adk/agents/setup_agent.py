"""Phase 2: Setup Agent for MCP-MD ADK.

This agent executes the 4-step MD setup workflow:
1. prepare_complex - Structure preparation and ligand parameterization
2. solvate - Add water box
3. build_topology - Generate Amber prmtop/rst7
4. run_simulation - Execute MD with OpenMM
"""

from google.adk.agents import LlmAgent
from google.adk.models.lite_llm import LiteLlm
from google.adk.tools import ToolContext
from google.adk.tools.function_tool import FunctionTool

from mcp_md_adk.config import get_litellm_model
from mcp_md_adk.prompts import get_setup_instruction
from mcp_md_adk.tools.mcp_setup import get_setup_tools
from mcp_md_adk.tools.custom_tools import get_workflow_status


# Create a wrapper that reads from session state
def get_workflow_status_tool(tool_context: ToolContext) -> dict:
    """Get current workflow progress and validate prerequisites. Call this before each step.

    Returns:
        dict: Current step info, completed steps, and validation status
    """
    state = tool_context.state
    completed_steps = list(state.get("completed_steps", []))
    outputs = dict(state.get("outputs", {}))

    return get_workflow_status(completed_steps, outputs)


def create_setup_agent() -> LlmAgent:
    """Create the Phase 2 setup agent.

    This agent:
    1. Reads SimulationBrief from session.state["simulation_brief"]
    2. Executes 4-step workflow using MCP tools
    3. Tracks progress via completed_steps in session.state
    4. Stores outputs in session.state["outputs"]

    Returns:
        Configured LlmAgent for setup phase
    """
    # Get all MCP tools for setup workflow
    mcp_tools = get_setup_tools()

    # Create FunctionTool for workflow status
    # Note: FunctionTool uses the function's __name__ and docstring automatically
    status_tool = FunctionTool(get_workflow_status_tool)

    # Combine all tools
    all_tools = mcp_tools + [status_tool]

    return LlmAgent(
        model=LiteLlm(model=get_litellm_model("setup")),
        name="setup_agent",
        description="Executes 4-step MD setup workflow",
        instruction=get_setup_instruction(),
        tools=all_tools,
        output_key="setup_result",  # Saves final result to session.state
    )


def update_state_after_tool(tool_name: str, result: dict, ctx) -> None:
    """Callback to update session state after tool execution.

    Updates completed_steps and outputs based on tool results.

    Args:
        tool_name: Name of the executed tool
        result: Tool result dictionary
        ctx: Invocation context with session access
    """
    from mcp_md_adk.utils import (
        canonical_tool_name,
        extract_output_paths,
        add_error_recovery_hints,
        TOOL_TO_STEP,
    )

    state = ctx.session.state

    # Add error recovery hints if failed
    if not result.get("success", True):
        result = add_error_recovery_hints(result)

    # Extract output paths and update outputs
    outputs = state.get("outputs", {})
    outputs.update(extract_output_paths(result))
    state["outputs"] = outputs

    # Track completed steps
    canonical_name = canonical_tool_name(tool_name)
    if canonical_name in TOOL_TO_STEP and result.get("success"):
        completed = state.get("completed_steps", [])
        step_name = TOOL_TO_STEP[canonical_name]
        if step_name not in completed:
            completed.append(step_name)
        state["completed_steps"] = completed

    # Log decision
    from datetime import datetime
    from mcp_md_adk.utils import compress_tool_result

    decision_log = state.get("decision_log", [])
    decision_log.append({
        "tool": tool_name,
        "result": compress_tool_result(tool_name, result),
        "timestamp": datetime.now().isoformat(),
    })
    state["decision_log"] = decision_log
