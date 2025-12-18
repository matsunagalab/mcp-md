"""Phase 3: Validation Agent for MCP-MD ADK.

This agent validates setup outputs and generates a comprehensive report.
"""

from google.adk.agents import LlmAgent
from google.adk.models.lite_llm import LiteLlm
from google.adk.tools import ToolContext
from google.adk.tools.function_tool import FunctionTool

from mcp_md_adk.config import get_litellm_model
from mcp_md_adk.prompts import get_validation_instruction
from mcp_md_adk.tools.custom_tools import run_validation
from mcp_md_adk.utils import safe_dict, safe_list


# Create a wrapper that reads from session state
def run_validation_tool(tool_context: ToolContext) -> dict:
    """Run validation and generate report. Reads parameters from session state.

    Returns:
        dict: validation_results and final_report
    """
    state = tool_context.state

    return run_validation(
        simulation_brief=safe_dict(state.get("simulation_brief", {})),
        session_dir=str(state.get("session_dir", "")),
        setup_outputs=safe_dict(state.get("outputs", {})),
        decision_log=safe_list(state.get("decision_log", [])),
        compressed_setup=str(state.get("compressed_setup", "")),
    )


def create_validation_agent() -> LlmAgent:
    """Create the Phase 3 validation agent.

    This agent:
    1. Reads setup outputs from session.state
    2. Validates required files exist
    3. Generates a comprehensive markdown report
    4. Saves results to session.state["validation_result"]

    Returns:
        Configured LlmAgent for validation phase
    """
    # Create FunctionTool for validation
    # Note: FunctionTool uses the function's __name__ and docstring automatically
    validation_tool = FunctionTool(run_validation_tool)

    return LlmAgent(
        model=LiteLlm(model=get_litellm_model("compress")),  # Use fast model for validation
        name="validation_agent",
        description="Validates setup outputs and generates report",
        instruction=get_validation_instruction(),
        tools=[validation_tool],
        output_key="validation_result",
    )
