"""Phase 1: Clarification Agent for MCP-MD ADK.

This agent handles user interaction to gather MD simulation requirements
and generates a structured SimulationBrief.
"""

from google.adk.agents import LlmAgent
from google.adk.models.lite_llm import LiteLlm
from google.adk.tools.function_tool import FunctionTool

from mcp_md_adk.config import get_litellm_model
from mcp_md_adk.prompts import get_clarification_instruction
from mcp_md_adk.tools.mcp_setup import get_clarification_tools
from mcp_md_adk.tools.custom_tools import generate_simulation_brief


def create_clarification_agent() -> LlmAgent:
    """Create the Phase 1 clarification agent.

    This agent:
    1. Uses fetch_molecules/inspect_molecules to analyze structures
    2. Asks clarification questions based on inspection
    3. Generates SimulationBrief via generate_simulation_brief tool
    4. Saves result to session.state["simulation_brief"] via output_key

    Returns:
        Configured LlmAgent for clarification phase
    """
    # Get MCP tools for structure inspection
    mcp_tools = get_clarification_tools()

    # Create FunctionTool for generating SimulationBrief
    generate_brief_tool = FunctionTool(generate_simulation_brief)

    # Combine all tools
    all_tools = mcp_tools + [generate_brief_tool]

    return LlmAgent(
        model=LiteLlm(model=get_litellm_model("clarification")),
        name="clarification_agent",
        description="Gathers MD simulation requirements and generates SimulationBrief",
        instruction=get_clarification_instruction(),
        tools=all_tools,
        output_key="simulation_brief",  # Saves to session.state["simulation_brief"]
    )
