"""Full 3-Phase MD Setup Agent for MCP-MD ADK.

This module integrates all three phases using ADK's SequentialAgent:
- Phase 1: Clarification (requirements gathering)
- Phase 2: Setup (MCP tool orchestration)
- Phase 3: Validation (QC and report generation)
"""

from google.adk.agents import SequentialAgent

from mcp_md_adk.agents.clarification_agent import create_clarification_agent
from mcp_md_adk.agents.setup_agent import create_setup_agent
from mcp_md_adk.agents.validation_agent import create_validation_agent


def create_full_agent() -> SequentialAgent:
    """Create the full 3-phase MD setup agent.

    The SequentialAgent orchestrates:
    1. clarification_agent → outputs simulation_brief
    2. setup_agent → executes 4-step workflow
    3. validation_agent → validates and generates report

    State flows through session.state between agents.

    Returns:
        Configured SequentialAgent for complete MD workflow
    """
    return SequentialAgent(
        name="full_md_agent",
        description="Complete MD simulation setup workflow with 3 phases",
        sub_agents=[
            create_clarification_agent(),
            create_setup_agent(),
            create_validation_agent(),
        ],
    )


def create_clarification_only_agent() -> SequentialAgent:
    """Create an agent that only runs Phase 1 (clarification).

    Useful for interactive mode where user reviews SimulationBrief
    before proceeding to setup.

    Returns:
        SequentialAgent with only clarification phase
    """
    return SequentialAgent(
        name="clarification_only_agent",
        description="Phase 1: Gather MD simulation requirements",
        sub_agents=[
            create_clarification_agent(),
        ],
    )


def create_setup_validation_agent() -> SequentialAgent:
    """Create an agent that runs Phase 2-3 (setup + validation).

    Useful for continuing after interactive clarification review.

    Returns:
        SequentialAgent with setup and validation phases
    """
    return SequentialAgent(
        name="setup_validation_agent",
        description="Phase 2-3: Execute setup and validate outputs",
        sub_agents=[
            create_setup_agent(),
            create_validation_agent(),
        ],
    )
