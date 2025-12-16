"""Full MD Setup Agent - 3-Phase Integration.

This module integrates all three phases of MD setup into a single LangGraph workflow:
- Phase 1: Clarification (requirements gathering)
- Phase 2: Setup (MCP tool orchestration)
- Phase 3: Validation (QC and report generation)

Uses SqliteSaver for checkpoint persistence and human-in-the-loop feedback.
"""

import sqlite3
from pathlib import Path
from typing import Optional

from langchain_core.messages import HumanMessage
from langgraph.checkpoint.sqlite import SqliteSaver
from langgraph.graph import END, START, StateGraph

# Phase 1
from mcp_md.clarification_agent import create_clarification_graph

# Phase 2
from mcp_md.setup_agent import create_setup_graph

# Phase 3
from mcp_md.validation_agent import create_validation_graph

from mcp_md.state_full import (
    FullAgentInputState,
    FullAgentOutputState,
    FullAgentState,
)


# =============================================================================
# STATE TRANSFORMATION NODES
# =============================================================================


def prepare_setup_input(state: FullAgentState) -> dict:
    """Transform clarification output to setup input.

    Converts SimulationBrief from Phase 1 to the format expected by
    the setup agent (Phase 2).
    """
    simulation_brief = state.get("simulation_brief", {})

    # Handle Pydantic model if present
    if hasattr(simulation_brief, "model_dump"):
        brief_dict = simulation_brief.model_dump()
    elif hasattr(simulation_brief, "dict"):
        brief_dict = simulation_brief.dict()
    else:
        brief_dict = simulation_brief

    return {
        "simulation_brief": brief_dict,
        "setup_messages": [HumanMessage(content="Starting MD setup")],
        "completed_steps": [],
        "current_step_index": 0,
        "outputs": {},
    }


def prepare_validation_input(state: FullAgentState) -> dict:
    """Transform setup output to validation input.

    Maps the output fields from Phase 2 to the input format
    expected by the validation agent (Phase 3).
    """
    return {
        "simulation_brief": state.get("simulation_brief", {}),
        "setup_outputs": state.get("outputs", {}),
        "decision_log": state.get("decision_log", []),
        "compressed_setup": state.get("compressed_setup", ""),
    }


# =============================================================================
# GRAPH CONSTRUCTION
# =============================================================================


def create_full_agent(
    checkpoint_path: Optional[Path] = None, interrupt_after_clarification: bool = True
):
    """Create the full 3-phase MD setup agent.

    Args:
        checkpoint_path: Path to SQLite checkpoint file.
                        Defaults to "checkpoints/full_workflow.db"
        interrupt_after_clarification: Whether to pause after Phase 1 for
                                      user review of SimulationBrief.
                                      Defaults to True.

    Returns:
        Compiled LangGraph with checkpointing enabled
    """
    # Build main graph
    builder = StateGraph(
        FullAgentState, input=FullAgentInputState, output=FullAgentOutputState
    )

    # Phase subgraphs
    clarification_graph = create_clarification_graph()
    setup_graph = create_setup_graph()
    validation_graph = create_validation_graph()

    # Add nodes
    builder.add_node("clarification_phase", clarification_graph)
    builder.add_node("prepare_setup_input", prepare_setup_input)
    builder.add_node("setup_phase", setup_graph)
    builder.add_node("prepare_validation_input", prepare_validation_input)
    builder.add_node("validation_phase", validation_graph)

    # Connect phases sequentially
    builder.add_edge(START, "clarification_phase")
    builder.add_edge("clarification_phase", "prepare_setup_input")
    builder.add_edge("prepare_setup_input", "setup_phase")
    builder.add_edge("setup_phase", "prepare_validation_input")
    builder.add_edge("prepare_validation_input", "validation_phase")
    builder.add_edge("validation_phase", END)

    # Setup checkpointing
    if checkpoint_path is None:
        checkpoint_path = Path("checkpoints/full_workflow.db")
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(str(checkpoint_path), check_same_thread=False)
    memory = SqliteSaver(conn)

    # Compile with optional interrupt for human-in-the-loop
    interrupt_nodes = ["clarification_phase"] if interrupt_after_clarification else []

    return builder.compile(checkpointer=memory, interrupt_after=interrupt_nodes)


def create_full_agent_no_interrupt(checkpoint_path: Optional[Path] = None):
    """Create full agent without interrupts (for automated testing).

    This is a convenience function that creates the full agent
    without pausing after the clarification phase.

    Args:
        checkpoint_path: Path to SQLite checkpoint file

    Returns:
        Compiled LangGraph without interrupts
    """
    return create_full_agent(checkpoint_path, interrupt_after_clarification=False)


# Pre-compiled agent for import convenience (with default settings)
full_agent = create_full_agent()
