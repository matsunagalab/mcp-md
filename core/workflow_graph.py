"""
Workflow Graph - LangGraph StateGraph construction and compilation.

Defines the fixed workflow skeleton with conditional edges for error handling.
"""

import logging
from typing import Dict
from functools import partial

from langgraph.graph import StateGraph, END
from langgraph.checkpoint.sqlite import SqliteSaver
from langchain_core.tools import BaseTool

from .workflow_state import WorkflowState
from .workflow_nodes import (
    planner_node,
    structure_fetch_node,
    repair_node,
    ligand_param_node,
    complex_node,
    assemble_node,
    qc_node,
    should_retry
)
from .mcp_integration import load_all_mcp_tools

logger = logging.getLogger(__name__)


async def create_workflow_graph(checkpoint_path: str = "checkpoints/workflow.db"):
    """Build and compile the MD workflow graph
    
    Args:
        checkpoint_path: Path to SQLite database for checkpointing
    
    Returns:
        Compiled LangGraph application with checkpointing enabled
    """
    logger.info("Creating workflow graph...")
    
    # Load all MCP tools
    mcp_tools = await load_all_mcp_tools()
    logger.info(f"Loaded {len(mcp_tools)} MCP tools")
    
    # Create state graph
    graph = StateGraph(WorkflowState)
    
    # Add nodes (using partial to bind tools to async functions)
    graph.add_node("planner", planner_node)
    graph.add_node("fetch", partial(structure_fetch_node, tools=mcp_tools))
    graph.add_node("repair", partial(repair_node, tools=mcp_tools))
    graph.add_node("ligand_param", partial(ligand_param_node, tools=mcp_tools))
    graph.add_node("complex", partial(complex_node, tools=mcp_tools))
    graph.add_node("assemble", partial(assemble_node, tools=mcp_tools))
    graph.add_node("qc", partial(qc_node, tools=mcp_tools))
    
    # Define edges (fixed skeleton)
    graph.set_entry_point("planner")
    graph.add_edge("planner", "fetch")
    graph.add_edge("fetch", "repair")
    graph.add_edge("repair", "ligand_param")
    graph.add_edge("ligand_param", "complex")
    
    # Conditional edges for error handling on complex step
    graph.add_conditional_edges(
        "complex",
        should_retry,
        {
            "retry": "complex",  # Retry same step
            "human_feedback": END,  # Stop and wait for human
            "continue": "assemble"  # Proceed to next step
        }
    )
    
    graph.add_edge("assemble", "qc")
    
    # Conditional edges for error handling on QC step
    graph.add_conditional_edges(
        "qc",
        should_retry,
        {
            "retry": "qc",
            "human_feedback": END,
            "continue": END  # Success!
        }
    )
    
    # Setup checkpointer for persistence
    memory = SqliteSaver.from_conn_string(checkpoint_path)
    logger.info(f"Checkpointer initialized: {checkpoint_path}")
    
    # Compile graph with checkpointing
    app = graph.compile(checkpointer=memory)
    
    logger.info("Workflow graph compiled successfully")
    return app


async def create_workflow_graph_with_interrupts(
    checkpoint_path: str = "checkpoints/workflow.db",
    interrupt_before: list[str] | None = None,
    interrupt_after: list[str] | None = None
):
    """Build workflow graph with human-in-the-loop interrupts
    
    Args:
        checkpoint_path: Path to SQLite database for checkpointing
        interrupt_before: List of node names to interrupt before
        interrupt_after: List of node names to interrupt after
    
    Returns:
        Compiled LangGraph application with interrupts enabled
        
    Example:
        ```python
        # Interrupt before complex generation for human approval
        app = await create_workflow_graph_with_interrupts(
            interrupt_before=["complex"]
        )
        ```
    """
    logger.info("Creating workflow graph with interrupts...")
    
    # Load all MCP tools
    mcp_tools = await load_all_mcp_tools()
    
    # Create state graph
    graph = StateGraph(WorkflowState)
    
    # Add nodes
    graph.add_node("planner", planner_node)
    graph.add_node("fetch", partial(structure_fetch_node, tools=mcp_tools))
    graph.add_node("repair", partial(repair_node, tools=mcp_tools))
    graph.add_node("ligand_param", partial(ligand_param_node, tools=mcp_tools))
    graph.add_node("complex", partial(complex_node, tools=mcp_tools))
    graph.add_node("assemble", partial(assemble_node, tools=mcp_tools))
    graph.add_node("qc", partial(qc_node, tools=mcp_tools))
    
    # Define edges
    graph.set_entry_point("planner")
    graph.add_edge("planner", "fetch")
    graph.add_edge("fetch", "repair")
    graph.add_edge("repair", "ligand_param")
    graph.add_edge("ligand_param", "complex")
    
    graph.add_conditional_edges(
        "complex",
        should_retry,
        {
            "retry": "complex",
            "human_feedback": END,
            "continue": "assemble"
        }
    )
    
    graph.add_edge("assemble", "qc")
    
    graph.add_conditional_edges(
        "qc",
        should_retry,
        {
            "retry": "qc",
            "human_feedback": END,
            "continue": END
        }
    )
    
    # Setup checkpointer
    memory = SqliteSaver.from_conn_string(checkpoint_path)
    
    # Compile with interrupts
    app = graph.compile(
        checkpointer=memory,
        interrupt_before=interrupt_before,
        interrupt_after=interrupt_after
    )
    
    logger.info(f"Workflow graph compiled with interrupts: "
                f"before={interrupt_before}, after={interrupt_after}")
    return app


def get_workflow_steps() -> list[str]:
    """Get list of all workflow step names
    
    Returns:
        List of step names in execution order
    """
    return [
        "planner",
        "fetch",
        "repair",
        "ligand_param",
        "complex",
        "assemble",
        "qc"
    ]

