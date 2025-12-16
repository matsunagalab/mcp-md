
"""User Clarification and Simulation Brief Generation.

This module implements the clarification phase of MD setup using LangGraph 1.0+ patterns:
- Command API for node routing
- Structured outputs with Pydantic
- MessagesState-based state management
"""

from datetime import datetime
from typing import Literal

from langchain.chat_models import init_chat_model
from langchain_core.messages import AIMessage, HumanMessage, get_buffer_string
from langgraph.graph import END, START, StateGraph
from langgraph.types import Command

from mcp_md.prompts import (
    clarify_requirements_prompt,
    generate_simulation_brief_prompt,
)
from mcp_md.state_scope import (
    AgentInputState,
    AgentState,
    ClarifyWithUser,
    SimulationBrief,
)

def get_today_str() -> str:
    """Get current date formatted for prompts."""
    return datetime.now().strftime("%a %b %-d, %Y")

# Initialize model (LangGraph 1.0+ compatible)
# Using Anthropic Claude Haiku 4.5
model = init_chat_model(model="anthropic:claude-haiku-4-5-20251001", temperature=0.0)

# Alternative models (uncomment to use):
# Option 1: OpenAI
# model = init_chat_model(model="openai:gpt-4o", temperature=0.0)
# Option 2: Ollama (local)
# from langchain_ollama import ChatOllama
# model = ChatOllama(model="gemma2:9b", temperature=0.0)

def clarify_requirements(
    state: AgentState,
) -> Command[Literal["generate_simulation_brief", "__end__"]]:
    """Determine if sufficient information exists to proceed with MD setup.

    Returns:
        Command: Routes to END if clarification needed, otherwise to generate_simulation_brief
    """
    structured_model = model.with_structured_output(ClarifyWithUser)
    response = structured_model.invoke(
        [
            HumanMessage(
                content=clarify_requirements_prompt.format(
                    messages=get_buffer_string(messages=state["messages"]),
                    date=get_today_str(),
                )
            )
        ]
    )

    if response.need_clarification:
        # Need more info - return to user
        return Command(
            goto=END, update={"messages": [AIMessage(content=response.question)]}
        )
    else:
        # Have enough info - proceed to brief generation
        return Command(
            goto="generate_simulation_brief",
            update={"messages": [AIMessage(content=response.verification)]},
        )

def generate_simulation_brief(state: AgentState) -> dict:
    """Generate structured simulation brief from conversation history.

    Returns:
        dict: Updated state with simulation_brief, research_brief, and setup_messages
    """
    structured_model = model.with_structured_output(SimulationBrief)
    response = structured_model.invoke(
        [
            HumanMessage(
                content=generate_simulation_brief_prompt.format(
                    messages=get_buffer_string(state.get("messages", [])),
                    date=get_today_str(),
                )
            )
        ]
    )

    return {
        "simulation_brief": response,
        "research_brief": str(response.model_dump()),
        "setup_messages": [
            HumanMessage(content=f"Starting MD setup with: {response.model_dump_json()}")
        ],
    }

# Build clarification graph (LangGraph 1.0+ pattern)
clarification_builder = StateGraph(AgentState, input_schema=AgentInputState)
clarification_builder.add_node("clarify_requirements", clarify_requirements)
clarification_builder.add_node("generate_simulation_brief", generate_simulation_brief)
clarification_builder.add_edge(START, "clarify_requirements")
clarification_builder.add_edge("generate_simulation_brief", END)
clarification_graph = clarification_builder.compile()

def create_clarification_graph():
    """Create and return the clarification graph.

    Factory function for creating the clarification phase graph.
    Useful for integration into larger workflows.

    Returns:
        Compiled StateGraph for clarification phase
    """
    builder = StateGraph(AgentState, input_schema=AgentInputState)
    builder.add_node("clarify_requirements", clarify_requirements)
    builder.add_node("generate_simulation_brief", generate_simulation_brief)
    builder.add_edge(START, "clarify_requirements")
    builder.add_edge("generate_simulation_brief", END)
    return builder.compile()
