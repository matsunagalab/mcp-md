
"""State Definitions and Pydantic Schemas for Clarification Phase.

This module defines state objects and structured schemas for the MD setup workflow
using LangGraph 1.0+ patterns:
- MessagesState for conversation tracking
- Annotated fields with proper reducers
- Pydantic models for structured outputs
"""

import operator
from typing import Annotated, Optional, Sequence

from langchain_core.messages import BaseMessage
from langgraph.graph import MessagesState
from langgraph.graph.message import add_messages
from pydantic import BaseModel, Field


# ===== STATE DEFINITIONS (LangGraph 1.0+) =====

class AgentInputState(MessagesState):
    """Input state for the full agent - only contains messages from user input.

    Used as input_schema for the main graph to define the public input interface.
    """

    pass


class AgentState(MessagesState):
    """Main state for the full multi-phase MD setup system.

    Extends MessagesState with additional fields for MD setup coordination.
    All fields use proper reducers for state updates.
    """

    # Phase 1: Clarification
    research_brief: Optional[str] = None  # Compatibility with deep_research pattern
    simulation_brief: Optional["SimulationBrief"] = None

    # Phase 2: Setup (will be added in Notebook 2)
    setup_messages: Annotated[Sequence[BaseMessage], add_messages] = []
    decision_log: Annotated[list[dict], operator.add] = []
    outputs: dict = {}

    # Phase 3: Validation & Export (will be added in Notebook 4)
    qc_results: dict = {}
    exports: dict = {}
    final_report: str = ""


# ===== STRUCTURED OUTPUT SCHEMAS (Pydantic) =====

class ClarifyWithUser(BaseModel):
    """Schema for user clarification decision and questions.

    Used with .with_structured_output() for LLM decision making.
    """

    need_clarification: bool = Field(
        description="Whether the user needs to be asked a clarifying question.",
    )
    question: str = Field(
        description="A question to ask the user to clarify the simulation requirements.",
    )
    verification: str = Field(
        description="Verification message that we will start setup after user provides information.",
    )


class SimulationBrief(BaseModel):
    """Schema for structured simulation brief generation.

    Transforms conversation into structured MD setup parameters.
    """

    # Structure
    pdb_id: Optional[str] = Field(default=None, description="PDB ID (e.g., 1ABC)")
    fasta_sequence: Optional[str] = Field(
        default=None, description="FASTA sequence for de novo generation"
    )
    ligand_smiles: Optional[str] = Field(default=None, description="Ligand SMILES string")

    # Simulation parameters
    ph: float = Field(default=7.0, description="pH value")
    salt_concentration: float = Field(default=0.15, description="Salt concentration (M)")
    water_model: str = Field(default="TIP3P", description="Water model")
    box_padding: float = Field(default=12.0, description="Box padding (Å)")
    force_field: str = Field(default="ff19SB", description="Protein force field")

    # Workflow preferences
    use_boltz2_docking: bool = Field(default=True, description="Use Boltz-2 for docking")
    refine_with_smina: bool = Field(default=False, description="Refine with Smina")
    output_formats: list[str] = Field(default=["amber"], description="Output formats")
