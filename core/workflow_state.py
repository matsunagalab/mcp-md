"""
Workflow State - TypedDict for LangGraph state management.

Defines the state structure for the MD workflow graph.
"""

from typing import TypedDict, Annotated, Sequence


class WorkflowState(TypedDict, total=False):
    """State for MD workflow graph
    
    Attributes:
        query: User's natural language query
        pdb_id: PDB ID for structure fetching
        ligand_smiles: SMILES string for ligand
        current_step: Current workflow step name
        outputs: Dictionary of outputs from each step
        user_preferences: User configuration (pH, salt, forcefield, etc.)
        decision_log: Log of all decisions made (append-only)
        retry_count: Number of retries for current step
        error: Error message if step failed
        thread_id: Session identifier for checkpointing
    """
    
    # Input parameters
    query: str
    pdb_id: str | None
    ligand_smiles: str | None
    
    # Execution state
    current_step: str
    outputs: dict
    
    # User configuration
    user_preferences: dict  # {ph: 7.4, salt: 0.15, water_model: "TIP3P", ...}
    
    # Decision logging (append-only with Annotated)
    decision_log: Annotated[Sequence[dict], "append"]
    
    # Error handling
    retry_count: int
    error: str | None
    
    # Session management
    thread_id: str


def create_initial_state(
    query: str,
    pdb_id: str | None = None,
    ligand_smiles: str | None = None,
    user_preferences: dict | None = None,
    thread_id: str = "default"
) -> WorkflowState:
    """Create initial workflow state
    
    Args:
        query: User's natural language query
        pdb_id: Optional PDB ID
        ligand_smiles: Optional SMILES string
        user_preferences: Optional user configuration
        thread_id: Session identifier
    
    Returns:
        Initial workflow state dictionary
    """
    return WorkflowState(
        query=query,
        pdb_id=pdb_id,
        ligand_smiles=ligand_smiles,
        current_step="planner",
        outputs={},
        user_preferences=user_preferences or {
            "ph": 7.4,
            "salt_concentration": 0.15,
            "water_model": "TIP3P",
            "force_field": "ff19SB",
            "box_padding": 12.0,
        },
        decision_log=[],
        retry_count=0,
        error=None,
        thread_id=thread_id,
    )

