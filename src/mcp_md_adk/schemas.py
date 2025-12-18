"""Pydantic schemas for MCP-MD ADK.

Re-exports schemas from mcp_md for compatibility and adds ADK-specific schemas.
"""

# Re-export all schemas from mcp_md
from mcp_md.state_scope import (
    SimulationBrief,
    ClarifyWithUser,
)

__all__ = [
    "SimulationBrief",
    "ClarifyWithUser",
]
