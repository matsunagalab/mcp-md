"""Configuration settings for MCP-MD.

This module centralizes all configuration settings, loading from
environment variables with sensible defaults.

Usage:
    from mcp_md.config import settings

    # Access settings
    model = settings.clarification_model
    output_dir = settings.output_dir
"""

import os
from pathlib import Path

# Try to use pydantic-settings if available, otherwise use dataclass
try:
    from pydantic_settings import BaseSettings

    class MCPMDSettings(BaseSettings):
        """MCP-MD configuration settings.

        Settings are loaded from environment variables with MCPMD_ prefix.
        For example: MCPMD_OUTPUT_DIR=/custom/output
        """

        # Output directory
        output_dir: str = "output"

        # LLM Models
        clarification_model: str = "anthropic:claude-haiku-4-5-20251001"
        setup_model: str = "anthropic:claude-sonnet-4-20250514"
        compress_model: str = "anthropic:claude-haiku-4-5-20251001"

        # Timeouts (seconds)
        default_timeout: int = 300
        md_simulation_timeout: int = 3600  # 1 hour for long simulations
        boltz_timeout: int = 1800  # 30 min for structure prediction

        # Token limits
        max_message_history: int = 6  # Number of messages to keep in context

        # MCP Server settings
        structure_server_path: str = "servers/structure_server.py"
        solvation_server_path: str = "servers/solvation_server.py"
        amber_server_path: str = "servers/amber_server.py"
        genesis_server_path: str = "servers/genesis_server.py"
        md_simulation_server_path: str = "servers/md_simulation_server.py"

        # Default simulation parameters
        default_temperature: float = 300.0  # Kelvin
        default_pressure: float = 1.0  # bar
        default_timestep: float = 2.0  # fs
        default_box_padding: float = 12.0  # Angstrom

        class Config:
            env_prefix = "MCPMD_"
            env_file = ".env"
            extra = "ignore"

except ImportError:
    # Fallback to simple dataclass-like class
    class MCPMDSettings:
        """MCP-MD configuration settings (fallback without pydantic-settings)."""

        def __init__(self):
            # Output directory
            self.output_dir = os.getenv("MCPMD_OUTPUT_DIR", "output")

            # LLM Models
            self.clarification_model = os.getenv(
                "MCPMD_CLARIFICATION_MODEL", "anthropic:claude-haiku-4-5-20251001"
            )
            self.setup_model = os.getenv(
                "MCPMD_SETUP_MODEL", "anthropic:claude-sonnet-4-20250514"
            )
            self.compress_model = os.getenv(
                "MCPMD_COMPRESS_MODEL", "anthropic:claude-haiku-4-5-20251001"
            )

            # Timeouts (seconds)
            self.default_timeout = int(os.getenv("MCPMD_DEFAULT_TIMEOUT", "300"))
            self.md_simulation_timeout = int(os.getenv("MCPMD_MD_SIMULATION_TIMEOUT", "3600"))
            self.boltz_timeout = int(os.getenv("MCPMD_BOLTZ_TIMEOUT", "1800"))

            # Token limits
            self.max_message_history = int(os.getenv("MCPMD_MAX_MESSAGE_HISTORY", "6"))

            # MCP Server paths
            self.structure_server_path = os.getenv(
                "MCPMD_STRUCTURE_SERVER_PATH", "servers/structure_server.py"
            )
            self.solvation_server_path = os.getenv(
                "MCPMD_SOLVATION_SERVER_PATH", "servers/solvation_server.py"
            )
            self.amber_server_path = os.getenv(
                "MCPMD_AMBER_SERVER_PATH", "servers/amber_server.py"
            )
            self.genesis_server_path = os.getenv(
                "MCPMD_GENESIS_SERVER_PATH", "servers/genesis_server.py"
            )
            self.md_simulation_server_path = os.getenv(
                "MCPMD_MD_SIMULATION_SERVER_PATH", "servers/md_simulation_server.py"
            )

            # Default simulation parameters
            self.default_temperature = float(os.getenv("MCPMD_DEFAULT_TEMPERATURE", "300.0"))
            self.default_pressure = float(os.getenv("MCPMD_DEFAULT_PRESSURE", "1.0"))
            self.default_timestep = float(os.getenv("MCPMD_DEFAULT_TIMESTEP", "2.0"))
            self.default_box_padding = float(os.getenv("MCPMD_DEFAULT_BOX_PADDING", "12.0"))


# Singleton instance
settings = MCPMDSettings()


def get_output_dir() -> Path:
    """Get resolved output directory path.

    Returns:
        Absolute Path to output directory
    """
    return Path(settings.output_dir).resolve()


def get_server_path(server_name: str) -> str:
    """Get path to an MCP server script.

    Args:
        server_name: Name of the server (structure, solvation, amber, genesis, md_simulation)

    Returns:
        Path to the server script

    Raises:
        ValueError: If server_name is not recognized
    """
    server_map = {
        "structure": settings.structure_server_path,
        "solvation": settings.solvation_server_path,
        "amber": settings.amber_server_path,
        "genesis": settings.genesis_server_path,
        "md_simulation": settings.md_simulation_server_path,
    }
    if server_name not in server_map:
        raise ValueError(f"Unknown server: {server_name}. Available: {list(server_map.keys())}")
    return server_map[server_name]
