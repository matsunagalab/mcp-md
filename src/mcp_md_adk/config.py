"""Configuration settings for MCP-MD ADK.

This module centralizes all configuration settings for the ADK implementation,
extending the base mcp_md.config with ADK-specific settings.

Usage:
    from mcp_md_adk.config import settings, get_litellm_model

    # Access settings
    model = get_litellm_model("clarification")
"""


# Re-export base config settings
from mcp_md.config import settings, get_output_dir, get_server_path


def get_litellm_model(model_type: str) -> str:
    """Convert mcp_md model string to LiteLLM format.

    mcp_md uses "anthropic:claude-xxx" format
    LiteLLM uses "anthropic/claude-xxx" format

    Args:
        model_type: Type of model ("clarification", "setup", or "compress")

    Returns:
        LiteLLM-compatible model string
    """
    model_map = {
        "clarification": settings.clarification_model,
        "setup": settings.setup_model,
        "compress": settings.compress_model,
    }

    model_str = model_map.get(model_type, settings.setup_model)

    # Convert "anthropic:model" to "anthropic/model" for LiteLLM
    if ":" in model_str:
        provider, model_name = model_str.split(":", 1)
        return f"{provider}/{model_name}"

    return model_str


def get_adk_config() -> dict:
    """Get ADK-specific configuration as a dictionary.

    Returns:
        Dictionary with ADK configuration settings
    """
    return {
        "clarification_model": get_litellm_model("clarification"),
        "setup_model": get_litellm_model("setup"),
        "compress_model": get_litellm_model("compress"),
        "output_dir": str(get_output_dir()),
        "max_message_history": settings.max_message_history,
        "default_timeout": settings.default_timeout,
        "md_simulation_timeout": settings.md_simulation_timeout,
        "servers": {
            "structure": get_server_path("structure"),
            "genesis": get_server_path("genesis"),
            "solvation": get_server_path("solvation"),
            "amber": get_server_path("amber"),
            "md_simulation": get_server_path("md_simulation"),
        },
    }
