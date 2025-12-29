"""ADK Tool configurations for MCP-MD workflow."""

try:
    from mcp_md_adk.tools.mcp_setup import create_mcp_toolsets
    __all__ = ["create_mcp_toolsets"]
except ImportError:
    # google-adk not installed
    __all__ = []
