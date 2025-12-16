
"""MCP Integration for MD Setup Tools."""

import sys
from langchain_mcp_adapters.client import MultiServerMCPClient

def create_mcp_client() -> MultiServerMCPClient:
    """Create MCP client with all 5 active servers configured."""
    python_exe = sys.executable

    server_config = {
        "structure": {
            "transport": "stdio",
            "command": python_exe,
            "args": ["-m", "servers.structure_server"]
        },
        "genesis": {
            "transport": "stdio",
            "command": python_exe,
            "args": ["-m", "servers.genesis_server"]
        },
        "solvation": {
            "transport": "stdio",
            "command": python_exe,
            "args": ["-m", "servers.solvation_server"]
        },
        "amber": {
            "transport": "stdio",
            "command": python_exe,
            "args": ["-m", "servers.amber_server"]
        },
        "md_simulation": {
            "transport": "stdio",
            "command": python_exe,
            "args": ["-m", "servers.md_simulation_server"]
        }
    }

    return MultiServerMCPClient(server_config)

async def load_mcp_tools():
    """Load all MCP tools as a dictionary."""
    client = create_mcp_client()
    tools = await client.get_tools()
    return {tool.name: tool for tool in tools}
