"""
MCP Integration - LangChain MCP Adapters for tool loading.

Provides MCP client setup and tool loading functionality using
the official langchain-mcp-adapters package.
"""

import sys
import logging
from typing import Dict
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain_core.tools import BaseTool

logger = logging.getLogger(__name__)


def create_mcp_client() -> MultiServerMCPClient:
    """Create MCP client with all FastMCP servers configured
    
    Returns:
        MultiServerMCPClient instance configured for all servers
    """
    # Get Python executable path for subprocess execution
    python_exe = sys.executable
    
    # Configure all MCP servers with stdio transport
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
        "complex": {
            "transport": "stdio",
            "command": python_exe,
            "args": ["-m", "servers.complex_server"]
        },
        "ligand": {
            "transport": "stdio",
            "command": python_exe,
            "args": ["-m", "servers.ligand_server"]
        },
        "assembly": {
            "transport": "stdio",
            "command": python_exe,
            "args": ["-m", "servers.assembly_server"]
        },
        "export": {
            "transport": "stdio",
            "command": python_exe,
            "args": ["-m", "servers.export_server"]
        },
        "qc_min": {
            "transport": "stdio",
            "command": python_exe,
            "args": ["-m", "servers.qc_min_server"]
        }
    }
    
    logger.info(f"Creating MCP client with {len(server_config)} servers")
    return MultiServerMCPClient(server_config)


async def load_all_mcp_tools() -> Dict[str, BaseTool]:
    """Load all MCP tools from configured servers
    
    Returns:
        Dictionary mapping tool names to LangChain Tool objects
        
    Note:
        MultiServerMCPClient is stateless by default. Each tool invocation
        creates a fresh ClientSession, executes the tool, and cleans up.
    """
    client = create_mcp_client()
    
    try:
        logger.info("Loading tools from all MCP servers...")
        tools = await client.get_tools()
        
        # Convert list to dict for easy lookup by name
        tools_dict = {tool.name: tool for tool in tools}
        
        logger.info(f"Loaded {len(tools_dict)} tools from MCP servers")
        logger.debug(f"Available tools: {list(tools_dict.keys())}")
        
        return tools_dict
    
    except Exception as e:
        logger.error(f"Failed to load MCP tools: {e}")
        raise


async def load_mcp_tools_stateful(server_name: str):
    """Load tools from a specific server with stateful session
    
    Args:
        server_name: Name of the MCP server (e.g., "structure", "genesis")
    
    Returns:
        List of tools from the specified server
        
    Note:
        Use this for servers that maintain state between tool calls.
        The session must be used within an async context manager.
    
    Example:
        ```python
        from langchain_mcp_adapters.tools import load_mcp_tools
        
        client = create_mcp_client()
        async with client.session(server_name) as session:
            tools = await load_mcp_tools(session)
            # Use tools within this context
        ```
    """
    from langchain_mcp_adapters.tools import load_mcp_tools
    
    client = create_mcp_client()
    
    logger.info(f"Creating stateful session for server: {server_name}")
    async with client.session(server_name) as session:
        tools = await load_mcp_tools(session)
        logger.info(f"Loaded {len(tools)} tools from {server_name}")
        return tools

