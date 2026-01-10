"""
Unified MCP Server - Combines all MDZen MCP tools into a single server.

This server is designed for Docker deployment with Streamable HTTP transport.
It imports and re-exports all tools from individual servers under a unified namespace.

Usage:
    python unified_server.py                    # stdio transport (default)
    python unified_server.py --http --port 3000 # Streamable HTTP transport
"""

import argparse
import os
import sys
from pathlib import Path

from fastmcp import FastMCP

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from common.utils import setup_logger  # noqa: E402

logger = setup_logger(__name__)

# Create unified FastMCP server
mcp = FastMCP(
    "MDZen MCP Server",
    instructions="""
MDZen is a molecular dynamics simulation setup assistant.

Available tool categories:
- Research: download_structure, get_alphafold_structure, inspect_molecules, search_proteins, get_protein_info
- Structure: clean_protein, clean_ligand, prepare_complex, split_molecules, merge_structures
- Genesis: boltz2_protein_from_seq, rdkit_validate_smiles, pubchem_search_similar
- Solvation: solvate_structure, embed_in_membrane, list_available_lipids
- Amber: build_amber_system
- MD Simulation: run_md_simulation, analyze_rmsd, analyze_rmsf, calculate_distance, analyze_hydrogen_bonds

Workflow:
1. Download or prepare structure (research/structure tools)
2. Clean and prepare complex (structure tools)
3. Solvate with water or membrane (solvation tools)
4. Build Amber topology (amber tools)
5. Run and analyze MD simulation (md_simulation tools)
""",
)


def import_server_tools():
    """Import and register tools from all individual servers."""
    # Import individual servers (this registers their tools with their mcp instances)
    from servers import research_server
    from servers import structure_server
    from servers import genesis_server
    from servers import solvation_server
    from servers import amber_server
    from servers import md_simulation_server

    # Get all tool functions from each server's mcp instance
    servers = [
        ("research", research_server.mcp),
        ("structure", structure_server.mcp),
        ("genesis", genesis_server.mcp),
        ("solvation", solvation_server.mcp),
        ("amber", amber_server.mcp),
        ("md_simulation", md_simulation_server.mcp),
    ]

    registered_tools = []

    for name, server_mcp in servers:
        # In FastMCP 2.x, tools are stored in _tool_manager
        if hasattr(server_mcp, "_tool_manager") and server_mcp._tool_manager:
            tool_manager = server_mcp._tool_manager
            if hasattr(tool_manager, "_tools"):
                for tool_name, tool in tool_manager._tools.items():
                    # Register tool with unified server
                    if hasattr(tool, "fn"):
                        mcp.tool()(tool.fn)
                        registered_tools.append(f"{name}:{tool_name}")
                        logger.debug(f"Registered tool: {name}:{tool_name}")

    logger.info(f"Registered {len(registered_tools)} tools from {len(servers)} servers")
    return registered_tools


# Health check endpoint for Docker
@mcp.tool()
def health_check() -> dict:
    """Check if the MCP server is running and healthy.

    Returns:
        dict: Status information including server name and available tools count.
    """
    return {
        "status": "healthy",
        "server": "MDZen MCP Server",
        "version": "0.2.0",
    }


def main():
    """Run the unified MCP server."""
    parser = argparse.ArgumentParser(description="MDZen Unified MCP Server")
    parser.add_argument(
        "--http",
        action="store_true",
        help="Use Streamable HTTP transport instead of stdio",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=3000,
        help="Port for HTTP transport (default: 3000)",
    )
    parser.add_argument(
        "--host",
        type=str,
        default="0.0.0.0",
        help="Host for HTTP transport (default: 0.0.0.0)",
    )
    args = parser.parse_args()

    # Import and register tools from all servers
    try:
        tools = import_server_tools()
        logger.info(f"Successfully loaded {len(tools)} tools")
    except Exception as e:
        logger.error(f"Failed to import server tools: {e}")
        # Continue anyway - at least health_check will be available

    # Run server with appropriate transport
    if args.http:
        logger.info(f"Starting HTTP server on {args.host}:{args.port}")
        mcp.run(transport="http", host=args.host, port=args.port)
    else:
        logger.info("Starting stdio server")
        mcp.run()


if __name__ == "__main__":
    main()
