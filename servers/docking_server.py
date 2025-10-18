"""
Docking Server - smina-based ligand docking.

Provides MCP tools for:
- Ligand docking with smina
- Pose clustering and analysis
"""

import logging
import json
from pathlib import Path
from typing import Tuple
from mcp.server import Server
from mcp.types import Tool, TextContent

from .base_server import BaseMCPServer
from tools.smina_wrapper import SminaWrapper
from core.utils import setup_logger

logger = setup_logger(__name__)


class DockingServer(BaseMCPServer):
    """MCP Server for docking operations"""
    
    def __init__(self):
        super().__init__("docking_server", "0.1.0")
        self.smina = SminaWrapper()
        self.setup_handlers()
    
    def setup_handlers(self):
        """Register MCP tool handlers"""
        
        @self.server.list_tools()
        async def list_tools() -> list[Tool]:
            return [
                Tool(
                    name="dock_ligand_smina",
                    description="Dock ligand using smina",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "receptor_pdb": {"type": "string"},
                            "ligand_pdb": {"type": "string"},
                            "center": {
                                "type": "array",
                                "items": {"type": "number"},
                                "minItems": 3,
                                "maxItems": 3
                            },
                            "size": {
                                "type": "array",
                                "items": {"type": "number"},
                                "minItems": 3,
                                "maxItems": 3
                            },
                            "scoring": {"type": "string", "default": "vinardo"},
                            "exhaustiveness": {"type": "integer", "default": 8}
                        },
                        "required": ["receptor_pdb", "ligand_pdb", "center", "size"]
                    }
                ),
            ]
        
        @self.server.call_tool()
        async def call_tool(name: str, arguments: dict) -> list[TextContent]:
            try:
                if name == "dock_ligand_smina":
                    result = await self.dock_ligand_smina(
                        receptor_pdb=arguments["receptor_pdb"],
                        ligand_pdb=arguments["ligand_pdb"],
                        center=tuple(arguments["center"]),
                        size=tuple(arguments["size"]),
                        scoring=arguments.get("scoring", "vinardo"),
                        exhaustiveness=arguments.get("exhaustiveness", 8)
                    )
                else:
                    raise ValueError(f"Unknown tool: {name}")
                
                return self.create_tool_response(json.dumps(result, indent=2))
            
            except Exception as e:
                logger.error(f"Tool {name} failed: {e}")
                return self.create_tool_response(
                    json.dumps({"error": str(e)}),
                    is_error=True
                )
    
    async def dock_ligand_smina(
        self,
        receptor_pdb: str,
        ligand_pdb: str,
        center: Tuple[float, float, float],
        size: Tuple[float, float, float],
        scoring: str = "vinardo",
        exhaustiveness: int = 8
    ) -> dict:
        """Dock ligand with smina"""
        logger.info(f"Docking {ligand_pdb} to {receptor_pdb}")
        
        output_dir = self.get_output_path("docking")
        
        result = self.smina.dock_ligand(
            receptor=receptor_pdb,
            ligand=ligand_pdb,
            center=center,
            size=size,
            output_dir=output_dir,
            scoring=scoring,
            exhaustiveness=exhaustiveness
        )
        
        return result


async def main():
    """Run docking server"""
    server = DockingServer()
    await server.run()


if __name__ == "__main__":
    import asyncio
    asyncio.run(main())

