"""
Assembly Server - tleap system building.

Provides MCP tools for:
- System solvation
- Ion addition
- Complete system assembly
"""

import logging
import json
from mcp.types import Tool, TextContent

from .base_server import BaseMCPServer
from tools.ambertools_wrapper import AmberToolsWrapper
from tools.packmol_wrapper import PackmolMemgenWrapper
from core.utils import setup_logger

logger = setup_logger(__name__)


class AssemblyServer(BaseMCPServer):
    """MCP Server for system assembly"""
    
    def __init__(self):
        super().__init__("assembly_server", "0.1.0")
        self.ambertools = AmberToolsWrapper()
        self.packmol_memgen = PackmolMemgenWrapper()
        self.setup_handlers()
    
    def setup_handlers(self):
        """Register MCP tool handlers"""
        
        @self.server.list_tools()
        async def list_tools() -> list[Tool]:
            return [
                Tool(
                    name="build_system_tleap",
                    description="Build complete MD system with tleap (solvation + ions)",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "protein_pdb": {"type": "string"},
                            "ligand_lib": {"type": "string"},
                            "forcefield": {"type": "string", "default": "leaprc.protein.ff19SB"},
                            "water_model": {"type": "string", "default": "tip3p"},
                            "box_padding": {"type": "number", "default": 10.0},
                            "salt_conc": {"type": "number", "default": 0.15}
                        }
                    }
                ),
                Tool(
                    name="build_membrane_system",
                    description="Build membrane protein system with Packmol-Memgen",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "protein_pdb": {"type": "string"},
                            "lipid_composition": {
                                "type": "object",
                                "description": "Lipid composition dict, e.g. {'POPC': 0.7, 'POPE': 0.3}"
                            },
                            "membrane_type": {"type": "string", "default": "bilayer"},
                            "dist_to_bilayer": {"type": "number", "default": 15.0}
                        },
                        "required": ["protein_pdb"]
                    }
                ),
            ]
        
        @self.server.call_tool()
        async def call_tool(name: str, arguments: dict) -> list[TextContent]:
            try:
                if name == "build_system_tleap":
                    result = await self.build_system_tleap(
                        protein_pdb=arguments.get("protein_pdb"),
                        ligand_lib=arguments.get("ligand_lib"),
                        forcefield=arguments.get("forcefield", "leaprc.protein.ff19SB"),
                        water_model=arguments.get("water_model", "tip3p"),
                        box_padding=arguments.get("box_padding", 10.0),
                        salt_conc=arguments.get("salt_conc", 0.15)
                    )
                elif name == "build_membrane_system":
                    result = await self.build_membrane_system(
                        protein_pdb=arguments["protein_pdb"],
                        lipid_composition=arguments.get("lipid_composition", {"POPC": 1.0}),
                        membrane_type=arguments.get("membrane_type", "bilayer"),
                        dist_to_bilayer=arguments.get("dist_to_bilayer", 15.0)
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
    
    async def build_system_tleap(
        self,
        protein_pdb: str = None,
        ligand_lib: str = None,
        forcefield: str = "leaprc.protein.ff19SB",
        water_model: str = "tip3p",
        box_padding: float = 10.0,
        salt_conc: float = 0.15
    ) -> dict:
        """Build system with tleap"""
        logger.info("Building system with tleap")
        
        output_dir = self.get_output_path("system")
        
        result = self.ambertools.build_system_tleap(
            protein_pdb=protein_pdb,
            ligand_lib=ligand_lib,
            output_dir=output_dir,
            forcefield=forcefield,
            water_model=water_model,
            box_padding=box_padding,
            salt_conc=salt_conc
        )
        
        return result
    
    async def build_membrane_system(
        self,
        protein_pdb: str,
        lipid_composition: dict = None,
        membrane_type: str = "bilayer",
        dist_to_bilayer: float = 15.0
    ) -> dict:
        """Build membrane protein system"""
        logger.info(f"Building membrane system: {membrane_type}")
        
        if lipid_composition is None:
            lipid_composition = {"POPC": 1.0}
        
        output_dir = self.get_output_path("membrane_system")
        
        result = self.packmol_memgen.build_membrane_system(
            protein_pdb=protein_pdb,
            output_dir=output_dir,
            lipid_composition=lipid_composition,
            membrane_type=membrane_type,
            dist_to_bilayer=dist_to_bilayer
        )
        
        return result


async def main():
    """Run assembly server"""
    server = AssemblyServer()
    await server.run()


if __name__ == "__main__":
    import asyncio
    asyncio.run(main())

