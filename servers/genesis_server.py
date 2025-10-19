"""
Genesis Server - Boltz-2 structure prediction from FASTA sequences.

Provides MCP tools for:
- Predicting protein structures from FASTA (de novo)
- Completing missing residues
- Generating structures without PDB data
"""

import logging
import json
from pathlib import Path
from typing import List, Dict, Any
from mcp.server import Server
from mcp.types import Tool, TextContent

from .base_server import BaseMCPServer
from tools.boltz2_wrapper import Boltz2Wrapper
from core.utils import setup_logger

logger = setup_logger(__name__)


class GenesisServer(BaseMCPServer):
    """MCP Server for Boltz-2 structure prediction"""
    
    def __init__(self):
        super().__init__("genesis_server", "0.1.0")
        self.boltz2 = Boltz2Wrapper()
        self.setup_handlers()
    
    def setup_handlers(self):
        """Register MCP tool handlers"""
        
        @self.server.list_tools()
        async def list_tools() -> list[Tool]:
            return [
                Tool(
                    name="boltz2_protein_from_seq",
                    description="Generate protein structure from FASTA sequence using Boltz-2",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "fasta": {"type": "string", "description": "Protein FASTA sequence"},
                            "use_msa": {"type": "boolean", "default": True, "description": "Use MSA for prediction"},
                            "num_models": {"type": "integer", "default": 5, "description": "Number of models to generate"}
                        },
                        "required": ["fasta"]
                    }
                ),
                Tool(
                    name="boltz2_complete_missing",
                    description="Complete missing residues in a structure using Boltz-2",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "pdb_file": {"type": "string", "description": "Input PDB file with missing residues"},
                            "fasta": {"type": "string", "description": "Full FASTA sequence"},
                            "use_msa": {"type": "boolean", "default": True}
                        },
                        "required": ["pdb_file", "fasta"]
                    }
                ),
            ]
        
        @self.server.call_tool()
        async def call_tool(name: str, arguments: dict) -> list[TextContent]:
            """Handle tool calls"""
            try:
                logger.info(f"Calling tool: {name}")
                
                if name == "boltz2_protein_from_seq":
                    result = await self.boltz2_protein_from_seq(
                        fasta=arguments["fasta"],
                        use_msa=arguments.get("use_msa", True),
                        num_models=arguments.get("num_models", 5)
                    )
                elif name == "boltz2_complete_missing":
                    result = await self.boltz2_complete_missing(
                        pdb_file=arguments["pdb_file"],
                        fasta=arguments["fasta"],
                        use_msa=arguments.get("use_msa", True)
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
    
    async def boltz2_protein_from_seq(
        self,
        fasta: str,
        use_msa: bool = True,
        num_models: int = 5
    ) -> dict:
        """Generate protein structure from FASTA sequence"""
        logger.info(f"Generating structure from FASTA (MSA={use_msa}, models={num_models})")
        
        output_dir = self.get_output_path("boltz2_genesis")
        
        sequences = [
            {"protein": {"id": "protein_A", "sequence": fasta}}
        ]
        
        result = self.boltz2.predict_structure(
            sequences=sequences,
            output_dir=output_dir,
            use_msa=use_msa,
            num_models=num_models
        )
        
        logger.info(f"Structure prediction completed: {result.get('output_pdb')}")
        return result
    
    async def boltz2_complete_missing(
        self,
        pdb_file: str,
        fasta: str,
        use_msa: bool = True
    ) -> dict:
        """Complete missing residues using Boltz-2"""
        logger.info(f"Completing missing residues in {pdb_file}")
        
        output_dir = self.get_output_path("boltz2_completion")
        
        # Use Boltz-2 to predict full structure with template
        sequences = [
            {
                "protein": {
                    "id": "protein_A",
                    "sequence": fasta
                }
            }
        ]
        
        result = self.boltz2.predict_structure(
            sequences=sequences,
            output_dir=output_dir,
            use_msa=use_msa,
            num_models=1,
            template_pdb=pdb_file
        )
        
        logger.info(f"Missing residue completion completed")
        return result


async def main():
    """Run genesis server"""
    server = GenesisServer()
    await server.run()


if __name__ == "__main__":
    import asyncio
    asyncio.run(main())

