"""
Complex Server - Boltz-2 complex prediction + Smina refinement.

Provides MCP tools for:
- Protein-ligand complex prediction with Boltz-2
- Binding affinity prediction
- Virtual screening
- Local refinement with Smina
"""

import logging
import json
from pathlib import Path
from typing import List, Tuple, Dict, Any
from mcp.server import Server
from mcp.types import Tool, TextContent

from .base_server import BaseMCPServer
from tools.boltz2_wrapper import Boltz2Wrapper
from tools.smina_wrapper import SminaWrapper
from core.utils import setup_logger

logger = setup_logger(__name__)


class ComplexServer(BaseMCPServer):
    """MCP Server for protein-ligand complex prediction"""
    
    def __init__(self):
        super().__init__("complex_server", "0.1.0")
        self.boltz2 = Boltz2Wrapper()
        self.smina = SminaWrapper()
        self.setup_handlers()
    
    def setup_handlers(self):
        """Register MCP tool handlers"""
        
        @self.server.list_tools()
        async def list_tools() -> list[Tool]:
            return [
                Tool(
                    name="boltz2_complex",
                    description="Predict protein-ligand complex with binding affinity using Boltz-2",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "protein_pdb": {"type": "string", "description": "Protein PDB file or FASTA sequence"},
                            "ligand_smiles": {
                                "type": "array",
                                "items": {"type": "string"},
                                "description": "List of ligand SMILES"
                            },
                            "use_msa": {"type": "boolean", "default": True},
                            "num_models": {"type": "integer", "default": 5},
                            "top_k": {"type": "integer", "default": 10, "description": "Number of top poses to return"}
                        },
                        "required": ["protein_pdb", "ligand_smiles"]
                    }
                ),
                Tool(
                    name="screen_ligands",
                    description="Screen multiple ligands for binding affinity (Boltz-2)",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "protein_pdb": {"type": "string"},
                            "ligand_smiles_list": {
                                "type": "array",
                                "items": {"type": "string"},
                                "description": "List of ligand SMILES to screen"
                            },
                            "screening_mode": {
                                "type": "string",
                                "enum": ["binary", "value"],
                                "default": "binary",
                                "description": "binary: hit/no-hit, value: IC50 prediction"
                            }
                        },
                        "required": ["protein_pdb", "ligand_smiles_list"]
                    }
                ),
                Tool(
                    name="smina_refine",
                    description="Refine Boltz-2 complex poses with Smina local search",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "receptor_pdb": {"type": "string"},
                            "ligand_pdb": {"type": "string"},
                            "center": {
                                "type": "array",
                                "items": {"type": "number"},
                                "minItems": 3,
                                "maxItems": 3,
                                "description": "Binding site center [x, y, z]"
                            },
                            "size": {
                                "type": "array",
                                "items": {"type": "number"},
                                "minItems": 3,
                                "maxItems": 3,
                                "description": "Search box size [x, y, z]"
                            },
                            "scoring": {"type": "string", "default": "vinardo"},
                            "exhaustiveness": {"type": "integer", "default": 16, "description": "Search exhaustiveness"}
                        },
                        "required": ["receptor_pdb", "ligand_pdb", "center", "size"]
                    }
                ),
            ]
        
        @self.server.call_tool()
        async def call_tool(name: str, arguments: dict) -> list[TextContent]:
            """Handle tool calls"""
            try:
                logger.info(f"Calling tool: {name}")
                
                if name == "boltz2_complex":
                    result = await self.boltz2_complex(
                        protein_pdb=arguments["protein_pdb"],
                        ligand_smiles=arguments["ligand_smiles"],
                        use_msa=arguments.get("use_msa", True),
                        num_models=arguments.get("num_models", 5),
                        top_k=arguments.get("top_k", 10)
                    )
                elif name == "screen_ligands":
                    result = await self.screen_ligands(
                        protein_pdb=arguments["protein_pdb"],
                        ligand_smiles_list=arguments["ligand_smiles_list"],
                        screening_mode=arguments.get("screening_mode", "binary")
                    )
                elif name == "smina_refine":
                    result = await self.smina_refine(
                        receptor_pdb=arguments["receptor_pdb"],
                        ligand_pdb=arguments["ligand_pdb"],
                        center=tuple(arguments["center"]),
                        size=tuple(arguments["size"]),
                        scoring=arguments.get("scoring", "vinardo"),
                        exhaustiveness=arguments.get("exhaustiveness", 16)
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
    
    async def boltz2_complex(
        self,
        protein_pdb: str,
        ligand_smiles: List[str],
        use_msa: bool = True,
        num_models: int = 5,
        top_k: int = 10
    ) -> dict:
        """Predict protein-ligand complex with Boltz-2"""
        logger.info(f"Predicting complex with {len(ligand_smiles)} ligands")
        
        output_dir = self.get_output_path("boltz2_complex")
        
        # Check if protein_pdb is a file or FASTA sequence
        if Path(protein_pdb).exists():
            # Load FASTA from PDB
            from Bio import SeqIO
            from Bio.PDB import PDBParser
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("protein", protein_pdb)
            # Extract sequence (simplified - assumes single chain)
            protein_fasta = str(list(structure.get_residues())[0].get_parent().get_parent().get_sequence())
        else:
            # Assume it's a FASTA sequence
            protein_fasta = protein_pdb
        
        result = self.boltz2.predict_complex_with_affinity(
            protein_fasta=protein_fasta,
            ligand_smiles=ligand_smiles,
            output_dir=output_dir,
            use_msa=use_msa,
            num_models=num_models
        )
        
        # Sort by affinity and return top_k
        if "poses" in result:
            sorted_poses = sorted(result["poses"], key=lambda x: x.get("affinity_score", 0), reverse=True)
            result["top_poses"] = sorted_poses[:top_k]
        
        logger.info(f"Complex prediction completed: {len(result.get('poses', []))} poses")
        return result
    
    async def screen_ligands(
        self,
        protein_pdb: str,
        ligand_smiles_list: List[str],
        screening_mode: str = "binary"
    ) -> dict:
        """Screen multiple ligands for binding"""
        logger.info(f"Screening {len(ligand_smiles_list)} ligands (mode={screening_mode})")
        
        output_dir = self.get_output_path("boltz2_screening")
        
        # Extract protein FASTA
        if Path(protein_pdb).exists():
            from Bio import SeqIO
            from Bio.PDB import PDBParser
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("protein", protein_pdb)
            protein_fasta = str(list(structure.get_residues())[0].get_parent().get_parent().get_sequence())
        else:
            protein_fasta = protein_pdb
        
        result = self.boltz2.screen_ligands(
            protein_fasta=protein_fasta,
            ligand_smiles_list=ligand_smiles_list,
            output_dir=output_dir,
            screening_mode=screening_mode
        )
        
        # Rank ligands by affinity
        if "results" in result:
            ranked = sorted(result["results"], key=lambda x: x.get("score", 0), reverse=True)
            result["ranked_ligands"] = ranked
        
        logger.info(f"Screening completed: {len(result.get('results', []))} ligands")
        return result
    
    async def smina_refine(
        self,
        receptor_pdb: str,
        ligand_pdb: str,
        center: Tuple[float, float, float],
        size: Tuple[float, float, float],
        scoring: str = "vinardo",
        exhaustiveness: int = 16
    ) -> dict:
        """Refine pose with Smina local search"""
        logger.info(f"Refining pose with Smina (exhaustiveness={exhaustiveness})")
        
        output_dir = self.get_output_path("smina_refinement")
        
        result = self.smina.dock(
            receptor=receptor_pdb,
            ligand=ligand_pdb,
            center=center,
            size=size,
            output_dir=output_dir,
            scoring=scoring,
            exhaustiveness=exhaustiveness,
            num_modes=5
        )
        
        logger.info(f"Smina refinement completed")
        return result


async def main():
    """Run complex server"""
    server = ComplexServer()
    await server.run()


if __name__ == "__main__":
    import asyncio
    asyncio.run(main())

