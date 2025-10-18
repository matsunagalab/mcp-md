"""
Ligand Server - RDKit 3D generation and AmberTools parameterization.

Provides MCP tools for:
- SMILES to 3D structure generation
- GAFF2 parameterization with AM1-BCC charges
- Ligand library creation for tleap
"""

import logging
import json
from pathlib import Path
from typing import Optional
from mcp.server import Server
from mcp.types import Tool, TextContent

from .base_server import BaseMCPServer
from tools.rdkit_wrapper import RDKitWrapper
from tools.ambertools_wrapper import AmberToolsWrapper
from core.utils import setup_logger

logger = setup_logger(__name__)


class LigandServer(BaseMCPServer):
    """MCP Server for ligand operations"""
    
    def __init__(self):
        super().__init__("ligand_server", "0.1.0")
        self.rdkit = RDKitWrapper()
        self.ambertools = AmberToolsWrapper()
        self.setup_handlers()
    
    def setup_handlers(self):
        """Register MCP tool handlers"""
        
        @self.server.list_tools()
        async def list_tools() -> list[Tool]:
            return [
                Tool(
                    name="smiles_to_3d",
                    description="Generate 3D structure from SMILES using RDKit",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "smiles": {"type": "string", "description": "SMILES string"},
                            "optimize": {"type": "boolean", "default": True},
                            "method": {
                                "type": "string",
                                "enum": ["mmff", "uff"],
                                "default": "mmff"
                            }
                        },
                        "required": ["smiles"]
                    }
                ),
                Tool(
                    name="generate_gaff_params",
                    description="Generate GAFF2 parameters with AM1-BCC charges (AmberTools)",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "ligand_file": {"type": "string", "description": "Input structure file"},
                            "net_charge": {"type": "integer", "default": 0},
                            "charge_method": {
                                "type": "string",
                                "enum": ["bcc", "gas", "resp"],
                                "default": "bcc",
                                "description": "bcc=AM1-BCC, gas=Gasteiger, resp=RESP"
                            },
                            "residue_name": {"type": "string", "default": "LIG"}
                        },
                        "required": ["ligand_file"]
                    }
                ),
                Tool(
                    name="create_ligand_lib",
                    description="Create tleap library file for ligand",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "mol2_file": {"type": "string"},
                            "frcmod_file": {"type": "string"},
                            "residue_name": {"type": "string", "default": "LIG"}
                        },
                        "required": ["mol2_file", "frcmod_file"]
                    }
                ),
                Tool(
                    name="parameterize_ligand_complete",
                    description="Complete ligand parameterization: SMILES -> 3D -> GAFF2 -> library",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "smiles": {"type": "string"},
                            "net_charge": {"type": "integer", "default": 0},
                            "residue_name": {"type": "string", "default": "LIG"},
                            "charge_method": {"type": "string", "default": "bcc"}
                        },
                        "required": ["smiles"]
                    }
                ),
            ]
        
        @self.server.call_tool()
        async def call_tool(name: str, arguments: dict) -> list[TextContent]:
            try:
                if name == "smiles_to_3d":
                    result = await self.smiles_to_3d(
                        smiles=arguments["smiles"],
                        optimize=arguments.get("optimize", True),
                        method=arguments.get("method", "mmff")
                    )
                elif name == "generate_gaff_params":
                    result = await self.generate_gaff_params(
                        ligand_file=arguments["ligand_file"],
                        net_charge=arguments.get("net_charge", 0),
                        charge_method=arguments.get("charge_method", "bcc"),
                        residue_name=arguments.get("residue_name", "LIG")
                    )
                elif name == "create_ligand_lib":
                    result = await self.create_ligand_lib(
                        mol2_file=arguments["mol2_file"],
                        frcmod_file=arguments["frcmod_file"],
                        residue_name=arguments.get("residue_name", "LIG")
                    )
                elif name == "parameterize_ligand_complete":
                    result = await self.parameterize_ligand_complete(
                        smiles=arguments["smiles"],
                        net_charge=arguments.get("net_charge", 0),
                        residue_name=arguments.get("residue_name", "LIG"),
                        charge_method=arguments.get("charge_method", "bcc")
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
    
    async def smiles_to_3d(
        self,
        smiles: str,
        optimize: bool = True,
        method: str = "mmff"
    ) -> dict:
        """Generate 3D structure from SMILES"""
        logger.info(f"Generating 3D structure from SMILES: {smiles}")
        
        # Validate SMILES
        if not self.rdkit.validate_smiles(smiles):
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        # Get formal charge
        formal_charge = self.rdkit.get_formal_charge(smiles)
        
        output_file = self.get_output_path("ligand_3d.mol2")
        
        result = self.rdkit.smiles_to_3d(
            smiles=smiles,
            output_file=output_file,
            optimize=optimize,
            method=method
        )
        
        result["formal_charge"] = formal_charge
        
        return result
    
    async def generate_gaff_params(
        self,
        ligand_file: str,
        net_charge: int = 0,
        charge_method: str = "bcc",
        residue_name: str = "LIG"
    ) -> dict:
        """Generate GAFF2 parameters"""
        logger.info(f"Generating GAFF2 parameters for {ligand_file}")
        
        if not self.validate_input_file(ligand_file):
            raise FileNotFoundError(f"Ligand file not found: {ligand_file}")
        
        output_dir = self.get_output_path("gaff_params")
        
        result = self.ambertools.generate_gaff_params(
            input_file=ligand_file,
            output_dir=output_dir,
            net_charge=net_charge,
            charge_method=charge_method,
            residue_name=residue_name
        )
        
        return result
    
    async def create_ligand_lib(
        self,
        mol2_file: str,
        frcmod_file: str,
        residue_name: str = "LIG"
    ) -> dict:
        """Create tleap library file"""
        logger.info(f"Creating ligand library for {residue_name}")
        
        if not self.validate_input_file(mol2_file):
            raise FileNotFoundError(f"MOL2 file not found: {mol2_file}")
        
        if not self.validate_input_file(frcmod_file):
            raise FileNotFoundError(f"FRCMOD file not found: {frcmod_file}")
        
        output_dir = self.get_output_path("ligand_lib")
        
        result = self.ambertools.create_ligand_lib(
            mol2_file=mol2_file,
            frcmod_file=frcmod_file,
            output_dir=output_dir,
            residue_name=residue_name
        )
        
        return result
    
    async def parameterize_ligand_complete(
        self,
        smiles: str,
        net_charge: int = 0,
        residue_name: str = "LIG",
        charge_method: str = "bcc"
    ) -> dict:
        """Complete ligand parameterization workflow"""
        logger.info(f"Complete parameterization workflow for: {smiles}")
        
        # Step 1: SMILES to 3D
        struct_3d = await self.smiles_to_3d(smiles=smiles, optimize=True)
        logger.info(f"Generated 3D structure: {struct_3d['output_file']}")
        
        # Use detected formal charge if net_charge not specified
        if net_charge == 0 and "formal_charge" in struct_3d:
            net_charge = struct_3d["formal_charge"]
            logger.info(f"Using detected formal charge: {net_charge}")
        
        # Step 2: Generate GAFF parameters
        gaff_params = await self.generate_gaff_params(
            ligand_file=struct_3d["output_file"],
            net_charge=net_charge,
            charge_method=charge_method,
            residue_name=residue_name
        )
        logger.info(f"Generated GAFF parameters: {gaff_params['mol2']}")
        
        # Step 3: Create library
        lib_result = await self.create_ligand_lib(
            mol2_file=gaff_params["mol2"],
            frcmod_file=gaff_params["frcmod"],
            residue_name=residue_name
        )
        logger.info(f"Created library: {lib_result['lib']}")
        
        # Combine results
        complete_result = {
            "smiles": smiles,
            "net_charge": net_charge,
            "residue_name": residue_name,
            "charge_method": charge_method,
            "initial_3d": struct_3d["output_file"],
            "gaff_mol2": gaff_params["mol2"],
            "frcmod": gaff_params["frcmod"],
            "library": lib_result["lib"],
            "pdb": lib_result["pdb"],
            "charges": gaff_params["charges"],
            "total_charge": gaff_params["total_charge"]
        }
        
        return complete_result


async def main():
    """Run ligand server"""
    server = LigandServer()
    await server.run()


if __name__ == "__main__":
    import asyncio
    asyncio.run(main())

