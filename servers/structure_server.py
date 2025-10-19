"""
Structure Server - PDB retrieval and structure cleaning.

Provides MCP tools for:
- Fetching structures from PDB/AlphaFold
- PDBFixer structure cleaning
- Protonation with PDB2PQR
- Structure validation
"""

import logging
import json
from pathlib import Path
from typing import Optional, List, Dict, Any
from mcp.server import Server
from mcp.types import Tool, TextContent
import httpx

from .base_server import BaseMCPServer
from tools.pdbfixer_wrapper import PDBFixerWrapper
from tools.pdb2pqr_wrapper import PDB2PQRWrapper
from core.utils import setup_logger, ensure_directory, count_atoms_in_pdb, get_pdb_chains

logger = setup_logger(__name__)


class StructureServer(BaseMCPServer):
    """MCP Server for structure retrieval and cleaning"""
    
    def __init__(self):
        super().__init__("structure_server", "0.1.0")
        self.pdbfixer = PDBFixerWrapper()
        self.pdb2pqr = PDB2PQRWrapper()
        self.setup_handlers()
    
    def setup_handlers(self):
        """Register MCP tool handlers"""
        
        @self.server.list_tools()
        async def list_tools() -> list[Tool]:
            return [
                Tool(
                    name="fetch_pdb",
                    description="Fetch PDB structure from RCSB PDB, AlphaFold DB, or PDB-REDO",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "pdb_id": {"type": "string", "description": "PDB ID (e.g., 1ABC)"},
                            "source": {
                                "type": "string",
                                "enum": ["pdb", "alphafold", "pdb-redo"],
                                "default": "pdb",
                                "description": "Source database"
                            }
                        },
                        "required": ["pdb_id"]
                    }
                ),
                Tool(
                    name="clean_structure",
                    description="Clean PDB structure with PDBFixer",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "pdb_file": {"type": "string"},
                            "remove_water": {"type": "boolean", "default": True},
                            "fix_missing": {"type": "boolean", "default": True}
                        },
                        "required": ["pdb_file"]
                    }
                ),
                Tool(
                    name="add_hydrogens",
                    description="Add hydrogens to structure",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "pdb_file": {"type": "string"},
                            "ph": {"type": "number", "default": 7.0}
                        },
                        "required": ["pdb_file"]
                    }
                ),
                Tool(
                    name="protonate_structure",
                    description="Protonate structure using PDB2PQR+PROPKA at specified pH",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "pdb_file": {"type": "string"},
                            "ph": {"type": "number", "default": 7.0},
                            "forcefield": {"type": "string", "default": "AMBER"}
                        },
                        "required": ["pdb_file"]
                    }
                ),
                Tool(
                    name="detect_modifications",
                    description="Detect disulfide bonds and modifications",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "pdb_file": {"type": "string"}
                        },
                        "required": ["pdb_file"]
                    }
                ),
                Tool(
                    name="validate_structure",
                    description="Validate PDB structure",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "pdb_file": {"type": "string"}
                        },
                        "required": ["pdb_file"]
                    }
                ),
            ]
        
        @self.server.call_tool()
        async def call_tool(name: str, arguments: dict) -> list[TextContent]:
            try:
                if name == "fetch_pdb":
                    result = await self.fetch_pdb(
                        pdb_id=arguments["pdb_id"],
                        source=arguments.get("source", "pdb")
                    )
                elif name == "clean_structure":
                    result = await self.clean_structure(
                        pdb_file=arguments["pdb_file"],
                        remove_water=arguments.get("remove_water", True),
                        fix_missing=arguments.get("fix_missing", True)
                    )
                elif name == "add_hydrogens":
                    result = await self.add_hydrogens(
                        pdb_file=arguments["pdb_file"],
                        ph=arguments.get("ph", 7.0)
                    )
                elif name == "protonate_structure":
                    result = await self.protonate_structure(
                        pdb_file=arguments["pdb_file"],
                        ph=arguments.get("ph", 7.0),
                        forcefield=arguments.get("forcefield", "AMBER")
                    )
                elif name == "detect_modifications":
                    result = await self.detect_modifications(
                        pdb_file=arguments["pdb_file"]
                    )
                elif name == "validate_structure":
                    result = await self.validate_structure(
                        pdb_file=arguments["pdb_file"]
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
    
    async def fetch_pdb(self, pdb_id: str, source: str = "pdb") -> dict:
        """Fetch PDB structure from database"""
        logger.info(f"Fetching {pdb_id} from {source}")
        
        pdb_id = pdb_id.upper()
        output_file = self.get_output_path(f"{pdb_id}.pdb")
        
        if source == "pdb":
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        elif source == "alphafold":
            # AlphaFold uses different ID format
            url = f"https://alphafold.ebi.ac.uk/files/AF-{pdb_id}-F1-model_v4.pdb"
        elif source == "pdb-redo":
            url = f"https://pdb-redo.eu/db/{pdb_id}/{pdb_id}_final.pdb"
        else:
            raise ValueError(f"Unknown source: {source}")
        
        # Download
        async with httpx.AsyncClient() as client:
            response = await client.get(url)
            response.raise_for_status()
            
            with open(output_file, 'wb') as f:
                f.write(response.content)
        
        logger.info(f"Downloaded {pdb_id} to {output_file}")
        
        # Get basic info
        num_atoms = count_atoms_in_pdb(output_file)
        chains = get_pdb_chains(output_file)
        
        return {
            "pdb_id": pdb_id,
            "source": source,
            "file_path": str(output_file),
            "num_atoms": num_atoms,
            "chains": chains
        }
    
    async def clean_structure(
        self,
        pdb_file: str,
        remove_water: bool = True,
        fix_missing: bool = True
    ) -> dict:
        """Clean PDB structure"""
        logger.info(f"Cleaning structure: {pdb_file}")
        
        if not self.validate_input_file(pdb_file):
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")
        
        output_file = self.get_output_path("cleaned.pdb")
        
        result = self.pdbfixer.clean_structure(
            input_pdb=pdb_file,
            output_pdb=output_file,
            remove_water=remove_water,
            fix_missing=fix_missing
        )
        
        return result
    
    async def add_hydrogens(self, pdb_file: str, ph: float = 7.0) -> dict:
        """Add hydrogens to structure"""
        logger.info(f"Adding hydrogens at pH {ph}")
        
        if not self.validate_input_file(pdb_file):
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")
        
        output_file = self.get_output_path("protonated.pdb")
        
        result = self.pdbfixer.add_hydrogens_only(
            input_pdb=pdb_file,
            output_pdb=output_file,
            ph=ph
        )
        
        return result
    
    async def protonate_structure(
        self,
        pdb_file: str,
        ph: float = 7.0,
        forcefield: str = "AMBER"
    ) -> dict:
        """Protonate structure using PDB2PQR+PROPKA"""
        logger.info(f"Protonating structure at pH {ph} with PDB2PQR")
        
        if not self.validate_input_file(pdb_file):
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")
        
        output_file = self.get_output_path("protonated_pdb2pqr.pdb")
        
        result = self.pdb2pqr.protonate_structure(
            input_pdb=pdb_file,
            output_pdb=output_file,
            ph=ph,
            forcefield=forcefield
        )
        
        return result
    
    async def detect_modifications(self, pdb_file: str) -> dict:
        """Detect disulfide bonds and modifications"""
        logger.info(f"Detecting modifications: {pdb_file}")
        
        if not self.validate_input_file(pdb_file):
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")
        
        modifications = {
            "file_path": pdb_file,
            "disulfide_bonds": [],
            "modified_residues": [],
            "metal_sites": []
        }
        
        # Parse PDB for modifications
        with open(pdb_file, 'r') as f:
            for line in f:
                # Detect SSBOND records
                if line.startswith('SSBOND'):
                    parts = line.split()
                    if len(parts) >= 7:
                        bond = {
                            "res1": f"{parts[2]}_{parts[3]}",
                            "res2": f"{parts[5]}_{parts[6]}"
                        }
                        modifications["disulfide_bonds"].append(bond)
                
                # Detect modified residues (MODRES)
                elif line.startswith('MODRES'):
                    parts = line.split()
                    if len(parts) >= 5:
                        mod = {
                            "residue": parts[2],
                            "chain": parts[3],
                            "resnum": parts[4],
                            "standard": parts[5] if len(parts) > 5 else None
                        }
                        modifications["modified_residues"].append(mod)
                
                # Detect metal ions
                elif line.startswith('HETATM'):
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    
                    # Common metal ions
                    metals = ['ZN', 'MG', 'CA', 'FE', 'CU', 'MN', 'NA', 'K']
                    if res_name in metals or atom_name in metals:
                        chain = line[21:22]
                        resnum = line[22:26].strip()
                        metal = {
                            "element": res_name,
                            "chain": chain,
                            "resnum": resnum
                        }
                        if metal not in modifications["metal_sites"]:
                            modifications["metal_sites"].append(metal)
        
        logger.info(f"Found {len(modifications['disulfide_bonds'])} disulfide bonds")
        logger.info(f"Found {len(modifications['modified_residues'])} modified residues")
        logger.info(f"Found {len(modifications['metal_sites'])} metal sites")
        
        return modifications
    
    async def validate_structure(self, pdb_file: str) -> dict:
        """Validate PDB structure"""
        logger.info(f"Validating structure: {pdb_file}")
        
        if not self.validate_input_file(pdb_file):
            return {"valid": False, "error": "File not found"}
        
        # Basic validation
        num_atoms = count_atoms_in_pdb(pdb_file)
        chains = get_pdb_chains(pdb_file)
        
        validation = {
            "valid": num_atoms > 0,
            "file_path": pdb_file,
            "num_atoms": num_atoms,
            "chains": chains
        }
        
        if num_atoms == 0:
            validation["error"] = "No atoms found in PDB file"
        
        return validation


async def main():
    """Run structure server"""
    server = StructureServer()
    await server.run()


if __name__ == "__main__":
    import asyncio
    asyncio.run(main())

