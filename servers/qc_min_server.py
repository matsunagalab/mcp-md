"""
QC/Min Server - Structure quality control and minimization.

Provides MCP tools for:
- Custom MolProbity checks (clashes, bonds, chirality)
- OpenMM minimization
- QC report generation
"""

import logging
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, Any
from mcp.types import Tool, TextContent

from .base_server import BaseMCPServer
from tools.molprobity_wrapper import MolProbityWrapper
from core.utils import setup_logger

logger = setup_logger(__name__)


class QCMinServer(BaseMCPServer):
    """MCP Server for quality control and minimization"""
    
    def __init__(self):
        super().__init__("qc_min_server", "0.1.0")
        self.molprobity = MolProbityWrapper()
        self.setup_handlers()
    
    def setup_handlers(self):
        """Register MCP tool handlers"""
        
        @self.server.list_tools()
        async def list_tools() -> list[Tool]:
            return [
                Tool(
                    name="molprobity_check",
                    description="Run custom MolProbity quality checks (clashes, bonds, chirality)",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "pdb_file": {"type": "string", "description": "PDB file to check"},
                            "clash_cutoff": {"type": "number", "default": 2.0, "description": "Clash distance cutoff (Å)"}
                        },
                        "required": ["pdb_file"]
                    }
                ),
                Tool(
                    name="openmm_minimize",
                    description="Minimize structure with OpenMM",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "prmtop": {"type": "string", "description": "Amber topology file"},
                            "inpcrd": {"type": "string", "description": "Amber coordinate file"},
                            "max_iterations": {"type": "integer", "default": 5000},
                            "tolerance": {"type": "number", "default": 10.0, "description": "Energy tolerance (kJ/mol)"}
                        },
                        "required": ["prmtop", "inpcrd"]
                    }
                ),
                Tool(
                    name="qc_report",
                    description="Generate comprehensive QC report (JSON + Markdown)",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "pdb_file": {"type": "string", "description": "PDB file to analyze"},
                            "prmtop": {"type": "string", "description": "Topology file (optional)"},
                            "minimize": {"type": "boolean", "default": False, "description": "Run minimization before QC"}
                        },
                        "required": ["pdb_file"]
                    }
                ),
            ]
        
        @self.server.call_tool()
        async def call_tool(name: str, arguments: dict) -> list[TextContent]:
            """Handle tool calls"""
            try:
                logger.info(f"Calling tool: {name}")
                
                if name == "molprobity_check":
                    result = await self.molprobity_check(
                        pdb_file=arguments["pdb_file"],
                        clash_cutoff=arguments.get("clash_cutoff", 2.0)
                    )
                elif name == "openmm_minimize":
                    result = await self.openmm_minimize(
                        prmtop=arguments["prmtop"],
                        inpcrd=arguments["inpcrd"],
                        max_iterations=arguments.get("max_iterations", 5000),
                        tolerance=arguments.get("tolerance", 10.0)
                    )
                elif name == "qc_report":
                    result = await self.qc_report(
                        pdb_file=arguments["pdb_file"],
                        prmtop=arguments.get("prmtop"),
                        minimize=arguments.get("minimize", False)
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
    
    async def molprobity_check(self, pdb_file: str, clash_cutoff: float = 2.0) -> dict:
        """Run MolProbity quality checks"""
        logger.info(f"Running MolProbity checks on {pdb_file}")
        
        result = self.molprobity.run_full_qc(pdb_file)
        result["clash_cutoff"] = clash_cutoff
        
        logger.info(f"MolProbity check completed: {result.get('overall_status')}")
        return result
    
    async def openmm_minimize(
        self,
        prmtop: str,
        inpcrd: str,
        max_iterations: int = 5000,
        tolerance: float = 10.0
    ) -> dict:
        """Minimize structure with OpenMM"""
        logger.info(f"Minimizing with OpenMM (max_iter={max_iterations})")
        
        output_dir = self.get_output_path("minimization")
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create minimization script
        script_path = output_dir / "minimize.py"
        pdb_output = output_dir / "minimized.pdb"
        
        script_content = f"""#!/usr/bin/env python
'''OpenMM Minimization Script'''

from openmm.app import *
from openmm import *
from openmm.unit import *
import sys

try:
    # Load system
    prmtop = AmberPrmtopFile('{prmtop}')
    inpcrd = AmberInpcrdFile('{inpcrd}')
    
    # Create system
    system = prmtop.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=1.0*nanometer,
        constraints=HBonds
    )
    
    # Setup minimization
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)
    
    # Get initial energy
    state = simulation.context.getState(getEnergy=True)
    initial_energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    print(f'Initial energy: {{initial_energy:.2f}} kJ/mol')
    
    # Minimize
    print('Minimizing...')
    simulation.minimizeEnergy(
        maxIterations={max_iterations},
        tolerance={tolerance}*kilojoules_per_mole
    )
    
    # Get final energy
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    final_energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    print(f'Final energy: {{final_energy:.2f}} kJ/mol')
    print(f'Energy change: {{final_energy - initial_energy:.2f}} kJ/mol')
    
    # Save
    positions = state.getPositions()
    PDBFile.writeFile(simulation.topology, positions, open('{pdb_output}', 'w'))
    
    # Write energy to file
    with open('{output_dir}/energy.txt', 'w') as f:
        f.write(f'Initial: {{initial_energy:.2f}}\\n')
        f.write(f'Final: {{final_energy:.2f}}\\n')
        f.write(f'Change: {{final_energy - initial_energy:.2f}}\\n')
    
    print('Minimization complete')
    sys.exit(0)

except Exception as e:
    print(f'Error: {{e}}', file=sys.stderr)
    sys.exit(1)
"""
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        script_path.chmod(0o755)
        
        # Run minimization
        import subprocess
        try:
            result = subprocess.run(
                ["python", str(script_path)],
                capture_output=True,
                text=True,
                timeout=300,
                check=True
            )
            
            # Read energy file
            energy_file = output_dir / "energy.txt"
            if energy_file.exists():
                with open(energy_file) as f:
                    energy_lines = f.readlines()
                    energies = {
                        "initial": float(energy_lines[0].split(":")[1].strip()),
                        "final": float(energy_lines[1].split(":")[1].strip()),
                        "change": float(energy_lines[2].split(":")[1].strip())
                    }
            else:
                energies = {}
            
            return {
                "status": "success",
                "minimized_pdb": str(pdb_output),
                "script": str(script_path),
                "energies": energies,
                "stdout": result.stdout
            }
        
        except subprocess.TimeoutExpired:
            logger.error("Minimization timeout")
            return {"error": "Minimization timeout", "status": "timeout"}
        except subprocess.CalledProcessError as e:
            logger.error(f"Minimization failed: {e.stderr}")
            return {"error": e.stderr, "status": "failed"}
    
    async def qc_report(
        self,
        pdb_file: str,
        prmtop: str = None,
        minimize: bool = False
    ) -> dict:
        """Generate comprehensive QC report"""
        logger.info(f"Generating QC report for {pdb_file}")
        
        output_dir = self.get_output_path("qc_report")
        output_dir.mkdir(parents=True, exist_ok=True)
        
        report = {
            "timestamp": datetime.utcnow().isoformat(),
            "pdb_file": pdb_file,
            "minimized": minimize
        }
        
        # Run MolProbity checks
        molprobity_result = await self.molprobity_check(pdb_file)
        report["molprobity"] = molprobity_result
        
        # Run minimization if requested
        if minimize and prmtop:
            # Need inpcrd from PDB
            inpcrd = pdb_file.replace(".pdb", ".inpcrd")
            if Path(inpcrd).exists():
                min_result = await self.openmm_minimize(prmtop, inpcrd)
                report["minimization"] = min_result
                
                # Re-run MolProbity on minimized structure
                if min_result.get("status") == "success":
                    minimized_pdb = min_result["minimized_pdb"]
                    molprobity_after = await self.molprobity_check(minimized_pdb)
                    report["molprobity_after_minimization"] = molprobity_after
        
        # Save JSON report
        json_path = output_dir / "qc_report.json"
        with open(json_path, 'w') as f:
            json.dump(report, f, indent=2)
        
        # Generate Markdown report
        md_path = output_dir / "qc_report.md"
        md_content = self._generate_markdown_report(report)
        with open(md_path, 'w') as f:
            f.write(md_content)
        
        report["json_report"] = str(json_path)
        report["markdown_report"] = str(md_path)
        
        logger.info(f"QC report generated: {json_path}")
        return report
    
    def _generate_markdown_report(self, report: dict) -> str:
        """Generate Markdown QC report"""
        md = f"""# QC Report

**Timestamp**: {report['timestamp']}
**PDB File**: {report['pdb_file']}
**Minimized**: {report['minimized']}

## MolProbity Analysis

**Overall Status**: {report['molprobity'].get('overall_status', 'unknown').upper()}

### Clashes
- **Number of clashes**: {report['molprobity']['clashes'].get('num_clashes', 0)}
- **Status**: {report['molprobity']['clashes'].get('status', 'unknown')}

### Bond Lengths
- **Number of outliers**: {report['molprobity']['bond_lengths'].get('num_outliers', 0)}
- **Status**: {report['molprobity']['bond_lengths'].get('status', 'unknown')}

### Chirality
- **Wrong chirality**: {report['molprobity']['chirality'].get('num_wrong_chirality', 0)}
- **Status**: {report['molprobity']['chirality'].get('status', 'unknown')}

"""
        
        # Add minimization section if present
        if "minimization" in report:
            min_data = report["minimization"]
            if "energies" in min_data:
                energies = min_data["energies"]
                md += f"""## Minimization

- **Initial Energy**: {energies.get('initial', 'N/A'):.2f} kJ/mol
- **Final Energy**: {energies.get('final', 'N/A'):.2f} kJ/mol
- **Energy Change**: {energies.get('change', 'N/A'):.2f} kJ/mol

"""
        
        return md


async def main():
    """Run QC/Min server"""
    server = QCMinServer()
    await server.run()


if __name__ == "__main__":
    import asyncio
    asyncio.run(main())

