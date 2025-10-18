"""
OpenMM wrapper for MD simulation script generation.
"""

import logging
from pathlib import Path
from typing import Union, Dict

logger = logging.getLogger(__name__)


class OpenMMWrapper:
    """Wrapper for OpenMM MD script generation"""
    
    def __init__(self):
        logger.info("OpenMM wrapper initialized")
    
    def create_workflow(
        self,
        prmtop: Union[str, Path],
        inpcrd: Union[str, Path],
        output_dir: Union[str, Path],
        protocol: str = "standard"
    ) -> Dict:
        """Create complete OpenMM MD workflow
        
        Args:
            prmtop: Topology file
            inpcrd: Coordinate file
            output_dir: Output directory
            protocol: Protocol type (standard, membrane, implicit)
        
        Returns:
            Dict with script paths
        """
        logger.info(f"Creating OpenMM workflow: {protocol}")
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate scripts
        minimize_script = self._create_minimization_script(prmtop, inpcrd, output_dir)
        equilibrate_script = self._create_equilibration_script(output_dir)
        production_script = self._create_production_script(output_dir)
        submit_script = self._create_submit_script(output_dir)
        
        return {
            "minimize_script": str(minimize_script),
            "equilibrate_script": str(equilibrate_script),
            "production_script": str(production_script),
            "submit_script": str(submit_script),
            "protocol": protocol
        }
    
    def _create_minimization_script(self, prmtop, inpcrd, output_dir) -> Path:
        """Create minimization script"""
        script_path = output_dir / "1_minimize.py"
        
        script = f"""#!/usr/bin/env python
'''OpenMM Minimization Script'''

from openmm.app import *
from openmm import *
from openmm.unit import *

# Load system
prmtop = AmberPrmtopFile('{prmtop}')
inpcrd = AmberInpcrdFile('{inpcrd}')

# Create system
system = prmtop.createSystem(
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*nanometer,
    constraints=HBonds
)

# Minimization
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)

print('Minimizing...')
simulation.minimizeEnergy(maxIterations=5000)

# Save
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('minimized.pdb', 'w'))

print('Minimization complete')
"""
        
        with open(script_path, 'w') as f:
            f.write(script)
        
        script_path.chmod(0o755)
        return script_path
    
    def _create_equilibration_script(self, output_dir) -> Path:
        """Create equilibration script"""
        script_path = output_dir / "2_equilibrate.py"
        
        script = """#!/usr/bin/env python
'''OpenMM NPT Equilibration Script'''

from openmm.app import *
from openmm import *
from openmm.unit import *

# Load minimized structure
pdb = PDBFile('minimized.pdb')
prmtop = AmberPrmtopFile('system.prmtop')

# Create system
system = prmtop.createSystem(
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*nanometer,
    constraints=HBonds
)

# Add barostat
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))

# Setup simulation
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

# Reporters
simulation.reporters.append(StateDataReporter('equilibration.log', 1000,
    step=True, time=True, temperature=True, volume=True))

# Run equilibration (500 ps)
print('Equilibrating...')
simulation.step(250000)

# Save
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('equilibrated.pdb', 'w'))

print('Equilibration complete')
"""
        
        with open(script_path, 'w') as f:
            f.write(script)
        
        script_path.chmod(0o755)
        return script_path
    
    def _create_production_script(self, output_dir) -> Path:
        """Create production MD script"""
        script_path = output_dir / "3_production.py"
        
        script = """#!/usr/bin/env python
'''OpenMM Production MD Script'''

from openmm.app import *
from openmm import *
from openmm.unit import *

# Load equilibrated structure
pdb = PDBFile('equilibrated.pdb')
prmtop = AmberPrmtopFile('system.prmtop')

# Create system
system = prmtop.createSystem(
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*nanometer,
    constraints=HBonds
)

# Add barostat
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))

# Setup simulation
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

# Reporters
simulation.reporters.append(DCDReporter('trajectory.dcd', 5000))
simulation.reporters.append(StateDataReporter('production.log', 1000,
    step=True, time=True, temperature=True, potentialEnergy=True, volume=True))

# Run production (10 ns)
print('Running production MD...')
simulation.step(5000000)

print('Production MD complete')
"""
        
        with open(script_path, 'w') as f:
            f.write(script)
        
        script_path.chmod(0o755)
        return script_path
    
    def _create_submit_script(self, output_dir) -> Path:
        """Create job submission script"""
        script_path = output_dir / "submit.sh"
        
        script = """#!/bin/bash
# OpenMM MD Workflow Submission Script

echo "Starting MD workflow..."

# 1. Minimization
echo "Step 1: Minimization"
python 1_minimize.py

# 2. Equilibration
echo "Step 2: Equilibration"
python 2_equilibrate.py

# 3. Production
echo "Step 3: Production MD"
python 3_production.py

echo "MD workflow complete!"
"""
        
        with open(script_path, 'w') as f:
            f.write(script)
        
        script_path.chmod(0o755)
        return script_path

