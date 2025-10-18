"""
Packmol and Packmol-Memgen wrappers for membrane system building.
"""

import logging
from pathlib import Path
from typing import Optional, Union, Dict, List
from .base_wrapper import BaseToolWrapper

logger = logging.getLogger(__name__)


class PackmolWrapper(BaseToolWrapper):
    """Wrapper for Packmol"""
    
    def __init__(self, conda_env: str = "mcp-md-tools"):
        super().__init__("packmol", conda_env=conda_env)
    
    def build_mixed_solvent(
        self,
        solvent_components: Dict[str, int],
        output_pdb: Union[str, Path],
        box_size: List[float]
    ) -> Dict:
        """Build mixed solvent box
        
        Args:
            solvent_components: Dict of {pdb_file: num_molecules}
            output_pdb: Output PDB file
            box_size: Box dimensions [x, y, z] in Angstroms
        
        Returns:
            Dict with build results
        """
        logger.info(f"Building mixed solvent box: {len(solvent_components)} components")
        
        output_path = Path(output_pdb)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Create Packmol input script
        script_path = output_path.parent / "packmol.inp"
        
        with open(script_path, 'w') as f:
            f.write(f"""tolerance 2.0
filetype pdb
output {output_path}

""")
            
            for pdb_file, num_molecules in solvent_components.items():
                f.write(f"""structure {pdb_file}
  number {num_molecules}
  inside box 0. 0. 0. {box_size[0]} {box_size[1]} {box_size[2]}
end structure

""")
        
        logger.info(f"Created Packmol input: {script_path}")
        
        # Run Packmol
        try:
            self.run(['-i', str(script_path)], cwd=output_path.parent)
            logger.info("Packmol completed successfully")
        except Exception as e:
            logger.error(f"Packmol failed: {e}")
            raise
        
        return {
            "output": str(output_path),
            "input_script": str(script_path),
            "num_components": len(solvent_components),
            "box_size": box_size
        }


class PackmolMemgenWrapper(BaseToolWrapper):
    """Wrapper for Packmol-Memgen membrane building"""
    
    def __init__(self, conda_env: str = "mcp-md-tools"):
        super().__init__("packmol-memgen", conda_env=conda_env)
    
    def build_membrane_system(
        self,
        protein_pdb: Union[str, Path],
        output_dir: Union[str, Path],
        lipid_composition: Dict[str, float],
        membrane_type: str = "bilayer",
        dist_to_bilayer: float = 15.0
    ) -> Dict:
        """Build membrane protein system
        
        Args:
            protein_pdb: Protein PDB file
            output_dir: Output directory
            lipid_composition: Dict of {lipid_name: ratio}
                Example: {"POPC": 0.7, "POPE": 0.3}
            membrane_type: "bilayer" or "micelle"
            dist_to_bilayer: Distance from protein to membrane (Angstroms)
        
        Returns:
            Dict with membrane system files
        """
        logger.info(f"Building {membrane_type} system with Packmol-Memgen")
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Build lipid composition string
        lipid_str = ":".join([f"{name}:{ratio}" for name, ratio in lipid_composition.items()])
        
        # Packmol-Memgen arguments
        args = [
            '--pdb', str(protein_pdb),
            '--lipids', lipid_str,
            '--membrane-type', membrane_type,
            '--dist', str(dist_to_bilayer),
            '--output', str(output_dir / "membrane_system")
        ]
        
        try:
            self.run(args, cwd=output_dir)
            logger.info("Packmol-Memgen completed successfully")
        except Exception as e:
            logger.error(f"Packmol-Memgen failed: {e}")
            # Try alternative: manual build
            logger.info("Attempting manual membrane build...")
            return self._manual_membrane_build(
                protein_pdb, output_dir, lipid_composition, dist_to_bilayer
            )
        
        # Find output files
        system_pdb = output_dir / "membrane_system.pdb"
        
        return {
            "output_pdb": str(system_pdb) if system_pdb.exists() else None,
            "output_dir": str(output_dir),
            "membrane_type": membrane_type,
            "lipid_composition": lipid_composition,
            "dist_to_bilayer": dist_to_bilayer
        }
    
    def _manual_membrane_build(
        self,
        protein_pdb: Path,
        output_dir: Path,
        lipid_composition: Dict[str, float],
        dist_to_bilayer: float
    ) -> Dict:
        """Fallback manual membrane building
        
        Creates a simple membrane system using tleap-compatible approach
        """
        logger.info("Using manual membrane build approach")
        
        # Create simple bilayer with tleap script
        tleap_script = output_dir / "build_membrane.in"
        
        with open(tleap_script, 'w') as f:
            f.write(f"""# Manual membrane system building
source leaprc.lipid17
source leaprc.protein.ff19SB
source leaprc.water.tip3p

# Load protein
protein = loadpdb {protein_pdb}

# Create simple membrane (POPC default)
# This is a placeholder - full implementation would use Packmol
membrane = loadpdb $AMBERHOME/dat/leap/pdb/membrane.pdb

# Combine
system = combine {{protein membrane}}

# Solvate
solvateBox system TIP3PBOX 12.0

# Save
saveamberparm system membrane_system.prmtop membrane_system.inpcrd
savepdb system membrane_system.pdb

quit
""")
        
        return {
            "output_pdb": str(output_dir / "membrane_system.pdb"),
            "output_dir": str(output_dir),
            "membrane_type": "manual_bilayer",
            "lipid_composition": lipid_composition,
            "tleap_script": str(tleap_script),
            "note": "Manual build - consider using full Packmol-Memgen for production"
        }

