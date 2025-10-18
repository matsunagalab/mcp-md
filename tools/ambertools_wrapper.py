"""
AmberTools wrapper for ligand parameterization and system building.

Provides interface to:
- antechamber: GAFF parameterization and charge calculation
- parmchk2: Missing parameter generation
- tleap: Topology building
"""

import logging
import re
from pathlib import Path
from typing import Optional, Union, List, Dict
from .base_wrapper import BaseToolWrapper

logger = logging.getLogger(__name__)


class AmberToolsWrapper(BaseToolWrapper):
    """Wrapper for AmberTools suite"""
    
    def __init__(self, conda_env: str = "mcp-md-tools"):
        # AmberTools is typically in conda environment
        super().__init__("antechamber", conda_env=conda_env)
        self.conda_env = conda_env
    
    def generate_gaff_params(
        self,
        input_file: Union[str, Path],
        output_dir: Union[str, Path],
        net_charge: int = 0,
        charge_method: str = "bcc",
        atom_type: str = "gaff2",
        residue_name: str = "LIG"
    ) -> Dict:
        """Generate GAFF parameters for ligand
        
        Args:
            input_file: Input structure file (mol2, pdb, sdf)
            output_dir: Output directory
            net_charge: Net molecular charge
            charge_method: Charge calculation method (bcc, gas, resp)
            atom_type: Atom type (gaff, gaff2)
            residue_name: Residue name (3-letter code)
        
        Returns:
            Dict with parameterization results
        """
        logger.info(f"Generating GAFF parameters for {input_file}")
        
        input_path = Path(input_file)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Determine input format
        input_format = input_path.suffix[1:].lower()  # Remove dot
        if input_format == 'sdf':
            input_format = 'mdl'
        
        # Output files
        output_mol2 = output_dir / f"{residue_name}_gaff.mol2"
        output_frcmod = output_dir / f"{residue_name}.frcmod"
        
        # Run antechamber
        antechamber_args = [
            '-i', str(input_path),
            '-fi', input_format,
            '-o', str(output_mol2),
            '-fo', 'mol2',
            '-c', charge_method,
            '-nc', str(net_charge),
            '-at', atom_type,
            '-rn', residue_name,
            '-pf', 'y'  # Remove temporary files
        ]
        
        logger.info(f"Running antechamber with charge method: {charge_method}")
        
        try:
            result = self.run(antechamber_args, cwd=output_dir)
            logger.info("Antechamber completed successfully")
        except Exception as e:
            logger.error(f"Antechamber failed: {e}")
            raise
        
        # Run parmchk2 to check for missing parameters
        parmchk2_wrapper = BaseToolWrapper("parmchk2", conda_env=self.conda_env)
        
        parmchk2_args = [
            '-i', str(output_mol2),
            '-f', 'mol2',
            '-o', str(output_frcmod),
            '-s', atom_type
        ]
        
        logger.info("Running parmchk2")
        
        try:
            parmchk2_wrapper.run(parmchk2_args, cwd=output_dir)
            logger.info("Parmchk2 completed successfully")
        except Exception as e:
            logger.error(f"Parmchk2 failed: {e}")
            raise
        
        # Parse charges from MOL2
        charges = self._parse_charges_from_mol2(output_mol2)
        
        return {
            "mol2": str(output_mol2),
            "frcmod": str(output_frcmod),
            "charges": charges,
            "total_charge": sum(charges) if charges else 0.0,
            "residue_name": residue_name,
            "charge_method": charge_method,
            "atom_type": atom_type
        }
    
    def create_ligand_lib(
        self,
        mol2_file: Union[str, Path],
        frcmod_file: Union[str, Path],
        output_dir: Union[str, Path],
        residue_name: str = "LIG"
    ) -> Dict:
        """Create tleap library file for ligand
        
        Args:
            mol2_file: GAFF-typed MOL2 file
            frcmod_file: Force modification file
            output_dir: Output directory
            residue_name: Residue name
        
        Returns:
            Dict with library file paths
        """
        logger.info(f"Creating tleap library for {residue_name}")
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create tleap input script
        leap_in = output_dir / f"{residue_name}_lib.in"
        lib_file = output_dir / f"{residue_name}.lib"
        pdb_file = output_dir / f"{residue_name}.pdb"
        
        leap_script = f"""
source leaprc.gaff2
{residue_name} = loadmol2 {mol2_file}
loadamberparams {frcmod_file}
check {residue_name}
saveoff {residue_name} {lib_file}
savepdb {residue_name} {pdb_file}
quit
"""
        
        with open(leap_in, 'w') as f:
            f.write(leap_script)
        
        logger.info(f"Created tleap script: {leap_in}")
        
        # Run tleap
        tleap_wrapper = BaseToolWrapper("tleap", conda_env=self.conda_env)
        
        try:
            tleap_wrapper.run(['-f', str(leap_in)], cwd=output_dir)
            logger.info("tleap completed successfully")
        except Exception as e:
            logger.error(f"tleap failed: {e}")
            raise
        
        return {
            "lib": str(lib_file),
            "pdb": str(pdb_file),
            "leap_in": str(leap_in),
            "mol2": str(mol2_file),
            "frcmod": str(frcmod_file),
            "residue_name": residue_name
        }
    
    def _parse_charges_from_mol2(self, mol2_file: Union[str, Path]) -> List[float]:
        """Parse atomic charges from MOL2 file
        
        Args:
            mol2_file: MOL2 file path
        
        Returns:
            List of atomic charges
        """
        charges = []
        
        with open(mol2_file, 'r') as f:
            in_atom_section = False
            
            for line in f:
                if line.startswith('@<TRIPOS>ATOM'):
                    in_atom_section = True
                    continue
                elif line.startswith('@<TRIPOS>'):
                    in_atom_section = False
                
                if in_atom_section and line.strip():
                    parts = line.split()
                    if len(parts) >= 9:
                        try:
                            charge = float(parts[8])
                            charges.append(charge)
                        except ValueError:
                            continue
        
        return charges
    
    def build_system_tleap(
        self,
        protein_pdb: Optional[Union[str, Path]],
        ligand_lib: Optional[Union[str, Path]],
        output_dir: Union[str, Path],
        forcefield: str = "leaprc.protein.ff19SB",
        water_model: str = "tip3p",
        box_padding: float = 10.0,
        neutralize: bool = True,
        salt_conc: float = 0.15
    ) -> Dict:
        """Build MD system with tleap
        
        Args:
            protein_pdb: Protein PDB file (optional)
            ligand_lib: Ligand library file (optional)
            output_dir: Output directory
            forcefield: Protein force field
            water_model: Water model
            box_padding: Box padding in Angstroms
            neutralize: Neutralize system with ions
            salt_conc: Salt concentration (M)
        
        Returns:
            Dict with topology and coordinate files
        """
        logger.info("Building system with tleap")
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Output files
        prmtop = output_dir / "system.prmtop"
        inpcrd = output_dir / "system.inpcrd"
        leap_in = output_dir / "tleap.in"
        
        # Build tleap script
        script_lines = [
            f"source {forcefield}",
            f"source leaprc.water.{water_model}",
        ]
        
        # Load ligand if provided
        if ligand_lib:
            script_lines.append(f"loadoff {ligand_lib}")
        
        # Load protein if provided
        if protein_pdb:
            script_lines.append(f"protein = loadpdb {protein_pdb}")
            complex_name = "protein"
        
        # Create complex if both protein and ligand
        if protein_pdb and ligand_lib:
            ligand_name = Path(ligand_lib).stem
            script_lines.append(f"ligand = loadpdb {ligand_name}.pdb")
            script_lines.append("complex = combine {protein ligand}")
            complex_name = "complex"
        
        # Solvate
        water_box = "TIP3PBOX" if water_model == "tip3p" else f"{water_model.upper()}BOX"
        script_lines.append(f"solvatebox {complex_name} {water_box} {box_padding}")
        
        # Add ions
        if neutralize:
            script_lines.append(f"addionsrand {complex_name} Na+ 0")
            script_lines.append(f"addionsrand {complex_name} Cl- 0")
        
        if salt_conc > 0:
            script_lines.append(f"addionsrand {complex_name} Na+ {int(salt_conc * 100)}")
            script_lines.append(f"addionsrand {complex_name} Cl- {int(salt_conc * 100)}")
        
        # Save system
        script_lines.extend([
            f"check {complex_name}",
            f"saveamberparm {complex_name} {prmtop} {inpcrd}",
            "quit"
        ])
        
        # Write script
        with open(leap_in, 'w') as f:
            f.write('\n'.join(script_lines))
        
        logger.info(f"Created tleap script: {leap_in}")
        
        # Run tleap
        tleap_wrapper = BaseToolWrapper("tleap", conda_env=self.conda_env)
        
        try:
            tleap_wrapper.run(['-f', str(leap_in)], cwd=output_dir)
            logger.info("tleap system building completed")
        except Exception as e:
            logger.error(f"tleap failed: {e}")
            raise
        
        return {
            "prmtop": str(prmtop),
            "inpcrd": str(inpcrd),
            "leap_in": str(leap_in)
        }

