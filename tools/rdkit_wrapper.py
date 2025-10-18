"""
RDKit wrapper for ligand 3D generation and manipulation.
"""

import logging
from pathlib import Path
from typing import Optional, Union
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

logger = logging.getLogger(__name__)


class RDKitWrapper:
    """Wrapper for RDKit ligand operations"""
    
    def __init__(self):
        logger.info("RDKit wrapper initialized")
    
    def smiles_to_3d(
        self,
        smiles: str,
        output_file: Union[str, Path],
        optimize: bool = True,
        num_confs: int = 1,
        method: str = "mmff"
    ) -> dict:
        """Generate 3D structure from SMILES
        
        Args:
            smiles: SMILES string
            output_file: Output file path (.mol2 or .sdf)
            optimize: Perform energy minimization
            num_confs: Number of conformers to generate
            method: Force field for optimization (mmff or uff)
        
        Returns:
            Dict with structure info
        """
        logger.info(f"Converting SMILES to 3D: {smiles}")
        
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates
        if num_confs > 1:
            # Multiple conformers
            conf_ids = AllChem.EmbedMultipleConfs(
                mol,
                numConfs=num_confs,
                randomSeed=42
            )
            if len(conf_ids) == 0:
                raise RuntimeError("Failed to generate conformers")
        else:
            # Single conformer
            result = AllChem.EmbedMolecule(mol, randomSeed=42)
            if result == -1:
                raise RuntimeError("Failed to generate 3D coordinates")
        
        # Optimize geometry
        if optimize:
            if method.lower() == "mmff":
                if num_confs > 1:
                    # Optimize all conformers
                    results = AllChem.MMFFOptimizeMoleculeConfs(mol)
                    energies = [r[1] for r in results]
                    best_conf = energies.index(min(energies))
                else:
                    AllChem.MMFFOptimizeMolecule(mol)
                    best_conf = 0
            elif method.lower() == "uff":
                if num_confs > 1:
                    results = AllChem.UFFOptimizeMoleculeConfs(mol)
                    energies = [r[1] for r in results]
                    best_conf = energies.index(min(energies))
                else:
                    AllChem.UFFOptimizeMolecule(mol)
                    best_conf = 0
            else:
                raise ValueError(f"Unknown method: {method}")
        else:
            best_conf = 0
        
        # Write output
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        if output_path.suffix.lower() == '.mol2':
            # Write MOL2
            writer = Chem.MolToMolBlock(mol, confId=best_conf)
            with open(output_path, 'w') as f:
                f.write(writer)
            logger.info(f"MOL2 file written (using MolBlock format)")
        elif output_path.suffix.lower() in ['.sdf', '.sd']:
            # Write SDF
            writer = Chem.SDWriter(str(output_path))
            writer.write(mol, confId=best_conf)
            writer.close()
        else:
            # Default to PDB
            Chem.MolToPDBFile(mol, str(output_path), confId=best_conf)
        
        logger.info(f"3D structure written to: {output_path}")
        
        # Get molecular properties
        mol_props = {
            "smiles": smiles,
            "output_file": str(output_path),
            "num_atoms": mol.GetNumAtoms(),
            "num_heavy_atoms": mol.GetNumHeavyAtoms(),
            "molecular_weight": Descriptors.MolWt(mol),
            "num_conformers": num_confs,
            "best_conformer": best_conf if num_confs > 1 else 0,
            "optimized": optimize,
            "method": method if optimize else None
        }
        
        return mol_props
    
    def mol2_to_pdb(
        self,
        mol2_file: Union[str, Path],
        pdb_file: Union[str, Path]
    ) -> dict:
        """Convert MOL2 to PDB
        
        Args:
            mol2_file: Input MOL2 file
            pdb_file: Output PDB file
        
        Returns:
            Conversion info
        """
        logger.info(f"Converting MOL2 to PDB: {mol2_file}")
        
        # Read MOL2 (try as MOL block first)
        with open(mol2_file, 'r') as f:
            mol_block = f.read()
        
        mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
        
        if mol is None:
            raise ValueError(f"Failed to read MOL2 file: {mol2_file}")
        
        # Write PDB
        Chem.MolToPDBFile(mol, str(pdb_file))
        
        logger.info(f"PDB written: {pdb_file}")
        
        return {
            "input": str(mol2_file),
            "output": str(pdb_file),
            "num_atoms": mol.GetNumAtoms()
        }
    
    def get_formal_charge(self, smiles: str) -> int:
        """Get formal charge from SMILES
        
        Args:
            smiles: SMILES string
        
        Returns:
            Total formal charge
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        return Chem.GetFormalCharge(mol)
    
    def validate_smiles(self, smiles: str) -> bool:
        """Validate SMILES string
        
        Args:
            smiles: SMILES string
        
        Returns:
            True if valid
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except:
            return False

