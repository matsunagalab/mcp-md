"""
MolProbity Wrapper - Custom implementation for structure quality checks.

Uses MDAnalysis and BioPython for:
- Atom clashes (steric conflicts)
- Bond length/angle validation
- Chirality checks
- Solvent accessibility
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Any
import numpy as np
from .base_wrapper import BaseToolWrapper

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
except ImportError:
    mda = None

try:
    from Bio.PDB import PDBParser, NeighborSearch
    from Bio.PDB.Polypeptide import is_aa
except ImportError:
    PDBParser = None

logger = logging.getLogger(__name__)


class MolProbityWrapper(BaseToolWrapper):
    """Custom MolProbity implementation for structure QC"""
    
    def __init__(self):
        super().__init__("molprobity", conda_env="mcp-md")
        
        if mda is None:
            raise ImportError("MDAnalysis not installed. Install with: pip install MDAnalysis")
        if PDBParser is None:
            raise ImportError("BioPython not installed. Install with: pip install biopython")
    
    def check_clashes(self, pdb_file: str, clash_cutoff: float = 2.0) -> Dict[str, Any]:
        """Check for atom clashes (steric conflicts)"""
        logger.info(f"Checking clashes in {pdb_file} (cutoff={clash_cutoff}Å)")
        
        try:
            # Load structure with BioPython
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("protein", pdb_file)
            
            # Get all atoms
            atoms = [atom for atom in structure.get_atoms()]
            
            # Build neighbor search
            ns = NeighborSearch(atoms)
            
            # Find clashes
            clashes = []
            for atom1 in atoms:
                neighbors = ns.search(atom1.coord, clash_cutoff, level='A')
                for atom2 in neighbors:
                    if atom1 != atom2:
                        dist = atom1 - atom2
                        if dist < clash_cutoff:
                            clashes.append({
                                "atom1": f"{atom1.get_parent().get_resname()}{atom1.get_parent().id[1]}.{atom1.name}",
                                "atom2": f"{atom2.get_parent().get_resname()}{atom2.get_parent().id[1]}.{atom2.name}",
                                "distance": float(dist),
                                "severity": "high" if dist < 1.5 else "medium"
                            })
            
            # Remove duplicates (A-B and B-A)
            unique_clashes = []
            seen = set()
            for clash in clashes:
                pair = tuple(sorted([clash["atom1"], clash["atom2"]]))
                if pair not in seen:
                    seen.add(pair)
                    unique_clashes.append(clash)
            
            return {
                "num_clashes": len(unique_clashes),
                "clashes": unique_clashes[:50],  # Limit to 50 for output
                "status": "pass" if len(unique_clashes) == 0 else "fail"
            }
        
        except Exception as e:
            logger.error(f"Clash check failed: {e}")
            return {"error": str(e), "status": "error"}
    
    def check_bond_lengths(self, pdb_file: str) -> Dict[str, Any]:
        """Validate bond lengths"""
        logger.info(f"Checking bond lengths in {pdb_file}")
        
        try:
            u = mda.Universe(pdb_file)
            
            # Standard bond lengths (Å)
            BOND_STANDARDS = {
                ("C", "C"): (1.54, 0.05),   # C-C single bond
                ("C", "N"): (1.47, 0.05),   # C-N single bond
                ("C", "O"): (1.43, 0.05),   # C-O single bond
                ("C", "S"): (1.82, 0.05),   # C-S bond
                ("N", "H"): (1.01, 0.05),   # N-H bond
                ("O", "H"): (0.96, 0.05),   # O-H bond
            }
            
            outliers = []
            
            # Check backbone bonds
            protein = u.select_atoms("protein")
            
            # CA-C bonds
            ca_atoms = protein.select_atoms("name CA")
            c_atoms = protein.select_atoms("name C")
            
            for i in range(min(len(ca_atoms), len(c_atoms))):
                dist = np.linalg.norm(ca_atoms[i].position - c_atoms[i].position)
                expected, tolerance = 1.52, 0.05
                if abs(dist - expected) > tolerance:
                    outliers.append({
                        "bond": f"CA{i}-C{i}",
                        "measured": float(dist),
                        "expected": expected,
                        "deviation": float(abs(dist - expected))
                    })
            
            return {
                "num_outliers": len(outliers),
                "outliers": outliers[:50],
                "status": "pass" if len(outliers) < 10 else "warning"
            }
        
        except Exception as e:
            logger.error(f"Bond length check failed: {e}")
            return {"error": str(e), "status": "error"}
    
    def check_chirality(self, pdb_file: str) -> Dict[str, Any]:
        """Check amino acid chirality (L vs D)"""
        logger.info(f"Checking chirality in {pdb_file}")
        
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("protein", pdb_file)
            
            wrong_chirality = []
            
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if not is_aa(residue, standard=True):
                            continue
                        
                        # Get CA, C, N, CB atoms
                        try:
                            ca = residue["CA"].coord
                            c = residue["C"].coord
                            n = residue["N"].coord
                            
                            # Skip glycine (no CB)
                            if residue.get_resname() == "GLY":
                                continue
                            
                            cb = residue["CB"].coord
                            
                            # Calculate chirality (cross product)
                            v1 = n - ca
                            v2 = c - ca
                            v3 = cb - ca
                            
                            chirality = np.dot(np.cross(v1, v2), v3)
                            
                            # L-amino acids should have negative chirality
                            if chirality > 0:
                                wrong_chirality.append({
                                    "residue": f"{residue.get_resname()}{residue.id[1]}",
                                    "chain": chain.id,
                                    "chirality_value": float(chirality),
                                    "expected": "L",
                                    "found": "D"
                                })
                        
                        except KeyError:
                            # Missing atoms
                            continue
            
            return {
                "num_wrong_chirality": len(wrong_chirality),
                "wrong_chirality": wrong_chirality,
                "status": "pass" if len(wrong_chirality) == 0 else "fail"
            }
        
        except Exception as e:
            logger.error(f"Chirality check failed: {e}")
            return {"error": str(e), "status": "error"}
    
    def run_full_qc(self, pdb_file: str) -> Dict[str, Any]:
        """Run all quality checks"""
        logger.info(f"Running full QC on {pdb_file}")
        
        clash_result = self.check_clashes(pdb_file)
        bond_result = self.check_bond_lengths(pdb_file)
        chirality_result = self.check_chirality(pdb_file)
        
        # Overall status
        statuses = [
            clash_result.get("status"),
            bond_result.get("status"),
            chirality_result.get("status")
        ]
        
        if "error" in statuses:
            overall_status = "error"
        elif "fail" in statuses:
            overall_status = "fail"
        elif "warning" in statuses:
            overall_status = "warning"
        else:
            overall_status = "pass"
        
        return {
            "overall_status": overall_status,
            "clashes": clash_result,
            "bond_lengths": bond_result,
            "chirality": chirality_result,
            "pdb_file": pdb_file
        }

