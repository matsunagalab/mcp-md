"""
PDB2PQR wrapper for pH-dependent protonation.

Uses PDB2PQR and PROPKA for accurate protonation state assignment.
"""

import logging
from pathlib import Path
from typing import Optional, Union, Dict
from .base_wrapper import BaseToolWrapper

logger = logging.getLogger(__name__)


class PDB2PQRWrapper(BaseToolWrapper):
    """Wrapper for PDB2PQR protonation"""
    
    def __init__(self):
        super().__init__("pdb2pqr", conda_env=None)
    
    def protonate_structure(
        self,
        input_pdb: Union[str, Path],
        output_pdb: Union[str, Path],
        ph: float = 7.0,
        forcefield: str = "AMBER"
    ) -> Dict:
        """Protonate structure at specified pH
        
        Args:
            input_pdb: Input PDB file
            output_pdb: Output PDB file
            ph: pH value (default 7.0)
            forcefield: Force field (AMBER, CHARMM, PARSE)
        
        Returns:
            Dict with protonation results
        """
        logger.info(f"Protonating structure at pH {ph}")
        
        output_path = Path(output_pdb)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # PDB2PQR arguments
        args = [
            '--ff', forcefield.lower(),
            '--with-ph', str(ph),
            '--titration-state-method', 'propka',
            '--drop-water',
            str(input_pdb),
            str(output_pdb)
        ]
        
        try:
            result = self.run(args)
            logger.info("Protonation completed successfully")
        except Exception as e:
            logger.error(f"PDB2PQR failed: {e}")
            raise
        
        return {
            "input": str(input_pdb),
            "output": str(output_pdb),
            "ph": ph,
            "forcefield": forcefield,
            "success": True
        }
    
    def get_pka_values(
        self,
        input_pdb: Union[str, Path],
        output_dir: Union[str, Path]
    ) -> Dict:
        """Calculate pKa values using PROPKA
        
        Args:
            input_pdb: Input PDB file
            output_dir: Output directory
        
        Returns:
            Dict with pKa analysis
        """
        logger.info("Calculating pKa values with PROPKA")
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Run PROPKA through PDB2PQR
        output_pqr = output_dir / "structure.pqr"
        
        args = [
            '--ff', 'amber',
            '--titration-state-method', 'propka',
            '--pka-output', str(output_dir / "pka.out"),
            str(input_pdb),
            str(output_pqr)
        ]
        
        try:
            result = self.run(args, cwd=output_dir)
            
            # Parse pKa output
            pka_results = self._parse_pka_output(output_dir / "pka.out")
            
            return {
                "input": str(input_pdb),
                "pka_file": str(output_dir / "pka.out"),
                "pka_values": pka_results,
                "success": True
            }
        
        except Exception as e:
            logger.error(f"pKa calculation failed: {e}")
            raise
    
    def _parse_pka_output(self, pka_file: Path) -> Dict:
        """Parse PROPKA pKa output
        
        Args:
            pka_file: PROPKA output file
        
        Returns:
            Parsed pKa values
        """
        pka_values = {}
        
        if not pka_file.exists():
            return pka_values
        
        with open(pka_file, 'r') as f:
            for line in f:
                # Parse pKa lines
                # Format: residue chain_id resnum pKa
                if line.strip() and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            residue = parts[0]
                            chain = parts[1]
                            resnum = parts[2]
                            pka = float(parts[3])
                            
                            key = f"{residue}_{chain}_{resnum}"
                            pka_values[key] = pka
                        except (ValueError, IndexError):
                            continue
        
        return pka_values

