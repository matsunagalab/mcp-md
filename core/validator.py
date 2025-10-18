"""
Validator - QC checks and error detection.
"""

import logging
from pathlib import Path
from typing import Dict, List

logger = logging.getLogger(__name__)


class WorkflowValidator:
    """Validate workflow steps and outputs"""
    
    def __init__(self):
        logger.info("Workflow Validator initialized")
    
    def validate_step_output(self, step_result: Dict) -> bool:
        """Validate step output
        
        Args:
            step_result: Result dictionary from step execution
        
        Returns:
            True if valid
        """
        if "error" in step_result:
            logger.error(f"Step has error: {step_result['error']}")
            return False
        
        return True
    
    def validate_structure(self, pdb_file: str) -> Dict:
        """Validate structure file
        
        Args:
            pdb_file: PDB file path
        
        Returns:
            Validation result dict
        """
        if not Path(pdb_file).exists():
            return {"valid": False, "error": "File not found"}
        
        # Basic validation
        with open(pdb_file, 'r') as f:
            content = f.read()
            if not content.strip():
                return {"valid": False, "error": "Empty file"}
        
        return {"valid": True}

