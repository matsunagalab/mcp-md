"""
Boltz-2 wrapper for structure and affinity prediction.

Provides high-level interface to Boltz-2 for:
- Structure prediction from FASTA
- Protein-ligand complex prediction from FASTA + SMILES
- Binding affinity prediction
- Missing residue completion
"""

import json
import logging
import yaml
from pathlib import Path
from typing import Optional, Union, List, Dict, Any
from .base_wrapper import BaseToolWrapper

logger = logging.getLogger(__name__)


class Boltz2Wrapper(BaseToolWrapper):
    """Wrapper for Boltz-2 biomolecular prediction"""
    
    def __init__(self):
        super().__init__("boltz", conda_env=None)
    
    def predict_structure(
        self,
        sequences: List[Dict[str, Any]],
        output_dir: Union[str, Path],
        use_msa: bool = True,
        num_models: int = 5
    ) -> Dict[str, Any]:
        """Predict structure from sequences
        
        Args:
            sequences: List of sequence dicts with type and sequence
            output_dir: Output directory
            use_msa: Use MSA server (requires internet)
            num_models: Number of models to generate
        
        Returns:
            Dict with predicted structures and confidence scores
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create YAML input
        yaml_input = self._create_yaml_input(sequences)
        yaml_path = output_dir / "boltz_input.yaml"
        
        with open(yaml_path, 'w') as f:
            yaml.dump(yaml_input, f)
        
        logger.info(f"Created Boltz-2 input YAML: {yaml_path}")
        
        # Run Boltz-2
        args = ["predict", str(yaml_path)]
        if use_msa:
            args.append("--use_msa_server")
        
        try:
            self.run(args, cwd=output_dir)
        except Exception as e:
            logger.error(f"Boltz-2 prediction failed: {e}")
            raise
        
        # Parse results
        results = self._parse_results(output_dir)
        results["yaml_input"] = str(yaml_path)
        
        return results
    
    def predict_complex_with_affinity(
        self,
        protein_fasta: str,
        ligand_smiles: List[str],
        output_dir: Union[str, Path],
        use_msa: bool = True,
        num_models: int = 5
    ) -> Dict[str, Any]:
        """Predict protein-ligand complex with affinity
        
        Args:
            protein_fasta: Protein FASTA sequence
            ligand_smiles: List of ligand SMILES strings
            output_dir: Output directory
            use_msa: Use MSA server
            num_models: Number of models
        
        Returns:
            Dict with structures, affinity predictions
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create sequences list
        sequences = [
            {"protein": {"id": "protein_A", "sequence": protein_fasta}}
        ]
        
        for i, smiles in enumerate(ligand_smiles):
            sequences.append({
                "ligand": {"id": f"ligand_{i}", "smiles": smiles}
            })
        
        # Create YAML with affinity prediction
        yaml_input = {
            "sequences": sequences,
            "affinity": {
                "enabled": True,
                "target_chain": "protein_A",
                "ligand_chain": "ligand_0"
            }
        }
        
        yaml_path = output_dir / "complex_affinity.yaml"
        with open(yaml_path, 'w') as f:
            yaml.dump(yaml_input, f)
        
        logger.info(f"Created complex+affinity YAML: {yaml_path}")
        
        # Run Boltz-2
        args = ["predict", str(yaml_path)]
        if use_msa:
            args.append("--use_msa_server")
        
        try:
            self.run(args, cwd=output_dir)
        except Exception as e:
            logger.error(f"Boltz-2 complex prediction failed: {e}")
            raise
        
        # Parse results
        results = self._parse_results(output_dir)
        results["yaml_input"] = str(yaml_path)
        
        # Parse affinity if available
        affinity_json = output_dir / "affinity.json"
        if affinity_json.exists():
            with open(affinity_json, 'r') as f:
                affinity_data = json.load(f)
            results["affinity"] = self._parse_affinity(affinity_data)
        
        return results
    
    def screen_ligands(
        self,
        protein_fasta: str,
        ligand_smiles_list: List[str],
        output_dir: Union[str, Path],
        screening_mode: str = "binary"
    ) -> Dict[str, Any]:
        """Screen multiple ligands for binding
        
        Args:
            protein_fasta: Protein sequence
            ligand_smiles_list: List of SMILES to screen
            output_dir: Output directory
            screening_mode: "binary" (hit discovery) or "quantitative" (optimization)
        
        Returns:
            Ranked screening results
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        all_results = []
        
        for i, smiles in enumerate(ligand_smiles_list):
            logger.info(f"Screening ligand {i+1}/{len(ligand_smiles_list)}: {smiles}")
            
            ligand_dir = output_dir / f"ligand_{i}"
            result = self.predict_complex_with_affinity(
                protein_fasta=protein_fasta,
                ligand_smiles=[smiles],
                output_dir=ligand_dir,
                use_msa=False,  # Skip MSA for screening speed
                num_models=1
            )
            
            # Extract relevant affinity score
            if "affinity" in result:
                if screening_mode == "binary":
                    score = result["affinity"].get("probability_binary", 0.0)
                else:
                    score = result["affinity"].get("pred_value", 0.0)
                
                all_results.append({
                    "smiles": smiles,
                    "affinity_score": score,
                    "structure": result["structures"][0] if result["structures"] else None,
                    "ligand_id": i
                })
        
        # Rank results
        reverse = (screening_mode == "binary")  # Higher is better for binary
        ranked = sorted(all_results, key=lambda x: x["affinity_score"], reverse=reverse)
        
        for rank, result in enumerate(ranked, 1):
            result["rank"] = rank
        
        return {
            "results": ranked,
            "ranked_by": "affinity_probability_binary" if screening_mode == "binary" else "affinity_pred_value",
            "num_ligands": len(ligand_smiles_list)
        }
    
    def _create_yaml_input(self, sequences: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Create Boltz-2 YAML input
        
        Args:
            sequences: List of sequence dictionaries
        
        Returns:
            YAML-compatible dict
        """
        return {"sequences": sequences}
    
    def _parse_results(self, output_dir: Path) -> Dict[str, Any]:
        """Parse Boltz-2 output files
        
        Args:
            output_dir: Output directory
        
        Returns:
            Parsed results dict
        """
        results = {
            "structures": [],
            "confidence": {}
        }
        
        # Find PDB structures
        pdb_files = sorted(output_dir.glob("*.pdb"))
        results["structures"] = [str(f) for f in pdb_files]
        
        # Parse confidence scores if available
        confidence_json = output_dir / "confidence.json"
        if confidence_json.exists():
            with open(confidence_json, 'r') as f:
                results["confidence"] = json.load(f)
        
        return results
    
    def _parse_affinity(self, affinity_data: Dict) -> Dict[str, float]:
        """Parse affinity prediction data
        
        Args:
            affinity_data: Raw affinity JSON data
        
        Returns:
            Parsed affinity dict with user-friendly keys
        """
        parsed = {}
        
        # Extract key metrics
        if "affinity_probability_binary" in affinity_data:
            parsed["probability_binary"] = affinity_data["affinity_probability_binary"]
        
        if "affinity_pred_value" in affinity_data:
            pred_value = affinity_data["affinity_pred_value"]
            parsed["pred_value"] = pred_value
            
            # Convert log10(IC50) to IC50 in μM
            parsed["ic50_um"] = 10 ** pred_value
        
        return parsed

