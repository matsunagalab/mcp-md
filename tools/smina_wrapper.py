"""
smina wrapper for molecular docking.

smina is a fork of AutoDock Vina with additional scoring functions.
"""

import logging
from pathlib import Path
from typing import Optional, Union, List, Dict, Tuple
from .base_wrapper import BaseToolWrapper

logger = logging.getLogger(__name__)


class SminaWrapper(BaseToolWrapper):
    """Wrapper for smina docking"""
    
    def __init__(self, conda_env: str = "mcp-md-tools"):
        super().__init__("smina", conda_env=conda_env)
    
    def dock_ligand(
        self,
        receptor: Union[str, Path],
        ligand: Union[str, Path],
        center: Tuple[float, float, float],
        size: Tuple[float, float, float],
        output_dir: Union[str, Path],
        scoring: str = "vinardo",
        exhaustiveness: int = 8,
        num_modes: int = 9
    ) -> Dict:
        """Dock ligand to receptor
        
        Args:
            receptor: Receptor PDBQT file
            ligand: Ligand PDBQT file
            center: Center of search box (x, y, z)
            size: Size of search box (x, y, z)
            output_dir: Output directory
            scoring: Scoring function (vina, vinardo, ad4_scoring)
            exhaustiveness: Search exhaustiveness
            num_modes: Number of binding modes
        
        Returns:
            Dict with docking results
        """
        logger.info(f"Docking {ligand} to {receptor}")
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        output_pdbqt = output_dir / "docked_poses.pdbqt"
        log_file = output_dir / "docking.log"
        
        # Build command
        args = [
            '-r', str(receptor),
            '-l', str(ligand),
            '-o', str(output_pdbqt),
            '--center_x', str(center[0]),
            '--center_y', str(center[1]),
            '--center_z', str(center[2]),
            '--size_x', str(size[0]),
            '--size_y', str(size[1]),
            '--size_z', str(size[2]),
            '--scoring', scoring,
            '--exhaustiveness', str(exhaustiveness),
            '--num_modes', str(num_modes),
            '--log', str(log_file)
        ]
        
        # Run smina
        try:
            result = self.run(args, cwd=output_dir)
            logger.info("Docking completed successfully")
        except Exception as e:
            logger.error(f"smina docking failed: {e}")
            raise
        
        # Parse results
        poses, scores = self._parse_docking_results(output_pdbqt)
        
        return {
            "output_pdbqt": str(output_pdbqt),
            "poses": poses,
            "scores": scores,
            "num_poses": len(poses),
            "best_score": scores[0] if scores else None,
            "log_file": str(log_file),
            "scoring": scoring
        }
    
    def _parse_docking_results(self, pdbqt_file: Path) -> Tuple[List[str], List[float]]:
        """Parse docking output PDBQT
        
        Args:
            pdbqt_file: Output PDBQT file
        
        Returns:
            Tuple of (pose_files, scores)
        """
        poses = []
        scores = []
        
        if not pdbqt_file.exists():
            logger.warning(f"Output PDBQT not found: {pdbqt_file}")
            return poses, scores
        
        # Split PDBQT into individual poses
        current_pose = []
        current_score = None
        pose_idx = 0
        
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith('REMARK VINA RESULT:'):
                    # Extract score
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            current_score = float(parts[3])
                        except ValueError:
                            pass
                
                if line.startswith('MODEL'):
                    current_pose = [line]
                elif line.startswith('ENDMDL'):
                    current_pose.append(line)
                    
                    # Write individual pose file
                    pose_file = pdbqt_file.parent / f"pose_{pose_idx}.pdb"
                    with open(pose_file, 'w') as pf:
                        pf.write(''.join(current_pose))
                    
                    poses.append(str(pose_file))
                    if current_score is not None:
                        scores.append(current_score)
                    
                    current_pose = []
                    current_score = None
                    pose_idx += 1
                else:
                    if current_pose:
                        current_pose.append(line)
        
        return poses, scores

