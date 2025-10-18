"""
Common data models for MCP-MD workflow.

Defines Pydantic models for various molecular dynamics components.
"""

from typing import Optional, Any
from enum import Enum
from pathlib import Path
from pydantic import BaseModel, Field


class ForceField(str, Enum):
    """Supported force fields"""
    FF19SB = "ff19SB"
    FF14SB = "ff14SB"
    OL15 = "OL15"
    OL3 = "OL3"
    LIPID17 = "lipid17"
    CHARMM36 = "CHARMM36"
    GLYCAM06 = "GLYCAM06"
    GAFF2 = "gaff2"


class WaterModel(str, Enum):
    """Supported water models"""
    TIP3P = "tip3p"
    TIP4PEW = "tip4pew"
    OPC = "opc"
    SPCE = "spce"


class ChargeMethod(str, Enum):
    """Ligand charge calculation methods"""
    AM1_BCC = "bcc"
    GASTEIGER = "gas"
    RESP = "resp"


class SystemType(str, Enum):
    """MD system types"""
    PROTEIN_ONLY = "protein_only"
    PROTEIN_LIGAND_COMPLEX = "protein_ligand_complex"
    MEMBRANE = "membrane_system"
    NUCLEIC = "nucleic_system"
    BOLTZ2_DENOVO = "boltz2_denovo"
    BOLTZ2_COMPLEX = "boltz2_complex"


class StructureInfo(BaseModel):
    """Structure information"""
    file_path: str = Field(..., description="Path to structure file")
    format: str = Field(..., description="File format (pdb, mol2, sdf, etc.)")
    num_atoms: Optional[int] = Field(None, description="Number of atoms")
    num_residues: Optional[int] = Field(None, description="Number of residues")
    chains: Optional[list[str]] = Field(None, description="Chain IDs")
    
    class Config:
        json_schema_extra = {
            "example": {
                "file_path": "/path/to/protein.pdb",
                "format": "pdb",
                "num_atoms": 1234,
                "num_residues": 150,
                "chains": ["A"]
            }
        }


class Boltz2Prediction(BaseModel):
    """Boltz-2 prediction results"""
    structures: list[str] = Field(..., description="Paths to predicted structures")
    confidence: dict[str, Any] = Field(..., description="pLDDT and pAE scores")
    yaml_input: Optional[str] = Field(None, description="Input YAML path")
    affinity: Optional[dict[str, float]] = Field(None, description="Affinity prediction")
    
    class Config:
        json_schema_extra = {
            "example": {
                "structures": ["model_0.pdb", "model_1.pdb"],
                "confidence": {
                    "plddt": [85.2, 83.1],
                    "pae": [[0.5, 1.2], [1.2, 0.4]]
                },
                "affinity": {
                    "probability_binary": 0.85,
                    "pred_value": -6.2,
                    "ic50_um": 0.63
                }
            }
        }


class LigandParams(BaseModel):
    """Ligand parameterization results"""
    mol2: str = Field(..., description="Path to GAFF-typed MOL2 file")
    frcmod: str = Field(..., description="Path to force modification file")
    charges: list[float] = Field(..., description="Atomic partial charges")
    total_charge: float = Field(..., description="Total charge")
    residue_name: str = Field(default="LIG", description="Residue name")
    
    class Config:
        json_schema_extra = {
            "example": {
                "mol2": "ligand_gaff.mol2",
                "frcmod": "ligand.frcmod",
                "charges": [-0.5, 0.1, 0.2],
                "total_charge": 0.0012,
                "residue_name": "LIG"
            }
        }


class DockingResult(BaseModel):
    """Docking results"""
    poses: list[str] = Field(..., description="Paths to docked poses")
    scores: list[float] = Field(..., description="Docking scores")
    best_pose_idx: int = Field(..., description="Index of best pose")
    
    class Config:
        json_schema_extra = {
            "example": {
                "poses": ["pose_0.pdb", "pose_1.pdb"],
                "scores": [-8.5, -7.2],
                "best_pose_idx": 0
            }
        }


class AmberSystem(BaseModel):
    """Amber topology and coordinate files"""
    prmtop: str = Field(..., description="Path to topology file (.prmtop)")
    inpcrd: str = Field(..., description="Path to coordinate file (.inpcrd or .rst7)")
    leap_in: Optional[str] = Field(None, description="Path to tleap input script")
    
    class Config:
        json_schema_extra = {
            "example": {
                "prmtop": "system.prmtop",
                "inpcrd": "system.rst7",
                "leap_in": "tleap.in"
            }
        }


class OpenMMWorkflow(BaseModel):
    """OpenMM MD workflow scripts"""
    minimize_script: str = Field(..., description="Minimization script")
    equilibrate_script: str = Field(..., description="Equilibration script")
    production_script: str = Field(..., description="Production MD script")
    submit_script: Optional[str] = Field(None, description="Job submission script")
    
    class Config:
        json_schema_extra = {
            "example": {
                "minimize_script": "1_minimize.py",
                "equilibrate_script": "2_equilibrate.py",
                "production_script": "3_production.py",
                "submit_script": "submit.sh"
            }
        }


class WorkflowStep(BaseModel):
    """Single workflow step"""
    step_id: str = Field(..., description="Unique step identifier")
    name: str = Field(..., description="Step name")
    mcp_tool: str = Field(..., description="MCP tool function name")
    mcp_server: str = Field(..., description="MCP server name")
    params: dict[str, Any] = Field(default_factory=dict, description="Step parameters")
    inputs: list[str] = Field(default_factory=list, description="Input file paths")
    outputs: list[str] = Field(default_factory=list, description="Output file paths")
    dependencies: list[str] = Field(default_factory=list, description="Dependent step IDs")
    optional: bool = Field(default=False, description="Whether step is optional")
    status: str = Field(default="pending", description="Step status")
    
    class Config:
        json_schema_extra = {
            "example": {
                "step_id": "step_001",
                "name": "fetch_pdb",
                "mcp_tool": "fetch_pdb",
                "mcp_server": "structure_server",
                "params": {"pdb_id": "1ABC"},
                "outputs": ["1ABC.pdb"],
                "status": "pending"
            }
        }


class WorkflowPlan(BaseModel):
    """Complete workflow plan"""
    plan_id: str = Field(..., description="Unique plan identifier")
    system_type: SystemType = Field(..., description="System type")
    steps: list[WorkflowStep] = Field(..., description="Workflow steps")
    description: Optional[str] = Field(None, description="Plan description")
    created_at: str = Field(..., description="Creation timestamp")
    
    class Config:
        json_schema_extra = {
            "example": {
                "plan_id": "plan_20250118_001",
                "system_type": "protein_ligand_complex",
                "steps": [],
                "description": "Build protein-ligand complex for MD",
                "created_at": "2025-01-18T10:30:00"
            }
        }

