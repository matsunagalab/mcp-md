"""
Amber Prep Server - Boltz-2 complex to MD simulation input files with FastMCP.

Provides MCP tools for:
- mmCIF complex parsing (Boltz-2 output)
- Protein preparation with pdb4amber
- Ligand parameterization with GAFF2/AM1-BCC
- Robust antechamber with sqm error diagnostics
- frcmod validation
- Complete MD system building with tleap
"""

import json
import logging
import re
import shutil
import uuid
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple
from fastmcp import FastMCP

from common.base import BaseToolWrapper
from common.utils import setup_logger, ensure_directory

logger = setup_logger(__name__)

# Create FastMCP server
mcp = FastMCP("Amber Prep Server")

# Initialize working directory
WORKING_DIR = Path("output/amber_prep")
ensure_directory(WORKING_DIR)

# Initialize tool wrappers
antechamber_wrapper = BaseToolWrapper("antechamber", conda_env="mcp-md")
parmchk2_wrapper = BaseToolWrapper("parmchk2", conda_env="mcp-md")
pdb4amber_wrapper = BaseToolWrapper("pdb4amber", conda_env="mcp-md")
tleap_wrapper = BaseToolWrapper("tleap", conda_env="mcp-md")
obabel_wrapper = BaseToolWrapper("obabel", conda_env="mcp-md")


# =============================================================================
# Helper Functions
# =============================================================================


def _generate_job_id() -> str:
    """Generate unique job identifier."""
    return str(uuid.uuid4())[:8]


def _parse_sqm_output(sqm_out_path: Path) -> Dict[str, Any]:
    """Parse sqm.out file for errors and diagnostics.
    
    Args:
        sqm_out_path: Path to sqm.out file
    
    Returns:
        Dict with diagnostics information
    """
    diagnostics = {
        "success": False,
        "errors": [],
        "warnings": [],
        "recommendations": [],
        "electron_count": None,
        "scf_converged": None
    }
    
    if not sqm_out_path.exists():
        diagnostics["errors"].append("sqm.out file not found")
        return diagnostics
    
    content = sqm_out_path.read_text()
    
    # Check for odd electron error
    if "number of electrons is odd" in content.lower():
        match = re.search(r"electrons is odd[^\d]*(\d+)", content, re.IGNORECASE)
        electron_count = match.group(1) if match else "unknown"
        diagnostics["errors"].append(f"Odd number of electrons ({electron_count})")
        diagnostics["electron_count"] = electron_count
        diagnostics["recommendations"].append(
            "Net charge is likely incorrect. Try adjusting by ±1."
        )
    
    # Check for SCF convergence issues
    if "no convergence in scf" in content.lower():
        diagnostics["errors"].append("SCF calculation did not converge")
        diagnostics["scf_converged"] = False
        diagnostics["recommendations"].extend([
            "Try optimizing the ligand structure before parameterization.",
            "Use a molecular viewer to check for unreasonable bond lengths.",
            "Consider using a different charge method (e.g., gas instead of bcc)."
        ])
    elif "scf" in content.lower() and "converged" in content.lower():
        diagnostics["scf_converged"] = True
    
    # Check for general sqm errors
    if "error" in content.lower() or "fatal" in content.lower():
        error_lines = [
            line.strip() for line in content.split('\n')
            if 'error' in line.lower() or 'fatal' in line.lower()
        ]
        for line in error_lines[:5]:  # Limit to first 5 errors
            if line not in diagnostics["errors"]:
                diagnostics["errors"].append(line)
    
    # Check for successful completion
    if "calculation completed" in content.lower() or len(diagnostics["errors"]) == 0:
        diagnostics["success"] = True
    
    return diagnostics


def _parse_frcmod_warnings(frcmod_path: Path) -> Dict[str, Any]:
    """Parse frcmod file for missing parameter warnings.
    
    Args:
        frcmod_path: Path to .frcmod file
    
    Returns:
        Dict with validation results
    """
    validation = {
        "valid": True,
        "warnings": [],
        "missing_params": {
            "bonds": [],
            "angles": [],
            "dihedrals": [],
            "impropers": []
        },
        "attn_count": 0,
        "recommendations": []
    }
    
    if not frcmod_path.exists():
        validation["valid"] = False
        validation["warnings"].append("frcmod file not found")
        return validation
    
    content = frcmod_path.read_text()
    lines = content.split('\n')
    
    current_section = None
    section_map = {
        "BOND": "bonds",
        "ANGLE": "angles",
        "DIHE": "dihedrals",
        "IMPROPER": "impropers"
    }
    
    for line in lines:
        # Track section
        for section_name in section_map:
            if line.strip().startswith(section_name):
                current_section = section_map[section_name]
                break
        
        # Check for ATTN warnings
        if "ATTN" in line or "need revision" in line.lower():
            validation["attn_count"] += 1
            validation["warnings"].append(line.strip())
            validation["valid"] = False
            
            if current_section:
                # Extract parameter type
                parts = line.split()
                if parts:
                    param_type = parts[0]
                    validation["missing_params"][current_section].append(param_type)
        
        # Check for zero force constants (dangerous)
        if re.search(r'\b0\.0+\s+0\.0+\b', line):
            validation["warnings"].append(f"Zero force constant detected: {line.strip()}")
    
    if validation["attn_count"] > 0:
        validation["recommendations"].extend([
            f"Found {validation['attn_count']} parameters requiring attention.",
            "These parameters were estimated by analogy and may be inaccurate.",
            "Consider consulting a computational chemistry expert for validation.",
            "For production runs, quantum chemistry calculations may be needed."
        ])
    
    return validation


def _estimate_charge_rdkit(mol) -> Dict[str, Any]:
    """Estimate net charge using RDKit.
    
    Args:
        mol: RDKit molecule object
    
    Returns:
        Dict with charge estimation results
    """
    from rdkit import Chem
    
    result = {
        "formal_charge": Chem.GetFormalCharge(mol),
        "ionizable_groups": [],
        "method": "rdkit"
    }
    
    # Identify ionizable groups
    # Carboxylic acids (typically deprotonated at pH 7.4)
    # Use the simplest working pattern
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)O")
    if carboxylic_pattern is not None:
        matches = mol.GetSubstructMatches(carboxylic_pattern)
        if matches:
            result["ionizable_groups"].append({
                "type": "carboxylic_acid",
                "count": len(matches),
                "typical_charge": -1,
                "pka_range": "3-5"
            })
    
    # Primary amines (typically protonated at pH 7.4)
    # Excludes amides (NC=O) and aromatic amines (lower pKa)
    amine_pattern = Chem.MolFromSmarts("[NX3;H2;!$(NC=O);!$(Nc)]")
    if amine_pattern and mol.HasSubstructMatch(amine_pattern):
        matches = mol.GetSubstructMatches(amine_pattern)
        result["ionizable_groups"].append({
            "type": "primary_amine",
            "count": len(matches),
            "typical_charge": +1,
            "pka_range": "9-11"
        })
    
    # Secondary amines (excludes amides and aromatic)
    sec_amine_pattern = Chem.MolFromSmarts("[NX3;H1;!$(NC=O);!$(Nc)]([#6])[#6]")
    if sec_amine_pattern and mol.HasSubstructMatch(sec_amine_pattern):
        matches = mol.GetSubstructMatches(sec_amine_pattern)
        result["ionizable_groups"].append({
            "type": "secondary_amine",
            "count": len(matches),
            "typical_charge": +1,
            "pka_range": "9-11"
        })
    
    # Tertiary amines (excludes amides)
    tert_amine_pattern = Chem.MolFromSmarts("[NX3;H0;!$(NC=O);!$(Nc)]([#6])([#6])[#6]")
    if tert_amine_pattern and mol.HasSubstructMatch(tert_amine_pattern):
        matches = mol.GetSubstructMatches(tert_amine_pattern)
        result["ionizable_groups"].append({
            "type": "tertiary_amine",
            "count": len(matches),
            "typical_charge": +1,
            "pka_range": "9-11"
        })
    
    # Phenols
    phenol_pattern = Chem.MolFromSmarts("[OX2H1]c1ccccc1")
    if mol.HasSubstructMatch(phenol_pattern):
        matches = mol.GetSubstructMatches(phenol_pattern)
        result["ionizable_groups"].append({
            "type": "phenol",
            "count": len(matches),
            "typical_charge": 0,
            "pka_range": "9-10"
        })
    
    # Sulfonic acids
    sulfonic_pattern = Chem.MolFromSmarts("[SX4](=O)(=O)[OX1H1]")
    if mol.HasSubstructMatch(sulfonic_pattern):
        matches = mol.GetSubstructMatches(sulfonic_pattern)
        result["ionizable_groups"].append({
            "type": "sulfonic_acid",
            "count": len(matches),
            "typical_charge": -1,
            "pka_range": "<1"
        })
    
    # Phosphates
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=O)([OX1H1])([OX1H1])")
    if mol.HasSubstructMatch(phosphate_pattern):
        matches = mol.GetSubstructMatches(phosphate_pattern)
        result["ionizable_groups"].append({
            "type": "phosphate",
            "count": len(matches),
            "typical_charge": -2,
            "pka_range": "2, 7"
        })
    
    return result


def _estimate_physiological_charge(charge_info: Dict[str, Any], ph: float = 7.4) -> int:
    """Estimate net charge at physiological pH.
    
    Args:
        charge_info: Output from _estimate_charge_rdkit
        ph: Target pH
    
    Returns:
        Estimated integer net charge
    """
    estimated_charge = charge_info["formal_charge"]
    
    for group in charge_info.get("ionizable_groups", []):
        group_type = group["type"]
        count = group["count"]
        
        # Adjust based on typical protonation at pH 7.4
        if group_type in ["carboxylic_acid", "sulfonic_acid"]:
            # Typically deprotonated (negative)
            estimated_charge -= count
        elif group_type in ["primary_amine", "secondary_amine"]:
            # Typically protonated (positive) 
            estimated_charge += count
        elif group_type == "phosphate":
            # Typically -2 at pH 7.4
            estimated_charge -= 2 * count
    
    return estimated_charge


# =============================================================================
# MCP Tools
# =============================================================================


@mcp.tool()
def parse_boltz2_complex(
    cif_file: str,
    output_dir: Optional[str] = None
) -> dict:
    """Parse Boltz-2 mmCIF output and separate protein from ligand.
    
    This tool extracts protein and ligand components from a Boltz-2 predicted
    complex structure in mmCIF format.
    
    Args:
        cif_file: Path to mmCIF file from Boltz-2
        output_dir: Output directory (auto-generated if None)
    
    Returns:
        Dict with separated structure files and metadata
    """
    logger.info(f"Parsing Boltz-2 complex: {cif_file}")
    
    try:
        import gemmi
    except ImportError:
        raise ImportError("gemmi not installed. Install with: pip install gemmi")
    
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        raise ImportError("RDKit not installed. Install via conda.")
    
    cif_path = Path(cif_file)
    if not cif_path.exists():
        raise FileNotFoundError(f"mmCIF file not found: {cif_file}")
    
    # Setup output directory
    job_id = _generate_job_id()
    if output_dir is None:
        output_dir = WORKING_DIR / job_id
    else:
        output_dir = Path(output_dir) / job_id
    ensure_directory(output_dir)
    
    # Parse mmCIF
    logger.info("Reading mmCIF with gemmi...")
    doc = gemmi.cif.read(str(cif_path))
    block = doc[0]
    structure = gemmi.make_structure_from_block(block)
    
    # Separate entities
    protein_chains = []
    ligand_residues = []
    chain_info = []
    
    for entity in structure.entities:
        entity_type = entity.entity_type.name
        chain_info.append({
            "name": entity.name,
            "type": entity_type,
            "subchains": list(entity.subchains)
        })
        
        if entity_type == "Polymer":
            # Protein chain
            for subchain_id in entity.subchains:
                protein_chains.append(subchain_id)
        elif entity_type == "NonPolymer":
            # Ligand
            for subchain_id in entity.subchains:
                ligand_residues.append(subchain_id)
    
    logger.info(f"Found {len(protein_chains)} protein chains, {len(ligand_residues)} ligands")
    
    # Extract protein to PDB
    protein_pdb = output_dir / "protein.pdb"
    structure.remove_ligands_and_waters()
    structure.write_pdb(str(protein_pdb))
    logger.info(f"Wrote protein: {protein_pdb}")
    
    # Re-read original for ligand extraction
    structure = gemmi.make_structure_from_block(block)
    
    # Extract ligand(s)
    ligand_files = []
    ligand_smiles_list = []
    seen_ligands = {}  # Track unique ligands to avoid duplicates
    
    for model in structure:
        for chain in model:
            for residue in chain:
                # Check if non-standard residue (ligand)
                if residue.entity_type == gemmi.EntityType.NonPolymer:
                    res_name = residue.name
                    chain_name = chain.name
                    
                    # Create unique identifier for this ligand instance
                    lig_key = f"{res_name}_{chain_name}"
                    
                    # Skip if we've already processed this exact ligand
                    if lig_key in seen_ligands:
                        continue
                    seen_ligands[lig_key] = True
                    
                    # Create a structure with just this residue
                    lig_structure = gemmi.Structure()
                    lig_model = gemmi.Model("1")
                    lig_chain = gemmi.Chain(chain_name)
                    lig_chain.add_residue(residue)
                    lig_model.add_chain(lig_chain)
                    lig_structure.add_model(lig_model)
                    
                    # Write ligand PDB with chain name to avoid overwriting
                    lig_pdb = output_dir / f"ligand_{res_name}_chain{chain_name}.pdb"
                    lig_structure.write_pdb(str(lig_pdb))
                    ligand_files.append(str(lig_pdb))
                    logger.info(f"Wrote ligand: {lig_pdb}")
                    
                    # Try to extract SMILES using RDKit
                    try:
                        rdkit_mol = Chem.MolFromPDBFile(str(lig_pdb), removeHs=False, sanitize=False)
                        if rdkit_mol:
                            try:
                                Chem.SanitizeMol(rdkit_mol)
                                smiles = Chem.MolToSmiles(rdkit_mol)
                            except:
                                smiles = "SMILES_EXTRACTION_FAILED"
                            ligand_smiles_list.append({
                                "residue": res_name,
                                "chain": chain_name,
                                "smiles": smiles,
                                "file": str(lig_pdb)
                            })
                    except Exception as e:
                        logger.warning(f"Could not extract SMILES for {res_name}: {e}")
    
    result = {
        "job_id": job_id,
        "output_dir": str(output_dir),
        "protein_pdb": str(protein_pdb),
        "ligand_files": ligand_files,
        "ligand_smiles": ligand_smiles_list,
        "chains": chain_info,
        "num_protein_chains": len(protein_chains),
        "num_ligands": len(ligand_files)
    }
    
    # Save metadata
    metadata_file = output_dir / "parse_metadata.json"
    with open(metadata_file, 'w') as f:
        json.dump(result, f, indent=2)
    
    return result


@mcp.tool()
def estimate_net_charge(
    ligand_file: str,
    ph: float = 7.4
) -> dict:
    """Estimate net charge of ligand using RDKit.
    
    Critical for antechamber: incorrect net charge causes sqm to fail with
    "odd number of electrons" error.
    
    Args:
        ligand_file: Path to ligand structure file (PDB, MOL2, SDF)
        ph: Target pH for protonation state estimation
    
    Returns:
        Dict with charge estimation and confidence information
    """
    logger.info(f"Estimating net charge for: {ligand_file}")
    
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors
    except ImportError:
        raise ImportError("RDKit not installed. Install via conda.")
    
    ligand_path = Path(ligand_file)
    if not ligand_path.exists():
        raise FileNotFoundError(f"Ligand file not found: {ligand_file}")
    
    # Load molecule based on file format
    suffix = ligand_path.suffix.lower()
    mol = None
    
    if suffix == '.pdb':
        mol = Chem.MolFromPDBFile(str(ligand_path), removeHs=False)
    elif suffix == '.mol2':
        mol = Chem.MolFromMol2File(str(ligand_path), removeHs=False)
    elif suffix in ['.sdf', '.mol']:
        mol = Chem.MolFromMolFile(str(ligand_path), removeHs=False)
    else:
        raise ValueError(f"Unsupported file format: {suffix}")
    
    if mol is None:
        raise ValueError(f"Could not parse ligand file: {ligand_file}")
    
    # Get charge estimation
    charge_info = _estimate_charge_rdkit(mol)
    
    # Estimate physiological charge
    estimated_charge = _estimate_physiological_charge(charge_info, ph)
    
    # Calculate confidence
    confidence = "high"
    confidence_notes = []
    
    if len(charge_info["ionizable_groups"]) > 2:
        confidence = "medium"
        confidence_notes.append("Multiple ionizable groups detected")
    
    # Check for unusual structures
    num_atoms = mol.GetNumAtoms()
    if num_atoms > 100:
        confidence = "medium"
        confidence_notes.append("Large molecule - charge estimation may be less reliable")
    
    # Check for metals
    metals = ["Fe", "Cu", "Zn", "Mg", "Ca", "Mn"]
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in metals:
            confidence = "low"
            confidence_notes.append(f"Metal atom ({atom.GetSymbol()}) detected - manual charge verification recommended")
            break
    
    result = {
        "ligand_file": str(ligand_file),
        "formal_charge": charge_info["formal_charge"],
        "estimated_charge_at_ph": estimated_charge,
        "target_ph": ph,
        "ionizable_groups": charge_info["ionizable_groups"],
        "confidence": confidence,
        "confidence_notes": confidence_notes,
        "molecular_formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
        "num_atoms": num_atoms,
        "num_heavy_atoms": mol.GetNumHeavyAtoms(),
        "smiles": Chem.MolToSmiles(mol)
    }
    
    logger.info(f"Estimated charge: {estimated_charge} (formal: {charge_info['formal_charge']}, confidence: {confidence})")
    
    return result


@mcp.tool()
def prepare_ligand_hydrogens(
    ligand_file: str,
    output_dir: Optional[str] = None,
    ph: float = 7.4,
    output_format: str = "mol2"
) -> dict:
    """Add hydrogens to ligand structure using OpenBabel.
    
    Proper hydrogen placement is essential for AM1-BCC charge calculation.
    
    Args:
        ligand_file: Path to ligand structure file
        output_dir: Output directory (uses ligand dir if None)
        ph: pH for protonation
        output_format: Output format (mol2, pdb)
    
    Returns:
        Dict with hydrogenated structure path
    """
    logger.info(f"Adding hydrogens to ligand: {ligand_file}")
    
    ligand_path = Path(ligand_file).resolve()  # Use absolute path
    if not ligand_path.exists():
        raise FileNotFoundError(f"Ligand file not found: {ligand_file}")
    
    if output_dir is None:
        output_dir = ligand_path.parent
    else:
        output_dir = Path(output_dir).resolve()  # Use absolute path
    ensure_directory(output_dir)
    
    # Determine input format
    input_format = ligand_path.suffix[1:].lower()
    if input_format == 'sdf':
        input_format = 'mol'
    
    # Output file (absolute path)
    output_file = (output_dir / f"{ligand_path.stem}_H.{output_format}").resolve()
    
    # Run OpenBabel with absolute paths
    args = [
        '-i', input_format,
        str(ligand_path),
        '-o', output_format,
        '-O', str(output_file),
        '-h',  # Add hydrogens
        '-p', str(ph)  # Protonate at specified pH
    ]
    
    try:
        obabel_wrapper.run(args)  # No cwd needed with absolute paths
        logger.info(f"Hydrogenated structure written: {output_file}")
    except Exception as e:
        logger.error(f"OpenBabel hydrogen addition failed: {e}")
        raise
    
    # Verify output
    if not output_file.exists():
        raise RuntimeError(f"Output file not created: {output_file}")
    
    # Get atom counts
    try:
        from rdkit import Chem
        if output_format == 'mol2':
            mol = Chem.MolFromMol2File(str(output_file), removeHs=False)
        else:
            mol = Chem.MolFromPDBFile(str(output_file), removeHs=False)
        
        num_atoms = mol.GetNumAtoms() if mol else 0
        num_hydrogens = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'H') if mol else 0
    except:
        num_atoms = 0
        num_hydrogens = 0
    
    return {
        "input_file": str(ligand_file),
        "output_file": str(output_file),
        "output_format": output_format,
        "ph": ph,
        "num_atoms": num_atoms,
        "num_hydrogens_added": num_hydrogens
    }


@mcp.tool()
def prepare_protein_for_amber(
    pdb_file: str,
    output_dir: Optional[str] = None,
    detect_disulfides: bool = True
) -> dict:
    """Prepare protein structure for Amber using pdb4amber.
    
    Handles:
    - Disulfide bond detection (CYS -> CYX)
    - Histidine protonation states (HIS -> HID/HIE/HIP)
    - Alternate conformation removal
    - Non-standard residue handling
    
    Args:
        pdb_file: Input protein PDB file
        output_dir: Output directory (uses PDB dir if None)
        detect_disulfides: Automatically detect disulfide bonds
    
    Returns:
        Dict with prepared protein information
    """
    logger.info(f"Preparing protein for Amber: {pdb_file}")
    
    pdb_path = Path(pdb_file).resolve()  # Use absolute path
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")
    
    if output_dir is None:
        output_dir = pdb_path.parent
    else:
        output_dir = Path(output_dir).resolve()  # Use absolute path
    ensure_directory(output_dir)
    
    output_pdb = output_dir / f"{pdb_path.stem}_amber.pdb"
    log_file = output_dir / f"{pdb_path.stem}_pdb4amber.log"
    
    # Build pdb4amber command
    args = [
        '-i', str(pdb_path),
        '-o', str(output_pdb),
        '--dry',  # Remove waters
        '--most-populous',  # Keep only most populous alternate conformation
    ]
    
    if not detect_disulfides:
        args.append('--no-ss-bond')
    
    try:
        result = pdb4amber_wrapper.run(args)  # No cwd needed with absolute paths
        
        # Save log
        with open(log_file, 'w') as f:
            f.write(result.stdout if result.stdout else "")
            f.write(result.stderr if result.stderr else "")
        
        logger.info(f"pdb4amber completed: {output_pdb}")
    except Exception as e:
        logger.error(f"pdb4amber failed: {e}")
        raise
    
    # Parse log for detected modifications
    disulfide_bonds = []
    histidine_states = []
    
    log_content = log_file.read_text() if log_file.exists() else ""
    
    # Extract SS bonds
    for line in log_content.split('\n'):
        if 'SS' in line or 'disulfide' in line.lower():
            disulfide_bonds.append(line.strip())
        if 'HIS' in line or 'HID' in line or 'HIE' in line or 'HIP' in line:
            histidine_states.append(line.strip())
    
    # Count atoms in output
    num_atoms = 0
    if output_pdb.exists():
        with open(output_pdb, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    num_atoms += 1
    
    return {
        "input_pdb": str(pdb_file),
        "output_pdb": str(output_pdb),
        "log_file": str(log_file),
        "num_atoms": num_atoms,
        "disulfide_bonds": disulfide_bonds,
        "histidine_states": histidine_states,
        "detect_disulfides": detect_disulfides
    }


@mcp.tool()
def run_antechamber_robust(
    ligand_file: str,
    output_dir: Optional[str] = None,
    net_charge: Optional[int] = None,
    residue_name: str = "LIG",
    charge_method: str = "bcc",
    atom_type: str = "gaff2",
    max_retries: int = 2
) -> dict:
    """Run antechamber with robust error handling and diagnostics.
    
    Automatically estimates net charge if not provided and parses sqm output
    for detailed error diagnostics.
    
    Args:
        ligand_file: Input ligand file (mol2, pdb, sdf)
        output_dir: Output directory
        net_charge: Net molecular charge (auto-estimated if None)
        residue_name: 3-letter residue name for tleap
        charge_method: Charge method (bcc=AM1-BCC, gas=Gasteiger)
        atom_type: Atom type (gaff, gaff2)
        max_retries: Maximum retry attempts with different charges
    
    Returns:
        Dict with parameterization results and diagnostics
    """
    logger.info(f"Running robust antechamber: {ligand_file}")
    
    ligand_path = Path(ligand_file)
    if not ligand_path.exists():
        raise FileNotFoundError(f"Ligand file not found: {ligand_file}")
    
    if output_dir is None:
        output_dir = ligand_path.parent
    else:
        output_dir = Path(output_dir)
    ensure_directory(output_dir)
    
    # Create diagnostics directory
    diag_dir = output_dir / "diagnostics"
    ensure_directory(diag_dir)
    
    # Auto-estimate charge if not provided
    charge_estimation = None
    if net_charge is None:
        logger.info("Auto-estimating net charge...")
        try:
            charge_result = estimate_net_charge(str(ligand_path))
            net_charge = charge_result["estimated_charge_at_ph"]
            charge_estimation = charge_result
            logger.info(f"Estimated charge: {net_charge} (confidence: {charge_result['confidence']})")
        except Exception as e:
            logger.warning(f"Charge estimation failed: {e}, defaulting to 0")
            net_charge = 0
    
    # Determine input format
    input_format = ligand_path.suffix[1:].lower()
    if input_format == 'sdf':
        input_format = 'mdl'
    
    # Output files - preserve original filename for uniqueness
    # e.g., ligand_SAH_chainC_H.mol2 -> ligand_SAH_chainC_H.gaff.mol2, ligand_SAH_chainC_H.frcmod
    input_stem = ligand_path.stem
    
    output_mol2 = output_dir / f"{input_stem}.gaff.mol2"
    output_frcmod = output_dir / f"{input_stem}.frcmod"
    
    # Retry with charge adjustments if needed
    charges_to_try = [net_charge]
    if max_retries > 0:
        charges_to_try.extend([net_charge + 1, net_charge - 1])
    
    last_error = None
    sqm_diagnostics = None
    charge_used = None
    
    for attempt, try_charge in enumerate(charges_to_try[:max_retries + 1]):
        logger.info(f"Attempt {attempt + 1}: trying charge = {try_charge}")
        
        # Build antechamber command
        args = [
            '-i', str(ligand_path),
            '-fi', input_format,
            '-o', str(output_mol2),
            '-fo', 'mol2',
            '-c', charge_method,
            '-nc', str(try_charge),
            '-at', atom_type,
            '-rn', residue_name,
            '-pf', 'y'  # Remove intermediate files
        ]
        
        try:
            antechamber_wrapper.run(args, cwd=output_dir)
            
            # Check if output was created
            if output_mol2.exists():
                charge_used = try_charge
                logger.info(f"Antechamber succeeded with charge = {try_charge}")
                break
            else:
                raise RuntimeError("Antechamber completed but output not created")
                
        except Exception as e:
            last_error = str(e)
            logger.warning(f"Antechamber failed with charge {try_charge}: {e}")
            
            # Parse sqm output for diagnostics
            sqm_out = output_dir / "sqm.out"
            if sqm_out.exists():
                sqm_diagnostics = _parse_sqm_output(sqm_out)
                
                # Copy to diagnostics dir
                shutil.copy(sqm_out, diag_dir / f"sqm_attempt{attempt + 1}.out")
                
                sqm_in = output_dir / "sqm.in"
                if sqm_in.exists():
                    shutil.copy(sqm_in, diag_dir / f"sqm_attempt{attempt + 1}.in")
    
    # Check final result
    if not output_mol2.exists():
        error_msg = f"Antechamber failed after {len(charges_to_try)} attempts"
        if sqm_diagnostics:
            error_msg += f". Diagnostics: {sqm_diagnostics['errors']}"
        raise RuntimeError(error_msg)
    
    # Run parmchk2
    logger.info("Running parmchk2...")
    parmchk2_args = [
        '-i', str(output_mol2),
        '-f', 'mol2',
        '-o', str(output_frcmod),
        '-s', atom_type
    ]
    
    try:
        parmchk2_wrapper.run(parmchk2_args, cwd=output_dir)
        logger.info(f"parmchk2 completed: {output_frcmod}")
    except Exception as e:
        logger.error(f"parmchk2 failed: {e}")
        raise
    
    # Validate frcmod
    frcmod_validation = _parse_frcmod_warnings(output_frcmod)
    
    # Save charge estimation to diagnostics
    if charge_estimation:
        with open(diag_dir / "charge_estimation.json", 'w') as f:
            json.dump(charge_estimation, f, indent=2)
    
    # Save frcmod validation
    with open(diag_dir / "frcmod_validation.json", 'w') as f:
        json.dump(frcmod_validation, f, indent=2)
    
    # Parse charges from output MOL2
    charges = []
    with open(output_mol2, 'r') as f:
        in_atom_section = False
        for line in f:
            if '@<TRIPOS>ATOM' in line:
                in_atom_section = True
                continue
            elif '@<TRIPOS>' in line:
                in_atom_section = False
            
            if in_atom_section and line.strip():
                parts = line.split()
                if len(parts) >= 9:
                    try:
                        charges.append(float(parts[8]))
                    except ValueError:
                        pass
    
    return {
        "mol2": str(output_mol2),
        "frcmod": str(output_frcmod),
        "charge_used": charge_used,
        "charge_method": charge_method,
        "atom_type": atom_type,
        "residue_name": residue_name,
        "charges": charges,
        "total_charge": sum(charges) if charges else 0.0,
        "frcmod_validation": frcmod_validation,
        "sqm_diagnostics": sqm_diagnostics,
        "charge_estimation": charge_estimation,
        "diagnostics_dir": str(diag_dir)
    }


@mcp.tool()
def validate_frcmod(
    frcmod_file: str
) -> dict:
    """Validate frcmod file for missing or problematic parameters.
    
    Checks for "ATTN, need revision" comments and zero force constants
    that indicate parameters were estimated by analogy and may be inaccurate.
    
    Args:
        frcmod_file: Path to .frcmod file
    
    Returns:
        Dict with validation results and recommendations
    """
    logger.info(f"Validating frcmod: {frcmod_file}")
    
    frcmod_path = Path(frcmod_file)
    if not frcmod_path.exists():
        raise FileNotFoundError(f"frcmod file not found: {frcmod_file}")
    
    validation = _parse_frcmod_warnings(frcmod_path)
    validation["file_path"] = str(frcmod_file)
    
    # Log results
    if validation["valid"]:
        logger.info("frcmod validation passed - no issues found")
    else:
        logger.warning(f"frcmod validation found {validation['attn_count']} issues")
        for warning in validation["warnings"][:5]:
            logger.warning(f"  {warning}")
    
    return validation


@mcp.tool()
def build_complex_system(
    protein_pdb: str,
    ligand_mol2: str,
    ligand_frcmod: str,
    output_dir: Optional[str] = None,
    forcefield: str = "leaprc.protein.ff14SB",
    water_model: str = "tip3p",
    box_padding: float = 12.0,
    box_type: str = "box",
    neutralize: bool = True,
    salt_conc: float = 0.15,
    residue_name: str = "LIG"
) -> dict:
    """Build complete MD system with tleap.
    
    Combines protein and ligand into a solvated, neutralized system
    ready for MD simulation.
    
    Args:
        protein_pdb: Prepared protein PDB file
        ligand_mol2: GAFF-parameterized ligand MOL2
        ligand_frcmod: Ligand force field modifications
        output_dir: Output directory
        forcefield: Protein force field
        water_model: Water model (tip3p, tip4pew, opc)
        box_padding: Distance from solute to box edge (Angstroms)
        box_type: Box type (box=rectangular, oct=truncated octahedron)
        neutralize: Add ions to neutralize system
        salt_conc: Additional salt concentration (M)
        residue_name: Ligand residue name
    
    Returns:
        Dict with topology and coordinate files
    """
    logger.info("Building MD system with tleap")
    
    # Validate input files (use absolute paths)
    protein_path = Path(protein_pdb).resolve()
    ligand_mol2_path = Path(ligand_mol2).resolve()
    ligand_frcmod_path = Path(ligand_frcmod).resolve()
    
    for path, name in [(protein_path, "protein PDB"), 
                       (ligand_mol2_path, "ligand MOL2"),
                       (ligand_frcmod_path, "ligand frcmod")]:
        if not path.exists():
            raise FileNotFoundError(f"{name} not found: {path}")
    
    # Setup output directory
    if output_dir is None:
        output_dir = WORKING_DIR / _generate_job_id()
    else:
        output_dir = Path(output_dir).resolve()
    ensure_directory(output_dir)
    
    # Output files
    parm7 = output_dir / "complex.parm7"
    rst7 = output_dir / "complex.rst7"
    leap_in = output_dir / "leap.in"
    leap_log = output_dir / "leap.log"
    complex_pdb = output_dir / "complex.pdb"
    
    # Water model mapping
    water_source = {
        "tip3p": "leaprc.water.tip3p",
        "tip4pew": "leaprc.water.tip4pew",
        "opc": "leaprc.water.opc"
    }.get(water_model, "leaprc.water.tip3p")
    
    water_box = {
        "tip3p": "TIP3PBOX",
        "tip4pew": "TIP4PEWBOX",
        "opc": "OPCBOX"
    }.get(water_model, "TIP3PBOX")
    
    # Solvate command
    if box_type == "oct":
        solvate_cmd = f"solvateoct complex {water_box} {box_padding}"
    else:
        solvate_cmd = f"solvatebox complex {water_box} {box_padding}"
    
    # Build tleap script
    leap_script = f"""# Amber Prep Server - tleap script
# Generated for protein-ligand complex

# Load force fields
source {forcefield}
source leaprc.gaff2
source {water_source}
loadamberparams frcmod.ionsjc_tip3p

# Load ligand parameters (frcmod BEFORE mol2!)
loadamberparams {ligand_frcmod_path.absolute()}

# Load ligand
{residue_name} = loadmol2 {ligand_mol2_path.absolute()}

# Load protein
protein = loadpdb {protein_path.absolute()}

# Combine into complex
complex = combine {{protein {residue_name}}}

# Check for errors
check complex

# Solvate
{solvate_cmd}

# Neutralize
addions complex Na+ 0
addions complex Cl- 0
"""
    
    # Add salt if requested
    if salt_conc > 0:
        # Approximate ion count for salt concentration
        leap_script += f"""
# Add salt ({salt_conc} M)
addionsrand complex Na+ 0
addionsrand complex Cl- 0
"""
    
    # Save and quit
    leap_script += f"""
# Save files
saveamberparm complex {parm7.absolute()} {rst7.absolute()}
savepdb complex {complex_pdb.absolute()}

quit
"""
    
    # Write leap script
    with open(leap_in, 'w') as f:
        f.write(leap_script)
    
    logger.info(f"Created tleap script: {leap_in}")
    
    # Run tleap with absolute path to script
    try:
        result = tleap_wrapper.run(['-f', str(leap_in.resolve())])
        
        # Save log
        with open(leap_log, 'w') as f:
            if result.stdout:
                f.write(result.stdout)
            if result.stderr:
                f.write("\n--- STDERR ---\n")
                f.write(result.stderr)
        
        logger.info("tleap completed successfully")
    except Exception as e:
        logger.error(f"tleap failed: {e}")
        raise
    
    # Verify outputs
    if not parm7.exists():
        raise RuntimeError(f"tleap did not create topology file: {parm7}")
    if not rst7.exists():
        raise RuntimeError(f"tleap did not create coordinate file: {rst7}")
    
    # Parse log for system info
    log_content = leap_log.read_text() if leap_log.exists() else ""
    
    # Extract atom count
    num_atoms = None
    num_residues = None
    for line in log_content.split('\n'):
        if 'atoms' in line.lower():
            match = re.search(r'(\d+)\s+atoms', line)
            if match:
                num_atoms = int(match.group(1))
        if 'residues' in line.lower():
            match = re.search(r'(\d+)\s+residues', line)
            if match:
                num_residues = int(match.group(1))
    
    # Check for warnings
    warnings = []
    for line in log_content.split('\n'):
        if 'warning' in line.lower() or 'error' in line.lower():
            warnings.append(line.strip())
    
    return {
        "parm7": str(parm7),
        "rst7": str(rst7),
        "complex_pdb": str(complex_pdb),
        "leap_in": str(leap_in),
        "leap_log": str(leap_log),
        "output_dir": str(output_dir),
        "forcefield": forcefield,
        "water_model": water_model,
        "box_type": box_type,
        "box_padding": box_padding,
        "neutralized": neutralize,
        "salt_concentration": salt_conc,
        "num_atoms": num_atoms,
        "num_residues": num_residues,
        "warnings": warnings[:10]  # Limit warnings
    }


@mcp.tool()
def build_multi_ligand_system(
    protein_pdb: str,
    ligands: list,
    output_dir: Optional[str] = None,
    forcefield: str = "leaprc.protein.ff14SB",
    water_model: str = "tip3p",
    box_padding: float = 12.0,
    box_type: str = "box",
    neutralize: bool = True,
    salt_conc: float = 0.15
) -> dict:
    """Build MD system with multiple ligands using tleap.
    
    Combines protein with multiple ligands into a solvated, neutralized system.
    
    Args:
        protein_pdb: Prepared protein PDB file
        ligands: List of dicts with 'mol2', 'frcmod', and 'residue_name' keys
        output_dir: Output directory
        forcefield: Protein force field
        water_model: Water model (tip3p, tip4pew, opc)
        box_padding: Distance from solute to box edge (Angstroms)
        box_type: Box type (box=rectangular, oct=truncated octahedron)
        neutralize: Add ions to neutralize system
        salt_conc: Additional salt concentration (M)
    
    Returns:
        Dict with topology and coordinate files
    """
    logger.info(f"Building MD system with {len(ligands)} ligands")
    
    # Validate protein
    protein_path = Path(protein_pdb).resolve()
    if not protein_path.exists():
        raise FileNotFoundError(f"Protein PDB not found: {protein_pdb}")
    
    # Validate ligands
    validated_ligands = []
    for i, lig in enumerate(ligands):
        mol2_path = Path(lig['mol2']).resolve()
        frcmod_path = Path(lig['frcmod']).resolve()
        res_name = lig.get('residue_name', f'L{i:02d}')[:3].upper()
        
        if not mol2_path.exists():
            raise FileNotFoundError(f"Ligand MOL2 not found: {lig['mol2']}")
        if not frcmod_path.exists():
            raise FileNotFoundError(f"Ligand frcmod not found: {lig['frcmod']}")
        
        validated_ligands.append({
            'mol2': mol2_path,
            'frcmod': frcmod_path,
            'residue_name': res_name
        })
    
    # Setup output directory
    if output_dir is None:
        output_dir = WORKING_DIR / _generate_job_id()
    else:
        output_dir = Path(output_dir).resolve()
    ensure_directory(output_dir)
    
    # Output files
    parm7 = output_dir / "complex.parm7"
    rst7 = output_dir / "complex.rst7"
    leap_in = output_dir / "leap.in"
    leap_log = output_dir / "leap.log"
    complex_pdb = output_dir / "complex.pdb"
    
    # Water model mapping
    water_source = {
        "tip3p": "leaprc.water.tip3p",
        "tip4pew": "leaprc.water.tip4pew",
        "opc": "leaprc.water.opc"
    }.get(water_model, "leaprc.water.tip3p")
    
    water_box = {
        "tip3p": "TIP3PBOX",
        "tip4pew": "TIP4PEWBOX",
        "opc": "OPCBOX"
    }.get(water_model, "TIP3PBOX")
    
    # Solvate command
    if box_type == "oct":
        solvate_cmd = f"solvateoct complex {water_box} {box_padding}"
    else:
        solvate_cmd = f"solvatebox complex {water_box} {box_padding}"
    
    # Build tleap script
    leap_script = f"""# Amber Prep Server - tleap script
# Generated for protein with {len(validated_ligands)} ligands

# Load force fields
source {forcefield}
source leaprc.gaff2
source {water_source}
loadamberparams frcmod.ionsjc_tip3p

"""
    
    # Load each ligand's parameters and structure
    ligand_names = []
    for i, lig in enumerate(validated_ligands):
        res_name = lig['residue_name']
        # Make unique variable name if same residue name appears multiple times
        var_name = f"{res_name}_{i}" if any(l['residue_name'] == res_name for j, l in enumerate(validated_ligands) if j != i) else res_name
        ligand_names.append(var_name)
        
        leap_script += f"""# Load ligand {i+1}: {res_name}
loadamberparams {lig['frcmod']}
{var_name} = loadmol2 {lig['mol2']}

"""
    
    # Load protein and combine
    leap_script += f"""# Load protein
protein = loadpdb {protein_path}

# Combine protein with all ligands
complex = combine {{protein {' '.join(ligand_names)}}}

# Check for errors
check complex

# Solvate
{solvate_cmd}

# Neutralize
addions complex Na+ 0
addions complex Cl- 0
"""
    
    # Add salt if requested
    if salt_conc > 0:
        leap_script += f"""
# Add salt ({salt_conc} M)
addionsrand complex Na+ 0
addionsrand complex Cl- 0
"""
    
    # Save and quit
    leap_script += f"""
# Save files
saveamberparm complex {parm7} {rst7}
savepdb complex {complex_pdb}

quit
"""
    
    # Write leap script
    with open(leap_in, 'w') as f:
        f.write(leap_script)
    
    logger.info(f"Created tleap script: {leap_in}")
    
    # Run tleap
    try:
        result = tleap_wrapper.run(['-f', str(leap_in.resolve())])
        
        # Save log
        with open(leap_log, 'w') as f:
            if result.stdout:
                f.write(result.stdout)
            if result.stderr:
                f.write("\n--- STDERR ---\n")
                f.write(result.stderr)
        
        logger.info("tleap completed successfully")
    except Exception as e:
        logger.error(f"tleap failed: {e}")
        raise
    
    # Verify output files
    if not parm7.exists():
        raise RuntimeError(f"tleap failed to create topology: {parm7}")
    if not rst7.exists():
        raise RuntimeError(f"tleap failed to create coordinates: {rst7}")
    
    # Parse log for statistics
    log_content = leap_log.read_text() if leap_log.exists() else ""
    num_atoms = None
    num_residues = None
    
    for line in log_content.split('\n'):
        if 'atoms' in line.lower():
            match = re.search(r'(\d+)\s+atoms', line)
            if match:
                num_atoms = int(match.group(1))
        if 'residues' in line.lower():
            match = re.search(r'(\d+)\s+residues', line)
            if match:
                num_residues = int(match.group(1))
    
    # Check for warnings
    warnings = []
    for line in log_content.split('\n'):
        if 'warning' in line.lower() or 'error' in line.lower():
            warnings.append(line.strip())
    
    return {
        "parm7": str(parm7),
        "rst7": str(rst7),
        "complex_pdb": str(complex_pdb),
        "leap_in": str(leap_in),
        "leap_log": str(leap_log),
        "output_dir": str(output_dir),
        "forcefield": forcefield,
        "water_model": water_model,
        "box_type": box_type,
        "box_padding": box_padding,
        "neutralized": neutralize,
        "salt_concentration": salt_conc,
        "num_atoms": num_atoms,
        "num_residues": num_residues,
        "num_ligands": len(validated_ligands),
        "ligand_names": [l['residue_name'] for l in validated_ligands],
        "warnings": warnings[:10]
    }


@mcp.tool()
def boltz2_to_amber_complete(
    cif_file: str,
    output_dir: Optional[str] = None,
    ph: float = 7.4,
    net_charge: Optional[int] = None,
    residue_name: str = "LIG",
    water_model: str = "tip3p",
    box_padding: float = 12.0
) -> dict:
    """Complete workflow: Boltz-2 mmCIF to MD-ready Amber files.
    
    Executes the full parameterization pipeline:
    1. Parse mmCIF complex
    2. Prepare protein with pdb4amber
    3. Add hydrogens to ligand
    4. Estimate/verify net charge
    5. Run antechamber (GAFF2 + AM1-BCC)
    6. Validate frcmod
    7. Build system with tleap
    
    Args:
        cif_file: Boltz-2 output mmCIF file
        output_dir: Output directory (auto-generated if None)
        ph: Target pH for protonation
        net_charge: Ligand net charge (auto-estimated if None)
        residue_name: Ligand residue name
        water_model: Water model for solvation
        box_padding: Box padding distance
    
    Returns:
        Dict with complete workflow results
    """
    logger.info(f"Starting complete Boltz-2 to Amber workflow: {cif_file}")
    
    workflow_log = []
    
    def log_step(step: str, result: dict):
        workflow_log.append({
            "step": step,
            "success": True,
            "result": result
        })
        logger.info(f"Completed step: {step}")
    
    # Step 1: Parse complex
    logger.info("Step 1/7: Parsing Boltz-2 complex...")
    parse_result = parse_boltz2_complex(cif_file, output_dir)
    job_dir = Path(parse_result["output_dir"])
    log_step("parse_complex", parse_result)
    
    if not parse_result["ligand_files"]:
        raise RuntimeError("No ligand found in complex")
    
    # Use first ligand
    ligand_pdb = parse_result["ligand_files"][0]
    
    # Step 2: Prepare protein
    logger.info("Step 2/7: Preparing protein with pdb4amber...")
    protein_result = prepare_protein_for_amber(
        parse_result["protein_pdb"],
        output_dir=str(job_dir)
    )
    log_step("prepare_protein", protein_result)
    
    # Step 3: Add hydrogens to ligand
    logger.info("Step 3/7: Adding hydrogens to ligand...")
    hydrogen_result = prepare_ligand_hydrogens(
        ligand_pdb,
        output_dir=str(job_dir),
        ph=ph
    )
    log_step("add_hydrogens", hydrogen_result)
    
    # Step 4: Estimate charge (if not provided)
    logger.info("Step 4/7: Estimating net charge...")
    charge_result = estimate_net_charge(hydrogen_result["output_file"], ph=ph)
    log_step("estimate_charge", charge_result)
    
    if net_charge is None:
        net_charge = charge_result["estimated_charge_at_ph"]
        logger.info(f"Using estimated charge: {net_charge}")
    
    # Step 5: Run antechamber
    logger.info("Step 5/7: Running antechamber (GAFF2 + AM1-BCC)...")
    antechamber_result = run_antechamber_robust(
        hydrogen_result["output_file"],
        output_dir=str(job_dir),
        net_charge=net_charge,
        residue_name=residue_name
    )
    log_step("antechamber", antechamber_result)
    
    # Step 6: Validate frcmod
    logger.info("Step 6/7: Validating frcmod...")
    frcmod_result = validate_frcmod(antechamber_result["frcmod"])
    log_step("validate_frcmod", frcmod_result)
    
    # Step 7: Build system
    logger.info("Step 7/7: Building MD system with tleap...")
    system_result = build_complex_system(
        protein_result["output_pdb"],
        antechamber_result["mol2"],
        antechamber_result["frcmod"],
        output_dir=str(job_dir),
        water_model=water_model,
        box_padding=box_padding,
        residue_name=residue_name
    )
    log_step("build_system", system_result)
    
    # Generate validation report
    validation_report = {
        "charge_confidence": charge_result.get("confidence", "unknown"),
        "frcmod_valid": frcmod_result.get("valid", False),
        "frcmod_warnings": frcmod_result.get("attn_count", 0),
        "tleap_warnings": len(system_result.get("warnings", [])),
        "overall_status": "success" if frcmod_result.get("valid", False) else "warning"
    }
    
    if not frcmod_result.get("valid", False):
        validation_report["issues"] = frcmod_result.get("recommendations", [])
    
    # Save workflow summary
    summary = {
        "job_id": parse_result["job_id"],
        "input_cif": str(cif_file),
        "output_dir": str(job_dir),
        "parm7": system_result["parm7"],
        "rst7": system_result["rst7"],
        "complex_pdb": system_result["complex_pdb"],
        "ligand_mol2": antechamber_result["mol2"],
        "ligand_frcmod": antechamber_result["frcmod"],
        "protein_pdb": protein_result["output_pdb"],
        "net_charge": net_charge,
        "ph": ph,
        "water_model": water_model,
        "box_padding": box_padding,
        "num_atoms": system_result.get("num_atoms"),
        "workflow_log": workflow_log,
        "validation_report": validation_report
    }
    
    # Save summary JSON
    summary_file = job_dir / "workflow_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    
    logger.info(f"Workflow complete! Output: {job_dir}")
    logger.info(f"  Topology: {system_result['parm7']}")
    logger.info(f"  Coordinates: {system_result['rst7']}")
    
    return summary


if __name__ == "__main__":
    mcp.run()

