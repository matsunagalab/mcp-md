"""
Workflow Nodes - Individual step implementations for LangGraph.

Each node represents a step in the MD workflow pipeline.
Nodes receive WorkflowState and return updates to merge into state.
"""

import logging
from typing import Dict
from datetime import datetime
from pathlib import Path

from langchain_core.tools import BaseTool

from .workflow_state import WorkflowState

logger = logging.getLogger(__name__)


def planner_node(state: WorkflowState) -> Dict:
    """Planning node: Initialize workflow with fixed skeleton
    
    Args:
        state: Current workflow state
    
    Returns:
        State updates (current_step, user_preferences)
    """
    logger.info(f"Planning workflow for query: {state.get('query', 'N/A')}")
    
    # Set default user preferences if not provided
    user_prefs = state.get("user_preferences", {})
    if not user_prefs:
        user_prefs = {
            "ph": 7.4,
            "salt_concentration": 0.15,
            "water_model": "TIP3P",
            "force_field": "ff19SB",
            "box_padding": 12.0,
        }
    
    return {
        "current_step": "fetch",
        "user_preferences": user_prefs,
        "decision_log": [{
            "timestamp": datetime.utcnow().isoformat(),
            "step": "planner",
            "action": "initialize_workflow",
            "params": user_prefs
        }]
    }


async def structure_fetch_node(state: WorkflowState, tools: Dict[str, BaseTool]) -> Dict:
    """Structure fetching node: Retrieve PDB structure
    
    Args:
        state: Current workflow state
        tools: Dictionary of available MCP tools
    
    Returns:
        State updates (outputs, current_step, decision_log)
    """
    logger.info(f"Fetching structure: PDB ID = {state.get('pdb_id')}")
    
    pdb_id = state.get("pdb_id")
    if not pdb_id:
        return {
            "error": "No PDB ID provided for structure fetch",
            "retry_count": state.get("retry_count", 0) + 1
        }
    
    try:
        fetch_pdb_tool = tools.get("fetch_pdb")
        if not fetch_pdb_tool:
            raise ValueError("fetch_pdb tool not available")
        
        result = await fetch_pdb_tool.ainvoke({"pdb_id": pdb_id})
        
        return {
            "outputs": {**state.get("outputs", {}), "structure": result},
            "current_step": "repair",
            "decision_log": [{
                "timestamp": datetime.utcnow().isoformat(),
                "step": "fetch",
                "tool": "fetch_pdb",
                "params": {"pdb_id": pdb_id},
                "result": "success"
            }],
            "error": None,
            "retry_count": 0
        }
    
    except Exception as e:
        logger.error(f"Structure fetch failed: {e}")
        return {
            "error": f"Failed to fetch PDB {pdb_id}: {str(e)}",
            "retry_count": state.get("retry_count", 0) + 1
        }


async def repair_node(state: WorkflowState, tools: Dict[str, BaseTool]) -> Dict:
    """Structure repair node: Clean and fix PDB structure
    
    Args:
        state: Current workflow state
        tools: Dictionary of available MCP tools
    
    Returns:
        State updates (outputs, current_step, decision_log)
    """
    logger.info("Repairing structure with PDBFixer")
    
    structure_data = state.get("outputs", {}).get("structure")
    if not structure_data:
        return {
            "error": "No structure data available for repair",
            "retry_count": state.get("retry_count", 0) + 1
        }
    
    try:
        clean_structure_tool = tools.get("clean_structure")
        if not clean_structure_tool:
            raise ValueError("clean_structure tool not available")
        
        result = await clean_structure_tool.ainvoke({
            "pdb_file": structure_data.get("file_path")
        })
        
        return {
            "outputs": {**state.get("outputs", {}), "repaired_structure": result},
            "current_step": "ligand_param",
            "decision_log": [{
                "timestamp": datetime.utcnow().isoformat(),
                "step": "repair",
                "tool": "clean_structure",
                "result": "success"
            }],
            "error": None,
            "retry_count": 0
        }
    
    except Exception as e:
        logger.error(f"Structure repair failed: {e}")
        return {
            "error": f"Failed to repair structure: {str(e)}",
            "retry_count": state.get("retry_count", 0) + 1
        }


async def ligand_param_node(state: WorkflowState, tools: Dict[str, BaseTool]) -> Dict:
    """Ligand parameterization node: Generate GAFF2 parameters
    
    Args:
        state: Current workflow state
        tools: Dictionary of available MCP tools
    
    Returns:
        State updates (outputs, current_step, decision_log)
    """
    logger.info("Parameterizing ligand with GAFF2/AM1-BCC")
    
    ligand_smiles = state.get("ligand_smiles")
    if not ligand_smiles:
        # No ligand, skip this step
        return {
            "current_step": "assemble",
            "decision_log": [{
                "timestamp": datetime.utcnow().isoformat(),
                "step": "ligand_param",
                "action": "skipped",
                "reason": "No ligand provided"
            }]
        }
    
    try:
        param_tool = tools.get("parameterize_ligand_complete")
        if not param_tool:
            raise ValueError("parameterize_ligand_complete tool not available")
        
        result = await param_tool.ainvoke({
            "smiles": ligand_smiles,
            "charge_method": "bcc"  # AM1-BCC
        })
        
        return {
            "outputs": {**state.get("outputs", {}), "ligand_params": result},
            "current_step": "complex",
            "decision_log": [{
                "timestamp": datetime.utcnow().isoformat(),
                "step": "ligand_param",
                "tool": "parameterize_ligand_complete",
                "params": {"smiles": ligand_smiles, "method": "AM1-BCC"},
                "result": "success"
            }],
            "error": None,
            "retry_count": 0
        }
    
    except Exception as e:
        logger.error(f"Ligand parameterization failed: {e}")
        return {
            "error": f"Failed to parameterize ligand: {str(e)}",
            "retry_count": state.get("retry_count", 0) + 1
        }


async def complex_node(state: WorkflowState, tools: Dict[str, BaseTool]) -> Dict:
    """Complex generation node: Generate protein-ligand complex
    
    Args:
        state: Current workflow state
        tools: Dictionary of available MCP tools
    
    Returns:
        State updates (outputs, current_step, decision_log)
    """
    logger.info("Generating protein-ligand complex")
    
    # Check if we have both protein and ligand
    has_ligand = state.get("ligand_smiles") is not None
    
    if not has_ligand:
        # No complex needed, proceed to assembly
        return {
            "current_step": "assemble",
            "decision_log": [{
                "timestamp": datetime.utcnow().isoformat(),
                "step": "complex",
                "action": "skipped",
                "reason": "No ligand for complex generation"
            }]
        }
    
    try:
        # Use Boltz-2 for complex prediction
        boltz2_complex_tool = tools.get("boltz2_complex")
        if not boltz2_complex_tool:
            raise ValueError("boltz2_complex tool not available")
        
        protein_data = state.get("outputs", {}).get("repaired_structure")
        ligand_smiles = state.get("ligand_smiles")
        
        result = await boltz2_complex_tool.ainvoke({
            "protein_pdb": protein_data.get("output"),
            "ligand_smiles": ligand_smiles,
            "top_k": 5
        })
        
        return {
            "outputs": {**state.get("outputs", {}), "complex": result},
            "current_step": "assemble",
            "decision_log": [{
                "timestamp": datetime.utcnow().isoformat(),
                "step": "complex",
                "tool": "boltz2_complex",
                "params": {"top_k": 5},
                "result": "success"
            }],
            "error": None,
            "retry_count": 0
        }
    
    except Exception as e:
        logger.error(f"Complex generation failed: {e}")
        return {
            "error": f"Failed to generate complex: {str(e)}",
            "retry_count": state.get("retry_count", 0) + 1
        }


async def assemble_node(state: WorkflowState, tools: Dict[str, BaseTool]) -> Dict:
    """System assembly node: Build complete MD system with tleap
    
    Args:
        state: Current workflow state
        tools: Dictionary of available MCP tools
    
    Returns:
        State updates (outputs, current_step, decision_log)
    """
    logger.info("Assembling MD system with tleap")
    
    try:
        build_system_tool = tools.get("build_system_tleap")
        if not build_system_tool:
            raise ValueError("build_system_tleap tool not available")
        
        outputs = state.get("outputs", {})
        user_prefs = state.get("user_preferences", {})
        
        # Determine input structure (complex or protein only)
        if "complex" in outputs:
            protein_pdb = outputs["complex"]["structures"][0]
        else:
            protein_pdb = outputs.get("repaired_structure", {}).get("output")
        
        # Determine ligand parameters
        ligand_lib = outputs.get("ligand_params", {}).get("library")
        ligand_frcmod = outputs.get("ligand_params", {}).get("frcmod")
        
        result = await build_system_tool.ainvoke({
            "protein_pdb": protein_pdb,
            "ligand_lib": ligand_lib,
            "ligand_frcmod": ligand_frcmod,
            "forcefield": user_prefs.get("force_field", "ff19SB"),
            "water_model": user_prefs.get("water_model", "TIP3P"),
            "salt_concentration": user_prefs.get("salt_concentration", 0.15),
            "box_padding": user_prefs.get("box_padding", 12.0)
        })
        
        return {
            "outputs": {**outputs, "system": result},
            "current_step": "qc",
            "decision_log": [{
                "timestamp": datetime.utcnow().isoformat(),
                "step": "assemble",
                "tool": "build_system_tleap",
                "params": user_prefs,
                "result": "success"
            }],
            "error": None,
            "retry_count": 0
        }
    
    except Exception as e:
        logger.error(f"System assembly failed: {e}")
        return {
            "error": f"Failed to assemble system: {str(e)}",
            "retry_count": state.get("retry_count", 0) + 1
        }


async def qc_node(state: WorkflowState, tools: Dict[str, BaseTool]) -> Dict:
    """Quality control node: Run QC checks and minimization
    
    Args:
        state: Current workflow state
        tools: Dictionary of available MCP tools
    
    Returns:
        State updates (outputs, current_step, decision_log)
    """
    logger.info("Running QC checks and minimization")
    
    try:
        minimize_tool = tools.get("openmm_minimize")
        if not minimize_tool:
            raise ValueError("openmm_minimize tool not available")
        
        system_data = state.get("outputs", {}).get("system", {})
        
        result = await minimize_tool.ainvoke({
            "prmtop": system_data.get("prmtop"),
            "inpcrd": system_data.get("inpcrd"),
            "max_iterations": 5000
        })
        
        return {
            "outputs": {**state.get("outputs", {}), "minimized": result},
            "current_step": "complete",
            "decision_log": [{
                "timestamp": datetime.utcnow().isoformat(),
                "step": "qc",
                "tool": "openmm_minimize",
                "params": {"max_iterations": 5000},
                "result": "success"
            }],
            "error": None,
            "retry_count": 0
        }
    
    except Exception as e:
        logger.error(f"QC/minimization failed: {e}")
        return {
            "error": f"Failed QC/minimization: {str(e)}",
            "retry_count": state.get("retry_count", 0) + 1
        }


def should_retry(state: WorkflowState) -> str:
    """Conditional edge function: Determine if step should retry
    
    Args:
        state: Current workflow state
    
    Returns:
        Edge name to follow: "retry", "human_feedback", or "continue"
    """
    error = state.get("error")
    retry_count = state.get("retry_count", 0)
    
    if error and retry_count < 1:
        logger.warning(f"Retrying after error: {error}")
        return "retry"
    elif error:
        logger.error(f"Max retries exceeded, requesting human feedback: {error}")
        return "human_feedback"
    
    return "continue"

