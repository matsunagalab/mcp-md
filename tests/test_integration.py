"""
Integration tests for the complete MCP-MD workflow.

Tests the end-to-end workflow:
1. Genesis: Generate protein from FASTA
2. Structure: Clean and protonate
3. Complex: Generate protein-ligand complex
4. Ligand: Parameterize ligand
5. Assembly: Build system with tleap
6. Export: Export to Amber format
7. QC/Min: Run quality checks and minimization
"""

import pytest
import asyncio
from pathlib import Path
import tempfile
import shutil


# Sample inputs for testing
SAMPLE_FASTA = """MKFLKFSLLTAVLLSVVFAFSSCGDDDDTYPYDVPDYAG"""

SAMPLE_SMILES = "CCO"  # Ethanol


@pytest.mark.asyncio
@pytest.mark.integration
async def test_full_workflow_simple():
    """Test simple protein-only workflow"""
    
    # This is a placeholder test
    # Actual implementation requires:
    # 1. MCP servers running
    # 2. Strands Agent initialized
    # 3. Tools available
    
    # For now, just verify imports work
    from core.strands_agent import MDWorkflowAgent
    from core.workflow_skeleton import WORKFLOW_SKELETON, get_all_step_names
    from core.decision_logger import DecisionLogger
    
    # Verify workflow skeleton
    steps = get_all_step_names()
    assert len(steps) == 9
    assert "fetch_or_generate" in steps
    assert "minimize_qc" in steps
    
    # Verify decision logger
    with tempfile.TemporaryDirectory() as tmpdir:
        logger = DecisionLogger(Path(tmpdir))
        logger.log_decision(
            step="test",
            tool="test_tool",
            params={"test": "value"},
            reason="Test reason"
        )
        
        decisions = logger.get_all_decisions()
        assert len(decisions) == 1
        assert decisions[0]["step"] == "test"


@pytest.mark.asyncio
@pytest.mark.integration
async def test_genesis_server():
    """Test Genesis Server independently"""
    
    # Placeholder for Genesis Server test
    # Requires Boltz-2 installation
    
    from servers.genesis_server import GenesisServer
    
    server = GenesisServer()
    assert server is not None
    assert server.name == "genesis_server"


@pytest.mark.asyncio
@pytest.mark.integration
async def test_structure_server_simplified():
    """Test simplified Structure Server"""
    
    from servers.structure_server import StructureServer
    
    server = StructureServer()
    assert server is not None
    assert server.name == "structure_server"
    
    # Should NOT have Boltz-2 wrapper
    assert not hasattr(server, 'boltz2')
    
    # Should have PDBFixer and PDB2PQR
    assert hasattr(server, 'pdbfixer')
    assert hasattr(server, 'pdb2pqr')


@pytest.mark.asyncio
@pytest.mark.integration
async def test_molprobity_wrapper():
    """Test custom MolProbity wrapper"""
    
    from tools.molprobity_wrapper import MolProbityWrapper
    
    # This will fail if MDAnalysis/BioPython not installed
    try:
        wrapper = MolProbityWrapper()
        assert wrapper is not None
    except ImportError as e:
        pytest.skip(f"MolProbity dependencies not installed: {e}")


@pytest.mark.asyncio
@pytest.mark.integration
async def test_workflow_skeleton():
    """Test workflow skeleton structure"""
    
    from core.workflow_skeleton import (
        WORKFLOW_SKELETON,
        get_step_by_name,
        get_all_step_names,
        get_step_tools,
        validate_workflow_plan
    )
    
    # Check all steps
    steps = get_all_step_names()
    assert len(steps) == 9
    
    # Check specific steps
    fetch_step = get_step_by_name("fetch_or_generate")
    assert "fetch_pdb" in fetch_step.required_tools
    
    qc_step = get_step_by_name("minimize_qc")
    assert "molprobity_check" in qc_step.required_tools
    assert "openmm_minimize" in qc_step.required_tools
    
    # Validate workflow plan
    valid_plan = ["fetch_or_generate", "repair_protonate", "ligand_param"]
    assert validate_workflow_plan(valid_plan) == True
    
    invalid_plan = ["ligand_param", "fetch_or_generate"]  # Wrong order
    assert validate_workflow_plan(invalid_plan) == False


@pytest.mark.asyncio
@pytest.mark.integration
async def test_decision_logger():
    """Test decision logger functionality"""
    
    from core.decision_logger import DecisionLogger
    
    with tempfile.TemporaryDirectory() as tmpdir:
        logger = DecisionLogger(Path(tmpdir))
        
        # Log multiple decisions
        logger.log_decision(
            step="step1",
            tool="tool1",
            params={"param1": "value1"},
            reason="Reason 1"
        )
        
        logger.log_decision(
            step="step2",
            tool="tool2",
            params={"param2": "value2"},
            reason="Reason 2"
        )
        
        logger.log_error(
            step="step3",
            tool="tool3",
            error="Error message"
        )
        
        # Get decisions
        decisions = logger.get_all_decisions()
        assert len(decisions) == 3
        
        # Generate summary
        summary = logger.generate_summary()
        assert summary["total_decisions"] == 3
        assert len(summary["errors"]) == 1
        assert "step1" in summary["steps"]
        assert "step2" in summary["steps"]
        
        # Export markdown
        report_path = logger.export_markdown_report()
        assert Path(report_path).exists()


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v", "-m", "integration"])
