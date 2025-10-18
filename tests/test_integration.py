"""
Integration tests for MCP-MD workflow.
"""

import pytest
from pathlib import Path
from core.planner import MDWorkflowPlanner
from core.workflow import WorkflowEngine
from core.models import SystemType


def test_planner_creates_workflow():
    """Test workflow planner"""
    planner = MDWorkflowPlanner()
    
    query = "Create MD system for PDB 1ABC with water and ions"
    plan = planner.plan_from_query(query)
    
    assert plan is not None
    assert len(plan.steps) > 0
    assert plan.system_type == SystemType.PROTEIN_ONLY


def test_workflow_execution():
    """Test workflow execution"""
    planner = MDWorkflowPlanner()
    engine = WorkflowEngine()
    
    query = "Build protein-ligand complex MD system"
    plan = planner.plan_from_query(query)
    
    # This would actually execute the workflow
    # For now, just validate plan structure
    assert plan is not None
    assert all(step.mcp_server for step in plan.steps)


def test_plan_saving():
    """Test plan markdown saving"""
    planner = MDWorkflowPlanner()
    
    query = "Predict structure from FASTA and create MD system"
    plan = planner.plan_from_query(query)
    
    output_path = Path("test_plan.md")
    planner.save_plan(plan, str(output_path))
    
    assert output_path.exists()
    
    # Cleanup
    output_path.unlink()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

