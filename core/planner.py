"""
Planner - Natural language to MD workflow conversion.

Uses LM Studio to convert user queries into structured MD workflows.
"""

import logging
import json
from pathlib import Path
from typing import Optional, Dict, List
from datetime import datetime

from .llm_client import LMStudioClient
from .models import WorkflowPlan, WorkflowStep, SystemType
from .utils import generate_unique_id

logger = logging.getLogger(__name__)


class MDWorkflowPlanner:
    """MD Workflow Planner with LLM integration"""
    
    WORKFLOW_TEMPLATES = {
        "protein_only": [
            ("fetch_pdb", "structure_server"),
            ("clean_structure", "structure_server"),
            ("build_system_tleap", "assembly_server"),
            ("create_openmm_workflow", "protocol_server"),
            ("package_system", "export_server"),
        ],
        "protein_ligand_complex": [
            ("fetch_pdb", "structure_server"),
            ("clean_structure", "structure_server"),
            ("smiles_to_3d", "ligand_server"),
            ("generate_gaff_params", "ligand_server"),
            ("create_ligand_lib", "ligand_server"),
            ("dock_ligand_smina", "docking_server"),
            ("build_system_tleap", "assembly_server"),
            ("create_openmm_workflow", "protocol_server"),
            ("package_system", "export_server"),
        ],
        "boltz2_denovo": [
            ("predict_structure_boltz2", "structure_server"),
            ("validate_structure", "structure_server"),
            ("build_system_tleap", "assembly_server"),
            ("create_openmm_workflow", "protocol_server"),
            ("package_system", "export_server"),
        ],
        "boltz2_complex": [
            ("predict_complex_with_affinity", "structure_server"),
            ("parameterize_ligand_complete", "ligand_server"),
            ("build_system_tleap", "assembly_server"),
            ("create_openmm_workflow", "protocol_server"),
            ("package_system", "export_server"),
        ],
    }
    
    def __init__(self, llm_client: Optional[LMStudioClient] = None):
        """Initialize planner
        
        Args:
            llm_client: LM Studio client (creates default if None)
        """
        self.llm = llm_client or LMStudioClient()
        logger.info("MD Workflow Planner initialized")
    
    def plan_from_query(self, user_query: str) -> WorkflowPlan:
        """Create workflow plan from natural language query
        
        Args:
            user_query: User's natural language request
        
        Returns:
            WorkflowPlan object
        """
        logger.info(f"Planning workflow from query: {user_query}")
        
        # Analyze query with LLM
        analyzed = self._analyze_query(user_query)
        
        # Select workflow template
        system_type = analyzed.get("system_type", "protein_only")
        template = self.WORKFLOW_TEMPLATES.get(system_type, self.WORKFLOW_TEMPLATES["protein_only"])
        
        # Create workflow steps
        steps = []
        for idx, (tool_name, server_name) in enumerate(template):
            step = WorkflowStep(
                step_id=f"step_{idx:03d}",
                name=tool_name,
                mcp_tool=tool_name,
                mcp_server=server_name,
                params={},
                dependencies=[f"step_{idx-1:03d}"] if idx > 0 else []
            )
            steps.append(step)
        
        # Create plan
        plan = WorkflowPlan(
            plan_id=generate_unique_id("plan"),
            system_type=SystemType(system_type),
            steps=steps,
            description=user_query,
            created_at=datetime.now().isoformat()
        )
        
        return plan
    
    def _analyze_query(self, query: str) -> Dict:
        """Analyze query with LLM
        
        Args:
            query: User query
        
        Returns:
            Dict with extracted information
        """
        system_prompt = """You are an MD workflow expert. Analyze the user's request and extract:
- system_type: protein_only, protein_ligand_complex, boltz2_denovo, boltz2_complex
- components: list of components (protein, ligand, water, ions)
- input_source: where structures come from

Respond with JSON only."""
        
        user_prompt = f"Analyze this MD workflow request:\n{query}"
        
        try:
            if not self.llm.is_available():
                logger.warning("LM Studio not available, using default analysis")
                return {"system_type": "protein_only", "components": ["protein"]}
            
            response = self.llm.complete_sync(user_prompt, system=system_prompt, temperature=0.1)
            
            # Parse JSON response
            import re
            json_match = re.search(r'\{.*\}', response, re.DOTALL)
            if json_match:
                return json.loads(json_match.group())
            
        except Exception as e:
            logger.warning(f"LLM analysis failed: {e}")
        
        # Fallback
        return {"system_type": "protein_only", "components": ["protein"]}
    
    def save_plan(self, plan: WorkflowPlan, output_path: str):
        """Save plan to markdown file
        
        Args:
            plan: WorkflowPlan object
            output_path: Output file path
        """
        logger.info(f"Saving plan to: {output_path}")
        
        md_content = self._generate_markdown(plan)
        
        with open(output_path, 'w') as f:
            f.write(md_content)
        
        logger.info("Plan saved successfully")
    
    def _generate_markdown(self, plan: WorkflowPlan) -> str:
        """Generate markdown representation of plan
        
        Args:
            plan: WorkflowPlan object
        
        Returns:
            Markdown string
        """
        md = f"""# MD Workflow Plan: {plan.plan_id}

## Description
{plan.description}

## System Type
{plan.system_type.value}

## Workflow Steps

"""
        
        for step in plan.steps:
            md += f"""### {step.step_id}: {step.name}

- **Server**: {step.mcp_server}
- **Tool**: {step.mcp_tool}
- **Dependencies**: {', '.join(step.dependencies) if step.dependencies else 'None'}
- **Status**: {step.status}

"""
        
        return md


class QueryAnalyzer:
    """Analyze natural language queries"""
    pass


class StepSelector:
    """Select workflow steps based on system type"""
    pass


class ParameterInferencer:
    """Infer missing parameters"""
    pass


class DAGBuilder:
    """Build dependency graph"""
    pass

