"""
Strands Agent - Core workflow orchestrator with LM Studio integration.

Provides:
- MCP client management for all servers
- Interactive CLI interface
- Workflow orchestration
- Decision logging
"""

import logging
import os
from pathlib import Path
from typing import List, Optional, Dict, Any

try:
    from mcp import StdioServerParameters, stdio_client
    from strands import Agent
    from strands.models.openai import OpenAIModel
    from strands.tools.mcp import MCPClient
    from prompt_toolkit import prompt
except ImportError as e:
    raise ImportError(
        f"Required packages not installed: {e}\n"
        "Install with: pip install strands-ai prompt-toolkit mcp"
    )

from .utils import setup_logger
from .decision_logger import DecisionLogger
from .workflow_skeleton import WORKFLOW_SKELETON

logger = setup_logger(__name__)


class MDWorkflowAgent:
    """MD Workflow orchestration agent using Strands + LM Studio"""
    
    def __init__(
        self,
        lm_studio_url: str = None,
        model_id: str = None,
        run_dir: Path = None
    ):
        """Initialize MD Workflow Agent
        
        Args:
            lm_studio_url: LM Studio API URL (default: http://localhost:1234/v1)
            model_id: Model ID (default: gemma-3-12b)
            run_dir: Run directory for outputs (default: runs/<timestamp>)
        """
        # Get config from env or defaults
        self.lm_studio_url = lm_studio_url or os.getenv(
            "LM_STUDIO_BASE_URL",
            "http://localhost:1234/v1"
        )
        self.model_id = model_id or os.getenv(
            "LM_STUDIO_MODEL",
            "gemma-3-12b"
        )
        
        # Setup run directory
        if run_dir is None:
            from datetime import datetime
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            run_dir = Path(f"runs/{timestamp}")
        self.run_dir = Path(run_dir)
        self.run_dir.mkdir(parents=True, exist_ok=True)
        
        # Decision logger
        self.decision_logger = DecisionLogger(self.run_dir)
        
        # Initialize LM Studio model
        logger.info(f"Initializing LM Studio: {self.lm_studio_url}, model={self.model_id}")
        self.model = OpenAIModel(
            client_args={
                "api_key": "lm-studio",  # Dummy key for local LM Studio
                "base_url": self.lm_studio_url,
            },
            model_id=self.model_id,
        )
        
        # MCP clients
        self.mcp_clients = self._create_mcp_clients()
        
        logger.info("MD Workflow Agent initialized")
    
    def _create_mcp_clients(self) -> Dict[str, MCPClient]:
        """Create MCP clients for all servers"""
        conda_env = os.getenv("CONDA_DEFAULT_ENV", "mcp-md")
        
        servers = [
            "structure_server",
            "genesis_server",
            "complex_server",
            "ligand_server",
            "assembly_server",
            "export_server",
            "qc_min_server",
        ]
        
        clients = {}
        for server_name in servers:
            logger.info(f"Creating MCP client for {server_name}")
            clients[server_name] = MCPClient(
                lambda s=server_name: stdio_client(
                    StdioServerParameters(
                        command="conda",
                        args=["run", "-n", conda_env, "python", "-m", f"servers.{s}"],
                    )
                )
            )
        
        return clients
    
    def run_interactive(self):
        """Run interactive CLI"""
        logger.info("Starting interactive CLI")
        
        # Enter MCP client context
        clients_list = list(self.mcp_clients.values())
        
        with clients_list[0], clients_list[1], clients_list[2], clients_list[3], \
             clients_list[4], clients_list[5], clients_list[6]:
            
            # Collect all tools from all servers
            tools = []
            for server_name, client in self.mcp_clients.items():
                logger.info(f"Loading tools from {server_name}")
                server_tools = client.list_tools_sync()
                tools.extend(server_tools)
                logger.info(f"  Loaded {len(server_tools)} tools")
            
            logger.info(f"Total tools available: {len(tools)}")
            
            # Create Strands agent
            agent = Agent(model=self.model, tools=tools)
            
            # Print banner
            print("=" * 60)
            print("MD Workflow Agent")
            print(f"LM Studio: {self.lm_studio_url}")
            print(f"Model: {self.model_id}")
            print(f"Run directory: {self.run_dir}")
            print("=" * 60)
            print("Type your query or 'exit' to quit")
            print()
            
            # Interactive loop
            while True:
                try:
                    user_input = prompt("> ")
                    
                    if user_input.lower().strip() in ["exit", "quit", "q"]:
                        print("Goodbye!")
                        break
                    
                    if not user_input.strip():
                        continue
                    
                    # Log user query
                    self.decision_logger.log_decision(
                        step="user_query",
                        tool="user_input",
                        params={"query": user_input},
                        reason="User provided natural language query"
                    )
                    
                    # Agent processes query
                    response = agent(user_input)
                    print()
                    
                except KeyboardInterrupt:
                    print("\nInterrupted. Type 'exit' to quit.")
                    continue
                except Exception as e:
                    logger.error(f"Error: {e}")
                    print(f"Error: {e}")
                    continue
    
    def run_workflow_from_plan(self, plan_file: str):
        """Run workflow from a plan file
        
        Args:
            plan_file: Path to workflow plan (Markdown or JSON)
        """
        logger.info(f"Running workflow from plan: {plan_file}")
        
        # TODO: Implement plan parsing and execution
        # 1. Parse plan file (Markdown/JSON)
        # 2. Extract workflow steps
        # 3. Map steps to MCP tools
        # 4. Execute in sequence with validation
        # 5. Log all decisions
        
        raise NotImplementedError("Workflow from plan not yet implemented")
    
    def suggest_workflow(self, query: str) -> Dict[str, Any]:
        """Suggest workflow steps for a given query
        
        Args:
            query: Natural language description of desired MD system
        
        Returns:
            Dict with suggested workflow steps
        """
        logger.info(f"Suggesting workflow for: {query}")
        
        # Use LLM to suggest workflow based on skeleton
        system_prompt = f"""You are an expert in Molecular Dynamics simulations.
Given a user query, suggest a workflow using the following skeleton:
{WORKFLOW_SKELETON}

For each step, specify:
1. Which MCP tool to use
2. Required parameters
3. Rationale

Return as structured JSON.
"""
        
        # TODO: Implement LLM-based workflow suggestion
        raise NotImplementedError("Workflow suggestion not yet implemented")

