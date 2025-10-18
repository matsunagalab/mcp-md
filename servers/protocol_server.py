"""
Protocol Server - OpenMM MD script generation.

Provides MCP tools for:
- OpenMM minimization scripts
- OpenMM equilibration scripts
- OpenMM production scripts
- Complete MD workflow generation
"""

import logging
import json
from pathlib import Path
from mcp.types import Tool, TextContent

from .base_server import BaseMCPServer
from tools.openmm_wrapper import OpenMMWrapper
from core.utils import setup_logger

logger = setup_logger(__name__)


class ProtocolServer(BaseMCPServer):
    """MCP Server for MD protocol generation"""
    
    def __init__(self):
        super().__init__("protocol_server", "0.1.0")
        self.openmm = OpenMMWrapper()
        self.setup_handlers()
    
    def setup_handlers(self):
        """Register MCP tool handlers"""
        
        @self.server.list_tools()
        async def list_tools() -> list[Tool]:
            return [
                Tool(
                    name="create_openmm_workflow",
                    description="Create complete OpenMM MD workflow scripts",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "prmtop": {"type": "string"},
                            "inpcrd": {"type": "string"},
                            "protocol": {"type": "string", "default": "standard"}
                        },
                        "required": ["prmtop", "inpcrd"]
                    }
                ),
            ]
        
        @self.server.call_tool()
        async def call_tool(name: str, arguments: dict) -> list[TextContent]:
            try:
                if name == "create_openmm_workflow":
                    result = await self.create_openmm_workflow(
                        prmtop=arguments["prmtop"],
                        inpcrd=arguments["inpcrd"],
                        protocol=arguments.get("protocol", "standard")
                    )
                else:
                    raise ValueError(f"Unknown tool: {name}")
                
                return self.create_tool_response(json.dumps(result, indent=2))
            
            except Exception as e:
                logger.error(f"Tool {name} failed: {e}")
                return self.create_tool_response(
                    json.dumps({"error": str(e)}),
                    is_error=True
                )
    
    async def create_openmm_workflow(
        self,
        prmtop: str,
        inpcrd: str,
        protocol: str = "standard"
    ) -> dict:
        """Create OpenMM workflow scripts"""
        logger.info(f"Creating OpenMM workflow ({protocol})")
        
        output_dir = self.get_output_path("openmm_workflow")
        
        result = self.openmm.create_workflow(
            prmtop=prmtop,
            inpcrd=inpcrd,
            output_dir=output_dir,
            protocol=protocol
        )
        
        return result


async def main():
    """Run protocol server"""
    server = ProtocolServer()
    await server.run()


if __name__ == "__main__":
    import asyncio
    asyncio.run(main())

