"""
Export Server - Format conversion and packaging.

Provides MCP tools for:
- Amber format export
- Package system files
"""

import logging
import json
import shutil
from pathlib import Path
from mcp.types import Tool, TextContent

from .base_server import BaseMCPServer
from core.utils import setup_logger

logger = setup_logger(__name__)


class ExportServer(BaseMCPServer):
    """MCP Server for export and packaging"""
    
    def __init__(self):
        super().__init__("export_server", "0.1.0")
        self.setup_handlers()
    
    def setup_handlers(self):
        """Register MCP tool handlers"""
        
        @self.server.list_tools()
        async def list_tools() -> list[Tool]:
            return [
                Tool(
                    name="package_system",
                    description="Package all system files into archive",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "system_dir": {"type": "string"},
                            "output_name": {"type": "string", "default": "md_system"}
                        },
                        "required": ["system_dir"]
                    }
                ),
            ]
        
        @self.server.call_tool()
        async def call_tool(name: str, arguments: dict) -> list[TextContent]:
            try:
                if name == "package_system":
                    result = await self.package_system(
                        system_dir=arguments["system_dir"],
                        output_name=arguments.get("output_name", "md_system")
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
    
    async def package_system(
        self,
        system_dir: str,
        output_name: str = "md_system"
    ) -> dict:
        """Package system files"""
        logger.info(f"Packaging system: {system_dir}")
        
        system_path = Path(system_dir)
        if not system_path.exists():
            raise FileNotFoundError(f"System directory not found: {system_dir}")
        
        output_dir = self.get_output_path("packages")
        output_dir.mkdir(exist_ok=True)
        
        # Create archive
        archive_path = output_dir / output_name
        shutil.make_archive(str(archive_path), 'zip', system_path)
        
        final_path = archive_path.with_suffix('.zip')
        
        logger.info(f"Package created: {final_path}")
        
        return {
            "archive_path": str(final_path),
            "format": "zip",
            "source_dir": str(system_path)
        }


async def main():
    """Run export server"""
    server = ExportServer()
    await server.run()


if __name__ == "__main__":
    import asyncio
    asyncio.run(main())

