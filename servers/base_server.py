"""
Base MCP server class with common functionality.
"""

import logging
from pathlib import Path
from typing import Optional, Any
from mcp.server import Server
from mcp.types import Tool, TextContent

logger = logging.getLogger(__name__)


class BaseMCPServer:
    """Base class for MCP servers
    
    Provides common functionality for all MCP-MD servers.
    """
    
    def __init__(self, name: str, version: str = "0.1.0"):
        """Initialize base server
        
        Args:
            name: Server name
            version: Server version
        """
        self.name = name
        self.version = version
        self.server = Server(name)
        self.working_dir = Path.cwd() / "output"
        self.working_dir.mkdir(exist_ok=True, parents=True)
        
        logger.info(f"{self.name} v{self.version} initialized")
        logger.info(f"Working directory: {self.working_dir}")
    
    def setup_handlers(self):
        """Set up MCP protocol handlers
        
        Should be overridden by subclasses to register tools.
        """
        raise NotImplementedError("Subclasses must implement setup_handlers()")
    
    def get_server(self) -> Server:
        """Get MCP server instance
        
        Returns:
            MCP Server instance
        """
        return self.server
    
    def create_tool_response(
        self,
        content: Any,
        is_error: bool = False
    ) -> list[TextContent]:
        """Create tool response
        
        Args:
            content: Response content (will be converted to string)
            is_error: Whether this is an error response
        
        Returns:
            List of TextContent for MCP response
        """
        text = str(content)
        return [TextContent(type="text", text=text)]
    
    def get_output_path(self, filename: str) -> Path:
        """Get output file path
        
        Args:
            filename: Output filename
        
        Returns:
            Full path to output file
        """
        return self.working_dir / filename
    
    def validate_input_file(self, file_path: str) -> bool:
        """Validate input file exists
        
        Args:
            file_path: Path to input file
        
        Returns:
            True if file exists and is readable
        """
        path = Path(file_path)
        if not path.is_file():
            logger.error(f"Input file not found: {file_path}")
            return False
        return True
    
    async def run(self):
        """Run the MCP server
        
        Should be called after setup_handlers().
        """
        from mcp.server.stdio import stdio_server
        
        async with stdio_server() as (read_stream, write_stream):
            await self.server.run(
                read_stream,
                write_stream,
                self.server.create_initialization_options()
            )

