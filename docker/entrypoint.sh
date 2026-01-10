#!/bin/bash
# MDZen MCP Server Entrypoint
# Starts the unified MCP server with Streamable HTTP transport

set -e

# Default values
HOST="${MDZEN_HOST:-0.0.0.0}"
PORT="${MDZEN_PORT:-3000}"
LOG_LEVEL="${MDZEN_LOG_LEVEL:-INFO}"

echo "=========================================="
echo "MDZen MCP Server"
echo "=========================================="
echo "Host: $HOST"
echo "Port: $PORT"
echo "Log Level: $LOG_LEVEL"
echo "Working Directory: $(pwd)"
echo "=========================================="

# Export log level for Python servers
export MDZEN_LOG_LEVEL="$LOG_LEVEL"

# Change to app directory for imports
cd /app

# Start unified MCP server with HTTP transport
exec python -u servers/unified_server.py --http --host "$HOST" --port "$PORT"
