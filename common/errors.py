"""Standardized error handling for MCP tools.

This module provides a unified error format that helps LLMs
understand and recover from errors during MD workflow execution.
"""

from typing import Optional


def create_validation_error(
    field: str,
    message: str,
    expected: Optional[str] = None,
    actual: Optional[str] = None,
) -> dict:
    """Create error result for input validation failures.

    Args:
        field: Name of the field that failed validation
        message: Description of the validation error
        expected: What was expected (optional)
        actual: What was received (optional)

    Returns:
        Standardized validation error dictionary
    """
    hints = [f"Check the '{field}' parameter"]
    if expected:
        hints.append(f"Expected: {expected}")
    if actual:
        hints.append(f"Received: {actual}")

    return {
        "success": False,
        "error_type": "ValidationError",
        "message": f"Validation failed for '{field}': {message}",
        "hints": hints,
        "context": {"field": field, "expected": expected, "actual": actual},
        "recoverable": True,
        "errors": [f"{field}: {message}"],
        "warnings": [],
    }


def create_file_not_found_error(
    file_path: str,
    file_type: str = "file",
) -> dict:
    """Create error result for missing files.

    Args:
        file_path: Path to the missing file
        file_type: Type of file (e.g., "PDB", "topology", "trajectory")

    Returns:
        Standardized file not found error dictionary
    """
    error_msg = f"{file_type} not found: {file_path}"
    return {
        "success": False,
        "error_type": "FileNotFoundError",
        "message": error_msg,
        "hints": [
            f"Verify the {file_type} path is correct",
            "Check that the file exists",
            "Ensure the path is absolute, not relative",
        ],
        "context": {"file_path": file_path, "file_type": file_type},
        "recoverable": True,
        "errors": [error_msg],
        "warnings": [],
    }


def create_tool_not_available_error(
    tool_name: str,
    install_hint: Optional[str] = None,
) -> dict:
    """Create error result for missing external tools.

    Args:
        tool_name: Name of the missing tool (e.g., "tleap", "antechamber")
        install_hint: How to install the tool

    Returns:
        Standardized tool not available error dictionary
    """
    hints = [f"Tool '{tool_name}' is not available in PATH"]
    if install_hint:
        hints.append(f"Install hint: {install_hint}")
    else:
        hints.append("Ensure AmberTools is installed and conda environment is activated")

    return {
        "success": False,
        "error_type": "ToolNotAvailableError",
        "message": f"Required tool '{tool_name}' not found",
        "hints": hints,
        "context": {"tool_name": tool_name},
        "recoverable": False,  # Usually requires environment setup
        "errors": [f"Tool not found: {tool_name}"],
        "warnings": [],
    }
