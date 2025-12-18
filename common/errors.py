"""Standardized error handling for MCP tools.

This module provides a unified error format that helps LLMs
understand and recover from errors during MD workflow execution.
"""

from typing import List, Optional


def create_error_result(
    error: Exception,
    hints: Optional[List[str]] = None,
    context: Optional[dict] = None,
    recoverable: bool = True,
) -> dict:
    """Create standardized error result for MCP tools.

    Returns a dictionary that:
    1. Sets success=False
    2. Provides error_type and message for LLM understanding
    3. Includes hints for potential recovery actions
    4. Preserves context for debugging

    Args:
        error: The exception that occurred
        hints: Recovery hints for the LLM (e.g., "Check file path exists")
        context: Additional context (e.g., file paths, parameters)
        recoverable: Whether the error can be recovered from

    Returns:
        Standardized error dictionary compatible with MCP tool responses

    Example:
        >>> try:
        ...     # some operation
        ... except FileNotFoundError as e:
        ...     return create_error_result(
        ...         e,
        ...         hints=["Check that the PDB file exists", "Verify the path is absolute"],
        ...         context={"file_path": "/path/to/missing.pdb"}
        ...     )
    """
    return {
        "success": False,
        "error_type": type(error).__name__,
        "message": str(error),
        "hints": hints or [],
        "context": context or {},
        "recoverable": recoverable,
        "errors": [str(error)],
        "warnings": [],
    }


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
    return create_error_result(
        FileNotFoundError(f"{file_type} not found: {file_path}"),
        hints=[
            f"Verify the {file_type} path is correct",
            "Check that the file exists",
            "Ensure the path is absolute, not relative",
        ],
        context={"file_path": file_path, "file_type": file_type},
        recoverable=True,
    )


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


def merge_errors(base_result: dict, *additional_errors: str) -> dict:
    """Merge additional error messages into an existing result.

    Args:
        base_result: Existing result dictionary
        additional_errors: Error messages to add

    Returns:
        Updated result with merged errors
    """
    result = base_result.copy()
    errors = list(result.get("errors", []))
    errors.extend(additional_errors)
    result["errors"] = errors
    if additional_errors:
        result["success"] = False
    return result
