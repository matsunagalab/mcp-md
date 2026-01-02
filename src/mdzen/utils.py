"""Common utilities for MCP-MD ADK.

Contains all utility functions needed by the ADK agents.
"""

import json
import logging
import warnings
from contextlib import contextmanager
from datetime import datetime


# =============================================================================
# DATE AND TIME UTILITIES
# =============================================================================


def get_today_str() -> str:
    """Get today's date as a formatted string.

    Returns:
        Date string in YYYY-MM-DD format
    """
    return datetime.now().strftime("%Y-%m-%d")


def format_duration(seconds: float) -> str:
    """Format duration in human-readable format.

    Args:
        seconds: Duration in seconds

    Returns:
        Human-readable duration string (e.g., "2m 30s")
    """
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes = int(seconds // 60)
    remaining_seconds = seconds % 60
    if minutes < 60:
        return f"{minutes}m {remaining_seconds:.0f}s"
    hours = minutes // 60
    remaining_minutes = minutes % 60
    return f"{hours}h {remaining_minutes}m"


# =============================================================================
# ADK STATE HELPERS
# =============================================================================


def safe_dict(value, default: dict | None = None) -> dict:
    """Safely convert a value to dict, handling JSON strings.

    ADK state may serialize values as JSON strings. This function handles
    both dict and string inputs.

    Args:
        value: Value to convert (dict, str, or other)
        default: Default value if conversion fails

    Returns:
        Dictionary representation of the value
    """
    if default is None:
        default = {}

    if value is None:
        return default
    if isinstance(value, dict):
        return value
    if isinstance(value, str):
        if not value:
            return default
        try:
            parsed = json.loads(value)
            if isinstance(parsed, dict):
                return parsed
            return default
        except json.JSONDecodeError:
            return default
    return default


def safe_list(value, default: list | None = None) -> list:
    """Safely convert a value to list, handling JSON strings.

    ADK state may serialize values as JSON strings. This function handles
    both list and string inputs.

    Args:
        value: Value to convert (list, str, or other)
        default: Default value if conversion fails

    Returns:
        List representation of the value
    """
    if default is None:
        default = []

    if value is None:
        return default
    if isinstance(value, list):
        return value
    if isinstance(value, str):
        if not value:
            return default
        try:
            parsed = json.loads(value)
            if isinstance(parsed, list):
                return parsed
            return default
        except json.JSONDecodeError:
            return default
    return default


@contextmanager
def suppress_adk_unknown_agent_warnings():
    """Context manager to suppress ADK 'unknown agent' warnings.

    The ADK runner sometimes emits "Event from an unknown agent" warnings
    when using nested SequentialAgents. This context manager temporarily
    suppresses these warnings.

    Usage:
        with suppress_adk_unknown_agent_warnings():
            async for event in runner.run_async(...):
                ...
    """

    class UnknownAgentFilter(logging.Filter):
        """Filter out 'Event from an unknown agent' warnings."""

        def filter(self, record):
            return "unknown agent" not in record.getMessage()

    unknown_filter = UnknownAgentFilter()
    original_settings = {}

    # Target ADK loggers (both naming conventions)
    logger_names = [
        "google_adk.runners",
        "google_adk",
        "google.adk.runners",
        "google.adk",
    ]

    # Apply filter to all relevant loggers
    for logger_name in logger_names:
        logger = logging.getLogger(logger_name)
        original_settings[logger_name] = {
            "level": logger.level,
            "propagate": logger.propagate,
        }
        logger.addFilter(unknown_filter)
        logger.setLevel(logging.CRITICAL)  # Only show CRITICAL
        logger.propagate = False  # Don't propagate to parent loggers

    # Also apply to root logger handlers
    root_logger = logging.getLogger()
    for handler in root_logger.handlers:
        handler.addFilter(unknown_filter)

    # Suppress warnings module messages
    warnings.filterwarnings("ignore", message=".*unknown agent.*")

    try:
        yield
    finally:
        # Restore original settings and remove filters
        for logger_name in logger_names:
            logger = logging.getLogger(logger_name)
            logger.removeFilter(unknown_filter)
            if logger_name in original_settings:
                logger.setLevel(original_settings[logger_name]["level"])
                logger.propagate = original_settings[logger_name]["propagate"]
        for handler in root_logger.handlers:
            handler.removeFilter(unknown_filter)


__all__ = [
    # Date/time
    "get_today_str",
    "format_duration",
    # ADK state helpers
    "safe_dict",
    "safe_list",
    "suppress_adk_unknown_agent_warnings",
]
