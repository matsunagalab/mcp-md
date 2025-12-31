"""Session management for MDZen.

This module provides session service configuration and state initialization
for ADK-based workflows.
"""

import uuid
from pathlib import Path
from typing import Optional

from google.adk.sessions import InMemorySessionService, DatabaseSessionService

from mdzen.config import get_output_dir


def _generate_job_id(length: int = 8) -> str:
    """Generate unique job identifier using UUID.

    Args:
        length: Length of ID (default: 8 characters)

    Returns:
        Unique job ID string (without prefix)
    """
    return uuid.uuid4().hex[:length]


def create_session_service(
    checkpoint_dir: str = "checkpoints",
    in_memory: bool = False,
) -> InMemorySessionService | DatabaseSessionService:
    """Create a session service for state persistence.

    Args:
        checkpoint_dir: Directory for SQLite database (if persistent)
        in_memory: Use in-memory storage (no persistence)

    Returns:
        Configured session service instance
    """
    if in_memory:
        return InMemorySessionService()

    # Ensure checkpoint directory exists
    Path(checkpoint_dir).mkdir(parents=True, exist_ok=True)

    # Create database session service
    # Note: Use sqlite+aiosqlite for async compatibility
    db_path = Path(checkpoint_dir) / "adk_sessions.db"
    return DatabaseSessionService(
        db_url=f"sqlite+aiosqlite:///{db_path}"
    )


async def initialize_session_state(
    session_service,
    app_name: str,
    user_id: str,
    session_id: Optional[str] = None,
) -> str:
    """Create and initialize a new session with required state.

    Args:
        session_service: SessionService instance
        app_name: Application name
        user_id: User identifier
        session_id: Optional session ID (generated if not provided)

    Returns:
        Session ID (format: job_XXXXXXXX)
    """
    # Generate session ID if not provided (use same ID for session and directory)
    if session_id is None:
        job_id = _generate_job_id()
        session_id = f"job_{job_id}"
    else:
        # Extract job_id from session_id if it's in job_XXXXXXXX format
        job_id = session_id.replace("job_", "") if session_id.startswith("job_") else session_id

    # Create session directory using same job_id for consistency
    session_dir = create_session_directory(job_id)

    # Initialize state with required fields
    # IMPORTANT: State must be passed during create_session, not modified after
    # InMemorySessionService doesn't persist modifications to session.state
    initial_state = {
        "session_dir": session_dir,
        "completed_steps": [],
        "outputs": {"session_dir": session_dir},
        "decision_log": [],
        "simulation_brief": None,
        "compressed_setup": "",
        "validation_result": None,
    }

    # Create session with initial state (async method)
    await session_service.create_session(
        app_name=app_name,
        user_id=user_id,
        session_id=session_id,
        state=initial_state,
    )

    return session_id


def create_session_directory(job_id: Optional[str] = None) -> str:
    """Create a unique session directory for workflow outputs.

    The directory is created under the configured output directory.
    Uses consistent job_XXXXXXXX naming format.

    Args:
        job_id: Optional job ID (generated if not provided)

    Returns:
        Absolute path to the session directory
    """
    output_base = get_output_dir()
    if job_id is None:
        job_id = _generate_job_id()
    session_dir = output_base / f"job_{job_id}"
    session_dir.mkdir(parents=True, exist_ok=True)

    return str(session_dir.resolve())


async def get_session_state(
    session_service,
    app_name: str,
    user_id: str,
    session_id: str,
) -> dict:
    """Get current state from a session.

    Args:
        session_service: SessionService instance
        app_name: Application name
        user_id: User identifier
        session_id: Session ID

    Returns:
        Current session state dictionary
    """
    session = await session_service.get_session(
        app_name=app_name,
        user_id=user_id,
        session_id=session_id,
    )

    if session is None:
        return {}

    return dict(session.state)


async def update_session_state(
    session_service,
    app_name: str,
    user_id: str,
    session_id: str,
    updates: dict,
) -> None:
    """Update session state with new values.

    Args:
        session_service: SessionService instance
        app_name: Application name
        user_id: User identifier
        session_id: Session ID
        updates: Dictionary of state updates
    """
    session = await session_service.get_session(
        app_name=app_name,
        user_id=user_id,
        session_id=session_id,
    )

    if session is not None:
        for key, value in updates.items():
            session.state[key] = value
