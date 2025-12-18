"""CLI commands for MCP-MD ADK.

This module provides Typer CLI commands for running the ADK-based
MD setup workflow.
"""

import asyncio
from datetime import datetime
from typing import Optional

import typer
from rich.console import Console
from prompt_toolkit import PromptSession

app = typer.Typer(
    name="adk",
    help="Google ADK-based MD setup commands",
)
console = Console()

# Constants
APP_NAME = "mcp_md_adk"
DEFAULT_USER = "default"


@app.command()
def run(
    request: Optional[str] = typer.Argument(
        None,
        help="MD setup request (optional, prompts if not provided)",
    ),
    batch: bool = typer.Option(
        False,
        "--batch",
        "-b",
        help="Run in batch mode (no human-in-the-loop)",
    ),
    checkpoint_dir: str = typer.Option(
        "checkpoints",
        help="Directory for session persistence",
    ),
    session_id: Optional[str] = typer.Option(
        None,
        help="Session ID for resuming",
    ),
):
    """Run MD setup using Google ADK.

    Examples:
        # Interactive mode
        python main.py adk run "Setup MD for PDB 1AKE"

        # Batch mode
        python main.py adk run --batch "Setup MD for PDB 1AKE, 1ns at 300K"

        # Resume session
        python main.py adk run --session-id md_session_xxx
    """
    asyncio.run(_run_async(request, batch, checkpoint_dir, session_id))


async def _run_async(
    request: Optional[str],
    batch: bool,
    checkpoint_dir: str,
    session_id: Optional[str],
):
    """Async implementation of the run command."""

    from mcp_md_adk.state.session_manager import (
        create_session_service,
    )

    # Create session service
    session_service = create_session_service(
        checkpoint_dir=checkpoint_dir,
        in_memory=batch,  # Use in-memory for batch mode
    )

    # Generate or use provided session ID
    if session_id is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        session_id = f"md_session_{timestamp}"

    console.print("=" * 60)
    console.print("[bold cyan]MCP-MD ADK Mode[/bold cyan]")
    console.print(f"Session ID: {session_id}")
    console.print(f"Mode: {'Batch' if batch else 'Interactive'}")
    console.print("=" * 60)

    # Get initial request if not provided
    if not request:
        console.print("\nDescribe your MD simulation setup:")
        console.print("(e.g., 'Setup MD for PDB 1AKE in water, 1 ns at 300K')")
        request = input("\n> ").strip()

        if request.lower() in ["quit", "exit", "q"]:
            console.print("[yellow]Session ended.[/yellow]")
            return

    if batch:
        await _run_batch(session_service, session_id, request)
    else:
        await _run_interactive(session_service, session_id, request)


async def _run_batch(session_service, session_id: str, request: str):
    """Run in batch mode (no interrupts)."""
    from google.adk.runners import Runner
    from google.genai import types
    from mcp_md_adk.agents.full_agent import create_full_agent
    from mcp_md_adk.state.session_manager import (
        initialize_session_state,
        get_session_state,
    )

    # Initialize session (async)
    await initialize_session_state(
        session_service=session_service,
        app_name=APP_NAME,
        user_id=DEFAULT_USER,
        session_id=session_id,
    )

    # Create full agent and runner
    agent = create_full_agent()
    runner = Runner(
        app_name=APP_NAME,
        agent=agent,
        session_service=session_service,
    )

    console.print("[dim]Running full workflow...[/dim]\n")

    # Create message content
    user_message = types.Content(
        role="user",
        parts=[types.Part(text=request)],
    )

    # Run agent
    async for event in runner.run_async(
        user_id=DEFAULT_USER,
        session_id=session_id,
        new_message=user_message,
    ):
        if event.is_final_response():
            pass  # Result captured in session state

    # Show results
    state = await get_session_state(session_service, APP_NAME, DEFAULT_USER, session_id)

    if state.get("validation_result"):
        validation = state["validation_result"]
        if "final_report" in validation:
            console.print("\n[bold green]Workflow Complete![/bold green]")
            console.print(validation["final_report"])
        else:
            console.print("\n[bold green]Batch Complete![/bold green]")
            console.print(f"Session directory: {state.get('session_dir')}")

    # Show generated files
    outputs = state.get("outputs", {})
    if outputs:
        console.print("\n[bold]Generated Files:[/bold]")
        for key, value in outputs.items():
            if key != "session_dir":
                console.print(f"  {key}: {value}")


async def _run_interactive(session_service, session_id: str, request: str):
    """Run in interactive mode with human-in-the-loop."""
    import json
    from google.adk.runners import Runner
    from google.genai import types
    from mcp_md_adk.agents.full_agent import (
        create_clarification_only_agent,
        create_setup_validation_agent,
    )
    from mcp_md_adk.state.session_manager import (
        initialize_session_state,
        get_session_state,
    )

    # Create async prompt session
    prompt_session = PromptSession()

    async def async_prompt(message: str) -> str:
        return await prompt_session.prompt_async(message)

    def create_message(text: str) -> types.Content:
        """Create a user message content object."""
        return types.Content(role="user", parts=[types.Part(text=text)])

    # Initialize session (async)
    await initialize_session_state(
        session_service=session_service,
        app_name=APP_NAME,
        user_id=DEFAULT_USER,
        session_id=session_id,
    )

    # Phase 1: Clarification
    console.print("\n[bold]Phase 1: Clarification[/bold]")
    console.print("[dim]Analyzing your request...[/dim]\n")

    clarification_agent = create_clarification_only_agent()
    runner = Runner(
        app_name=APP_NAME,
        agent=clarification_agent,
        session_service=session_service,
    )

    # Run clarification
    async for event in runner.run_async(
        user_id=DEFAULT_USER,
        session_id=session_id,
        new_message=create_message(request),
    ):
        if event.is_final_response():
            console.print(f"[blue]Agent:[/blue] {event.content}\n")

    # Interactive clarification loop
    state = await get_session_state(session_service, APP_NAME, DEFAULT_USER, session_id)

    while True:
        simulation_brief = state.get("simulation_brief")

        if simulation_brief:
            console.print("\n[bold green]SimulationBrief Generated:[/bold green]")
            console.print(json.dumps(simulation_brief, indent=2))

            console.print("\n[yellow]Options:[/yellow]")
            console.print("  - Type 'continue' or 'yes' to proceed to Setup phase")
            console.print("  - Type 'quit' to exit")
            console.print("  - Or provide feedback to modify the brief\n")

            user_input = (await async_prompt(">> ")).strip()

            if user_input.lower() in ["quit", "exit", "q"]:
                console.print("[yellow]Session ended.[/yellow]")
                return
            elif user_input.lower() in ["continue", "yes", "y", "ok", "proceed"]:
                break
            else:
                # User wants to modify - send feedback
                async for event in runner.run_async(
                    user_id=DEFAULT_USER,
                    session_id=session_id,
                    new_message=create_message(user_input),
                ):
                    if event.is_final_response():
                        console.print(f"\n[blue]Agent:[/blue] {event.content}\n")

                state = await get_session_state(session_service, APP_NAME, DEFAULT_USER, session_id)
        else:
            # No brief yet, ask for more info
            user_input = (await async_prompt(">> ")).strip()

            if user_input.lower() in ["quit", "exit", "q"]:
                console.print("[yellow]Session ended.[/yellow]")
                return

            async for event in runner.run_async(
                user_id=DEFAULT_USER,
                session_id=session_id,
                new_message=create_message(user_input),
            ):
                if event.is_final_response():
                    console.print(f"\n[blue]Agent:[/blue] {event.content}\n")

            state = await get_session_state(session_service, APP_NAME, DEFAULT_USER, session_id)

    # Phase 2-3: Setup and Validation
    console.print("\n[bold]Phase 2-3: Setup & Validation[/bold]")
    console.print("[dim]Executing MD setup workflow... This may take a few minutes.[/dim]\n")

    setup_agent = create_setup_validation_agent()
    runner = Runner(
        app_name=APP_NAME,
        agent=setup_agent,
        session_service=session_service,
    )

    async for event in runner.run_async(
        user_id=DEFAULT_USER,
        session_id=session_id,
        new_message=create_message("continue"),
    ):
        if event.is_final_response():
            pass  # Final response handled below

    # Show results
    state = await get_session_state(session_service, APP_NAME, DEFAULT_USER, session_id)

    if state.get("validation_result"):
        validation = state["validation_result"]
        if "final_report" in validation:
            console.print("\n[bold green]Setup Complete![/bold green]")
            console.print(validation["final_report"])

    # Show generated files
    outputs = state.get("outputs", {})
    if outputs:
        console.print("\n[bold]Generated Files:[/bold]")
        for key, value in outputs.items():
            if key != "session_dir":
                console.print(f"  {key}: {value}")

    console.print(f"\n[green]Session complete! Session ID: {session_id}[/green]")
    console.print(f"[dim]Session directory: {state.get('session_dir')}[/dim]")


@app.command()
def info():
    """Show ADK implementation information."""
    console.print("[bold]MCP-MD ADK Implementation[/bold]")
    console.print()
    console.print("This is the Google ADK-based implementation of MCP-MD.")
    console.print()
    console.print("Features:")
    console.print("  - Google Agent Development Kit (ADK)")
    console.print("  - LiteLLM for Anthropic Claude models")
    console.print("  - McpToolset for MCP server integration")
    console.print("  - SequentialAgent for 3-phase workflow")
    console.print("  - SessionService for state persistence")
    console.print()
    console.print("Usage:")
    console.print("  python main.py adk run [request]")
    console.print("  python main.py adk run --batch [request]")
    console.print()
