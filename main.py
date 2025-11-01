"""
MCP-MD - Molecular Dynamics Input File Generation Agent

Main entry point for the MCP-MD workflow system.
"""

import typer
import asyncio
from pathlib import Path
from rich.console import Console
from rich.table import Table

from core.langgraph_agent import MDWorkflowAgent

app = typer.Typer(help="MD Input File Generation Agent with Boltz-2, AmberTools, and OpenMM")
console = Console()


@app.command()
def chat(
    lm_studio_url: str = typer.Option(
        "http://localhost:1234/v1",
        help="LM Studio API URL"
    ),
    model: str = typer.Option(
        "gemma-3-12b",
        help="Model ID"
    ),
    run_dir: str = typer.Option(
        None,
        help="Run directory (default: runs/<timestamp>)"
    )
):
    """Start interactive chat with MD Workflow Agent"""
    try:
        agent = MDWorkflowAgent(
            lm_studio_url=lm_studio_url,
            model_id=model,
            run_dir=Path(run_dir) if run_dir else None
        )
        agent.run_interactive()
    except ImportError as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        console.print("\nInstall required packages:")
        console.print("  uv pip install -e \".[openai]\"")
        raise typer.Exit(1)
    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        raise typer.Exit(1)


@app.command()
def list_servers():
    """List available MCP servers"""
    table = Table(title="Available MCP Servers")
    table.add_column("Server", style="cyan")
    table.add_column("Description", style="green")
    
    servers = [
        ("structure_server", "PDB retrieval and structure cleaning"),
        ("genesis_server", "Boltz-2 structure generation from FASTA"),
        ("complex_server", "Boltz-2 complex prediction + Smina refinement"),
        ("ligand_server", "RDKit 3D generation, AmberTools GAFF2 parameterization"),
        ("assembly_server", "tleap system building (solvation + ions)"),
        ("export_server", "Format conversion and packaging"),
        ("qc_min_server", "MolProbity QC checks + OpenMM minimization"),
    ]
    
    for server, desc in servers:
        table.add_row(server, desc)
    
    console.print(table)


@app.command()
def info():
    """Show system information"""
    console.print("[bold]MCP-MD: Molecular Dynamics Input File Generation Agent[/bold]")
    console.print()
    console.print("Features:")
    console.print("  • Boltz-2 structure and affinity prediction")
    console.print("  • AmberTools ligand parameterization (AM1-BCC)")
    console.print("  • smina molecular docking")
    console.print("  • OpenMM MD script generation")
    console.print("  • LM Studio LLM integration")
    console.print()
    console.print("For usage, run: [cyan]mcp-md --help[/cyan]")


def main():
    """Main entry point"""
    app()


if __name__ == "__main__":
    main()
