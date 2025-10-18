"""
MCP-MD - Molecular Dynamics Input File Generation Agent

Main entry point for the MCP-MD workflow system.
"""

import typer
import asyncio
from pathlib import Path
from rich.console import Console
from rich.table import Table

from core.planner import MDWorkflowPlanner
from core.workflow import WorkflowEngine

app = typer.Typer(help="MD Input File Generation Agent with Boltz-2, AmberTools, and OpenMM")
console = Console()


@app.command()
def plan(
    query: str = typer.Argument(..., help="Natural language workflow description"),
    output: str = typer.Option("workflow_plan.md", help="Output plan file")
):
    """Create workflow plan from natural language query"""
    console.print(f"[bold blue]Creating workflow plan...[/bold blue]")
    console.print(f"Query: {query}")
    
    planner = MDWorkflowPlanner()
    plan = planner.plan_from_query(query)
    
    planner.save_plan(plan, output)
    
    console.print(f"[bold green]✓ Plan saved to: {output}[/bold green]")
    console.print(f"System type: {plan.system_type.value}")
    console.print(f"Steps: {len(plan.steps)}")


@app.command()
def list_servers():
    """List available MCP servers"""
    table = Table(title="Available MCP Servers")
    table.add_column("Server", style="cyan")
    table.add_column("Description", style="green")
    
    servers = [
        ("structure_server", "PDB retrieval, Boltz-2 prediction, structure cleaning"),
        ("ligand_server", "RDKit 3D generation, AmberTools GAFF2 parameterization"),
        ("docking_server", "smina-based ligand docking"),
        ("assembly_server", "tleap system building (solvation + ions)"),
        ("protocol_server", "OpenMM MD script generation"),
        ("export_server", "Format conversion and packaging"),
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
