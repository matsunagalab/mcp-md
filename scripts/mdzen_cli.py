#!/usr/bin/env python
"""
MDZen CLI - Command-line interface for MD simulation setup.

This CLI wraps the FastMCP server functions for direct command-line usage.
Use this instead of MCP when you have a local conda environment with dependencies.

Usage:
    python scripts/mdzen_cli.py download --pdb-id 1AKE
    python scripts/mdzen_cli.py prepare --structure-file 1ake.cif --chains A
    python scripts/mdzen_cli.py solvate --pdb-file merged.pdb
    python scripts/mdzen_cli.py topology --pdb-file solvated.pdb --box-dims box.json
"""

import argparse
import asyncio
import json
import os
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from common.utils import setup_logger

logger = setup_logger(__name__)


def download_structure(args):
    """Download structure from PDB."""
    from servers.research_server import download_structure as _download

    async def run():
        result = await _download(
            pdb_id=args.pdb_id,
            format=args.format,
            output_dir=args.output_dir,
        )
        return result

    result = asyncio.run(run())
    print(json.dumps(result, indent=2, default=str))
    return result


def prepare_complex(args):
    """Prepare protein-ligand complex."""
    from servers.structure_server import prepare_complex as _prepare

    # Parse chains if provided
    chains = args.chains.split(",") if args.chains else None

    # Parse ligand_smiles if provided
    ligand_smiles = None
    if args.ligand_smiles:
        ligand_smiles = json.loads(args.ligand_smiles)

    # Parse include_types if provided
    include_types = args.include_types.split(",") if args.include_types else None

    # In FastMCP 2.x, decorated functions need .fn to access underlying function
    func = _prepare.fn if hasattr(_prepare, "fn") else _prepare

    result = func(
        structure_file=args.structure_file,
        output_dir=args.output_dir,
        select_chains=chains,
        ph=args.ph,
        cap_termini=args.cap_termini,
        process_proteins=not args.no_proteins,
        process_ligands=not args.no_ligands,
        run_parameterization=not args.no_parameterization,
        ligand_smiles=ligand_smiles,
        include_types=include_types,
        optimize_ligands=not args.no_optimize,
        charge_method=args.charge_method,
        atom_type=args.atom_type,
    )
    print(json.dumps(result, indent=2, default=str))
    return result


def solvate_structure(args):
    """Solvate structure with water box."""
    from servers.solvation_server import solvate_structure as _solvate

    func = _solvate.fn if hasattr(_solvate, "fn") else _solvate

    result = func(
        pdb_file=args.pdb_file,
        output_dir=args.output_dir,
        output_name=args.output_name,
        dist=args.distance,
        cubic=not args.rectangular,
        salt=not args.no_salt,
        salt_c=args.cation,
        salt_a=args.anion,
        saltcon=args.salt_concentration,
    )
    print(json.dumps(result, indent=2, default=str))
    return result


def build_topology(args):
    """Build Amber topology files."""
    from servers.amber_server import build_amber_system as _build

    # Parse box_dimensions if provided
    box_dimensions = None
    if args.box_dimensions:
        if os.path.isfile(args.box_dimensions):
            with open(args.box_dimensions) as f:
                box_dimensions = json.load(f)
        else:
            box_dimensions = json.loads(args.box_dimensions)

    # Parse ligand_params if provided
    ligand_params = None
    if args.ligand_params:
        if os.path.isfile(args.ligand_params):
            with open(args.ligand_params) as f:
                ligand_params = json.load(f)
        else:
            ligand_params = json.loads(args.ligand_params)

    func = _build.fn if hasattr(_build, "fn") else _build

    result = func(
        pdb_file=args.pdb_file,
        ligand_params=ligand_params,
        box_dimensions=box_dimensions,
        forcefield=args.forcefield,
        water_model=args.water_model,
        is_membrane=args.membrane,
        output_name=args.output_name,
        output_dir=args.output_dir,
    )
    print(json.dumps(result, indent=2, default=str))
    return result


def run_simulation(args):
    """Run MD simulation."""
    from servers.md_simulation_server import run_md_simulation as _run

    func = _run.fn if hasattr(_run, "fn") else _run

    result = func(
        prmtop_file=args.prmtop,
        inpcrd_file=args.inpcrd,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        simulation_time_ns=args.time,
        temperature_k=args.temperature,
        pressure_bar=args.pressure if not args.nvt else None,
        timestep_fs=args.timestep,
        report_interval_ps=args.report_interval,
        trajectory_interval_ps=args.trajectory_interval,
        platform=args.platform,
        minimize_first=not args.no_minimize,
        minimize_steps=args.minimize_steps,
    )
    print(json.dumps(result, indent=2, default=str))
    return result


def main():
    parser = argparse.ArgumentParser(
        description="MDZen CLI - Molecular Dynamics Setup Tools",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download structure
  python scripts/mdzen_cli.py download --pdb-id 1AKE --output-dir ./workdir

  # Prepare complex
  python scripts/mdzen_cli.py prepare --structure-file ./workdir/1ake.cif --chains A --output-dir ./workdir

  # Solvate
  python scripts/mdzen_cli.py solvate --pdb-file ./workdir/merged.pdb --output-dir ./workdir

  # Build topology
  python scripts/mdzen_cli.py topology --pdb-file ./workdir/solvated.pdb --box-dimensions ./workdir/box.json --output-dir ./workdir

  # Run simulation
  python scripts/mdzen_cli.py simulate --prmtop ./workdir/system.parm7 --inpcrd ./workdir/system.rst7 --time 0.1
""",
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Download command
    download_parser = subparsers.add_parser("download", help="Download structure from PDB")
    download_parser.add_argument("--pdb-id", required=True, help="PDB ID (e.g., 1AKE)")
    download_parser.add_argument("--format", default="cif", choices=["pdb", "cif"], help="Output format")
    download_parser.add_argument("--output-dir", default="./workdir", help="Output directory")
    download_parser.set_defaults(func=download_structure)

    # Prepare command
    prepare_parser = subparsers.add_parser("prepare", help="Prepare protein-ligand complex")
    prepare_parser.add_argument("--structure-file", required=True, help="Input structure file (PDB/CIF)")
    prepare_parser.add_argument("--chains", help="Chains to include (comma-separated, e.g., A,B)")
    prepare_parser.add_argument("--ph", type=float, default=7.4, help="Target pH for protonation")
    prepare_parser.add_argument("--cap-termini", action="store_true", help="Cap termini with ACE/NME")
    prepare_parser.add_argument("--no-proteins", action="store_true", help="Skip protein processing")
    prepare_parser.add_argument("--no-ligands", action="store_true", help="Skip ligand processing")
    prepare_parser.add_argument("--no-parameterization", action="store_true", help="Skip ligand parameterization")
    prepare_parser.add_argument("--ligand-smiles", help="JSON dict of ligand SMILES (e.g., '{\"ATP\": \"...\"}')")
    prepare_parser.add_argument("--include-types", help="Types to include (e.g., protein,ligand)")
    prepare_parser.add_argument("--no-optimize", action="store_true", help="Skip ligand geometry optimization")
    prepare_parser.add_argument("--charge-method", default="bcc", help="Charge method for antechamber")
    prepare_parser.add_argument("--atom-type", default="gaff2", help="Atom type for antechamber")
    prepare_parser.add_argument("--output-dir", default="./workdir", help="Output directory")
    prepare_parser.set_defaults(func=prepare_complex)

    # Solvate command
    solvate_parser = subparsers.add_parser("solvate", help="Solvate structure with water box")
    solvate_parser.add_argument("--pdb-file", required=True, help="Input PDB file")
    solvate_parser.add_argument("--distance", type=float, default=15.0, help="Distance to box edge (Å)")
    solvate_parser.add_argument("--rectangular", action="store_true", help="Use rectangular box (default: cubic)")
    solvate_parser.add_argument("--no-salt", action="store_true", help="Skip adding salt")
    solvate_parser.add_argument("--salt-concentration", type=float, default=0.15, help="Salt concentration (M)")
    solvate_parser.add_argument("--cation", default="Na+", help="Cation type")
    solvate_parser.add_argument("--anion", default="Cl-", help="Anion type")
    solvate_parser.add_argument("--output-name", default="solvated", help="Output file name (without extension)")
    solvate_parser.add_argument("--output-dir", default="./workdir", help="Output directory")
    solvate_parser.set_defaults(func=solvate_structure)

    # Topology command
    topology_parser = subparsers.add_parser("topology", help="Build Amber topology files")
    topology_parser.add_argument("--pdb-file", required=True, help="Input PDB file (solvated)")
    topology_parser.add_argument("--box-dimensions", help="Box dimensions (JSON string or file path)")
    topology_parser.add_argument("--ligand-params", help="Ligand parameters (JSON string or file path)")
    topology_parser.add_argument("--forcefield", default="ff14SB", help="Force field")
    topology_parser.add_argument("--water-model", default="tip3p", help="Water model")
    topology_parser.add_argument("--membrane", action="store_true", help="Membrane system")
    topology_parser.add_argument("--output-name", default="system", help="Output file name prefix")
    topology_parser.add_argument("--output-dir", default="./workdir", help="Output directory")
    topology_parser.set_defaults(func=build_topology)

    # Simulate command
    simulate_parser = subparsers.add_parser("simulate", help="Run MD simulation")
    simulate_parser.add_argument("--prmtop", required=True, help="Amber prmtop file")
    simulate_parser.add_argument("--inpcrd", required=True, help="Amber inpcrd/rst7 file")
    simulate_parser.add_argument("--time", type=float, default=0.1, help="Simulation time (ns)")
    simulate_parser.add_argument("--temperature", type=float, default=300.0, help="Temperature (K)")
    simulate_parser.add_argument("--pressure", type=float, default=1.0, help="Pressure (bar)")
    simulate_parser.add_argument("--nvt", action="store_true", help="Use NVT instead of NPT")
    simulate_parser.add_argument("--timestep", type=float, default=2.0, help="Timestep (fs)")
    simulate_parser.add_argument("--report-interval", type=float, default=1.0, help="Report interval (ps)")
    simulate_parser.add_argument("--trajectory-interval", type=float, default=10.0, help="Trajectory save interval (ps)")
    simulate_parser.add_argument("--platform", default="CPU", choices=["CPU", "CUDA", "OpenCL"], help="OpenMM platform")
    simulate_parser.add_argument("--no-minimize", action="store_true", help="Skip energy minimization")
    simulate_parser.add_argument("--minimize-steps", type=int, default=1000, help="Minimization steps")
    simulate_parser.add_argument("--output-prefix", default="md", help="Output file prefix")
    simulate_parser.add_argument("--output-dir", default="./workdir", help="Output directory")
    simulate_parser.set_defaults(func=run_simulation)

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    try:
        args.func(args)
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
