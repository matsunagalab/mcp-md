"""
MD Simulation Server - Molecular dynamics simulation & analysis with OpenMM and MDTraj.

Provides MCP tools for:
- OpenMM MD simulation (NVT/NPT equilibration, production)
- MDTraj trajectory analysis (RMSD, RMSF, distances, hydrogen bonds, etc.)
- Energy analysis
- Secondary structure analysis
"""

import logging
import numpy as np
from pathlib import Path
from typing import Optional
from fastmcp import FastMCP

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from common.utils import setup_logger, ensure_directory

logger = setup_logger(__name__)

# Create FastMCP server
mcp = FastMCP("MD Simulation Server")

# Initialize working directory
WORKING_DIR = Path("output/md_simulation")
ensure_directory(WORKING_DIR)


@mcp.tool()
def run_md_simulation(
    prmtop_file: str,
    inpcrd_file: str,
    simulation_time_ns: float = 1.0,
    temperature_kelvin: float = 300.0,
    pressure_bar: Optional[float] = None,
    timestep_fs: float = 2.0,
    output_frequency_ps: float = 10.0,
    trajectory_format: str = "dcd",
    restraint_file: Optional[str] = None
) -> dict:
    """Run MD simulation using OpenMM
    
    Args:
        prmtop_file: Amber topology file
        inpcrd_file: Amber coordinate file
        simulation_time_ns: Simulation time in nanoseconds
        temperature_kelvin: Temperature in Kelvin
        pressure_bar: Pressure in bar (None for NVT, set for NPT)
        timestep_fs: Integration timestep in femtoseconds
        output_frequency_ps: Output frequency in picoseconds
        trajectory_format: Trajectory format (dcd or pdb)
        restraint_file: Optional file with restraint definitions
    
    Returns:
        Dict with simulation results and file paths
    """
    logger.info(f"Starting MD simulation: {simulation_time_ns}ns at {temperature_kelvin}K")
    
    try:
        from openmm.app import AmberPrmtopFile, AmberInpcrdFile, PDBFile, DCDReporter, StateDataReporter
        from openmm import LangevinMiddleIntegrator, MonteCarloBarostat
        from openmm.app import Simulation, PME, HBonds
        from openmm.unit import (
            nanometer, kelvin, picosecond, picoseconds, femtosecond, 
            femtoseconds, bar, atmosphere, kilojoule_per_mole
        )
    except ImportError:
        raise ImportError("OpenMM not installed. Install with: conda install -c conda-forge openmm")
    
    prmtop_path = Path(prmtop_file)
    inpcrd_path = Path(inpcrd_file)
    
    if not prmtop_path.is_file():
        raise FileNotFoundError(f"Topology file not found: {prmtop_file}")
    if not inpcrd_path.is_file():
        raise FileNotFoundError(f"Coordinate file not found: {inpcrd_file}")
    
    output_dir = WORKING_DIR / "simulations"
    ensure_directory(output_dir)
    
    # Load system
    logger.info("Loading Amber files")
    prmtop = AmberPrmtopFile(str(prmtop_path))
    inpcrd = AmberInpcrdFile(str(inpcrd_path))
    
    # Create system
    logger.info("Creating OpenMM system")
    system = prmtop.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=1.0*nanometer,
        constraints=HBonds
    )
    
    # Add barostat if NPT
    if pressure_bar is not None:
        barostat = MonteCarloBarostat(
            pressure_bar * bar,
            temperature_kelvin * kelvin
        )
        system.addForce(barostat)
        ensemble = "NPT"
    else:
        ensemble = "NVT"
    
    # Create integrator
    integrator = LangevinMiddleIntegrator(
        temperature_kelvin * kelvin,
        1.0 / picosecond,
        timestep_fs * femtoseconds
    )
    
    # Create simulation
    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)
    
    # Apply restraints if provided
    if restraint_file and Path(restraint_file).is_file():
        logger.info(f"Applying restraints from {restraint_file}")
        # TODO: Implement restraint file parsing
    
    # Setup reporters
    trajectory_file = output_dir / f"trajectory.{trajectory_format}"
    log_file = output_dir / "simulation.log"
    energy_file = output_dir / "energy.dat"
    
    if trajectory_format.lower() == "dcd":
        simulation.reporters.append(DCDReporter(str(trajectory_file), int(output_frequency_ps / timestep_fs * 1000)))
    else:
        # For PDB format, save less frequently
        from openmm.app import PDBReporter
        simulation.reporters.append(PDBReporter(str(trajectory_file), int(output_frequency_ps / timestep_fs * 1000)))
    
    simulation.reporters.append(StateDataReporter(
        str(energy_file),
        int(output_frequency_ps / timestep_fs * 1000),
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=(ensemble == "NPT"),
        density=(ensemble == "NPT")
    ))
    
    # Get initial energy
    state = simulation.context.getState(getEnergy=True)
    initial_energy = state.getPotentialEnergy()
    logger.info(f"Initial energy: {initial_energy}")
    
    # Run simulation
    simulation_steps = int(simulation_time_ns * 1000000 / timestep_fs)  # Convert ns to steps
    logger.info(f"Running {simulation_steps} steps ({simulation_time_ns}ns)")
    
    simulation.step(simulation_steps)
    
    # Get final energy
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    final_energy = state.getPotentialEnergy()
    logger.info(f"Final energy: {final_energy}")
    
    # Save final structure
    final_pdb = output_dir / "final_structure.pdb"
    positions = state.getPositions()
    with open(final_pdb, 'w') as f:
        PDBFile.writeFile(simulation.topology, positions, f)
    
    logger.info(f"Simulation complete. Trajectory saved: {trajectory_file}")
    
    return {
        "ensemble": ensemble,
        "simulation_time_ns": simulation_time_ns,
        "temperature_kelvin": temperature_kelvin,
        "pressure_bar": pressure_bar,
        "timestep_fs": timestep_fs,
        "initial_energy_kj_mol": float(initial_energy._value),
        "final_energy_kj_mol": float(final_energy._value),
        "trajectory_file": str(trajectory_file),
        "final_structure": str(final_pdb),
        "energy_file": str(energy_file),
        "num_steps": simulation_steps
    }


@mcp.tool()
def analyze_rmsd(
    trajectory_file: str,
    topology_file: str,
    reference_file: Optional[str] = None,
    selection: str = "protein and name CA",
    start_frame: int = 0,
    end_frame: Optional[int] = None
) -> dict:
    """Calculate RMSD using MDTraj
    
    Args:
        trajectory_file: Trajectory file (DCD or PDB)
        topology_file: Topology file (PDB or PRMTOP)
        reference_file: Reference structure (default: first frame)
        selection: Atom selection string
        start_frame: Starting frame index
        end_frame: Ending frame index (None for all)
    
    Returns:
        Dict with RMSD analysis results
    """
    logger.info(f"Calculating RMSD: {trajectory_file}")
    
    try:
        import mdtraj as mdt
    except ImportError:
        raise ImportError("MDTraj not installed. Install with: conda install -c conda-forge mdtraj")
    
    traj_path = Path(trajectory_file)
    topo_path = Path(topology_file)
    
    if not traj_path.is_file():
        raise FileNotFoundError(f"Trajectory file not found: {trajectory_file}")
    if not topo_path.is_file():
        raise FileNotFoundError(f"Topology file not found: {topology_file}")
    
    # Load trajectory
    logger.info("Loading trajectory")
    if traj_path.suffix.lower() == ".dcd":
        traj = mdt.load_dcd(str(traj_path), top=str(topo_path))
    else:
        traj = mdt.load(str(traj_path))
    
    # Apply frame selection
    if end_frame is None:
        traj = traj[start_frame:]
    else:
        traj = traj[start_frame:end_frame]
    
    # Load reference
    if reference_file and Path(reference_file).is_file():
        ref = mdt.load(str(reference_file))
    else:
        ref = traj[0]
    
    # Select atoms
    selection_indices = traj.topology.select(selection)
    ref_selection = ref.atom_slice(selection_indices)
    
    # Calculate RMSD
    rmsd = mdt.rmsd(traj, ref_selection, atom_indices=selection_indices)
    
    # Calculate statistics
    mean_rmsd = float(np.mean(rmsd))
    std_rmsd = float(np.std(rmsd))
    min_rmsd = float(np.min(rmsd))
    max_rmsd = float(np.max(rmsd))
    
    logger.info(f"RMSD: mean={mean_rmsd:.2f}Å, std={std_rmsd:.2f}Å")
    
    return {
        "mean_rmsd_angstrom": mean_rmsd,
        "std_rmsd_angstrom": std_rmsd,
        "min_rmsd_angstrom": min_rmsd,
        "max_rmsd_angstrom": max_rmsd,
        "rmsd_values": rmsd.tolist(),
        "num_frames": len(rmsd),
        "selection": selection
    }


@mcp.tool()
def analyze_rmsf(
    trajectory_file: str,
    topology_file: str,
    selection: str = "protein and name CA",
    start_frame: int = 0,
    end_frame: Optional[int] = None
) -> dict:
    """Calculate RMSF (Root Mean Square Fluctuation) using MDTraj
    
    Args:
        trajectory_file: Trajectory file (DCD or PDB)
        topology_file: Topology file (PDB or PRMTOP)
        selection: Atom selection string
        start_frame: Starting frame index
        end_frame: Ending frame index (None for all)
    
    Returns:
        Dict with RMSF analysis results
    """
    logger.info(f"Calculating RMSF: {trajectory_file}")
    
    try:
        import mdtraj as mdt
    except ImportError:
        raise ImportError("MDTraj not installed. Install with: conda install -c conda-forge mdtraj")
    
    traj_path = Path(trajectory_file)
    topo_path = Path(topology_file)
    
    if not traj_path.is_file():
        raise FileNotFoundError(f"Trajectory file not found: {trajectory_file}")
    if not topo_path.is_file():
        raise FileNotFoundError(f"Topology file not found: {topology_file}")
    
    # Load trajectory
    logger.info("Loading trajectory")
    if traj_path.suffix.lower() == ".dcd":
        traj = mdt.load_dcd(str(traj_path), top=str(topo_path))
    else:
        traj = mdt.load(str(traj_path))
    
    # Apply frame selection
    if end_frame is None:
        traj = traj[start_frame:]
    else:
        traj = traj[start_frame:end_frame]
    
    # Select atoms
    selection_indices = traj.topology.select(selection)
    traj_selection = traj.atom_slice(selection_indices)
    
    # Calculate RMSF
    rmsf = mdt.rmsf(traj_selection, traj_selection)
    
    # Get residue information
    residues = [atom.residue for atom in traj_selection.topology.atoms]
    residue_names = [f"{res.name}{res.index}" for res in residues]
    
    mean_rmsf = float(np.mean(rmsf))
    max_rmsf = float(np.max(rmsf))
    max_rmsf_residue = residue_names[int(np.argmax(rmsf))]
    
    logger.info(f"RMSF: mean={mean_rmsf:.2f}Å, max={max_rmsf:.2f}Å ({max_rmsf_residue})")
    
    return {
        "mean_rmsf_angstrom": mean_rmsf,
        "max_rmsf_angstrom": max_rmsf,
        "max_rmsf_residue": max_rmsf_residue,
        "rmsf_values": rmsf.tolist(),
        "residue_names": residue_names,
        "num_residues": len(rmsf)
    }


@mcp.tool()
def calculate_distance(
    trajectory_file: str,
    topology_file: str,
    atom1_selection: str,
    atom2_selection: str,
    start_frame: int = 0,
    end_frame: Optional[int] = None
) -> dict:
    """Calculate distance between two atom selections over trajectory
    
    Args:
        trajectory_file: Trajectory file (DCD or PDB)
        topology_file: Topology file (PDB or PRMTOP)
        atom1_selection: First atom selection string
        atom2_selection: Second atom selection string
        start_frame: Starting frame index
        end_frame: Ending frame index (None for all)
    
    Returns:
        Dict with distance analysis results
    """
    logger.info(f"Calculating distance: {trajectory_file}")
    
    try:
        import mdtraj as mdt
    except ImportError:
        raise ImportError("MDTraj not installed. Install with: conda install -c conda-forge mdtraj")
    
    traj_path = Path(trajectory_file)
    topo_path = Path(topology_file)
    
    if not traj_path.is_file():
        raise FileNotFoundError(f"Trajectory file not found: {trajectory_file}")
    if not topo_path.is_file():
        raise FileNotFoundError(f"Topology file not found: {topology_file}")
    
    # Load trajectory
    logger.info("Loading trajectory")
    if traj_path.suffix.lower() == ".dcd":
        traj = mdt.load_dcd(str(traj_path), top=str(topo_path))
    else:
        traj = mdt.load(str(traj_path))
    
    # Apply frame selection
    if end_frame is None:
        traj = traj[start_frame:]
    else:
        traj = traj[start_frame:end_frame]
    
    # Select atoms
    indices1 = traj.topology.select(atom1_selection)
    indices2 = traj.topology.select(atom2_selection)
    
    if len(indices1) == 0:
        raise ValueError(f"No atoms selected: {atom1_selection}")
    if len(indices2) == 0:
        raise ValueError(f"No atoms selected: {atom2_selection}")
    
    # Calculate distances (between centers of mass if multiple atoms)
    if len(indices1) == 1 and len(indices2) == 1:
        # Single atom to single atom
        distances = mdt.compute_distances(traj, [[indices1[0], indices2[0]]])
        distances = distances.flatten()
    else:
        # Center of mass calculation
        distances = []
        for frame in traj:
            com1 = np.mean(frame.xyz[0, indices1, :], axis=0)
            com2 = np.mean(frame.xyz[0, indices2, :], axis=0)
            dist = np.linalg.norm(com1 - com2)
            distances.append(dist)
        distances = np.array(distances)
    
    # Calculate statistics
    mean_dist = float(np.mean(distances))
    std_dist = float(np.std(distances))
    min_dist = float(np.min(distances))
    max_dist = float(np.max(distances))
    
    logger.info(f"Distance: mean={mean_dist:.2f}Å, std={std_dist:.2f}Å")
    
    return {
        "mean_distance_angstrom": mean_dist,
        "std_distance_angstrom": std_dist,
        "min_distance_angstrom": min_dist,
        "max_distance_angstrom": max_dist,
        "distances": distances.tolist(),
        "num_frames": len(distances),
        "atom1_selection": atom1_selection,
        "atom2_selection": atom2_selection
    }


@mcp.tool()
def analyze_hydrogen_bonds(
    trajectory_file: str,
    topology_file: str,
    donor_selection: str = "protein",
    acceptor_selection: str = "protein",
    distance_cutoff: float = 3.0,
    angle_cutoff: float = 120.0,
    start_frame: int = 0,
    end_frame: Optional[int] = None
) -> dict:
    """Analyze hydrogen bonds using MDTraj
    
    Args:
        trajectory_file: Trajectory file (DCD or PDB)
        topology_file: Topology file (PDB or PRMTOP)
        donor_selection: Donor atom selection
        acceptor_selection: Acceptor atom selection
        distance_cutoff: Distance cutoff in Angstroms
        angle_cutoff: Angle cutoff in degrees
        start_frame: Starting frame index
        end_frame: Ending frame index (None for all)
    
    Returns:
        Dict with hydrogen bond analysis results
    """
    logger.info(f"Analyzing hydrogen bonds: {trajectory_file}")
    
    try:
        import mdtraj as mdt
    except ImportError:
        raise ImportError("MDTraj not installed. Install with: conda install -c conda-forge mdtraj")
    
    traj_path = Path(trajectory_file)
    topo_path = Path(topology_file)
    
    if not traj_path.is_file():
        raise FileNotFoundError(f"Trajectory file not found: {trajectory_file}")
    if not topo_path.is_file():
        raise FileNotFoundError(f"Topology file not found: {topology_file}")
    
    # Load trajectory
    logger.info("Loading trajectory")
    if traj_path.suffix.lower() == ".dcd":
        traj = mdt.load_dcd(str(traj_path), top=str(topo_path))
    else:
        traj = mdt.load(str(traj_path))
    
    # Apply frame selection
    if end_frame is None:
        traj = traj[start_frame:]
    else:
        traj = traj[start_frame:end_frame]
    
    # Select atoms
    donor_indices = traj.topology.select(donor_selection)
    acceptor_indices = traj.topology.select(acceptor_selection)
    
    # Calculate hydrogen bonds
    hbonds = mdt.baker_hubbard(
        traj,
        distance_cutoff=distance_cutoff / 10.0,  # Convert to nm
        angle_cutoff=np.deg2rad(angle_cutoff)
    )
    
    # Analyze hydrogen bonds
    hbond_list = []
    for frame_idx, frame_hbonds in enumerate(hbonds):
        for hbond in frame_hbonds:
            donor_idx = hbond[0]
            h_idx = hbond[1]
            acceptor_idx = hbond[2]
            
            # Check if donor and acceptor are in selections
            if donor_idx in donor_indices and acceptor_idx in acceptor_indices:
                donor_atom = traj.topology.atom(donor_idx)
                acceptor_atom = traj.topology.atom(acceptor_idx)
                
                hbond_list.append({
                    "frame": int(frame_idx),
                    "donor": f"{donor_atom.residue.name}{donor_atom.residue.index}.{donor_atom.name}",
                    "acceptor": f"{acceptor_atom.residue.name}{acceptor_atom.residue.index}.{acceptor_atom.name}"
                })
    
    # Calculate frequency
    hbond_freq = {}
    for hbond in hbond_list:
        key = f"{hbond['donor']}-{hbond['acceptor']}"
        hbond_freq[key] = hbond_freq.get(key, 0) + 1
    
    total_hbonds = len(hbond_list)
    num_frames = len(hbonds)
    avg_hbonds_per_frame = total_hbonds / num_frames if num_frames > 0 else 0
    
    logger.info(f"Found {total_hbonds} hydrogen bonds (avg {avg_hbonds_per_frame:.1f} per frame)")
    
    return {
        "total_hbonds": total_hbonds,
        "num_frames": num_frames,
        "avg_hbonds_per_frame": avg_hbonds_per_frame,
        "hbond_frequency": hbond_freq,
        "distance_cutoff_angstrom": distance_cutoff,
        "angle_cutoff_degrees": angle_cutoff
    }


@mcp.tool()
def analyze_secondary_structure(
    trajectory_file: str,
    topology_file: str,
    selection: str = "protein",
    start_frame: int = 0,
    end_frame: Optional[int] = None
) -> dict:
    """Analyze secondary structure using MDTraj
    
    Args:
        trajectory_file: Trajectory file (DCD or PDB)
        topology_file: Topology file (PDB or PRMTOP)
        selection: Atom selection string
        start_frame: Starting frame index
        end_frame: Ending frame index (None for all)
    
    Returns:
        Dict with secondary structure analysis results
    """
    logger.info(f"Analyzing secondary structure: {trajectory_file}")
    
    try:
        import mdtraj as mdt
    except ImportError:
        raise ImportError("MDTraj not installed. Install with: conda install -c conda-forge mdtraj")
    
    traj_path = Path(trajectory_file)
    topo_path = Path(topology_file)
    
    if not traj_path.is_file():
        raise FileNotFoundError(f"Trajectory file not found: {trajectory_file}")
    if not topo_path.is_file():
        raise FileNotFoundError(f"Topology file not found: {topology_file}")
    
    # Load trajectory
    logger.info("Loading trajectory")
    if traj_path.suffix.lower() == ".dcd":
        traj = mdt.load_dcd(str(traj_path), top=str(topo_path))
    else:
        traj = mdt.load(str(traj_path))
    
    # Apply frame selection
    if end_frame is None:
        traj = traj[start_frame:]
    else:
        traj = traj[start_frame:end_frame]
    
    # Select atoms
    selection_indices = traj.topology.select(selection)
    traj_selection = traj.atom_slice(selection_indices)
    
    # Calculate secondary structure
    ss = mdt.compute_dssp(traj_selection, simplified=True)
    
    # Calculate statistics
    helix_fraction = np.mean(ss == 'H')
    sheet_fraction = np.mean(ss == 'E')
    coil_fraction = np.mean(ss == 'C')
    
    logger.info(f"Secondary structure: Helix={helix_fraction:.2%}, Sheet={sheet_fraction:.2%}, Coil={coil_fraction:.2%}")
    
    return {
        "helix_fraction": float(helix_fraction),
        "sheet_fraction": float(sheet_fraction),
        "coil_fraction": float(coil_fraction),
        "secondary_structure": ss.tolist(),
        "num_residues": ss.shape[1],
        "num_frames": ss.shape[0]
    }


@mcp.tool()
def analyze_contacts(
    trajectory_file: str,
    topology_file: str,
    group1_selection: str,
    group2_selection: str,
    cutoff: float = 5.0,
    start_frame: int = 0,
    end_frame: Optional[int] = None
) -> dict:
    """Analyze contacts between two groups using MDTraj
    
    Args:
        trajectory_file: Trajectory file (DCD or PDB)
        topology_file: Topology file (PDB or PRMTOP)
        group1_selection: First group atom selection
        group2_selection: Second group atom selection
        cutoff: Distance cutoff in Angstroms
        start_frame: Starting frame index
        end_frame: Ending frame index (None for all)
    
    Returns:
        Dict with contact analysis results
    """
    logger.info(f"Analyzing contacts: {trajectory_file}")
    
    try:
        import mdtraj as mdt
    except ImportError:
        raise ImportError("MDTraj not installed. Install with: conda install -c conda-forge mdtraj")
    
    traj_path = Path(trajectory_file)
    topo_path = Path(topology_file)
    
    if not traj_path.is_file():
        raise FileNotFoundError(f"Trajectory file not found: {trajectory_file}")
    if not topo_path.is_file():
        raise FileNotFoundError(f"Topology file not found: {topology_file}")
    
    # Load trajectory
    logger.info("Loading trajectory")
    if traj_path.suffix.lower() == ".dcd":
        traj = mdt.load_dcd(str(traj_path), top=str(topo_path))
    else:
        traj = mdt.load(str(traj_path))
    
    # Apply frame selection
    if end_frame is None:
        traj = traj[start_frame:]
    else:
        traj = traj[start_frame:end_frame]
    
    # Select atoms
    indices1 = traj.topology.select(group1_selection)
    indices2 = traj.topology.select(group2_selection)
    
    if len(indices1) == 0:
        raise ValueError(f"No atoms selected: {group1_selection}")
    if len(indices2) == 0:
        raise ValueError(f"No atoms selected: {group2_selection}")
    
    # Calculate contacts
    contacts = []
    for frame in traj:
        frame_contacts = []
        for i in indices1:
            for j in indices2:
                dist = np.linalg.norm(frame.xyz[0, i, :] - frame.xyz[0, j, :]) * 10.0  # Convert to Angstroms
                if dist < cutoff:
                    atom_i = traj.topology.atom(i)
                    atom_j = traj.topology.atom(j)
                    frame_contacts.append({
                        "atom1": f"{atom_i.residue.name}{atom_i.residue.index}.{atom_i.name}",
                        "atom2": f"{atom_j.residue.name}{atom_j.residue.index}.{atom_j.name}",
                        "distance": float(dist)
                    })
        contacts.append(frame_contacts)
    
    # Calculate statistics
    num_contacts_per_frame = [len(c) for c in contacts]
    avg_contacts = float(np.mean(num_contacts_per_frame))
    max_contacts = int(np.max(num_contacts_per_frame))
    
    logger.info(f"Contacts: avg={avg_contacts:.1f} per frame, max={max_contacts}")
    
    return {
        "avg_contacts_per_frame": avg_contacts,
        "max_contacts_per_frame": max_contacts,
        "contact_frames": contacts,
        "cutoff_angstrom": cutoff,
        "group1_selection": group1_selection,
        "group2_selection": group2_selection,
        "num_frames": len(contacts)
    }


@mcp.tool()
def analyze_energy_timeseries(
    energy_file: str
) -> dict:
    """Analyze energy timeseries from simulation log
    
    Args:
        energy_file: Energy log file from OpenMM simulation
    
    Returns:
        Dict with energy analysis results
    """
    logger.info(f"Analyzing energy timeseries: {energy_file}")
    
    energy_path = Path(energy_file)
    if not energy_path.is_file():
        raise FileNotFoundError(f"Energy file not found: {energy_file}")
    
    # Parse energy file
    import pandas as pd
    
    try:
        # Try to read as CSV/TSV
        df = pd.read_csv(energy_path, sep='\s+', comment='#')
    except:
        # Fallback: manual parsing
        data = []
        with open(energy_path, 'r') as f:
            header = None
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split()
                if len(parts) > 0:
                    if header is None:
                        header = parts
                    else:
                        if len(parts) == len(header):
                            data.append([float(p) for p in parts])
        
        if not data:
            raise ValueError("Could not parse energy file")
        df = pd.DataFrame(data, columns=header)
    
    # Extract energy columns
    energy_columns = {}
    for col in df.columns:
        col_lower = col.lower()
        if 'potential' in col_lower or 'pe' in col_lower:
            energy_columns['potential'] = col
        elif 'kinetic' in col_lower or 'ke' in col_lower:
            energy_columns['kinetic'] = col
        elif 'total' in col_lower or 'te' in col_lower:
            energy_columns['total'] = col
        elif 'temperature' in col_lower or 'temp' in col_lower:
            energy_columns['temperature'] = col
    
    # Calculate statistics
    results = {
        "num_frames": len(df)
    }
    
    for energy_type, col in energy_columns.items():
        if col in df.columns:
            values = df[col].values
            results[f"{energy_type}_mean"] = float(np.mean(values))
            results[f"{energy_type}_std"] = float(np.std(values))
            results[f"{energy_type}_min"] = float(np.min(values))
            results[f"{energy_type}_max"] = float(np.max(values))
    
    logger.info(f"Analyzed {len(df)} energy frames")
    
    return results


if __name__ == "__main__":
    mcp.run()

