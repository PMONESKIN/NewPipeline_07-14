"""
Module 03 — Step 3: MD Stability Simulations (OpenMM)
=====================================================
Runs molecular dynamics simulations on top HADDOCK3-ranked candidates
to validate binding stability. Uses CHARMM36 force field with
physiological conditions (310K, 150mM NaCl).

Selects candidates based on HADDOCK3 scores (fast or full).
Uses folded peptide PDBs from Module 02.

Features:
  - CHARMM36 force field
  - Platform cascade: CUDA > OpenCL > Metal > CPU
  - Checkpoint/resume support
  - Graceful SIGINT handling
  - Full energy logging

Usage:
    python3 modules/03_docking/md_stability.py --run-dir path [--n 10] [--ns 50]

Inputs:
    {run_dir}/candidates/candidate_pool.json (with haddock_fast_score or haddock_full_score)
    {run_dir}/docking/receptor/{PDB_ID}_receptor.pdb
    data/candidates/folded_structures/ (peptide PDBs)

Outputs:
    {run_dir}/md_results/{candidate_id}/
      - solvated.pdb
      - trajectory.dcd
      - energy.csv
      - checkpoint.chk
      - summary.json
"""

import argparse
import json
import logging
import signal
import sys
from pathlib import Path

import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]


def select_platform():
    """Select the best available compute platform."""
    from openmm import Platform

    preferences = ["CUDA", "OpenCL", "Metal", "CPU"]
    for name in preferences:
        try:
            platform = Platform.getPlatformByName(name)
            log.info(f"  Platform: {name}")
            return platform
        except Exception:
            continue

    log.warning("  No GPU platform found, using CPU")
    return Platform.getPlatformByName("CPU")


def build_complex(receptor_pdb: Path, ligand_pdb: Path, output_path: Path) -> Path:
    """Combine receptor and ligand PDBs into a single complex."""
    receptor_lines = []
    ligand_lines = []

    with open(receptor_pdb) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                receptor_lines.append(line)

    with open(ligand_pdb) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                ligand_lines.append(line)

    with open(output_path, "w") as f:
        for line in receptor_lines:
            f.write(line)
        f.write("TER\n")
        for line in ligand_lines:
            f.write(line)
        f.write("TER\nEND\n")

    return output_path


def run_md_simulation(
    complex_pdb: Path,
    output_dir: Path,
    ns: float = 50.0,
    temperature_K: float = 310.15,
    ionic_strength: float = 0.15,
    checkpoint_interval_ps: float = 500,
    report_interval_ps: float = 50,
    resume: bool = False,
) -> dict:
    """
    Run MD simulation on a protein-peptide complex.

    Returns:
        Dict with energy stats and completion status
    """
    from openmm import unit, LangevinMiddleIntegrator, MonteCarloBarostat
    from openmm.app import (
        ForceField, Modeller, Simulation, PDBFile,
        DCDReporter, StateDataReporter, CheckpointReporter,
    )
    from pdbfixer import PDBFixer

    output_dir.mkdir(parents=True, exist_ok=True)

    solvated_pdb = output_dir / "solvated.pdb"
    traj_path = output_dir / "trajectory.dcd"
    energy_path = output_dir / "energy.csv"
    checkpoint_path = output_dir / "checkpoint.chk"
    final_pdb = output_dir / "final.pdb"

    dt = 0.002  # ps
    total_steps = int(ns * 1000 / dt)
    checkpoint_steps = int(checkpoint_interval_ps / dt)
    report_steps = int(report_interval_ps / dt)

    forcefield = ForceField("charmm36.xml", "charmm36/water.xml")

    # Load or build solvated system
    if resume and solvated_pdb.exists():
        log.info("  Loading saved solvated system for resume...")
        pdb_loaded = PDBFile(str(solvated_pdb))
        topology = pdb_loaded.topology
        positions = pdb_loaded.positions
    else:
        log.info("  Fixing structure...")
        fixer = PDBFixer(filename=str(complex_pdb))
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=7.4)
        fixer.removeHeterogens(keepWater=False)

        log.info("  Setting up system (CHARMM36, 150mM NaCl)...")
        modeller = Modeller(fixer.topology, fixer.positions)
        modeller.addSolvent(
            forcefield,
            model="tip3p",
            padding=1.0 * unit.nanometers,
            ionicStrength=ionic_strength * unit.molar,
        )

        topology = modeller.topology
        positions = modeller.positions

        with open(solvated_pdb, "w") as f:
            PDBFile.writeFile(topology, positions, f)

        log.info(f"  Total atoms: {topology.getNumAtoms()}")

    # Create system
    system = forcefield.createSystem(
        topology,
        nonbondedMethod=2,  # PME
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=2,  # HBonds
    )
    system.addForce(
        MonteCarloBarostat(1.0 * unit.atmospheres, temperature_K * unit.kelvin)
    )

    integrator = LangevinMiddleIntegrator(
        temperature_K * unit.kelvin,
        1.0 / unit.picoseconds,
        dt * unit.picoseconds,
    )

    platform = select_platform()
    simulation = Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(positions)

    # Resume or minimize + equilibrate
    if resume and checkpoint_path.exists():
        log.info("  Resuming from checkpoint...")
        simulation.loadCheckpoint(str(checkpoint_path))
        current_step = simulation.context.getState().getStepCount()
        remaining_steps = total_steps - current_step
        log.info(f"  Resumed at step {current_step} ({current_step * dt / 1000:.1f} ns)")

        if remaining_steps <= 0:
            log.info("  Simulation already complete!")
            return {"completed": True, "steps_run": current_step, "total_steps": total_steps}

        # New files for resumed segment (avoids DCD corruption)
        traj_path = output_dir / f"trajectory_resumed_{current_step}.dcd"
        energy_path = output_dir / f"energy_resumed_{current_step}.csv"
    else:
        remaining_steps = total_steps

        log.info("  Minimizing energy...")
        simulation.minimizeEnergy()

        log.info("  Equilibration (NVT 100ps + NPT 100ps)...")
        simulation.context.setVelocitiesToTemperature(temperature_K * unit.kelvin)
        simulation.step(int(100 / dt))  # NVT
        simulation.step(int(100 / dt))  # NPT

    # Reporters
    simulation.reporters.append(DCDReporter(str(traj_path), report_steps))
    simulation.reporters.append(StateDataReporter(
        str(energy_path), report_steps,
        step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        totalEnergy=True, temperature=True, volume=True, speed=True,
        separator=",",
    ))
    simulation.reporters.append(StateDataReporter(
        sys.stdout, report_steps * 4,
        step=True, time=True, temperature=True, speed=True,
        remainingTime=True, totalSteps=total_steps,
    ))
    simulation.reporters.append(CheckpointReporter(str(checkpoint_path), checkpoint_steps))

    # Graceful shutdown
    stop_requested = [False]

    def signal_handler(sig, frame):
        log.info("\n  SIGINT received — saving checkpoint...")
        stop_requested[0] = True

    original_handler = signal.signal(signal.SIGINT, signal_handler)

    # Production MD
    log.info(f"  Production MD ({ns} ns, {remaining_steps} steps)...")
    steps_done = 0
    chunk_size = checkpoint_steps

    while steps_done < remaining_steps:
        if stop_requested[0]:
            simulation.saveCheckpoint(str(checkpoint_path))
            log.info(f"  Stopped at step {steps_done}/{remaining_steps}")
            break

        batch = min(chunk_size, remaining_steps - steps_done)
        simulation.step(batch)
        steps_done += batch

    signal.signal(signal.SIGINT, original_handler)

    # Save final structure
    positions = simulation.context.getState(getPositions=True).getPositions()
    with open(final_pdb, "w") as f:
        PDBFile.writeFile(simulation.topology, positions, f)

    # Analyze
    results = {
        "completed": not stop_requested[0],
        "steps_run": steps_done,
        "total_steps": total_steps,
        "ns_completed": round(steps_done * dt / 1000, 2),
    }

    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(results, indent=2))

    return results


def select_top_candidates(candidates: list[dict], n: int = 10) -> list[dict]:
    """Select top N candidates by best available HADDOCK score."""
    # Prefer full score, fall back to fast score
    for score_field in ["haddock_full_score", "haddock_fast_score"]:
        scored = [c for c in candidates if c.get(score_field) is not None]
        if scored:
            scored.sort(key=lambda c: c[score_field])  # lower = better
            log.info(f"  Selecting top {n} by {score_field}")
            return scored[:n]

    return []


def run(run_dir: str, n: int = 10, ns: float = 50.0, resume: bool = False) -> list[dict]:
    """
    Run MD stability simulations for top N HADDOCK3-ranked candidates.
    """
    sys.path.insert(0, str(ROOT))
    from modules.run_manager import RunManager

    rm = RunManager(run_dir=run_dir)
    config = rm.config
    candidates = rm.load_candidates()

    md_cfg = config.get("md_stability", {})
    temperature = md_cfg.get("temperature_K", 310.15)
    ionic = md_cfg.get("ionic_strength_molar", 0.15)
    checkpoint_ps = md_cfg.get("checkpoint_interval_ps", 500)
    report_ps = md_cfg.get("report_interval_ps", 50)

    top = select_top_candidates(candidates, n)

    if not top:
        log.error("No HADDOCK3-scored candidates. Run haddock3_dock.py first.")
        return candidates

    log.info(f"Running MD stability for top {len(top)} candidates ({ns} ns each)...")

    # Get receptor PDB
    target_cfg = config["targets"][0]
    interface = rm.load_interface(target_cfg["name"])
    receptor_pdb = rm.docking_dir / "receptor" / f"{interface['pdb_id']}_receptor.pdb"
    folded_dir = rm.folded_dir

    if not receptor_pdb.exists():
        log.error(f"Receptor PDB not found: {receptor_pdb}")
        log.error("Run haddock3_dock.py first (it extracts the receptor).")
        return candidates

    for i, candidate in enumerate(top):
        cid = candidate["id"]
        score_field = "haddock_full_score" if candidate.get("haddock_full_score") is not None else "haddock_fast_score"
        score = candidate[score_field]

        log.info(f"\n[{i+1}/{len(top)}] {cid} ({score_field}={score:.1f})")
        log.info(f"  Sequence: {candidate['sequence']}")

        ligand_pdb = folded_dir / f"{cid}.pdb"
        if not ligand_pdb.exists():
            log.warning(f"  Peptide PDB not found: {ligand_pdb}")
            continue

        try:
            cand_md_dir = rm.md_dir / cid
            complex_pdb = cand_md_dir / "complex.pdb"
            cand_md_dir.mkdir(parents=True, exist_ok=True)
            build_complex(receptor_pdb, ligand_pdb, complex_pdb)

            results = run_md_simulation(
                complex_pdb=complex_pdb,
                output_dir=cand_md_dir,
                ns=ns,
                temperature_K=temperature,
                ionic_strength=ionic,
                checkpoint_interval_ps=checkpoint_ps,
                report_interval_ps=report_ps,
                resume=resume,
            )

            log.info(f"  Completed: {results['ns_completed']} ns")
            log.info(f"  Output: {cand_md_dir}")

        except Exception as e:
            log.error(f"  MD failed: {e}")

    rm.update_status("md_complete")
    log.info(f"\nMD stability simulations complete.")

    return candidates


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True, help="Run directory")
    parser.add_argument("--n", type=int, default=10, help="Number of top candidates")
    parser.add_argument("--ns", type=float, default=50.0, help="Simulation length (ns)")
    parser.add_argument("--resume", action="store_true", help="Resume from checkpoint")
    args = parser.parse_args()
    run(run_dir=args.run_dir, n=args.n, ns=args.ns, resume=args.resume)
