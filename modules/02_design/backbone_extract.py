"""
Module 02 — Backbone Extraction
================================
Extracts the binding region backbone from a co-crystal PDB at a
specified length, then energy minimizes it with the receptor present.

This gives ProteinMPNN an optimized backbone template for the
target peptide length.

Usage:
    python3 modules/02_design/backbone_extract.py --run-dir path [--length 8]
"""

import argparse
import json
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]


def extract_binding_region(pdb_path: Path, ligand_chain: str, receptor_chain: str,
                           target_length: int, fixed_motif: str = None,
                           output_path: Path = None) -> Path:
    """
    Extract a backbone region of target_length centered on the binding motif.
    Keeps the full receptor chain intact.
    """
    from Bio.PDB import PDBParser, PDBIO, Select
    from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))

    # Get ligand residues and sequence
    ligand_residues = []
    for model in structure:
        for chain in model:
            if chain.id == ligand_chain:
                ligand_residues = [r for r in chain if is_aa(r, standard=True)]

    full_length = len(ligand_residues)
    if target_length >= full_length:
        log.info(f"  Target length {target_length} >= chain length {full_length}, using full chain")
        if output_path:
            io = PDBIO()
            io.set_structure(structure)
            io.save(str(output_path))
        return output_path or pdb_path

    # Find the motif center if specified
    seq = "".join(protein_letters_3to1.get(r.get_resname(), "X") for r in ligand_residues)
    residue_ids = [r.get_id()[1] for r in ligand_residues]

    if fixed_motif and fixed_motif in seq:
        motif_start = seq.index(fixed_motif)
        motif_center = motif_start + len(fixed_motif) // 2
    else:
        # Center on middle of chain
        motif_center = full_length // 2

    # Calculate extraction window centered on motif
    half = target_length // 2
    start_idx = max(0, motif_center - half)
    end_idx = start_idx + target_length

    # Adjust if we go past the end
    if end_idx > full_length:
        end_idx = full_length
        start_idx = end_idx - target_length

    keep_res_ids = set(residue_ids[start_idx:end_idx])
    extracted_seq = seq[start_idx:end_idx]

    log.info(f"  Full chain: {seq} ({full_length} aa)")
    log.info(f"  Extracted:  {extracted_seq} ({target_length} aa, positions {start_idx+1}-{end_idx})")

    if fixed_motif:
        if fixed_motif in extracted_seq:
            log.info(f"  Motif '{fixed_motif}' preserved in extraction")
        else:
            log.warning(f"  Motif '{fixed_motif}' NOT in extracted region!")

    # Write new PDB with only the extracted region + full receptor
    class RegionSelect(Select):
        def accept_residue(self, residue):
            chain = residue.get_parent()
            if chain.id == ligand_chain:
                return residue.get_id()[1] in keep_res_ids
            elif chain.id == receptor_chain:
                return residue.get_id()[0] == " "  # skip heterogens
            return False

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path), RegionSelect())

    log.info(f"  Saved extracted complex → {output_path.name}")
    return output_path


def energy_minimize(pdb_path: Path, output_path: Path, steps: int = 500) -> Path:
    """
    Energy minimize the extracted backbone with receptor present.
    Uses OpenMM if available, otherwise skips minimization.
    """
    try:
        from openmm.app import PDBFile, ForceField, Modeller, Simulation
        from openmm import LangevinMiddleIntegrator
        from openmm import unit
        from pdbfixer import PDBFixer

        log.info("  Energy minimizing extracted backbone...")

        fixer = PDBFixer(filename=str(pdb_path))
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=7.4)
        fixer.removeHeterogens(keepWater=False)

        forcefield = ForceField("charmm36.xml", "charmm36/water.xml")
        modeller = Modeller(fixer.topology, fixer.positions)

        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=0,  # NoCutoff for in-vacuum minimization
            constraints=2,  # HBonds
        )

        integrator = LangevinMiddleIntegrator(
            300 * unit.kelvin, 1.0 / unit.picoseconds, 0.002 * unit.picoseconds
        )

        simulation = Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)
        simulation.minimizeEnergy(maxIterations=steps)

        positions = simulation.context.getState(getPositions=True).getPositions()
        with open(output_path, "w") as f:
            PDBFile.writeFile(simulation.topology, positions, f)

        log.info(f"  Minimized → {output_path.name}")
        return output_path

    except ImportError:
        log.warning("  OpenMM not available — skipping energy minimization")
        import shutil
        shutil.copy(pdb_path, output_path)
        return output_path


def run(run_dir: str, length: int = None) -> dict:
    """Extract and minimize backbone for each target."""
    sys.path.insert(0, str(ROOT))
    from modules.run_manager import RunManager

    rm = RunManager(run_dir=run_dir)
    config = rm.config

    results = {}

    for target_cfg in config.get("targets", []):
        name = target_cfg["name"]
        target_length = length or target_cfg.get("design_chain_length", 16)

        log.info(f"\n--- Extracting backbone for {name} (length={target_length}) ---")

        try:
            interface = rm.load_interface(name)
            pdb_path = rm.structures_dir / f"{interface['pdb_id']}.pdb"

            # Extract region
            extracted_pdb = rm.processed_dir / f"{interface['pdb_id']}_extracted_{target_length}aa.pdb"
            extract_binding_region(
                pdb_path=pdb_path,
                ligand_chain=interface["ligand_chain"],
                receptor_chain=interface["receptor_chain"],
                target_length=target_length,
                fixed_motif=interface.get("fixed_motif"),
                output_path=extracted_pdb,
            )

            # Energy minimize
            minimized_pdb = rm.processed_dir / f"{interface['pdb_id']}_minimized_{target_length}aa.pdb"
            energy_minimize(extracted_pdb, minimized_pdb)

            # Update interface JSON with backbone info
            interface["extracted_pdb"] = str(extracted_pdb)
            interface["minimized_pdb"] = str(minimized_pdb)
            interface["design_chain_length"] = target_length
            rm.save_interface(name, interface)

            results[name] = str(minimized_pdb)

        except Exception as e:
            log.error(f"{name}: backbone extraction failed — {e}")

    log.info(f"\nDone. {len(results)} backbone(s) extracted.")
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--length", type=int, default=None, help="Target peptide length")
    args = parser.parse_args()
    run(run_dir=args.run_dir, length=args.length)
