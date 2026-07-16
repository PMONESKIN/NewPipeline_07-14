"""
Module 01 — Step 2: Interface Analysis
=======================================
Analyzes the binding interface between receptor and ligand chains
in a co-crystal PDB structure. Produces a JSON summary used by
all downstream modules.

If active_residues are specified in config.yaml, validates they exist.
If not specified, auto-detects interface residues using a distance cutoff.

Usage:
    python modules/01_targets/analyze_interface.py

Inputs:
    data/structures/{PDB_ID}.pdb
    config.yaml (target definitions)

Outputs:
    data/processed/{target_name}_interface.json
"""

import json
import logging
from pathlib import Path

import numpy as np
import yaml
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]

INTERFACE_DISTANCE_CUTOFF = 4.5  # Angstroms


def load_config(config_path: str = None) -> dict:
    path = Path(config_path) if config_path else ROOT / "config.yaml"
    with open(path) as f:
        return yaml.safe_load(f)


def get_chain_residues(structure, chain_id: str) -> list:
    """Get all standard amino acid residues from a chain."""
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                return [r for r in chain if is_aa(r, standard=True)]
    return []


def get_chain_sequence(residues: list) -> str:
    """Extract one-letter amino acid sequence from a list of residues."""
    return "".join(
        protein_letters_3to1.get(r.get_resname(), "X") for r in residues
    )


def get_residue_numbers(residues: list) -> list[int]:
    """Get PDB residue numbers from a list of residues."""
    return [r.get_id()[1] for r in residues]


def validate_active_residues(
    active_residues: list[int], receptor_residues: list
) -> list[int]:
    """
    Validate that all user-specified active residues exist in the receptor chain.
    Returns the validated list. Raises ValueError if any are missing.
    """
    existing_numbers = set(get_residue_numbers(receptor_residues))
    missing = [r for r in active_residues if r not in existing_numbers]

    if missing:
        raise ValueError(
            f"Active residues not found in receptor chain: {missing}. "
            f"Available residue range: {min(existing_numbers)}-{max(existing_numbers)}"
        )

    log.info(f"  Validated {len(active_residues)} active residues from config")
    return active_residues


def detect_interface_residues(
    receptor_residues: list,
    ligand_residues: list,
    cutoff: float = INTERFACE_DISTANCE_CUTOFF,
) -> list[int]:
    """
    Auto-detect receptor residues at the binding interface.
    Finds all receptor residues with any heavy atom within cutoff
    distance of any ligand heavy atom.

    Args:
        receptor_residues: list of receptor Bio.PDB Residue objects
        ligand_residues: list of ligand Bio.PDB Residue objects
        cutoff: distance threshold in Angstroms

    Returns:
        Sorted list of receptor residue numbers at the interface
    """
    # Collect all receptor heavy atoms
    receptor_atoms = []
    for res in receptor_residues:
        for atom in res:
            if atom.element != "H":
                receptor_atoms.append(atom)

    if not receptor_atoms:
        raise ValueError("No heavy atoms found in receptor chain")

    # Build neighbor search from receptor atoms
    ns = NeighborSearch(receptor_atoms)

    # Find receptor residues near any ligand atom
    interface_residues = set()
    for res in ligand_residues:
        for atom in res:
            if atom.element == "H":
                continue
            nearby = ns.search(atom.get_vector().get_array(), cutoff, level="R")
            for neighbor_res in nearby:
                if is_aa(neighbor_res, standard=True):
                    interface_residues.add(neighbor_res.get_id()[1])

    result = sorted(interface_residues)
    log.info(f"  Auto-detected {len(result)} interface residues (cutoff={cutoff} Å)")
    log.info(f"  Residues: {result}")
    return result


def compute_binding_box(ligand_residues: list, padding: float = 10.0) -> dict:
    """
    Compute a binding box centered on the ligand chain.
    Used by AutoDock Vina for docking.

    Args:
        ligand_residues: list of ligand Bio.PDB Residue objects
        padding: extra space around the ligand in each direction (Angstroms)

    Returns:
        Dict with center_x/y/z and size_x/y/z
    """
    coords = []
    for res in ligand_residues:
        for atom in res:
            if atom.element != "H":
                coords.append(atom.get_vector().get_array())

    if not coords:
        raise ValueError("No heavy atoms found in ligand chain")

    coords = np.array(coords)
    center = coords.mean(axis=0)
    span = coords.max(axis=0) - coords.min(axis=0)
    size = span + 2 * padding

    return {
        "center_x": round(float(center[0]), 2),
        "center_y": round(float(center[1]), 2),
        "center_z": round(float(center[2]), 2),
        "size_x": round(float(size[0]), 1),
        "size_y": round(float(size[1]), 1),
        "size_z": round(float(size[2]), 1),
    }


def analyze_target(target: dict, config: dict) -> dict:
    """
    Run interface analysis for a single target.

    Returns:
        Interface summary dict to be saved as JSON
    """
    name = target["name"]
    pdb_id = target["pdb_id"].upper()
    receptor_chain = target["receptor_chain"]
    ligand_chain = target["ligand_chain"]

    structures_dir = Path(config["outputs"]["structures"])
    pdb_path = structures_dir / f"{pdb_id}.pdb"

    if not pdb_path.exists():
        raise FileNotFoundError(
            f"{pdb_path} not found. Run fetch_structures.py first."
        )

    log.info(f"\n--- Analyzing {name} ({pdb_id}) ---")
    log.info(f"  Receptor: chain {receptor_chain} | Ligand: chain {ligand_chain}")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, str(pdb_path))

    # Get residues for both chains
    receptor_residues = get_chain_residues(structure, receptor_chain)
    ligand_residues = get_chain_residues(structure, ligand_chain)

    if not receptor_residues:
        raise ValueError(f"No residues found in receptor chain {receptor_chain}")
    if not ligand_residues:
        raise ValueError(f"No residues found in ligand chain {ligand_chain}")

    log.info(f"  Receptor: {len(receptor_residues)} residues")
    log.info(f"  Ligand: {len(ligand_residues)} residues")

    # Active residues: validate from config or auto-detect
    config_residues = target.get("active_residues")
    if config_residues:
        active_residues = validate_active_residues(config_residues, receptor_residues)
        residues_source = "config"
    else:
        active_residues = detect_interface_residues(receptor_residues, ligand_residues)
        residues_source = "auto_detected"

    # Ligand sequence
    ligand_sequence = get_chain_sequence(ligand_residues)
    log.info(f"  Ligand sequence: {ligand_sequence}")

    # Design chain length: from config or auto-detect
    design_chain_length = target.get("design_chain_length") or len(ligand_residues)
    log.info(f"  Design chain length: {design_chain_length}")

    # Binding box: from config or auto-compute from ligand
    config_box = target.get("binding_box", {})
    if config_box and config_box.get("center_x") is not None:
        binding_box = {
            "center_x": config_box["center_x"],
            "center_y": config_box["center_y"],
            "center_z": config_box["center_z"],
            "size_x": config_box.get("size_x", 40.0),
            "size_y": config_box.get("size_y", 40.0),
            "size_z": config_box.get("size_z", 40.0),
        }
        log.info(f"  Binding box: from config")
    else:
        binding_box = compute_binding_box(ligand_residues)
        log.info(
            f"  Binding box: auto-computed | center=({binding_box['center_x']}, "
            f"{binding_box['center_y']}, {binding_box['center_z']})"
        )

    # Fixed motif
    fixed_motif = target.get("fixed_motif")
    if fixed_motif:
        if fixed_motif in ligand_sequence:
            log.info(f"  Fixed motif '{fixed_motif}' found in ligand sequence")
        else:
            log.warning(
                f"  Fixed motif '{fixed_motif}' NOT found in ligand sequence "
                f"'{ligand_sequence}' — ProteinMPNN will redesign all positions"
            )
            fixed_motif = None

    # Build result
    result = {
        "target_name": name,
        "pdb_id": pdb_id,
        "pdb_path": str(pdb_path),
        "receptor_chain": receptor_chain,
        "ligand_chain": ligand_chain,
        "receptor_residue_count": len(receptor_residues),
        "ligand_residue_count": len(ligand_residues),
        "active_residues": active_residues,
        "active_residues_source": residues_source,
        "ligand_sequence": ligand_sequence,
        "design_chain_length": design_chain_length,
        "binding_box": binding_box,
        "fixed_motif": fixed_motif,
        "mpnn_temperatures": target.get("mpnn_temperatures", [0.1, 0.2, 0.3]),
        "seed_sequences": target.get("seed_sequences", []),
    }

    return result


def run(config_path: str = None) -> list[dict]:
    """
    Analyze interfaces for all targets in config.yaml.

    Returns:
        List of interface result dicts
    """
    config = load_config(config_path)
    output_dir = Path(config["outputs"]["processed"])
    output_dir.mkdir(parents=True, exist_ok=True)

    targets = config.get("targets", [])
    if not targets:
        log.error("No targets found in config.yaml.")
        return []

    results = []

    for target in targets:
        name = target.get("name", "UNKNOWN")
        try:
            result = analyze_target(target, config)

            # Save per-target JSON
            safe_name = name.lower().replace("/", "_").replace(" ", "_")
            output_path = output_dir / f"{safe_name}_interface.json"
            with open(output_path, "w") as f:
                json.dump(result, f, indent=2)
            log.info(f"  Saved → {output_path}")

            results.append(result)

        except Exception as e:
            log.error(f"{name}: analysis failed — {e}")

    log.info(f"\nDone. {len(results)} target(s) analyzed.")
    return results


if __name__ == "__main__":
    run()
