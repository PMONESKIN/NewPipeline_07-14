"""
Module 01 — Step 2: Interface Analysis
=======================================
Analyzes the binding interface between receptor and ligand chains.
Produces a JSON summary used by all downstream modules.

Usage:
    python3 modules/01_targets/analyze_interface.py --run-dir path
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import yaml
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
INTERFACE_DISTANCE_CUTOFF = 4.5


def get_chain_residues(structure, chain_id):
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                return [r for r in chain if is_aa(r, standard=True)]
    return []


def get_chain_sequence(residues):
    return "".join(protein_letters_3to1.get(r.get_resname(), "X") for r in residues)


def detect_interface_residues(receptor_residues, ligand_residues, cutoff=INTERFACE_DISTANCE_CUTOFF):
    receptor_atoms = [a for r in receptor_residues for a in r if a.element != "H"]
    if not receptor_atoms:
        raise ValueError("No heavy atoms in receptor chain")

    ns = NeighborSearch(receptor_atoms)
    interface = set()
    for res in ligand_residues:
        for atom in res:
            if atom.element == "H":
                continue
            nearby = ns.search(atom.get_vector().get_array(), cutoff, level="R")
            for nr in nearby:
                if is_aa(nr, standard=True):
                    interface.add(nr.get_id()[1])

    return sorted(interface)


def validate_active_residues(active_residues, receptor_residues):
    existing = set(r.get_id()[1] for r in receptor_residues)
    missing = [r for r in active_residues if r not in existing]
    if missing:
        raise ValueError(f"Active residues not found: {missing}")
    return active_residues


def compute_binding_box(ligand_residues, padding=10.0):
    coords = []
    for res in ligand_residues:
        for atom in res:
            if atom.element != "H":
                coords.append(atom.get_vector().get_array())
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


def analyze_target(target, rm):
    name = target["name"]
    pdb_id = target["pdb_id"].upper()
    receptor_chain = target["receptor_chain"]
    ligand_chain = target["ligand_chain"]

    pdb_path = rm.structures_dir / f"{pdb_id}.pdb"
    if not pdb_path.exists():
        raise FileNotFoundError(f"{pdb_path} not found. Run fetch_structures.py first.")

    log.info(f"\n--- Analyzing {name} ({pdb_id}) ---")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, str(pdb_path))

    receptor_residues = get_chain_residues(structure, receptor_chain)
    ligand_residues = get_chain_residues(structure, ligand_chain)

    if not receptor_residues:
        raise ValueError(f"No residues in receptor chain {receptor_chain}")
    if not ligand_residues:
        raise ValueError(f"No residues in ligand chain {ligand_chain}")

    log.info(f"  Receptor: {len(receptor_residues)} residues | Ligand: {len(ligand_residues)} residues")

    config_residues = target.get("active_residues")
    if config_residues:
        active_residues = validate_active_residues(config_residues, receptor_residues)
        residues_source = "config"
        log.info(f"  Validated {len(active_residues)} active residues from config")
    else:
        active_residues = detect_interface_residues(receptor_residues, ligand_residues)
        residues_source = "auto_detected"
        log.info(f"  Auto-detected {len(active_residues)} interface residues")

    ligand_sequence = get_chain_sequence(ligand_residues)
    design_chain_length = target.get("design_chain_length") or len(ligand_residues)
    binding_box = compute_binding_box(ligand_residues)

    fixed_motif = target.get("fixed_motif")
    if fixed_motif and fixed_motif not in ligand_sequence:
        log.warning(f"  Fixed motif '{fixed_motif}' not found in ligand — ignoring")
        fixed_motif = None

    return {
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


def run(run_dir: str) -> list[dict]:
    sys.path.insert(0, str(ROOT))
    from modules.run_manager import RunManager

    rm = RunManager(run_dir=run_dir)
    config = rm.config
    results = []

    for target in config.get("targets", []):
        name = target.get("name", "UNKNOWN")
        try:
            result = analyze_target(target, rm)
            rm.save_interface(name, result)
            log.info(f"  Saved → {rm.processed_dir}")
            results.append(result)
        except Exception as e:
            log.error(f"{name}: analysis failed — {e}")

    log.info(f"\nDone. {len(results)} target(s) analyzed.")
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    args = parser.parse_args()
    run(run_dir=args.run_dir)
