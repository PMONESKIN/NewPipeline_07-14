"""
Module 01 — Step 1: Structure Retrieval
=======================================
Downloads co-crystal PDB structures for all targets.

Usage:
    python3 modules/01_targets/fetch_structures.py --run-dir path
"""

import argparse
import logging
import sys
from pathlib import Path

import requests
import yaml
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"


def fetch_pdb(pdb_id: str, output_dir: Path) -> Path:
    pdb_id = pdb_id.upper()
    output_path = output_dir / f"{pdb_id}.pdb"

    if output_path.exists():
        log.info(f"{pdb_id}: already downloaded, skipping.")
        return output_path

    url = RCSB_URL.format(pdb_id=pdb_id)
    log.info(f"Downloading PDB {pdb_id} from RCSB...")
    response = requests.get(url, timeout=30)
    if response.status_code != 200:
        raise RuntimeError(f"Failed to download PDB {pdb_id} (HTTP {response.status_code})")

    output_path.write_text(response.text)
    log.info(f"Saved → {output_path}")
    return output_path


def inspect_chains(pdb_path: Path) -> dict:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))

    chains = {}
    for model in structure:
        for chain in model:
            residues = [r for r in chain if is_aa(r, standard=True)]
            if not residues:
                continue
            seq = "".join(protein_letters_3to1.get(r.get_resname(), "X") for r in residues)
            chains[chain.id] = {"residue_count": len(residues), "sequence": seq}

    log.info(f"\nChain summary for {pdb_path.name}:")
    for chain_id, info in chains.items():
        preview = info["sequence"][:20] + ("..." if len(info["sequence"]) > 20 else "")
        log.info(f"  Chain {chain_id}: {info['residue_count']} residues | {preview}")

    return chains


def run(run_dir: str) -> dict:
    sys.path.insert(0, str(ROOT))
    from modules.run_manager import RunManager

    rm = RunManager(run_dir=run_dir)
    config = rm.config
    targets = config.get("targets", [])
    fetched = {}

    for target in targets:
        name = target["name"]
        pdb_id = target.get("pdb_id")
        if not pdb_id:
            log.error(f"{name}: no pdb_id specified. Skipping.")
            continue

        log.info(f"\n--- {name} ---")
        try:
            pdb_path = fetch_pdb(pdb_id, rm.structures_dir)
            inspect_chains(pdb_path)

            receptor_chain = target.get("receptor_chain")
            ligand_chain = target.get("ligand_chain")
            chain_info = inspect_chains(pdb_path) if False else {}  # already printed

            fetched[name] = str(pdb_path)
        except Exception as e:
            log.error(f"{name}: fetch failed — {e}")

    log.info(f"\nDone. {len(fetched)} structure(s) fetched.")
    return fetched


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    args = parser.parse_args()
    run(run_dir=args.run_dir)
