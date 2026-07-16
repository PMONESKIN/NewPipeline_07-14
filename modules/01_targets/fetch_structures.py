"""
Module 01 — Step 1: Structure Retrieval
=======================================
Downloads co-crystal PDB structures for all targets defined in config.yaml.

Each target must specify a pdb_id of a co-crystal structure containing
both the receptor and ligand (peptide) chains in the bound conformation.

Usage:
    python modules/01_targets/fetch_structures.py

Outputs:
    data/structures/{PDB_ID}.pdb
"""

import logging
from pathlib import Path

import requests
import yaml
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"


def load_config(config_path: str = None) -> dict:
    path = Path(config_path) if config_path else ROOT / "config.yaml"
    with open(path) as f:
        return yaml.safe_load(f)


def fetch_pdb(pdb_id: str, output_dir: Path) -> Path:
    """
    Download a PDB file from RCSB.

    Args:
        pdb_id: 4-character PDB ID (e.g., "2FLU")
        output_dir: Directory to save the file

    Returns:
        Path to the downloaded PDB file
    """
    pdb_id = pdb_id.upper()
    output_path = output_dir / f"{pdb_id}.pdb"

    if output_path.exists():
        log.info(f"{pdb_id}: already downloaded, skipping.")
        return output_path

    url = RCSB_URL.format(pdb_id=pdb_id)
    log.info(f"Downloading PDB {pdb_id} from RCSB...")

    response = requests.get(url, timeout=30)
    if response.status_code != 200:
        raise RuntimeError(
            f"Failed to download PDB {pdb_id} (HTTP {response.status_code}). "
            f"Check the PDB ID at https://www.rcsb.org/structure/{pdb_id}"
        )

    output_path.write_text(response.text)
    log.info(f"Saved → {output_path}")
    return output_path


def inspect_chains(pdb_path: Path) -> dict:
    """
    Print chain summary for a PDB file so the user can verify
    which chain is the receptor vs the peptide.

    Returns:
        Dict mapping chain_id -> {residue_count, sequence}
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))

    chains = {}
    for model in structure:
        for chain in model:
            residues = [r for r in chain if is_aa(r, standard=True)]
            if not residues:
                continue
            seq = "".join(
                protein_letters_3to1.get(r.get_resname(), "X") for r in residues
            )
            chains[chain.id] = {
                "residue_count": len(residues),
                "sequence": seq,
            }

    log.info(f"\nChain summary for {pdb_path.name}:")
    for chain_id, info in chains.items():
        preview = info["sequence"][:20]
        suffix = "..." if len(info["sequence"]) > 20 else ""
        log.info(
            f"  Chain {chain_id}: {info['residue_count']} residues | "
            f"{preview}{suffix}"
        )

    return chains


def run(config_path: str = None) -> dict:
    """
    Fetch PDB structures for all targets in config.yaml.

    Returns:
        Dict mapping target_name -> pdb_path
    """
    config = load_config(config_path)
    output_dir = Path(config["outputs"]["structures"])
    output_dir.mkdir(parents=True, exist_ok=True)

    targets = config.get("targets", [])
    if not targets:
        log.error("No targets found in config.yaml.")
        return {}

    fetched = {}

    for target in targets:
        name = target["name"]
        pdb_id = target.get("pdb_id")

        if not pdb_id:
            log.error(f"{name}: no pdb_id specified. Skipping.")
            continue

        log.info(f"\n--- {name} ---")

        try:
            pdb_path = fetch_pdb(pdb_id, output_dir)
            chain_info = inspect_chains(pdb_path)

            # Validate that the specified chains exist
            receptor_chain = target.get("receptor_chain")
            ligand_chain = target.get("ligand_chain")

            if receptor_chain and receptor_chain not in chain_info:
                log.error(
                    f"  receptor_chain '{receptor_chain}' not found in {pdb_id}. "
                    f"Available chains: {list(chain_info.keys())}"
                )
            if ligand_chain and ligand_chain not in chain_info:
                log.error(
                    f"  ligand_chain '{ligand_chain}' not found in {pdb_id}. "
                    f"Available chains: {list(chain_info.keys())}"
                )

            fetched[name] = pdb_path

        except Exception as e:
            log.error(f"{name}: fetch failed — {e}")

    log.info(f"\nDone. {len(fetched)} structure(s) saved to {output_dir}")
    return fetched


if __name__ == "__main__":
    run()
