"""
Module 02 — Step 3: Peptide 3D Structure Generation
====================================================
Generates 3D PDB structures for each designed peptide sequence.

Currently uses PeptideBuilder (extended conformation).

TODO: Integrate ESMFold for predicted folded structures.
  ESMFold gives realistic folded conformations + pLDDT confidence scores.
  Options:
    - ESMFold API (free, no install): https://api.esmatlas.com/foldSequence/v1/pdb/
    - Local ESMFold: pip install fair-esm (requires openfold + GPU torch)
  Both output the same PDB format — drop-in replacement for PeptideBuilder.
  When available, ESMFold structures will improve docking accuracy since
  HADDOCK3 gets a realistic starting conformation instead of an extended chain.

Usage:
    python3 modules/02_design/fold_peptides.py

Inputs:
    data/candidates/candidate_pool.json

Outputs:
    data/candidates/folded_structures/{candidate_id}.pdb
"""

import json
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]


def load_config(config_path: str = None) -> dict:
    import yaml
    path = Path(config_path) if config_path else ROOT / "config.yaml"
    with open(path) as f:
        return yaml.safe_load(f)


def build_structure(sequence: str) -> str:
    """
    Build a 3D peptide structure from sequence using PeptideBuilder.
    Returns the structure as a PDB-format string.
    """
    from PeptideBuilder import make_structure
    from Bio.PDB import PDBIO
    import io

    n = len(sequence)
    phi = [-180.0] * (n - 1)
    psi_im1 = [180.0] * (n - 1)

    structure = make_structure(sequence, phi, psi_im1)

    string_io = io.StringIO()
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(string_io)

    return string_io.getvalue()


def run(config_path: str = None) -> list[dict]:
    """
    Generate 3D structures for all candidates in the pool.
    """
    config = load_config(config_path)
    candidates_dir = Path(config["outputs"]["candidates"])
    pool_path = candidates_dir / "candidate_pool.json"

    if not pool_path.exists():
        log.error("candidate_pool.json not found. Run candidate_pool.py first.")
        return []

    with open(pool_path) as f:
        candidates = json.load(f)

    folded_dir = candidates_dir / "folded_structures"
    folded_dir.mkdir(parents=True, exist_ok=True)

    total = len(candidates)
    built = 0
    failed = 0

    log.info(f"Generating 3D structures for {total} peptides (PeptideBuilder)...")

    for i, candidate in enumerate(candidates):
        cid = candidate["id"]
        seq = candidate["sequence"]

        pdb_path = folded_dir / f"{cid}.pdb"

        # Skip if already built
        if pdb_path.exists() and candidate.get("folded_pdb") is not None:
            built += 1
            continue

        try:
            pdb_string = build_structure(seq)
            pdb_path.write_text(pdb_string)
            candidate["folded_pdb"] = str(pdb_path)
            built += 1

            if (i + 1) % 20 == 0:
                log.info(f"  [{i+1}/{total}] {built} structures built...")

        except Exception as e:
            candidate["folded_pdb"] = None
            failed += 1
            log.error(f"  [{i+1}/{total}] {cid} failed: {e}")

    # Save updated pool
    with open(pool_path, "w") as f:
        json.dump(candidates, f, indent=2)

    log.info(f"\nDone. {built} structures built, {failed} failed.")
    log.info(f"Structures saved to: {folded_dir}")

    return candidates


if __name__ == "__main__":
    run()
