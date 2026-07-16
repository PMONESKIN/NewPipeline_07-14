"""
Module 02 — Step 3: Peptide 3D Structure Generation
====================================================
Generates 3D PDB structures for each candidate using PeptideBuilder.

TODO: Integrate ESMFold for predicted folded structures.
  ESMFold gives realistic folded conformations + pLDDT confidence scores.
  Both output the same PDB format — drop-in replacement for PeptideBuilder.

Usage:
    python3 modules/02_design/fold_peptides.py --run-dir path
"""

import argparse
import json
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]


def build_structure(sequence: str) -> str:
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


def run(run_dir: str) -> list[dict]:
    sys.path.insert(0, str(ROOT))
    from modules.run_manager import RunManager

    rm = RunManager(run_dir=run_dir)
    candidates = rm.load_candidates()

    if not candidates:
        log.error("No candidates found. Run candidate_pool.py first.")
        return []

    folded_dir = rm.folded_dir
    total = len(candidates)
    built = 0

    log.info(f"Generating 3D structures for {total} peptides...")

    for i, candidate in enumerate(candidates):
        cid = candidate["id"]
        seq = candidate["sequence"]
        pdb_path = folded_dir / f"{cid}.pdb"

        if pdb_path.exists() and candidate.get("folded_pdb"):
            built += 1
            continue

        try:
            pdb_string = build_structure(seq)
            pdb_path.write_text(pdb_string)
            candidate["folded_pdb"] = str(pdb_path)
            built += 1
        except Exception as e:
            candidate["folded_pdb"] = None
            log.error(f"  {cid} failed: {e}")

        if (i + 1) % 20 == 0:
            log.info(f"  [{i+1}/{total}] {built} structures built...")

    rm.save_candidates(candidates)
    log.info(f"\nDone. {built} structures built.")
    return candidates


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    args = parser.parse_args()
    run(run_dir=args.run_dir)
