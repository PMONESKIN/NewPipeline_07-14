"""
Module 04 — Step 2: Re-dock Modified Peptides
===============================================
Runs HADDOCK3 docking on modified peptide candidates to verify
that modifications (CPP tags, D-amino acids, etc.) don't disrupt binding.

Usage:
    python3 modules/04_modification/redock_modified.py --run-dir path [--mode fast]

Inputs:
    {run_dir}/candidates/modified_candidates.json
    data/candidates/folded_structures/{modified_id}.pdb

Outputs:
    {run_dir}/docking/haddock3_modified/ (per-candidate results)
    Updates modified_candidates.json with HADDOCK scores
"""

import argparse
import importlib.util
import json
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]


def load_haddock_module():
    """Load haddock3_dock.py using importlib (numbered directory)."""
    module_path = ROOT / "modules" / "03_docking" / "haddock3_dock.py"
    spec = importlib.util.spec_from_file_location("haddock3_dock", str(module_path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run(run_dir: str, mode: str = "fast") -> list[dict]:
    """Re-dock all modified candidates through HADDOCK3."""
    sys.path.insert(0, str(ROOT))
    from modules.run_manager import RunManager

    rm = RunManager(run_dir=run_dir)
    config = rm.config

    mod_path = rm.candidates_dir / "modified_candidates.json"
    if not mod_path.exists():
        log.error("No modified candidates found. Run modify_peptides.py first.")
        return []

    with open(mod_path) as f:
        modified = json.load(f)

    if not modified:
        log.error("Modified candidates file is empty.")
        return []

    log.info(f"Re-docking {len(modified)} modified candidate(s) in {mode} mode...")

    # Load HADDOCK3 docking functions
    h3 = load_haddock_module()

    # Prepare receptor
    target_cfg = config["targets"][0]
    interface = rm.load_interface(target_cfg["name"])
    pdb_path = rm.structures_dir / f"{interface['pdb_id']}.pdb"

    receptor_dir = rm.docking_dir / "receptor"
    receptor_dir.mkdir(parents=True, exist_ok=True)
    receptor_pdb = receptor_dir / f"{interface['pdb_id']}_receptor.pdb"

    if not receptor_pdb.exists():
        log.info(f"Extracting receptor chain {interface['receptor_chain']}...")
        h3.extract_receptor_pdb(pdb_path, interface["receptor_chain"], receptor_pdb)

    folded_dir = rm.folded_dir
    haddock_dir = rm.docking_dir / "haddock3_modified"
    haddock_dir.mkdir(parents=True, exist_ok=True)

    haddock_cfg = config.get("docking", {}).get("haddock3", {})
    sampling = haddock_cfg.get("sampling", {})
    score_field = "haddock_fast_score" if mode == "fast" else "haddock_full_score"

    scored = 0
    failed = 0

    for i, candidate in enumerate(modified):
        mid = candidate["id"]
        seq = candidate["sequence"]

        if candidate.get(score_field) is not None:
            scored += 1
            continue

        log.info(f"\n[{i+1}/{len(modified)}] {mid} ({seq})")

        score = h3.dock_candidate(
            candidate=candidate,
            receptor_pdb=receptor_pdb,
            interface=interface,
            folded_dir=folded_dir,
            output_dir=haddock_dir,
            mode=mode,
            sampling=sampling,
        )

        if score is not None:
            candidate[score_field] = score
            scored += 1
            log.info(f"  HADDOCK score: {score:.1f}")

            # Compare to original
            orig_id = candidate.get("original_id")
            if orig_id:
                log.info(f"  Original: {candidate['original_sequence']}")
                log.info(f"  Modified: {seq}")
                log.info(f"  Modifications: {', '.join(candidate.get('modifications', []))}")
        else:
            failed += 1
            log.error(f"  Re-docking failed")

    # Save updated modified candidates
    with open(mod_path, "w") as f:
        json.dump(modified, f, indent=2)

    log.info(f"\nRe-docking complete.")
    log.info(f"  Scored: {scored}/{len(modified)}")
    log.info(f"  Failed: {failed}/{len(modified)}")

    # Print comparison table
    if scored > 0:
        log.info(f"\n{'ID':<20} {'Original':<20} {'Modified':<25} {'Score':<10} {'Mods'}")
        log.info(f"{'-'*20} {'-'*20} {'-'*25} {'-'*10} {'-'*30}")
        for c in modified:
            if c.get(score_field) is not None:
                orig = c.get("original_sequence", "")[:18]
                mod_seq = c["sequence"][:23]
                score = f"{c[score_field]:.1f}"
                mods = ", ".join(c.get("modifications", []))[:28]
                log.info(f"{c['id']:<20} {orig:<20} {mod_seq:<25} {score:<10} {mods}")

    rm.update_status("modified_redocked")

    return modified


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True, help="Run directory")
    parser.add_argument("--mode", choices=["fast", "full"], default="fast", help="Docking mode")
    args = parser.parse_args()
    run(run_dir=args.run_dir, mode=args.mode)
