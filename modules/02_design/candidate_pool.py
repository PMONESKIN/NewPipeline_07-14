"""
Module 02 — Step 2: Candidate Pool Assembly
============================================
Merges seed sequences + ProteinMPNN designs, deduplicates, assigns IDs.
All reads/writes go through the run directory.

Usage:
    python3 modules/02_design/candidate_pool.py --run-dir path
"""

import argparse
import json
import logging
import re
import sys
from collections import Counter
from datetime import date
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
STANDARD_AAS = set("ACDEFGHIKLMNPQRSTVWY")


def make_prefix(target_name):
    clean = re.sub(r"[^a-zA-Z0-9]", "", target_name).upper()
    return clean[:8] if clean else "UNK"


def is_valid_sequence(seq, min_len, max_len):
    if not seq:
        return False, "empty"
    seq = seq.upper().strip()
    if not min_len <= len(seq) <= max_len:
        return False, f"length {len(seq)}"
    non_standard = set(seq) - STANDARD_AAS
    if non_standard:
        return False, f"non-standard: {non_standard}"
    return True, ""


def run(run_dir: str) -> list[dict]:
    sys.path.insert(0, str(ROOT))
    from modules.run_manager import RunManager

    rm = RunManager(run_dir=run_dir)
    config = rm.config

    min_len = config["candidates"]["length"]["min"]
    max_len = config["candidates"]["length"]["max"]

    all_candidates = []

    # Load seeds from config
    for target_cfg in config.get("targets", []):
        for seq in target_cfg.get("seed_sequences", []):
            all_candidates.append({
                "sequence": seq.upper().strip(),
                "length": len(seq),
                "target_name": target_cfg["name"],
                "design_source": "seed",
                "mpnn_score": None, "sampling_temperature": None, "sequence_recovery": None,
            })

    # Load MPNN designs from run
    mpnn_path = rm.candidates_dir / "mpnn_designs.json"
    if mpnn_path.exists():
        with open(mpnn_path) as f:
            all_candidates.extend(json.load(f))
        log.info(f"Loaded MPNN designs from {mpnn_path}")

    if not all_candidates:
        log.error("No candidates found.")
        return []

    # Filter + dedup
    seen = set()
    filtered = []
    rejected = {"length": 0, "non_standard": 0, "duplicate": 0}

    for c in all_candidates:
        seq = c.get("sequence", "").upper().strip()
        valid, reason = is_valid_sequence(seq, min_len, max_len)
        if not valid:
            if "length" in reason: rejected["length"] += 1
            else: rejected["non_standard"] += 1
            continue
        if seq in seen:
            rejected["duplicate"] += 1
            continue
        seen.add(seq)
        filtered.append(c)

    # Assign IDs
    pool = []
    counters = {}
    for c in filtered:
        target = c.get("target_name", "Unknown")
        counters[target] = counters.get(target, 0) + 1
        prefix = make_prefix(target)
        uid = f"{prefix}-{counters[target]:03d}"

        pool.append({
            "id": uid,
            "sequence": c["sequence"].upper(),
            "length": len(c["sequence"]),
            "target_name": target,
            "design_source": c.get("design_source", "unknown"),
            "mpnn_score": c.get("mpnn_score"),
            "sampling_temperature": c.get("sampling_temperature"),
            "sequence_recovery": c.get("sequence_recovery"),
            "haddock_fast_score": None,
            "haddock_full_score": None,
        })

    rm.save_candidates(pool)

    # Summary report
    summary = f"# Candidate Pool Summary\n**Date:** {date.today()}\n**Total:** {len(pool)}\n"
    summary += f"\nRejected: {rejected}\n"
    rm.reports_dir.mkdir(parents=True, exist_ok=True)
    (rm.reports_dir / "candidate_pool_summary.md").write_text(summary)

    log.info(f"\nPool: {len(pool)} candidates (rejected {sum(rejected.values())})")
    for target, count in counters.items():
        log.info(f"  {target}: {count}")

    return pool


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    args = parser.parse_args()
    run(run_dir=args.run_dir)
