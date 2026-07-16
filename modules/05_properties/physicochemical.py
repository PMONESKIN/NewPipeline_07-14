"""
Module 05 — Physicochemical Properties
========================================
Calculates MW, charge, pI, GRAVY, instability index for each candidate.
All reads/writes go through the run directory.

Usage:
    python3 modules/05_properties/physicochemical.py --run-dir path
"""

import argparse
import json
import logging
import sys
from pathlib import Path

from Bio.SeqUtils.ProtParam import ProteinAnalysis

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
INSTABILITY_THRESHOLD = 40.0


def compute_properties(sequence):
    analysis = ProteinAnalysis(sequence)
    mw = analysis.molecular_weight()
    n = len(sequence)
    ala_pct = sequence.count("A") / n * 100
    val_pct = sequence.count("V") / n * 100
    ile_leu_pct = (sequence.count("I") + sequence.count("L")) / n * 100
    aliphatic = ala_pct + 2.9 * val_pct + 3.9 * ile_leu_pct
    instability = analysis.instability_index()

    return {
        "molecular_weight_da": round(mw, 1),
        "net_charge_ph74": round(analysis.charge_at_pH(7.4), 2),
        "isoelectric_point": round(analysis.isoelectric_point(), 2),
        "gravy": round(analysis.gravy(), 3),
        "instability_index": round(instability, 1),
        "is_unstable": instability > INSTABILITY_THRESHOLD,
        "aliphatic_index": round(aliphatic, 1),
    }


def run(run_dir: str) -> list[dict]:
    sys.path.insert(0, str(ROOT))
    from modules.run_manager import RunManager

    rm = RunManager(run_dir=run_dir)
    candidates = rm.load_candidates()

    if not candidates:
        log.error("No candidates found.")
        return []

    log.info(f"Computing physicochemical properties for {len(candidates)} candidates...")

    for c in candidates:
        try:
            c["physicochemical"] = compute_properties(c["sequence"])
        except Exception as e:
            log.error(f"  {c['id']}: {e}")
            c["physicochemical"] = None

    rm.save_candidates(candidates)

    mws = [c["physicochemical"]["molecular_weight_da"] for c in candidates if c.get("physicochemical")]
    unstable = sum(1 for c in candidates if c.get("physicochemical", {}).get("is_unstable"))
    log.info(f"Done. MW range: {min(mws):.0f}-{max(mws):.0f} Da. Unstable: {unstable}/{len(candidates)}")
    return candidates


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    args = parser.parse_args()
    run(run_dir=args.run_dir)
