"""
Module 05 — Permeability Prediction
=====================================
Rule-based permeability estimates. All reads/writes through run directory.

Usage:
    python3 modules/05_properties/permeability.py --run-dir path
"""

import argparse
import json
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]

MW_CATEGORIES = [
    (500, "possible_passive", "May passively permeate"),
    (1000, "requires_enhancer", "Needs permeation enhancer"),
    (2000, "unlikely_passive", "Needs CPP or delivery system"),
    (float("inf"), "impermeable_passive", "Requires active delivery"),
]


def predict_by_rules(sequence, mw=None):
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    analysis = ProteinAnalysis(sequence)
    if mw is None:
        mw = analysis.molecular_weight()

    category = description = ""
    for threshold, cat, desc in MW_CATEGORIES:
        if mw < threshold:
            category, description = cat, desc
            break

    charge = abs(analysis.charge_at_pH(7.4))
    gravy = analysis.gravy()
    hbd_count = sum(1 for aa in sequence if aa in "STNQKRHYW")

    if category == "possible_passive" and charge <= 3:
        overall = "potentially_permeable"
    elif category == "requires_enhancer":
        overall = "needs_enhancement"
    else:
        overall = "needs_delivery_system"

    return {
        "method": "rule_based",
        "molecular_weight_da": round(mw, 1),
        "mw_category": category,
        "mw_description": description,
        "net_charge_magnitude": round(charge, 1),
        "gravy": round(gravy, 3),
        "hbd_count": hbd_count,
        "overall_permeability": overall,
        "confidence": "low",
        "note": "Rule-based only. Experimental PAMPA/Caco-2 needed.",
    }


def run(run_dir: str) -> list[dict]:
    sys.path.insert(0, str(ROOT))
    from modules.run_manager import RunManager

    rm = RunManager(run_dir=run_dir)
    candidates = rm.load_candidates()

    if not candidates:
        log.error("No candidates found.")
        return []

    log.info(f"Predicting permeability for {len(candidates)} candidates...")
    cats = {}

    for c in candidates:
        mw = c.get("physicochemical", {}).get("molecular_weight_da")
        c["permeability"] = predict_by_rules(c["sequence"], mw)
        cat = c["permeability"]["mw_category"]
        cats[cat] = cats.get(cat, 0) + 1

    rm.save_candidates(candidates)
    log.info(f"Done. Categories: {cats}")
    return candidates


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    args = parser.parse_args()
    run(run_dir=args.run_dir)
