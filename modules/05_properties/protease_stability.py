"""
Module 05 — Protease Stability Prediction
===========================================
Predicts protease cleavage sites. All reads/writes through run directory.

Usage:
    python3 modules/05_properties/protease_stability.py --run-dir path
"""

import argparse
import json
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]

PROTEASES = {
    "trypsin": {"cuts_after": set("KR"), "blocked_by_next": "P", "relevance": "serum, intracellular"},
    "chymotrypsin": {"cuts_after": set("FYW"), "blocked_by_next": "P", "relevance": "serum, GI tract"},
    "elastase": {"cuts_after": set("AVSGT"), "blocked_by_next": "P", "relevance": "skin, neutrophils"},
}


def find_cleavage_sites(sequence):
    results = {}
    total_sites = 0
    for name, rules in PROTEASES.items():
        sites = []
        for i, aa in enumerate(sequence):
            if aa in rules["cuts_after"]:
                if i + 1 < len(sequence) and sequence[i + 1] == rules["blocked_by_next"]:
                    continue
                sites.append({"position": i + 1, "residue": aa, "context": sequence[max(0,i-2):i+3]})
        results[name] = {"count": len(sites), "sites": sites, "relevance": rules["relevance"]}
        total_sites += len(sites)

    risk = "low" if total_sites == 0 else ("moderate" if total_sites <= 3 else "high")
    return {"total_cleavage_sites": total_sites, "protease_risk": risk, "proteases": results}


def suggest_fixes(sequence, cleavage_data):
    suggestions = []
    for protease, data in cleavage_data["proteases"].items():
        for site in data["sites"]:
            suggestions.append(f"Position {site['position']} ({site['residue']}): D-{site['residue']} blocks {protease}")
    return suggestions


def run(run_dir: str) -> list[dict]:
    sys.path.insert(0, str(ROOT))
    from modules.run_manager import RunManager

    rm = RunManager(run_dir=run_dir)
    candidates = rm.load_candidates()

    if not candidates:
        log.error("No candidates found.")
        return []

    log.info(f"Analyzing protease stability for {len(candidates)} candidates...")
    risk_counts = {"low": 0, "moderate": 0, "high": 0}

    for c in candidates:
        cleavage = find_cleavage_sites(c["sequence"])
        c["protease_cleavage"] = cleavage
        c["protease_suggestions"] = suggest_fixes(c["sequence"], cleavage)
        risk_counts[cleavage["protease_risk"]] += 1

    rm.save_candidates(candidates)
    log.info(f"Done. Low: {risk_counts['low']}, Moderate: {risk_counts['moderate']}, High: {risk_counts['high']}")
    return candidates


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    args = parser.parse_args()
    run(run_dir=args.run_dir)
