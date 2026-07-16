"""
Module 05 — Step 2: Protease Stability Prediction
===================================================
Predicts protease cleavage sites in each candidate peptide.
Identifies where trypsin, chymotrypsin, and elastase would cut,
and assigns a risk level.

This is rule-based (regex on P1 position) — does not consider
3D structure or kinetics. Use as a flag for review, not definitive.

Usage:
    python3 modules/05_properties/protease_stability.py

Inputs:
    data/candidates/candidate_pool.json

Outputs:
    Updates candidate_pool.json with protease cleavage data
"""

import json
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]

# Protease cleavage rules (P1 position, not followed by P)
PROTEASES = {
    "trypsin": {
        "cuts_after": set("KR"),
        "blocked_by_next": "P",
        "relevance": "serum, intracellular",
    },
    "chymotrypsin": {
        "cuts_after": set("FYW"),
        "blocked_by_next": "P",
        "relevance": "serum, GI tract",
    },
    "elastase": {
        "cuts_after": set("AVSGT"),
        "blocked_by_next": "P",
        "relevance": "skin, neutrophils",
    },
}


def load_config(config_path: str = None) -> dict:
    import yaml
    path = Path(config_path) if config_path else ROOT / "config.yaml"
    with open(path) as f:
        return yaml.safe_load(f)


def find_cleavage_sites(sequence: str) -> dict:
    """
    Find protease cleavage sites in a peptide sequence.

    Returns dict with per-protease sites and overall risk.
    """
    results = {}
    total_sites = 0

    for protease_name, rules in PROTEASES.items():
        sites = []
        for i, aa in enumerate(sequence):
            if aa in rules["cuts_after"]:
                # Check if next residue blocks cleavage
                if i + 1 < len(sequence) and sequence[i + 1] == rules["blocked_by_next"]:
                    continue
                sites.append({
                    "position": i + 1,  # 1-indexed
                    "residue": aa,
                    "context": sequence[max(0, i-2):i+3],
                })

        results[protease_name] = {
            "count": len(sites),
            "sites": sites,
            "relevance": rules["relevance"],
        }
        total_sites += len(sites)

    # Risk assessment
    if total_sites == 0:
        risk = "low"
    elif total_sites <= 3:
        risk = "moderate"
    else:
        risk = "high"

    return {
        "total_cleavage_sites": total_sites,
        "protease_risk": risk,
        "proteases": results,
    }


def suggest_fixes(sequence: str, cleavage_data: dict) -> list[str]:
    """
    Suggest D-amino acid substitutions to improve protease resistance.
    """
    suggestions = []
    for protease, data in cleavage_data["proteases"].items():
        for site in data["sites"]:
            pos = site["position"]
            aa = site["residue"]
            suggestions.append(
                f"Position {pos} ({aa}): D-{aa} substitution blocks {protease} cleavage"
            )
    return suggestions


def run(config_path: str = None) -> list[dict]:
    """Compute protease stability for all candidates."""
    config = load_config(config_path)
    candidates_dir = Path(config["outputs"]["candidates"])
    pool_path = candidates_dir / "candidate_pool.json"

    if not pool_path.exists():
        log.error("candidate_pool.json not found.")
        return []

    with open(pool_path) as f:
        candidates = json.load(f)

    total = len(candidates)
    log.info(f"Analyzing protease stability for {total} candidates...")

    risk_counts = {"low": 0, "moderate": 0, "high": 0}

    for candidate in candidates:
        seq = candidate["sequence"]

        cleavage = find_cleavage_sites(seq)
        suggestions = suggest_fixes(seq, cleavage)

        candidate["protease_cleavage"] = cleavage
        candidate["protease_suggestions"] = suggestions

        risk_counts[cleavage["protease_risk"]] += 1

    # Save
    with open(pool_path, "w") as f:
        json.dump(candidates, f, indent=2)

    log.info(f"\nDone. Protease analysis for {total} candidates.")
    log.info(f"  Low risk:      {risk_counts['low']}")
    log.info(f"  Moderate risk: {risk_counts['moderate']}")
    log.info(f"  High risk:     {risk_counts['high']}")

    # Show details for first 10
    log.info(f"\n{'ID':<16} {'Sequence':<22} {'Sites':<7} {'Risk':<10} {'Trypsin':<8} {'Chymo':<8} {'Elastase'}")
    log.info(f"{'-'*16} {'-'*22} {'-'*7} {'-'*10} {'-'*8} {'-'*8} {'-'*8}")

    for c in candidates[:10]:
        cl = c.get("protease_cleavage", {})
        seq = c["sequence"][:20] + ".." if len(c["sequence"]) > 20 else c["sequence"]
        total_s = cl.get("total_cleavage_sites", 0)
        risk = cl.get("protease_risk", "?")
        tryp = cl.get("proteases", {}).get("trypsin", {}).get("count", 0)
        chymo = cl.get("proteases", {}).get("chymotrypsin", {}).get("count", 0)
        elast = cl.get("proteases", {}).get("elastase", {}).get("count", 0)
        log.info(f"{c['id']:<16} {seq:<22} {total_s:<7} {risk:<10} {tryp:<8} {chymo:<8} {elast}")

    return candidates


if __name__ == "__main__":
    run()
