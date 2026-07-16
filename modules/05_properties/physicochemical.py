"""
Module 05 — Step 1: Physicochemical Properties
================================================
Calculates core physicochemical properties for each candidate
using Biopython's ProtParam. All deterministic, no external calls.

Properties computed:
  - Molecular weight (Da)
  - Net charge at pH 7.4
  - Isoelectric point (pI)
  - GRAVY index (hydrophobicity, Kyte-Doolittle)
  - Instability index (>40 = likely unstable)
  - Aliphatic index (thermostability indicator)

Usage:
    python3 modules/05_properties/physicochemical.py

Inputs:
    data/candidates/candidate_pool.json

Outputs:
    Updates candidate_pool.json with physicochemical properties
"""

import json
import logging
from pathlib import Path

from Bio.SeqUtils.ProtParam import ProteinAnalysis

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]

INSTABILITY_THRESHOLD = 40.0


def load_config(config_path: str = None) -> dict:
    import yaml
    path = Path(config_path) if config_path else ROOT / "config.yaml"
    with open(path) as f:
        return yaml.safe_load(f)


def compute_properties(sequence: str) -> dict:
    """
    Compute all physicochemical properties for a peptide sequence.

    Returns dict with all computed values.
    """
    analysis = ProteinAnalysis(sequence)

    mw = analysis.molecular_weight()
    charge = analysis.charge_at_pH(7.4)
    pi = analysis.isoelectric_point()
    gravy = analysis.gravy()
    instability = analysis.instability_index()

    # Aliphatic index (not in ProtParam, compute manually)
    # AI = X(Ala) + a*X(Val) + b*X(Ile+Leu)
    # where a=2.9, b=3.9, X = mole percent
    n = len(sequence)
    ala_pct = sequence.count("A") / n * 100
    val_pct = sequence.count("V") / n * 100
    ile_leu_pct = (sequence.count("I") + sequence.count("L")) / n * 100
    aliphatic = ala_pct + 2.9 * val_pct + 3.9 * ile_leu_pct

    return {
        "molecular_weight_da": round(mw, 1),
        "net_charge_ph74": round(charge, 2),
        "isoelectric_point": round(pi, 2),
        "gravy": round(gravy, 3),
        "instability_index": round(instability, 1),
        "is_unstable": instability > INSTABILITY_THRESHOLD,
        "aliphatic_index": round(aliphatic, 1),
    }


def run(config_path: str = None) -> list[dict]:
    """Compute physicochemical properties for all candidates."""
    config = load_config(config_path)
    candidates_dir = Path(config["outputs"]["candidates"])
    pool_path = candidates_dir / "candidate_pool.json"

    if not pool_path.exists():
        log.error("candidate_pool.json not found.")
        return []

    with open(pool_path) as f:
        candidates = json.load(f)

    total = len(candidates)
    log.info(f"Computing physicochemical properties for {total} candidates...")

    for i, candidate in enumerate(candidates):
        seq = candidate["sequence"]

        try:
            props = compute_properties(seq)
            candidate["physicochemical"] = props
        except Exception as e:
            log.error(f"  {candidate['id']}: failed — {e}")
            candidate["physicochemical"] = None

    # Save
    with open(pool_path, "w") as f:
        json.dump(candidates, f, indent=2)

    # Summary stats
    mws = [c["physicochemical"]["molecular_weight_da"] for c in candidates if c.get("physicochemical")]
    unstable = sum(1 for c in candidates if c.get("physicochemical", {}).get("is_unstable"))

    log.info(f"\nDone. Properties computed for {total} candidates.")
    log.info(f"  MW range: {min(mws):.0f} - {max(mws):.0f} Da")
    log.info(f"  Unstable (index > 40): {unstable}/{total}")

    # Print top-level summary table
    log.info(f"\n{'ID':<16} {'MW':<10} {'Charge':<8} {'pI':<6} {'GRAVY':<8} {'Instab':<8} {'Stable?'}")
    log.info(f"{'-'*16} {'-'*10} {'-'*8} {'-'*6} {'-'*8} {'-'*8} {'-'*7}")

    for c in candidates[:15]:
        p = c.get("physicochemical")
        if not p:
            continue
        stable = "yes" if not p["is_unstable"] else "NO"
        log.info(
            f"{c['id']:<16} {p['molecular_weight_da']:<10.0f} {p['net_charge_ph74']:<8.1f} "
            f"{p['isoelectric_point']:<6.1f} {p['gravy']:<8.3f} {p['instability_index']:<8.1f} {stable}"
        )

    if total > 15:
        log.info(f"  ... and {total - 15} more")

    return candidates


if __name__ == "__main__":
    run()
