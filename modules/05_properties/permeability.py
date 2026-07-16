"""
Module 05 — Step 3: Permeability Prediction
=============================================
Estimates membrane permeability for each candidate peptide.

Methods (in order of preference):
  1. PerMM web API — physics-based, uses 3D structure (when available)
  2. Rule-based — MW/TPSA heuristics (always available)

TODO: Integrate local PerMM (PyPerMM) when available.

Usage:
    python3 modules/05_properties/permeability.py

Inputs:
    data/candidates/candidate_pool.json

Outputs:
    Updates candidate_pool.json with permeability predictions
"""

import json
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]

# MW-based permeability categories (Bos & Meinardi "500 Da rule")
MW_CATEGORIES = [
    (500, "possible_passive", "May passively permeate membranes"),
    (1000, "requires_enhancer", "Needs permeation enhancer or delivery vehicle"),
    (2000, "unlikely_passive", "Unlikely to cross membranes passively, needs CPP or delivery system"),
    (float("inf"), "impermeable_passive", "Too large for passive permeation, requires active delivery"),
]


def load_config(config_path: str = None) -> dict:
    import yaml
    path = Path(config_path) if config_path else ROOT / "config.yaml"
    with open(path) as f:
        return yaml.safe_load(f)


def predict_by_rules(sequence: str, mw: float = None) -> dict:
    """
    Rule-based permeability prediction using MW and sequence properties.

    This is a rough screen — not a substitute for experimental PAMPA/Caco-2.
    """
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    analysis = ProteinAnalysis(sequence)
    if mw is None:
        mw = analysis.molecular_weight()

    # MW category
    category = "unknown"
    description = ""
    for threshold, cat, desc in MW_CATEGORIES:
        if mw < threshold:
            category = cat
            description = desc
            break

    # Charge — highly charged peptides have poor passive permeability
    charge = abs(analysis.charge_at_pH(7.4))
    charge_penalty = "high_charge" if charge > 3 else "acceptable_charge"

    # Hydrophobicity — more hydrophobic = better passive permeation (generally)
    gravy = analysis.gravy()
    if gravy > 0:
        hydrophobicity = "hydrophobic_favorable"
    elif gravy > -1:
        hydrophobicity = "moderate"
    else:
        hydrophobicity = "hydrophilic_unfavorable"

    # Hydrogen bond donors — more HBD = worse permeability (Lipinski)
    # Rough estimate: count N-H and O-H capable sidechains
    hbd_residues = set("STNQKRHYW")
    hbd_count = sum(1 for aa in sequence if aa in hbd_residues)
    hbd_penalty = "excessive_hbd" if hbd_count > 5 else "acceptable_hbd"

    # Overall assessment
    if category == "possible_passive" and charge_penalty == "acceptable_charge":
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
        "charge_assessment": charge_penalty,
        "gravy": round(gravy, 3),
        "hydrophobicity_assessment": hydrophobicity,
        "hbd_count": hbd_count,
        "hbd_assessment": hbd_penalty,
        "overall_permeability": overall,
        "confidence": "low",
        "note": "Rule-based prediction only. Experimental PAMPA/Caco-2 needed for validation.",
    }


def run(config_path: str = None) -> list[dict]:
    """Predict permeability for all candidates."""
    config = load_config(config_path)
    candidates_dir = Path(config["outputs"]["candidates"])
    pool_path = candidates_dir / "candidate_pool.json"

    if not pool_path.exists():
        log.error("candidate_pool.json not found.")
        return []

    with open(pool_path) as f:
        candidates = json.load(f)

    total = len(candidates)
    log.info(f"Predicting permeability for {total} candidates...")

    category_counts = {}

    for candidate in candidates:
        seq = candidate["sequence"]

        # Use physicochemical MW if already computed
        mw = None
        if candidate.get("physicochemical"):
            mw = candidate["physicochemical"].get("molecular_weight_da")

        perm = predict_by_rules(seq, mw)
        candidate["permeability"] = perm

        cat = perm["mw_category"]
        category_counts[cat] = category_counts.get(cat, 0) + 1

    # Save
    with open(pool_path, "w") as f:
        json.dump(candidates, f, indent=2)

    log.info(f"\nDone. Permeability predicted for {total} candidates.")
    log.info(f"\n  Category breakdown:")
    for cat, count in sorted(category_counts.items()):
        log.info(f"    {cat}: {count}")

    # Show details for first 10
    log.info(f"\n{'ID':<16} {'MW':<10} {'Category':<22} {'Charge':<8} {'GRAVY':<8} {'Overall'}")
    log.info(f"{'-'*16} {'-'*10} {'-'*22} {'-'*8} {'-'*8} {'-'*20}")

    for c in candidates[:10]:
        p = c.get("permeability", {})
        log.info(
            f"{c['id']:<16} {p.get('molecular_weight_da', 0):<10.0f} "
            f"{p.get('mw_category', '?'):<22} {p.get('net_charge_magnitude', 0):<8.1f} "
            f"{p.get('gravy', 0):<8.3f} {p.get('overall_permeability', '?')}"
        )

    log.info(f"\n  Note: These are rule-based estimates (low confidence).")
    log.info(f"  Experimental PAMPA or Caco-2 assays required for validation.")

    return candidates


if __name__ == "__main__":
    run()
