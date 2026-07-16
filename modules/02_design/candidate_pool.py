"""
Module 02 — Step 3: Candidate Pool Assembly
============================================
Merges seed sequences + ProteinMPNN designed sequences into a single
standardized candidate pool for Module 03 (docking).

Each candidate gets:
  - A unique ID (e.g., NRF2KEAP1-001, MYTARGET-001)
  - Standardized fields for docking + property prediction
  - Source tag (seed vs. designed)

Applies pre-filters before docking:
  - Length within config range
  - Standard amino acids only
  - No duplicate sequences

Usage:
    python3 modules/02_design/candidate_pool.py

Outputs:
    data/candidates/candidate_pool.json
    outputs/reports/candidate_pool_summary.md
"""

import json
import logging
import re
from collections import Counter
from datetime import date
from pathlib import Path

import yaml

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
STANDARD_AAS = set("ACDEFGHIKLMNPQRSTVWY")


def load_config(config_path: str = None) -> dict:
    path = Path(config_path) if config_path else ROOT / "config.yaml"
    with open(path) as f:
        return yaml.safe_load(f)


def make_prefix(target_name: str) -> str:
    """Generate a short ID prefix from the target name."""
    # Remove special chars, take first 4 chars uppercase
    clean = re.sub(r"[^a-zA-Z0-9]", "", target_name).upper()
    return clean[:8] if clean else "UNK"


def is_valid_sequence(seq: str, min_len: int, max_len: int) -> tuple[bool, str]:
    """Validate a peptide sequence."""
    if not seq:
        return False, "empty sequence"
    seq = seq.upper().strip()
    if not min_len <= len(seq) <= max_len:
        return False, f"length {len(seq)} outside [{min_len}, {max_len}]"
    non_standard = set(seq) - STANDARD_AAS
    if non_standard:
        return False, f"non-standard residues: {non_standard}"
    return True, ""


def load_seed_sequences(config: dict) -> list[dict]:
    """Load seed_sequences from each target's config entry."""
    seeds = []
    for target_cfg in config.get("targets", []):
        target_name = target_cfg.get("name", "Unknown")
        for seq in target_cfg.get("seed_sequences", []):
            seeds.append({
                "sequence": seq.upper().strip(),
                "length": len(seq),
                "target_name": target_name,
                "design_source": "seed",
                "mpnn_score": None,
                "sampling_temperature": None,
                "sequence_recovery": None,
            })
    if seeds:
        log.info(f"Loaded {len(seeds)} seed sequences from config.yaml")
    return seeds


def load_candidates(candidates_dir: Path, config: dict) -> list[dict]:
    """Load seed sequences + MPNN designs."""
    all_candidates = []

    seeds = load_seed_sequences(config)
    all_candidates.extend(seeds)

    mpnn_path = candidates_dir / "mpnn_designs.json"
    if mpnn_path.exists():
        with open(mpnn_path) as f:
            mpnn = json.load(f)
        log.info(f"Loaded {len(mpnn)} ProteinMPNN designed candidates")
        all_candidates.extend(mpnn)
    else:
        log.warning("mpnn_designs.json not found — run design_peptides.py first")

    return all_candidates


def build_pool(config_path: str = None) -> list[dict]:
    """
    Load all candidates, filter, deduplicate, assign IDs, and save.
    """
    config = load_config(config_path)
    candidates_dir = Path(config["outputs"]["candidates"])
    reports_dir = Path(config["outputs"]["reports"])
    candidates_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    min_len = config["candidates"]["length"]["min"]
    max_len = config["candidates"]["length"]["max"]

    raw_candidates = load_candidates(candidates_dir, config)

    if not raw_candidates:
        log.error("No candidates found. Run design_peptides.py first.")
        return []

    # Filter and deduplicate
    seen_sequences = set()
    filtered = []
    rejected_counts = {"length": 0, "non_standard": 0, "duplicate": 0}

    for c in raw_candidates:
        seq = c.get("sequence", "").upper().strip()
        valid, reason = is_valid_sequence(seq, min_len, max_len)

        if not valid:
            if "length" in reason:
                rejected_counts["length"] += 1
            elif "non-standard" in reason:
                rejected_counts["non_standard"] += 1
            continue

        if seq in seen_sequences:
            rejected_counts["duplicate"] += 1
            continue

        seen_sequences.add(seq)
        filtered.append(c)

    log.info(f"\nFiltering results:")
    log.info(f"  Input:     {len(raw_candidates)} candidates")
    log.info(f"  Rejected:  {sum(rejected_counts.values())} ({rejected_counts})")
    log.info(f"  Remaining: {len(filtered)} candidates")

    # Assign IDs per target
    pool = []
    target_counters = {}

    for c in filtered:
        target = c.get("target_name", "Unknown")
        target_counters[target] = target_counters.get(target, 0) + 1
        prefix = make_prefix(target)
        uid = f"{prefix}-{target_counters[target]:03d}"

        pool.append({
            "id": uid,
            "sequence": c["sequence"].upper(),
            "length": len(c["sequence"]),
            "target_name": target,
            "design_source": c.get("design_source", "unknown"),
            "mpnn_score": c.get("mpnn_score"),
            "sampling_temperature": c.get("sampling_temperature"),
            "sequence_recovery": c.get("sequence_recovery"),
            # Fields filled by Module 03 (docking)
            "vina_score": None,
            "haddock_score": None,
            "ensemble_score": None,
            "uncertainty_flag": None,
            # Fields filled by Module 04 (properties)
            "predicted_stability": None,
            "predicted_solubility": None,
            "predicted_permeability": None,
            "property_flags": [],
        })

    # Save pool
    pool_path = candidates_dir / "candidate_pool.json"
    with open(pool_path, "w") as f:
        json.dump(pool, f, indent=2)

    # Generate summary
    summary = generate_summary(pool, rejected_counts, len(raw_candidates))
    summary_path = reports_dir / "candidate_pool_summary.md"
    summary_path.write_text(summary)

    log.info(f"\nCandidate pool saved to {pool_path}")
    log.info(f"Summary report saved to {summary_path}")
    log.info(f"\nPool breakdown:")
    for target, count in target_counters.items():
        des = sum(1 for c in pool if c["target_name"] == target and c["design_source"] == "proteinmpnn")
        seed = sum(1 for c in pool if c["target_name"] == target and c["design_source"] == "seed")
        log.info(f"  {target}: {count} total ({seed} seeds, {des} designed)")

    return pool


def generate_summary(pool: list[dict], rejected: dict, total_input: int) -> str:
    target_counts = Counter(c["target_name"] for c in pool)
    source_counts = Counter(c["design_source"] for c in pool)
    length_dist = Counter(c["length"] for c in pool)

    lines = [
        f"# PeptideScreen — Candidate Pool Summary",
        f"**Date:** {date.today().isoformat()}  ",
        f"**Total candidates:** {len(pool)}",
        "",
        "---",
        "",
        "## Filtering",
        f"- Raw candidates: {total_input}",
        f"- Rejected (length): {rejected.get('length', 0)}",
        f"- Rejected (non-standard): {rejected.get('non_standard', 0)}",
        f"- Rejected (duplicates): {rejected.get('duplicate', 0)}",
        f"- **Final pool: {len(pool)}**",
        "",
        "## By Target",
        "",
        "| Target | Count |",
        "|--------|-------|",
    ]
    for target, count in sorted(target_counts.items()):
        lines.append(f"| {target} | {count} |")

    lines += [
        "",
        "## By Source",
        "",
        "| Source | Count |",
        "|--------|-------|",
    ]
    for source, count in sorted(source_counts.items()):
        lines.append(f"| {source} | {count} |")

    lines += [
        "",
        "## Length Distribution",
        "",
        "| Length (aa) | Count |",
        "|------------|-------|",
    ]
    for length in sorted(length_dist.keys()):
        lines.append(f"| {length} | {length_dist[length]} |")

    lines += [
        "",
        "## Next Step",
        "Run Module 03 (docking) to score all candidates.",
        "",
        "*Generated by PeptideScreen Module 02*",
    ]

    return "\n".join(lines)


if __name__ == "__main__":
    build_pool()
