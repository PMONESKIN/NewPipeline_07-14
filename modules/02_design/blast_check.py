"""
Module 02 — Step 4: BLAST Sequence Check
=========================================
Checks designed peptides against UniProt/NCBI for known matches.

Usage:
    python3 modules/02_design/blast_check.py --run-dir path
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]


def blast_sequence(sequence, max_hits=5):
    from Bio.Blast import NCBIWWW, NCBIXML
    try:
        result_handle = NCBIWWW.qblast("blastp", "swissprot", sequence,
                                        hitlist_size=max_hits, expect=10.0,
                                        word_size=2, matrix_name="PAM30")
        blast_record = next(NCBIXML.parse(result_handle))

        hits = []
        for alignment in blast_record.alignments[:max_hits]:
            for hsp in alignment.hsps[:1]:
                identity = hsp.identities / hsp.align_length * 100 if hsp.align_length > 0 else 0
                hits.append({
                    "title": alignment.title[:100],
                    "identity_pct": round(identity, 1),
                    "evalue": hsp.expect,
                    "coverage_pct": round(hsp.align_length / len(sequence) * 100, 1),
                })
        return hits
    except Exception as e:
        log.error(f"  BLAST failed: {e}")
        return []


def classify_hits(hits):
    if not hits:
        return {"known_match": False, "best_identity": 0, "flags": [],
                "summary": "No significant matches — novel sequence"}

    best = hits[0]
    flags = []
    if best["identity_pct"] >= 95: flags.append("KNOWN_SEQUENCE")
    elif best["identity_pct"] >= 80: flags.append("HIGH_SIMILARITY")
    elif best["identity_pct"] >= 60: flags.append("MODERATE_SIMILARITY")

    title_lower = best["title"].lower()
    if any(w in title_lower for w in ["toxin", "venom"]): flags.append("TOXIN_MATCH")
    if any(w in title_lower for w in ["allergen"]): flags.append("ALLERGEN_MATCH")

    return {
        "known_match": best["identity_pct"] >= 80,
        "best_identity": best["identity_pct"],
        "flags": flags,
        "summary": f"Best: {best['identity_pct']:.0f}% — {best['title'][:60]}",
    }


def run(run_dir: str) -> list[dict]:
    sys.path.insert(0, str(ROOT))
    from modules.run_manager import RunManager

    rm = RunManager(run_dir=run_dir)
    candidates = rm.load_candidates()

    if not candidates:
        log.error("No candidates found.")
        return []

    total = len(candidates)
    log.info(f"BLAST check on {total} candidates (~30 sec each, ~{total * 30 // 60} min total)...")

    for i, c in enumerate(candidates):
        if c.get("blast_checked"):
            continue

        log.info(f"[{i+1}/{total}] {c['id']} ({c['sequence']})")
        hits = blast_sequence(c["sequence"])
        classification = classify_hits(hits)

        c["blast_checked"] = True
        c["blast_best_identity"] = classification["best_identity"]
        c["blast_known_match"] = classification["known_match"]
        c["blast_summary"] = classification["summary"]

        if classification["known_match"]:
            log.info(f"  KNOWN: {classification['summary']}")
        else:
            log.info(f"  Novel")

        time.sleep(3)

        if (i + 1) % 10 == 0:
            rm.save_candidates(candidates)

    rm.save_candidates(candidates)
    known = sum(1 for c in candidates if c.get("blast_known_match"))
    log.info(f"\nDone. Known: {known}, Novel: {total - known}")
    return candidates


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    args = parser.parse_args()
    run(run_dir=args.run_dir)
