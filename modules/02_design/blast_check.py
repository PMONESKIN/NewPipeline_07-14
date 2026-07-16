"""
Module 02 — Step 4: BLAST Sequence Check
=========================================
Checks designed peptide sequences against UniProt/NCBI databases
to identify matches with known proteins or peptides.

Catches:
  - "You just redesigned a known protein fragment"
  - High similarity to known bioactive peptides
  - Potential toxin/allergen matches

Uses NCBI BLAST web API via Biopython (free, no key needed).

Usage:
    python3 modules/02_design/blast_check.py

Inputs:
    data/candidates/candidate_pool.json

Outputs:
    Updates candidate_pool.json with blast_hits field
    data/candidates/blast_results.json (detailed results)
"""

import json
import logging
import time
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]


def load_config(config_path: str = None) -> dict:
    import yaml
    path = Path(config_path) if config_path else ROOT / "config.yaml"
    with open(path) as f:
        return yaml.safe_load(f)


def blast_sequence(sequence: str, max_hits: int = 5) -> list[dict]:
    """
    Run BLASTP against NCBI nr database for a peptide sequence.

    Args:
        sequence: amino acid sequence
        max_hits: maximum number of hits to return

    Returns:
        List of hit dicts with: title, organism, identity, evalue, coverage
    """
    from Bio.Blast import NCBIWWW, NCBIXML

    try:
        # Submit BLAST query
        result_handle = NCBIWWW.qblast(
            "blastp",           # protein BLAST
            "swissprot",        # SwissProt database (curated, faster than nr)
            sequence,
            hitlist_size=max_hits,
            expect=10.0,        # e-value threshold (relaxed for short peptides)
            word_size=2,        # short word size for short sequences
            matrix_name="PAM30", # better for short sequences than BLOSUM62
        )

        blast_records = NCBIXML.parse(result_handle)
        blast_record = next(blast_records)

        hits = []
        for alignment in blast_record.alignments[:max_hits]:
            for hsp in alignment.hsps[:1]:  # top HSP per hit
                identity = hsp.identities / hsp.align_length * 100 if hsp.align_length > 0 else 0
                coverage = hsp.align_length / len(sequence) * 100

                hits.append({
                    "title": alignment.title[:100],
                    "identity_pct": round(identity, 1),
                    "evalue": hsp.expect,
                    "coverage_pct": round(coverage, 1),
                    "aligned_query": str(hsp.query),
                    "aligned_subject": str(hsp.sbjct),
                    "score": hsp.score,
                })

        return hits

    except Exception as e:
        log.error(f"  BLAST failed: {e}")
        return []


def classify_hits(hits: list[dict], sequence: str) -> dict:
    """
    Classify BLAST results into categories.

    Returns:
        Dict with: known_match (bool), best_identity, flags, summary
    """
    if not hits:
        return {
            "known_match": False,
            "best_identity": 0,
            "flags": [],
            "summary": "No significant matches found — novel sequence",
        }

    best = hits[0]
    best_identity = best["identity_pct"]

    flags = []
    if best_identity >= 95:
        flags.append("KNOWN_SEQUENCE")
    elif best_identity >= 80:
        flags.append("HIGH_SIMILARITY")
    elif best_identity >= 60:
        flags.append("MODERATE_SIMILARITY")

    # Check for concerning matches
    title_lower = best["title"].lower()
    if any(word in title_lower for word in ["toxin", "venom", "poison"]):
        flags.append("TOXIN_MATCH")
    if any(word in title_lower for word in ["allergen", "allergenic"]):
        flags.append("ALLERGEN_MATCH")

    summary = f"Best match: {best_identity:.0f}% identity — {best['title'][:60]}"

    return {
        "known_match": best_identity >= 80,
        "best_identity": best_identity,
        "flags": flags,
        "summary": summary,
    }


def run(config_path: str = None) -> list[dict]:
    """
    Run BLAST check on all candidates.

    Note: NCBI BLAST has rate limits. This will take ~30 seconds per
    candidate due to server-side processing. For 121 candidates,
    expect ~1 hour total.
    """
    config = load_config(config_path)
    candidates_dir = Path(config["outputs"]["candidates"])
    pool_path = candidates_dir / "candidate_pool.json"

    if not pool_path.exists():
        log.error("candidate_pool.json not found. Run candidate_pool.py first.")
        return []

    with open(pool_path) as f:
        candidates = json.load(f)

    total = len(candidates)
    checked = 0
    known_matches = 0

    log.info(f"Running BLAST check on {total} candidates...")
    log.info(f"(~30 seconds per candidate, total ~{total * 30 // 60} minutes)")

    all_blast_results = {}

    for i, candidate in enumerate(candidates):
        cid = candidate["id"]
        seq = candidate["sequence"]

        # Skip if already checked
        if candidate.get("blast_checked"):
            checked += 1
            continue

        log.info(f"[{i+1}/{total}] {cid} ({seq})")

        hits = blast_sequence(seq)
        classification = classify_hits(hits, seq)

        candidate["blast_checked"] = True
        candidate["blast_best_identity"] = classification["best_identity"]
        candidate["blast_known_match"] = classification["known_match"]
        candidate["blast_flags"] = classification["flags"]
        candidate["blast_summary"] = classification["summary"]

        all_blast_results[cid] = {
            "sequence": seq,
            "hits": hits,
            "classification": classification,
        }

        checked += 1
        if classification["known_match"]:
            known_matches += 1
            log.info(f"  KNOWN MATCH: {classification['summary']}")
        elif classification["flags"]:
            log.info(f"  {', '.join(classification['flags'])}: {classification['summary']}")
        else:
            log.info(f"  Novel — no significant matches")

        # Rate limit: NCBI asks for max 3 requests per second
        time.sleep(3)

        # Checkpoint every 10 candidates
        if (i + 1) % 10 == 0:
            with open(pool_path, "w") as f:
                json.dump(candidates, f, indent=2)
            log.info(f"  Checkpoint saved ({checked} checked)")

    # Final save
    with open(pool_path, "w") as f:
        json.dump(candidates, f, indent=2)

    # Save detailed results
    blast_path = candidates_dir / "blast_results.json"
    with open(blast_path, "w") as f:
        json.dump(all_blast_results, f, indent=2)

    log.info(f"\nBLAST check complete.")
    log.info(f"  Checked: {checked}/{total}")
    log.info(f"  Known matches: {known_matches}")
    log.info(f"  Novel sequences: {checked - known_matches}")

    return candidates


if __name__ == "__main__":
    run()
