"""
Module 03 — Step 2: Score Ranking
==================================
Ranks candidates by HADDOCK3 docking score.
Lower HADDOCK score = better predicted binding.

Usage:
    python3 modules/03_docking/score_ranking.py --run-dir path [--mode fast]

Inputs:
    {run_dir}/candidates/candidate_pool.json (with haddock scores)

Outputs:
    Updates candidate_pool.json with rank
    Prints ranked table to console
"""

import argparse
import json
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]


def run(run_dir: str, mode: str = "fast") -> list[dict]:
    """
    Rank candidates by HADDOCK3 score.
    """
    sys.path.insert(0, str(ROOT))
    from modules.run_manager import RunManager

    rm = RunManager(run_dir=run_dir)
    candidates = rm.load_candidates()

    score_field = "haddock_fast_score" if mode == "fast" else "haddock_full_score"
    rank_field = f"rank_{mode}"

    scored = [c for c in candidates if c.get(score_field) is not None]
    unscored = [c for c in candidates if c.get(score_field) is None]

    if not scored:
        log.error(f"No candidates with {score_field}. Run haddock3_dock.py --mode {mode} first.")
        return candidates

    # Sort by score (lower = better)
    scored.sort(key=lambda c: c[score_field])

    # Assign ranks
    for rank, c in enumerate(scored, 1):
        for candidate in candidates:
            if candidate["id"] == c["id"]:
                candidate[rank_field] = rank
                break

    rm.save_candidates(candidates)

    # Print results
    log.info(f"\n{'='*70}")
    log.info(f"HADDOCK3 {mode.upper()} screen ranking — {len(scored)} candidates scored")
    log.info(f"{'='*70}")
    log.info(f"")
    log.info(f"{'Rank':<6} {'ID':<16} {'Sequence':<22} {'Score':<10} {'MPNN':<8}")
    log.info(f"{'-'*6} {'-'*16} {'-'*22} {'-'*10} {'-'*8}")

    for rank, c in enumerate(scored, 1):
        seq = c["sequence"][:20] + ".." if len(c["sequence"]) > 20 else c["sequence"]
        score = f"{c[score_field]:.1f}"
        mpnn = f"{c['mpnn_score']:.2f}" if c.get("mpnn_score") is not None else "N/A"
        marker = " ***" if rank <= 10 else ""
        log.info(f"{rank:<6} {c['id']:<16} {seq:<22} {score:<10} {mpnn:<8}{marker}")

        if rank == 10:
            log.info(f"{'-'*6} {'-'*16} {'-'*22} {'-'*10} {'-'*8}")

    log.info(f"\n  Total scored: {len(scored)}")
    log.info(f"  Unscored:    {len(unscored)}")
    log.info(f"  Best score:  {scored[0][score_field]:.1f} ({scored[0]['id']})")
    log.info(f"  Worst score: {scored[-1][score_field]:.1f} ({scored[-1]['id']})")

    if mode == "fast":
        log.info(f"\nNext step: review top candidates, then run full validation:")
        log.info(f"  python3 modules/03_docking/haddock3_dock.py --run-dir {run_dir} --mode full --n 20")

    rm.update_status(f"ranked_{mode}")

    return candidates


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True, help="Run directory")
    parser.add_argument("--mode", choices=["fast", "full"], default="fast", help="Which scores to rank")
    args = parser.parse_args()
    run(run_dir=args.run_dir, mode=args.mode)
