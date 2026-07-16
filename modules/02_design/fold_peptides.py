"""
Module 02 — Step 3: Peptide 3D Structure Generation
====================================================
Generates 3D PDB structures for each candidate.

Tries ESMFold first (realistic predicted fold + pLDDT scores).
Falls back to PeptideBuilder (extended conformation) if ESMFold
is not available.

ESMFold requires: pip install fair-esm openfold torch (with GPU)

Usage:
    python3 modules/02_design/fold_peptides.py --run-dir path
    python3 modules/02_design/fold_peptides.py --run-dir path --method esmfold
    python3 modules/02_design/fold_peptides.py --run-dir path --method peptidebuilder
"""

import argparse
import json
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]

# Global ESMFold model cache
_esmfold_model = None


def esmfold_available() -> bool:
    """Check if ESMFold can be loaded."""
    try:
        import torch
        if not torch.cuda.is_available():
            return False
        import esm
        import esm.esmfold.v1.pretrained
        return True
    except ImportError:
        return False


def load_esmfold():
    """Load ESMFold model (downloads ~3GB on first run)."""
    global _esmfold_model
    if _esmfold_model is not None:
        return _esmfold_model

    import torch
    import esm

    log.info("Loading ESMFold model (~3GB download on first run)...")
    _esmfold_model = esm.pretrained.esmfold_v1()
    _esmfold_model = _esmfold_model.eval().cuda()
    log.info("ESMFold loaded on GPU.")
    return _esmfold_model


def fold_esmfold(sequence: str) -> tuple[str, float]:
    """
    Fold with ESMFold (GPU). Returns (pdb_string, mean_pLDDT).
    """
    import torch

    model = load_esmfold()
    with torch.no_grad():
        pdb_string = model.infer_pdb(sequence)

    # Parse pLDDT from B-factor column
    plddt_values = []
    for line in pdb_string.split("\n"):
        if line.startswith("ATOM"):
            try:
                plddt_values.append(float(line[60:66].strip()))
            except (ValueError, IndexError):
                pass

    mean_plddt = sum(plddt_values) / len(plddt_values) if plddt_values else 0.0
    return pdb_string, round(mean_plddt, 1)


def fold_peptidebuilder(sequence: str) -> tuple[str, float]:
    """
    Fallback: extended conformation with PeptideBuilder. Returns (pdb_string, 0.0).
    """
    from PeptideBuilder import make_structure
    from Bio.PDB import PDBIO
    import io

    n = len(sequence)
    phi = [-180.0] * (n - 1)
    psi_im1 = [180.0] * (n - 1)

    structure = make_structure(sequence, phi, psi_im1)
    string_io = io.StringIO()
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(string_io)
    return string_io.getvalue(), 0.0


def run(run_dir: str, method: str = "auto") -> list[dict]:
    """
    Generate 3D structures for all candidates.

    method: "auto" (try ESMFold, fall back to PeptideBuilder),
            "esmfold" (require ESMFold, fail if unavailable),
            "peptidebuilder" (always use PeptideBuilder)
    """
    sys.path.insert(0, str(ROOT))
    from modules.run_manager import RunManager

    rm = RunManager(run_dir=run_dir)
    candidates = rm.load_candidates()

    if not candidates:
        log.error("No candidates found. Run candidate_pool.py first.")
        return []

    # Decide which method to use
    use_esmfold = False
    if method == "esmfold":
        if esmfold_available():
            use_esmfold = True
        else:
            log.error("ESMFold requested but not available. Needs: pip install fair-esm openfold + GPU with CUDA.")
            return []
    elif method == "auto":
        use_esmfold = esmfold_available()

    if use_esmfold:
        log.info("Using ESMFold (GPU) — predicted folded structures with pLDDT scores")
        fold_fn = fold_esmfold
    else:
        log.info("Using PeptideBuilder — extended conformation (no GPU/ESMFold detected)")
        fold_fn = fold_peptidebuilder

    folded_dir = rm.folded_dir
    total = len(candidates)
    built = 0

    log.info(f"Generating 3D structures for {total} peptides...")

    for i, candidate in enumerate(candidates):
        cid = candidate["id"]
        seq = candidate["sequence"]
        pdb_path = folded_dir / f"{cid}.pdb"

        if pdb_path.exists() and candidate.get("folded_pdb"):
            built += 1
            continue

        try:
            pdb_string, plddt = fold_fn(seq)
            pdb_path.write_text(pdb_string)
            candidate["folded_pdb"] = str(pdb_path)
            candidate["plddt"] = plddt if plddt > 0 else None
            candidate["fold_method"] = "esmfold" if use_esmfold else "peptidebuilder"
            built += 1

            if use_esmfold and (i + 1) % 10 == 0:
                log.info(f"  [{i+1}/{total}] {built} folded (pLDDT avg so far: "
                        f"{sum(c['plddt'] for c in candidates[:i+1] if c.get('plddt')) / max(built,1):.1f})")
        except Exception as e:
            candidate["folded_pdb"] = None
            candidate["plddt"] = None
            log.error(f"  {cid} failed: {e}")

        if not use_esmfold and (i + 1) % 20 == 0:
            log.info(f"  [{i+1}/{total}] {built} structures built...")

    rm.save_candidates(candidates)

    log.info(f"\nDone. {built}/{total} structures generated.")
    log.info(f"Method: {'ESMFold (GPU)' if use_esmfold else 'PeptideBuilder (extended)'}")

    if use_esmfold:
        plddts = [c["plddt"] for c in candidates if c.get("plddt")]
        if plddts:
            log.info(f"pLDDT range: {min(plddts):.1f} - {max(plddts):.1f}, mean: {sum(plddts)/len(plddts):.1f}")

    return candidates


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--method", choices=["auto", "esmfold", "peptidebuilder"], default="auto",
                       help="auto: try ESMFold then fallback. esmfold: require GPU. peptidebuilder: always extended.")
    args = parser.parse_args()
    run(run_dir=args.run_dir, method=args.method)
