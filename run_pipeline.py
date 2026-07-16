#!/usr/bin/env python3
"""
PeptideScreen — Master Pipeline Runner
========================================
Runs the complete peptide discovery pipeline with human decision
points between phases.

Each run gets its own timestamped folder under runs/.

Usage:
    python3 run_pipeline.py                         # full pipeline
    python3 run_pipeline.py --skip-md               # skip MD (run on Colab later)
    python3 run_pipeline.py --resume <run_dir>      # resume existing run
    python3 run_pipeline.py --fast-only              # fast screen only, no full validation

Pipeline phases:
    Phase 1: Design        — ProteinMPNN + fold + BLAST
    Phase 2: Fast Screen   — HADDOCK3 rigidbody (all candidates)
    ⟹ HUMAN: review scores, select top N for full validation
    Phase 3: Full Dock     — HADDOCK3 full (top N)
    ⟹ HUMAN: review, select candidates to modify
    Phase 4: Modification  — add CPP tags, D-amino acids, etc.
    Phase 5: Properties    — physicochemical + protease + permeability
    Phase 6: MD Stability  — OpenMM (optional, can run on Colab)
    Phase 7: Report        — final ranked report
"""

import argparse
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))


def load_module(module_path: str):
    """Load a module from a file path with numbered directories."""
    full_path = ROOT / module_path
    spec = importlib.util.spec_from_file_location(full_path.stem, str(full_path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def pause(message: str) -> str:
    """Pause and ask the user what to do."""
    print(f"\n{'='*60}")
    print(f"  {message}")
    print(f"{'='*60}")
    response = input("\n  Press Enter to continue, or type 'skip' to skip: ").strip()
    return response


def main():
    parser = argparse.ArgumentParser(description="PeptideScreen — Full Pipeline")
    parser.add_argument("--resume", help="Resume an existing run directory")
    parser.add_argument("--skip-md", action="store_true", help="Skip MD stability (run on Colab later)")
    parser.add_argument("--skip-blast", action="store_true", help="Skip BLAST check (slow)")
    parser.add_argument("--fast-only", action="store_true", help="Fast screen only, no full validation")
    parser.add_argument("--haddock-n", type=int, default=20, help="Top N for full HADDOCK3 validation")
    parser.add_argument("--md-n", type=int, default=10, help="Top N for MD stability")
    parser.add_argument("--md-ns", type=float, default=50.0, help="MD simulation length (ns)")
    args = parser.parse_args()

    from modules.run_manager import RunManager

    # ================================================================
    # CREATE OR RESUME RUN
    # ================================================================
    if args.resume:
        rm = RunManager(run_dir=args.resume)
    else:
        rm = RunManager()

    run_dir = str(rm.run_dir)
    print(f"\n{'='*60}")
    print(f"  PeptideScreen Pipeline")
    print(f"  Run: {rm.run_dir.name}")
    print(f"  Directory: {run_dir}")
    print(f"{'='*60}\n")

    # ================================================================
    # PHASE 1: DESIGN
    # ================================================================
    print("\n" + "=" * 60)
    print("  PHASE 1: Peptide Design")
    print("=" * 60)

    # 1a. Structure generation
    print("\n--- Step 1/3: Generating 3D peptide structures ---")
    fold = load_module("modules/02_design/fold_peptides.py")
    fold.run()

    # 1b. BLAST check (optional)
    if not args.skip_blast:
        response = pause("Run BLAST sequence check? (~1 hour for all candidates)")
        if response.lower() != "skip":
            print("\n--- Step 2/3: BLAST sequence check ---")
            blast = load_module("modules/02_design/blast_check.py")
            blast.run()
        else:
            print("  Skipping BLAST check.")
    else:
        print("  Skipping BLAST check (--skip-blast).")

    print("\n--- Step 3/3: Design phase complete ---")

    # ================================================================
    # PHASE 2: FAST SCREEN
    # ================================================================
    response = pause("Start HADDOCK3 fast screen (all candidates)?")
    if response.lower() == "skip":
        print("  Skipping fast screen.")
    else:
        print("\n" + "=" * 60)
        print("  PHASE 2: HADDOCK3 Fast Screen")
        print("=" * 60)

        haddock = load_module("modules/03_docking/haddock3_dock.py")
        haddock.run(run_dir=run_dir, mode="fast")

        # Show ranking
        ranking = load_module("modules/03_docking/score_ranking.py")
        ranking.run(run_dir=run_dir, mode="fast")

    # ================================================================
    # PHASE 3: FULL VALIDATION
    # ================================================================
    if not args.fast_only:
        response = pause(f"Run HADDOCK3 full validation on top {args.haddock_n} candidates?")
        if response.lower() == "skip":
            print("  Skipping full validation.")
        else:
            print("\n" + "=" * 60)
            print("  PHASE 3: HADDOCK3 Full Validation")
            print("=" * 60)

            haddock = load_module("modules/03_docking/haddock3_dock.py")
            haddock.run(run_dir=run_dir, mode="full", n=args.haddock_n)

            ranking = load_module("modules/03_docking/score_ranking.py")
            ranking.run(run_dir=run_dir, mode="full")

    # ================================================================
    # PHASE 4: MODIFICATION
    # ================================================================
    response = pause("Enter modification phase? (add CPP tags, D-amino acids, etc.)")
    if response.lower() == "skip":
        print("  Skipping modification.")
    else:
        print("\n" + "=" * 60)
        print("  PHASE 4: Peptide Modification")
        print("=" * 60)

        modify = load_module("modules/04_modification/modify_peptides.py")
        modified = modify.run(run_dir=run_dir)

        if modified:
            response = pause("Re-dock modified candidates through HADDOCK3?")
            if response.lower() != "skip":
                redock = load_module("modules/04_modification/redock_modified.py")
                redock.run(run_dir=run_dir, mode="fast")

    # ================================================================
    # PHASE 5: PROPERTIES
    # ================================================================
    print("\n" + "=" * 60)
    print("  PHASE 5: Property Prediction")
    print("=" * 60)

    print("\n--- Physicochemical properties ---")
    physchem = load_module("modules/05_properties/physicochemical.py")
    physchem.run()

    print("\n--- Protease stability ---")
    protease = load_module("modules/05_properties/protease_stability.py")
    protease.run()

    print("\n--- Permeability prediction ---")
    perm = load_module("modules/05_properties/permeability.py")
    perm.run()

    print("\n--- Properties report ---")
    prop_report = load_module("modules/05_properties/properties_report.py")
    prop_report.run()

    # ================================================================
    # PHASE 6: MD STABILITY (optional)
    # ================================================================
    if not args.skip_md:
        response = pause(f"Run MD stability simulations? (top {args.md_n}, {args.md_ns} ns each)")
        if response.lower() == "skip":
            print("  Skipping MD. Run on Colab for GPU acceleration:")
            print(f"  Use the OPENMM_V2.ipynb notebook with candidates from {run_dir}")
        else:
            print("\n" + "=" * 60)
            print("  PHASE 6: MD Stability (OpenMM)")
            print("=" * 60)

            md = load_module("modules/03_docking/md_stability.py")
            md.run(run_dir=run_dir, n=args.md_n, ns=args.md_ns)
    else:
        print("\n  MD stability skipped (--skip-md).")
        print(f"  Run on Colab with OPENMM_V2.ipynb using candidates from {run_dir}")

    # ================================================================
    # PHASE 7: FINAL REPORT
    # ================================================================
    print("\n" + "=" * 60)
    print("  PHASE 7: Final Report")
    print("=" * 60)

    report = load_module("modules/03_docking/docking_report.py")
    report.run(run_dir=run_dir)

    # ================================================================
    # DONE
    # ================================================================
    print(f"\n{'='*60}")
    print(f"  PIPELINE COMPLETE")
    print(f"{'='*60}")
    print(f"\n  Run directory: {run_dir}")
    print(f"  Docking report: {run_dir}/reports/docking_report.md")
    print(f"  Properties report: outputs/reports/properties_report.md")
    print(f"  Shortlist: {run_dir}/reports/shortlist.json")
    print(f"\n  To view the report:")
    print(f"    cat {run_dir}/reports/docking_report.md")
    print(f"\n  To resume this run:")
    print(f"    python3 run_pipeline.py --resume {run_dir}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
