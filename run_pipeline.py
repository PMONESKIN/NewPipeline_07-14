#!/usr/bin/env python3
"""
PeptideScreen — Master Pipeline Runner
========================================
Runs the COMPLETE pipeline from PDB fetch through final report.
Every step reads/writes through the run directory — single data path.

Usage:
    python3 run_pipeline.py                         # full pipeline
    python3 run_pipeline.py --skip-md               # skip MD (run on Colab)
    python3 run_pipeline.py --skip-blast            # skip BLAST (slow)
    python3 run_pipeline.py --resume <run_dir>      # resume existing run
    python3 run_pipeline.py --fast-only             # fast screen only

Pipeline:
    Phase 0: Setup       — fetch PDB, analyze interface
    Phase 1: Design      — ProteinMPNN + fold + BLAST
    Phase 2: Fast Screen — HADDOCK3 rigidbody (all candidates)
    Phase 3: Full Dock   — HADDOCK3 full (top N)
    Phase 4: Modify      — CPP tags, D-amino acids, etc. + redock
    Phase 5: Properties  — physicochemical + protease + permeability
    Phase 6: MD          — OpenMM stability (optional)
    Phase 7: Report      — final ranked report
"""

import argparse
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))


def load_module(module_path: str):
    full_path = ROOT / module_path
    spec = importlib.util.spec_from_file_location(full_path.stem, str(full_path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def pause(message: str) -> str:
    print(f"\n{'='*60}")
    print(f"  {message}")
    print(f"{'='*60}")
    response = input("\n  Press Enter to continue, or type 'skip' to skip: ").strip()
    return response


def phase_header(phase_num: int, title: str):
    print(f"\n{'='*60}")
    print(f"  PHASE {phase_num}: {title}")
    print(f"{'='*60}\n")


def main():
    parser = argparse.ArgumentParser(description="PeptideScreen — Full Pipeline")
    parser.add_argument("--resume", help="Resume an existing run directory")
    parser.add_argument("--skip-md", action="store_true", help="Skip MD stability")
    parser.add_argument("--skip-blast", action="store_true", help="Skip BLAST check")
    parser.add_argument("--fast-only", action="store_true", help="Fast screen only")
    parser.add_argument("--haddock-n", type=int, default=20, help="Top N for full validation")
    parser.add_argument("--md-n", type=int, default=10, help="Top N for MD")
    parser.add_argument("--md-ns", type=float, default=50.0, help="MD length (ns)")
    args = parser.parse_args()

    from modules.run_manager import RunManager

    if args.resume:
        rm = RunManager(run_dir=args.resume)
    else:
        rm = RunManager()

    run_dir = str(rm.run_dir)

    print(f"\n{'='*60}")
    print(f"  PeptideScreen Pipeline")
    print(f"  Run: {rm.run_dir.name}")
    print(f"{'='*60}")

    # ================================================================
    # PHASE 0: TARGET SETUP
    # ================================================================
    phase_header(0, "Target Setup")

    print("--- Fetching PDB structures ---")
    fetch = load_module("modules/01_targets/fetch_structures.py")
    fetch.run(run_dir=run_dir)

    print("\n--- Analyzing binding interface ---")
    analyze = load_module("modules/01_targets/analyze_interface.py")
    analyze.run(run_dir=run_dir)

    rm.update_status("targets_ready")

    # ================================================================
    # PHASE 1: DESIGN
    # ================================================================
    phase_header(1, "Peptide Design")

    # Check if ProteinMPNN is set up
    mpnn_dir = ROOT / "tools" / "ProteinMPNN"
    if not mpnn_dir.exists():
        print("--- Setting up ProteinMPNN (one-time) ---")
        setup_mpnn = load_module("modules/02_design/setup_proteinmpnn.py")
        setup_mpnn.run()

    print("--- Running ProteinMPNN ---")
    design = load_module("modules/02_design/design_peptides.py")
    design.run(run_dir=run_dir)

    print("\n--- Building candidate pool ---")
    pool = load_module("modules/02_design/candidate_pool.py")
    pool.run(run_dir=run_dir)

    print("\n--- Generating 3D structures ---")
    fold = load_module("modules/02_design/fold_peptides.py")
    fold.run(run_dir=run_dir)

    # BLAST (optional)
    if not args.skip_blast:
        response = pause("Run BLAST sequence check? (~1 hour)")
        if response.lower() != "skip":
            print("\n--- BLAST check ---")
            blast = load_module("modules/02_design/blast_check.py")
            blast.run(run_dir=run_dir)
    else:
        print("  Skipping BLAST (--skip-blast)")

    rm.update_status("design_complete")

    # ================================================================
    # PHASE 2: FAST SCREEN
    # ================================================================
    response = pause("Start HADDOCK3 fast screen? (all candidates, ~3 min each)")
    if response.lower() != "skip":
        phase_header(2, "HADDOCK3 Fast Screen")

        haddock = load_module("modules/03_docking/haddock3_dock.py")
        haddock.run(run_dir=run_dir, mode="fast")

        ranking = load_module("modules/03_docking/score_ranking.py")
        ranking.run(run_dir=run_dir, mode="fast")

    # ================================================================
    # PHASE 3: FULL VALIDATION
    # ================================================================
    if not args.fast_only:
        response = pause(f"Run HADDOCK3 full validation on top {args.haddock_n}?")
        if response.lower() != "skip":
            phase_header(3, "HADDOCK3 Full Validation")

            haddock = load_module("modules/03_docking/haddock3_dock.py")
            haddock.run(run_dir=run_dir, mode="full", n=args.haddock_n)

            ranking = load_module("modules/03_docking/score_ranking.py")
            ranking.run(run_dir=run_dir, mode="full")

    # ================================================================
    # PHASE 4: MODIFICATION
    # ================================================================
    response = pause("Enter modification phase? (add CPP tags, D-amino acids, etc.)")
    if response.lower() != "skip":
        phase_header(4, "Peptide Modification")

        modify = load_module("modules/04_modification/modify_peptides.py")
        modified = modify.run(run_dir=run_dir)

        if modified:
            response = pause("Re-dock modified candidates?")
            if response.lower() != "skip":
                redock = load_module("modules/04_modification/redock_modified.py")
                redock.run(run_dir=run_dir, mode="fast")

    # ================================================================
    # PHASE 5: PROPERTIES
    # ================================================================
    phase_header(5, "Property Prediction")

    print("--- Physicochemical ---")
    physchem = load_module("modules/05_properties/physicochemical.py")
    physchem.run(run_dir=run_dir)

    print("\n--- Protease stability ---")
    protease = load_module("modules/05_properties/protease_stability.py")
    protease.run(run_dir=run_dir)

    print("\n--- Permeability ---")
    perm = load_module("modules/05_properties/permeability.py")
    perm.run(run_dir=run_dir)

    print("\n--- Properties report ---")
    prop_report = load_module("modules/05_properties/properties_report.py")
    prop_report.run(run_dir=run_dir)

    # ================================================================
    # PHASE 6: MD STABILITY
    # ================================================================
    if not args.skip_md:
        response = pause(f"Run MD stability? (top {args.md_n}, {args.md_ns} ns each)")
        if response.lower() != "skip":
            phase_header(6, "MD Stability (OpenMM)")
            try:
                md = load_module("modules/03_docking/md_stability.py")
                md.run(run_dir=run_dir, n=args.md_n, ns=args.md_ns)
            except Exception as e:
                print(f"  MD skipped: {e}")
                print("  Run on Colab with OPENMM_V2.ipynb instead.")
    else:
        print(f"\n  MD skipped (--skip-md). Run on Colab with OPENMM_V2.ipynb")
        print(f"  Candidates at: {run_dir}/candidates/candidate_pool.json")

    # ================================================================
    # PHASE 7: FINAL REPORT
    # ================================================================
    phase_header(7, "Final Report")

    report = load_module("modules/03_docking/docking_report.py")
    report.run(run_dir=run_dir)

    # ================================================================
    # DONE
    # ================================================================
    print(f"\n{'='*60}")
    print(f"  PIPELINE COMPLETE")
    print(f"{'='*60}")
    print(f"  Run: {rm.run_dir.name}")
    print(f"  Directory: {run_dir}")
    print(f"  Reports:")
    print(f"    {run_dir}/reports/docking_report.md")
    print(f"    {run_dir}/reports/properties_report.md")
    print(f"    {run_dir}/reports/shortlist.json")
    print(f"  Resume: python3 run_pipeline.py --resume {run_dir}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
