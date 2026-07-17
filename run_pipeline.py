#!/usr/bin/env python3
"""
PeptideScreen — Master Pipeline Runner
========================================
State-of-the-art peptide discovery pipeline.

Flow:
  Phase 0: Target Setup (fetch PDB, analyze interface)
  Phase 1: Backbone (extract from PDB or RFdiffusion)
  Phase 2: Design (ProteinMPNN + candidate pool)
  Phase 3: AF2-Multimer Filter (predict complex, filter by ipTM/pLDDT/PAE)
  Phase 4: HADDOCK3 Full Docking (on AF2-validated candidates only)
  Phase 5: Modification (CPP tags, D-amino acids, redock)
  Phase 6: Properties (physicochemical, protease, permeability)
  Phase 7: MD Stability (OpenMM, optional)
  Phase 8: Final Report

Usage:
    python3 run_pipeline.py --no-pause
    python3 run_pipeline.py --no-pause --use-rfdiffusion --peptide-length 8
    python3 run_pipeline.py --resume <run_dir>
"""

import argparse
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

_no_pause = False


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
    if _no_pause:
        print("  [--no-pause] Continuing automatically...")
        return ""
    response = input("\n  Press Enter to continue, or type 'skip' to skip: ").strip()
    return response


def phase_header(phase_num: int, title: str):
    print(f"\n{'='*60}")
    print(f"  PHASE {phase_num}: {title}")
    print(f"{'='*60}\n")


def main():
    parser = argparse.ArgumentParser(description="PeptideScreen — Full Pipeline")
    parser.add_argument("--resume", help="Resume an existing run directory")
    parser.add_argument("--no-pause", action="store_true", help="Run without pause points")
    parser.add_argument("--skip-md", action="store_true", help="Skip MD stability")
    parser.add_argument("--skip-blast", action="store_true", help="Skip BLAST check")
    parser.add_argument("--use-rfdiffusion", action="store_true", help="Use RFdiffusion for backbone generation")
    parser.add_argument("--peptide-length", type=int, default=None, help="Target peptide length (for backbone extraction)")
    parser.add_argument("--n-backbones", type=int, default=100, help="Number of RFdiffusion backbones")
    parser.add_argument("--haddock-n", type=int, default=None, help="Limit HADDOCK3 to top N AF2-passed candidates")
    parser.add_argument("--md-n", type=int, default=5, help="Top N for MD stability")
    parser.add_argument("--md-ns", type=float, default=50.0, help="MD simulation length (ns)")
    parser.add_argument("--iptm-cutoff", type=float, default=0.7, help="AF2 ipTM filter threshold")
    parser.add_argument("--plddt-cutoff", type=float, default=80.0, help="AF2 pLDDT filter threshold")
    args = parser.parse_args()

    global _no_pause
    _no_pause = args.no_pause

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
    print(f"  PeptideScreen Pipeline v2.0")
    print(f"  Run: {rm.run_dir.name}")
    print(f"  Directory: {run_dir}")
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
    # PHASE 1: BACKBONE GENERATION
    # ================================================================
    phase_header(1, "Backbone Generation")

    if args.use_rfdiffusion:
        print("--- RFdiffusion backbone generation ---")
        response = pause(f"Generate {args.n_backbones} diverse backbones with RFdiffusion?")
        if response.lower() != "skip":
            try:
                rfdiff = load_module("modules/02_design/rfdiffusion_backbones.py")
                length = args.peptide_length or 8
                rfdiff.run(run_dir=run_dir, n_backbones=args.n_backbones,
                          peptide_length=length)
            except Exception as e:
                print(f"  RFdiffusion failed: {e}")
                print("  Falling back to co-crystal backbone extraction...")
                backbone = load_module("modules/02_design/backbone_extract.py")
                backbone.run(run_dir=run_dir, length=args.peptide_length)
    else:
        print("--- Extracting backbone from co-crystal structure ---")
        backbone = load_module("modules/02_design/backbone_extract.py")
        backbone.run(run_dir=run_dir, length=args.peptide_length)

    rm.update_status("backbone_ready")

    # ================================================================
    # PHASE 2: SEQUENCE DESIGN
    # ================================================================
    phase_header(2, "Sequence Design (ProteinMPNN)")

    # Check ProteinMPNN setup
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

    # BLAST check (optional)
    if not args.skip_blast:
        response = pause("Run BLAST sequence check? (~30 sec per candidate)")
        if response.lower() != "skip":
            print("\n--- BLAST check ---")
            blast = load_module("modules/02_design/blast_check.py")
            blast.run(run_dir=run_dir)
    else:
        print("  Skipping BLAST (--skip-blast)")

    rm.update_status("design_complete")

    # ================================================================
    # PHASE 3: AF2-MULTIMER FILTER
    # ================================================================
    phase_header(3, "AF2-Multimer Structure Prediction + Filter")

    print("--- Running ColabFold (AlphaFold2-Multimer) ---")
    print("--- This predicts complex structure and filters by ipTM/pLDDT/PAE ---")

    af2 = load_module("modules/03_docking/af2_filter.py")
    af2.run(run_dir=run_dir, iptm_cutoff=args.iptm_cutoff,
            plddt_cutoff=args.plddt_cutoff)

    # Check how many passed
    candidates = rm.load_candidates()
    passed = [c for c in candidates if c.get("af2_pass")]
    print(f"\n  AF2 filter: {len(passed)}/{len(candidates)} candidates passed")

    if not passed:
        print("\n  WARNING: No candidates passed AF2 filter!")
        print("  Consider lowering thresholds: --iptm-cutoff 0.5 --plddt-cutoff 60")
        response = pause("Continue anyway with all candidates?")
        if response.lower() == "skip":
            print("  Pipeline stopped.")
            return

    # ================================================================
    # PHASE 4: HADDOCK3 FULL DOCKING
    # ================================================================
    phase_header(4, "HADDOCK3 Full Docking (AF2-validated candidates)")

    # Only dock candidates that passed AF2 filter
    n_to_dock = args.haddock_n or len(passed)

    response = pause(f"Run HADDOCK3 full docking on {min(n_to_dock, len(passed))} AF2-validated candidates?")
    if response.lower() != "skip":
        haddock = load_module("modules/03_docking/haddock3_dock.py")

        # Mark non-AF2-passed candidates so HADDOCK skips them
        for c in candidates:
            if not c.get("af2_pass"):
                c["haddock_full_score"] = None  # ensure they're skipped

        rm.save_candidates(candidates)
        haddock.run(run_dir=run_dir, mode="full", n=n_to_dock)

        ranking = load_module("modules/03_docking/score_ranking.py")
        ranking.run(run_dir=run_dir, mode="full")

    # ================================================================
    # PHASE 5: MODIFICATION
    # ================================================================
    response = pause("Enter modification phase? (add CPP tags, D-amino acids, etc.)")
    if response.lower() != "skip":
        phase_header(5, "Peptide Modification")

        modify = load_module("modules/04_modification/modify_peptides.py")
        modified = modify.run(run_dir=run_dir)

        if modified:
            response = pause("Re-dock modified candidates through HADDOCK3?")
            if response.lower() != "skip":
                redock = load_module("modules/04_modification/redock_modified.py")
                redock.run(run_dir=run_dir, mode="full")

    # ================================================================
    # PHASE 6: PROPERTIES
    # ================================================================
    phase_header(6, "Property Prediction")

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
    # PHASE 7: MD STABILITY
    # ================================================================
    if not args.skip_md:
        response = pause(f"Run MD stability? (top {args.md_n}, {args.md_ns} ns each)")
        if response.lower() != "skip":
            phase_header(7, "MD Stability (OpenMM)")
            try:
                md = load_module("modules/03_docking/md_stability.py")
                md.run(run_dir=run_dir, n=args.md_n, ns=args.md_ns)
            except Exception as e:
                print(f"  MD skipped: {e}")
                print("  Run on Colab with OPENMM_V2.ipynb for GPU acceleration.")
    else:
        print(f"\n  MD skipped (--skip-md). Run separately if needed.")

    # ================================================================
    # PHASE 8: FINAL REPORT
    # ================================================================
    phase_header(8, "Final Report")

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
    print(f"  Resume: python3 run_pipeline.py --resume {run_dir}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
