"""
Module 04 — Step 1: Peptide Modification
=========================================
Interactive CLI for applying modifications to top docking candidates.

Supported modifications:
  - CPP tags: prepend cell-penetrating peptide sequences (R4, R8, TAT, penetratin)
  - D-amino acids: swap specific positions to D-form (protease resistance)
  - Truncation: shorten the peptide from either end
  - Cyclization: flag for head-to-tail cyclization
  - Custom: user types any modification

Usage:
    python3 modules/04_modification/modify_peptides.py --run-dir path

Inputs:
    {run_dir}/candidates/candidate_pool.json (with HADDOCK scores)

Outputs:
    {run_dir}/candidates/modified_candidates.json
    data/candidates/folded_structures/{modified_id}.pdb
"""

import json
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]

CPP_TAGS = {
    "R4": "RRRR",
    "R8": "RRRRRRRR",
    "R9": "RRRRRRRRR",
    "TAT": "GRKKRRQRRRPQ",
    "penetratin": "RQIKIWFQNRRMKWKK",
}


def load_config(config_path: str = None) -> dict:
    import yaml
    path = Path(config_path) if config_path else ROOT / "config.yaml"
    with open(path) as f:
        return yaml.safe_load(f)


def build_structure(sequence: str, output_path: Path):
    """Build a PDB for the modified peptide."""
    from PeptideBuilder import make_structure
    from Bio.PDB import PDBIO

    n = len(sequence)
    phi = [-180.0] * (n - 1)
    psi_im1 = [180.0] * (n - 1)

    structure = make_structure(sequence, phi, psi_im1)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path))


def show_top_candidates(candidates: list[dict], n: int = 20):
    """Display top candidates for user selection."""
    # Try full score first, then fast score
    score_field = None
    for field in ["haddock_full_score", "haddock_fast_score"]:
        if any(c.get(field) is not None for c in candidates):
            score_field = field
            break

    if score_field:
        scored = [c for c in candidates if c.get(score_field) is not None]
        scored.sort(key=lambda c: c[score_field])
    else:
        scored = candidates

    print(f"\n{'Rank':<6} {'ID':<16} {'Sequence':<25} {'Score':<10}")
    print(f"{'-'*6} {'-'*16} {'-'*25} {'-'*10}")

    for rank, c in enumerate(scored[:n], 1):
        seq = c["sequence"]
        score = f"{c[score_field]:.1f}" if score_field and c.get(score_field) is not None else "N/A"
        print(f"{rank:<6} {c['id']:<16} {seq:<25} {score:<10}")


def apply_cpp_tag(sequence: str, tag_name: str, position: str = "N") -> str:
    """Prepend or append a CPP tag."""
    tag_seq = CPP_TAGS.get(tag_name)
    if not tag_seq:
        print(f"  Unknown tag '{tag_name}'. Available: {list(CPP_TAGS.keys())}")
        return sequence

    if position.upper() == "N":
        return tag_seq + sequence
    else:
        return sequence + tag_seq


def apply_d_amino_acids(sequence: str, positions: list[int]) -> tuple[str, str]:
    """
    Mark positions as D-amino acids. Returns modified sequence and a note.
    The sequence itself doesn't change (same letters), but we track which
    positions are D-form for downstream tools.
    """
    note = f"D-amino acids at positions: {positions}"
    return sequence, note


def apply_truncation(sequence: str, start: int, end: int) -> str:
    """Truncate peptide to positions start:end (1-indexed, inclusive)."""
    return sequence[start - 1:end]


def interactive_modify(candidates: list[dict], run_dir: Path) -> list[dict]:
    """Interactive modification session."""
    modified = []
    folded_dir = ROOT / "data" / "candidates" / "folded_structures"
    folded_dir.mkdir(parents=True, exist_ok=True)

    print("\n=== Peptide Modification ===")
    show_top_candidates(candidates)

    while True:
        print("\n" + "=" * 50)
        candidate_id = input("\nEnter candidate ID to modify (or 'done' to finish): ").strip()

        if candidate_id.lower() == "done":
            break

        # Find the candidate
        candidate = None
        for c in candidates:
            if c["id"] == candidate_id:
                candidate = c
                break

        if not candidate:
            print(f"  Candidate '{candidate_id}' not found.")
            continue

        original_seq = candidate["sequence"]
        print(f"\n  Original: {original_seq} ({len(original_seq)} aa)")

        modified_seq = original_seq
        modifications = []

        # Modification loop for this candidate
        while True:
            print(f"\n  Current:  {modified_seq} ({len(modified_seq)} aa)")
            print(f"\n  Modifications available:")
            print(f"    1. Add CPP tag (R4, R8, R9, TAT, penetratin)")
            print(f"    2. D-amino acid substitution")
            print(f"    3. Truncate")
            print(f"    4. Custom prepend/append")
            print(f"    5. Done with this candidate")

            choice = input("\n  Choice [1-5]: ").strip()

            if choice == "1":
                print(f"\n  Available tags:")
                for name, seq in CPP_TAGS.items():
                    print(f"    {name}: {seq} ({len(seq)} aa)")
                tag = input("  Tag name: ").strip()
                pos = input("  Position [N-term/C-term] (N): ").strip() or "N"
                modified_seq = apply_cpp_tag(modified_seq, tag, pos)
                modifications.append(f"CPP:{tag}:{pos}-term")
                print(f"  → {modified_seq}")

            elif choice == "2":
                positions_str = input("  Positions to swap to D-form (comma-separated): ").strip()
                try:
                    positions = [int(p.strip()) for p in positions_str.split(",")]
                    modified_seq, note = apply_d_amino_acids(modified_seq, positions)
                    modifications.append(f"D-aa:{positions}")
                    print(f"  → {note}")
                except ValueError:
                    print("  Invalid positions.")

            elif choice == "3":
                start = input(f"  Start position [1-{len(modified_seq)}] (1): ").strip()
                end = input(f"  End position [1-{len(modified_seq)}] ({len(modified_seq)}): ").strip()
                start = int(start) if start else 1
                end = int(end) if end else len(modified_seq)
                modified_seq = apply_truncation(modified_seq, start, end)
                modifications.append(f"truncate:{start}-{end}")
                print(f"  → {modified_seq}")

            elif choice == "4":
                custom = input("  Sequence to add: ").strip().upper()
                pos = input("  Position [N-term/C-term] (N): ").strip() or "N"
                if pos.upper() == "N":
                    modified_seq = custom + modified_seq
                else:
                    modified_seq = modified_seq + custom
                modifications.append(f"custom:{custom}:{pos}-term")
                print(f"  → {modified_seq}")

            elif choice == "5":
                break

        if modifications:
            mod_id = f"{candidate_id}_mod"
            mod_entry = {
                "id": mod_id,
                "original_id": candidate_id,
                "original_sequence": original_seq,
                "sequence": modified_seq,
                "length": len(modified_seq),
                "target_name": candidate["target_name"],
                "design_source": "modified",
                "modifications": modifications,
                "haddock_fast_score": None,
                "haddock_full_score": None,
            }

            # Build PDB for modified peptide
            pdb_path = folded_dir / f"{mod_id}.pdb"
            try:
                build_structure(modified_seq, pdb_path)
                mod_entry["folded_pdb"] = str(pdb_path)
                print(f"\n  PDB generated: {pdb_path.name}")
            except Exception as e:
                print(f"\n  Warning: PDB generation failed: {e}")
                mod_entry["folded_pdb"] = None

            modified.append(mod_entry)
            print(f"\n  Saved: {mod_id}")
            print(f"  {original_seq} → {modified_seq}")
            print(f"  Modifications: {', '.join(modifications)}")

            another = input("\n  Modify another candidate? [y/N]: ").strip()
            if another.lower() != "y":
                break

    return modified


def run(run_dir: str) -> list[dict]:
    """Run interactive modification session."""
    sys.path.insert(0, str(ROOT))
    from modules.run_manager import RunManager

    rm = RunManager(run_dir=run_dir)
    candidates = rm.load_candidates()

    modified = interactive_modify(candidates, rm.run_dir)

    if modified:
        # Save modified candidates
        mod_path = rm.candidates_dir / "modified_candidates.json"
        with open(mod_path, "w") as f:
            json.dump(modified, f, indent=2)

        print(f"\n{len(modified)} modified candidate(s) saved to {mod_path}")
        print(f"\nNext step: re-dock modified candidates:")
        print(f"  python3 modules/04_modification/redock_modified.py --run-dir {run_dir}")
    else:
        print("\nNo modifications made.")

    return modified


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True, help="Run directory")
    args = parser.parse_args()
    run(run_dir=args.run_dir)
