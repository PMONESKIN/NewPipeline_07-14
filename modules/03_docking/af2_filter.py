"""
Module 03 — AF2-Multimer Filter
=================================
Runs ColabFold (AlphaFold2-Multimer) on each candidate to predict
the peptide-receptor complex structure. Filters candidates by
confidence metrics: ipTM, pLDDT, and interface PAE.

This replaces the HADDOCK3 fast screen with a more accurate,
structure-based filter. Output PDB structures are used as input
for HADDOCK3 full docking.

Usage:
    python3 modules/03_docking/af2_filter.py --run-dir path

Requirements:
    pip install colabfold[alphafold]  (on Colab with GPU)
"""

import argparse
import json
import logging
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]

# Filter thresholds (based on published binder design benchmarks)
IPTM_THRESHOLD = 0.7
PLDDT_THRESHOLD = 80.0
PAE_THRESHOLD = 7.0


def get_receptor_sequence(interface: dict, structures_dir: Path) -> str:
    """Extract receptor amino acid sequence from PDB."""
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1

    pdb_path = structures_dir / f"{interface['pdb_id']}.pdb"
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))

    for model in structure:
        for chain in model:
            if chain.id == interface["receptor_chain"]:
                seq = "".join(
                    protein_letters_3to1.get(r.get_resname(), "X")
                    for r in chain if is_aa(r, standard=True)
                )
                return seq
    return ""


def write_colabfold_input(candidates: list[dict], receptor_seq: str,
                          output_dir: Path) -> Path:
    """
    Write a FASTA file with receptor:peptide pairs for ColabFold batch.
    ColabFold expects sequences separated by ':' for multimer prediction.
    """
    fasta_path = output_dir / "af2_input.fasta"
    with open(fasta_path, "w") as f:
        for c in candidates:
            # ColabFold format: >name\nchain1:chain2
            f.write(f">{c['id']}\n")
            f.write(f"{receptor_seq}:{c['sequence']}\n")

    log.info(f"  Wrote {len(candidates)} sequences to {fasta_path.name}")
    return fasta_path


def run_colabfold(fasta_path: Path, output_dir: Path,
                  num_models: int = 1, num_recycles: int = 3) -> bool:
    """Run ColabFold batch on the input FASTA."""
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "colabfold_batch",
        str(fasta_path),
        str(output_dir),
        "--num-models", str(num_models),
        "--num-recycle", str(num_recycles),
        "--msa-mode", "single_sequence",  # faster, sufficient for short peptides
        "--model-type", "alphafold2_multimer_v3",
    ]

    log.info(f"  Running ColabFold ({len(open(fasta_path).readlines())//2} complexes)...")
    log.info(f"  This may take 1-2 minutes per candidate on GPU...")

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)

    if result.returncode != 0:
        log.error(f"  ColabFold failed: {result.stderr[:500]}")
        return False

    return True


def parse_colabfold_results(output_dir: Path, candidates: list[dict]) -> list[dict]:
    """
    Parse ColabFold output files for each candidate.
    Extracts ipTM, pLDDT, PAE, and saves PDB paths.
    """
    results = []

    for c in candidates:
        cid = c["id"]

        # ColabFold output naming: {id}_unrelaxed_rank_001_*.pdb
        pdb_files = sorted(output_dir.glob(f"{cid}_unrelaxed_rank_001_*.pdb"))
        json_files = sorted(output_dir.glob(f"{cid}_scores_rank_001_*.json"))

        if not pdb_files:
            # Try alternative naming
            pdb_files = sorted(output_dir.glob(f"{cid}*.pdb"))
            json_files = sorted(output_dir.glob(f"{cid}*scores*.json"))

        if not pdb_files:
            log.warning(f"  {cid}: no PDB output found")
            c["af2_pdb"] = None
            c["af2_iptm"] = None
            c["af2_plddt"] = None
            c["af2_pae"] = None
            c["af2_pass"] = False
            continue

        # Best PDB
        best_pdb = pdb_files[0]
        c["af2_pdb"] = str(best_pdb)

        # Parse pLDDT from PDB B-factor column (peptide chain only)
        plddt_values = []
        with open(best_pdb) as f:
            for line in f:
                if line.startswith("ATOM"):
                    try:
                        plddt_values.append(float(line[60:66].strip()))
                    except (ValueError, IndexError):
                        pass

        c["af2_plddt"] = round(np.mean(plddt_values), 1) if plddt_values else None

        # Parse scores JSON if available
        if json_files:
            try:
                with open(json_files[0]) as f:
                    scores = json.load(f)
                c["af2_iptm"] = round(float(scores.get("iptm", 0)), 3)
                c["af2_ptm"] = round(float(scores.get("ptm", 0)), 3)

                # Calculate interface PAE (mean PAE between chains)
                if "pae" in scores:
                    pae_matrix = np.array(scores["pae"])
                    # Receptor is the longer chain, peptide is shorter
                    # PAE between chains = off-diagonal blocks
                    receptor_len = pae_matrix.shape[0] - len(c["sequence"])
                    peptide_len = len(c["sequence"])
                    if receptor_len > 0:
                        # Extract cross-chain PAE block
                        cross_pae = pae_matrix[:receptor_len, receptor_len:]
                        c["af2_pae"] = round(float(np.mean(cross_pae)), 2)
                    else:
                        c["af2_pae"] = None
                else:
                    c["af2_pae"] = None

            except (json.JSONDecodeError, KeyError, ValueError) as e:
                log.warning(f"  {cid}: score parsing failed — {e}")
                c["af2_iptm"] = None
                c["af2_pae"] = None
        else:
            c["af2_iptm"] = None
            c["af2_pae"] = None

        # Apply filter
        iptm_pass = c.get("af2_iptm") is not None and c["af2_iptm"] >= IPTM_THRESHOLD
        plddt_pass = c.get("af2_plddt") is not None and c["af2_plddt"] >= PLDDT_THRESHOLD
        pae_pass = c.get("af2_pae") is None or c["af2_pae"] <= PAE_THRESHOLD

        c["af2_pass"] = iptm_pass and plddt_pass and pae_pass

        status = "PASS" if c["af2_pass"] else "FAIL"
        iptm_str = f"{c['af2_iptm']:.3f}" if c.get("af2_iptm") is not None else "N/A"
        plddt_str = f"{c['af2_plddt']:.1f}" if c.get("af2_plddt") is not None else "N/A"
        pae_str = f"{c['af2_pae']:.1f}" if c.get("af2_pae") is not None else "N/A"

        log.info(f"  {cid}: ipTM={iptm_str} pLDDT={plddt_str} PAE={pae_str} → {status}")

    return candidates


def extract_peptide_pdb(af2_complex_pdb: Path, peptide_len: int,
                        output_path: Path) -> Path:
    """
    Extract the peptide chain from an AF2 complex prediction.
    The peptide is the last chain (shorter one).
    """
    from Bio.PDB import PDBParser, PDBIO, Select
    from Bio.PDB.Polypeptide import is_aa

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", str(af2_complex_pdb))

    # Find the shortest chain (that's the peptide)
    chains = []
    for model in structure:
        for chain in model:
            residues = [r for r in chain if is_aa(r, standard=True)]
            chains.append((chain.id, len(residues)))

    chains.sort(key=lambda x: x[1])
    peptide_chain_id = chains[0][0]  # shortest chain

    class PeptideSelect(Select):
        def accept_chain(self, chain):
            return chain.id == peptide_chain_id

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path), PeptideSelect())

    return output_path


def run(run_dir: str, iptm_cutoff: float = None, plddt_cutoff: float = None,
        pae_cutoff: float = None) -> list[dict]:
    """
    Run AF2-Multimer filter on all candidates.
    """
    sys.path.insert(0, str(ROOT))
    from modules.run_manager import RunManager

    rm = RunManager(run_dir=run_dir)
    config = rm.config
    candidates = rm.load_candidates()

    if not candidates:
        log.error("No candidates found.")
        return []

    # Override thresholds if specified
    global IPTM_THRESHOLD, PLDDT_THRESHOLD, PAE_THRESHOLD
    if iptm_cutoff is not None:
        IPTM_THRESHOLD = iptm_cutoff
    if plddt_cutoff is not None:
        PLDDT_THRESHOLD = plddt_cutoff
    if pae_cutoff is not None:
        PAE_THRESHOLD = pae_cutoff

    # Get receptor sequence
    target_cfg = config["targets"][0]
    interface = rm.load_interface(target_cfg["name"])
    receptor_seq = get_receptor_sequence(interface, rm.structures_dir)

    if not receptor_seq:
        log.error("Could not extract receptor sequence.")
        return candidates

    log.info(f"Receptor: {len(receptor_seq)} aa")
    log.info(f"Candidates: {len(candidates)}")
    log.info(f"Thresholds: ipTM>{IPTM_THRESHOLD} pLDDT>{PLDDT_THRESHOLD} PAE<{PAE_THRESHOLD}")

    # Set up AF2 output directory
    af2_dir = rm.docking_dir / "af2_predictions"
    af2_dir.mkdir(parents=True, exist_ok=True)

    # Write input FASTA
    fasta_path = write_colabfold_input(candidates, receptor_seq, af2_dir)

    # Run ColabFold
    af2_output = af2_dir / "output"
    success = run_colabfold(fasta_path, af2_output)

    if not success:
        log.error("ColabFold failed. Check installation: pip install colabfold[alphafold]")
        return candidates

    # Parse results and filter
    candidates = parse_colabfold_results(af2_output, candidates)

    # Save AF2-predicted peptide PDBs for HADDOCK input
    for c in candidates:
        if c.get("af2_pass") and c.get("af2_pdb"):
            peptide_pdb = rm.folded_dir / f"{c['id']}.pdb"
            try:
                extract_peptide_pdb(
                    Path(c["af2_pdb"]),
                    len(c["sequence"]),
                    peptide_pdb,
                )
                c["folded_pdb"] = str(peptide_pdb)
                c["fold_method"] = "af2_multimer"
            except Exception as e:
                log.warning(f"  {c['id']}: peptide extraction failed — {e}")

    # Save updated candidates
    rm.save_candidates(candidates)

    # Summary
    passed = sum(1 for c in candidates if c.get("af2_pass"))
    failed = len(candidates) - passed

    log.info(f"\n{'='*60}")
    log.info(f"AF2-Multimer Filter Results")
    log.info(f"{'='*60}")
    log.info(f"  Total:  {len(candidates)}")
    log.info(f"  Passed: {passed}")
    log.info(f"  Failed: {failed}")
    log.info(f"  Pass rate: {100*passed/max(len(candidates),1):.0f}%")

    if passed > 0:
        log.info(f"\n  Passed candidates:")
        passed_list = [c for c in candidates if c.get("af2_pass")]
        passed_list.sort(key=lambda c: c.get("af2_iptm", 0), reverse=True)
        for c in passed_list:
            log.info(f"    {c['id']}  ipTM={c.get('af2_iptm',0):.3f}  "
                    f"pLDDT={c.get('af2_plddt',0):.1f}  {c['sequence']}")

    rm.update_status("af2_filtered")
    return candidates


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--iptm-cutoff", type=float, default=None)
    parser.add_argument("--plddt-cutoff", type=float, default=None)
    parser.add_argument("--pae-cutoff", type=float, default=None)
    args = parser.parse_args()
    run(run_dir=args.run_dir, iptm_cutoff=args.iptm_cutoff,
        plddt_cutoff=args.plddt_cutoff, pae_cutoff=args.pae_cutoff)
