"""
Module 02 — Step 1: ProteinMPNN Peptide Design
===============================================
Runs ProteinMPNN to generate novel peptide sequences.
Reads interface JSONs from the run directory.

Usage:
    python3 modules/02_design/design_peptides.py --run-dir path
"""

import argparse
import json
import logging
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
MPNN_DIR = ROOT / "tools" / "ProteinMPNN"


def get_device_flag():
    try:
        import torch
        if torch.cuda.is_available():
            return "--use_gpu"
    except ImportError:
        pass
    return ""


def find_motif_positions(pdb_path, chain_id, motif):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                seq = "".join(protein_letters_3to1.get(r.get_resname(), "X")
                             for r in chain if is_aa(r, standard=True))
                idx = seq.find(motif)
                if idx >= 0:
                    positions = list(range(idx + 1, idx + 1 + len(motif)))
                    log.info(f"  Motif '{motif}' at positions {positions}")
                    return positions
                else:
                    log.warning(f"  Motif '{motif}' NOT found in chain {chain_id}")
    return []


def run_proteinmpnn(pdb_path, output_dir, fixed_chains, design_chains,
                    num_sequences=30, temperatures=None, seed=42,
                    fixed_positions=None):
    if not MPNN_DIR.exists():
        raise FileNotFoundError(f"ProteinMPNN not found at {MPNN_DIR}. Run setup_proteinmpnn.py first.")

    temperatures = temperatures or [0.1, 0.2, 0.3]
    output_dir.mkdir(parents=True, exist_ok=True)

    parsed_dir = output_dir / "parsed"
    parsed_dir.mkdir(exist_ok=True)
    parsed_jsonl = parsed_dir / "parsed.jsonl"

    with tempfile.TemporaryDirectory() as tmp:
        import shutil
        shutil.copy(pdb_path, Path(tmp) / pdb_path.name)
        subprocess.run([sys.executable, str(MPNN_DIR / "helper_scripts" / "parse_multiple_chains.py"),
                       "--input_path", tmp, "--output_path", str(parsed_jsonl)],
                      capture_output=True, text=True, check=True)

    assigned_jsonl = parsed_dir / "assigned.jsonl"
    subprocess.run([sys.executable, str(MPNN_DIR / "helper_scripts" / "assign_fixed_chains.py"),
                   "--input_path", str(parsed_jsonl), "--output_path", str(assigned_jsonl),
                   "--chain_list", design_chains], capture_output=True, text=True, check=True)

    fixed_pos_jsonl = None
    if fixed_positions:
        fixed_pos_jsonl = parsed_dir / "fixed_positions.jsonl"
        pos_dict = {pdb_path.stem: {design_chains: fixed_positions}}
        for fc in fixed_chains.split(","):
            fc = fc.strip()
            if fc:
                pos_dict[pdb_path.stem][fc] = []
        with open(fixed_pos_jsonl, "w") as f:
            f.write(json.dumps(pos_dict) + "\n")

    temp_str = " ".join(str(t) for t in temperatures)
    cmd = [sys.executable, str(MPNN_DIR / "protein_mpnn_run.py"),
           "--jsonl_path", str(parsed_jsonl), "--chain_id_jsonl", str(assigned_jsonl),
           "--out_folder", str(output_dir), "--num_seq_per_target", str(num_sequences),
           "--sampling_temp", temp_str, "--seed", str(seed), "--batch_size", "1"]
    if fixed_pos_jsonl:
        cmd.extend(["--fixed_positions_jsonl", str(fixed_pos_jsonl)])
    device_flag = get_device_flag()
    if device_flag:
        cmd.append(device_flag)

    log.info(f"  Running ProteinMPNN ({num_sequences} seqs x {len(temperatures)} temps)...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        # If GPU failed, retry on CPU
        if device_flag and "--use_gpu" in cmd:
            log.warning(f"  GPU failed, retrying on CPU...")
            log.warning(f"  Error: {result.stderr[:200]}")
            cmd.remove("--use_gpu")
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                raise RuntimeError(f"ProteinMPNN failed on CPU too:\n{result.stderr[:500]}")
        else:
            raise RuntimeError(f"ProteinMPNN failed:\n{result.stderr[:500]}")

    fasta_files = list((output_dir / "seqs").glob("*.fa"))
    if not fasta_files:
        raise FileNotFoundError("No FASTA output found")
    return fasta_files[0]


def parse_mpnn_fasta(fasta_path, target_name, design_chain_len=None):
    entries = []
    with open(fasta_path) as f:
        lines = f.read().strip().split("\n")
    i = 0
    while i < len(lines):
        if lines[i].startswith(">"):
            entries.append((lines[i][1:], lines[i+1].strip() if i+1 < len(lines) else ""))
            i += 2
        else:
            i += 1

    if not entries:
        return []

    _, ref_seq = entries[0]
    if design_chain_len is None:
        design_chain_len = len(ref_seq)

    candidates = []
    for header, seq in entries[1:]:
        score = temp = seq_recovery = None
        for part in header.split(","):
            part = part.strip()
            if part.startswith("score="):
                try: score = float(part.split("=")[1])
                except: pass
            elif part.startswith("T="):
                try: temp = float(part.split("=")[1])
                except: pass
            elif part.startswith("seq_recovery="):
                try: seq_recovery = float(part.split("=")[1])
                except: pass

        candidates.append({
            "sequence": seq[:design_chain_len].upper(),
            "length": min(len(seq), design_chain_len),
            "mpnn_score": score,
            "sampling_temperature": temp,
            "sequence_recovery": seq_recovery,
            "design_source": "proteinmpnn",
            "target_name": target_name,
        })
    return candidates


def run(run_dir: str) -> list[dict]:
    sys.path.insert(0, str(ROOT))
    from modules.run_manager import RunManager

    rm = RunManager(run_dir=run_dir)
    config = rm.config
    all_candidates = []

    mpnn_output_dir = rm.candidates_dir / "mpnn_output"
    mpnn_output_dir.mkdir(parents=True, exist_ok=True)

    for target_cfg in config.get("targets", []):
        name = target_cfg["name"]
        try:
            interface = rm.load_interface(name)
            pdb_path = Path(interface["pdb_path"])
            if not pdb_path.exists():
                pdb_path = rm.structures_dir / f"{interface['pdb_id']}.pdb"

            fixed_motif = interface.get("fixed_motif")
            fixed_positions = None
            if fixed_motif:
                fixed_positions = find_motif_positions(pdb_path, interface["ligand_chain"], fixed_motif)

            safe_name = name.lower().replace("/", "_").replace(" ", "_")
            fasta = run_proteinmpnn(
                pdb_path=pdb_path,
                output_dir=mpnn_output_dir / safe_name,
                fixed_chains=interface["receptor_chain"],
                design_chains=interface["ligand_chain"],
                num_sequences=config["candidates"]["initial_pool"],
                temperatures=interface.get("mpnn_temperatures", [0.1, 0.2, 0.3]),
                fixed_positions=fixed_positions,
            )

            candidates = parse_mpnn_fasta(fasta, name, interface["design_chain_length"])
            all_candidates.extend(candidates)
            log.info(f"  {name}: {len(candidates)} candidates designed")
        except Exception as e:
            log.error(f"{name} design failed: {e}")

    # Save designs
    output_path = rm.candidates_dir / "mpnn_designs.json"
    with open(output_path, "w") as f:
        json.dump(all_candidates, f, indent=2)

    log.info(f"\nDesign complete: {len(all_candidates)} candidates → {output_path}")
    return all_candidates


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    args = parser.parse_args()
    run(run_dir=args.run_dir)
