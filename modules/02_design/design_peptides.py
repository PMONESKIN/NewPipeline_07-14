"""
Module 02 — Step 2: ProteinMPNN Peptide Design
===============================================
Runs ProteinMPNN to generate novel peptide sequences for each target,
using co-crystal PDB structures as design templates.

Reads target parameters from the interface JSON files produced by
Module 01 (analyze_interface.py), including:
  - PDB path, chain assignments
  - Fixed motif positions
  - Design chain length
  - MPNN temperatures

Usage:
    python3 modules/02_design/design_peptides.py

Run after setup_proteinmpnn.py and Module 01.
Outputs: data/candidates/mpnn_designs.json
"""

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


def load_config(config_path: str = None) -> dict:
    path = Path(config_path) if config_path else ROOT / "config.yaml"
    with open(path) as f:
        return yaml.safe_load(f)


def load_interface(target_name: str, processed_dir: Path) -> dict:
    """Load the interface JSON produced by analyze_interface.py."""
    safe_name = target_name.lower().replace("/", "_").replace(" ", "_")
    path = processed_dir / f"{safe_name}_interface.json"
    if not path.exists():
        raise FileNotFoundError(
            f"Interface file not found: {path}. Run analyze_interface.py first."
        )
    with open(path) as f:
        return json.load(f)


def get_device_flag() -> str:
    """Return the appropriate device flag for ProteinMPNN."""
    try:
        import torch
        if torch.cuda.is_available():
            return "--use_gpu"
    except ImportError:
        pass
    return ""


def find_motif_positions(pdb_path: Path, chain_id: str, motif: str) -> list[int]:
    """
    Find 1-indexed sequential positions of a motif within a PDB chain.
    Returns positions in ProteinMPNN's format: 1-indexed, sequential.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))

    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                seq = ""
                for res in chain:
                    if is_aa(res, standard=True):
                        aa = protein_letters_3to1.get(res.get_resname(), "X")
                        seq += aa

                idx = seq.find(motif)
                if idx >= 0:
                    positions = list(range(idx + 1, idx + 1 + len(motif)))
                    log.info(
                        f"  Motif '{motif}' found at positions {positions} "
                        f"in chain {chain_id}"
                    )
                    return positions
                else:
                    log.warning(
                        f"  Motif '{motif}' NOT found in chain {chain_id} "
                        f"({len(seq)} aa). No positions locked."
                    )
    return []


def run_proteinmpnn(
    pdb_path: Path,
    output_dir: Path,
    fixed_chains: str,
    design_chains: str,
    num_sequences: int = 30,
    temperatures: list[float] = None,
    seed: int = 42,
    fixed_positions: list[int] | None = None,
) -> Path:
    """
    Run ProteinMPNN on a PDB file.

    Returns:
        Path to output FASTA file with designed sequences
    """
    if not MPNN_DIR.exists():
        raise FileNotFoundError(
            f"ProteinMPNN not found at {MPNN_DIR}. "
            "Run: python3 modules/02_design/setup_proteinmpnn.py"
        )

    temperatures = temperatures or [0.1, 0.2, 0.3]
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Parse PDB chains into MPNN's JSON format
    parsed_dir = output_dir / "parsed"
    parsed_dir.mkdir(exist_ok=True)
    parsed_jsonl = parsed_dir / "parsed.jsonl"

    with tempfile.TemporaryDirectory() as tmp_pdb_dir:
        import shutil
        shutil.copy(pdb_path, Path(tmp_pdb_dir) / pdb_path.name)

        parse_result = subprocess.run(
            [
                sys.executable,
                str(MPNN_DIR / "helper_scripts" / "parse_multiple_chains.py"),
                "--input_path", tmp_pdb_dir,
                "--output_path", str(parsed_jsonl),
            ],
            capture_output=True, text=True,
        )

    if parse_result.returncode != 0:
        raise RuntimeError(f"MPNN chain parsing failed:\n{parse_result.stderr}")
    log.info("  Chains parsed")

    # Step 2: Assign fixed vs designable chains
    assigned_jsonl = parsed_dir / "assigned.jsonl"
    assign_result = subprocess.run(
        [
            sys.executable,
            str(MPNN_DIR / "helper_scripts" / "assign_fixed_chains.py"),
            "--input_path", str(parsed_jsonl),
            "--output_path", str(assigned_jsonl),
            "--chain_list", design_chains,
        ],
        capture_output=True, text=True,
    )

    if assign_result.returncode != 0:
        raise RuntimeError(f"MPNN chain assignment failed:\n{assign_result.stderr}")
    log.info(f"  Chains assigned (fixed: {fixed_chains}, design: {design_chains})")

    # Step 3: Lock specific positions if fixed_motif was specified
    fixed_pos_jsonl = None
    if fixed_positions:
        fixed_pos_jsonl = parsed_dir / "fixed_positions.jsonl"
        pdb_name = pdb_path.stem

        pos_dict = {pdb_name: {}}
        pos_dict[pdb_name][design_chains] = fixed_positions
        for fc in fixed_chains.split(","):
            fc = fc.strip()
            if fc:
                pos_dict[pdb_name][fc] = []

        with open(fixed_pos_jsonl, "w") as f:
            f.write(json.dumps(pos_dict) + "\n")

        log.info(
            f"  Fixed positions: {len(fixed_positions)} residues locked "
            f"in chain {design_chains}"
        )

    # Step 4: Run ProteinMPNN
    temp_str = " ".join(str(t) for t in temperatures)
    device_flag = get_device_flag()

    mpnn_cmd = [
        sys.executable,
        str(MPNN_DIR / "protein_mpnn_run.py"),
        "--jsonl_path", str(parsed_jsonl),
        "--chain_id_jsonl", str(assigned_jsonl),
        "--out_folder", str(output_dir),
        "--num_seq_per_target", str(num_sequences),
        "--sampling_temp", temp_str,
        "--seed", str(seed),
        "--batch_size", "1",
    ]

    if fixed_pos_jsonl:
        mpnn_cmd.extend(["--fixed_positions_jsonl", str(fixed_pos_jsonl)])

    if device_flag:
        mpnn_cmd.append(device_flag)

    log.info(f"  Running ProteinMPNN ({num_sequences} seqs x {len(temperatures)} temps)...")
    mpnn_result = subprocess.run(mpnn_cmd, capture_output=True, text=True)

    if mpnn_result.returncode != 0:
        raise RuntimeError(f"ProteinMPNN failed:\n{mpnn_result.stderr}")

    fasta_files = list((output_dir / "seqs").glob("*.fa"))
    if not fasta_files:
        raise FileNotFoundError(
            f"No FASTA output in {output_dir}/seqs/. Check ProteinMPNN output."
        )

    log.info(f"  Design complete: {fasta_files[0]}")
    return fasta_files[0]


def parse_mpnn_fasta(
    fasta_path: Path,
    target_name: str,
    design_chain_len: int | None = None,
) -> list[dict]:
    """
    Parse ProteinMPNN output FASTA into candidate dicts.

    MPNN outputs full multi-chain sequences. The designed chain portion
    is placed first. We extract only that portion.
    """
    entries = []

    with open(fasta_path) as f:
        lines = f.read().strip().split("\n")

    i = 0
    while i < len(lines):
        if lines[i].startswith(">"):
            header = lines[i][1:]
            seq = lines[i + 1].strip() if i + 1 < len(lines) else ""
            entries.append((header, seq))
            i += 2
        else:
            i += 1

    if not entries:
        return []

    # First entry is the reference
    _, ref_seq = entries[0]
    if design_chain_len is None:
        design_chain_len = len(ref_seq)

    candidates = []
    for header, seq in entries[1:]:  # Skip reference
        score = None
        temp = None
        seq_recovery = None

        for part in header.split(","):
            part = part.strip()
            if part.startswith("score="):
                try:
                    score = float(part.split("=")[1])
                except ValueError:
                    pass
            elif part.startswith("T="):
                try:
                    temp = float(part.split("=")[1])
                except ValueError:
                    pass
            elif part.startswith("seq_recovery="):
                try:
                    seq_recovery = float(part.split("=")[1])
                except ValueError:
                    pass

        designed_seq = seq[:design_chain_len].upper()

        candidates.append({
            "sequence": designed_seq,
            "length": len(designed_seq),
            "mpnn_score": score,
            "sampling_temperature": temp,
            "sequence_recovery": seq_recovery,
            "design_source": "proteinmpnn",
            "target_name": target_name,
        })

    return candidates


def design_target(interface: dict, config: dict, output_dir: Path) -> list[dict]:
    """
    Run ProteinMPNN design for a single target using its interface JSON.
    """
    target_name = interface["target_name"]
    pdb_path = Path(interface["pdb_path"])
    receptor_chain = interface["receptor_chain"]
    ligand_chain = interface["ligand_chain"]
    design_chain_length = interface["design_chain_length"]
    temperatures = interface.get("mpnn_temperatures", [0.1, 0.2, 0.3])
    fixed_motif = interface.get("fixed_motif")

    log.info(f"\n--- Designing peptides for {target_name} ---")
    log.info(
        f"  PDB: {pdb_path.name} | receptor: chain {receptor_chain} | "
        f"design: chain {ligand_chain} ({design_chain_length} aa)"
    )

    if not pdb_path.exists():
        raise FileNotFoundError(
            f"{pdb_path} not found. Run fetch_structures.py first."
        )

    # Find positions to lock if fixed_motif is specified
    fixed_positions = None
    if fixed_motif:
        fixed_positions = find_motif_positions(pdb_path, ligand_chain, fixed_motif)
        if not fixed_positions:
            log.warning(f"  fixed_motif '{fixed_motif}' not found — all positions designable")

    # Output directory for this target
    safe_name = target_name.lower().replace("/", "_").replace(" ", "_").replace("-", "_")
    mpnn_output = output_dir / safe_name

    fasta = run_proteinmpnn(
        pdb_path=pdb_path,
        output_dir=mpnn_output,
        fixed_chains=receptor_chain,
        design_chains=ligand_chain,
        num_sequences=config["candidates"]["initial_pool"],
        temperatures=temperatures,
        fixed_positions=fixed_positions,
    )

    candidates = parse_mpnn_fasta(fasta, target_name, design_chain_len=design_chain_length)
    log.info(f"  Generated {len(candidates)} candidates")
    return candidates


def run(config_path: str = None) -> list[dict]:
    """
    Run ProteinMPNN design for all targets.
    Reads interface JSONs from Module 01 output.
    """
    config = load_config(config_path)
    processed_dir = Path(config["outputs"]["processed"])
    candidates_dir = Path(config["outputs"]["candidates"])
    mpnn_output_dir = candidates_dir / "mpnn_output"
    mpnn_output_dir.mkdir(parents=True, exist_ok=True)

    targets = config.get("targets", [])
    if not targets:
        log.error("No targets found in config.yaml.")
        return []

    log.info(f"Designing peptides for {len(targets)} target(s)...")
    all_candidates = []

    for target_cfg in targets:
        target_name = target_cfg.get("name", "UNKNOWN")
        try:
            interface = load_interface(target_name, processed_dir)
            candidates = design_target(interface, config, mpnn_output_dir)
            all_candidates.extend(candidates)
        except Exception as e:
            log.error(f"{target_name} design failed: {e}")
            log.error("  Continuing with remaining targets...")

    # Save all candidates
    output_path = candidates_dir / "mpnn_designs.json"
    with open(output_path, "w") as f:
        json.dump(all_candidates, f, indent=2)

    log.info(f"\nDesign complete: {len(all_candidates)} total candidates -> {output_path}")

    for target_cfg in targets:
        name = target_cfg["name"]
        count = sum(1 for c in all_candidates if c["target_name"] == name)
        log.info(f"  {name}: {count} candidates")

    return all_candidates


if __name__ == "__main__":
    run()
