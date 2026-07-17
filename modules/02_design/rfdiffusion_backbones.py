"""
Module 02 — RFdiffusion Backbone Generation
=============================================
Generates diverse peptide backbone geometries that fit the target
binding pocket using RFdiffusion. Each backbone is then fed to
ProteinMPNN for sequence design.

This is OPTIONAL — only used when you want to explore backbone
diversity beyond the co-crystal structure.

Usage:
    python3 modules/02_design/rfdiffusion_backbones.py --run-dir path --n-backbones 100 --peptide-length 8

Requirements:
    pip install rfdiffusion (or clone from GitHub)
    RFdiffusion model weights (~1GB download on first run)
"""

import argparse
import json
import logging
import os
import subprocess
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]


def check_rfdiffusion():
    """Check if RFdiffusion is installed."""
    try:
        result = subprocess.run(
            ["python", "-c", "import rfdiffusion; print('ok')"],
            capture_output=True, text=True, timeout=30
        )
        return result.returncode == 0
    except Exception:
        return False


def setup_rfdiffusion():
    """Install RFdiffusion if not present."""
    log.info("Installing RFdiffusion...")
    cmds = [
        "pip install -q dgl -f https://data.dgl.ai/wheels/torch-2.1/cu121/repo.html",
        "pip install -q e3nn",
        "pip install -q git+https://github.com/RosettaCommons/RFdiffusion.git",
    ]
    for cmd in cmds:
        result = subprocess.run(cmd.split(), capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            log.error(f"  Install failed: {result.stderr[:200]}")
            return False

    # Download model weights
    log.info("Downloading RFdiffusion model weights...")
    weights_dir = ROOT / "tools" / "rfdiffusion_weights"
    weights_dir.mkdir(parents=True, exist_ok=True)

    weight_urls = [
        "http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt",
        "http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt",
    ]

    for url in weight_urls:
        fname = url.split("/")[-1]
        target = weights_dir / fname
        if not target.exists():
            subprocess.run(["wget", "-q", url, "-O", str(target)], timeout=300)

    return True


def generate_backbones(
    receptor_pdb: Path,
    hotspot_residues: list[int],
    receptor_chain: str,
    peptide_length: int,
    n_backbones: int,
    output_dir: Path,
) -> list[Path]:
    """
    Generate diverse peptide backbones using RFdiffusion.

    Args:
        receptor_pdb: Path to receptor PDB
        hotspot_residues: Receptor residues to target (binding pocket)
        receptor_chain: Chain ID of receptor
        peptide_length: Length of peptide to generate
        n_backbones: Number of backbones to generate
        output_dir: Where to save backbone PDBs

    Returns:
        List of paths to generated backbone PDBs
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Format hotspot residues for RFdiffusion
    hotspots = ",".join(f"{receptor_chain}{r}" for r in hotspot_residues[:10])

    # RFdiffusion config
    cmd = [
        sys.executable, "-m", "rfdiffusion.inference.model_runner",
        f"inference.output_prefix={output_dir}/backbone",
        f"inference.input_pdb={receptor_pdb}",
        f"inference.num_designs={n_backbones}",
        f"contigmap.contigs=[{receptor_chain}1-1000/0 {peptide_length}-{peptide_length}]",
        f"ppi.hotspot_res=[{hotspots}]",
        "inference.ckpt_override_path=tools/rfdiffusion_weights/Complex_base_ckpt.pt",
        "denoiser.noise_scale_ca=0",
        "denoiser.noise_scale_frame=0",
    ]

    log.info(f"  Generating {n_backbones} backbones ({peptide_length} aa)...")
    log.info(f"  Hotspots: {hotspots}")

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

    if result.returncode != 0:
        log.error(f"  RFdiffusion failed: {result.stderr[:300]}")
        return []

    # Collect output PDBs
    backbone_pdbs = sorted(output_dir.glob("backbone_*.pdb"))
    log.info(f"  Generated {len(backbone_pdbs)} backbones")

    return backbone_pdbs


def filter_backbones_by_pocket(
    backbone_pdbs: list[Path],
    receptor_pdb: Path,
    hotspot_residues: list[int],
    receptor_chain: str,
    contact_cutoff: float = 8.0,
    min_contacts: int = 3,
) -> list[Path]:
    """
    Filter backbones: keep only those that land in the binding pocket.
    A backbone passes if its peptide chain makes at least min_contacts
    contacts with the hotspot residues.
    """
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import is_aa

    parser = PDBParser(QUIET=True)
    passed = []

    for pdb_path in backbone_pdbs:
        structure = parser.get_structure(pdb_path.stem, str(pdb_path))

        # Find peptide and receptor chains
        chains = []
        for model in structure:
            for chain in model:
                residues = [r for r in chain if is_aa(r, standard=True)]
                chains.append((chain.id, len(residues), residues))

        if len(chains) < 2:
            continue

        # Sort by length — receptor is longer, peptide is shorter
        chains.sort(key=lambda x: x[1], reverse=True)
        receptor_residues = chains[0][2]
        peptide_residues = chains[-1][2]

        # Count contacts between peptide and hotspot residues
        hotspot_set = set(hotspot_residues)
        contacts = 0

        for pep_res in peptide_residues:
            pep_ca = pep_res["CA"] if "CA" in pep_res else None
            if pep_ca is None:
                continue

            for rec_res in receptor_residues:
                if rec_res.get_id()[1] not in hotspot_set:
                    continue
                rec_ca = rec_res["CA"] if "CA" in rec_res else None
                if rec_ca is None:
                    continue

                dist = pep_ca - rec_ca
                if dist < contact_cutoff:
                    contacts += 1
                    break

        if contacts >= min_contacts:
            passed.append(pdb_path)

    log.info(f"  Pocket filter: {len(passed)}/{len(backbone_pdbs)} backbones in binding pocket")
    return passed


def run(run_dir: str, n_backbones: int = 100, peptide_length: int = 8) -> list[Path]:
    """Generate and filter RFdiffusion backbones."""
    sys.path.insert(0, str(ROOT))
    from modules.run_manager import RunManager

    rm = RunManager(run_dir=run_dir)
    config = rm.config

    # Check/install RFdiffusion
    if not check_rfdiffusion():
        log.info("RFdiffusion not found. Installing...")
        if not setup_rfdiffusion():
            log.error("RFdiffusion installation failed.")
            log.error("Install manually: pip install git+https://github.com/RosettaCommons/RFdiffusion.git")
            return []

    all_backbones = []

    for target_cfg in config.get("targets", []):
        name = target_cfg["name"]
        log.info(f"\n--- Generating backbones for {name} ---")

        try:
            interface = rm.load_interface(name)
            receptor_pdb = rm.structures_dir / f"{interface['pdb_id']}.pdb"

            backbone_dir = rm.processed_dir / "rfdiffusion_backbones"

            backbones = generate_backbones(
                receptor_pdb=receptor_pdb,
                hotspot_residues=interface["active_residues"],
                receptor_chain=interface["receptor_chain"],
                peptide_length=peptide_length,
                n_backbones=n_backbones,
                output_dir=backbone_dir,
            )

            if backbones:
                filtered = filter_backbones_by_pocket(
                    backbones,
                    receptor_pdb,
                    interface["active_residues"],
                    interface["receptor_chain"],
                )
                all_backbones.extend(filtered)

                # Save backbone list
                backbone_info = {
                    "target": name,
                    "total_generated": len(backbones),
                    "passed_filter": len(filtered),
                    "backbone_paths": [str(p) for p in filtered],
                }
                info_path = backbone_dir / "backbone_info.json"
                with open(info_path, "w") as f:
                    json.dump(backbone_info, f, indent=2)

        except Exception as e:
            log.error(f"{name}: backbone generation failed — {e}")

    log.info(f"\nDone. {len(all_backbones)} filtered backbones ready for ProteinMPNN.")
    rm.update_status("backbones_generated")
    return all_backbones


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--n-backbones", type=int, default=100)
    parser.add_argument("--peptide-length", type=int, default=8)
    args = parser.parse_args()
    run(run_dir=args.run_dir, n_backbones=args.n_backbones,
        peptide_length=args.peptide_length)
