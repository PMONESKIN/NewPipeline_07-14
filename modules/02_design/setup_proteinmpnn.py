"""
Module 02 — Step 1: ProteinMPNN Setup
======================================
One-time setup: clones the official ProteinMPNN repository and verifies
the environment is ready for peptide design.

ProteinMPNN (Dauparas et al. 2022) is not on PyPI — it runs as a set of
Python scripts. This script clones it into tools/ProteinMPNN/ and checks
that all dependencies are available.

Usage:
    python3 modules/02_design/setup_proteinmpnn.py

Run once before design_peptides.py.
"""

import logging
import subprocess
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
PROTEINMPNN_REPO = "https://github.com/dauparas/ProteinMPNN.git"
MPNN_DIR = ROOT / "tools" / "ProteinMPNN"


def clone_proteinmpnn() -> None:
    """Clone ProteinMPNN into tools/ProteinMPNN/ if not already present."""
    MPNN_DIR.parent.mkdir(parents=True, exist_ok=True)

    if MPNN_DIR.exists():
        log.info(f"ProteinMPNN already cloned at {MPNN_DIR}")
        return

    log.info(f"Cloning ProteinMPNN from {PROTEINMPNN_REPO}...")
    result = subprocess.run(
        ["git", "clone", PROTEINMPNN_REPO, str(MPNN_DIR)],
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        raise RuntimeError(f"git clone failed:\n{result.stderr}")

    log.info(f"ProteinMPNN cloned to {MPNN_DIR}")


def verify_environment() -> None:
    """Check all required packages are importable."""
    required = {"torch": "torch", "numpy": "numpy", "Bio": "biopython"}
    missing = []

    for pkg, pip_name in required.items():
        try:
            __import__(pkg)
            log.info(f"  [OK] {pkg}")
        except ImportError:
            log.error(f"  [MISSING] {pkg} — pip install {pip_name}")
            missing.append(pip_name)

    if missing:
        raise RuntimeError(
            f"Missing packages. Run: pip install {' '.join(missing)}"
        )


def check_scripts() -> None:
    """Verify key ProteinMPNN scripts are present."""
    required = [
        MPNN_DIR / "protein_mpnn_run.py",
        MPNN_DIR / "helper_scripts" / "parse_multiple_chains.py",
        MPNN_DIR / "helper_scripts" / "assign_fixed_chains.py",
    ]

    for script in required:
        if script.exists():
            log.info(f"  [OK] {script.relative_to(MPNN_DIR)}")
        else:
            raise FileNotFoundError(
                f"Expected script not found: {script}. "
                "Try deleting tools/ProteinMPNN/ and re-running this setup."
            )


def check_device() -> None:
    """Report available compute device."""
    import torch

    if torch.cuda.is_available():
        device = f"CUDA ({torch.cuda.get_device_name(0)})"
    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        device = "MPS (Apple Silicon) — ProteinMPNN will use CPU instead"
    else:
        device = "CPU (no GPU — design will be slower but works fine)"

    log.info(f"  Compute device: {device}")


def run() -> None:
    log.info("=== ProteinMPNN Setup ===\n")

    log.info("Step 1: Checking Python environment...")
    verify_environment()

    log.info("\nStep 2: Cloning ProteinMPNN...")
    clone_proteinmpnn()

    log.info("\nStep 3: Verifying scripts...")
    check_scripts()

    log.info("\nStep 4: Checking compute device...")
    check_device()

    log.info("\nSetup complete. Ready to run design_peptides.py")


if __name__ == "__main__":
    run()
