#!/bin/bash
# PeptideScreen — Linux/AWS Setup Script
# Run once on a fresh clone.
#
# Prerequisites:
#   - Python 3.11+
#   - Git
#   - (Optional) CUDA GPU for ProteinMPNN + OpenMM
#
# Usage: bash setup.sh

set -e
echo ""
echo "=== PeptideScreen Setup ==="
echo ""

# Check Python
if ! command -v python3 &> /dev/null; then
    echo "ERROR: python3 not found."
    echo "  Ubuntu/Debian: sudo apt install python3 python3-venv python3-pip"
    echo "  Amazon Linux:  sudo yum install python3"
    exit 1
fi
echo "[OK] $(python3 --version)"

# Virtual environment
echo ""
echo "[1/4] Creating virtual environment..."
if [ -d ".venv" ]; then
    echo "      .venv already exists, skipping."
else
    python3 -m venv .venv
fi
source .venv/bin/activate
echo "      Activated."

# Install core dependencies
echo ""
echo "[2/4] Installing core dependencies..."
pip install --upgrade pip --quiet
pip install -r requirements.txt --quiet
echo "      Done."

# ProteinMPNN
echo ""
echo "[3/4] Setting up ProteinMPNN..."
if [ -f "tools/ProteinMPNN/protein_mpnn_run.py" ]; then
    echo "      Already present, skipping."
else
    python3 modules/02_design/setup_proteinmpnn.py
fi

# AutoDock Vina
echo ""
echo "[4/4] Checking AutoDock Vina..."
if command -v vina &> /dev/null; then
    echo "      [OK] $(vina --version 2>&1 | head -1)"
elif [ -f "tools/vina" ]; then
    echo "      [OK] Found at tools/vina"
else
    echo "      Vina not found — downloading Linux binary..."
    mkdir -p tools
    ARCH=$(uname -m)
    if [ "$ARCH" = "x86_64" ]; then
        curl -L "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_linux_x86_64" -o tools/vina && chmod +x tools/vina
        echo "      [OK] $(tools/vina --version 2>&1 | head -1)"
    else
        echo "      Unknown architecture ($ARCH). Download manually:"
        echo "      https://github.com/ccsb-scripps/AutoDock-Vina/releases"
    fi
fi

echo ""
echo "=== Core setup complete ==="
echo ""
echo "Optional components (install separately):"
echo ""
echo "  HADDOCK3 (docking validation):"
echo "    pip install haddock3"
echo "    # Requires CNS solver: http://cns-online.org (academic license)"
echo "    export CNS_SOLVE=/path/to/cns_solve"
echo ""
echo "  OpenMM (MD stability simulations):"
echo "    conda install -c conda-forge openmm pdbfixer cudatoolkit"
echo ""
echo "Activate environment in future sessions:"
echo "  source .venv/bin/activate"
echo ""
echo "Next step: Edit config.yaml with your target, then run:"
echo "  python3 modules/01_targets/fetch_structures.py"
