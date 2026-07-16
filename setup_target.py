#!/usr/bin/env python3
"""
PeptideScreen — Interactive Config Setup
=========================================
CLI walkthrough that generates config.yaml with all pipeline settings.
Every field has a sensible default — press Enter to accept.

Uses the RCSB REST API to look up chain info — does NOT download PDB files.
The actual download happens in fetch_structures.py.

Usage:
    python3 setup_target.py

Outputs:
    config.yaml
"""

import sys
print("Loading...", flush=True)
from pathlib import Path

import requests
import yaml
print("Ready.\n", flush=True)

ROOT = Path(__file__).resolve().parent
RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"


def prompt(text: str, default=None) -> str:
    """Prompt user with optional default. Returns stripped input or default."""
    if default is not None:
        display = f"{text} [{default}]: "
    else:
        display = f"{text}: "
    sys.stdout.flush()
    val = input(display).strip()
    if not val and default is not None:
        return str(default)
    return val


def prompt_list(text: str, default: list = None) -> list:
    """Prompt for comma-separated list. Returns list of strings."""
    default_str = ",".join(str(x) for x in default) if default else ""
    raw = prompt(text, default_str if default else None)
    if not raw:
        return []
    return [x.strip() for x in raw.split(",") if x.strip()]


def prompt_int(text: str, default: int) -> int:
    val = prompt(text, default)
    try:
        return int(val)
    except ValueError:
        print(f"  Invalid number, using default: {default}")
        return default


def prompt_float(text: str, default: float) -> float:
    val = prompt(text, default)
    try:
        return float(val)
    except ValueError:
        print(f"  Invalid number, using default: {default}")
        return default


def lookup_pdb_chains(pdb_id: str) -> dict:
    """
    Query RCSB REST API for chain info. No file download needed.

    Returns:
        Dict mapping chain_id -> {residue_count, sequence, description}
    """
    pdb_id = pdb_id.upper()

    # Get entry metadata
    print(f"  Looking up {pdb_id} on RCSB...", end=" ", flush=True)
    resp = requests.get(RCSB_ENTRY_URL.format(pdb_id=pdb_id), timeout=15)
    if resp.status_code != 200:
        print("FAILED")
        print(f"  PDB ID '{pdb_id}' not found. Check at https://www.rcsb.org/structure/{pdb_id}")
        sys.exit(1)

    entry = resp.json()
    entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])

    if not entity_ids:
        print("FAILED")
        print(f"  No polymer entities found in {pdb_id}.")
        sys.exit(1)

    # Get each entity's chain info and sequence
    chains = {}
    for eid in entity_ids:
        eresp = requests.get(RCSB_ENTITY_URL.format(pdb_id=pdb_id, entity_id=eid), timeout=15)
        if eresp.status_code != 200:
            continue
        edata = eresp.json()

        # Get chain IDs mapped to this entity
        chain_ids = (
            edata.get("rcsb_polymer_entity_container_identifiers", {})
            .get("auth_asym_ids", [])
        )
        seq = edata.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can", "")
        description = edata.get("rcsb_polymer_entity", {}).get("pdbx_description", "")
        seq_len = len(seq) if seq else 0

        for cid in chain_ids:
            chains[cid] = {
                "residue_count": seq_len,
                "sequence": seq,
                "description": description,
            }

    print("done.", flush=True)
    return chains


def setup_target() -> dict:
    """Walk through target setup. Returns target config dict."""
    print("\n── Target Setup ──\n", flush=True)

    name = prompt("Target name (e.g., Nrf2/KEAP1)")
    if not name:
        print("  Target name is required.")
        sys.exit(1)

    pdb_id = prompt("PDB ID (e.g., 2FLU)")
    if not pdb_id:
        print("  PDB ID is required.")
        sys.exit(1)
    pdb_id = pdb_id.upper()

    # Look up chain info from RCSB API (no download)
    chains = lookup_pdb_chains(pdb_id)

    if not chains:
        print("  No protein chains found in this PDB.")
        sys.exit(1)

    print(f"\n  Chain summary for {pdb_id}:", flush=True)
    chain_ids = list(chains.keys())
    for cid, info in chains.items():
        desc = info["description"][:40] if info["description"] else ""
        preview = info["sequence"][:25]
        suffix = "..." if len(info["sequence"]) > 25 else ""
        print(f"    Chain {cid}: {info['residue_count']} aa | {desc} | {preview}{suffix}")

    # Chain selection
    print(flush=True)
    receptor_chain = prompt(f"Receptor chain [{'/'.join(chain_ids)}]")
    if receptor_chain not in chains:
        print(f"  Chain '{receptor_chain}' not found. Available: {chain_ids}")
        sys.exit(1)

    ligand_chain = prompt(f"Ligand chain [{'/'.join(chain_ids)}]")
    if ligand_chain not in chains:
        print(f"  Chain '{ligand_chain}' not found. Available: {chain_ids}")
        sys.exit(1)

    if receptor_chain == ligand_chain:
        print("  Receptor and ligand chains must be different.")
        sys.exit(1)

    ligand_seq = chains[ligand_chain]["sequence"]
    ligand_len = chains[ligand_chain]["residue_count"]
    print(f"\n  Ligand sequence: {ligand_seq} ({ligand_len} aa)", flush=True)

    # Optional settings
    print("\n── Optional Settings (press Enter to use defaults) ──\n", flush=True)

    # Active residues
    raw_residues = prompt("Active residues (comma-separated, Enter to auto-detect later)")
    if raw_residues:
        try:
            active_residues = [int(x.strip()) for x in raw_residues.split(",")]
        except ValueError:
            print("  Invalid residue numbers. Will auto-detect in analyze_interface.py.")
            active_residues = None
    else:
        active_residues = None
        print("  Will auto-detect in analyze_interface.py.", flush=True)

    # Fixed motif
    fixed_motif = prompt("Fixed motif to lock during design (Enter to skip)")
    if fixed_motif:
        if fixed_motif.upper() in ligand_seq.upper():
            idx = ligand_seq.upper().index(fixed_motif.upper())
            print(f"  Found '{fixed_motif}' in ligand at position {idx + 1}-{idx + len(fixed_motif)}")
        else:
            print(f"  WARNING: '{fixed_motif}' not found in ligand sequence.")
            print(f"  Motif will be ignored — all positions will be redesigned.")
            fixed_motif = None

    # Design chain length
    design_chain_length = prompt_int("Design chain length", ligand_len)

    # MPNN temperatures
    raw_temps = prompt_list("MPNN temperatures (comma-separated)", [0.1, 0.2, 0.3])
    try:
        temperatures = [float(t) for t in raw_temps]
    except ValueError:
        print("  Invalid temperatures, using defaults.")
        temperatures = [0.1, 0.2, 0.3]

    # Seed sequences
    seeds = prompt_list("Seed sequences (comma-separated, Enter to skip)")

    # Build target dict
    target = {
        "name": name,
        "pdb_id": pdb_id,
        "receptor_chain": receptor_chain,
        "ligand_chain": ligand_chain,
        "design_chain_length": design_chain_length,
        "mpnn_temperatures": temperatures,
    }

    if active_residues:
        target["active_residues"] = active_residues
    if fixed_motif:
        target["fixed_motif"] = fixed_motif
    if seeds:
        target["seed_sequences"] = seeds

    # Summary
    print(f"\n── Target Summary ──", flush=True)
    print(f"  Name:             {name}")
    print(f"  PDB:              {pdb_id}")
    print(f"  Receptor:         chain {receptor_chain} ({chains[receptor_chain]['residue_count']} aa)")
    print(f"  Ligand:           chain {ligand_chain} ({ligand_len} aa)")
    print(f"  Active residues:  {len(active_residues) if active_residues else 'auto-detect'}")
    print(f"  Fixed motif:      {fixed_motif or 'none'}")
    print(f"  Design length:    {design_chain_length}")
    print(f"  Temperatures:     {temperatures}")
    print(f"  Seed sequences:   {len(seeds)}")

    return target


def setup_docking() -> dict:
    """Walk through docking settings."""
    print("\n── Docking Settings ──\n", flush=True)

    vina_exhaust = prompt_int("Vina exhaustiveness", 8)
    vina_modes = prompt_int("Vina num modes", 9)
    haddock_max = prompt_int("HADDOCK3 max candidates", 20)
    rigidbody = prompt_int("HADDOCK3 rigidbody sampling", 1000)
    flexref = prompt_int("HADDOCK3 flexref sampling", 200)
    emref = prompt_int("HADDOCK3 emref sampling", 200)

    return {
        "methods": ["vina", "haddock3"],
        "vina": {
            "exhaustiveness": vina_exhaust,
            "num_modes": vina_modes,
        },
        "haddock3": {
            "max_candidates": haddock_max,
            "sampling": {
                "rigidbody": rigidbody,
                "flexref": flexref,
                "emref": emref,
            },
        },
    }


def setup_md() -> dict:
    """Walk through MD stability settings."""
    print("\n── MD Stability Settings ──\n", flush=True)

    force_field = prompt("Force field", "charmm36")
    temp = prompt_float("Temperature (K)", 310.15)
    ionic = prompt_float("Ionic strength (M)", 0.15)
    ns = prompt_float("Simulation length (ns)", 10.0)
    checkpoint = prompt_float("Checkpoint interval (ps)", 500)
    report = prompt_float("Report interval (ps)", 50)

    return {
        "force_field": force_field,
        "temperature_K": temp,
        "ionic_strength_molar": ionic,
        "default_ns": ns,
        "checkpoint_interval_ps": checkpoint,
        "report_interval_ps": report,
    }


def setup_candidates() -> dict:
    """Walk through candidate pipeline settings."""
    print("\n── Candidate Pipeline Settings ──\n", flush=True)

    pool = prompt_int("Initial pool size", 50)
    shortlist = prompt_int("Final shortlist size", 10)
    min_len = prompt_int("Min peptide length (aa)", 5)
    max_len = prompt_int("Max peptide length (aa)", 30)
    max_mw = prompt_int("Max molecular weight (Da)", 2500)

    return {
        "initial_pool": pool,
        "final_shortlist": shortlist,
        "length": {
            "min": min_len,
            "max": max_len,
        },
        "filters": {
            "max_molecular_weight": max_mw,
            "flag_aggregation_prone": True,
        },
    }


def main():
    print("\n=== PeptideScreen — Config Setup ===", flush=True)

    # Collect targets
    targets = []
    while True:
        target = setup_target()
        targets.append(target)

        another = prompt("\nAdd another target? [y/N]", "N")
        if another.lower() != "y":
            break

    # Pipeline settings
    docking = setup_docking()
    md = setup_md()
    candidates = setup_candidates()

    # Build full config
    config = {
        "project": {
            "name": "PeptideScreen",
            "version": "1.0",
            "description": "Computational peptide discovery pipeline.",
        },
        "targets": targets,
        "docking": docking,
        "md_stability": md,
        "candidates": candidates,
        "outputs": {
            "structures": "data/structures/",
            "candidates": "data/candidates/",
            "docking": "data/docking/",
            "processed": "data/processed/",
            "reports": "outputs/reports/",
            "figures": "outputs/figures/",
        },
    }

    # Save
    config_path = ROOT / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)

    print(f"\nConfig saved → {config_path}")
    print("\nNext steps:")
    print("  python3 modules/01_targets/fetch_structures.py")
    print("  python3 modules/01_targets/analyze_interface.py")


if __name__ == "__main__":
    main()
