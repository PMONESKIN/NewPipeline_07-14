"""
Module 03 — Step 4: Docking Report Generator
=============================================
Generates a comprehensive markdown report with HADDOCK3-ranked candidates,
score tables, property summaries, and methodology notes.
No LLM calls — all deterministic.

Usage:
    python3 modules/03_docking/docking_report.py --run-dir path

Inputs:
    {run_dir}/candidates/candidate_pool.json (with HADDOCK scores + properties)

Outputs:
    {run_dir}/reports/docking_report.md
    {run_dir}/reports/shortlist.json
"""

import argparse
import json
import logging
from datetime import datetime
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]


def fmt(val, decimals=1, suffix=""):
    if val is None:
        return "N/A"
    if isinstance(val, float):
        return f"{val:.{decimals}f}{suffix}"
    return str(val)


def generate_report(candidates: list[dict], config: dict) -> str:
    """Generate the full docking report as markdown."""
    today = datetime.now().strftime("%Y-%m-%d")

    # Determine which score field is available
    score_field = None
    for field in ["haddock_full_score", "haddock_fast_score"]:
        if any(c.get(field) is not None for c in candidates):
            score_field = field
            break

    scored = [c for c in candidates if c.get(score_field) is not None] if score_field else []
    scored.sort(key=lambda c: c[score_field])  # lower = better

    targets = sorted(set(c["target_name"] for c in candidates))

    sections = []

    # Header
    sections.append(f"""# PeptideScreen — Docking Report

**Date:** {today}
**Pipeline:** PeptideScreen v1.0
**Docking method:** HADDOCK3 (local)
**Score field:** {score_field or 'none'}

---

## Executive Summary

| Metric | Value |
|--------|-------|
| Total candidates | {len(candidates)} |
| HADDOCK3 scored | {len(scored)} |
| Scoring mode | {'Full (rigidbody + flexref + emref)' if 'full' in (score_field or '') else 'Fast (rigidbody only)'} |

---
""")

    # Per-target ranked table
    for target in targets:
        target_scored = [c for c in scored if c["target_name"] == target]

        sections.append(f"## {target} — Ranked Candidates\n")

        if not target_scored:
            sections.append(f"*No scored candidates for {target}.*\n")
            continue

        # Summary table
        header = "| Rank | ID | Sequence | Len | HADDOCK | MPNN | Stable? | Protease | Permeability |"
        divider = "|------|-----|----------|-----|---------|------|---------|----------|-------------|"
        sections.append(header)
        sections.append(divider)

        for rank, c in enumerate(target_scored[:20], 1):
            seq = c["sequence"][:18] + ".." if len(c["sequence"]) > 18 else c["sequence"]
            hadd = fmt(c.get(score_field), 1)
            mpnn = fmt(c.get("mpnn_score"), 2)

            # Properties (if computed)
            pc = c.get("physicochemical", {})
            stable = "yes" if pc and not pc.get("is_unstable") else ("NO" if pc else "?")

            pr = c.get("protease_cleavage", {})
            prot_risk = pr.get("protease_risk", "?")

            pm = c.get("permeability", {})
            perm = pm.get("mw_category", "?")[:15]

            sections.append(
                f"| {rank} | {c['id']} | `{seq}` | {c['length']} | {hadd} | {mpnn} | {stable} | {prot_risk} | {perm} |"
            )

        sections.append("")

        # Top 10 detailed
        sections.append(f"### Top 10 — Detailed\n")
        for rank, c in enumerate(target_scored[:10], 1):
            lines = [
                f"#### {rank}. {c['id']}",
                "",
                "| Property | Value |",
                "|----------|-------|",
                f"| Sequence | `{c['sequence']}` |",
                f"| Length | {c['length']} aa |",
                f"| Design source | {c.get('design_source', '?')} |",
                f"| HADDOCK score | {fmt(c.get(score_field), 1)} |",
                f"| MPNN score | {fmt(c.get('mpnn_score'), 2)} |",
            ]

            # Physicochemical
            pc = c.get("physicochemical", {})
            if pc:
                lines.extend([
                    f"| MW | {fmt(pc.get('molecular_weight_da'), 0)} Da |",
                    f"| Net charge (pH 7.4) | {fmt(pc.get('net_charge_ph74'), 1)} |",
                    f"| pI | {fmt(pc.get('isoelectric_point'), 1)} |",
                    f"| GRAVY | {fmt(pc.get('gravy'), 3)} |",
                    f"| Instability index | {fmt(pc.get('instability_index'), 1)} ({'unstable' if pc.get('is_unstable') else 'stable'}) |",
                ])

            # Protease
            pr = c.get("protease_cleavage", {})
            if pr:
                lines.append(f"| Protease risk | {pr.get('protease_risk', '?')} ({pr.get('total_cleavage_sites', 0)} sites) |")

            # Permeability
            pm = c.get("permeability", {})
            if pm:
                lines.append(f"| Permeability | {pm.get('overall_permeability', '?')} |")

            lines.extend(["", "---", ""])
            sections.append("\n".join(lines))

    # Methodology
    haddock_cfg = config.get("docking", {}).get("haddock3", {})
    sampling = haddock_cfg.get("sampling", {})

    sections.append(f"""## Methodology

**HADDOCK3 (local execution):**
- AIR restraints generated from interface analysis (active residues)
- Fast screen: rigidbody docking ({sampling.get('rigidbody_fast', 50)} models)
- Full validation: rigidbody ({sampling.get('rigidbody', 1000)}) + flexref ({sampling.get('flexref', 200)}) + emref ({sampling.get('emref', 200)})
- Scoring from caprieval module (lower = better binding)

**Physicochemical:** Biopython ProtParam (deterministic)

**Protease stability:** Rule-based P1 position matching (trypsin, chymotrypsin, elastase)

**Permeability:** MW-based rules (Bos & Meinardi 500 Da rule) — low confidence

### Caveats

1. HADDOCK scores are relative rankings, not absolute binding affinities.
2. Peptide starting structures use extended conformation (PeptideBuilder). ESMFold integration for folded structures is recommended.
3. Protease predictions are sequence-based only — 3D structure not considered.
4. Permeability estimates are rule-based. Experimental PAMPA/Caco-2 required for validation.
5. All designs are novel sequences — no experimental validation has been performed.

---

*Generated by PeptideScreen Module 03*
""")

    return "\n".join(sections)


def run(run_dir: str) -> None:
    """Generate docking report for a run."""
    import sys
    sys.path.insert(0, str(ROOT))
    from modules.run_manager import RunManager

    rm = RunManager(run_dir=run_dir)
    config = rm.config
    candidates = rm.load_candidates()

    # Check for any HADDOCK scores
    has_scores = any(
        c.get("haddock_fast_score") is not None or c.get("haddock_full_score") is not None
        for c in candidates
    )
    if not has_scores:
        log.error("No HADDOCK scores found. Run haddock3_dock.py first.")
        return

    scored_count = sum(
        1 for c in candidates
        if c.get("haddock_fast_score") is not None or c.get("haddock_full_score") is not None
    )
    log.info(f"Generating report for {scored_count} scored candidates...")

    report_md = generate_report(candidates, config)
    report_path = rm.reports_dir / "docking_report.md"
    report_path.write_text(report_md)
    log.info(f"Report saved -> {report_path}")

    # Save shortlist
    score_field = "haddock_full_score" if any(c.get("haddock_full_score") is not None for c in candidates) else "haddock_fast_score"
    scored = [c for c in candidates if c.get(score_field) is not None]
    scored.sort(key=lambda c: c[score_field])

    shortlist = []
    for target in set(c["target_name"] for c in candidates):
        group = [c for c in scored if c["target_name"] == target]
        shortlist.extend(group[:10])

    shortlist_path = rm.reports_dir / "shortlist.json"
    shortlist_path.write_text(json.dumps(shortlist, indent=2))
    log.info(f"Shortlist ({len(shortlist)} candidates) -> {shortlist_path}")

    rm.update_status("report_generated")
    log.info("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True, help="Run directory")
    args = parser.parse_args()
    run(run_dir=args.run_dir)
