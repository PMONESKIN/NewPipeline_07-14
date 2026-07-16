"""
Run Manager — Creates and manages independent run directories.

Each pipeline execution gets its own timestamped folder under runs/.
ALL modules read/write through this manager — single source of truth.

Usage:
    from modules.run_manager import RunManager

    rm = RunManager()           # creates new run dir
    rm = RunManager("path")     # resume existing run
"""

import json
import shutil
from datetime import datetime
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
RUNS_DIR = ROOT / "runs"


class RunManager:
    def __init__(self, run_dir: str = None, config_path: str = None):
        if run_dir:
            self.run_dir = Path(run_dir)
            if not self.run_dir.is_absolute():
                self.run_dir = ROOT / self.run_dir
            if not self.run_dir.exists():
                raise FileNotFoundError(f"Run directory not found: {self.run_dir}")
        else:
            self.run_dir = self._create_run(config_path)

        self._setup_subdirs()

    def _create_run(self, config_path: str = None) -> Path:
        config_path = Path(config_path) if config_path else ROOT / "config.yaml"
        config = yaml.safe_load(config_path.read_text())

        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M")
        targets = config.get("targets", [])
        target_names = "_".join(
            t.get("name", "unknown").lower().replace("/", "-").replace(" ", "-")
            for t in targets
        )
        if not target_names:
            target_names = "unnamed"

        run_name = f"{timestamp}_{target_names}"
        run_dir = RUNS_DIR / run_name
        run_dir.mkdir(parents=True, exist_ok=True)

        # Copy config snapshot
        shutil.copy(config_path, run_dir / "run_config.yaml")

        # Save run metadata
        run_info = {
            "run_name": run_name,
            "created": datetime.now().isoformat(),
            "targets": [t.get("name") for t in targets],
            "status": "created",
        }
        (run_dir / "run_info.json").write_text(json.dumps(run_info, indent=2))

        print(f"\n=== New run: {run_name} ===")
        print(f"  Directory: {run_dir}")

        return run_dir

    def _setup_subdirs(self):
        self.candidates_dir = self.run_dir / "candidates"
        self.docking_dir = self.run_dir / "docking"
        self.reports_dir = self.run_dir / "reports"
        self.md_dir = self.run_dir / "md_results"
        self.processed_dir = self.run_dir / "processed"
        self.structures_dir = self.run_dir / "structures"
        self.folded_dir = self.run_dir / "candidates" / "folded_structures"

        for d in [self.candidates_dir, self.docking_dir, self.reports_dir,
                  self.md_dir, self.processed_dir, self.structures_dir,
                  self.folded_dir]:
            d.mkdir(exist_ok=True)

    @property
    def config(self) -> dict:
        return yaml.safe_load((self.run_dir / "run_config.yaml").read_text())

    @property
    def candidate_pool_path(self) -> Path:
        return self.candidates_dir / "candidate_pool.json"

    def load_candidates(self) -> list[dict]:
        if not self.candidate_pool_path.exists():
            return []
        with open(self.candidate_pool_path) as f:
            return json.load(f)

    def save_candidates(self, candidates: list[dict]):
        with open(self.candidate_pool_path, "w") as f:
            json.dump(candidates, f, indent=2)

    def load_interface(self, target_name: str) -> dict:
        safe_name = target_name.lower().replace("/", "_").replace(" ", "_")
        path = self.processed_dir / f"{safe_name}_interface.json"
        if not path.exists():
            raise FileNotFoundError(f"Interface file not found: {path}")
        with open(path) as f:
            return json.load(f)

    def save_interface(self, target_name: str, data: dict):
        safe_name = target_name.lower().replace("/", "_").replace(" ", "_")
        path = self.processed_dir / f"{safe_name}_interface.json"
        with open(path, "w") as f:
            json.dump(data, f, indent=2)

    def update_status(self, status: str):
        info_path = self.run_dir / "run_info.json"
        info = json.loads(info_path.read_text())
        info["status"] = status
        info["last_updated"] = datetime.now().isoformat()
        info_path.write_text(json.dumps(info, indent=2))

    @staticmethod
    def list_runs() -> list[Path]:
        if not RUNS_DIR.exists():
            return []
        runs = sorted(RUNS_DIR.iterdir(), reverse=True)
        return [r for r in runs if r.is_dir() and (r / "run_info.json").exists()]
