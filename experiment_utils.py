"""
Experiment utilities for TACTICS research.

NOTE: This module is NOT part of the distributed TACTICS package.
It lives at the project root and is used only for personal research
and manuscript preparation.

Usage:
    # From project root or experiments/ folder
    import sys
    sys.path.insert(0, str(Path(__file__).parent))  # if needed
    
    from experiment_utils import ExperimentRun, load_dataset
    
    run = ExperimentRun("rws_benchmark", dataset="thrombin")
    run.save_config(config.model_dump())
    
    # ... run experiment ...
    
    fig.savefig(run.figure_path("convergence"))
    run.log.info("Experiment complete")
"""

from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, Any, Union
import yaml
import logging
import json
import sys


# =============================================================================
# Path Configuration
# =============================================================================

def get_project_root() -> Path:
    """Find the TACTICS project root directory."""
    # Start from this file's location and walk up
    current = Path(__file__).resolve().parent
    
    while current != current.parent:
        if (current / "pyproject.toml").exists():
            return current
        current = current.parent
    
    # Fallback: assume we're running from project root
    return Path.cwd()


PROJECT_ROOT = get_project_root()
DATA_DIR = PROJECT_ROOT / "data"
OUTPUTS_DIR = PROJECT_ROOT / "outputs" / "runs"
EXPERIMENTS_DIR = PROJECT_ROOT / "experiments"


# =============================================================================
# Dataset Management
# =============================================================================

class Dataset:
    """Represents a dataset with reagents, scores, and optional enumerated library."""
    
    def __init__(self, name: str):
        self.name = name
        self.reagents_dir = DATA_DIR / "reagents" / name
        self.scores_dir = DATA_DIR / "scores" / name
        self.library_dir = DATA_DIR / "libraries" / name
        
        if not self.reagents_dir.exists():
            available = [d.name for d in (DATA_DIR / "reagents").iterdir() if d.is_dir()]
            raise ValueError(
                f"Dataset '{name}' not found. Available datasets: {available}"
            )
    
    @property
    def reagent_files(self) -> Dict[str, Path]:
        """Return dict of reagent name -> file path."""
        return {
            f.stem: f 
            for f in self.reagents_dir.glob("*.smi")
        }
    
    @property
    def score_file(self) -> Optional[Path]:
        """Return the primary score file if it exists."""
        candidates = [
            self.scores_dir / "scores.csv",
            self.scores_dir / "docking_scores.csv",
            self.scores_dir / "rocs_scores.csv",
            self.scores_dir / "smiles_scores.csv",
        ]
        for c in candidates:
            if c.exists():
                return c
        
        # Try any CSV in the directory
        csvs = list(self.scores_dir.glob("*.csv"))
        return csvs[0] if csvs else None
    
    @property
    def has_scores(self) -> bool:
        return self.score_file is not None
    
    @property
    def has_enumerated_library(self) -> bool:
        return self.library_dir.exists() and any(self.library_dir.glob("*.csv"))
    
    def get_smarts(self) -> Optional[str]:
        """Read SMARTS pattern from dataset README if available."""
        readme = self.reagents_dir / "README.md"
        if not readme.exists():
            return None
        
        content = readme.read_text()
        # Simple extraction - looks for SMARTS: `pattern`
        import re
        match = re.search(r"SMARTS:\s*`([^`]+)`", content)
        return match.group(1) if match else None
    
    def __repr__(self) -> str:
        return (
            f"Dataset('{self.name}', "
            f"reagents={list(self.reagent_files.keys())}, "
            f"has_scores={self.has_scores})"
        )


def list_datasets() -> list[str]:
    """List all available datasets."""
    reagents_dir = DATA_DIR / "reagents"
    if not reagents_dir.exists():
        return []
    return [d.name for d in reagents_dir.iterdir() if d.is_dir()]


def load_dataset(name: str) -> Dataset:
    """Load a dataset by name."""
    return Dataset(name)


# =============================================================================
# Experiment Run Management
# =============================================================================

class ExperimentRun:
    """
    Manages an experiment run with standardized output organization.
    
    Creates a timestamped directory structure:
        outputs/runs/{YYYY-MM-DD}_{name}/
        ├── config.yaml
        ├── logs/
        │   └── experiment.log
        └── figures/
    
    Usage:
        run = ExperimentRun("rws_benchmark")
        run.save_config({"iterations": 500, "seed": 42})
        
        run.log.info("Starting experiment")
        
        fig.savefig(run.figure_path("convergence"))
        
        run.save_results({"best_score": 9.5, "recovery": 0.85})
    """
    
    def __init__(
        self,
        name: str,
        dataset: Optional[str] = None,
        base_dir: Optional[Path] = None,
        timestamp: Optional[datetime] = None,
    ):
        """
        Initialize an experiment run.
        
        Args:
            name: Descriptive experiment name (e.g., "rws_thermal_sweep")
            dataset: Optional dataset name to associate with run
            base_dir: Override default outputs directory
            timestamp: Override timestamp (useful for resuming runs)
        """
        self.name = name
        self.dataset_name = dataset
        self.timestamp = timestamp or datetime.now()
        
        # Set up directories
        base = base_dir or OUTPUTS_DIR
        self.run_dir = base / f"{self.timestamp:%Y-%m-%d}_{name}"
        self.logs_dir = self.run_dir / "logs"
        self.figures_dir = self.run_dir / "figures"
        
        # Create directories
        self.logs_dir.mkdir(parents=True, exist_ok=True)
        self.figures_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize logger
        self._setup_logging()
        
        # Load dataset if specified
        self.dataset = Dataset(dataset) if dataset else None
        
        # Track metadata
        self._metadata = {
            "experiment_name": name,
            "dataset": dataset,
            "timestamp": self.timestamp.isoformat(),
            "run_dir": str(self.run_dir),
        }
        
        # Announce
        print(f"Experiment outputs: {self.run_dir}")
        self.log.info(f"Initialized experiment: {name}")
        if dataset:
            self.log.info(f"Dataset: {dataset}")
    
    def _setup_logging(self):
        """Configure logging to file and console."""
        self.log = logging.getLogger(f"TACTICS.{self.name}")
        self.log.setLevel(logging.DEBUG)
        
        # Avoid duplicate handlers if re-initialized
        if self.log.handlers:
            self.log.handlers.clear()
        
        # File handler (DEBUG level - everything)
        log_file = self.logs_dir / "experiment.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s'
        ))
        self.log.addHandler(file_handler)
        
        # Console handler (INFO level)
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(logging.Formatter(
            '%(levelname)s: %(message)s'
        ))
        self.log.addHandler(console_handler)
    
    def save_config(self, config: Union[Dict[str, Any], Any]):
        """
        Save experiment configuration to config.yaml.
        
        Args:
            config: Configuration dict or Pydantic model (with .model_dump())
        """
        if hasattr(config, 'model_dump'):
            config_dict = config.model_dump()
        else:
            config_dict = config
        
        full_config = {
            **self._metadata,
            "parameters": config_dict,
        }
        
        config_path = self.run_dir / "config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(full_config, f, default_flow_style=False, sort_keys=False)
        
        self.log.debug(f"Saved config to {config_path}")
    
    def save_results(self, results: Dict[str, Any], filename: str = "results.yaml"):
        """Save experiment results."""
        results_path = self.run_dir / filename
        with open(results_path, "w") as f:
            yaml.dump(results, f, default_flow_style=False)
        self.log.info(f"Saved results to {results_path}")
    
    def save_dataframe(self, df, filename: str):
        """Save a pandas/polars DataFrame to CSV."""
        path = self.run_dir / filename
        df.to_csv(path, index=False)
        self.log.info(f"Saved DataFrame to {path}")
        return path
    
    def figure_path(self, name: str, ext: str = "png") -> Path:
        """Get path for saving a figure."""
        return self.figures_dir / f"{name}.{ext}"
    
    def log_path(self, name: str = "experiment") -> Path:
        """Get path for a log file."""
        return self.logs_dir / f"{name}.log"
    
    def add_detailed_logger(self, name: str) -> logging.Logger:
        """
        Create an additional logger for detailed output (e.g., TS iterations).
        
        Returns a logger that writes to logs/{name}.log
        """
        logger = logging.getLogger(f"TACTICS.{self.name}.{name}")
        logger.setLevel(logging.DEBUG)
        
        handler = logging.FileHandler(self.log_path(name))
        handler.setFormatter(logging.Formatter(
            '%(asctime)s - %(message)s'
        ))
        logger.addHandler(handler)
        
        return logger
    
    @classmethod
    def from_existing(cls, run_dir: Union[str, Path]) -> "ExperimentRun":
        """Load an existing run for analysis."""
        run_dir = Path(run_dir)
        if not run_dir.exists():
            raise ValueError(f"Run directory not found: {run_dir}")
        
        # Parse name and timestamp from directory name
        dir_name = run_dir.name
        date_str, name = dir_name.split("_", 1)
        timestamp = datetime.strptime(date_str, "%Y-%m-%d")
        
        # Load config to get dataset
        config_path = run_dir / "config.yaml"
        dataset = None
        if config_path.exists():
            with open(config_path) as f:
                config = yaml.safe_load(f)
                dataset = config.get("dataset")
        
        return cls(
            name=name,
            dataset=dataset,
            base_dir=run_dir.parent,
            timestamp=timestamp,
        )
    
    def __repr__(self) -> str:
        return f"ExperimentRun('{self.name}', dir='{self.run_dir}')"


# =============================================================================
# Convenience Functions
# =============================================================================

def list_runs(n: int = 10) -> list[Path]:
    """List recent experiment runs."""
    if not OUTPUTS_DIR.exists():
        return []
    
    runs = sorted(OUTPUTS_DIR.iterdir(), reverse=True)
    return runs[:n]


def load_run(run_dir: Union[str, Path]) -> ExperimentRun:
    """Load an existing run for analysis."""
    return ExperimentRun.from_existing(run_dir)


def compare_configs(run1: Union[str, Path], run2: Union[str, Path]) -> Dict[str, Any]:
    """Compare configurations between two runs."""
    def load_config(path):
        config_file = Path(path) / "config.yaml"
        with open(config_file) as f:
            return yaml.safe_load(f)
    
    c1 = load_config(run1)
    c2 = load_config(run2)
    
    # Find differences
    diffs = {}
    all_keys = set(c1.get("parameters", {}).keys()) | set(c2.get("parameters", {}).keys())
    
    for key in all_keys:
        v1 = c1.get("parameters", {}).get(key)
        v2 = c2.get("parameters", {}).get(key)
        if v1 != v2:
            diffs[key] = {"run1": v1, "run2": v2}
    
    return diffs