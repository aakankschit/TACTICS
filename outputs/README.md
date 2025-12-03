# TACTICS Experiment Outputs

Generated outputs from experiments and examples (GITIGNORED, regenerable).

## Structure

```
outputs/
└── runs/
    └── {YYYY-MM-DD}_{experiment_name}/
        ├── config.yaml
        ├── logs/
        ├── figures/
        └── results.csv
```

## Usage

All outputs are automatically organized by `experiment_utils.py`:

```python
from experiment_utils import ExperimentRun

# Automatically creates timestamped output directory
run = ExperimentRun("strategy_comparison", dataset="thrombin")

# Outputs saved to: outputs/runs/2025-12-03_strategy_comparison/
```

## Notes

- **Everything in `outputs/` is gitignored** (regenerable)
- Run experiments to regenerate any needed outputs
- Organized by date and experiment name for easy tracking
- Each run directory is self-contained with config and results
