# TACTICS Experiments

Research scripts for manuscript preparation and benchmarking (NOT distributed with package).

## Structure

- **`strategy_benchmarks/`**: Compare different selection strategies
  - RWS vs Bayes-UCB comparisons
  - Strategy recovery analysis
  - Legacy vs current validation

- **`warmup_analysis/`**: Test warmup configurations
  - Warmup variance analysis
  - Balanced sampling studies

- **`validation/`**: Validation and verification scripts
  - Legacy system validation
  - ROCS comparison studies

- **`_templates/`**: Template scripts for new experiments

- **`_debug/`**: Debug and diagnostic scripts (temporary)

## Usage

Use `experiment_utils.py` at the repository root for output management:

```python
from experiment_utils import ExperimentRun

# Create organized experiment run
run = ExperimentRun("my_experiment", dataset="thrombin")

# Outputs automatically organized in outputs/runs/{date}_{name}/
```

## Adding New Experiments

1. Copy template from `_templates/experiment_template.py`
2. Update paths to use `data/` structure
3. Use `experiment_utils` for output management
4. Categorize into appropriate subdirectory

## Notes

- All experiment scripts are **research-only**, not distributed
- Use `experiment_utils.py` at root for consistent output handling
- Outputs go to `outputs/runs/` (gitignored)
