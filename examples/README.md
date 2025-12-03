# TACTICS Examples

Clean, user-facing examples demonstrating TACTICS functionality.

## Quick Start Examples

- **`quickstart_thrombin.py`** - Minimal working example with thrombin dataset
- **`bayes_ucb_example.py`** - Demonstrates Bayes-UCB selection strategy
- **`simple_multi_smarts_demo.py`** - Multiple reaction patterns example

## Jupyter Notebooks

Interactive tutorials in `notebooks/`:

- **`adenine_dataset.ipynb`** - Adenine receptor dataset walkthrough
- **`quinazoline_dataset.ipynb`** - Quinazoline dataset analysis
- **`seh_analysis.ipynb`** - sEH (DEL) dataset exploration
- **`smarts_toolkit.ipynb`** - SMARTS pattern validation tools

## Usage

All examples use the reorganized data structure:

```python
# Data files are in data/ directory
reagent_files = [
    "data/reagents/thrombin/acids.smi",
    "data/reagents/thrombin/coupled_aa_sub.smi"
]

score_file = "data/scores/thrombin/product_scores.csv"
```

## For Research/Benchmarking

Research scripts and comparisons have been moved to `experiments/` directory:
- Strategy benchmarks → `experiments/strategy_benchmarks/`
- Warmup analysis → `experiments/warmup_analysis/`
- Validation studies → `experiments/validation/`

## Notes

- Examples are **user-facing** and **distributed** with the package
- Focus on demonstrating TACTICS functionality clearly
- Keep examples simple and well-documented
- Research code belongs in `experiments/`, not here
