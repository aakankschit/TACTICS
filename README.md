# TACTICS: Thompson Sampling-Assisted Chemical Targeting and Iterative Compound Selection for Drug Discovery

<p align="center">
  <img src="https://raw.githubusercontent.com/aakankschit/TACTICS/main/docs/source/_static/images/TACTICS_logo.png" alt="TACTICS Logo" width="600">
</p>

A comprehensive library for Thompson Sampling-based optimization of chemical combinatorial libraries, featuring a unified architecture with flexible strategy selection, modern Pydantic configuration, and preset configurations for out-of-the-box usage.

## Quick Start with Interactive Tutorials

TACTICS includes interactive [marimo](https://marimo.io) notebooks for learning and exploration.
For full documentation, see the [TACTICS Documentation](https://aakankschit.github.io/TACTICS/).

### Installation

```bash
pip install chem-tactics[tutorials]  # Includes marimo
```

### Running Tutorials

**As an interactive app** (recommended for exploration):
```bash
marimo run tutorials/thompson_sampling_tutorial.py
```

**In edit mode** (for learning/modification):
```bash
marimo edit tutorials/thompson_sampling_tutorial.py
```

### Available Tutorials

| Tutorial | Description |
|----------|-------------|
| `library_enumeration_tutorial.py` | SynthesisPipeline and enumeration |
| `thompson_sampling_tutorial.py` | Selection strategies comparison |
| `reaction_config_builder.py` | ReactionConfig builder |
| `library_component_comparison.py` | Library component analysis |
| `legacy_vs_current_comparison.py` | Legacy vs current benchmark |

> **Note**: Tutorials default to the bundled Thrombin dataset. Select "Local Data" mode to use your own files.

## Key Features

- **Unified Thompson Sampling Framework**: Single `ThompsonSampler` with pluggable selection strategies
- **Multiple Selection Strategies**:
  - Greedy (pure exploitation)
  - Roulette Wheel (adaptive thermal cycling)
  - UCB (Upper Confidence Bound)
  - Epsilon-Greedy (balanced exploration/exploitation)
  - Bayes-UCB (Bayesian upper confidence bound)
  - Boltzmann (temperature-based selection)
- **Warmup Strategies**: Balanced (recommended), Standard, Enhanced
- **Preset Configurations**: 5 ready-to-use presets for common use cases
- **Modern Pydantic Configuration**: Type-safe configuration with full validation
- **Parallel Processing**: Batch mode with multiprocessing for expensive evaluators
- **Multiple Evaluators**: Lookup, Database, ROCS, Fred, ML classifiers, and more
- **Synthesis Pipeline**: `SynthesisPipeline` architecture for single-step, alternative SMARTS, and multi-step reactions
- **SMARTS Toolkit**: `ReactionDef` with built-in validation, visualization, and protecting group support
- **Library Enumeration**: Efficient generation of combinatorial reaction products with `write_enumerated_library()`
- **Library Analysis**: Comprehensive analysis and visualization tools
- **Polars DataFrames**: Fast, efficient data handling throughout

## Package Structure

```
TACTICS/
├── thompson_sampling/
│   ├── config.py              # ThompsonSamplingConfig (Pydantic v2)
│   ├── presets.py             # Preset configurations
│   ├── factories.py           # Factory functions for component creation
│   ├── core/                  # Core unified sampler
│   │   ├── sampler.py         # ThompsonSampler (unified)
│   │   ├── evaluators.py      # All evaluator classes
│   │   └── evaluator_config.py # Evaluator Pydantic configs
│   ├── strategies/            # Selection strategies
│   │   ├── greedy.py
│   │   ├── roulette_wheel.py
│   │   ├── ucb.py
│   │   ├── epsilon_greedy.py
│   │   ├── bayes_ucb.py
│   │   └── config.py          # Strategy Pydantic configs
│   ├── warmup/                # Warmup strategies
│   │   └── config.py          # Warmup Pydantic configs (Balanced, Standard, Enhanced)
│   └── baseline.py            # Random baseline sampling
├── library_enumeration/       # Library generation tools
│   ├── synthesis_pipeline.py  # SynthesisPipeline - main entry point
│   ├── enumeration_utils.py   # EnumerationResult, EnumerationError
│   ├── file_writer.py         # write_enumerated_library()
│   ├── generate_products.py   # Product generation utilities
│   └── smarts_toolkit/        # SMARTS validation and configuration
│       ├── config.py          # ReactionDef, ReactionConfig, StepInput, DeprotectionSpec
│       ├── _validator.py      # ValidationResult, internal validation
│       └── constants.py       # Protecting groups, salt fragments
└── library_analysis/          # Analysis and visualization
```

## Repository Structure

```
TACTICS/
├── src/TACTICS/              # Core package (pip installable)
│   ├── thompson_sampling/    # Thompson Sampling algorithms
│   ├── library_enumeration/  # Library generation tools
│   ├── library_analysis/     # Analysis and visualization
│   └── data/                 # Bundled tutorial datasets
│       └── thrombin/         # Thrombin inhibitor dataset
│
├── tutorials/                # Interactive marimo tutorials
├── tests/                    # Unit and integration tests
└── docs/                     # Sphinx documentation
```

## Quick Start

### Simple Out-of-the-Box Usage with Presets (Recommended)

The easiest way to get started is using presets with `SynthesisPipeline`:

```python
from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef
from TACTICS.thompson_sampling import ThompsonSampler, get_preset
from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

# 1. Create synthesis pipeline (single source of truth for reactions)
rxn_config = ReactionConfig(
    reactions=[ReactionDef(
        reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
        step_index=0,
        description="Amide coupling"
    )],
    reagent_file_list=["acids.smi", "amines.smi"]
)
pipeline = SynthesisPipeline(rxn_config)

# 2. Create evaluator config
evaluator = LookupEvaluatorConfig(ref_filename="scores.csv")

# 3. Get a preset configuration
config = get_preset(
    "fast_exploration",  # Quick screening with epsilon-greedy
    synthesis_pipeline=pipeline,
    evaluator_config=evaluator,
    mode="minimize",  # Use "minimize" for docking scores
    num_iterations=1000
)

# 4. Create sampler from config and run
sampler = ThompsonSampler.from_config(config)
warmup_df = sampler.warm_up(num_warmup_trials=config.num_warmup_trials)
results_df = sampler.search(num_cycles=config.num_ts_iterations)
sampler.close()

# 5. Analyze top results
print(results_df.sort("score").head(10))
```

**Available Presets:**
- `"fast_exploration"` - Epsilon-greedy strategy, quick screening
- `"parallel_batch"` - Batch processing with multiprocessing (for slow evaluators)
- `"conservative_exploit"` - Greedy strategy, focus on best reagents
- `"balanced_sampling"` - UCB strategy with theoretical guarantees
- `"diverse_coverage"` - Maximum diversity exploration

### Parallel Batch Processing (for Expensive Evaluators)

For slow evaluators (docking, ML models), use batch mode with multiprocessing:

```python
from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef
from TACTICS.thompson_sampling import ThompsonSampler, get_preset
from TACTICS.thompson_sampling.core.evaluator_config import FredEvaluatorConfig

# Create synthesis pipeline
rxn_config = ReactionConfig(
    reactions=[ReactionDef(
        reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
        step_index=0
    )],
    reagent_file_list=["acids.smi", "amines.smi"]
)
pipeline = SynthesisPipeline(rxn_config)

# Configure slow evaluator (molecular docking)
evaluator = FredEvaluatorConfig(design_unit_file="receptor.oedu")

# Get parallel batch preset
config = get_preset(
    "parallel_batch",
    synthesis_pipeline=pipeline,
    evaluator_config=evaluator,
    mode="minimize",  # Docking scores (lower is better)
    batch_size=100,   # Sample 100 compounds per cycle
)

# Create sampler and run
sampler = ThompsonSampler.from_config(config)
warmup_df = sampler.warm_up(num_warmup_trials=config.num_warmup_trials)
results_df = sampler.search(num_cycles=config.num_ts_iterations)
sampler.close()
```

### Custom Configuration (Advanced)

For full control, create custom configurations:

```python
from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef
from TACTICS.thompson_sampling import ThompsonSampler, ThompsonSamplingConfig
from TACTICS.thompson_sampling.strategies.config import RouletteWheelConfig
from TACTICS.thompson_sampling.warmup.config import BalancedWarmupConfig
from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

# Create synthesis pipeline
rxn_config = ReactionConfig(
    reactions=[ReactionDef(
        reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
        step_index=0
    )],
    reagent_file_list=["acids.smi", "amines.smi"]
)
pipeline = SynthesisPipeline(rxn_config)

# Create fully customized configuration
config = ThompsonSamplingConfig(
    synthesis_pipeline=pipeline,
    num_ts_iterations=5000,
    num_warmup_trials=5,
    strategy_config=RouletteWheelConfig(
        mode="maximize",
        alpha=0.1,  # Initial heating temperature
        beta=0.1,   # Initial cooling temperature
    ),
    warmup_config=BalancedWarmupConfig(
        observations_per_reagent=5,
        use_per_reagent_variance=True,
    ),
    evaluator_config=LookupEvaluatorConfig(
        ref_filename="scores.csv",
        score_col="binding_affinity"
    ),
    batch_size=10,
    log_filename="optimization.log"
)

# Create sampler and run
sampler = ThompsonSampler.from_config(config)
warmup_df = sampler.warm_up(num_warmup_trials=config.num_warmup_trials)
results_df = sampler.search(num_cycles=config.num_ts_iterations)
sampler.close()

# Save results
results_df.write_csv("my_results.csv")
```

### Random Baseline Sampling

```python
from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef
from TACTICS.thompson_sampling import RandomBaselineConfig, run_random_baseline
from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

# Create synthesis pipeline
rxn_config = ReactionConfig(
    reactions=[ReactionDef(
        reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
        step_index=0
    )],
    reagent_file_list=["acids.smi", "amines.smi"]
)
pipeline = SynthesisPipeline(rxn_config)

config = RandomBaselineConfig(
    synthesis_pipeline=pipeline,
    evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv"),
    num_trials=1000,
    num_to_save=100,
    ascending_output=False,
    outfile_name="random_results.csv"
)

results_df = run_random_baseline(config)
```

## Configuration

### Pydantic Configuration Models

The package uses Pydantic v2 for robust configuration validation:

```python
from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef
from TACTICS.thompson_sampling import ThompsonSamplingConfig
from TACTICS.thompson_sampling.strategies.config import EpsilonGreedyConfig
from TACTICS.thompson_sampling.warmup.config import BalancedWarmupConfig
from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

# Create synthesis pipeline
rxn_config = ReactionConfig(
    reactions=[ReactionDef(
        reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
        step_index=0
    )],
    reagent_file_list=["acids.smi", "amines.smi"]
)
pipeline = SynthesisPipeline(rxn_config)

# Automatic validation and type checking
config = ThompsonSamplingConfig(
    synthesis_pipeline=pipeline,  # Required: single source of truth
    num_ts_iterations=1000,
    strategy_config=EpsilonGreedyConfig(mode="maximize", epsilon=0.2),
    warmup_config=BalancedWarmupConfig(),
    evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv"),
)
```

### Configuration Validation

```python
from pydantic import ValidationError

# Invalid configuration raises ValidationError
try:
    rxn = ReactionDef(
        reaction_smarts="invalid-smarts",  # ValidationError: Invalid SMARTS
        step_index=0
    )
except ValidationError as e:
    print(f"Configuration error: {e}")
```

## Testing

The package includes comprehensive tests for configuration validation:

```bash
# Run all tests
pytest tests/

# Run configuration tests
pytest tests/test_config_validation.py -v

# Run with coverage
pytest tests/ --cov=TACTICS --cov-report=html
```

## Documentation

- **Full Documentation**: [TACTICS Documentation](https://aakankschit.github.io/TACTICS/)
- **Interactive Tutorials**: See `tutorials/` for marimo notebooks
- **API Reference**: Build locally with `cd docs && make html`

## Installation

```bash
# Clone repository and install package in development mode
git clone https://github.com/aakankschit/TACTICS.git
cd TACTICS
pip install -e .

# With interactive tutorials (marimo):
pip install -e ".[tutorials]"

# With test dependencies:
pip install -e ".[test]"
```

## Requirements

- Python 3.11+
- Multiprocessing support

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Ensure all tests pass
6. Submit a pull request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use TACTICS in your research, please cite:

```bibtex
@software{tactics,
    title={TACTICS: Thompson Sampling-Assisted Chemical Targeting and Iterative Compound Selection for Drug Discovery},
    author={Aakankschit Nandkeolyar},
    year={2025},
    url={https://github.com/your-org/TACTICS}
}
```

## Support

For questions and support:
- Open an issue on GitHub
- Contact: anandkeo@uci.edu

---

**This work is based on [previous work](https://pubs.acs.org/doi/10.1021/acs.jcim.3c01790) by Patrick Walters.**
**This project is a collaboration between the University of California Irvine, Leiden University and Groningen University.**
