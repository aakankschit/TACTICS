# TACTICS: Thompson Sampling-Assisted Chemical Targeting and Iterative Compound Selection for Drug Discovery

TACTICS is a comprehensive library for Thompson Sampling-based optimization of chemical reaction libraries, featuring a unified architecture with flexible strategy selection, modern Pydantic configuration, and preset configurations for out-of-the-box usage.

## 🚀 Key Features

- **Unified Thompson Sampling Framework**: Single `ThompsonSampler` with pluggable selection strategies
- **Multiple Selection Strategies**:
  - Greedy (pure exploitation)
  - Roulette Wheel (adaptive thermal cycling)
  - UCB (Upper Confidence Bound)
  - Epsilon-Greedy (balanced exploration/exploitation)
  - Boltzmann (temperature-based selection)
- **Advanced Warmup Strategies**: Standard, Stratified, Enhanced, Latin Hypercube
- **Preset Configurations**: 5 ready-to-use presets for common use cases
- **Modern Pydantic Configuration**: Type-safe configuration with full validation
- **Parallel Processing**: Batch mode with multiprocessing for expensive evaluators
- **Multiple Evaluators**: Lookup, Database, ROCS, Fred, ML classifiers, and more
- **Library Enumeration**: Efficient generation of combinatorial reaction products
- **Library Analysis**: Comprehensive analysis and visualization tools
- **Polars DataFrames**: Fast, efficient data handling throughout

## 📦 Package Structure

```
TACTICS/
├── thompson_sampling/
│   ├── config.py          # Pydantic configuration models
│   ├── main.py            # run_ts() convenience wrapper
│   ├── presets.py         # Preset configurations
│   ├── core/              # Core unified sampler
│   │   ├── sampler.py         # ThompsonSampler (unified)
│   │   ├── evaluators.py      # All evaluator classes
│   │   └── evaluator_config.py # Evaluator Pydantic configs
│   ├── strategies/        # Selection strategies
│   │   ├── greedy.py
│   │   ├── roulette_wheel.py
│   │   ├── ucb.py
│   │   ├── epsilon_greedy.py
│   │   └── config.py      # Strategy Pydantic configs
│   ├── warmup/            # Warmup strategies
│   │   └── config.py      # Warmup Pydantic configs
│   ├── baseline.py        # Random baseline sampling
│   └── legacy/            # Legacy code (deprecated)
├── library_enumeration/   # Library generation tools
└── library_analysis/      # Analysis and visualization
```

## 🎯 Quick Start

### Simple Out-of-the-Box Usage with Presets (Recommended)

The easiest way to get started is using presets - just provide your reaction and reagents:

```python
from TACTICS.thompson_sampling import run_ts, get_preset
from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

# 1. Create evaluator config
evaluator = LookupEvaluatorConfig(ref_filename="scores.csv")

# 2. Get a preset configuration
config = get_preset(
    "fast_exploration",  # Quick screening with epsilon-greedy
    reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
    reagent_file_list=["aldehydes.smi", "amines.smi"],
    evaluator_config=evaluator,
    mode="minimize",  # Use "minimize" for docking scores
    num_iterations=1000
)

# 3. Run and get results
results_df = run_ts(config)
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
from TACTICS.thompson_sampling import run_ts, get_preset
from TACTICS.thompson_sampling.core.evaluator_config import FredEvaluatorConfig

# Configure slow evaluator
evaluator = FredEvaluatorConfig(design_unit_file="receptor.oedu")

# Get parallel batch preset
config = get_preset(
    "parallel_batch",
    reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
    reagent_file_list=["aldehydes.smi", "amines.smi"],
    evaluator_config=evaluator,
    mode="minimize",  # Docking scores
    num_iterations=1000,
    batch_size=100,   # Sample 100 compounds per cycle
    processes=8       # Use 8 CPU cores
)

results_df = run_ts(config)
```

### Custom Configuration (Advanced)

For full control, create custom configurations:

```python
from TACTICS.thompson_sampling import ThompsonSamplingConfig, run_ts
from TACTICS.thompson_sampling.strategies.config import EpsilonGreedyConfig
from TACTICS.thompson_sampling.warmup.config import StratifiedWarmupConfig
from TACTICS.thompson_sampling.core.evaluator_config import DBEvaluatorConfig

config = ThompsonSamplingConfig(
    reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
    reagent_file_list=["aldehydes.smi", "amines.smi"],
    num_ts_iterations=1000,
    num_warmup_trials=5,
    strategy_config=EpsilonGreedyConfig(mode="maximize", epsilon=0.2, decay=0.995),
    warmup_config=StratifiedWarmupConfig(),
    evaluator_config=DBEvaluatorConfig(db_filename="scores.db"),
    results_filename="results.csv",
    log_filename="run.log"
)

results_df = run_ts(config)
```

## 🔧 Configuration

### Pydantic Configuration Models

The package uses Pydantic for robust configuration validation:

```python
from TACTICS.thompson_sampling import StandardSamplerConfig, EnhancedSamplerConfig

# Automatic validation and type checking
config = StandardSamplerConfig(
    sampler_type="standard",  # Must be "standard"
    ts_mode="maximize",       # Must be one of: maximize, minimize, maximize_boltzmann, minimize_boltzmann
    processes=4,              # Must be positive integer
    # ... other parameters
)
```

### Configuration Validation

```python
# Invalid configuration raises ValidationError
try:
    config = StandardSamplerConfig(
        sampler_type="invalid",  # ❌ ValidationError
        ts_mode="maximize",
        # ...
    )
except ValidationError as e:
    print(f"Configuration error: {e}")
```

## 🧪 Testing

The package includes comprehensive tests for configuration validation:

```bash
# Run all tests
pytest tests/

# Run configuration tests
pytest tests/test_config_validation.py -v

# Run with coverage
pytest tests/ --cov=TACTICS --cov-report=html
```

## 📚 Documentation

- **API Documentation**: See `docs/` for detailed API documentation
- **Examples**: Check `examples/` for usage examples
- **Configuration Guide**: See `TACTICS/README.md` for detailed configuration options

## 🛠️ Installation

```bash
# Clone repository
git clone https://github.com/your-org/TACTICS.git
cd TACTICS

# Install dependencies
pip install -r requirements.txt

# Install package in development mode
pip install -e .
```

## 📋 Requirements

- Python 3.8+
- Pydantic
- RDKit
- NumPy
- Pandas
- Polars
- tqdm
- Multiprocessing support

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Ensure all tests pass
6. Submit a pull request

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 📖 Citation

If you use TACTICS in your research, please cite:

```bibtex
@software{tactics,
    title={TACTICS: Thompson Sampling-Assisted Chemical Targeting and Iterative Compound Selection for Drug Discovery},
    author={Aakankschit Nandkeolyar},
    year={2025},
    url={https://github.com/your-org/TACTICS}
}
```

## 🆘 Support

For questions and support:
- Open an issue on GitHub
- Contact: anandkeo@uci.edu

---

**Key Recommendation**: Use the modern Pydantic configuration approach for new projects. The legacy interface is maintained for backward compatibility.
