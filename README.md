# PRISMS: Probabilistic Rapid Intelligent Screening of Molecular Space

A comprehensive library for Thompson Sampling-based optimization of chemical combinatorial libraries, featuring modern Pydantic configuration and improved package organization.

## 🚀 Key Features

- **Modern Configuration**: Pydantic-based configuration with type safety and validation
- **Two Thompson Sampling Strategies**: Standard (greedy selection) and Enhanced (thermal cycling)
- **Random Baseline Sampling**: Configurable random sampling for comparison studies
- **Library Enumeration**: Efficient generation of combinatorial reaction products
- **Library Analysis**: Comprehensive analysis and visualization tools
- **Advanced Search Methods**: Including DisallowTracker-enhanced strategies for guaranteed uniqueness
- **Clean Package Structure**: Well-organized modules with clear separation of concerns

## 📦 Package Structure

```
PRISMS/
├── thompson_sampling/
│   ├── config.py          # Pydantic configuration models
│   ├── main.py            # Main execution interface
│   ├── core/              # Core sampling algorithms
│   │   ├── standard_sampler.py
│   │   ├── enhanced_sampler.py
│   │   └── evaluators.py
│   ├── utils/             # Utility functions
│   │   ├── ts_logger.py
│   │   └── ts_utils.py
│   └── legacy/            # Original modules (backward compatibility)
├── library_enumeration/   # Library generation tools
└── library_analysis/      # Analysis and visualization
```

## 🎯 Quick Start

### Standard Thompson Sampling (Greedy Selection)

```python
from PRISMS.thompson_sampling import StandardSamplerConfig, run_ts

# Create configuration using Pydantic models
config = StandardSamplerConfig(
    sampler_type="standard",
    ts_mode="maximize",
    evaluator_class_name="DBEvaluator",
    evaluator_arg="scores.csv",
    reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
    num_ts_iterations=100,
    reagent_file_list=["aldehydes.smi", "amines.smi"],
    num_warmup_trials=3,
    results_filename="results.csv"
)

# Run Thompson Sampling with greedy selection
results_df = run_ts(config)
```

### Random Baseline Sampling

```python
from PRISMS.thompson_sampling import RandomBaselineConfig, run_random_baseline

# Create configuration for random baseline
config = RandomBaselineConfig(
    evaluator_class_name="LookupEvaluator",
    evaluator_arg={
        "ref_filename": "scores.csv",
        "ref_colname": "Score"
    },
    reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
    reagent_file_list=["aldehydes.smi", "amines.smi"],
    num_trials=1000,
    num_to_save=100,
    ascending_output=False,
    outfile_name="random_results.csv"
)

# Run random baseline sampling
results_df = run_random_baseline(config)
```

### Enhanced Thompson Sampling (Thermal Cycling)

```python
from PRISMS.thompson_sampling import EnhancedSamplerConfig

config = EnhancedSamplerConfig(
    sampler_type="enhanced",
    processes=4,
    scaling=1.0,
    percent_of_library=0.1,
    minimum_no_of_compounds_per_core=10,
    stopping_criteria=1000,
    evaluator_class_name="DBEvaluator",
    evaluator_arg="scores.csv",
    reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
    num_ts_iterations=100,
    reagent_file_list=["aldehydes.smi", "amines.smi"],
    num_warmup_trials=3
)

# Run Thompson Sampling with thermal cycling
results_df = run_ts(config)
```

## 🔧 Configuration

### Pydantic Configuration Models

The package uses Pydantic for robust configuration validation:

```python
from PRISMS.thompson_sampling import StandardSamplerConfig, EnhancedSamplerConfig, RandomBaselineConfig

# Automatic validation and type checking
config = StandardSamplerConfig(
    sampler_type="standard",  # Must be "standard"
    ts_mode="maximize",       # Must be one of: maximize, minimize, maximize_boltzmann, minimize_boltzmann
    processes=4,              # Must be positive integer
    # ... other parameters
)

# Random baseline configuration
random_config = RandomBaselineConfig(
    evaluator_class_name="LookupEvaluator",
    evaluator_arg={"ref_filename": "scores.csv", "ref_colname": "Score"},
    num_trials=1000,         # Must be positive integer
    num_to_save=100,         # Must be positive integer
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
pytest tests/ --cov=PRISMS --cov-report=html
```

## 📚 Documentation

- **API Documentation**: See `docs/` for detailed API documentation
- **Examples**: Check `examples/` for usage examples
- **Configuration Guide**: See `PRISMS/README.md` for detailed configuration options

## 🛠️ Installation

```bash
# Clone repository
git clone https://github.com/your-org/PRISMS.git
cd PRISMS

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

If you use PRISMS in your research, please cite:

```bibtex
@software{prisms,
    title={PRISMS: Probabilistic Rapid Intelligent Screening of Molecular Space},
    author={Your Name},
    year={2024},
    url={https://github.com/your-org/PRISMS}
}
```

## 🆘 Support

For questions and support:
- Open an issue on GitHub
- Contact: your-email@institution.edu

---

**This work is based on [previous work](https://pubs.acs.org/doi/10.1021/acs.jcim.3c01790) by Patrick Walters.**
**This project is a collaboration between the University of California Irvine and Leiden University.**