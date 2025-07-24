# TACTICS: Thompson Sampling-Assisted Chemical Targeting and Iterative Compound Selection for Drug Discovery

TACTICS is a comprehensive library for Thompson Sampling-based optimization of chemical reaction libraries, featuring modern Pydantic configuration and improved package organization.

## 🚀 Key Features

- **Modern Configuration**: Pydantic-based configuration with type safety and validation
- **Two Thompson Sampling Strategies**: Standard (greedy selection) and Enhanced (thermal cycling)
- **Library Enumeration**: Efficient generation of combinatorial reaction products
- **Library Analysis**: Comprehensive analysis and visualization tools
- **Advanced Search Methods**: Including DisallowTracker-enhanced strategies for guaranteed uniqueness
- **Clean Package Structure**: Well-organized modules with clear separation of concerns

## 📦 Package Structure

```
TACTICS/
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
from TACTICS.thompson_sampling import StandardSamplerConfig, run_ts

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

### Enhanced Thompson Sampling (Thermal Cycling)

```python
from TACTICS.thompson_sampling import EnhancedSamplerConfig

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
    author={Your Name},
    year={2024},
    url={https://github.com/your-org/TACTICS}
}
```

## 🆘 Support

For questions and support:
- Open an issue on GitHub
- Contact: your-email@institution.edu

---

**Key Recommendation**: Use the modern Pydantic configuration approach for new projects. The legacy interface is maintained for backward compatibility.
