"""
Library Enumeration - Chemical library generation and enumeration.

This module provides:
- SynthesisPipeline: Main class for library enumeration
- File writing utilities for output formats
- Utility functions for product generation

Quick Start
-----------

**Create a pipeline and enumerate:**

    from TACTICS.library_enumeration import (
        SynthesisPipeline,
        ReactionDef,
        ReactionConfig,
        write_enumerated_library,
    )

    # Define reaction
    config = ReactionConfig(
        reactions=[ReactionDef(reaction_smarts="...", step_index=0)],
        reagent_file_list=["acids.smi", "amines.smi"]
    )

    # Create pipeline
    pipeline = SynthesisPipeline(config)

    # Enumerate library
    results = pipeline.enumerate_library(n_jobs=4, show_progress=True)

    # Write results
    write_enumerated_library(results, "products.csv", format="csv")

**For Thompson Sampling:**

    from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig
    from TACTICS.thompson_sampling import ThompsonSamplingConfig

    pipeline = SynthesisPipeline(config)
    ts_config = ThompsonSamplingConfig(
        synthesis_pipeline=pipeline,
        ...
    )
"""

# Main pipeline class
from .synthesis_pipeline import SynthesisPipeline

# Result dataclasses and analysis utilities
from .enumeration_utils import (
    EnumerationResult,
    EnumerationError,
    AutoDetectionResult,
    results_to_dataframe,
    failures_to_dataframe,
    summarize_failures,
)

# File writing utilities
from .file_writer import (
    write_enumerated_library,
    write_products_chunked,
    write_products_to_files,
)

# Product generation utilities
from .generate_products import (
    enumerate_products,
    generate_all_combinations,
    LibraryEnumerator,  # Legacy, deprecated
)

# Multiprocessing utilities
from .multiprocessing_utils import initializer

# Re-export smarts_toolkit classes for convenience
from .smarts_toolkit import (
    ReactionDef,
    ReactionConfig,
    DeprotectionSpec,
    StepInput,
    InputSource,
    ValidationResult,
    ProtectingGroupInfo,
)

__all__ = [
    # Main class
    "SynthesisPipeline",
    # Results
    "EnumerationResult",
    "EnumerationError",
    "AutoDetectionResult",
    # File utilities
    "write_enumerated_library",
    "write_products_chunked",
    "write_products_to_files",
    "results_to_dataframe",
    "failures_to_dataframe",
    "summarize_failures",
    # Product generation
    "enumerate_products",
    "generate_all_combinations",
    # Multiprocessing
    "initializer",
    # Legacy
    "LibraryEnumerator",
    # Config (re-exported from smarts_toolkit)
    "ReactionDef",
    "ReactionConfig",
    "DeprotectionSpec",
    "StepInput",
    "InputSource",
    "ValidationResult",
    "ProtectingGroupInfo",
]
