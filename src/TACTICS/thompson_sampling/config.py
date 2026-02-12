"""
Thompson Sampling configuration.

Provides Pydantic configuration classes for Thompson Sampling optimization.
"""

from pydantic import BaseModel, Field, field_validator, model_validator
from typing import Literal, Optional, Union, List, TYPE_CHECKING, Any
import logging

from .strategies.config import (
    GreedyConfig,
    RouletteWheelConfig,
    UCBConfig,
    EpsilonGreedyConfig,
    BoltzmannConfig,
    BayesUCBConfig,
)
from .warmup.config import (
    StandardWarmupConfig,
    EnhancedWarmupConfig,
    BalancedWarmupConfig,
)
from .core.evaluator_config import (
    LookupEvaluatorConfig,
    DBEvaluatorConfig,
    FPEvaluatorConfig,
    MWEvaluatorConfig,
    ROCSEvaluatorConfig,
    FredEvaluatorConfig,
    MLClassifierEvaluatorConfig,
    CustomEvaluatorConfig,
)

if TYPE_CHECKING:
    from ..library_enumeration import SynthesisPipeline


# Type aliases for nested configs
StrategyConfigType = Union[
    GreedyConfig,
    RouletteWheelConfig,
    UCBConfig,
    EpsilonGreedyConfig,
    BoltzmannConfig,
    BayesUCBConfig,
]

WarmupConfigType = Union[
    BalancedWarmupConfig, StandardWarmupConfig, EnhancedWarmupConfig
]

EvaluatorConfigType = Union[
    LookupEvaluatorConfig,
    DBEvaluatorConfig,
    FPEvaluatorConfig,
    MWEvaluatorConfig,
    ROCSEvaluatorConfig,
    FredEvaluatorConfig,
    MLClassifierEvaluatorConfig,
    CustomEvaluatorConfig,
]


class ThompsonSamplingConfig(BaseModel):
    """
    Configuration for Thompson Sampling optimization.

    The SynthesisPipeline is the single source of truth for:
    - Reaction SMARTS patterns
    - Reagent file paths
    - Multi-step synthesis configuration

    Example:
        >>> from TACTICS.library_enumeration import (
        ...     SynthesisPipeline, ReactionDef, ReactionConfig
        ... )
        >>> from TACTICS.thompson_sampling import ThompsonSamplingConfig
        >>> from TACTICS.thompson_sampling.strategies.config import GreedyConfig
        >>> from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig
        >>>
        >>> # Create reaction config
        >>> config = ReactionConfig(
        ...     reactions=[ReactionDef(reaction_smarts="...", step_index=0)],
        ...     reagent_file_list=["acids.smi", "amines.smi"]
        ... )
        >>>
        >>> # Create pipeline
        >>> pipeline = SynthesisPipeline(config)
        >>>
        >>> # Create Thompson Sampling config
        >>> ts_config = ThompsonSamplingConfig(
        ...     synthesis_pipeline=pipeline,
        ...     num_ts_iterations=1000,
        ...     strategy_config=GreedyConfig(mode="maximize"),
        ...     evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv"),
        ... )
    """

    # Required: SynthesisPipeline contains all reaction and reagent info
    synthesis_pipeline: Any = Field(
        ...,
        description="SynthesisPipeline instance containing reaction config and reagent files",
    )

    # Thompson Sampling parameters
    num_ts_iterations: int = Field(
        ..., gt=0, description="Number of Thompson Sampling iterations"
    )
    num_warmup_trials: int = Field(
        default=3, gt=0, description="Number of warmup trials per reagent"
    )

    # Strategy and evaluator configs (required)
    strategy_config: StrategyConfigType = Field(
        ..., description="Selection strategy configuration"
    )
    evaluator_config: EvaluatorConfigType = Field(
        ..., description="Evaluator configuration"
    )

    # Optional warmup config
    warmup_config: Optional[WarmupConfigType] = Field(
        default=None, description="Warmup strategy configuration"
    )

    # Batch sampling parameters
    batch_size: int = Field(
        default=1, gt=0, description="Compounds to sample per iteration"
    )
    max_resamples: Optional[int] = Field(
        default=None, gt=0, description="Max resampling attempts for duplicates"
    )

    # Output
    results_filename: Optional[str] = Field(
        default="results.csv", description="Results output file"
    )
    log_filename: Optional[str] = Field(default=None, description="Log file path")

    # Performance (Note: multiprocessing typically not needed for compound generation)
    hide_progress: bool = Field(default=False, description="Hide progress bars")

    # Pre-enumerated library (optional, for testing)
    product_library_file: Optional[str] = Field(
        default=None,
        description="Path to CSV with pre-enumerated products (Product_Code, SMILES)",
    )

    # Bayesian update method
    use_boltzmann_weighting: bool = Field(
        default=False,
        description="Use Boltzmann-weighted Bayesian updates (legacy RWS)",
    )

    # SMARTS compatibility detection
    auto_detect_smarts_compatibility: bool = Field(
        default=True,
        description="Auto-detect reagent-pattern compatibility",
    )
    deprotect_for_compatibility: bool = Field(
        default=False,
        description="Apply deprotection during compatibility detection",
    )
    desalt_for_compatibility: bool = Field(
        default=False,
        description="Apply desalting during compatibility detection",
    )

    @field_validator("synthesis_pipeline")
    @classmethod
    def validate_pipeline(cls, v):
        """Ensure synthesis_pipeline is valid."""
        if v is None:
            raise ValueError("synthesis_pipeline is required")
        # Check it has required attributes
        if not hasattr(v, "reagent_file_list"):
            raise ValueError("synthesis_pipeline must have reagent_file_list attribute")
        if not hasattr(v, "enumerate_single"):
            raise ValueError("synthesis_pipeline must have enumerate_single method")
        return v

    def model_post_init(self, __context):
        """Set defaults and validate configuration."""
        # Set default warmup if not provided
        if self.warmup_config is None:
            object.__setattr__(self, "warmup_config", BalancedWarmupConfig())

    # === Convenience Properties ===

    @property
    def reagent_file_list(self) -> List[str]:
        """Get reagent file list from pipeline."""
        return self.synthesis_pipeline.reagent_file_list

    @property
    def num_components(self) -> int:
        """Get number of reagent components."""
        return self.synthesis_pipeline.num_components

    @property
    def num_steps(self) -> int:
        """Get number of reaction steps."""
        return self.synthesis_pipeline.num_steps

    model_config = {"arbitrary_types_allowed": True}


class RandomBaselineConfig(BaseModel):
    """
    Configuration for random baseline sampling.

    Used for comparison against Thompson Sampling.
    """

    synthesis_pipeline: Any = Field(..., description="SynthesisPipeline instance")
    evaluator_config: EvaluatorConfigType = Field(
        ..., description="Evaluator configuration"
    )
    num_trials: int = Field(..., gt=0, description="Number of random trials")
    num_to_save: int = Field(..., gt=0, description="Number of top results to save")
    ascending_output: bool = Field(default=False, description="Sort output ascending")
    outfile_name: Optional[str] = Field(default=None, description="Output file")
    log_filename: Optional[str] = Field(default=None, description="Log file")

    @field_validator("synthesis_pipeline")
    @classmethod
    def validate_pipeline(cls, v):
        """Ensure synthesis_pipeline is valid."""
        if v is None:
            raise ValueError("synthesis_pipeline is required")
        return v

    model_config = {"arbitrary_types_allowed": True}
