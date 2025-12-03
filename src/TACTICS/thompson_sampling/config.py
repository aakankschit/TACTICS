from pydantic import BaseModel, Field, field_validator
from typing import Literal, Optional, Union
import logging
from .strategies.config import (
    GreedyConfig,
    RouletteWheelConfig,
    UCBConfig,
    EpsilonGreedyConfig,
    BoltzmannConfig,
    BayesUCBConfig
)
from .warmup.config import (
    StandardWarmupConfig,
    EnhancedWarmupConfig,
    BalancedWarmupConfig
)
from .core.evaluator_config import (
    LookupEvaluatorConfig,
    DBEvaluatorConfig,
    FPEvaluatorConfig,
    MWEvaluatorConfig,
    ROCSEvaluatorConfig,
    FredEvaluatorConfig,
    MLClassifierEvaluatorConfig
)

# Import MultiStepSynthesisConfig for advanced reaction handling
try:
    from ..library_enumeration.smarts_toolkit.reaction_sequence import MultiStepSynthesisConfig
    REACTION_SEQUENCE_AVAILABLE = True
except ImportError:
    REACTION_SEQUENCE_AVAILABLE = False
    MultiStepSynthesisConfig = None

# Type aliases for nested configs
StrategyConfigType = Union[
    GreedyConfig,
    RouletteWheelConfig,
    UCBConfig,
    EpsilonGreedyConfig,
    BoltzmannConfig,
    BayesUCBConfig
]

WarmupConfigType = Union[
    BalancedWarmupConfig,
    StandardWarmupConfig,
    EnhancedWarmupConfig
]

EvaluatorConfigType = Union[
    LookupEvaluatorConfig,
    DBEvaluatorConfig,
    FPEvaluatorConfig,
    MWEvaluatorConfig,
    ROCSEvaluatorConfig,
    FredEvaluatorConfig,
    MLClassifierEvaluatorConfig
]

class ThompsonSamplingConfig(BaseModel):
    """
    Unified configuration for Thompson Sampling with nested component configs.

    This config supports multiple usage patterns:
    1. Legacy: Specify evaluator_class_name, evaluator_arg, selection_strategy (string)
    2. Modern: Use nested configs (strategy_config, warmup_config, evaluator_config)
    3. Advanced: Use reaction_sequence_config for multi-SMARTS or multi-step synthesis

    The modern approach with nested configs provides better validation and type safety.
    """

    # Core parameters - reaction specification
    reaction_smarts: Optional[str] = None  # Simple single SMARTS (backward compatible)
    reaction_sequence_config: Optional['MultiStepSynthesisConfig'] = Field(
        default=None,
        description="Advanced: Multi-SMARTS or multi-step synthesis configuration. "
                   "Takes precedence over reaction_smarts if both are provided."
    )
    reagent_file_list: list[str]
    num_ts_iterations: int
    num_warmup_trials: int = 3

    # Modern approach: Nested configs (recommended)
    strategy_config: Optional[StrategyConfigType] = None
    warmup_config: Optional[WarmupConfigType] = Field(default=None)
    evaluator_config: Optional[EvaluatorConfigType] = None

    # Legacy approach: String-based config (backward compatibility)
    evaluator_class_name: Optional[str] = None
    evaluator_arg: Optional[str | dict] = None
    selection_strategy: Optional[Literal[
        "greedy",
        "roulette_wheel",
        "ucb",
        "epsilon_greedy",
        "boltzmann",
        "bayes_ucb"
    ]] = None
    mode: Optional[Literal["maximize", "minimize"]] = None
    strategy_params: Optional[dict] = Field(default_factory=dict)

    # Batch sampling parameters
    batch_size: int = Field(default=1, gt=0)
    max_resamples: Optional[int] = Field(default=None, gt=0)

    # Output
    results_filename: Optional[str] = "results.csv"
    log_filename: Optional[str] = None

    # Performance / Multiprocessing
    processes: int = Field(default=1, gt=0)
    min_cpds_per_core: int = Field(default=10, gt=0)
    hide_progress: bool = False

    # Pre-enumerated library (optional, for testing)
    product_library_file: Optional[str] = Field(
        default=None,
        description="Path to CSV with Product_Code and SMILES columns for pre-enumerated products"
    )

    # Bayesian update method
    use_boltzmann_weighting: bool = Field(
        default=False,
        description="If True, use Boltzmann-weighted Bayesian updates (legacy RWS algorithm). "
                   "If False, use standard uniform-weighted Bayesian updates (default)."
    )

    @field_validator('batch_size')
    @classmethod
    def validate_batch_size(cls, v):
        """Ensure batch_size is positive."""
        if v < 1:
            raise ValueError('batch_size must be >= 1')
        return v

    @field_validator('processes')
    @classmethod
    def validate_processes(cls, v):
        """Ensure processes is positive."""
        if v < 1:
            raise ValueError('processes must be >= 1')
        return v

    def model_post_init(self, __context):
        """
        Validate that either modern or legacy config is provided (not both).
        Also validate reaction configuration.
        """
        modern_provided = (
            self.strategy_config is not None or
            self.evaluator_config is not None
        )
        legacy_provided = (
            self.evaluator_class_name is not None or
            self.selection_strategy is not None
        )

        if modern_provided and legacy_provided:
            raise ValueError(
                "Cannot mix modern (strategy_config, evaluator_config) and "
                "legacy (evaluator_class_name, selection_strategy) configuration. "
                "Please use one approach consistently."
            )

        if not modern_provided and not legacy_provided:
            raise ValueError(
                "Must provide either modern configs (strategy_config, evaluator_config) "
                "or legacy configs (evaluator_class_name, selection_strategy)"
            )

        # If legacy, ensure required fields are present
        if legacy_provided:
            if self.evaluator_class_name is None:
                raise ValueError("evaluator_class_name is required for legacy config")
            if self.evaluator_arg is None:
                raise ValueError("evaluator_arg is required for legacy config")

        # Validate reaction configuration
        has_simple_smarts = self.reaction_smarts is not None
        has_advanced_config = self.reaction_sequence_config is not None

        if not has_simple_smarts and not has_advanced_config:
            raise ValueError(
                "Must provide either reaction_smarts or reaction_sequence_config"
            )

        # Warn if both are provided (advanced takes precedence)
        if has_simple_smarts and has_advanced_config:
            logger = logging.getLogger(__name__)
            logger.warning(
                "Both reaction_smarts and reaction_sequence_config provided. "
                "Using reaction_sequence_config (advanced mode)."
            )

        # Set default warmup if not provided
        if self.warmup_config is None and modern_provided:
            self.warmup_config = BalancedWarmupConfig()


class RandomBaselineConfig(BaseModel):
    """Configuration for random baseline sampling"""

    evaluator_class_name: str
    evaluator_arg: str | dict
    reaction_smarts: str
    reagent_file_list: list[str]
    num_trials: int
    num_to_save: int
    ascending_output: bool = False
    outfile_name: Optional[str] = None
    log_filename: Optional[str] = None