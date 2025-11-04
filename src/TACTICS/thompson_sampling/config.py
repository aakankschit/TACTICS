from pydantic import BaseModel, Field, field_validator
from typing import Literal, Optional, Union
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
    StratifiedWarmupConfig,
    LatinHypercubeWarmupConfig
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
    StandardWarmupConfig,
    EnhancedWarmupConfig,
    StratifiedWarmupConfig,
    LatinHypercubeWarmupConfig
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

    This config supports two usage patterns:
    1. Legacy: Specify evaluator_class_name, evaluator_arg, selection_strategy (string)
    2. Modern: Use nested configs (strategy_config, warmup_config, evaluator_config)

    The modern approach with nested configs provides better validation and type safety.
    """

    # Core parameters
    reaction_smarts: str
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

        # Set default warmup if not provided
        if self.warmup_config is None and modern_provided:
            self.warmup_config = StratifiedWarmupConfig()


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