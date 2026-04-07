"""Factory functions for creating Thompson Sampling components from configs."""

from typing import Union
import json

from .strategies.config import (
    GreedyConfig,
    RouletteWheelConfig,
    UCBConfig,
    EpsilonGreedyConfig,
    BayesUCBConfig,
    TopTwoConfig,
)
from .strategies import (
    GreedySelection,
    RouletteWheelSelection,
    UCBSelection,
    EpsilonGreedySelection,
    BayesUCBSelection,
    TopTwoSelection,
)
from .strategies.base_strategy import SelectionStrategy

from .warmup.config import (
    StandardWarmupConfig,
    EnhancedWarmupConfig,
    BalancedWarmupConfig
)
from .warmup import (
    StandardWarmup,
    EnhancedWarmup,
    BalancedWarmup
)
from .warmup.base import WarmupStrategy

from .core.evaluator_config import (
    LookupEvaluatorConfig,
    DBEvaluatorConfig,
    FPEvaluatorConfig,
    MWEvaluatorConfig,
    ROCSEvaluatorConfig,
    FredEvaluatorConfig,
    MLClassifierEvaluatorConfig,
    CustomEvaluatorConfig
)
from .core.evaluators import (
    LookupEvaluator,
    DBEvaluator,
    FPEvaluator,
    MWEvaluator,
    ROCSEvaluator,
    FredEvaluator,
    MLClassifierEvaluator,
    Evaluator,
    CustomEvaluator
)


# Type aliases for config unions
StrategyConfig = Union[
    GreedyConfig,
    RouletteWheelConfig,
    UCBConfig,
    EpsilonGreedyConfig,
    BayesUCBConfig,
    TopTwoConfig,
]

WarmupConfig = Union[
    StandardWarmupConfig,
    EnhancedWarmupConfig,
    BalancedWarmupConfig
]

EvaluatorConfig = Union[
    LookupEvaluatorConfig,
    DBEvaluatorConfig,
    FPEvaluatorConfig,
    MWEvaluatorConfig,
    ROCSEvaluatorConfig,
    FredEvaluatorConfig,
    MLClassifierEvaluatorConfig,
    CustomEvaluatorConfig
]


def create_strategy(config: StrategyConfig) -> SelectionStrategy:
    """
    Create a selection strategy from a Pydantic config.

    Args:
        config: Strategy configuration (GreedyConfig, RouletteWheelConfig, etc.)

    Returns:
        SelectionStrategy: Instantiated strategy object

    Example:
        >>> config = RouletteWheelConfig(mode="maximize", alpha=0.1, beta=0.1)
        >>> strategy = create_strategy(config)
        >>> isinstance(strategy, RouletteWheelSelection)
        True
    """
    if isinstance(config, GreedyConfig):
        return GreedySelection(mode=config.mode)

    elif isinstance(config, RouletteWheelConfig):
        return RouletteWheelSelection(
            mode=config.mode,
            alpha=config.alpha,
            beta=config.beta,
            exploration_phase_end=config.exploration_phase_end,
            transition_phase_end=config.transition_phase_end,
            min_observations=config.min_observations,
            adaptive_temperature=config.adaptive_temperature,
            alpha_increment=config.alpha_increment,
            beta_increment=config.beta_increment,
            efficiency_threshold=config.efficiency_threshold,
            alpha_max=config.alpha_max,
            cats_exploration_fraction=config.cats_exploration_fraction,
            cats_range=config.cats_range,
            divergence_threshold=config.divergence_threshold,
            cats_ema_decay=config.cats_ema_decay,
        )

    elif isinstance(config, UCBConfig):
        return UCBSelection(mode=config.mode, c=config.c)

    elif isinstance(config, EpsilonGreedyConfig):
        return EpsilonGreedySelection(
            mode=config.mode,
            epsilon=config.epsilon,
            decay=config.decay
        )

    elif isinstance(config, BayesUCBConfig):
        return BayesUCBSelection(
            mode=config.mode,
            initial_p_high=config.initial_p_high,
            initial_p_low=config.initial_p_low,
            exploration_phase_end=config.exploration_phase_end,
            transition_phase_end=config.transition_phase_end,
            min_observations=config.min_observations,
            cats_exploration_fraction=config.cats_exploration_fraction,
            criticality_metric=config.criticality_metric,
            n_adaptive_sharpening=config.n_adaptive_sharpening,
        )

    elif isinstance(config, TopTwoConfig):
        return TopTwoSelection(
            mode=config.mode,
            beta=config.beta,
            heated_scale=config.heated_scale,
            cooled_scale=config.cooled_scale,
            min_observations=config.min_observations,
            adaptive_temperature=config.adaptive_temperature,
            scale_increment=config.scale_increment,
            cooled_scale_increment=config.cooled_scale_increment,
            efficiency_threshold=config.efficiency_threshold,
            heated_scale_max=config.heated_scale_max,
            adaptive_disagreement=config.adaptive_disagreement,
            disagreement_window=config.disagreement_window,
            disagreement_high_threshold=config.disagreement_high_threshold,
            disagreement_low_threshold=config.disagreement_low_threshold,
            disagreement_decay_rate=config.disagreement_decay_rate,
            ema_alpha=config.ema_alpha,
            heated_scale_min=config.heated_scale_min,
            gmic_convergence_gate=config.gmic_convergence_gate,
            max_growth_per_step=config.max_growth_per_step,
        )

    else:
        raise ValueError(f"Unknown strategy config type: {type(config)}")


def create_warmup(config: WarmupConfig) -> WarmupStrategy:
    """
    Create a warmup strategy from a Pydantic config.

    Args:
        config: Warmup configuration (BalancedWarmupConfig, StandardWarmupConfig, etc.)

    Returns:
        WarmupStrategy: Instantiated warmup strategy object

    Example:
        >>> config = BalancedWarmupConfig(observations_per_reagent=5)
        >>> warmup = create_warmup(config)
        >>> isinstance(warmup, BalancedWarmup)
        True
    """
    if isinstance(config, StandardWarmupConfig):
        return StandardWarmup()

    elif isinstance(config, EnhancedWarmupConfig):
        return EnhancedWarmup()

    elif isinstance(config, BalancedWarmupConfig):
        return BalancedWarmup(
            observations_per_reagent=config.observations_per_reagent,
            seed=config.seed,
            use_per_reagent_variance=config.use_per_reagent_variance,
            shrinkage_strength=config.shrinkage_strength
        )

    else:
        raise ValueError(f"Unknown warmup config type: {type(config)}")


def create_evaluator(config: EvaluatorConfig) -> Evaluator:
    """
    Create an evaluator from a Pydantic config.

    Args:
        config: Evaluator configuration (LookupEvaluatorConfig, DBEvaluatorConfig, etc.)

    Returns:
        Evaluator: Instantiated evaluator object

    Example:
        >>> config = LookupEvaluatorConfig(
        ...     ref_filename="scores.csv",
        ...     ref_colname="Score"
        ... )
        >>> evaluator = create_evaluator(config)
        >>> isinstance(evaluator, LookupEvaluator)
        True
    """
    if isinstance(config, LookupEvaluatorConfig):
        # LookupEvaluator expects a dict with ref_filename, compound_col, and score_col
        input_dict = {
            "ref_filename": config.ref_filename,
            "compound_col": config.compound_col,
            "score_col": config.score_col,
        }
        return LookupEvaluator(json.dumps(input_dict))

    elif isinstance(config, DBEvaluatorConfig):
        input_dict = {
            "db_filename": config.db_filename,
            "db_prefix": config.db_prefix
        }
        return DBEvaluator(input_dict)

    elif isinstance(config, FPEvaluatorConfig):
        input_dict = {"query_smiles": config.query_smiles}
        return FPEvaluator(input_dict)

    elif isinstance(config, MWEvaluatorConfig):
        return MWEvaluator()

    elif isinstance(config, ROCSEvaluatorConfig):
        input_dict = {"query_molfile": config.query_molfile}
        evaluator = ROCSEvaluator(input_dict)
        evaluator.set_max_confs(config.max_confs)
        return evaluator

    elif isinstance(config, FredEvaluatorConfig):
        input_dict = {"design_unit_file": config.design_unit_file}
        evaluator = FredEvaluator(input_dict)
        evaluator.set_max_confs(config.max_confs)
        return evaluator

    elif isinstance(config, MLClassifierEvaluatorConfig):
        input_dict = {"model_filename": config.model_filename}
        return MLClassifierEvaluator(input_dict)

    elif isinstance(config, CustomEvaluatorConfig):
        scoring_function = config.scoring_function
        return CustomEvaluator(scoring_function)

    else:
        raise ValueError(f"Unknown evaluator config type: {type(config)}")
