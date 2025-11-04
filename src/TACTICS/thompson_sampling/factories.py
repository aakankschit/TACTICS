"""Factory functions for creating Thompson Sampling components from configs."""

from typing import Union
import json

from .strategies.config import (
    GreedyConfig,
    RouletteWheelConfig,
    UCBConfig,
    EpsilonGreedyConfig,
    BoltzmannConfig,
    BayesUCBConfig
)
from .strategies import (
    GreedySelection,
    RouletteWheelSelection,
    UCBSelection,
    EpsilonGreedySelection,
    BayesUCBSelection
)
from .strategies.base_strategy import SelectionStrategy

from .warmup.config import (
    StandardWarmupConfig,
    EnhancedWarmupConfig,
    StratifiedWarmupConfig,
    LatinHypercubeWarmupConfig
)
from .warmup import (
    StandardWarmup,
    EnhancedWarmup,
    StratifiedWarmup,
    LatinHypercubeWarmup
)
from .warmup.base import WarmupStrategy

from .core.evaluator_config import (
    LookupEvaluatorConfig,
    DBEvaluatorConfig,
    FPEvaluatorConfig,
    MWEvaluatorConfig,
    ROCSEvaluatorConfig,
    FredEvaluatorConfig,
    MLClassifierEvaluatorConfig
)
from .core.evaluators import (
    LookupEvaluator,
    DBEvaluator,
    FPEvaluator,
    MWEvaluator,
    ROCSEvaluator,
    FredEvaluator,
    MLClassifierEvaluator,
    Evaluator
)


# Type aliases for config unions
StrategyConfig = Union[
    GreedyConfig,
    RouletteWheelConfig,
    UCBConfig,
    EpsilonGreedyConfig,
    BoltzmannConfig,
    BayesUCBConfig
]

WarmupConfig = Union[
    StandardWarmupConfig,
    EnhancedWarmupConfig,
    StratifiedWarmupConfig,
    LatinHypercubeWarmupConfig
]

EvaluatorConfig = Union[
    LookupEvaluatorConfig,
    DBEvaluatorConfig,
    FPEvaluatorConfig,
    MWEvaluatorConfig,
    ROCSEvaluatorConfig,
    FredEvaluatorConfig,
    MLClassifierEvaluatorConfig
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
            scaling=config.scaling,
            alpha_increment=config.alpha_increment,
            beta_increment=config.beta_increment,
            efficiency_threshold=config.efficiency_threshold
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
            efficiency_threshold=config.efficiency_threshold,
            p_high_bounds=config.p_high_bounds,
            p_low_bounds=config.p_low_bounds,
            delta_high=config.delta_high,
            delta_low=config.delta_low
        )

    elif isinstance(config, BoltzmannConfig):
        # Note: BoltzmannSelection not yet implemented in strategies/
        # This would need to be created similar to other strategies
        raise NotImplementedError(
            "BoltzmannSelection strategy not yet implemented. "
            "Use RouletteWheelSelection as an alternative."
        )

    else:
        raise ValueError(f"Unknown strategy config type: {type(config)}")


def create_warmup(config: WarmupConfig) -> WarmupStrategy:
    """
    Create a warmup strategy from a Pydantic config.

    Args:
        config: Warmup configuration (StandardWarmupConfig, StratifiedWarmupConfig, etc.)

    Returns:
        WarmupStrategy: Instantiated warmup strategy object

    Example:
        >>> config = StratifiedWarmupConfig()
        >>> warmup = create_warmup(config)
        >>> isinstance(warmup, StratifiedWarmup)
        True
    """
    if isinstance(config, StandardWarmupConfig):
        return StandardWarmup()

    elif isinstance(config, EnhancedWarmupConfig):
        return EnhancedWarmup()

    elif isinstance(config, StratifiedWarmupConfig):
        return StratifiedWarmup()

    elif isinstance(config, LatinHypercubeWarmupConfig):
        return LatinHypercubeWarmup()

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
        # LookupEvaluator expects a JSON string or dict
        input_dict = {
            "ref_filename": config.ref_filename,
            "ref_colname": config.ref_colname
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

    else:
        raise ValueError(f"Unknown evaluator config type: {type(config)}")
