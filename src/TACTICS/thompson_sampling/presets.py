"""Configuration presets for common Thompson Sampling use cases."""

from typing import Literal
from .config import ThompsonSamplingConfig
from .strategies.config import (
    GreedyConfig,
    RouletteWheelConfig,
    UCBConfig,
    EpsilonGreedyConfig
)
from .warmup.config import (
    StandardWarmupConfig,
    StratifiedWarmupConfig,
    LatinHypercubeWarmupConfig
)
from .core.evaluator_config import LookupEvaluatorConfig, DBEvaluatorConfig


class ConfigPresets:
    """
    Pre-configured Thompson Sampling configurations for common use cases.

    These presets provide sensible defaults for different scenarios,
    reducing configuration complexity while maintaining flexibility.

    All presets support both maximize and minimize modes via the `mode` parameter.
    """

    @staticmethod
    def fast_exploration(
        reaction_smarts: str,
        reagent_file_list: list[str],
        evaluator_config,
        num_iterations: int = 1000,
        mode: Literal["maximize", "minimize"] = "maximize"
    ) -> ThompsonSamplingConfig:
        """
        Fast exploration with epsilon-greedy strategy.

        Best for:
        - Quick initial screening
        - Balanced exploration/exploitation
        - Fast evaluators (LookupEvaluator, DBEvaluator)

        Configuration:
        - Strategy: Epsilon-Greedy (ε=0.2, decay=0.995)
        - Warmup: Stratified (balanced, no duplicates)
        - Batch: Single mode (batch_size=1)
        - Processes: 1 (sequential)

        Args:
            reaction_smarts: SMARTS string for reaction
            reagent_file_list: List of reagent file paths
            evaluator_config: Evaluator configuration
            num_iterations: Number of Thompson sampling iterations
            mode: "maximize" for highest scores, "minimize" for lowest scores (e.g., docking)
        """
        return ThompsonSamplingConfig(
            reaction_smarts=reaction_smarts,
            reagent_file_list=reagent_file_list,
            num_ts_iterations=num_iterations,
            num_warmup_trials=3,
            strategy_config=EpsilonGreedyConfig(
                mode=mode,
                epsilon=0.2,
                decay=0.995
            ),
            warmup_config=StratifiedWarmupConfig(),
            evaluator_config=evaluator_config,
            batch_size=1,
            processes=1,
            min_cpds_per_core=10
        )

    @staticmethod
    def parallel_batch(
        reaction_smarts: str,
        reagent_file_list: list[str],
        evaluator_config,
        num_iterations: int = 1000,
        batch_size: int = 100,
        processes: int = 4,
        mode: Literal["maximize", "minimize"] = "maximize"
    ) -> ThompsonSamplingConfig:
        """
        Parallel batch processing for computationally expensive evaluators.

        Best for:
        - Slow evaluators (ROCS, Fred, ML models, docking)
        - Large-scale screening
        - Multi-core systems

        Configuration:
        - Strategy: Roulette Wheel (adaptive thermal cycling)
        - Warmup: Latin Hypercube (optimal space coverage)
        - Batch: Batch mode (configurable batch_size)
        - Processes: Configurable (default=4)

        Args:
            reaction_smarts: SMARTS string for reaction
            reagent_file_list: List of reagent file paths
            evaluator_config: Evaluator configuration
            num_iterations: Number of Thompson sampling iterations
            batch_size: Number of compounds to sample per cycle
            processes: Number of CPU cores for parallel evaluation
            mode: "maximize" for highest scores, "minimize" for lowest scores (e.g., docking)
        """
        return ThompsonSamplingConfig(
            reaction_smarts=reaction_smarts,
            reagent_file_list=reagent_file_list,
            num_ts_iterations=num_iterations,
            num_warmup_trials=5,
            strategy_config=RouletteWheelConfig(
                mode=mode,
                alpha=0.1,
                beta=0.1,
                scaling=1.0,
                alpha_increment=0.01,
                beta_increment=0.001,
                efficiency_threshold=0.1
            ),
            warmup_config=LatinHypercubeWarmupConfig(),
            evaluator_config=evaluator_config,
            batch_size=batch_size,
            max_resamples=1000,
            processes=processes,
            min_cpds_per_core=10
        )

    @staticmethod
    def conservative_exploit(
        reaction_smarts: str,
        reagent_file_list: list[str],
        evaluator_config,
        num_iterations: int = 1000,
        mode: Literal["maximize", "minimize"] = "maximize"
    ) -> ThompsonSamplingConfig:
        """
        Conservative exploitation with greedy strategy.

        Best for:
        - Focus on best-performing reagents
        - Hit optimization (find the absolute best)
        - Well-characterized chemical space

        Configuration:
        - Strategy: Greedy (pure exploitation)
        - Warmup: Stratified (balanced initialization)
        - Batch: Single mode (batch_size=1)
        - Processes: 1 (sequential)

        Args:
            reaction_smarts: SMARTS string for reaction
            reagent_file_list: List of reagent file paths
            evaluator_config: Evaluator configuration
            num_iterations: Number of Thompson sampling iterations
            mode: "maximize" for highest scores, "minimize" for lowest scores (e.g., docking)
        """
        return ThompsonSamplingConfig(
            reaction_smarts=reaction_smarts,
            reagent_file_list=reagent_file_list,
            num_ts_iterations=num_iterations,
            num_warmup_trials=5,
            strategy_config=GreedyConfig(mode=mode),
            warmup_config=StratifiedWarmupConfig(),
            evaluator_config=evaluator_config,
            batch_size=1,
            processes=1,
            min_cpds_per_core=10
        )

    @staticmethod
    def balanced_sampling(
        reaction_smarts: str,
        reagent_file_list: list[str],
        evaluator_config,
        num_iterations: int = 1000,
        mode: Literal["maximize", "minimize"] = "maximize"
    ) -> ThompsonSamplingConfig:
        """
        Balanced exploration and exploitation with UCB strategy.

        Best for:
        - General-purpose screening
        - Theoretical guarantees
        - Diverse chemical space exploration

        Configuration:
        - Strategy: UCB (upper confidence bound, c=2.0)
        - Warmup: Stratified (balanced initialization)
        - Batch: Single mode (batch_size=1)
        - Processes: 1 (sequential)

        Args:
            reaction_smarts: SMARTS string for reaction
            reagent_file_list: List of reagent file paths
            evaluator_config: Evaluator configuration
            num_iterations: Number of Thompson sampling iterations
            mode: "maximize" for highest scores, "minimize" for lowest scores (e.g., docking)
        """
        return ThompsonSamplingConfig(
            reaction_smarts=reaction_smarts,
            reagent_file_list=reagent_file_list,
            num_ts_iterations=num_iterations,
            num_warmup_trials=3,
            strategy_config=UCBConfig(mode=mode, c=2.0),
            warmup_config=StratifiedWarmupConfig(),
            evaluator_config=evaluator_config,
            batch_size=1,
            processes=1,
            min_cpds_per_core=10
        )

    @staticmethod
    def diverse_coverage(
        reaction_smarts: str,
        reagent_file_list: list[str],
        evaluator_config,
        num_iterations: int = 1000,
        mode: Literal["maximize", "minimize"] = "maximize"
    ) -> ThompsonSamplingConfig:
        """
        Maximum diversity with roulette wheel and Latin Hypercube warmup.

        Best for:
        - Reagents ordered by properties
        - Maximum chemical diversity is critical
        - Exploration-heavy applications

        Configuration:
        - Strategy: Roulette Wheel (high exploration)
        - Warmup: Latin Hypercube (maximum space coverage)
        - Batch: Single mode (batch_size=1)
        - Processes: 1 (sequential)

        Args:
            reaction_smarts: SMARTS string for reaction
            reagent_file_list: List of reagent file paths
            evaluator_config: Evaluator configuration
            num_iterations: Number of Thompson sampling iterations
            mode: "maximize" for highest scores, "minimize" for lowest scores (e.g., docking)
        """
        return ThompsonSamplingConfig(
            reaction_smarts=reaction_smarts,
            reagent_file_list=reagent_file_list,
            num_ts_iterations=num_iterations,
            num_warmup_trials=5,
            strategy_config=RouletteWheelConfig(
                mode=mode,
                alpha=0.2,  # Higher alpha = more exploration
                beta=0.2,
                scaling=1.0,
                alpha_increment=0.02,
                beta_increment=0.002,
                efficiency_threshold=0.2
            ),
            warmup_config=LatinHypercubeWarmupConfig(),
            evaluator_config=evaluator_config,
            batch_size=1,
            processes=1,
            min_cpds_per_core=10
        )


# Convenience function for quick preset access
def get_preset(
    preset_name: str,
    reaction_smarts: str,
    reagent_file_list: list[str],
    evaluator_config,
    **kwargs
) -> ThompsonSamplingConfig:
    """
    Get a configuration preset by name.

    Args:
        preset_name: Name of preset ("fast_exploration", "parallel_batch", etc.)
        reaction_smarts: SMARTS string for reaction
        reagent_file_list: List of reagent file paths
        evaluator_config: Evaluator configuration
        **kwargs: Additional arguments passed to preset function, including:
            - mode: "maximize" or "minimize" (controls optimization direction)
            - num_iterations: Number of iterations
            - batch_size: (parallel_batch only)
            - processes: (parallel_batch only)

    Returns:
        ThompsonSamplingConfig: Configured preset

    Example:
        >>> from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig
        >>> evaluator = LookupEvaluatorConfig(ref_filename="scores.csv")
        >>>
        >>> # Maximize mode (default)
        >>> config = get_preset(
        ...     "fast_exploration",
        ...     reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        ...     reagent_file_list=["acids.smi", "amines.smi"],
        ...     evaluator_config=evaluator
        ... )
        >>>
        >>> # Minimize mode for docking
        >>> config = get_preset(
        ...     "parallel_batch",
        ...     reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        ...     reagent_file_list=["acids.smi", "amines.smi"],
        ...     evaluator_config=FredEvaluatorConfig(design_unit_file="receptor.oedu"),
        ...     mode="minimize",
        ...     batch_size=100,
        ...     processes=8
        ... )
    """
    presets = {
        "fast_exploration": ConfigPresets.fast_exploration,
        "parallel_batch": ConfigPresets.parallel_batch,
        "conservative_exploit": ConfigPresets.conservative_exploit,
        "balanced_sampling": ConfigPresets.balanced_sampling,
        "diverse_coverage": ConfigPresets.diverse_coverage,
    }

    if preset_name not in presets:
        raise ValueError(
            f"Unknown preset: {preset_name}. "
            f"Available presets: {list(presets.keys())}"
        )

    return presets[preset_name](
        reaction_smarts=reaction_smarts,
        reagent_file_list=reagent_file_list,
        evaluator_config=evaluator_config,
        **kwargs
    )
