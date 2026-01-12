"""Configuration presets for common Thompson Sampling use cases."""

from typing import Literal, Optional, TYPE_CHECKING, Union, List
from pathlib import Path
from .config import ThompsonSamplingConfig
from .strategies.config import (
    GreedyConfig,
    RouletteWheelConfig,
    UCBConfig,
    EpsilonGreedyConfig,
)
from .warmup.config import (
    StandardWarmupConfig,
    BalancedWarmupConfig,
    EnhancedWarmupConfig,
)
from .core.evaluator_config import LookupEvaluatorConfig, DBEvaluatorConfig

if TYPE_CHECKING:
    from ..library_enumeration import SynthesisPipeline


class ConfigPresets:
    """
    Pre-configured Thompson Sampling configurations for common use cases.

    These presets provide sensible defaults for different scenarios,
    reducing configuration complexity while maintaining flexibility.

    All presets support both maximize and minimize modes via the `mode` parameter.
    """

    @staticmethod
    def fast_exploration(
        synthesis_pipeline: "SynthesisPipeline",
        evaluator_config,
        num_iterations: int = 1000,
        mode: Literal["maximize", "minimize"] = "maximize",
        output_dir: Optional[str] = None,
    ) -> ThompsonSamplingConfig:
        """
        Fast exploration with epsilon-greedy strategy.

        Best for:
        - Quick initial screening
        - Balanced exploration/exploitation
        - Fast evaluators (LookupEvaluator, DBEvaluator)

        Configuration:
        - Strategy: Epsilon-Greedy (ε=0.2, decay=0.995)
        - Warmup: Balanced (K=3 observations per reagent, per-reagent variance)
        - Batch: Single mode (batch_size=1)

        Args:
            synthesis_pipeline: SynthesisPipeline with reaction config and reagent files
            evaluator_config: Evaluator configuration
            num_iterations: Number of Thompson sampling iterations
            mode: "maximize" for highest scores, "minimize" for lowest scores (e.g., docking)
            output_dir: Directory to save output files (optional). If specified, results and logs
                       will be saved to this directory with preset-specific names.
        """
        # Set output filenames if output_dir is specified
        results_filename = None
        log_filename = None
        if output_dir is not None:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            results_filename = str(output_path / "fast_exploration_results.csv")
            log_filename = str(output_path / "fast_exploration.log")

        return ThompsonSamplingConfig(
            synthesis_pipeline=synthesis_pipeline,
            num_ts_iterations=num_iterations,
            num_warmup_trials=3,
            strategy_config=EpsilonGreedyConfig(mode=mode, epsilon=0.2, decay=0.995),
            warmup_config=BalancedWarmupConfig(observations_per_reagent=3),
            evaluator_config=evaluator_config,
            batch_size=1,
            results_filename=results_filename,
            log_filename=log_filename,
        )

    @staticmethod
    def parallel_batch(
        synthesis_pipeline: "SynthesisPipeline",
        evaluator_config,
        num_iterations: int = 1000,
        batch_size: int = 100,
        mode: Literal["maximize", "minimize"] = "maximize",
        output_dir: Optional[str] = None,
    ) -> ThompsonSamplingConfig:
        """
        Parallel batch processing for computationally expensive evaluators.

        Best for:
        - Slow evaluators (ROCS, Fred, ML models, docking)
        - Large-scale screening
        - Multi-core systems

        Configuration:
        - Strategy: Roulette Wheel (adaptive thermal cycling)
        - Warmup: Balanced (K=5 observations per reagent, per-reagent variance)
        - Batch: Batch mode (configurable batch_size)

        Args:
            synthesis_pipeline: SynthesisPipeline with reaction config and reagent files
            evaluator_config: Evaluator configuration
            num_iterations: Number of Thompson sampling iterations
            batch_size: Number of compounds to sample per cycle
            mode: "maximize" for highest scores, "minimize" for lowest scores (e.g., docking)
            output_dir: Directory to save output files (optional). If specified, results and logs
                       will be saved to this directory with preset-specific names.
        """
        # Set output filenames if output_dir is specified
        results_filename = None
        log_filename = None
        if output_dir is not None:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            results_filename = str(output_path / "parallel_batch_results.csv")
            log_filename = str(output_path / "parallel_batch.log")

        return ThompsonSamplingConfig(
            synthesis_pipeline=synthesis_pipeline,
            num_ts_iterations=num_iterations,
            num_warmup_trials=5,
            strategy_config=RouletteWheelConfig(
                mode=mode,
                alpha=0.1,
                beta=0.1,
            ),
            warmup_config=BalancedWarmupConfig(observations_per_reagent=5),
            evaluator_config=evaluator_config,
            batch_size=batch_size,
            max_resamples=1000,
            results_filename=results_filename,
            log_filename=log_filename,
        )

    @staticmethod
    def conservative_exploit(
        synthesis_pipeline: "SynthesisPipeline",
        evaluator_config,
        num_iterations: int = 1000,
        mode: Literal["maximize", "minimize"] = "maximize",
        output_dir: Optional[str] = None,
    ) -> ThompsonSamplingConfig:
        """
        Conservative exploitation with greedy strategy.

        Best for:
        - Focus on best-performing reagents
        - Hit optimization (find the absolute best)
        - Well-characterized chemical space

        Configuration:
        - Strategy: Greedy (pure exploitation)
        - Warmup: Balanced (K=5 observations per reagent, per-reagent variance)
        - Batch: Single mode (batch_size=1)

        Args:
            synthesis_pipeline: SynthesisPipeline with reaction config and reagent files
            evaluator_config: Evaluator configuration
            num_iterations: Number of Thompson sampling iterations
            mode: "maximize" for highest scores, "minimize" for lowest scores (e.g., docking)
            output_dir: Directory to save output files (optional). If specified, results and logs
                       will be saved to this directory with preset-specific names.
        """
        # Set output filenames if output_dir is specified
        results_filename = None
        log_filename = None
        if output_dir is not None:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            results_filename = str(output_path / "conservative_exploit_results.csv")
            log_filename = str(output_path / "conservative_exploit.log")

        return ThompsonSamplingConfig(
            synthesis_pipeline=synthesis_pipeline,
            num_ts_iterations=num_iterations,
            num_warmup_trials=5,
            strategy_config=GreedyConfig(mode=mode),
            warmup_config=BalancedWarmupConfig(observations_per_reagent=5),
            evaluator_config=evaluator_config,
            batch_size=1,
            results_filename=results_filename,
            log_filename=log_filename,
        )

    @staticmethod
    def balanced_sampling(
        synthesis_pipeline: "SynthesisPipeline",
        evaluator_config,
        num_iterations: int = 1000,
        mode: Literal["maximize", "minimize"] = "maximize",
        output_dir: Optional[str] = None,
    ) -> ThompsonSamplingConfig:
        """
        Balanced exploration and exploitation with UCB strategy.

        Best for:
        - General-purpose screening
        - Theoretical guarantees
        - Diverse chemical space exploration

        Configuration:
        - Strategy: UCB (upper confidence bound, c=2.0)
        - Warmup: Balanced (K=3 observations per reagent, per-reagent variance)
        - Batch: Single mode (batch_size=1)

        Args:
            synthesis_pipeline: SynthesisPipeline with reaction config and reagent files
            evaluator_config: Evaluator configuration
            num_iterations: Number of Thompson sampling iterations
            mode: "maximize" for highest scores, "minimize" for lowest scores (e.g., docking)
            output_dir: Directory to save output files (optional). If specified, results and logs
                       will be saved to this directory with preset-specific names.
        """
        # Set output filenames if output_dir is specified
        results_filename = None
        log_filename = None
        if output_dir is not None:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            results_filename = str(output_path / "balanced_sampling_results.csv")
            log_filename = str(output_path / "balanced_sampling.log")

        return ThompsonSamplingConfig(
            synthesis_pipeline=synthesis_pipeline,
            num_ts_iterations=num_iterations,
            num_warmup_trials=3,
            strategy_config=UCBConfig(mode=mode, c=2.0),
            warmup_config=BalancedWarmupConfig(observations_per_reagent=3),
            evaluator_config=evaluator_config,
            batch_size=1,
            results_filename=results_filename,
            log_filename=log_filename,
        )

    @staticmethod
    def diverse_coverage(
        synthesis_pipeline: "SynthesisPipeline",
        evaluator_config,
        num_iterations: int = 1000,
        mode: Literal["maximize", "minimize"] = "maximize",
        output_dir: Optional[str] = None,
    ) -> ThompsonSamplingConfig:
        """
        Maximum diversity with roulette wheel and balanced warmup.

        Best for:
        - Reagents ordered by properties
        - Maximum chemical diversity is critical
        - Exploration-heavy applications

        Configuration:
        - Strategy: Roulette Wheel (high exploration)
        - Warmup: Balanced (K=5 observations per reagent, per-reagent variance)
        - Batch: Single mode (batch_size=1)

        Args:
            synthesis_pipeline: SynthesisPipeline with reaction config and reagent files
            evaluator_config: Evaluator configuration
            num_iterations: Number of Thompson sampling iterations
            mode: "maximize" for highest scores, "minimize" for lowest scores (e.g., docking)
            output_dir: Directory to save output files (optional). If specified, results and logs
                       will be saved to this directory with preset-specific names.
        """
        # Set output filenames if output_dir is specified
        results_filename = None
        log_filename = None
        if output_dir is not None:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            results_filename = str(output_path / "diverse_coverage_results.csv")
            log_filename = str(output_path / "diverse_coverage.log")

        return ThompsonSamplingConfig(
            synthesis_pipeline=synthesis_pipeline,
            num_ts_iterations=num_iterations,
            num_warmup_trials=5,
            strategy_config=RouletteWheelConfig(
                mode=mode,
                alpha=0.2,  # Higher alpha = more exploration
                beta=0.2,
            ),
            warmup_config=BalancedWarmupConfig(observations_per_reagent=5),
            evaluator_config=evaluator_config,
            batch_size=1,
            results_filename=results_filename,
            log_filename=log_filename,
        )

    @staticmethod
    def legacy_rws_maximize(
        synthesis_pipeline: "SynthesisPipeline",
        evaluator_config,
        num_iterations: int = 18500,
        max_resamples: int = 6000,
        output_dir: Optional[str] = None,
    ) -> ThompsonSamplingConfig:
        """
        Legacy RWS with Boltzmann weighting for maximize mode.

        This preset replicates the original RWS algorithm with:
        - Boltzmann-weighted Bayesian updates (exponential weighting of good scores)
        - Enhanced warmup (stochastic parallel pairing)
        - Automatic score scaling for maximize mode
        - Adaptive thermal cycling (α=0.1, β=0.1)

        Best for:
        - Replicating published RWS results
        - Problems with clear structure (ROCS similarity, etc.)
        - Maximum performance on structured optimization

        Expected performance: ~76±5 recovery (matching original RWS paper)

        Args:
            synthesis_pipeline: SynthesisPipeline with reaction config and reagent files
            evaluator_config: Evaluator configuration
            num_iterations: Number of iterations (default: 18500, matching paper)
            max_resamples: Early stopping after consecutive duplicates (default: 6000)
            output_dir: Directory to save output files (optional)
        """
        # Set output filenames if output_dir is specified
        results_filename = None
        log_filename = None
        if output_dir is not None:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            results_filename = str(output_path / "legacy_rws_maximize_results.csv")
            log_filename = str(output_path / "legacy_rws_maximize.log")

        return ThompsonSamplingConfig(
            synthesis_pipeline=synthesis_pipeline,
            num_ts_iterations=num_iterations,
            num_warmup_trials=5,
            strategy_config=RouletteWheelConfig(mode="maximize", alpha=0.1, beta=0.1),
            warmup_config=EnhancedWarmupConfig(),  # Stochastic parallel pairing
            evaluator_config=evaluator_config,
            batch_size=1,
            max_resamples=max_resamples,
            use_boltzmann_weighting=True,  # Enable legacy RWS algorithm
            results_filename=results_filename,
            log_filename=log_filename,
        )

    @staticmethod
    def legacy_rws_minimize(
        synthesis_pipeline: "SynthesisPipeline",
        evaluator_config,
        num_iterations: int = 18500,
        max_resamples: int = 6000,
        output_dir: Optional[str] = None,
    ) -> ThompsonSamplingConfig:
        """
        Legacy RWS with Boltzmann weighting for minimize mode (e.g., docking).

        This preset replicates the original RWS algorithm with:
        - Boltzmann-weighted Bayesian updates (exponential weighting of good scores)
        - Enhanced warmup (stochastic parallel pairing)
        - Automatic score scaling for minimize mode (inverts scores)
        - Adaptive thermal cycling (α=0.1, β=0.1)

        Best for:
        - Molecular docking (minimize binding energy)
        - Free energy calculations (minimize ΔG)
        - Any minimize-mode structured optimization

        Expected performance: ~76±5 recovery (matching original RWS paper)

        Args:
            synthesis_pipeline: SynthesisPipeline with reaction config and reagent files
            evaluator_config: Evaluator configuration (e.g., FredEvaluatorConfig)
            num_iterations: Number of iterations (default: 18500, matching paper)
            max_resamples: Early stopping after consecutive duplicates (default: 6000)
            output_dir: Directory to save output files (optional)
        """
        # Set output filenames if output_dir is specified
        results_filename = None
        log_filename = None
        if output_dir is not None:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            results_filename = str(output_path / "legacy_rws_minimize_results.csv")
            log_filename = str(output_path / "legacy_rws_minimize.log")

        return ThompsonSamplingConfig(
            synthesis_pipeline=synthesis_pipeline,
            num_ts_iterations=num_iterations,
            num_warmup_trials=5,
            strategy_config=RouletteWheelConfig(mode="minimize", alpha=0.1, beta=0.1),
            warmup_config=EnhancedWarmupConfig(),  # Stochastic parallel pairing
            evaluator_config=evaluator_config,
            batch_size=1,
            max_resamples=max_resamples,
            use_boltzmann_weighting=True,  # Enable legacy RWS algorithm
            results_filename=results_filename,
            log_filename=log_filename,
        )


# Convenience function for quick preset access
def get_preset(
    preset_name: str,
    synthesis_pipeline: "SynthesisPipeline",
    evaluator_config,
    **kwargs,
) -> ThompsonSamplingConfig:
    """
    Get a configuration preset by name.

    Args:
        preset_name: Name of preset ("fast_exploration", "parallel_batch", etc.)
        synthesis_pipeline: SynthesisPipeline with reaction config and reagent files
        evaluator_config: Evaluator configuration
        **kwargs: Additional arguments passed to preset function, including:
            - mode: "maximize" or "minimize" (controls optimization direction)
            - num_iterations: Number of iterations
            - batch_size: (parallel_batch only)
            - output_dir: Directory to save results and logs (all presets)

    Returns:
        ThompsonSamplingConfig: Configured preset

    Example:
        >>> from TACTICS.library_enumeration import SynthesisPipeline
        >>> from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef
        >>> from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig
        >>>
        >>> # Create pipeline
        >>> config = ReactionConfig(
        ...     reactions=[ReactionDef(reaction_smarts="...", step_index=0)],
        ...     reagent_file_list=["acids.smi", "amines.smi"]
        ... )
        >>> pipeline = SynthesisPipeline(config)
        >>>
        >>> # Get preset
        >>> evaluator = LookupEvaluatorConfig(ref_filename="scores.csv")
        >>> ts_config = get_preset(
        ...     "fast_exploration",
        ...     synthesis_pipeline=pipeline,
        ...     evaluator_config=evaluator
        ... )
    """
    presets = {
        "fast_exploration": ConfigPresets.fast_exploration,
        "parallel_batch": ConfigPresets.parallel_batch,
        "conservative_exploit": ConfigPresets.conservative_exploit,
        "balanced_sampling": ConfigPresets.balanced_sampling,
        "diverse_coverage": ConfigPresets.diverse_coverage,
        "legacy_rws_maximize": ConfigPresets.legacy_rws_maximize,
        "legacy_rws_minimize": ConfigPresets.legacy_rws_minimize,
    }

    if preset_name not in presets:
        raise ValueError(
            f"Unknown preset: {preset_name}. Available presets: {list(presets.keys())}"
        )

    return presets[preset_name](
        synthesis_pipeline=synthesis_pipeline,
        evaluator_config=evaluator_config,
        **kwargs,
    )
