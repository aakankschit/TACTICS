"""Configuration presets for common Thompson Sampling use cases.

Presets map directly to the method configurations validated in large-scale
benchmarking (114,450+ trials, 21 libraries). Users who call ``get_preset()``
with no arguments get the best-performing method.

Preset hierarchy:
    ``"recommended"``
        Enhanced warmup + Top-Two TS + Boltzmann weighting.
        Best overall method (86.1% mean top-100 recovery).
    ``"recommended_rws"``
        Enhanced warmup + Roulette Wheel Selection (CATS) + Boltzmann weighting.
        Close second (85.5%). The original TACTICS method.
    ``"baseline"``
        Balanced warmup + Greedy selection. Isolates the framework's warmup
        contribution (+1.5 pts on 2-component libraries via warmup alone).
    ``"legacy_rws"``
        Reproduces the original RWS algorithm from Zhao et al. (2025).
        Pass ``mode="minimize"`` for docking.
"""

from typing import Literal, Optional, TYPE_CHECKING
from pathlib import Path

from .config import ThompsonSamplingConfig
from .strategies.config import (
    GreedyConfig,
    RouletteWheelConfig,
    TopTwoConfig,
)
from .warmup.config import (
    BalancedWarmupConfig,
    EnhancedWarmupConfig,
)
from .core.evaluator_config import LookupEvaluatorConfig, DBEvaluatorConfig

if TYPE_CHECKING:
    from ..library_enumeration import SynthesisPipeline


def _output_paths(output_dir: Optional[str], prefix: str):
    """Helper to build results/log filenames from an output directory."""
    if output_dir is None:
        return None, None
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    return (
        str(output_path / f"{prefix}_results.csv"),
        str(output_path / f"{prefix}.log"),
    )


class ConfigPresets:
    """
    Pre-configured Thompson Sampling configurations validated by benchmarking.

    Each preset corresponds to a method configuration tested across 21
    combinatorial libraries (114,450+ trials). All presets support both
    ``"maximize"`` and ``"minimize"`` modes.

    Recommended (use for new work):
        - ``recommended``: Enhanced + TT-TS + Boltzmann (best overall)
        - ``recommended_rws``: Enhanced + RWS/CATS + Boltzmann (close second)

    Baseline:
        - ``baseline``: Balanced + Greedy (isolates warmup contribution)

    Legacy (for reproducing published results):
        - ``legacy_rws``: Original RWS algorithm (pass ``mode`` for direction)
    """

    @staticmethod
    def recommended(
        synthesis_pipeline: "SynthesisPipeline",
        evaluator_config,
        num_iterations: int = 1000,
        batch_size: int = 100,
        mode: Literal["maximize", "minimize"] = "maximize",
        output_dir: Optional[str] = None,
    ) -> ThompsonSamplingConfig:
        """
        Best-performing method from large-scale benchmarking.

        Uses Top-Two Thompson Sampling with Enhanced warmup, Boltzmann
        weighting, and GMIC-weighted component rotation. Adaptive
        per-component thermal cycling auto-tunes exploration on both
        balanced and imbalanced libraries.

        Performance: 86.1% mean top-100 recovery (114,450 trials, 21 libraries).

        Args:
            synthesis_pipeline: SynthesisPipeline with reaction config and reagent files
            evaluator_config: Evaluator configuration
            num_iterations: Number of Thompson sampling iterations
            batch_size: Number of compounds to sample per cycle (default: 100)
            mode: "maximize" for highest scores, "minimize" for lowest (e.g., docking)
            output_dir: Directory to save output files (optional)
        """
        results_filename, log_filename = _output_paths(output_dir, "recommended")
        return ThompsonSamplingConfig(
            synthesis_pipeline=synthesis_pipeline,
            num_ts_iterations=num_iterations,
            num_warmup_trials=5,
            strategy_config=TopTwoConfig(mode=mode),
            warmup_config=EnhancedWarmupConfig(),
            evaluator_config=evaluator_config,
            batch_size=batch_size,
            max_resamples=1000,
            use_boltzmann_weighting=True,
            results_filename=results_filename,
            log_filename=log_filename,
        )

    @staticmethod
    def recommended_rws(
        synthesis_pipeline: "SynthesisPipeline",
        evaluator_config,
        num_iterations: int = 1000,
        batch_size: int = 100,
        mode: Literal["maximize", "minimize"] = "maximize",
        output_dir: Optional[str] = None,
    ) -> ThompsonSamplingConfig:
        """
        Roulette Wheel Selection with Enhanced warmup and Boltzmann weighting.

        The original TACTICS method (CATS + GMIC-weighted rotation). Close
        second to TT-TS (85.5% vs 86.1% mean recovery). Preferred when
        thermal cycling behavior is well-understood for the target library.

        Performance: 85.5% mean top-100 recovery (114,450 trials, 21 libraries).

        Args:
            synthesis_pipeline: SynthesisPipeline with reaction config and reagent files
            evaluator_config: Evaluator configuration
            num_iterations: Number of Thompson sampling iterations
            batch_size: Number of compounds to sample per cycle (default: 100)
            mode: "maximize" for highest scores, "minimize" for lowest (e.g., docking)
            output_dir: Directory to save output files (optional)
        """
        results_filename, log_filename = _output_paths(output_dir, "recommended_rws")
        return ThompsonSamplingConfig(
            synthesis_pipeline=synthesis_pipeline,
            num_ts_iterations=num_iterations,
            num_warmup_trials=5,
            strategy_config=RouletteWheelConfig(mode=mode, alpha=0.1, beta=0.05),
            warmup_config=EnhancedWarmupConfig(),
            evaluator_config=evaluator_config,
            batch_size=batch_size,
            max_resamples=1000,
            use_boltzmann_weighting=True,
            results_filename=results_filename,
            log_filename=log_filename,
        )

    @staticmethod
    def baseline(
        synthesis_pipeline: "SynthesisPipeline",
        evaluator_config,
        num_iterations: int = 1000,
        mode: Literal["maximize", "minimize"] = "maximize",
        output_dir: Optional[str] = None,
    ) -> ThompsonSamplingConfig:
        """
        Baseline method: Balanced warmup + Greedy selection.

        Isolates the contribution of TACTICS' stratified warmup without any
        advanced selection strategy. Useful as a reference point — the gap
        between this baseline and ``recommended`` measures the value added
        by TT-TS and GMIC-weighted rotation.

        Performance: Balanced-Greedy provides +1.5 pts over Legacy-Greedy
        on 2-component libraries via warmup alone (significant on 10/11).

        Args:
            synthesis_pipeline: SynthesisPipeline with reaction config and reagent files
            evaluator_config: Evaluator configuration
            num_iterations: Number of Thompson sampling iterations
            mode: "maximize" for highest scores, "minimize" for lowest (e.g., docking)
            output_dir: Directory to save output files (optional)
        """
        results_filename, log_filename = _output_paths(output_dir, "baseline")
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
    def legacy_rws(
        synthesis_pipeline: "SynthesisPipeline",
        evaluator_config,
        num_iterations: int = 18500,
        max_resamples: int = 6000,
        mode: Literal["maximize", "minimize"] = "maximize",
        output_dir: Optional[str] = None,
    ) -> ThompsonSamplingConfig:
        """
        Legacy RWS (reproduces Zhao et al. 2025).

        Replicates the original Enhanced Thompson Sampling algorithm:
        - Enhanced warmup (stochastic parallel pairing)
        - Roulette Wheel Selection with Boltzmann thermal cycling
        - Unweighted round-robin component rotation (no GMIC)

        Use this only for reproducing published results. For new work,
        use ``"recommended"`` or ``"recommended_rws"`` instead.

        The ``mode`` parameter controls optimization direction:
        - ``"maximize"``: higher scores are better (e.g., ROCS similarity)
        - ``"minimize"``: lower scores are better (e.g., docking scores)

        Internally, minimize mode negates scores before Boltzmann weighting
        so that lower raw scores receive higher Boltzmann weights.

        Args:
            synthesis_pipeline: SynthesisPipeline with reaction config and reagent files
            evaluator_config: Evaluator configuration
            num_iterations: Number of iterations (default: 18500, matching paper)
            max_resamples: Early stopping after consecutive duplicates (default: 6000)
            mode: "maximize" or "minimize" (default: "maximize")
            output_dir: Directory to save output files (optional)
        """
        results_filename, log_filename = _output_paths(output_dir, "legacy_rws")
        return ThompsonSamplingConfig(
            synthesis_pipeline=synthesis_pipeline,
            num_ts_iterations=num_iterations,
            num_warmup_trials=5,
            strategy_config=RouletteWheelConfig(mode=mode, alpha=0.1, beta=0.1),
            warmup_config=EnhancedWarmupConfig(),
            evaluator_config=evaluator_config,
            batch_size=1,
            max_resamples=max_resamples,
            use_boltzmann_weighting=True,
            results_filename=results_filename,
            log_filename=log_filename,
        )


def get_preset(
    preset_name: str = "recommended",
    synthesis_pipeline: "SynthesisPipeline" = None,
    evaluator_config=None,
    **kwargs,
) -> ThompsonSamplingConfig:
    """
    Get a configuration preset by name.

    Args:
        preset_name: Name of preset (default: ``"recommended"``). Available:

            - ``"recommended"``: Enhanced + TT-TS + Boltzmann (best overall, 86.1%)
            - ``"recommended_rws"``: Enhanced + RWS/CATS + Boltzmann (85.5%)
            - ``"baseline"``: Balanced + Greedy (isolates warmup contribution)
            - ``"legacy_rws"``: Original RWS algorithm (reproduces Zhao et al. 2025)

        synthesis_pipeline: SynthesisPipeline with reaction config and reagent files
        evaluator_config: Evaluator configuration
        **kwargs: Additional arguments passed to preset function, including:
            - mode: "maximize" or "minimize"
            - num_iterations: Number of iterations
            - batch_size: Compounds per cycle (recommended/recommended_rws, default 100)
            - output_dir: Directory to save results and logs

    Returns:
        ThompsonSamplingConfig: Configured preset

    Example:
        >>> from TACTICS.library_enumeration import SynthesisPipeline
        >>> from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef
        >>> from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig
        >>>
        >>> pipeline = SynthesisPipeline(ReactionConfig(
        ...     reactions=[ReactionDef(reaction_smarts="...", step_index=0)],
        ...     reagent_file_list=["acids.smi", "amines.smi"]
        ... ))
        >>> evaluator = LookupEvaluatorConfig(ref_filename="scores.csv")
        >>> ts_config = get_preset(
        ...     synthesis_pipeline=pipeline,
        ...     evaluator_config=evaluator
        ... )
    """
    presets = {
        "recommended": ConfigPresets.recommended,
        "recommended_rws": ConfigPresets.recommended_rws,
        "baseline": ConfigPresets.baseline,
        "legacy_rws": ConfigPresets.legacy_rws,
    }

    if preset_name not in presets:
        raise ValueError(
            f"Unknown preset: {preset_name}. "
            f"Available presets: {list(presets.keys())}"
        )

    return presets[preset_name](
        synthesis_pipeline=synthesis_pipeline,
        evaluator_config=evaluator_config,
        **kwargs,
    )
