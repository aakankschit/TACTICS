"""Pydantic configuration models for selection strategies."""

from pydantic import BaseModel, Field, field_validator
from typing import Literal, Optional, Tuple, Union, Annotated


class GreedyConfig(BaseModel):
    """Configuration for Greedy selection strategy.

    Pure argmax Thompson Sampling: sample posteriors, pick the best.
    """

    strategy_type: Literal["greedy"] = "greedy"
    mode: Literal["maximize", "minimize"] = "maximize"


class RouletteWheelConfig(BaseModel):
    """
    Configuration for Roulette Wheel selection with Component-Aware Thompson Sampling (CATS).

    Combines thermal cycling with component criticality analysis for efficient
    exploration of ultra-large combinatorial libraries. CATS automatically adjusts
    exploration based on Shannon entropy-based criticality.

    References:
        Zhao, H., Nittinger, E. & Tyrchan, C. Enhanced Thompson Sampling by Roulette
        Wheel Selection for Screening Ultra-Large Combinatorial Libraries.
        bioRxiv 2024.05.16.594622 (2024)
    """

    strategy_type: Literal["roulette_wheel"] = "roulette_wheel"
    mode: Literal["maximize", "minimize", "maximize_boltzmann", "minimize_boltzmann"] = "maximize"

    # Thermal cycling parameters
    alpha: float = Field(default=0.1, gt=0, description="Base temperature for heated component")
    beta: float = Field(default=0.05, gt=0, description="Base temperature for cooled components")

    # CATS parameters
    exploration_phase_end: float = Field(
        default=0.20,
        gt=0,
        le=1,
        description="Fraction of iterations before CATS starts (default: 0.20 = 20%)"
    )
    transition_phase_end: float = Field(
        default=0.60,
        gt=0,
        le=1,
        description="Fraction of iterations when CATS is fully applied (default: 0.60 = 60%)"
    )
    min_observations: int = Field(
        default=5,
        gt=0,
        description="Minimum observations per reagent before trusting criticality"
    )

    # CATS exploration decay
    cats_exploration_fraction: Optional[float] = Field(
        default=0.3,
        ge=0,
        le=1,
        description=(
            "Fraction of total cycles during which CATS explores at full strength. "
            "After this point, CATS influence decays linearly if criticality remains low. "
            "Set to None to disable decay. (default: 0.3 = first 30% of cycles)"
        ),
    )

    # CATS multiplier range override
    cats_range: Optional[float] = Field(
        default=None,
        gt=0,
        description=(
            "Override alpha/beta-derived CATS multiplier range. If set, "
            "cats_max = cats_range, cats_min = 1/cats_range. "
            "When None, range is derived from alpha/beta ratio."
        ),
    )

    # CATS EMA smoothing for relative GMIC
    cats_ema_decay: Optional[float] = Field(
        default=None,
        description=(
            "EMA decay factor for smoothing relative GMIC in the CATS multiplier. "
            "When set, relative_gmic is replaced by an EMA that accumulates "
            "directional signal across cycles. Smaller values = smoother. "
            "None disables smoothing (default, backward compatible). "
            "Must be in (0, 1) when set. Typical range: 0.05-0.15."
        ),
    )

    @field_validator('cats_ema_decay')
    @classmethod
    def validate_cats_ema_decay(cls, v):
        """Validate cats_ema_decay is in (0, 1) when set."""
        if v is not None and (v <= 0 or v >= 1):
            raise ValueError(f"cats_ema_decay must be in (0, 1), got {v}")
        return v

    # GMIC divergence gate
    divergence_threshold: float = Field(
        default=0.1,
        gt=0,
        description="KL divergence threshold for switching from diversity to GMIC criticality mode"
    )

    # Adaptive temperature parameters (legacy RWS-inspired)
    adaptive_temperature: bool = Field(
        default=False,
        description="Enable adaptive temperature control (increase alpha/beta when sampling efficiency drops)"
    )
    alpha_increment: float = Field(
        default=0.01, ge=0,
        description="Amount to increase alpha when efficiency drops below threshold"
    )
    beta_increment: float = Field(
        default=0.001, ge=0,
        description="Amount to increase beta when zero unique compounds found"
    )
    efficiency_threshold: float = Field(
        default=0.10, ge=0, le=1,
        description="Efficiency below which alpha is incremented"
    )
    alpha_max: float = Field(
        default=2.0, gt=0,
        description="Maximum alpha value"
    )

    @field_validator('transition_phase_end')
    @classmethod
    def validate_phase_progression(cls, v, info):
        """Ensure transition_phase_end > exploration_phase_end."""
        if 'exploration_phase_end' in info.data and v <= info.data['exploration_phase_end']:
            raise ValueError(
                f"transition_phase_end ({v}) must be > exploration_phase_end "
                f"({info.data['exploration_phase_end']})"
            )
        return v


class UCBConfig(BaseModel):
    """Configuration for Upper Confidence Bound selection."""

    strategy_type: Literal["ucb"] = "ucb"
    mode: Literal["maximize", "minimize"] = "maximize"
    c: float = Field(default=2.0, gt=0, description="Exploration parameter (higher = more exploration)")


class EpsilonGreedyConfig(BaseModel):
    """Configuration for Epsilon-Greedy selection with decaying epsilon."""

    strategy_type: Literal["epsilon_greedy"] = "epsilon_greedy"
    mode: Literal["maximize", "minimize"] = "maximize"
    epsilon: float = Field(default=0.1, ge=0, le=1, description="Initial exploration probability")
    decay: float = Field(default=0.995, gt=0, le=1, description="Decay rate for epsilon per iteration")


class TopTwoConfig(BaseModel):
    """Configuration for Top-Two Thompson Sampling (TT-TS).

    TT-TS targets best-arm identification rather than regret minimization.
    It draws two independent posterior samples and explores challengers when
    there is genuine uncertainty about which reagent is best.

    Supports asymmetric thermal cycling via posterior std scaling. Criticality
    is used strictly for weighted component rotation — it does NOT produce a
    temperature multiplier.

    Reference:
        Russo, D. (2020). Simple Bayesian Algorithms for Best-Arm Identification.
        Operations Research, 68(6), 1625-1647.
    """

    strategy_type: Literal["top_two"] = "top_two"
    mode: Literal["maximize", "minimize"] = "maximize"
    beta: float = Field(
        default=0.5,
        ge=0,
        le=1,
        description=(
            "Probability of selecting the challenger when two posterior samples "
            "disagree on the best reagent. β=0.5 gives equal weight to exploration "
            "and exploitation. Higher β → more exploration of uncertain reagents."
        ),
    )

    # Thermal cycling parameters
    heated_scale: float = Field(
        default=1.5,
        gt=0,
        description=(
            "Multiplier on posterior std for the heated component. "
            ">1 inflates uncertainty → more TT-TS disagreement → exploration. "
            "Set to 1.0 to disable thermal cycling."
        ),
    )
    cooled_scale: float = Field(
        default=0.75,
        gt=0,
        description=(
            "Multiplier on posterior std for cooled components. "
            "<1 deflates uncertainty → more TT-TS agreement → exploitation. "
            "Set to 1.0 to disable thermal cycling."
        ),
    )

    # GMIC criticality (for weighted rotation only — does NOT affect temperature)
    min_observations: int = Field(
        default=5,
        gt=0,
        description="Minimum observations per reagent before trusting GMIC criticality",
    )

    # Adaptive thermal cycling parameters
    adaptive_temperature: bool = Field(
        default=False,
        description=(
            "Enable adaptive thermal cycling. Progressively increases heated_scale "
            "and decreases cooled_scale when sampling efficiency drops, counteracting "
            "posterior tightening over long searches."
        ),
    )
    scale_increment: float = Field(
        default=0.01,
        ge=0,
        description="Amount to increase heated_scale when efficiency drops below threshold",
    )
    cooled_scale_increment: float = Field(
        default=0.001,
        ge=0,
        description="Amount to decrease cooled_scale when zero unique compounds found",
    )
    efficiency_threshold: float = Field(
        default=0.10,
        ge=0,
        le=1,
        description="Efficiency below which heated_scale is incremented",
    )
    heated_scale_max: float = Field(
        default=5.0,
        gt=0,
        description="Maximum heated_scale value",
    )

    # Bidirectional disagreement-rate adaptive thermal cycling
    adaptive_disagreement: bool = Field(
        default=True,
        description=(
            "Enable bidirectional disagreement-rate adaptation. Tracks per-component "
            "TT-TS disagreement rate via exponential moving average. Reduces "
            "heated_scale when disagreement is saturated (>high_threshold) and "
            "increases it when too low (<low_threshold). Ensures TT-TS works "
            "optimally on both balanced and imbalanced libraries."
        ),
    )
    disagreement_high_threshold: float = Field(
        default=0.8,
        gt=0,
        le=1,
        description="Disagreement rate above which heated_scale is reduced toward 1.0.",
    )
    disagreement_low_threshold: float = Field(
        default=0.3,
        gt=0,
        le=1,
        description="Disagreement rate below which heated_scale is increased.",
    )
    disagreement_decay_rate: float = Field(
        default=0.95,
        gt=0,
        lt=1,
        description=(
            "Multiplicative decay factor for heated_scale adaptation. "
            "Applied to the excess (heated_scale - 1.0) each adaptation step."
        ),
    )
    ema_alpha: float = Field(
        default=0.02,
        gt=0,
        lt=1,
        description=(
            "Exponential moving average smoothing factor for per-component "
            "disagreement tracking. Smaller = smoother (longer memory)."
        ),
    )
    heated_scale_min: float = Field(
        default=1.0,
        gt=0,
        description="Minimum heated_scale value (floor). Default 1.0 = raw posteriors.",
    )
    gmic_convergence_gate: Optional[float] = Field(
        default=None,
        description=(
            "GMIC threshold above which heated_scale inflation is suppressed. "
            "When a component's GMIC exceeds this value, the adaptive rule will "
            "not inflate its heated_scale even if disagreement is low — the "
            "component is already converged and low disagreement is correct. "
            "Set to None to disable (default). Typical values: 0.5–1.0."
        ),
    )
    max_growth_per_step: Optional[float] = Field(
        default=None,
        description=(
            "Maximum heated_scale increase per adaptation step. Caps the "
            "per-step growth when a component converges rapidly "
            "(disagreement drops below low_threshold). "
            "Set to None for uncapped growth (default). "
            "Disabled after diagnostic analysis showed the 'runaway' on "
            "balanced 3-comp libraries has no measurable recovery impact — "
            "the smaller TACTICS advantage on those libraries is due to "
            "easier SAR (higher min_GMIC), not mechanism malfunction."
        ),
    )
    # Legacy fields kept for backward compatibility of disagreement_window
    disagreement_window: int = Field(
        default=200,
        gt=1,
        description="Rolling window size for global disagreement rate (diagnostics).",
    )

    @field_validator("disagreement_low_threshold")
    @classmethod
    def validate_threshold_ordering(cls, v, info):
        high = info.data.get("disagreement_high_threshold")
        if high is not None and v >= high:
            raise ValueError(
                f"disagreement_low_threshold ({v}) must be < "
                f"disagreement_high_threshold ({high})"
            )
        return v


class BayesUCBConfig(BaseModel):
    """
    Configuration for Bayes-UCB selection with Component-Aware Thompson Sampling (CATS).

    Uses Bayesian Upper Confidence Bounds with Student-t quantiles and combines
    percentile-based thermal cycling with component criticality analysis for
    efficient exploration of ultra-large combinatorial libraries.

    The percentile parameters serve as an analog to temperature in thermal cycling:
    - Higher percentile → wider confidence bounds → more exploration
    - Lower percentile → tighter bounds → more exploitation

    CATS automatically adjusts exploration based on Shannon entropy-based criticality.

    References:
        Kaufmann, E., Cappé, O., & Garivier, A. (2012). On Bayesian upper confidence
        bounds for bandit problems. In AISTATS.
    """

    strategy_type: Literal["bayes_ucb"] = "bayes_ucb"
    mode: Literal["maximize", "minimize"] = "maximize"

    # Thermal cycling percentile parameters
    initial_p_high: float = Field(
        default=0.90,
        ge=0.5,
        le=0.999,
        description="Base percentile for heated component (more exploration)"
    )
    initial_p_low: float = Field(
        default=0.60,
        ge=0.5,
        le=0.999,
        description="Base percentile for cooled components (more exploitation)"
    )

    # CATS parameters
    exploration_phase_end: float = Field(
        default=0.20,
        gt=0,
        le=1,
        description="Fraction of iterations before CATS starts (default: 0.20 = 20%)"
    )
    transition_phase_end: float = Field(
        default=0.60,
        gt=0,
        le=1,
        description="Fraction of iterations when CATS is fully applied (default: 0.60 = 60%)"
    )
    min_observations: int = Field(
        default=5,
        gt=0,
        description="Minimum observations per reagent before trusting criticality"
    )

    # CATS exploration decay
    cats_exploration_fraction: Optional[float] = Field(
        default=0.3,
        ge=0,
        le=1,
        description=(
            "Fraction of total cycles during which CATS explores at full strength. "
            "After this point, CATS influence decays linearly if criticality remains low. "
            "Set to None to disable decay. (default: 0.3 = first 30% of cycles)"
        ),
    )

    # CATS criticality metric
    criticality_metric: Literal["ipr", "shannon"] = Field(
        default="ipr",
        description=(
            "Metric for computing component criticality. "
            "'ipr' uses Inverse Participation Ratio (sensitive to probability concentration). "
            "'shannon' uses Shannon entropy (legacy, insensitive at large N)."
        ),
    )
    n_adaptive_sharpening: bool = Field(
        default=True,
        description=(
            "Enable N-adaptive sharpening of z-scores before softmax. "
            "Counteracts softmax flattening for components with many reagents (large N). "
            "Only applies when criticality_metric='ipr'."
        ),
    )

    @field_validator('transition_phase_end')
    @classmethod
    def validate_phase_progression(cls, v, info):
        """Ensure transition_phase_end > exploration_phase_end."""
        if 'exploration_phase_end' in info.data and v <= info.data['exploration_phase_end']:
            raise ValueError(
                f"transition_phase_end ({v}) must be > exploration_phase_end "
                f"({info.data['exploration_phase_end']})"
            )
        return v
