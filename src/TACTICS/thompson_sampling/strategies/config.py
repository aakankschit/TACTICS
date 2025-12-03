"""Pydantic configuration models for selection strategies."""

from pydantic import BaseModel, Field, field_validator
from typing import Literal, Tuple


class GreedyConfig(BaseModel):
    """Configuration for Greedy selection strategy."""

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


class BoltzmannConfig(BaseModel):
    """Configuration for Boltzmann/Softmax selection."""

    strategy_type: Literal["boltzmann"] = "boltzmann"
    mode: Literal["maximize_boltzmann", "minimize_boltzmann"] = "maximize_boltzmann"
    temperature: float = Field(default=1.0, gt=0, description="Temperature parameter (lower = more exploitation)")


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
