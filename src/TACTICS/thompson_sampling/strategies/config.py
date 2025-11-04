"""Pydantic configuration models for selection strategies."""

from pydantic import BaseModel, Field
from typing import Literal, Tuple


class GreedyConfig(BaseModel):
    """Configuration for Greedy selection strategy."""

    strategy_type: Literal["greedy"] = "greedy"
    mode: Literal["maximize", "minimize"] = "maximize"


class RouletteWheelConfig(BaseModel):
    """
    Configuration for Roulette Wheel selection with adaptive thermal cycling.

    Implements the enhanced Thompson sampling approach from:
    Zhao, H., Nittinger, E. & Tyrchan, C. Enhanced Thompson Sampling by Roulette
    Wheel Selection for Screening Ultra-Large Combinatorial Libraries.
    bioRxiv 2024.05.16.594622 (2024)
    """

    strategy_type: Literal["roulette_wheel"] = "roulette_wheel"
    mode: Literal["maximize", "minimize", "maximize_boltzmann", "minimize_boltzmann"] = "maximize"
    alpha: float = Field(default=0.1, gt=0, description="Initial temperature for heated component")
    beta: float = Field(default=0.1, gt=0, description="Initial temperature for cooled components")
    scaling: float = Field(default=1.0, gt=0, description="Scaling factor for scores")
    alpha_increment: float = Field(default=0.01, ge=0, description="Alpha increase when efficiency drops")
    beta_increment: float = Field(default=0.001, ge=0, description="Beta increase when no unique compounds")
    efficiency_threshold: float = Field(default=0.1, gt=0, le=1, description="Threshold for triggering temperature adjustment")


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
    Configuration for Adaptive Bayes-UCB selection with thermal cycling.

    This strategy uses Bayesian Upper Confidence Bounds with Student-t quantiles
    and adapts exploration/exploitation dynamically via percentile parameters.

    The percentile parameters serve as an analog to temperature in thermal cycling:
    - Higher percentile → wider confidence bounds → more exploration
    - Lower percentile → tighter bounds → more exploitation

    References:
        Kaufmann, E., Cappé, O., & Garivier, A. (2012). On Bayesian upper confidence
        bounds for bandit problems. In AISTATS.
    """

    strategy_type: Literal["bayes_ucb"] = "bayes_ucb"
    mode: Literal["maximize", "minimize"] = "maximize"

    # Initial percentile parameters
    initial_p_high: float = Field(
        default=0.90,
        ge=0.5,
        le=0.999,
        description="Initial percentile for heated component (more exploration)"
    )
    initial_p_low: float = Field(
        default=0.90,
        ge=0.5,
        le=0.999,
        description="Initial percentile for cooled components (more exploitation)"
    )

    # Adaptation parameters
    efficiency_threshold: float = Field(
        default=0.10,
        gt=0,
        lt=1,
        description="Sampling efficiency threshold (e.g., 0.1 = 10% unique compounds)"
    )

    # Percentile bounds
    p_high_bounds: Tuple[float, float] = Field(
        default=(0.85, 0.995),
        description="Lower and upper bounds for p_high percentile"
    )
    p_low_bounds: Tuple[float, float] = Field(
        default=(0.50, 0.90),
        description="Lower and upper bounds for p_low percentile"
    )

    # Update step sizes
    delta_high: float = Field(
        default=0.01,
        gt=0,
        description="Step size for increasing p_high when stuck"
    )
    delta_low: float = Field(
        default=0.005,
        gt=0,
        description="Step size for decreasing p_low when making progress"
    )
