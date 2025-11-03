"""Pydantic configuration models for selection strategies."""

from pydantic import BaseModel, Field
from typing import Literal


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
