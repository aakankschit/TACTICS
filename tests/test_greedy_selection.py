"""Tests for Greedy selection strategy."""

import pytest
import numpy as np

from TACTICS.thompson_sampling.strategies.greedy_selection import GreedySelection
from TACTICS.thompson_sampling.strategies.config import GreedyConfig
from TACTICS.thompson_sampling.factories import create_strategy
from TACTICS.thompson_sampling.core.reagent import Reagent


def _make_reagents(means, stds, n_samples_list):
    """Helper to create reagent lists with specified posteriors."""
    reagents = []
    for i, (m, s, n) in enumerate(zip(means, stds, n_samples_list)):
        r = Reagent(f"reagent_{i}", f"C{'C' * i}")
        r.mean = m
        r.std = s
        r.n_samples = n
        r.current_phase = "search"
        reagents.append(r)
    return reagents


class TestPlainGreedy:
    """GreedySelection is pure argmax/argmin Thompson Sampling."""

    def test_default_config(self):
        config = GreedyConfig()
        assert config.strategy_type == "greedy"
        assert config.mode == "maximize"

    def test_greedy_select_reagent(self):
        strategy = GreedySelection(mode="maximize")
        reagents = _make_reagents(
            means=[1.0, 5.0, 3.0],
            stds=[0.01, 0.01, 0.01],
            n_samples_list=[10, 10, 10],
        )
        rng = np.random.default_rng(42)
        # With tiny std, should almost always pick index 1 (mean=5.0)
        selections = [strategy.select_reagent(reagents, rng=rng) for _ in range(20)]
        assert all(s == 1 for s in selections)

    def test_greedy_get_component_criticality_returns_none(self):
        strategy = GreedySelection(mode="maximize")
        reagents = _make_reagents([1.0, 2.0], [0.1, 0.1], [10, 10])
        assert strategy.get_component_criticality(reagents) is None

    def test_greedy_get_component_state_returns_none(self):
        strategy = GreedySelection(mode="maximize")
        reagents = _make_reagents([1.0, 2.0], [0.1, 0.1], [10, 10])
        assert strategy.get_component_state(reagents, 0, 0, 100) is None

    def test_minimize_mode(self):
        strategy = GreedySelection(mode="minimize")
        reagents = _make_reagents(
            means=[5.0, 1.0, 3.0],
            stds=[0.01, 0.01, 0.01],
            n_samples_list=[10, 10, 10],
        )
        rng = np.random.default_rng(42)
        selections = [strategy.select_reagent(reagents, rng=rng) for _ in range(20)]
        assert all(s == 1 for s in selections)

    def test_disallow_mask(self):
        strategy = GreedySelection(mode="maximize")
        reagents = _make_reagents(
            means=[1.0, 100.0, 3.0],
            stds=[0.01, 0.01, 0.01],
            n_samples_list=[10, 10, 10],
        )
        rng = np.random.default_rng(42)
        # Disallow index 1 (the best) — should pick index 2
        selections = [
            strategy.select_reagent(reagents, disallow_mask={1}, rng=rng)
            for _ in range(20)
        ]
        assert all(s == 2 for s in selections)


class TestFactoryRoundtrip:
    """Factory correctly wires config to strategy."""

    def test_plain_greedy_factory(self):
        config = GreedyConfig()
        strategy = create_strategy(config)
        assert isinstance(strategy, GreedySelection)
        assert strategy.mode == "maximize"

    def test_minimize_greedy_factory(self):
        config = GreedyConfig(mode="minimize")
        strategy = create_strategy(config)
        assert isinstance(strategy, GreedySelection)
        assert strategy.mode == "minimize"
