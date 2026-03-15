"""Tests for Top-Two Thompson Sampling (TT-TS) selection strategy."""

import pytest
import numpy as np

from TACTICS.thompson_sampling.strategies.top_two_selection import TopTwoSelection
from TACTICS.thompson_sampling.strategies.config import TopTwoConfig
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


class TestTopTwoConfig:
    """TopTwoConfig Pydantic model validation."""

    def test_default_values(self):
        config = TopTwoConfig()
        assert config.strategy_type == "top_two"
        assert config.mode == "maximize"
        assert config.beta == 0.5
        assert config.heated_scale == 1.5
        assert config.cooled_scale == 0.75
        assert config.adaptive_disagreement is False
        assert config.disagreement_window == 200
        assert config.disagreement_threshold == 0.3
        assert config.disagreement_scale_increment == 0.05

    def test_custom_disagreement_fields(self):
        config = TopTwoConfig(
            adaptive_disagreement=True,
            disagreement_window=100,
            disagreement_threshold=0.2,
            disagreement_scale_increment=0.1,
        )
        assert config.adaptive_disagreement is True
        assert config.disagreement_window == 100
        assert config.disagreement_threshold == 0.2
        assert config.disagreement_scale_increment == 0.1

    def test_disagreement_window_must_be_gt_1(self):
        with pytest.raises(Exception):
            TopTwoConfig(disagreement_window=1)

    def test_disagreement_threshold_bounds(self):
        with pytest.raises(Exception):
            TopTwoConfig(disagreement_threshold=0)
        with pytest.raises(Exception):
            TopTwoConfig(disagreement_threshold=1.5)

    def test_disagreement_scale_increment_non_negative(self):
        # 0 is fine (effectively disables adaptation)
        config = TopTwoConfig(disagreement_scale_increment=0)
        assert config.disagreement_scale_increment == 0
        with pytest.raises(Exception):
            TopTwoConfig(disagreement_scale_increment=-0.1)


class TestTopTwoSelection:
    """Core TT-TS selection behavior."""

    def test_wide_posteriors_cause_disagreement(self):
        """With wide posteriors (high std), two samples often disagree."""
        strategy = TopTwoSelection(mode="maximize", heated_scale=1.0, cooled_scale=1.0)
        reagents = _make_reagents(
            means=[5.0, 5.0, 5.0],
            stds=[10.0, 10.0, 10.0],
            n_samples_list=[10, 10, 10],
        )
        rng = np.random.default_rng(42)
        selections = [strategy.select_reagent(reagents, rng=rng) for _ in range(100)]
        # With equal means and wide stds, we should see multiple reagents selected
        unique = set(selections)
        assert len(unique) > 1

    def test_tight_posteriors_reduce_disagreement(self):
        """With one dominant mean and tight posteriors, selection converges."""
        strategy = TopTwoSelection(mode="maximize", heated_scale=1.0, cooled_scale=1.0)
        reagents = _make_reagents(
            means=[1.0, 10.0, 1.0],
            stds=[0.001, 0.001, 0.001],
            n_samples_list=[100, 100, 100],
        )
        rng = np.random.default_rng(42)
        selections = [strategy.select_reagent(reagents, rng=rng) for _ in range(50)]
        # With dominant mean at index 1, nearly all selections should be 1
        assert all(s == 1 for s in selections)

    def test_minimize_mode(self):
        strategy = TopTwoSelection(mode="minimize", heated_scale=1.0, cooled_scale=1.0)
        reagents = _make_reagents(
            means=[5.0, 1.0, 3.0],
            stds=[0.001, 0.001, 0.001],
            n_samples_list=[10, 10, 10],
        )
        rng = np.random.default_rng(42)
        selections = [strategy.select_reagent(reagents, rng=rng) for _ in range(50)]
        assert all(s == 1 for s in selections)

    def test_disallow_mask(self):
        strategy = TopTwoSelection(mode="maximize", heated_scale=1.0, cooled_scale=1.0)
        reagents = _make_reagents(
            means=[1.0, 100.0, 3.0],
            stds=[0.001, 0.001, 0.001],
            n_samples_list=[10, 10, 10],
        )
        rng = np.random.default_rng(42)
        selections = [
            strategy.select_reagent(reagents, disallow_mask={1}, rng=rng)
            for _ in range(20)
        ]
        assert all(s == 2 for s in selections)


class TestDisagreementAdaptation:
    """Tests for the disagreement-rate adaptive thermal cycling mechanism."""

    def test_disabled_no_scale_change(self):
        """adaptive_disagreement=False: heated_scale stays constant even with zero disagreement."""
        strategy = TopTwoSelection(
            mode="maximize",
            heated_scale=1.5,
            adaptive_disagreement=False,
            disagreement_window=10,
            disagreement_threshold=0.3,
            disagreement_scale_increment=0.1,
        )
        # Force zero disagreement by using tight posteriors with a clear winner
        reagents = _make_reagents(
            means=[1.0, 100.0, 1.0],
            stds=[0.001, 0.001, 0.001],
            n_samples_list=[100, 100, 100],
        )
        rng = np.random.default_rng(42)
        for _ in range(50):
            strategy.select_reagent(reagents, rng=rng)

        assert strategy.heated_scale == 1.5  # unchanged

    def test_enabled_increases_heated_scale(self):
        """adaptive_disagreement=True: heated_scale increases when disagreement is low."""
        strategy = TopTwoSelection(
            mode="maximize",
            heated_scale=1.5,
            heated_scale_max=5.0,
            adaptive_disagreement=True,
            disagreement_window=10,
            disagreement_threshold=0.3,
            disagreement_scale_increment=0.1,
        )
        # Force zero disagreement with a clear winner
        reagents = _make_reagents(
            means=[1.0, 100.0, 1.0],
            stds=[0.001, 0.001, 0.001],
            n_samples_list=[100, 100, 100],
        )
        rng = np.random.default_rng(42)
        for _ in range(50):
            strategy.select_reagent(reagents, rng=rng)

        # heated_scale should have increased
        assert strategy.heated_scale > 1.5
        # disagreement_rate should be near 0
        assert strategy.disagreement_rate < 0.3

    def test_heated_scale_capped_at_max(self):
        """heated_scale never exceeds heated_scale_max."""
        strategy = TopTwoSelection(
            mode="maximize",
            heated_scale=4.9,
            heated_scale_max=5.0,
            adaptive_disagreement=True,
            disagreement_window=5,
            disagreement_threshold=0.5,
            disagreement_scale_increment=0.5,
        )
        reagents = _make_reagents(
            means=[1.0, 100.0, 1.0],
            stds=[0.001, 0.001, 0.001],
            n_samples_list=[100, 100, 100],
        )
        rng = np.random.default_rng(42)
        for _ in range(100):
            strategy.select_reagent(reagents, rng=rng)

        assert strategy.heated_scale <= 5.0

    def test_buffer_not_active_before_window_fills(self):
        """Rate stays at 1.0 until the buffer reaches disagreement_window size."""
        strategy = TopTwoSelection(
            mode="maximize",
            heated_scale=1.5,
            adaptive_disagreement=True,
            disagreement_window=100,
            disagreement_threshold=0.3,
            disagreement_scale_increment=0.1,
        )
        reagents = _make_reagents(
            means=[1.0, 100.0, 1.0],
            stds=[0.001, 0.001, 0.001],
            n_samples_list=[100, 100, 100],
        )
        rng = np.random.default_rng(42)
        # Run fewer than window size
        for _ in range(50):
            strategy.select_reagent(reagents, rng=rng)

        # Rate should still be 1.0 (initial) and heated_scale unchanged
        assert strategy.disagreement_rate == 1.0
        assert strategy.heated_scale == 1.5

    def test_reset_clears_disagreement_state(self):
        """reset_temperature clears all disagreement tracking state."""
        strategy = TopTwoSelection(
            mode="maximize",
            heated_scale=1.5,
            adaptive_disagreement=True,
            disagreement_window=10,
            disagreement_threshold=0.3,
            disagreement_scale_increment=0.1,
        )
        reagents = _make_reagents(
            means=[1.0, 100.0, 1.0],
            stds=[0.001, 0.001, 0.001],
            n_samples_list=[100, 100, 100],
        )
        rng = np.random.default_rng(42)
        for _ in range(50):
            strategy.select_reagent(reagents, rng=rng)

        # State should be modified
        assert strategy._total_selections > 0
        assert strategy.heated_scale > 1.5

        # Reset
        strategy.reset_temperature()

        assert strategy.heated_scale == 1.5
        assert strategy._disagreement_buffer == []
        assert strategy._disagreement_rate == 1.0
        assert strategy._total_selections == 0

    def test_high_disagreement_no_adaptation(self):
        """When disagreement is above threshold, heated_scale stays unchanged."""
        strategy = TopTwoSelection(
            mode="maximize",
            heated_scale=1.5,
            adaptive_disagreement=True,
            disagreement_window=10,
            disagreement_threshold=0.3,
            disagreement_scale_increment=0.1,
        )
        # Equal means + wide stds → high disagreement
        reagents = _make_reagents(
            means=[5.0, 5.0, 5.0],
            stds=[10.0, 10.0, 10.0],
            n_samples_list=[10, 10, 10],
        )
        rng = np.random.default_rng(42)
        for _ in range(50):
            strategy.select_reagent(reagents, rng=rng)

        # With equal means and wide stds, disagreement should be high
        assert strategy.disagreement_rate > 0.3
        assert strategy.heated_scale == 1.5


class TestFactoryRoundtrip:
    """Config fields wire through factory correctly."""

    def test_default_factory(self):
        config = TopTwoConfig()
        strategy = create_strategy(config)
        assert isinstance(strategy, TopTwoSelection)
        assert strategy.mode == "maximize"
        assert strategy.beta == 0.5
        assert strategy.adaptive_disagreement is False

    def test_disagreement_fields_wired(self):
        config = TopTwoConfig(
            adaptive_disagreement=True,
            disagreement_window=100,
            disagreement_threshold=0.2,
            disagreement_scale_increment=0.1,
        )
        strategy = create_strategy(config)
        assert isinstance(strategy, TopTwoSelection)
        assert strategy.adaptive_disagreement is True
        assert strategy.disagreement_window == 100
        assert strategy.disagreement_threshold == 0.2
        assert strategy.disagreement_scale_increment == 0.1

    def test_minimize_factory(self):
        config = TopTwoConfig(mode="minimize")
        strategy = create_strategy(config)
        assert strategy.mode == "minimize"

    def test_all_fields_roundtrip(self):
        config = TopTwoConfig(
            mode="minimize",
            beta=0.7,
            heated_scale=2.0,
            cooled_scale=0.5,
            min_observations=10,
            adaptive_temperature=True,
            scale_increment=0.02,
            cooled_scale_increment=0.002,
            efficiency_threshold=0.15,
            heated_scale_max=8.0,
            adaptive_disagreement=True,
            disagreement_window=50,
            disagreement_threshold=0.4,
            disagreement_scale_increment=0.08,
        )
        strategy = create_strategy(config)
        assert strategy.mode == "minimize"
        assert strategy.beta == 0.7
        assert strategy.heated_scale == 2.0
        assert strategy.cooled_scale == 0.5
        assert strategy.min_observations == 10
        assert strategy.adaptive_temperature is True
        assert strategy.scale_increment == 0.02
        assert strategy.cooled_scale_increment == 0.002
        assert strategy.efficiency_threshold == 0.15
        assert strategy.heated_scale_max == 8.0
        assert strategy.adaptive_disagreement is True
        assert strategy.disagreement_window == 50
        assert strategy.disagreement_threshold == 0.4
        assert strategy.disagreement_scale_increment == 0.08


class TestOrthogonality:
    """GMIC rotation and weighted rotation work independently of disagreement adaptation."""

    def test_rotate_component_roundrobin(self):
        strategy = TopTwoSelection(adaptive_disagreement=True)
        assert strategy.current_component_idx == 0
        strategy.rotate_component(3)
        assert strategy.current_component_idx == 1
        strategy.rotate_component(3)
        assert strategy.current_component_idx == 2
        strategy.rotate_component(3)
        assert strategy.current_component_idx == 0

    def test_rotate_component_weighted(self):
        strategy = TopTwoSelection(adaptive_disagreement=True)
        rng = np.random.default_rng(42)
        # Component 0: critical (one dominant reagent) → high GMIC → less likely heated
        # Component 1: flexible (all similar) → low GMIC → more likely heated
        critical = _make_reagents([1.0, 1.0, 10.0], [0.5, 0.5, 0.5], [10, 10, 10])
        flexible = _make_reagents([5.0, 5.0, 5.0], [1.0, 1.0, 1.0], [10, 10, 10])
        heated_indices = []
        for _ in range(200):
            strategy.rotate_component_weighted(2, [critical, flexible], rng=rng)
            heated_indices.append(strategy.current_component_idx)
        # Flexible component (index 1) should be heated more often
        frac_1 = sum(1 for x in heated_indices if x == 1) / len(heated_indices)
        assert frac_1 > 0.7

    def test_gmic_criticality(self):
        strategy = TopTwoSelection(adaptive_disagreement=True, min_observations=3)
        # Reagents with diverse means → high GMIC
        reagents = _make_reagents(
            means=[1.0, 10.0, 5.0],
            stds=[0.5, 0.5, 0.5],
            n_samples_list=[10, 10, 10],
        )
        gmic = strategy.get_component_criticality(reagents)
        assert gmic > 0

        # Reagents with similar means → low GMIC
        reagents_similar = _make_reagents(
            means=[5.0, 5.01, 4.99],
            stds=[0.5, 0.5, 0.5],
            n_samples_list=[10, 10, 10],
        )
        gmic_similar = strategy.get_component_criticality(reagents_similar)
        assert gmic_similar < gmic
