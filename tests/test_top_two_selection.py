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
        assert config.adaptive_disagreement is True
        assert config.disagreement_high_threshold == 0.8
        assert config.disagreement_low_threshold == 0.3
        assert config.disagreement_decay_rate == 0.95
        assert config.ema_alpha == 0.02
        assert config.heated_scale_min == 1.0
        assert config.disagreement_window == 200

    def test_custom_disagreement_fields(self):
        config = TopTwoConfig(
            adaptive_disagreement=True,
            disagreement_high_threshold=0.9,
            disagreement_low_threshold=0.2,
            disagreement_decay_rate=0.9,
            ema_alpha=0.05,
            heated_scale_min=0.8,
        )
        assert config.adaptive_disagreement is True
        assert config.disagreement_high_threshold == 0.9
        assert config.disagreement_low_threshold == 0.2
        assert config.disagreement_decay_rate == 0.9
        assert config.ema_alpha == 0.05
        assert config.heated_scale_min == 0.8

    def test_disagreement_window_must_be_gt_1(self):
        with pytest.raises(Exception):
            TopTwoConfig(disagreement_window=1)

    def test_threshold_ordering_validation(self):
        """low_threshold must be < high_threshold."""
        with pytest.raises(Exception):
            TopTwoConfig(disagreement_high_threshold=0.3, disagreement_low_threshold=0.3)
        with pytest.raises(Exception):
            TopTwoConfig(disagreement_high_threshold=0.2, disagreement_low_threshold=0.5)
        # Valid ordering should not raise
        TopTwoConfig(disagreement_high_threshold=0.8, disagreement_low_threshold=0.2)

    def test_decay_rate_bounds(self):
        with pytest.raises(Exception):
            TopTwoConfig(disagreement_decay_rate=0)
        with pytest.raises(Exception):
            TopTwoConfig(disagreement_decay_rate=1.0)

    def test_ema_alpha_bounds(self):
        with pytest.raises(Exception):
            TopTwoConfig(ema_alpha=0)
        with pytest.raises(Exception):
            TopTwoConfig(ema_alpha=1.0)


class TestTopTwoSelection:
    """Core TT-TS selection behavior."""

    def test_wide_posteriors_cause_disagreement(self):
        """With wide posteriors (high std), two samples often disagree."""
        strategy = TopTwoSelection(
            mode="maximize", heated_scale=1.0, cooled_scale=1.0,
            adaptive_disagreement=False,
        )
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
        strategy = TopTwoSelection(
            mode="maximize", heated_scale=1.0, cooled_scale=1.0,
            adaptive_disagreement=False,
        )
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
        strategy = TopTwoSelection(
            mode="minimize", heated_scale=1.0, cooled_scale=1.0,
            adaptive_disagreement=False,
        )
        reagents = _make_reagents(
            means=[5.0, 1.0, 3.0],
            stds=[0.001, 0.001, 0.001],
            n_samples_list=[10, 10, 10],
        )
        rng = np.random.default_rng(42)
        selections = [strategy.select_reagent(reagents, rng=rng) for _ in range(50)]
        assert all(s == 1 for s in selections)

    def test_disallow_mask(self):
        strategy = TopTwoSelection(
            mode="maximize", heated_scale=1.0, cooled_scale=1.0,
            adaptive_disagreement=False,
        )
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

    def test_per_component_heated_scale_used(self):
        """select_reagent uses per-component heated_scale when set."""
        strategy = TopTwoSelection(
            mode="maximize", heated_scale=1.5, adaptive_disagreement=False,
        )
        # Set per-component scale for component 0 to 1.0
        strategy._heated_scale_per_component[0] = 1.0
        strategy.current_component_idx = 0  # component 0 is heated

        reagents = _make_reagents(
            means=[5.0, 5.0, 5.0],
            stds=[1.0, 1.0, 1.0],
            n_samples_list=[10, 10, 10],
        )
        rng = np.random.default_rng(42)
        # Just verify it doesn't crash and uses the per-component scale
        strategy.select_reagent(reagents, rng=rng, component_idx=0)
        assert strategy.effective_heated_scale == 1.0


class TestDisagreementAdaptation:
    """Tests for bidirectional disagreement-rate adaptive thermal cycling."""

    def test_disabled_no_scale_change(self):
        """adaptive_disagreement=False: no per-component adaptation occurs."""
        strategy = TopTwoSelection(
            mode="maximize",
            heated_scale=1.5,
            adaptive_disagreement=False,
        )
        reagents = _make_reagents(
            means=[1.0, 100.0, 1.0],
            stds=[0.001, 0.001, 0.001],
            n_samples_list=[100, 100, 100],
        )
        rng = np.random.default_rng(42)
        for _ in range(200):
            strategy.select_reagent(reagents, rng=rng)
            strategy.adapt_heated_scale()

        assert strategy.heated_scale == 1.5  # global unchanged
        assert len(strategy._heated_scale_per_component) == 0  # no per-comp changes

    def test_low_disagreement_increases_scale(self):
        """When disagreement is low, per-component heated_scale increases."""
        strategy = TopTwoSelection(
            mode="maximize",
            heated_scale=1.5,
            heated_scale_max=5.0,
            adaptive_disagreement=True,
            disagreement_low_threshold=0.3,
            disagreement_high_threshold=0.8,
            ema_alpha=0.1,  # fast response for testing
        )
        # Clear winner → zero disagreement
        reagents = _make_reagents(
            means=[1.0, 100.0, 1.0],
            stds=[0.001, 0.001, 0.001],
            n_samples_list=[100, 100, 100],
        )
        rng = np.random.default_rng(42)
        strategy.current_component_idx = 0
        for _ in range(100):
            strategy.select_reagent(reagents, rng=rng, component_idx=0)
        strategy.adapt_heated_scale()

        # Per-component scale for comp 0 should have increased
        comp_scale = strategy._heated_scale_per_component.get(0, strategy.heated_scale)
        assert comp_scale > 1.5

    def test_high_disagreement_reduces_scale(self):
        """When disagreement is saturated, per-component heated_scale reduces toward 1.0."""
        strategy = TopTwoSelection(
            mode="maximize",
            heated_scale=1.5,
            adaptive_disagreement=True,
            disagreement_low_threshold=0.3,
            disagreement_high_threshold=0.8,
            ema_alpha=0.1,  # fast response for testing
            heated_scale_min=1.0,
        )
        # Many reagents with equal means + wide stds → near-100% disagreement
        # (3 reagents only gives ~67% disagreement; need many to exceed 0.8)
        n = 50
        reagents = _make_reagents(
            means=[5.0] * n,
            stds=[10.0] * n,
            n_samples_list=[10] * n,
        )
        rng = np.random.default_rng(42)
        strategy.current_component_idx = 0
        for _ in range(100):
            strategy.select_reagent(reagents, rng=rng, component_idx=0)
        strategy.adapt_heated_scale()

        # Per-component scale for comp 0 should have decreased
        comp_scale = strategy._heated_scale_per_component.get(0, strategy.heated_scale)
        assert comp_scale < 1.5

    def test_healthy_disagreement_no_change(self):
        """Disagreement in the dead zone (0.3-0.8) produces no adaptation."""
        strategy = TopTwoSelection(
            mode="maximize",
            heated_scale=1.5,
            adaptive_disagreement=True,
            disagreement_low_threshold=0.3,
            disagreement_high_threshold=0.8,
            ema_alpha=0.05,
        )
        # Manually set EMA to middle of dead zone
        strategy._component_disagreement_ema[0] = 0.5
        strategy._component_disagreement_counts[0] = 200
        strategy.current_component_idx = 0

        result = strategy.adapt_heated_scale()
        assert result is False
        assert 0 not in strategy._heated_scale_per_component

    def test_scale_respects_floor(self):
        """heated_scale never drops below heated_scale_min."""
        strategy = TopTwoSelection(
            mode="maximize",
            heated_scale=1.05,
            adaptive_disagreement=True,
            disagreement_high_threshold=0.8,
            heated_scale_min=1.0,
            disagreement_decay_rate=0.5,  # aggressive decay for testing
            ema_alpha=0.1,
        )
        # Force high disagreement
        strategy._component_disagreement_ema[0] = 0.95
        strategy._component_disagreement_counts[0] = 200
        strategy.current_component_idx = 0

        # Decay multiple times
        for _ in range(20):
            strategy.adapt_heated_scale()

        comp_scale = strategy._heated_scale_per_component.get(0, strategy.heated_scale)
        assert comp_scale >= 1.0

    def test_scale_respects_ceiling(self):
        """heated_scale never exceeds heated_scale_max."""
        strategy = TopTwoSelection(
            mode="maximize",
            heated_scale=4.9,
            heated_scale_max=5.0,
            adaptive_disagreement=True,
            disagreement_low_threshold=0.3,
            disagreement_decay_rate=0.5,  # aggressive growth for testing
            ema_alpha=0.1,
        )
        # Force low disagreement
        strategy._component_disagreement_ema[0] = 0.1
        strategy._component_disagreement_counts[0] = 200
        strategy.current_component_idx = 0

        for _ in range(20):
            strategy.adapt_heated_scale()

        comp_scale = strategy._heated_scale_per_component.get(0, strategy.heated_scale)
        assert comp_scale <= 5.0

    def test_ema_warmup_no_adaptation(self):
        """No adaptation before sufficient EMA observations."""
        strategy = TopTwoSelection(
            mode="maximize",
            heated_scale=1.5,
            adaptive_disagreement=True,
            ema_alpha=0.02,  # requires ~50 observations minimum
        )
        # Manually set high disagreement but low count
        strategy._component_disagreement_ema[0] = 0.95
        strategy._component_disagreement_counts[0] = 10  # below min_obs
        strategy.current_component_idx = 0

        result = strategy.adapt_heated_scale()
        assert result is False

    def test_per_component_ema_independence(self):
        """Disagreement for component 0 does not affect component 1's EMA."""
        strategy = TopTwoSelection(
            mode="maximize",
            heated_scale=1.5,
            adaptive_disagreement=True,
            ema_alpha=0.1,
        )

        # Record high disagreement for comp 0
        for _ in range(100):
            strategy._record_disagreement(True, component_idx=0)
        # Record low disagreement for comp 1
        for _ in range(100):
            strategy._record_disagreement(False, component_idx=1)

        assert strategy._component_disagreement_ema[0] > 0.8
        assert strategy._component_disagreement_ema[1] < 0.2

    def test_reset_clears_all_state(self):
        """reset_temperature clears per-component state."""
        strategy = TopTwoSelection(
            mode="maximize",
            heated_scale=1.5,
            adaptive_disagreement=True,
            ema_alpha=0.1,
        )
        # Accumulate some state
        for _ in range(100):
            strategy._record_disagreement(True, component_idx=0)
        strategy._heated_scale_per_component[0] = 1.2

        assert strategy._total_selections > 0
        assert len(strategy._component_disagreement_ema) > 0
        assert len(strategy._heated_scale_per_component) > 0

        strategy.reset_temperature()

        assert strategy.heated_scale == 1.5
        assert strategy._heated_scale_per_component == {}
        assert strategy._component_disagreement_ema == {}
        assert strategy._component_disagreement_counts == {}
        assert strategy._disagreement_buffer == []
        assert strategy._disagreement_rate == 1.0
        assert strategy._total_selections == 0

    def test_convergence_on_saturated_component(self):
        """With sustained high disagreement, heated_scale converges near 1.0."""
        strategy = TopTwoSelection(
            mode="maximize",
            heated_scale=1.5,
            adaptive_disagreement=True,
            disagreement_high_threshold=0.8,
            disagreement_decay_rate=0.95,
            heated_scale_min=1.0,
            ema_alpha=0.1,
        )
        # Simulate thrombin-like scenario: always high disagreement
        strategy._component_disagreement_ema[0] = 0.95
        strategy._component_disagreement_counts[0] = 200
        strategy.current_component_idx = 0

        for _ in range(100):
            strategy.adapt_heated_scale()

        comp_scale = strategy._heated_scale_per_component[0]
        # Should have converged close to 1.0
        assert comp_scale < 1.05


class TestComponentState:
    """Tests for get_component_state diagnostics."""

    def test_returns_dict(self):
        strategy = TopTwoSelection(adaptive_disagreement=True)
        reagents = _make_reagents([5.0, 5.0], [1.0, 1.0], [10, 10])
        state = strategy.get_component_state(reagents, 0, 5, 100)
        assert isinstance(state, dict)
        assert "gmic" in state
        assert "heated_scale" in state
        assert "disagreement_ema" in state
        assert "is_heated" in state

    def test_heated_component_flag(self):
        strategy = TopTwoSelection(adaptive_disagreement=True)
        strategy.current_component_idx = 1
        reagents = _make_reagents([5.0, 5.0], [1.0, 1.0], [10, 10])

        state0 = strategy.get_component_state(reagents, 0, 0, 10)
        state1 = strategy.get_component_state(reagents, 1, 0, 10)
        assert state0["is_heated"] is False
        assert state1["is_heated"] is True


class TestFactoryRoundtrip:
    """Config fields wire through factory correctly."""

    def test_default_factory(self):
        config = TopTwoConfig()
        strategy = create_strategy(config)
        assert isinstance(strategy, TopTwoSelection)
        assert strategy.mode == "maximize"
        assert strategy.beta == 0.5
        assert strategy.adaptive_disagreement is True

    def test_disagreement_fields_wired(self):
        config = TopTwoConfig(
            adaptive_disagreement=True,
            disagreement_high_threshold=0.9,
            disagreement_low_threshold=0.2,
            disagreement_decay_rate=0.9,
            ema_alpha=0.05,
            heated_scale_min=0.8,
        )
        strategy = create_strategy(config)
        assert isinstance(strategy, TopTwoSelection)
        assert strategy.adaptive_disagreement is True
        assert strategy.disagreement_high_threshold == 0.9
        assert strategy.disagreement_low_threshold == 0.2
        assert strategy.disagreement_decay_rate == 0.9
        assert strategy.ema_alpha == 0.05
        assert strategy.heated_scale_min == 0.8

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
            disagreement_high_threshold=0.9,
            disagreement_low_threshold=0.2,
            disagreement_decay_rate=0.9,
            ema_alpha=0.05,
            heated_scale_min=0.8,
            disagreement_window=50,
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
        assert strategy.disagreement_high_threshold == 0.9
        assert strategy.disagreement_low_threshold == 0.2
        assert strategy.disagreement_decay_rate == 0.9
        assert strategy.ema_alpha == 0.05
        assert strategy.heated_scale_min == 0.8
        assert strategy.disagreement_window == 50

    def test_disabled_backward_compat(self):
        """Explicit adaptive_disagreement=False disables adaptation."""
        config = TopTwoConfig(adaptive_disagreement=False)
        strategy = create_strategy(config)
        assert strategy.adaptive_disagreement is False


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
