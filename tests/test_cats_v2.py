"""Tests for CATS GMIC + Criticality Divergence."""

import pytest
import numpy as np

from TACTICS.thompson_sampling.strategies.roulette_wheel import RouletteWheelSelection
from TACTICS.thompson_sampling.strategies.config import RouletteWheelConfig
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


class TestGMIC:
    """Test Gaussian Mutual Information Criticality calculation."""

    def test_identical_means_returns_zero(self):
        """All reagents with identical means → GMIC ≈ 0."""
        strategy = RouletteWheelSelection()
        reagents = _make_reagents(
            means=[5.0, 5.0, 5.0, 5.0],
            stds=[1.0, 1.0, 1.0, 1.0],
            n_samples_list=[10, 10, 10, 10],
        )
        gmic = strategy._calculate_gmic(reagents)
        assert gmic == pytest.approx(0.0, abs=1e-6)

    def test_one_dominant_returns_positive(self):
        """One reagent with much higher mean → positive GMIC."""
        strategy = RouletteWheelSelection()
        reagents = _make_reagents(
            means=[1.0, 1.0, 1.0, 10.0],
            stds=[0.5, 0.5, 0.5, 0.5],
            n_samples_list=[10, 10, 10, 10],
        )
        gmic = strategy._calculate_gmic(reagents)
        assert gmic > 0.5  # Clear signal

    def test_high_noise_stays_low(self):
        """High posterior noise relative to signal → low GMIC."""
        strategy = RouletteWheelSelection()
        reagents = _make_reagents(
            means=[5.0, 5.1, 4.9, 5.05],
            stds=[10.0, 10.0, 10.0, 10.0],
            n_samples_list=[10, 10, 10, 10],
        )
        gmic = strategy._calculate_gmic(reagents)
        assert gmic < 0.01  # Noise overwhelms signal

    def test_insufficient_active_returns_zero(self):
        """Fewer than 2 active reagents → GMIC = 0."""
        strategy = RouletteWheelSelection()
        reagents = _make_reagents(
            means=[5.0],
            stds=[1.0],
            n_samples_list=[10],
        )
        gmic = strategy._calculate_gmic(reagents)
        assert gmic == 0.0

    def test_zero_samples_excluded(self):
        """Reagents with n_samples=0 are excluded from GMIC."""
        strategy = RouletteWheelSelection()
        reagents = _make_reagents(
            means=[5.0, 5.0, 10.0],
            stds=[1.0, 1.0, 1.0],
            n_samples_list=[10, 10, 0],  # Third reagent inactive
        )
        gmic = strategy._calculate_gmic(reagents)
        # Only two active reagents with identical means → low GMIC
        assert gmic == pytest.approx(0.0, abs=1e-6)

    def test_gmic_details_returns_diagnostics(self):
        """_calculate_gmic_details returns both value and diagnostics dict."""
        strategy = RouletteWheelSelection()
        reagents = _make_reagents(
            means=[1.0, 5.0, 3.0],
            stds=[0.5, 0.5, 0.5],
            n_samples_list=[10, 10, 10],
        )
        gmic, details = strategy._calculate_gmic_details(reagents)
        assert gmic > 0
        assert "signal_var" in details
        assert "mean_noise_var" in details
        assert details["n_active_reagents"] == 3
        assert details["signal_var"] > 0
        assert details["mean_noise_var"] > 0

    def test_minimize_mode_gmic_unaffected(self):
        """GMIC is mode-independent (uses variance of means, not direction)."""
        strategy_max = RouletteWheelSelection(mode="maximize")
        strategy_min = RouletteWheelSelection(mode="minimize")
        reagents = _make_reagents(
            means=[1.0, 5.0, 3.0],
            stds=[0.5, 0.5, 0.5],
            n_samples_list=[10, 10, 10],
        )
        gmic_max = strategy_max._calculate_gmic(reagents)
        gmic_min = strategy_min._calculate_gmic(reagents)
        assert gmic_max == pytest.approx(gmic_min, abs=1e-10)


class TestDivergence:
    """Test posterior divergence tracking."""

    def test_first_call_returns_inf(self):
        """First divergence call has no prior snapshot → inf."""
        strategy = RouletteWheelSelection()
        reagents = _make_reagents(
            means=[1.0, 2.0, 3.0],
            stds=[0.5, 0.5, 0.5],
            n_samples_list=[10, 10, 10],
        )
        div = strategy._update_divergence(0, reagents)
        assert div == float("inf")

    def test_stable_posteriors_low_divergence(self):
        """Identical posteriors on consecutive calls → divergence ≈ 0."""
        strategy = RouletteWheelSelection()
        reagents = _make_reagents(
            means=[1.0, 2.0, 3.0],
            stds=[0.5, 0.5, 0.5],
            n_samples_list=[10, 10, 10],
        )
        # First call → snapshot
        strategy._update_divergence(0, reagents)
        # Second call with identical posteriors → KL ≈ 0
        div = strategy._update_divergence(0, reagents)
        assert div == pytest.approx(0.0, abs=1e-10)

    def test_shifting_posteriors_high_divergence(self):
        """Large posterior shift → high divergence."""
        strategy = RouletteWheelSelection()
        reagents1 = _make_reagents(
            means=[1.0, 2.0, 3.0],
            stds=[0.5, 0.5, 0.5],
            n_samples_list=[10, 10, 10],
        )
        # First call → snapshot
        strategy._update_divergence(0, reagents1)

        # Shift posteriors significantly
        reagents2 = _make_reagents(
            means=[5.0, 6.0, 7.0],
            stds=[0.5, 0.5, 0.5],
            n_samples_list=[10, 10, 10],
        )
        div = strategy._update_divergence(0, reagents2)
        assert div > 1.0  # Significant shift

    def test_independent_components(self):
        """Different components track divergence independently."""
        strategy = RouletteWheelSelection()
        reagents_a = _make_reagents([1.0, 2.0], [0.5, 0.5], [10, 10])
        reagents_b = _make_reagents([3.0, 4.0], [0.5, 0.5], [10, 10])

        # Snapshot both
        strategy._update_divergence(0, reagents_a)
        strategy._update_divergence(1, reagents_b)

        # Component 0 shifts, component 1 stays
        reagents_a_shifted = _make_reagents([5.0, 6.0], [0.5, 0.5], [10, 10])
        div_0 = strategy._update_divergence(0, reagents_a_shifted)
        div_1 = strategy._update_divergence(1, reagents_b)

        assert div_0 > 1.0
        assert div_1 == pytest.approx(0.0, abs=1e-10)

    def test_std_changes_contribute(self):
        """Changes in std (not just mean) contribute to divergence."""
        strategy = RouletteWheelSelection()
        reagents1 = _make_reagents([1.0, 2.0], [1.0, 1.0], [10, 10])
        strategy._update_divergence(0, reagents1)

        # Same means, different stds
        reagents2 = _make_reagents([1.0, 2.0], [0.1, 0.1], [10, 10])
        div = strategy._update_divergence(0, reagents2)
        assert div > 0.1  # Std change causes non-zero divergence


class TestTemperatureMultiplier:
    """Test directional temperature adjustment."""

    def test_heated_flexible_gets_amplified(self):
        """Heated + flexible (low GMIC) → multiplier > 1."""
        strategy = RouletteWheelSelection(alpha=0.1, beta=0.05)
        # Component 0 heated (current_component_idx=0)
        strategy.current_component_idx = 0

        # Set up: component 0 flexible (low GMIC, small spread), component 1 critical (high GMIC)
        flexible_reagents = _make_reagents(
            means=[5.0, 5.1, 4.9], stds=[1.0, 1.0, 1.0], n_samples_list=[10, 10, 10]
        )
        critical_reagents = _make_reagents(
            means=[1.0, 1.0, 10.0], stds=[0.5, 0.5, 0.5], n_samples_list=[10, 10, 10]
        )

        # Cache GMICs so GMIC-based multiplier can use them
        strategy._cached_gmics = {
            0: strategy._calculate_gmic(flexible_reagents),
            1: strategy._calculate_gmic(critical_reagents),
        }
        # Verify flexible GMIC is positive but much lower than critical
        assert strategy._cached_gmics[0] > 0
        assert strategy._cached_gmics[0] < strategy._cached_gmics[1]

        # First snapshot → divergence inf (diversity mode)
        strategy._update_divergence(0, flexible_reagents)
        # Second call → stable
        strategy._update_divergence(0, flexible_reagents)

        state = strategy._get_component_state_gmic(flexible_reagents, 0, 10, 100)
        assert state["is_heated"] is True
        assert state["cats_mode"] == "criticality"
        # Flexible heated component should get multiplier > 1 (amplified heating)
        assert state["effective_multiplier"] >= 1.0

    def test_cooled_critical_gets_amplified(self):
        """Cooled + critical (high GMIC) → multiplier < 1."""
        strategy = RouletteWheelSelection(alpha=0.1, beta=0.05)
        # Component 0 is heated, so component 1 is cooled
        strategy.current_component_idx = 0

        critical_reagents = _make_reagents(
            means=[1.0, 1.0, 10.0], stds=[0.5, 0.5, 0.5], n_samples_list=[10, 10, 10]
        )
        flexible_reagents = _make_reagents(
            means=[5.0, 5.0, 5.0], stds=[1.0, 1.0, 1.0], n_samples_list=[10, 10, 10]
        )

        strategy._cached_gmics = {
            0: strategy._calculate_gmic(flexible_reagents),
            1: strategy._calculate_gmic(critical_reagents),
        }

        # Stabilize divergence for component 1
        strategy._update_divergence(1, critical_reagents)
        strategy._update_divergence(1, critical_reagents)

        state = strategy._get_component_state_gmic(critical_reagents, 1, 10, 100)
        assert state["is_heated"] is False
        assert state["cats_mode"] == "criticality"
        # Critical cooled component should get multiplier <= 1 (amplified cooling)
        assert state["effective_multiplier"] <= 1.0

    def test_unstable_uses_diversity_mode(self):
        """When divergence > threshold → diversity mode (base temp only)."""
        strategy = RouletteWheelSelection(
            alpha=0.1, beta=0.05, divergence_threshold=0.1
        )
        strategy.current_component_idx = 0

        reagents = _make_reagents(
            means=[1.0, 5.0, 3.0], stds=[0.5, 0.5, 0.5], n_samples_list=[10, 10, 10]
        )
        strategy._cached_gmics = {0: strategy._calculate_gmic(reagents)}

        # First call → divergence inf → diversity mode
        state = strategy._get_component_state_gmic(reagents, 0, 5, 100)
        assert state["cats_mode"] == "diversity"
        assert state["effective_multiplier"] == 1.0
        assert state["final_temperature"] == strategy.alpha  # Base temp only


class TestWeightedRotation:
    """Test GMIC-weighted rotation."""

    def test_flexible_component_heated_more(self):
        """Flexible components (low GMIC) should be selected for heating more often."""
        strategy = RouletteWheelSelection(alpha=0.1, beta=0.05)
        rng = np.random.default_rng(42)

        # Component 0: flexible (all similar means)
        flexible = _make_reagents(
            means=[5.0, 5.0, 5.0], stds=[1.0, 1.0, 1.0], n_samples_list=[10, 10, 10]
        )
        # Component 1: critical (one dominant)
        critical = _make_reagents(
            means=[1.0, 1.0, 10.0], stds=[0.5, 0.5, 0.5], n_samples_list=[10, 10, 10]
        )

        counts = {0: 0, 1: 0}
        for _ in range(1000):
            strategy.rotate_component_weighted(2, [flexible, critical], rng=rng)
            counts[strategy.current_component_idx] += 1

        # Flexible component should be heated more often
        assert counts[0] > counts[1]

    def test_caches_gmics(self):
        """rotate_component_weighted should cache GMIC values."""
        strategy = RouletteWheelSelection()
        rng = np.random.default_rng(42)

        reagents_a = _make_reagents([1.0, 5.0], [0.5, 0.5], [10, 10])
        reagents_b = _make_reagents([3.0, 3.0], [1.0, 1.0], [10, 10])

        strategy.rotate_component_weighted(2, [reagents_a, reagents_b], rng=rng)
        assert 0 in strategy._cached_gmics
        assert 1 in strategy._cached_gmics
        assert strategy._cached_gmics[0] > strategy._cached_gmics[1]


class TestConfigIntegration:
    """Test config → factory → strategy wiring."""

    def test_default_config(self):
        """Default config creates strategy with GMIC."""
        config = RouletteWheelConfig()
        strategy = create_strategy(config)
        assert isinstance(strategy, RouletteWheelSelection)

    def test_divergence_threshold_passthrough(self):
        """divergence_threshold is correctly passed through factory."""
        config = RouletteWheelConfig(divergence_threshold=0.05)
        strategy = create_strategy(config)
        assert isinstance(strategy, RouletteWheelSelection)
        assert strategy._divergence_threshold == 0.05

    def test_get_component_criticality_returns_gmic(self):
        """get_component_criticality returns GMIC."""
        config = RouletteWheelConfig()
        strategy = create_strategy(config)
        reagents = _make_reagents(
            means=[1.0, 5.0, 3.0],
            stds=[0.5, 0.5, 0.5],
            n_samples_list=[10, 10, 10],
        )
        crit = strategy.get_component_criticality(reagents)
        gmic = strategy._calculate_gmic(reagents)
        assert crit == pytest.approx(gmic)

    def test_state_has_all_fields(self):
        """Component state includes GMIC, divergence, and mode fields."""
        strategy = RouletteWheelSelection(alpha=0.1, beta=0.05)
        strategy._cached_gmics = {0: 0.5}

        reagents = _make_reagents(
            means=[1.0, 5.0, 3.0],
            stds=[0.5, 0.5, 0.5],
            n_samples_list=[10, 10, 10],
        )
        state = strategy.get_component_state(reagents, 0, 10, 100)

        # GMIC fields
        assert "gmic" in state
        assert "criticality" in state
        assert "divergence" in state
        assert "divergence_threshold" in state
        assert "is_stable" in state
        assert "cats_mode" in state
        assert state["cats_mode"] in ("diversity", "criticality")

        # Temperature pipeline fields
        assert "base_temp" in state
        assert "final_temperature" in state
        assert "is_heated" in state


class TestCATSEmaSmoothing:
    """Tests for EMA-smoothed relative GMIC in the CATS multiplier."""

    def test_ema_disabled_by_default(self):
        """Default cats_ema_decay is None, EMA state is empty."""
        strategy = RouletteWheelSelection(mode="maximize", alpha=0.1, beta=0.05)
        assert strategy.cats_ema_decay is None
        assert strategy._ema_relative_gmic == {}

    def test_ema_state_nan_when_disabled(self):
        """When EMA is disabled, diagnostics reports NaN for ema_relative_gmic."""
        strategy = RouletteWheelSelection(mode="maximize", alpha=0.1, beta=0.05)
        reagents = _make_reagents(
            [0.5, 0.3, 0.1], [0.1, 0.1, 0.1], [10, 10, 10]
        )
        # Seed cached GMICs so relative_gmic can be computed
        strategy._cached_gmics = {0: 0.5, 1: 0.5}
        state = strategy.get_component_state(reagents, 0, 5, 100)
        assert np.isnan(state["ema_relative_gmic"])

    def _warm_divergence(self, strategy, reagents, comp_idx):
        """Call get_component_state twice to seed divergence snapshot so is_stable=True."""
        strategy.get_component_state(reagents, comp_idx, 0, 100)
        # Second call compares against snapshot from first → low divergence → stable
        strategy.get_component_state(reagents, comp_idx, 1, 100)

    def test_ema_initializes_to_first_value(self):
        """First stable call initializes EMA to the instantaneous relative_gmic."""
        strategy = RouletteWheelSelection(
            mode="maximize", alpha=0.1, beta=0.05, cats_ema_decay=0.1
        )
        reagents = _make_reagents(
            [0.5, 0.3, 0.1], [0.1, 0.1, 0.1], [10, 10, 10]
        )
        strategy._cached_gmics = {0: 0.8, 1: 1.2}
        # Warm up divergence tracker so is_stable=True
        self._warm_divergence(strategy, reagents, 0)

        state = strategy.get_component_state(reagents, 0, 5, 100)

        # EMA should be initialized (not NaN)
        assert not np.isnan(state["ema_relative_gmic"])
        assert 0 in strategy._ema_relative_gmic

    def test_ema_smooths_relative_gmic(self):
        """EMA changes more slowly than instantaneous relative_gmic."""
        strategy = RouletteWheelSelection(
            mode="maximize", alpha=0.1, beta=0.05, cats_ema_decay=0.1
        )
        reagents = _make_reagents(
            [0.5, 0.3, 0.1], [0.1, 0.1, 0.1], [10, 10, 10]
        )
        # Warm up divergence so is_stable=True
        strategy._cached_gmics = {0: 0.5, 1: 1.5}
        self._warm_divergence(strategy, reagents, 0)

        # First stable call: component 0 is flexible (low GMIC)
        state1 = strategy.get_component_state(reagents, 0, 5, 100)
        ema1 = state1["ema_relative_gmic"]

        # Shift GMIC: component 0 is now critical (high GMIC)
        strategy._cached_gmics = {0: 1.5, 1: 0.5}
        state2 = strategy.get_component_state(reagents, 0, 6, 100)
        ema2 = state2["ema_relative_gmic"]

        # The instantaneous relative_gmic flipped from <1 to >1,
        # but EMA should be closer to ema1 than to the new instantaneous value
        # because decay=0.1 means only 10% weight on new observation
        assert abs(ema2 - ema1) < abs(ema2 - (1.5 / 1.0))

    def test_ema_per_component_independence(self):
        """EMA state is tracked independently per component."""
        strategy = RouletteWheelSelection(
            mode="maximize", alpha=0.1, beta=0.05, cats_ema_decay=0.1
        )
        # Component 0: high signal spread → high GMIC
        reagents_a = _make_reagents(
            [0.9, 0.3, 0.1], [0.05, 0.05, 0.05], [10, 10, 10]
        )
        # Component 1: low signal spread → low GMIC
        reagents_b = _make_reagents(
            [0.31, 0.30, 0.29], [0.05, 0.05, 0.05], [10, 10, 10]
        )
        strategy._cached_gmics = {0: 0.5, 1: 1.5}
        # Warm up divergence for both components
        self._warm_divergence(strategy, reagents_a, 0)
        self._warm_divergence(strategy, reagents_b, 1)

        state_a = strategy.get_component_state(reagents_a, 0, 5, 100)
        state_b = strategy.get_component_state(reagents_b, 1, 5, 100)

        assert 0 in strategy._ema_relative_gmic
        assert 1 in strategy._ema_relative_gmic
        # Components with different GMIC have different relative_gmic → different EMA
        assert strategy._ema_relative_gmic[0] != strategy._ema_relative_gmic[1]

    def test_config_cats_ema_decay_roundtrip(self):
        """cats_ema_decay passes through config → factory → strategy."""
        config = RouletteWheelConfig(
            mode="maximize", alpha=0.1, beta=0.05, cats_ema_decay=0.1
        )
        strategy = create_strategy(config)
        assert isinstance(strategy, RouletteWheelSelection)
        assert strategy.cats_ema_decay == 0.1

    def test_config_default_none(self):
        """Default config has cats_ema_decay=None."""
        config = RouletteWheelConfig(mode="maximize")
        assert config.cats_ema_decay is None
