"""
Comprehensive test suite for Adaptive Bayes-UCB Selection Strategy.

Tests UCB index computation, percentile adaptation, thermal cycling,
and integration with the ThompsonSampler.
"""

import pytest
import numpy as np
from scipy import stats

from TACTICS.thompson_sampling.strategies.bayes_ucb_selection import BayesUCBSelection
from TACTICS.thompson_sampling.core.reagent import Reagent


class TestBayesUCBIndexComputation:
    """Test UCB index computation with various scenarios."""

    def test_ucb_index_with_sufficient_samples(self):
        """Test UCB computation for reagents with n >= 2."""
        strategy = BayesUCBSelection(mode="maximize", initial_p_high=0.95)

        # Create mock reagent with sufficient samples
        reagent = Reagent("test_reagent", "CC")
        reagent.mean = 5.0
        reagent.std = 1.0
        reagent.n_samples = 10
        reagent.current_phase = "search"

        # Compute UCB indices
        ucb_indices = strategy._compute_ucb_indices([reagent], percentile=0.95)

        # Verify UCB is greater than mean (exploration bonus)
        assert ucb_indices[0] > reagent.mean

        # Verify calculation matches expected Student-t formula
        df = reagent.n_samples - 1
        t_quantile = stats.t.ppf(0.95, df)
        expected_ucb = reagent.mean + reagent.std * t_quantile / np.sqrt(
            reagent.n_samples
        )
        assert np.isclose(ucb_indices[0], expected_ucb)

    def test_ucb_index_with_insufficient_samples(self):
        """Test UCB computation for reagents with n < 2 (conservative bonus)."""
        strategy = BayesUCBSelection(mode="maximize")

        # Create reagent with only 1 sample
        reagent = Reagent("test_reagent", "CC")
        reagent.mean = 5.0
        reagent.std = 1.0
        reagent.n_samples = 1
        reagent.current_phase = "search"

        ucb_indices = strategy._compute_ucb_indices([reagent], percentile=0.95)

        # Should use conservative bonus: mean + 3*std
        expected_ucb = reagent.mean + 3.0 * reagent.std
        assert np.isclose(ucb_indices[0], expected_ucb)

    def test_ucb_minimize_mode(self):
        """Test UCB computation in minimize mode (lower confidence bound)."""
        strategy = BayesUCBSelection(mode="minimize", initial_p_high=0.95)

        reagent = Reagent("test_reagent", "CC")
        reagent.mean = 5.0
        reagent.std = 1.0
        reagent.n_samples = 10
        reagent.current_phase = "search"

        ucb_indices = strategy._compute_ucb_indices([reagent], percentile=0.95)

        # In minimize mode, should subtract exploration bonus
        df = reagent.n_samples - 1
        t_quantile = stats.t.ppf(0.95, df)
        expected_ucb = reagent.mean - reagent.std * t_quantile / np.sqrt(
            reagent.n_samples
        )
        assert np.isclose(ucb_indices[0], expected_ucb)

    def test_higher_percentile_increases_exploration(self):
        """Test that higher percentile leads to larger UCB values."""
        strategy = BayesUCBSelection(mode="maximize")

        reagent = Reagent("test_reagent", "CC")
        reagent.mean = 5.0
        reagent.std = 1.0
        reagent.n_samples = 10
        reagent.current_phase = "search"

        ucb_low = strategy._compute_ucb_indices([reagent], percentile=0.80)[0]
        ucb_high = strategy._compute_ucb_indices([reagent], percentile=0.95)[0]

        # Higher percentile should give higher UCB (more exploration)
        assert ucb_high > ucb_low

    def test_ucb_multiple_reagents(self):
        """Test UCB computation across multiple reagents."""
        strategy = BayesUCBSelection(mode="maximize")

        reagents = []
        for i in range(5):
            r = Reagent(f"reagent_{i}", "CC")
            r.mean = float(i)
            r.std = 1.0
            r.n_samples = 5
            r.current_phase = "search"
            reagents.append(r)

        ucb_indices = strategy._compute_ucb_indices(reagents, percentile=0.95)

        # Should have one index per reagent
        assert len(ucb_indices) == 5

        # Reagent with higher mean should generally have higher UCB
        # (though uncertainty can affect this)
        assert ucb_indices[-1] > ucb_indices[0]


class TestCATSCriticality:
    """Test CATS criticality-based adjustment (replaces old percentile adaptation)."""

    def test_criticality_with_uniform_distribution(self):
        """Test criticality is low (flexible) when all reagents have similar means."""
        strategy = BayesUCBSelection(initial_p_high=0.90, initial_p_low=0.60)

        # Create reagents with very similar means
        reagents = []
        for i in range(5):
            r = Reagent(f"reagent_{i}", "CC")
            r.mean = 5.0 + i * 0.01  # Nearly identical means
            r.std = 1.0
            r.n_samples = 10
            reagents.append(r)

        criticality = strategy._calculate_criticality(reagents)
        # Low criticality = flexible component
        assert criticality < 0.3

    def test_criticality_with_peaked_distribution(self):
        """Test criticality is high (critical) when one reagent dominates."""
        strategy = BayesUCBSelection(initial_p_high=0.90, initial_p_low=0.60)

        reagents = []
        for i in range(5):
            r = Reagent(f"reagent_{i}", "CC")
            r.mean = 0.0 if i < 4 else 10.0  # One much better than others
            r.std = 1.0
            r.n_samples = 10
            reagents.append(r)

        criticality = strategy._calculate_criticality(reagents)
        # High criticality = critical component
        assert criticality > 0.5

    def test_criticality_returns_neutral_with_insufficient_data(self):
        """Test criticality is neutral when data is insufficient."""
        strategy = BayesUCBSelection(
            initial_p_high=0.90, initial_p_low=0.60, min_observations=5
        )

        reagents = []
        for i in range(5):
            r = Reagent(f"reagent_{i}", "CC")
            r.mean = float(i)
            r.std = 1.0
            r.n_samples = 2  # Below min_observations
            reagents.append(r)

        criticality = strategy._calculate_criticality(reagents)
        assert criticality == 0.5  # Neutral

    def test_reset_percentiles(self):
        """Test resetting percentiles to initial values."""
        strategy = BayesUCBSelection(initial_p_high=0.90, initial_p_low=0.60)

        # Change percentiles
        strategy.p_high = 0.95
        strategy.p_low = 0.50

        # Reset
        strategy.reset_percentiles()

        assert strategy.p_high == 0.90
        assert strategy.p_low == 0.60


class TestThermalCycling:
    """Test component rotation for thermal cycling."""

    def test_component_rotation(self):
        """Test cycling through components."""
        strategy = BayesUCBSelection()

        assert strategy.current_component_idx == 0

        strategy.rotate_component(n_components=3)
        assert strategy.current_component_idx == 1

        strategy.rotate_component(n_components=3)
        assert strategy.current_component_idx == 2

        strategy.rotate_component(n_components=3)
        assert strategy.current_component_idx == 0  # Should wrap around

    def test_heated_component_uses_high_percentile(self):
        """Test heated component gets higher percentile."""
        strategy = BayesUCBSelection(initial_p_high=0.95, initial_p_low=0.70)

        # Create mock reagents
        reagents = []
        for i in range(3):
            r = Reagent(f"reagent_{i}", "CC")
            r.mean = 5.0
            r.std = 1.0
            r.n_samples = 10
            r.current_phase = "search"
            reagents.append(r)

        # Set heated component to index 1
        strategy.current_component_idx = 1

        # Select from component 1 (heated)
        ucb_heated = strategy._compute_ucb_indices(reagents, strategy.p_high)

        # Select from component 0 (cooled)
        ucb_cooled = strategy._compute_ucb_indices(reagents, strategy.p_low)

        # Heated should have higher UCB values (more exploration)
        assert np.mean(ucb_heated) > np.mean(ucb_cooled)


class TestReagentSelection:
    """Test reagent selection logic."""

    def test_select_reagent_maximize(self):
        """Test selecting reagent with highest UCB in maximize mode."""
        strategy = BayesUCBSelection(mode="maximize")

        # Create reagents with different means
        reagents = []
        for i in range(5):
            r = Reagent(f"reagent_{i}", "CC")
            r.mean = float(i)  # Increasing means
            r.std = 0.1  # Small std to minimize uncertainty
            r.n_samples = 10
            r.current_phase = "search"
            reagents.append(r)

        # Select reagent
        selected_idx = strategy.select_reagent(reagents, component_idx=0)

        # Should select reagent with highest mean (reagent_4)
        assert selected_idx == 4

    def test_select_reagent_minimize(self):
        """Test selecting reagent with lowest UCB in minimize mode."""
        strategy = BayesUCBSelection(mode="minimize")

        # Create reagents with different means
        reagents = []
        for i in range(5):
            r = Reagent(f"reagent_{i}", "CC")
            r.mean = float(i)
            r.std = 0.1
            r.n_samples = 10
            r.current_phase = "search"
            reagents.append(r)

        # Select reagent
        selected_idx = strategy.select_reagent(reagents, component_idx=0)

        # Should select reagent with lowest mean (reagent_0)
        assert selected_idx == 0

    def test_select_reagent_with_disallow_mask(self):
        """Test that disallow mask is respected."""
        strategy = BayesUCBSelection(mode="maximize")

        reagents = []
        for i in range(5):
            r = Reagent(f"reagent_{i}", "CC")
            r.mean = float(i)
            r.std = 0.1
            r.n_samples = 10
            r.current_phase = "search"
            reagents.append(r)

        # Disallow the best reagent (index 4)
        disallow_mask = {4}
        selected_idx = strategy.select_reagent(reagents, disallow_mask=disallow_mask)

        # Should select next best (reagent_3)
        assert selected_idx == 3
        assert selected_idx not in disallow_mask

    def test_select_batch(self):
        """Test batch selection."""
        strategy = BayesUCBSelection(mode="maximize")

        reagents = []
        for i in range(5):
            r = Reagent(f"reagent_{i}", "CC")
            r.mean = float(i)
            r.std = 0.1
            r.n_samples = 10
            r.current_phase = "search"
            reagents.append(r)

        # Select batch of 3
        selected_indices = strategy.select_batch(reagents, batch_size=3)

        assert len(selected_indices) == 3
        # All should be reagent_4 (highest mean)
        assert all(idx == 4 for idx in selected_indices)


class TestBayesUCBConfig:
    """Test Pydantic configuration."""

    def test_config_defaults(self):
        """Test default configuration values."""
        from TACTICS.thompson_sampling.strategies.config import BayesUCBConfig

        config = BayesUCBConfig()

        assert config.strategy_type == "bayes_ucb"
        assert config.mode == "maximize"
        assert config.initial_p_high == 0.90
        assert config.initial_p_low == 0.60
        # CATS parameters (replaced old percentile adaptation parameters)
        assert config.exploration_phase_end == 0.20
        assert config.transition_phase_end == 0.60
        assert config.min_observations == 5

    def test_config_validation(self):
        """Test configuration validation."""
        from TACTICS.thompson_sampling.strategies.config import BayesUCBConfig

        # Valid config
        config = BayesUCBConfig(
            mode="minimize", initial_p_high=0.95, initial_p_low=0.80
        )
        assert config.mode == "minimize"

        # Invalid percentile (too high)
        with pytest.raises(Exception):  # Pydantic ValidationError
            BayesUCBConfig(initial_p_high=1.1)

        # Invalid percentile (too low)
        with pytest.raises(Exception):
            BayesUCBConfig(initial_p_low=0.3)

    def test_factory_creates_strategy(self):
        """Test factory function creates BayesUCBSelection from config."""
        from TACTICS.thompson_sampling.strategies.config import BayesUCBConfig
        from TACTICS.thompson_sampling.factories import create_strategy

        config = BayesUCBConfig(
            mode="maximize", initial_p_high=0.95, initial_p_low=0.80
        )

        strategy = create_strategy(config)

        assert isinstance(strategy, BayesUCBSelection)
        assert strategy.mode == "maximize"
        assert strategy.p_high == 0.95
        assert strategy.p_low == 0.80


class TestBayesUCBEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_zero_std_reagent(self):
        """Test handling reagent with zero standard deviation."""
        strategy = BayesUCBSelection(mode="maximize")

        reagent = Reagent("test", "CC")
        reagent.mean = 5.0
        reagent.std = 0.0  # Zero std
        reagent.n_samples = 10
        reagent.current_phase = "search"

        # Should not crash, UCB should equal mean
        ucb_indices = strategy._compute_ucb_indices([reagent], percentile=0.95)
        assert np.isclose(ucb_indices[0], reagent.mean)

    def test_very_high_uncertainty(self):
        """Test handling reagent with very high uncertainty."""
        strategy = BayesUCBSelection(mode="maximize")

        reagent = Reagent("test", "CC")
        reagent.mean = 5.0
        reagent.std = 100.0  # Very high std
        reagent.n_samples = 3
        reagent.current_phase = "search"

        ucb_indices = strategy._compute_ucb_indices([reagent], percentile=0.95)

        # UCB should be much larger than mean due to uncertainty
        assert ucb_indices[0] > reagent.mean + 10

    def test_all_reagents_disallowed(self):
        """Test selection when all reagents are disallowed."""
        strategy = BayesUCBSelection(mode="maximize")

        reagents = [Reagent(f"r{i}", "CC") for i in range(3)]
        for r in reagents:
            r.mean = 5.0
            r.std = 1.0
            r.n_samples = 5
            r.current_phase = "search"

        # Disallow all reagents
        disallow_mask = {0, 1, 2}

        # Should still select something (handle gracefully)
        # In this implementation, -inf is assigned to disallowed in maximize mode
        # so np.argmax will select first element (index 0)
        selected_idx = strategy.select_reagent(reagents, disallow_mask=disallow_mask)
        assert selected_idx in {0, 1, 2}
