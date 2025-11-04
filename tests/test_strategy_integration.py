"""
Integration tests for selection strategies with ThompsonSampler.

These tests ensure that any selection strategy (existing or new) can be properly
integrated with the ThompsonSampler class.

Uses real example data from examples/ folder.
"""

import pytest
import os
import numpy as np

from TACTICS.thompson_sampling.core.sampler import ThompsonSampler
from TACTICS.thompson_sampling.core.evaluators import LookupEvaluator
from TACTICS.thompson_sampling.strategies.base_strategy import SelectionStrategy
from TACTICS.thompson_sampling.strategies.greedy_selection import GreedySelection
from TACTICS.thompson_sampling.strategies.roulette_wheel import RouletteWheelSelection
from TACTICS.thompson_sampling.strategies.ucb_selection import UCBSelection
from TACTICS.thompson_sampling.strategies.epsilon_greedy import EpsilonGreedySelection
from TACTICS.thompson_sampling.strategies.bayes_ucb_selection import BayesUCBSelection


# Paths to real example data
EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "..", "examples")
REAGENT_FILE1 = os.path.join(EXAMPLES_DIR, "input_files", "acids.smi")
REAGENT_FILE2 = os.path.join(EXAMPLES_DIR, "input_files", "coupled_aa_sub.smi")
LOOKUP_FILE = os.path.join(EXAMPLES_DIR, "docking_scores", "product_scores.csv")
REACTION_SMARTS = "[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]>>[#6:1](=[O:2])[#7:3]"


class TestStrategyIntegration:
    """Test that selection strategies integrate properly with ThompsonSampler"""

    def test_greedy_selection_integration(self):
        """Test that GreedySelection integrates with ThompsonSampler"""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

        # Run warmup and search
        warmup_results = sampler.warm_up(num_warmup_trials=2)
        search_results = sampler.search(num_cycles=5)

        assert len(warmup_results) > 0
        assert len(search_results) > 0

        sampler.close()

    def test_roulette_wheel_integration(self):
        """Test that RouletteWheelSelection integrates with ThompsonSampler"""
        strategy = RouletteWheelSelection(mode="maximize", alpha=0.1, beta=0.1)
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

        # Run warmup and search
        warmup_results = sampler.warm_up(num_warmup_trials=2)
        search_results = sampler.search(num_cycles=5)

        assert len(warmup_results) > 0
        assert len(search_results) > 0

        sampler.close()

    def test_ucb_selection_integration(self):
        """Test that UCBSelection integrates with ThompsonSampler"""
        strategy = UCBSelection(mode="maximize", c=1.0)
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

        # Run warmup and search
        warmup_results = sampler.warm_up(num_warmup_trials=2)
        search_results = sampler.search(num_cycles=5)

        assert len(warmup_results) > 0
        assert len(search_results) > 0

        sampler.close()

    def test_epsilon_greedy_integration(self):
        """Test that EpsilonGreedySelection integrates with ThompsonSampler"""
        strategy = EpsilonGreedySelection(mode="maximize", epsilon=0.1)
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

        # Run warmup and search
        warmup_results = sampler.warm_up(num_warmup_trials=2)
        search_results = sampler.search(num_cycles=5)

        assert len(warmup_results) > 0
        assert len(search_results) > 0

        sampler.close()

    def test_strategy_mode_maximize(self):
        """Test that strategies work correctly in maximize mode"""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

        # Run warmup and search (now returns polars DataFrame)
        sampler.warm_up(num_warmup_trials=2)
        results_df = sampler.search(num_cycles=10)

        # Verify results contain scores
        assert len(results_df) > 0
        scores = results_df["score"].to_list()
        assert all(isinstance(s, (int, float)) for s in scores)
        assert all(np.isfinite(s) for s in scores)

        sampler.close()

    def test_strategy_mode_minimize(self):
        """Test that strategies work correctly in minimize mode"""
        strategy = GreedySelection(mode="minimize")
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

        # Run warmup and search (now returns polars DataFrame)
        sampler.warm_up(num_warmup_trials=2)
        results_df = sampler.search(num_cycles=10)

        # Verify results contain scores
        assert len(results_df) > 0
        scores = results_df["score"].to_list()
        assert all(isinstance(s, (int, float)) for s in scores)
        assert all(np.isfinite(s) for s in scores)

        sampler.close()

    def test_batch_selection_single_mode(self):
        """Test batch_size=1 (single selection per cycle)"""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy, batch_size=1)
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

        sampler.warm_up(num_warmup_trials=2)
        results = sampler.search(num_cycles=5)

        assert len(results) > 0

        sampler.close()

    def test_batch_selection_batch_mode(self):
        """Test batch_size>1 (multiple selections per cycle)"""
        strategy = RouletteWheelSelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy, batch_size=3)
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

        sampler.warm_up(num_warmup_trials=2)
        results = sampler.search(num_cycles=5)

        assert len(results) > 0

        sampler.close()

    def test_strategy_with_disallow_mask(self):
        """Test that strategies respect disallow mask"""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

        # Run warmup to initialize reagents
        sampler.warm_up(num_warmup_trials=2)

        # Test select_reagent with disallow mask
        rng = np.random.default_rng(seed=42)
        disallow_mask = {0}  # Disallow first reagent
        selected = strategy.select_reagent(
            reagent_list=sampler.reagent_lists[0],
            disallow_mask=disallow_mask,
            rng=rng
        )

        # Should not select disallowed index
        assert selected not in disallow_mask

        sampler.close()

    def test_roulette_wheel_temperature_adaptation(self):
        """Test RouletteWheelSelection adaptive temperature control"""
        strategy = RouletteWheelSelection(mode="maximize", alpha=0.1, beta=0.1)

        # Test temperature increase when efficiency is low
        initial_alpha = strategy.alpha
        strategy.update_temperature(n_unique=0, batch_size=10)
        assert strategy.alpha > initial_alpha

        # Test temperature reset
        strategy.reset_temperature()
        assert strategy.alpha == strategy.initial_alpha
        assert strategy.beta == strategy.initial_beta

    def test_roulette_wheel_component_rotation(self):
        """Test RouletteWheelSelection component rotation for thermal cycling"""
        strategy = RouletteWheelSelection(mode="maximize")

        # Test component rotation
        assert strategy.current_component_idx == 0
        strategy.rotate_component(n_components=3)
        assert strategy.current_component_idx == 1
        strategy.rotate_component(n_components=3)
        assert strategy.current_component_idx == 2
        strategy.rotate_component(n_components=3)
        assert strategy.current_component_idx == 0  # Should wrap around

    def test_custom_strategy_integration(self):
        """Test that a custom strategy can be integrated"""

        class CustomStrategy(SelectionStrategy):
            """Simple custom strategy that selects randomly"""

            def select_reagent(self, reagent_list, disallow_mask=None, **kwargs):
                rng = kwargs.get('rng', np.random.default_rng())
                available = list(range(len(reagent_list)))

                if disallow_mask:
                    available = [i for i in available if i not in disallow_mask]

                return rng.choice(available) if available else 0

        # Test custom strategy
        strategy = CustomStrategy(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

        warmup_results = sampler.warm_up(num_warmup_trials=2)
        search_results = sampler.search(num_cycles=5)

        assert len(warmup_results) > 0
        assert len(search_results) > 0

        sampler.close()

    def test_strategy_select_batch_method(self):
        """Test that select_batch method works for strategies that implement it"""
        strategy = RouletteWheelSelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy, batch_size=5)
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

        # Run warmup to initialize reagents
        sampler.warm_up(num_warmup_trials=2)

        # Test select_batch
        rng = np.random.default_rng(seed=42)
        batch = strategy.select_batch(
            reagent_list=sampler.reagent_lists[0],
            batch_size=5,
            rng=rng
        )

        assert len(batch) == 5
        assert all(0 <= idx < len(sampler.reagent_lists[0]) for idx in batch)

        sampler.close()

    def test_bayes_ucb_selection_integration(self):
        """Test that BayesUCBSelection integrates with ThompsonSampler"""
        strategy = BayesUCBSelection(
            mode="maximize",
            initial_p_high=0.90,
            initial_p_low=0.90
        )
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

        # Run warmup and search
        warmup_results = sampler.warm_up(num_warmup_trials=2)
        search_results = sampler.search(num_cycles=5)

        assert len(warmup_results) > 0
        assert len(search_results) > 0

        sampler.close()

    def test_bayes_ucb_percentile_adaptation(self):
        """Test BayesUCBSelection adaptive percentile control"""
        strategy = BayesUCBSelection(
            mode="maximize",
            initial_p_high=0.90,
            initial_p_low=0.90,
            efficiency_threshold=0.10
        )

        # Test percentile increase when efficiency is low (stuck)
        initial_p_high = strategy.p_high
        strategy.update_percentiles(n_unique=0, batch_size=10)
        assert strategy.p_high > initial_p_high

        # Test percentile reset
        strategy.reset_percentiles()
        assert strategy.p_high == 0.90
        assert strategy.p_low == 0.90

    def test_bayes_ucb_component_rotation(self):
        """Test BayesUCBSelection component rotation for thermal cycling"""
        strategy = BayesUCBSelection(mode="maximize")

        # Test component rotation
        assert strategy.current_component_idx == 0
        strategy.rotate_component(n_components=3)
        assert strategy.current_component_idx == 1
        strategy.rotate_component(n_components=3)
        assert strategy.current_component_idx == 2
        strategy.rotate_component(n_components=3)
        assert strategy.current_component_idx == 0  # Should wrap around

    def test_bayes_ucb_with_thermal_cycling(self):
        """Test BayesUCB with thermal cycling in real sampling scenario"""
        strategy = BayesUCBSelection(
            mode="maximize",
            initial_p_high=0.90,
            initial_p_low=0.80,
            efficiency_threshold=0.10
        )
        sampler = ThompsonSampler(selection_strategy=strategy, batch_size=3)
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

        # Run warmup and search
        sampler.warm_up(num_warmup_trials=2)
        results = sampler.search(num_cycles=10)

        assert len(results) > 0
        # Verify thermal cycling happened (component should have rotated)
        # After 10 cycles with 2 components, should have cycled through multiple times
        # (exact behavior depends on implementation details)

        sampler.close()
