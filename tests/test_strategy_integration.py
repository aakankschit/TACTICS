"""
Integration tests for selection strategies with ThompsonSampler.

These tests ensure that any selection strategy (existing or new) can be properly
integrated with the ThompsonSampler class.
"""

import pytest
import tempfile
import os
import numpy as np
import pandas as pd

from TACTICS.thompson_sampling.core.sampler import ThompsonSampler
from TACTICS.thompson_sampling.core.evaluators import LookupEvaluator
from TACTICS.thompson_sampling.strategies.base_strategy import SelectionStrategy
from TACTICS.thompson_sampling.strategies.greedy_selection import GreedySelection
from TACTICS.thompson_sampling.strategies.roulette_wheel import RouletteWheelSelection
from TACTICS.thompson_sampling.strategies.ucb_selection import UCBSelection
from TACTICS.thompson_sampling.strategies.epsilon_greedy import EpsilonGreedySelection


class TestStrategyIntegration:
    """Test that selection strategies integrate properly with ThompsonSampler"""

    def setup_method(self):
        """Set up test fixtures"""
        self.temp_dir = tempfile.mkdtemp()

        # Create test reagent files
        self.reagent_file1 = os.path.join(self.temp_dir, "reagents1.smi")
        self.reagent_file2 = os.path.join(self.temp_dir, "reagents2.smi")

        with open(self.reagent_file1, "w") as f:
            f.write("CCO\tethanol\n")
            f.write("CCCO\tpropanol\n")
            f.write("CCCCO\tbutanol\n")

        with open(self.reagent_file2, "w") as f:
            f.write("CC(=O)O\tacetic_acid\n")
            f.write("CCC(=O)O\tpropionic_acid\n")

        # Create test lookup file
        self.lookup_file = os.path.join(self.temp_dir, "scores.csv")
        lookup_data = pd.DataFrame({
            'Product_Code': [
                'ethanol_acetic_acid', 'ethanol_propionic_acid',
                'propanol_acetic_acid', 'propanol_propionic_acid',
                'butanol_acetic_acid', 'butanol_propionic_acid'
            ],
            'Scores': [0.5, 0.7, 0.6, 0.8, 0.4, 0.9]
        })
        lookup_data.to_csv(self.lookup_file, index=False)

        # Simple reaction SMARTS
        self.reaction_smarts = "[C:1][OH:2].[C:3](=[O:4])[OH:5]>>[C:1][O:2][C:3](=[O:4])"

    def teardown_method(self):
        """Clean up test fixtures"""
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_greedy_selection_integration(self):
        """Test that GreedySelection integrates with ThompsonSampler"""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

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
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

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
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

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
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

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
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

        # Run warmup and search
        sampler.warm_up(num_warmup_trials=2)
        results = sampler.search(num_cycles=10)

        # Verify results contain scores
        scores = [r[0] for r in results]
        assert len(scores) > 0
        assert all(isinstance(s, (int, float)) for s in scores)

        sampler.close()

    def test_strategy_mode_minimize(self):
        """Test that strategies work correctly in minimize mode"""
        strategy = GreedySelection(mode="minimize")
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

        # Run warmup and search
        sampler.warm_up(num_warmup_trials=2)
        results = sampler.search(num_cycles=10)

        # Verify results contain scores
        scores = [r[0] for r in results]
        assert len(scores) > 0
        assert all(isinstance(s, (int, float)) for s in scores)

        sampler.close()

    def test_batch_selection_single_mode(self):
        """Test batch_size=1 (single selection per cycle)"""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy, batch_size=1)
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

        sampler.warm_up(num_warmup_trials=2)
        results = sampler.search(num_cycles=5)

        assert len(results) > 0

        sampler.close()

    def test_batch_selection_batch_mode(self):
        """Test batch_size>1 (multiple selections per cycle)"""
        strategy = RouletteWheelSelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy, batch_size=3)
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

        sampler.warm_up(num_warmup_trials=2)
        results = sampler.search(num_cycles=5)

        assert len(results) > 0

        sampler.close()

    def test_strategy_with_disallow_mask(self):
        """Test that strategies respect disallow mask"""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

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
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

        warmup_results = sampler.warm_up(num_warmup_trials=2)
        search_results = sampler.search(num_cycles=5)

        assert len(warmup_results) > 0
        assert len(search_results) > 0

        sampler.close()

    def test_strategy_select_batch_method(self):
        """Test that select_batch method works for strategies that implement it"""
        strategy = RouletteWheelSelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy, batch_size=5)
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

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

    def test_all_strategies_with_parallel_evaluation(self):
        """Test that all strategies work with parallel evaluation"""
        strategies = [
            GreedySelection(mode="maximize"),
            RouletteWheelSelection(mode="maximize"),
            UCBSelection(mode="maximize"),
            EpsilonGreedySelection(mode="maximize")
        ]

        for strategy in strategies:
            sampler = ThompsonSampler(
                selection_strategy=strategy,
                processes=2,  # Enable parallel evaluation
                batch_size=2
            )
            sampler.read_reagents([self.reagent_file1, self.reagent_file2])
            sampler.set_reaction(self.reaction_smarts)
            sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

            warmup_results = sampler.warm_up(num_warmup_trials=2)
            search_results = sampler.search(num_cycles=3)

            assert len(warmup_results) > 0
            assert len(search_results) > 0

            sampler.close()

    def test_strategy_with_max_resamples(self):
        """Test that strategies respect max_resamples stopping criterion"""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            selection_strategy=strategy,
            batch_size=1,
            max_resamples=5  # Stop after 5 consecutive duplicate samples
        )
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

        sampler.warm_up(num_warmup_trials=2)

        # With greedy selection and small library, should hit max_resamples quickly
        results = sampler.search(num_cycles=100)

        # Should stop early due to max_resamples
        assert len(results) < 100

        sampler.close()
