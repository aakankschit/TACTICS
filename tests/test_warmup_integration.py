"""
Integration tests for warmup strategies with ThompsonSampler.

These tests ensure that any warmup strategy (existing or new) can be properly
integrated with the ThompsonSampler class.

Uses real example data from examples/ folder.
"""

import pytest
import os
import numpy as np

from TACTICS.thompson_sampling.core.sampler import ThompsonSampler
from TACTICS.thompson_sampling.core.evaluators import LookupEvaluator
from TACTICS.thompson_sampling.strategies.greedy_selection import GreedySelection
from TACTICS.thompson_sampling.warmup import (
    WarmupStrategy, StandardWarmup, BalancedWarmup,
    EnhancedWarmup
)


# Paths to real example data
EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "..", "examples")
REAGENT_FILE1 = os.path.join(EXAMPLES_DIR, "input_files", "acids.smi")
REAGENT_FILE2 = os.path.join(EXAMPLES_DIR, "input_files", "coupled_aa_sub.smi")
LOOKUP_FILE = os.path.join(EXAMPLES_DIR, "docking_scores", "product_scores.csv")
REACTION_SMARTS = "[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]>>[#6:1](=[O:2])[#7:3]"


class TestWarmupIntegration:
    """Test that warmup strategies integrate properly with ThompsonSampler"""

    def test_standard_warmup_integration(self):
        """Test that StandardWarmup integrates with ThompsonSampler"""
        warmup_strategy = StandardWarmup()
        selection_strategy = GreedySelection(mode="maximize")

        sampler = ThompsonSampler(
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy
        )
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

        # Run warmup
        warmup_results = sampler.warm_up(num_warmup_trials=2)

        assert len(warmup_results) > 0

        # Verify all reagents have been initialized
        for reagent_list in sampler.reagent_lists:
            for reagent in reagent_list:
                assert reagent.mean is not None
                assert reagent.std is not None

        sampler.close()

    def test_balanced_warmup_integration(self):
        """Test that BalancedWarmup integrates with ThompsonSampler"""
        warmup_strategy = BalancedWarmup(observations_per_reagent=3)
        selection_strategy = GreedySelection(mode="maximize")

        sampler = ThompsonSampler(
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy
        )
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

        # Run warmup
        warmup_results = sampler.warm_up(num_warmup_trials=3)

        assert len(warmup_results) > 0

        # Verify all reagents have been initialized
        for reagent_list in sampler.reagent_lists:
            for reagent in reagent_list:
                assert reagent.mean is not None
                assert reagent.std is not None

        sampler.close()

    def test_enhanced_warmup_integration(self):
        """Test that EnhancedWarmup integrates with ThompsonSampler"""
        warmup_strategy = EnhancedWarmup()
        selection_strategy = GreedySelection(mode="maximize")

        sampler = ThompsonSampler(
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy
        )
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

        # Run warmup
        warmup_results = sampler.warm_up(num_warmup_trials=2)

        assert len(warmup_results) > 0

        # Verify all reagents have been initialized
        for reagent_list in sampler.reagent_lists:
            for reagent in reagent_list:
                assert reagent.mean is not None
                assert reagent.std is not None

        sampler.close()

    def test_balanced_warmup_per_reagent_variance(self):
        """Test that BalancedWarmup with per-reagent variance integrates with ThompsonSampler"""
        warmup_strategy = BalancedWarmup(
            observations_per_reagent=5,
            use_per_reagent_variance=True,
            shrinkage_strength=3.0
        )
        selection_strategy = GreedySelection(mode="maximize")

        sampler = ThompsonSampler(
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy
        )
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

        # Run warmup
        warmup_results = sampler.warm_up(num_warmup_trials=5)

        assert len(warmup_results) > 0

        # Verify all reagents have been initialized
        for reagent_list in sampler.reagent_lists:
            for reagent in reagent_list:
                assert reagent.mean is not None
                assert reagent.std is not None

        sampler.close()

    def test_warmup_generates_correct_combinations(self):
        """Test that warmup generates valid reagent combinations"""
        warmup_strategy = StandardWarmup()
        selection_strategy = GreedySelection(mode="maximize")

        sampler = ThompsonSampler(
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy
        )
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

        # Generate warmup combinations
        combinations = warmup_strategy.generate_warmup_combinations(
            sampler.reagent_lists,
            num_warmup_trials=2,
            disallow_tracker=sampler._disallow_tracker
        )

        # Verify combinations are valid
        assert len(combinations) > 0

        for combo in combinations:
            assert len(combo) == len(sampler.reagent_lists)
            for component_idx, reagent_idx in enumerate(combo):
                assert 0 <= reagent_idx < len(sampler.reagent_lists[component_idx])

        sampler.close()

    def test_warmup_expected_evaluations(self):
        """Test that warmup strategies report correct expected evaluations"""
        strategies = [
            StandardWarmup(),
            BalancedWarmup(observations_per_reagent=3),
            EnhancedWarmup()
        ]

        for warmup_strategy in strategies:
            selection_strategy = GreedySelection(mode="maximize")
            sampler = ThompsonSampler(
                selection_strategy=selection_strategy,
                warmup_strategy=warmup_strategy
            )
            sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])

            # Get expected evaluations
            expected = warmup_strategy.get_expected_evaluations(
                sampler.reagent_lists,
                num_warmup_trials=2
            )

            assert isinstance(expected, int)
            assert expected > 0

            sampler.close()

    def test_warmup_strategy_names(self):
        """Test that warmup strategies have descriptive names"""
        strategies = [
            StandardWarmup(),
            BalancedWarmup(observations_per_reagent=3),
            EnhancedWarmup()
        ]

        for warmup_strategy in strategies:
            name = warmup_strategy.get_name()
            assert isinstance(name, str)
            assert len(name) > 0

            description = warmup_strategy.get_description()
            assert isinstance(description, str)
            assert len(description) > 0

    def test_custom_warmup_integration(self):
        """Test that a custom warmup strategy can be integrated"""

        class CustomWarmup(WarmupStrategy):
            """Simple custom warmup that randomly samples combinations"""

            def generate_warmup_combinations(self, reagent_lists, num_warmup_trials,
                                           disallow_tracker):
                import random
                combinations = []

                # Generate random combinations
                for _ in range(num_warmup_trials * 2):
                    combo = [
                        random.randint(0, len(reagent_list) - 1)
                        for reagent_list in reagent_lists
                    ]
                    combinations.append(combo)

                return combinations

            def get_name(self):
                return "Custom Random Warmup"

        # Test custom warmup
        warmup_strategy = CustomWarmup()
        selection_strategy = GreedySelection(mode="maximize")

        sampler = ThompsonSampler(
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy
        )
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

        warmup_results = sampler.warm_up(num_warmup_trials=2)

        assert len(warmup_results) > 0

        sampler.close()

    def test_warmup_initializes_priors(self):
        """Test that warmup correctly initializes reagent priors"""
        warmup_strategy = StandardWarmup()
        selection_strategy = GreedySelection(mode="maximize")

        sampler = ThompsonSampler(
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy
        )
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

        # Before warmup, reagents should not have initialized priors
        for reagent_list in sampler.reagent_lists:
            for reagent in reagent_list:
                assert reagent.mean is None or reagent.mean == 0

        # Run warmup
        warmup_results = sampler.warm_up(num_warmup_trials=2)

        # After warmup, all reagents should have valid priors
        for reagent_list in sampler.reagent_lists:
            for reagent in reagent_list:
                assert reagent.mean is not None
                assert reagent.std is not None
                assert reagent.std > 0  # Standard deviation should be positive

        sampler.close()

    def test_warmup_to_search_workflow(self):
        """Test complete workflow from warmup to search"""
        warmup_strategies = [
            StandardWarmup(),
            BalancedWarmup(observations_per_reagent=3),
            EnhancedWarmup()
        ]

        for warmup_strategy in warmup_strategies:
            selection_strategy = GreedySelection(mode="maximize")
            sampler = ThompsonSampler(
                selection_strategy=selection_strategy,
                warmup_strategy=warmup_strategy
            )
            sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
            sampler.set_reaction(REACTION_SMARTS)
            sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

            # Run warmup (returns polars DataFrame)
            warmup_df = sampler.warm_up(num_warmup_trials=2)
            assert len(warmup_df) > 0

            # Run search (returns polars DataFrame)
            search_df = sampler.search(num_cycles=5)
            assert len(search_df) > 0

            # Verify search results have valid scores
            scores = search_df["score"].to_list()
            smiles_list = search_df["SMILES"].to_list()
            names_list = search_df["Name"].to_list()

            for score, smiles, name in zip(scores, smiles_list, names_list):
                assert isinstance(score, (int, float))
                assert np.isfinite(score)

            sampler.close()

    def test_warmup_handles_failed_evaluations(self):
        """Test that warmup handles failed evaluations gracefully"""
        warmup_strategy = StandardWarmup()
        selection_strategy = GreedySelection(mode="maximize")

        sampler = ThompsonSampler(
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy
        )
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": LOOKUP_FILE}))

        # Run warmup - should handle any NaN scores gracefully (returns polars DataFrame)
        warmup_df = sampler.warm_up(num_warmup_trials=2)

        # Should have some valid results
        assert len(warmup_df) > 0

        # All returned results should have finite scores
        scores = warmup_df["score"].to_list()
        for score in scores:
            assert np.isfinite(score)

        sampler.close()
