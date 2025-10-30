"""
Integration tests for warmup strategies with ThompsonSampler.

These tests ensure that any warmup strategy (existing or new) can be properly
integrated with the ThompsonSampler class.
"""

import pytest
import tempfile
import os
import numpy as np
import pandas as pd

from TACTICS.thompson_sampling.core.sampler import ThompsonSampler
from TACTICS.thompson_sampling.core.evaluators import LookupEvaluator
from TACTICS.thompson_sampling.strategies.greedy_selection import GreedySelection
from TACTICS.thompson_sampling.warmup import (
    WarmupStrategy, StandardWarmup, StratifiedWarmup,
    EnhancedWarmup, LatinHypercubeWarmup
)


class TestWarmupIntegration:
    """Test that warmup strategies integrate properly with ThompsonSampler"""

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

    def test_standard_warmup_integration(self):
        """Test that StandardWarmup integrates with ThompsonSampler"""
        warmup_strategy = StandardWarmup()
        selection_strategy = GreedySelection(mode="maximize")

        sampler = ThompsonSampler(
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy
        )
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

        # Run warmup
        warmup_results = sampler.warm_up(num_warmup_trials=2)

        assert len(warmup_results) > 0
        # Standard warmup: each reagent tested num_trials times
        # Expected: (3 + 2) reagents × 2 trials = 10 evaluations
        assert len(warmup_results) >= 8  # Allow some to fail

        # Verify all reagents have been initialized
        for reagent_list in sampler.reagent_lists:
            for reagent in reagent_list:
                assert reagent.mean is not None
                assert reagent.std is not None

        sampler.close()

    def test_stratified_warmup_integration(self):
        """Test that StratifiedWarmup integrates with ThompsonSampler"""
        warmup_strategy = StratifiedWarmup()
        selection_strategy = GreedySelection(mode="maximize")

        sampler = ThompsonSampler(
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy
        )
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

        # Run warmup
        warmup_results = sampler.warm_up(num_warmup_trials=2)

        assert len(warmup_results) > 0

        # Verify all reagents have been initialized
        for reagent_list in sampler.reagent_lists:
            for reagent in reagent_list:
                assert reagent.mean is not None
                assert reagent.std is not None

        sampler.close()

    def test_enhanced_warmup_integration(self):
        """Test that EnhancedWarmup integrates with ThompsonSampler"""
        warmup_strategy = EnhancedWarmup(
            anchor_strategy="max_variance",
            num_anchors=2
        )
        selection_strategy = GreedySelection(mode="maximize")

        sampler = ThompsonSampler(
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy
        )
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

        # Run warmup
        warmup_results = sampler.warm_up(num_warmup_trials=2)

        assert len(warmup_results) > 0

        # Verify all reagents have been initialized
        for reagent_list in sampler.reagent_lists:
            for reagent in reagent_list:
                assert reagent.mean is not None
                assert reagent.std is not None

        sampler.close()

    def test_latin_hypercube_warmup_integration(self):
        """Test that LatinHypercubeWarmup integrates with ThompsonSampler"""
        warmup_strategy = LatinHypercubeWarmup()
        selection_strategy = GreedySelection(mode="maximize")

        sampler = ThompsonSampler(
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy
        )
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

        # Run warmup
        warmup_results = sampler.warm_up(num_warmup_trials=2)

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
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

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
            StratifiedWarmup(),
            EnhancedWarmup(),
            LatinHypercubeWarmup()
        ]

        for warmup_strategy in strategies:
            selection_strategy = GreedySelection(mode="maximize")
            sampler = ThompsonSampler(
                selection_strategy=selection_strategy,
                warmup_strategy=warmup_strategy
            )
            sampler.read_reagents([self.reagent_file1, self.reagent_file2])

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
            StratifiedWarmup(),
            EnhancedWarmup(),
            LatinHypercubeWarmup()
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
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

        warmup_results = sampler.warm_up(num_warmup_trials=2)

        assert len(warmup_results) > 0

        sampler.close()

    def test_warmup_with_parallel_evaluation(self):
        """Test that warmup strategies work with parallel evaluation"""
        warmup_strategies = [
            StandardWarmup(),
            StratifiedWarmup(),
            EnhancedWarmup(),
            LatinHypercubeWarmup()
        ]

        for warmup_strategy in warmup_strategies:
            selection_strategy = GreedySelection(mode="maximize")
            sampler = ThompsonSampler(
                selection_strategy=selection_strategy,
                warmup_strategy=warmup_strategy,
                processes=2  # Enable parallel evaluation
            )
            sampler.read_reagents([self.reagent_file1, self.reagent_file2])
            sampler.set_reaction(self.reaction_smarts)
            sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

            warmup_results = sampler.warm_up(num_warmup_trials=2)

            assert len(warmup_results) > 0

            # Verify all reagents initialized
            for reagent_list in sampler.reagent_lists:
                for reagent in reagent_list:
                    assert reagent.mean is not None

            sampler.close()

    def test_warmup_initializes_priors(self):
        """Test that warmup correctly initializes reagent priors"""
        warmup_strategy = StandardWarmup()
        selection_strategy = GreedySelection(mode="maximize")

        sampler = ThompsonSampler(
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy
        )
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

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
            StratifiedWarmup(),
            EnhancedWarmup(),
            LatinHypercubeWarmup()
        ]

        for warmup_strategy in warmup_strategies:
            selection_strategy = GreedySelection(mode="maximize")
            sampler = ThompsonSampler(
                selection_strategy=selection_strategy,
                warmup_strategy=warmup_strategy
            )
            sampler.read_reagents([self.reagent_file1, self.reagent_file2])
            sampler.set_reaction(self.reaction_smarts)
            sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

            # Run warmup
            warmup_results = sampler.warm_up(num_warmup_trials=2)
            assert len(warmup_results) > 0

            # Run search
            search_results = sampler.search(num_cycles=5)
            assert len(search_results) > 0

            # Verify search results have valid scores
            for result in search_results:
                score, smiles, name = result
                assert isinstance(score, (int, float))
                assert not np.isnan(score)

            sampler.close()

    def test_warmup_handles_failed_evaluations(self):
        """Test that warmup handles failed evaluations gracefully"""
        warmup_strategy = StandardWarmup()
        selection_strategy = GreedySelection(mode="maximize")

        sampler = ThompsonSampler(
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy
        )
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)
        sampler.set_evaluator(LookupEvaluator({"ref_filename": self.lookup_file}))

        # Run warmup - should handle any NaN scores gracefully
        warmup_results = sampler.warm_up(num_warmup_trials=2)

        # Should have some valid results
        assert len(warmup_results) > 0

        # All returned results should have finite scores
        for result in warmup_results:
            score = result[0]
            assert np.isfinite(score)

        sampler.close()
