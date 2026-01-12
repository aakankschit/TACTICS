"""
Integration tests for warmup strategies with ThompsonSampler.

These tests ensure that any warmup strategy (existing or new) can be properly
integrated with the ThompsonSampler class.

Uses temporary test files created in fixtures.
"""

import pytest
import os
import tempfile
import shutil
import numpy as np

from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef
from TACTICS.thompson_sampling.core.sampler import ThompsonSampler
from TACTICS.thompson_sampling.core.evaluators import MWEvaluator
from TACTICS.thompson_sampling.strategies.greedy_selection import GreedySelection
from TACTICS.thompson_sampling.warmup import (
    WarmupStrategy,
    StandardWarmup,
    BalancedWarmup,
    EnhancedWarmup,
)


# Amide coupling reaction SMARTS
REACTION_SMARTS = "[C:1](=O)[OH].[N:2]>>[C:1](=O)[N:2]"


class TestWarmupIntegration:
    """Test that warmup strategies integrate properly with ThompsonSampler"""

    @pytest.fixture(autouse=True)
    def setup_temp_files(self):
        """Create temporary reagent files for testing"""
        self.temp_dir = tempfile.mkdtemp()

        # Create acid reagents file (carboxylic acids)
        self.acids_file = os.path.join(self.temp_dir, "acids.smi")
        with open(self.acids_file, "w") as f:
            f.write("CC(=O)O\tCA0\n")  # acetic acid
            f.write("CCC(=O)O\tCA1\n")  # propionic acid
            f.write("CCCC(=O)O\tCA2\n")  # butyric acid
            f.write("CCCCC(=O)O\tCA3\n")  # valeric acid
            f.write("c1ccccc1C(=O)O\tCA4\n")  # benzoic acid
            f.write("CC(C)C(=O)O\tCA5\n")  # isobutyric acid
            f.write("c1ccc(C(=O)O)cc1\tCA6\n")  # para-benzoic
            f.write("OC(=O)CC(=O)O\tCA7\n")  # malonic acid

        # Create amine reagents file
        self.amines_file = os.path.join(self.temp_dir, "amines.smi")
        with open(self.amines_file, "w") as f:
            f.write("CN\tAA0\n")  # methylamine
            f.write("CCN\tAA1\n")  # ethylamine
            f.write("CCCN\tAA2\n")  # propylamine
            f.write("CCCCN\tAA3\n")  # butylamine
            f.write("c1ccccc1N\tAA4\n")  # aniline
            f.write("CC(C)N\tAA5\n")  # isopropylamine
            f.write("NCCN\tAA6\n")  # ethylenediamine
            f.write("c1ccc(N)cc1\tAA7\n")  # para-aniline

        # Create pipeline
        config = ReactionConfig(
            reactions=[ReactionDef(reaction_smarts=REACTION_SMARTS, step_index=0)],
            reagent_file_list=[self.acids_file, self.amines_file],
        )
        self.pipeline = SynthesisPipeline(config)

        yield

        # Cleanup
        shutil.rmtree(self.temp_dir)

    def test_standard_warmup_integration(self):
        """Test that StandardWarmup integrates with ThompsonSampler"""
        warmup_strategy = StandardWarmup()
        selection_strategy = GreedySelection(mode="maximize")

        sampler = ThompsonSampler(
            self.pipeline,
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy,
        )
        sampler.read_reagents(self.pipeline.reagent_file_list)
        sampler.set_evaluator(MWEvaluator())

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
        warmup_strategy = BalancedWarmup(observations_per_reagent=2)
        selection_strategy = GreedySelection(mode="maximize")

        sampler = ThompsonSampler(
            self.pipeline,
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy,
        )
        sampler.read_reagents(self.pipeline.reagent_file_list)
        sampler.set_evaluator(MWEvaluator())

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
        warmup_strategy = EnhancedWarmup()
        selection_strategy = GreedySelection(mode="maximize")

        sampler = ThompsonSampler(
            self.pipeline,
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy,
        )
        sampler.read_reagents(self.pipeline.reagent_file_list)
        sampler.set_evaluator(MWEvaluator())

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
            observations_per_reagent=3,
            use_per_reagent_variance=True,
            shrinkage_strength=3.0,
        )
        selection_strategy = GreedySelection(mode="maximize")

        sampler = ThompsonSampler(
            self.pipeline,
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy,
        )
        sampler.read_reagents(self.pipeline.reagent_file_list)
        sampler.set_evaluator(MWEvaluator())

        # Run warmup
        warmup_results = sampler.warm_up(num_warmup_trials=3)

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
            self.pipeline,
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy,
        )
        sampler.read_reagents(self.pipeline.reagent_file_list)
        sampler.set_evaluator(MWEvaluator())

        # Generate warmup combinations
        combinations = warmup_strategy.generate_warmup_combinations(
            sampler.reagent_lists,
            num_warmup_trials=2,
            disallow_tracker=sampler._disallow_tracker,
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
            BalancedWarmup(observations_per_reagent=2),
            EnhancedWarmup(),
        ]

        for warmup_strategy in strategies:
            selection_strategy = GreedySelection(mode="maximize")
            sampler = ThompsonSampler(
                self.pipeline,
                selection_strategy=selection_strategy,
                warmup_strategy=warmup_strategy,
            )
            sampler.read_reagents(self.pipeline.reagent_file_list)

            # Get expected evaluations
            expected = warmup_strategy.get_expected_evaluations(
                sampler.reagent_lists, num_warmup_trials=2
            )

            assert isinstance(expected, int)
            assert expected > 0

            sampler.close()

    def test_warmup_strategy_names(self):
        """Test that warmup strategies have descriptive names"""
        strategies = [
            StandardWarmup(),
            BalancedWarmup(observations_per_reagent=2),
            EnhancedWarmup(),
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

            def generate_warmup_combinations(
                self, reagent_lists, num_warmup_trials, disallow_tracker
            ):
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
            self.pipeline,
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy,
        )
        sampler.read_reagents(self.pipeline.reagent_file_list)
        sampler.set_evaluator(MWEvaluator())

        warmup_results = sampler.warm_up(num_warmup_trials=2)

        assert len(warmup_results) > 0

        sampler.close()

    def test_warmup_initializes_priors(self):
        """Test that warmup correctly initializes reagent priors"""
        warmup_strategy = StandardWarmup()
        selection_strategy = GreedySelection(mode="maximize")

        sampler = ThompsonSampler(
            self.pipeline,
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy,
        )
        sampler.read_reagents(self.pipeline.reagent_file_list)
        sampler.set_evaluator(MWEvaluator())

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
            BalancedWarmup(observations_per_reagent=2),
            EnhancedWarmup(),
        ]

        for warmup_strategy in warmup_strategies:
            selection_strategy = GreedySelection(mode="maximize")
            sampler = ThompsonSampler(
                self.pipeline,
                selection_strategy=selection_strategy,
                warmup_strategy=warmup_strategy,
            )
            sampler.read_reagents(self.pipeline.reagent_file_list)
            sampler.set_evaluator(MWEvaluator())

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
            self.pipeline,
            selection_strategy=selection_strategy,
            warmup_strategy=warmup_strategy,
        )
        sampler.read_reagents(self.pipeline.reagent_file_list)
        sampler.set_evaluator(MWEvaluator())

        # Run warmup - should handle any NaN scores gracefully (returns polars DataFrame)
        warmup_df = sampler.warm_up(num_warmup_trials=2)

        # Should have some valid results
        assert len(warmup_df) > 0

        # All returned results should have finite scores
        scores = warmup_df["score"].to_list()
        for score in scores:
            assert np.isfinite(score)

        sampler.close()
