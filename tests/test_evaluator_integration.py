"""
Integration tests for evaluator classes with ThompsonSampler.

These tests ensure that any evaluator (existing or new) can be properly
integrated with the ThompsonSampler class.

Uses real example data from examples/ folder.
"""

import pytest
import os
import numpy as np

from TACTICS.thompson_sampling.core.sampler import ThompsonSampler
from TACTICS.thompson_sampling.core.evaluators import (
    Evaluator, MWEvaluator, FPEvaluator, LookupEvaluator
)
from TACTICS.thompson_sampling.strategies.greedy_selection import GreedySelection


# Paths to real example data
EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "..", "examples")
REAGENT_FILE1 = os.path.join(EXAMPLES_DIR, "input_files", "acids.smi")
REAGENT_FILE2 = os.path.join(EXAMPLES_DIR, "input_files", "coupled_aa_sub.smi")
LOOKUP_FILE = os.path.join(EXAMPLES_DIR, "docking_scores", "product_scores.csv")
REACTION_SMARTS = "[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]>>[#6:1](=[O:2])[#7:3]"


class TestEvaluatorIntegration:
    """Test that evaluators integrate properly with ThompsonSampler"""

    def test_lookup_evaluator_integration(self):
        """Test that LookupEvaluator integrates with ThompsonSampler"""
        # Create sampler with LookupEvaluator
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)

        evaluator = LookupEvaluator({"ref_filename": LOOKUP_FILE})
        sampler.set_evaluator(evaluator)

        # Verify evaluator is set correctly
        assert sampler.evaluator is not None
        assert isinstance(sampler.evaluator, LookupEvaluator)
        assert sampler.evaluator.counter == 0

        # Test evaluation
        smiles, name, score = sampler.evaluate([0, 0])
        assert isinstance(score, (int, float))
        assert sampler.evaluator.counter == 1

        sampler.close()

    def test_evaluator_counter_increments(self):
        """Test that evaluator counter increments correctly"""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy, batch_size=5)
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)

        evaluator = LookupEvaluator({"ref_filename": LOOKUP_FILE})
        sampler.set_evaluator(evaluator)

        initial_count = evaluator.counter

        # Evaluate batch of compounds
        batch_results = sampler.evaluate_batch([[0, 0], [0, 1], [1, 0]])

        # Counter should increment by number of evaluations
        assert evaluator.counter == initial_count + 3
        assert len(batch_results) == 3

        sampler.close()

    def test_evaluator_with_warmup(self):
        """Test that evaluators work correctly during warmup phase"""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)

        evaluator = LookupEvaluator({"ref_filename": LOOKUP_FILE})
        sampler.set_evaluator(evaluator)

        # Run warmup (returns polars DataFrame)
        warmup_df = sampler.warm_up(num_warmup_trials=2)

        # Verify warmup ran and evaluator was called
        assert len(warmup_df) > 0
        assert evaluator.counter > 0

        # All warmup results should have valid scores
        scores = warmup_df["score"].to_list()
        smiles_list = warmup_df["SMILES"].to_list()
        names_list = warmup_df["Name"].to_list()

        for score, smiles, name in zip(scores, smiles_list, names_list):
            assert isinstance(score, (int, float))
            assert np.isfinite(score)

        sampler.close()

    def test_evaluator_with_full_workflow(self):
        """Test evaluator through complete warmup + search workflow"""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy, batch_size=1)
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)

        evaluator = LookupEvaluator({"ref_filename": LOOKUP_FILE})
        sampler.set_evaluator(evaluator)

        # Run warmup
        warmup_results = sampler.warm_up(num_warmup_trials=2)
        warmup_count = evaluator.counter

        # Run search
        search_results = sampler.search(num_cycles=5)

        # Verify evaluator was called during both phases
        assert warmup_count > 0
        assert evaluator.counter > warmup_count
        assert len(search_results) > 0

        sampler.close()

    def test_mw_evaluator_integration(self):
        """Test that MWEvaluator integrates with ThompsonSampler (no lookup needed)"""
        # Create sampler with MWEvaluator
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)

        evaluator = MWEvaluator()
        sampler.set_evaluator(evaluator)

        # Verify evaluator is set correctly
        assert sampler.evaluator is not None
        assert isinstance(sampler.evaluator, MWEvaluator)
        assert sampler.evaluator.counter == 0

        # Test single evaluation
        smiles, name, score = sampler.evaluate([0, 0])
        assert isinstance(score, (int, float))
        assert not np.isnan(score)
        assert sampler.evaluator.counter == 1

        sampler.close()

    def test_fp_evaluator_integration(self):
        """Test that FPEvaluator integrates with ThompsonSampler (no lookup needed)"""
        # Create sampler with FPEvaluator
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)

        # Use a simple target molecule
        evaluator = FPEvaluator({"query_smiles": "CC(=O)NC1CCCCC1"})
        sampler.set_evaluator(evaluator)

        # Verify evaluator is set correctly
        assert sampler.evaluator is not None
        assert isinstance(sampler.evaluator, FPEvaluator)
        assert sampler.evaluator.counter == 0

        # Test evaluation
        smiles, name, score = sampler.evaluate([0, 0])
        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0  # Tanimoto similarity
        assert sampler.evaluator.counter == 1

        sampler.close()

    def test_custom_evaluator_integration(self):
        """Test that a custom evaluator can be integrated"""

        class CustomEvaluator(Evaluator):
            """Simple custom evaluator for testing"""

            def __init__(self):
                self.num_evaluations = 0

            @property
            def counter(self):
                return self.num_evaluations

            def evaluate(self, mol):
                self.num_evaluations += 1
                # Return number of atoms as score
                return mol.GetNumAtoms()

        # Create sampler with custom evaluator
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([REAGENT_FILE1, REAGENT_FILE2])
        sampler.set_reaction(REACTION_SMARTS)

        evaluator = CustomEvaluator()
        sampler.set_evaluator(evaluator)

        # Test evaluation
        smiles, name, score = sampler.evaluate([0, 0])
        assert isinstance(score, (int, float))
        assert score > 0
        assert evaluator.counter == 1

        sampler.close()
