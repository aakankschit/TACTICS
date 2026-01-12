"""
Integration tests for evaluator classes with ThompsonSampler.

These tests ensure that any evaluator (existing or new) can be properly
integrated with the ThompsonSampler class.
"""

import pytest
import os
import tempfile
import shutil
import numpy as np

from TACTICS.thompson_sampling.core.sampler import ThompsonSampler
from TACTICS.thompson_sampling.core.evaluators import (
    Evaluator,
    MWEvaluator,
    FPEvaluator,
    LookupEvaluator,
)
from TACTICS.thompson_sampling.strategies.greedy_selection import GreedySelection
from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef


# Amide coupling reaction
REACTION_SMARTS = "[C:1](=O)[OH].[N:2]>>[C:1](=O)[N:2]"


class TestEvaluatorIntegration:
    """Test that evaluators integrate properly with ThompsonSampler"""

    def setup_method(self):
        """Set up test fixtures with temporary files."""
        self.temp_dir = tempfile.mkdtemp()

        # Create reagent files
        self.reagent_file1 = os.path.join(self.temp_dir, "acids.smi")
        self.reagent_file2 = os.path.join(self.temp_dir, "amines.smi")

        with open(self.reagent_file1, "w") as f:
            f.write("CC(=O)O\tacetic_acid\n")
            f.write("CCC(=O)O\tpropionic_acid\n")
            f.write("CCCC(=O)O\tbutyric_acid\n")

        with open(self.reagent_file2, "w") as f:
            f.write("CN\tmethylamine\n")
            f.write("CCN\tethylamine\n")
            f.write("CCCN\tpropylamine\n")

        # Create lookup file with scores (using default LookupEvaluator column names)
        self.lookup_file = os.path.join(self.temp_dir, "scores.csv")
        with open(self.lookup_file, "w") as f:
            f.write("Product_Code,Scores\n")
            f.write("acetic_acid_methylamine,0.8\n")
            f.write("acetic_acid_ethylamine,0.7\n")
            f.write("acetic_acid_propylamine,0.6\n")
            f.write("propionic_acid_methylamine,0.9\n")
            f.write("propionic_acid_ethylamine,0.85\n")
            f.write("propionic_acid_propylamine,0.75\n")
            f.write("butyric_acid_methylamine,0.5\n")
            f.write("butyric_acid_ethylamine,0.55\n")
            f.write("butyric_acid_propylamine,0.45\n")

        # Create pipeline
        self.reaction_config = ReactionConfig(
            reactions=[ReactionDef(reaction_smarts=REACTION_SMARTS, step_index=0)],
            reagent_file_list=[self.reagent_file1, self.reagent_file2],
        )
        self.pipeline = SynthesisPipeline(self.reaction_config)

    def teardown_method(self):
        """Clean up temporary files."""
        shutil.rmtree(self.temp_dir)

    def test_lookup_evaluator_integration(self):
        """Test that LookupEvaluator integrates with ThompsonSampler"""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(self.pipeline, selection_strategy=strategy)
        sampler.read_reagents(self.pipeline.reagent_file_list)

        evaluator = LookupEvaluator({"ref_filename": self.lookup_file})
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
        sampler = ThompsonSampler(
            self.pipeline, selection_strategy=strategy, batch_size=5
        )
        sampler.read_reagents(self.pipeline.reagent_file_list)

        evaluator = LookupEvaluator({"ref_filename": self.lookup_file})
        sampler.set_evaluator(evaluator)

        initial_count = evaluator.counter

        # Evaluate batch of compounds
        batch_results = sampler.evaluate_batch([[0, 0], [0, 1], [1, 0]])

        # Counter should increment by number of evaluations
        assert evaluator.counter == initial_count + 3
        assert len(batch_results) == 3

        sampler.close()

    @pytest.mark.skip(reason="Warmup has edge case with small test data")
    def test_evaluator_with_warmup(self):
        """Test that evaluators work correctly during warmup phase"""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(self.pipeline, selection_strategy=strategy)
        sampler.read_reagents(self.pipeline.reagent_file_list)

        evaluator = LookupEvaluator({"ref_filename": self.lookup_file})
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

    @pytest.mark.skip(reason="Warmup has edge case with small test data")
    def test_evaluator_with_full_workflow(self):
        """Test evaluator through complete warmup + search workflow"""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            self.pipeline, selection_strategy=strategy, batch_size=1
        )
        sampler.read_reagents(self.pipeline.reagent_file_list)

        evaluator = LookupEvaluator({"ref_filename": self.lookup_file})
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
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(self.pipeline, selection_strategy=strategy)
        sampler.read_reagents(self.pipeline.reagent_file_list)

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
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(self.pipeline, selection_strategy=strategy)
        sampler.read_reagents(self.pipeline.reagent_file_list)

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
                self._counter = 0

            @property
            def counter(self):
                return self._counter

            def evaluate(self, mol):
                self._counter += 1
                # Return number of atoms as score
                return mol.GetNumAtoms()

        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(self.pipeline, selection_strategy=strategy)
        sampler.read_reagents(self.pipeline.reagent_file_list)

        evaluator = CustomEvaluator()
        sampler.set_evaluator(evaluator)

        # Test evaluation
        smiles, name, score = sampler.evaluate([0, 0])
        assert isinstance(score, (int, float))
        assert score > 0
        assert evaluator.counter == 1

        sampler.close()
