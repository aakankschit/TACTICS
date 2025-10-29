"""
Integration tests for evaluator classes with ThompsonSampler.

These tests ensure that any evaluator (existing or new) can be properly
integrated with the ThompsonSampler class.
"""

import pytest
import tempfile
import os
import numpy as np
import pandas as pd
from rdkit import Chem
from sqlitedict import SqliteDict

from TACTICS.thompson_sampling.core.sampler import ThompsonSampler
from TACTICS.thompson_sampling.core.evaluators import (
    Evaluator, MWEvaluator, FPEvaluator, LookupEvaluator, DBEvaluator
)
from TACTICS.thompson_sampling.strategies.greedy_selection import GreedySelection


class TestEvaluatorIntegration:
    """Test that evaluators integrate properly with ThompsonSampler"""

    def setup_method(self):
        """Set up test fixtures"""
        self.temp_dir = tempfile.mkdtemp()

        # Create test reagent files
        self.reagent_file1 = os.path.join(self.temp_dir, "reagents1.smi")
        self.reagent_file2 = os.path.join(self.temp_dir, "reagents2.smi")

        with open(self.reagent_file1, "w") as f:
            f.write("CCO\tethanol\n")
            f.write("CCCO\tpropanol\n")

        with open(self.reagent_file2, "w") as f:
            f.write("CC(=O)O\tacetic_acid\n")
            f.write("CCC(=O)O\tpropionic_acid\n")

        # Create test lookup file
        self.lookup_file = os.path.join(self.temp_dir, "scores.csv")
        lookup_data = pd.DataFrame({
            'Product_Code': [
                'ethanol_acetic_acid',
                'ethanol_propionic_acid',
                'propanol_acetic_acid',
                'propanol_propionic_acid'
            ],
            'Scores': [0.5, 0.7, 0.6, 0.8]
        })
        lookup_data.to_csv(self.lookup_file, index=False)

        # Create test database
        self.db_file = os.path.join(self.temp_dir, "scores.db")
        db = SqliteDict(self.db_file)
        db['test_ethanol_acetic_acid'] = 0.5
        db['test_ethanol_propionic_acid'] = 0.7
        db['test_propanol_acetic_acid'] = 0.6
        db['test_propanol_propionic_acid'] = 0.8
        db.commit()
        db.close()

        # Simple reaction SMARTS (alcohol + acid -> ester)
        self.reaction_smarts = "[C:1][OH:2].[C:3](=[O:4])[OH:5]>>[C:1][O:2][C:3](=[O:4])"

    def teardown_method(self):
        """Clean up test fixtures"""
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_mw_evaluator_integration(self):
        """Test that MWEvaluator integrates with ThompsonSampler"""
        # Create sampler with MWEvaluator
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)

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

        # Test multiple evaluations
        for i in range(3):
            sampler.evaluate([0, 0])
        assert sampler.evaluator.counter == 4

    def test_fp_evaluator_integration(self):
        """Test that FPEvaluator integrates with ThompsonSampler"""
        # Create sampler with FPEvaluator
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)

        # Use a simple target molecule
        evaluator = FPEvaluator({"query_smiles": "CCOC(=O)C"})
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

    def test_lookup_evaluator_integration(self):
        """Test that LookupEvaluator integrates with ThompsonSampler"""
        # Create sampler with LookupEvaluator
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)

        evaluator = LookupEvaluator({"ref_filename": self.lookup_file})
        sampler.set_evaluator(evaluator)

        # Verify evaluator is set correctly
        assert sampler.evaluator is not None
        assert isinstance(sampler.evaluator, LookupEvaluator)
        assert sampler.evaluator.counter == 0

        # Test evaluation
        smiles, name, score = sampler.evaluate([0, 0])
        assert isinstance(score, (int, float))
        assert score == 0.5  # Known score from lookup file
        assert sampler.evaluator.counter == 1

    def test_db_evaluator_integration(self):
        """Test that DBEvaluator integrates with ThompsonSampler"""
        # Create sampler with DBEvaluator
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)

        evaluator = DBEvaluator({
            "db_filename": self.db_file,
            "db_prefix": "test_"
        })
        sampler.set_evaluator(evaluator)

        # Verify evaluator is set correctly
        assert sampler.evaluator is not None
        assert isinstance(sampler.evaluator, DBEvaluator)
        assert sampler.evaluator.counter == 0

        # Test evaluation
        smiles, name, score = sampler.evaluate([0, 0])
        assert isinstance(score, (int, float))
        assert score == 0.5  # Known score from database
        assert sampler.evaluator.counter == 1

    def test_evaluator_counter_increments(self):
        """Test that evaluator counter increments correctly"""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy, batch_size=5)
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)

        evaluator = MWEvaluator()
        sampler.set_evaluator(evaluator)

        initial_count = evaluator.counter

        # Evaluate batch of compounds
        batch_results = sampler.evaluate_batch([[0, 0], [0, 1], [1, 0]])

        # Counter should increment by number of evaluations
        assert evaluator.counter == initial_count + 3
        assert len(batch_results) == 3

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
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)

        evaluator = CustomEvaluator()
        sampler.set_evaluator(evaluator)

        # Test evaluation
        smiles, name, score = sampler.evaluate([0, 0])
        assert isinstance(score, (int, float))
        assert score > 0
        assert evaluator.counter == 1

    def test_evaluator_with_warmup(self):
        """Test that evaluators work correctly during warmup phase"""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy)
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)

        evaluator = LookupEvaluator({"ref_filename": self.lookup_file})
        sampler.set_evaluator(evaluator)

        # Run warmup
        warmup_results = sampler.warm_up(num_warmup_trials=2)

        # Verify warmup ran and evaluator was called
        assert len(warmup_results) > 0
        assert evaluator.counter > 0

        # All warmup results should have valid scores
        for result in warmup_results:
            score, smiles, name = result
            assert isinstance(score, (int, float))
            assert not np.isnan(score)

    def test_evaluator_with_full_workflow(self):
        """Test evaluator through complete warmup + search workflow"""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy, batch_size=1)
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)

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

        # Cleanup
        sampler.close()

    def test_parallel_evaluator_multiprocessing_warning(self):
        """Test that fast evaluators warn about multiprocessing inefficiency"""
        import logging

        # Capture log output
        logger = logging.getLogger('TACTICS.thompson_sampling.core.sampler')
        logger.setLevel(logging.WARNING)

        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            selection_strategy=strategy,
            processes=4  # Enable multiprocessing
        )
        sampler.read_reagents([self.reagent_file1, self.reagent_file2])
        sampler.set_reaction(self.reaction_smarts)

        # LookupEvaluator and DBEvaluator should trigger warning
        evaluator = LookupEvaluator({"ref_filename": self.lookup_file})

        # This should log a warning about multiprocessing inefficiency
        sampler.set_evaluator(evaluator)

        # The sampler should still work despite the warning
        assert sampler.evaluator is not None
        sampler.close()
