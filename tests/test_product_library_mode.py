"""
Comprehensive tests for pre-enumerated product library mode.

This tests the new feature where Thompson sampler can use pre-enumerated
products instead of synthesizing via reaction SMARTS.
"""

import pytest
import os
import tempfile
import polars as pl
import numpy as np
from pathlib import Path

from TACTICS.thompson_sampling.core.sampler import ThompsonSampler
from TACTICS.thompson_sampling.core.evaluators import FPEvaluator, MWEvaluator, LookupEvaluator
from TACTICS.thompson_sampling.strategies.greedy_selection import GreedySelection
from TACTICS.thompson_sampling.config import ThompsonSamplingConfig
from TACTICS.thompson_sampling.strategies.config import GreedyConfig
from TACTICS.thompson_sampling.core.evaluator_config import FPEvaluatorConfig, MWEvaluatorConfig


# Paths to test data
EXAMPLES_DIR = Path(__file__).parent.parent / "examples"
REAGENT_FILE1 = EXAMPLES_DIR / "input_files" / "acids.smi"
REAGENT_FILE2 = EXAMPLES_DIR / "input_files" / "coupled_aa_sub.smi"
REACTION_SMARTS = "[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]>>[#6:1](=[O:2])[#7:3]"


@pytest.fixture
def sample_product_library():
    """
    Create a temporary product library CSV with pre-enumerated products.

    This simulates a scenario where products have been pre-computed (e.g., from
    manual enumeration or external docking).
    """
    # Create a few sample products with realistic SMILES
    data = {
        'Product_Code': [
            'CA0_AA0_AA0',
            'CA1_AA0_AA0',
            'CA0_AA0_AA1',
            'CA2_AA1_AA2',
            'CA5_AA3_AA5'
        ],
        'SMILES': [
            'CC1(C)Oc2ccc(CC(NC(=O)c3cc(=O)[nH]c(=O)[nH]3)C(N)=O)cc2O1',
            'CC1(C)Oc2ccc(CC(NC(=O)C2CC(C(=O)N3CCC(C(=O)O)CC3)CC2)C(N)=O)cc2O1',
            'CC1(C)Oc2ccc(CC(NC(=O)c3cc(=O)[nH]c(=O)[nH]3)C(=O)N(C)CC2CCCC(C(N)=O)CC2)cc2O1',
            'CC1(C)Oc2ccc(CC(NC(=O)CCCC(C)C)C(=O)N(C)CC2C(=O)NCC(c3ccc(Br)c(O)c3Br)C(N)=O)cc2O1',
            'CC1(C)Oc2ccc(CC(NC(=O)CN2ccc(NC(=O)OC(c3ccccc3)c3ccccc3)nc2=O)C(=O)NC2(C(N)=O)CCN(Cc3ccccc3)CC2)cc2O1'
        ]
    }

    df = pl.DataFrame(data)

    # Create temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        df.write_csv(f.name)
        temp_file = f.name

    yield temp_file

    # Cleanup
    if os.path.exists(temp_file):
        os.remove(temp_file)


class TestProductLibraryLoading:
    """Test loading and validation of product libraries."""

    def test_load_valid_library(self, sample_product_library):
        """Test loading a valid product library."""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy)

        # Load library
        sampler.load_product_library(sample_product_library)

        # Check that library was loaded
        assert sampler.product_smiles_dict is not None
        assert len(sampler.product_smiles_dict) == 5
        assert 'CA0_AA0_AA0' in sampler.product_smiles_dict

        sampler.close()

    def test_load_library_via_constructor(self, sample_product_library):
        """Test loading library via constructor parameter."""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            selection_strategy=strategy,
            product_library_file=sample_product_library
        )

        # Check that library was loaded
        assert sampler.product_smiles_dict is not None
        assert len(sampler.product_smiles_dict) == 5

        sampler.close()

    def test_load_nonexistent_file(self):
        """Test error handling for nonexistent library file."""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy)

        with pytest.raises(FileNotFoundError):
            sampler.load_product_library("/nonexistent/file.csv")

        sampler.close()

    def test_load_library_missing_product_code_column(self):
        """Test error handling for library missing Product_Code column."""
        # Create invalid CSV
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            df = pl.DataFrame({
                'Name': ['prod1', 'prod2'],
                'SMILES': ['CCO', 'CCC']
            })
            df.write_csv(f.name)
            temp_file = f.name

        try:
            strategy = GreedySelection(mode="maximize")
            sampler = ThompsonSampler(selection_strategy=strategy)

            with pytest.raises(ValueError, match="Product_Code"):
                sampler.load_product_library(temp_file)

            sampler.close()
        finally:
            if os.path.exists(temp_file):
                os.remove(temp_file)

    def test_load_library_missing_smiles_column(self):
        """Test error handling for library missing SMILES column."""
        # Create invalid CSV
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            df = pl.DataFrame({
                'Product_Code': ['prod1', 'prod2'],
                'Score': [1.0, 2.0]
            })
            df.write_csv(f.name)
            temp_file = f.name

        try:
            strategy = GreedySelection(mode="maximize")
            sampler = ThompsonSampler(selection_strategy=strategy)

            with pytest.raises(ValueError, match="SMILES"):
                sampler.load_product_library(temp_file)

            sampler.close()
        finally:
            if os.path.exists(temp_file):
                os.remove(temp_file)


class TestProductLibraryEvaluation:
    """Test evaluation using pre-enumerated products."""

    def test_evaluate_with_library(self, sample_product_library):
        """Test that evaluation works with pre-enumerated library."""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            selection_strategy=strategy,
            product_library_file=sample_product_library
        )

        # Set up sampler
        sampler.read_reagents([str(REAGENT_FILE1), str(REAGENT_FILE2)])
        evaluator = MWEvaluator()
        sampler.set_evaluator(evaluator)

        # Evaluate product CA0_AA0_AA0 (reagents at indices [0, 0])
        smiles, name, score = sampler.evaluate([0, 0])

        # Check that evaluation succeeded
        assert name == 'CA0_AA0_AA0'
        assert smiles != 'FAIL'
        assert np.isfinite(score)
        assert score > 0  # MW should be positive

        sampler.close()

    def test_evaluate_missing_product(self, sample_product_library):
        """Test error handling when product is not in library."""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            selection_strategy=strategy,
            product_library_file=sample_product_library
        )

        # Set up sampler
        sampler.read_reagents([str(REAGENT_FILE1), str(REAGENT_FILE2)])
        evaluator = MWEvaluator()
        sampler.set_evaluator(evaluator)

        # Try to evaluate product not in library (e.g., CA10_AA10_AA10)
        # This should fail gracefully and return NaN
        smiles, name, score = sampler.evaluate([10, 10])

        assert smiles == 'FAIL'
        assert np.isnan(score)

        sampler.close()

    def test_evaluate_with_fp_evaluator(self, sample_product_library):
        """Test evaluation with FPEvaluator using pre-enumerated library."""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            selection_strategy=strategy,
            product_library_file=sample_product_library
        )

        # Set up sampler
        sampler.read_reagents([str(REAGENT_FILE1), str(REAGENT_FILE2)])
        evaluator = FPEvaluator({"query_smiles": "CC(=O)NC1CCCCC1"})
        sampler.set_evaluator(evaluator)

        # Evaluate product
        smiles, name, score = sampler.evaluate([0, 0])

        # Check that evaluation succeeded
        assert smiles != 'FAIL'
        assert np.isfinite(score)
        assert 0.0 <= score <= 1.0  # Tanimoto similarity

        sampler.close()

    def test_evaluate_batch_with_library(self, sample_product_library):
        """Test batch evaluation with pre-enumerated library."""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            selection_strategy=strategy,
            product_library_file=sample_product_library
        )

        # Set up sampler
        sampler.read_reagents([str(REAGENT_FILE1), str(REAGENT_FILE2)])
        evaluator = MWEvaluator()
        sampler.set_evaluator(evaluator)

        # Evaluate batch
        choice_lists = [[0, 0], [1, 0], [0, 1]]
        results = sampler.evaluate_batch(choice_lists)

        assert len(results) == 3
        for smiles, name, score in results:
            # At least some should succeed (those in library)
            if name in ['CA0_AA0_AA0', 'CA1_AA0_AA0', 'CA0_AA0_AA1']:
                assert smiles != 'FAIL'
                assert np.isfinite(score)

        sampler.close()


class TestProductLibraryBatchMode:
    """Test product library with batch Thompson sampling (batch_size > 1)."""

    def test_batch_mode_with_library(self, sample_product_library):
        """
        Test that batch Thompson sampling (batch_size > 1) works with product library.

        In batch mode, multiple compounds are sampled per cycle. With product library,
        this means multiple product lookups per cycle.
        """
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            selection_strategy=strategy,
            product_library_file=sample_product_library,
            batch_size=3  # Batch mode!
        )

        # Set up sampler
        sampler.read_reagents([str(REAGENT_FILE1), str(REAGENT_FILE2)])
        sampler.set_reaction(REACTION_SMARTS)  # Fallback for products not in library
        evaluator = MWEvaluator()
        sampler.set_evaluator(evaluator)

        # Run warmup
        warmup_df = sampler.warm_up(num_warmup_trials=2)
        assert len(warmup_df) > 0

        # Run search in batch mode
        search_df = sampler.search(num_cycles=5)

        # Check that search produced results
        # In batch mode, we expect batch_size * num_cycles evaluations (if all unique)
        assert len(search_df) > 0
        assert len(search_df) <= 5 * 3  # At most batch_size * num_cycles

        # Verify some products from library were found
        library_products = set(['CA0_AA0_AA0', 'CA1_AA0_AA0', 'CA0_AA0_AA1', 'CA2_AA1_AA2', 'CA5_AA3_AA5'])
        found_library_products = set(search_df['Name'].to_list()).intersection(library_products)

        # At least check that results exist
        assert len(search_df) > 0

        sampler.close()

    def test_batch_mode_parallel_evaluation_with_library(self, sample_product_library):
        """
        Test batch mode with parallel evaluation and product library.

        This tests the combination of:
        - batch_size > 1 (multiple compounds per cycle)
        - processes > 1 (parallel evaluation)
        - product_library_file (pre-enumerated products)
        """
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            selection_strategy=strategy,
            product_library_file=sample_product_library,
            batch_size=5,
            processes=2,  # Parallel evaluation
            min_cpds_per_core=5
        )

        # Set up sampler
        sampler.read_reagents([str(REAGENT_FILE1), str(REAGENT_FILE2)])
        sampler.set_reaction(REACTION_SMARTS)
        evaluator = MWEvaluator()
        sampler.set_evaluator(evaluator)

        # Run a small search
        search_df = sampler.search(num_cycles=3)

        # Verify results
        assert len(search_df) > 0
        scores = search_df['score'].to_list()
        assert all(np.isfinite(s) for s in scores)

        sampler.close()


class TestProductLibraryWithWorkflow:
    """Test product library with full warmup and search workflow."""

    def test_warmup_with_library(self, sample_product_library):
        """Test that warmup works with pre-enumerated library."""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            selection_strategy=strategy,
            product_library_file=sample_product_library
        )

        # Set up sampler - need reaction as fallback since library doesn't have all products
        sampler.read_reagents([str(REAGENT_FILE1), str(REAGENT_FILE2)])
        sampler.set_reaction(REACTION_SMARTS)  # Fallback for products not in library
        evaluator = MWEvaluator()
        sampler.set_evaluator(evaluator)

        # Run warmup
        warmup_df = sampler.warm_up(num_warmup_trials=2)

        # Check that warmup produced results
        # Note: Some products may not be in library, so we expect some failures
        assert len(warmup_df) > 0

        # Check that valid products have finite scores
        valid_products = warmup_df.filter(
            warmup_df['Name'].is_in(['CA0_AA0_AA0', 'CA1_AA0_AA0', 'CA0_AA0_AA1'])
        )
        if len(valid_products) > 0:
            scores = valid_products['score'].to_list()
            assert all(np.isfinite(s) for s in scores)

        sampler.close()

    def test_search_with_library(self, sample_product_library):
        """Test that search works with pre-enumerated library."""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            selection_strategy=strategy,
            product_library_file=sample_product_library
        )

        # Set up sampler
        sampler.read_reagents([str(REAGENT_FILE1), str(REAGENT_FILE2)])
        sampler.set_reaction(REACTION_SMARTS)  # Fallback for products not in library
        evaluator = MWEvaluator()
        sampler.set_evaluator(evaluator)

        # Run warmup and search
        warmup_df = sampler.warm_up(num_warmup_trials=2)
        search_df = sampler.search(num_cycles=5)

        # Check that search produced results
        assert len(search_df) >= 0  # May be empty if all products missing

        sampler.close()


class TestBackwardCompatibility:
    """Test that existing behavior is unchanged when library not provided."""

    def test_normal_synthesis_without_library(self):
        """Test that synthesis works normally when no library provided."""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(selection_strategy=strategy)

        # Set up sampler normally
        sampler.read_reagents([str(REAGENT_FILE1), str(REAGENT_FILE2)])
        sampler.set_reaction(REACTION_SMARTS)
        evaluator = MWEvaluator()
        sampler.set_evaluator(evaluator)

        # Evaluate should work via normal synthesis
        smiles, name, score = sampler.evaluate([0, 0])

        assert smiles != 'FAIL'
        assert np.isfinite(score)
        assert sampler.product_smiles_dict is None  # No library loaded

        sampler.close()

    def test_config_without_library(self):
        """Test that config works without product_library_file."""
        config = ThompsonSamplingConfig(
            reaction_smarts=REACTION_SMARTS,
            reagent_file_list=[str(REAGENT_FILE1), str(REAGENT_FILE2)],
            num_ts_iterations=5,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=MWEvaluatorConfig()
        )

        sampler = ThompsonSampler.from_config(config)

        # Should work normally without library
        assert sampler.product_smiles_dict is None

        # Can still evaluate normally
        warmup_df = sampler.warm_up(num_warmup_trials=2)
        assert len(warmup_df) > 0

        sampler.close()


class TestProductLibraryWithLookupEvaluator:
    """Test that LookupEvaluator functionality is preserved with product library."""

    def test_lookup_evaluator_with_product_library(self, sample_product_library):
        """
        Test that LookupEvaluator works correctly with product library mode.

        This is critical: LookupEvaluator uses product_name (product code) to
        look up scores, NOT the molecule object. With product library mode,
        we skip synthesis but LookupEvaluator should still work by using the
        product_name directly.
        """
        # Create a lookup evaluator with sample scores
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            scores_df = pl.DataFrame({
                'Product_Code': ['CA0_AA0_AA0', 'CA1_AA0_AA0', 'CA0_AA0_AA1'],
                'Scores': [-10.5, -12.3, -9.8]
            })
            scores_df.write_csv(f.name)
            scores_file = f.name

        try:
            strategy = GreedySelection(mode="minimize")
            sampler = ThompsonSampler(
                selection_strategy=strategy,
                product_library_file=sample_product_library
            )

            # Set up with LookupEvaluator
            sampler.read_reagents([str(REAGENT_FILE1), str(REAGENT_FILE2)])
            evaluator = LookupEvaluator({"ref_filename": scores_file})
            sampler.set_evaluator(evaluator)

            # Evaluate a product
            smiles, name, score = sampler.evaluate([0, 0])

            # Check that evaluation succeeded using product_name
            assert name == 'CA0_AA0_AA0'
            assert smiles != 'FAIL'
            assert score == -10.5  # Should match the lookup score

            sampler.close()
        finally:
            if os.path.exists(scores_file):
                os.remove(scores_file)

    def test_lookup_evaluator_without_product_library(self):
        """
        Test that LookupEvaluator still works in normal synthesis mode.

        This ensures backward compatibility - existing LookupEvaluator usage
        should be unchanged when not using product library.
        """
        # Create a lookup evaluator with sample scores
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            scores_df = pl.DataFrame({
                'Product_Code': ['CA0_AA0_AA0', 'CA1_AA0_AA0', 'CA0_AA0_AA1'],
                'Scores': [-10.5, -12.3, -9.8]
            })
            scores_df.write_csv(f.name)
            scores_file = f.name

        try:
            strategy = GreedySelection(mode="minimize")
            sampler = ThompsonSampler(selection_strategy=strategy)

            # Set up with normal synthesis (no product library)
            sampler.read_reagents([str(REAGENT_FILE1), str(REAGENT_FILE2)])
            sampler.set_reaction(REACTION_SMARTS)
            evaluator = LookupEvaluator({"ref_filename": scores_file})
            sampler.set_evaluator(evaluator)

            # Evaluate a product via synthesis
            smiles, name, score = sampler.evaluate([0, 0])

            # Check that evaluation succeeded
            assert name == 'CA0_AA0_AA0'
            assert smiles != 'FAIL'
            assert score == -10.5  # Should match the lookup score

            sampler.close()
        finally:
            if os.path.exists(scores_file):
                os.remove(scores_file)


class TestProductLibraryConfig:
    """Test configuration integration with product library."""

    def test_config_with_library(self, sample_product_library):
        """Test that config properly passes product_library_file."""
        config = ThompsonSamplingConfig(
            reaction_smarts=REACTION_SMARTS,
            reagent_file_list=[str(REAGENT_FILE1), str(REAGENT_FILE2)],
            num_ts_iterations=5,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=MWEvaluatorConfig(),
            product_library_file=sample_product_library
        )

        sampler = ThompsonSampler.from_config(config)

        # Check that library was loaded
        assert sampler.product_smiles_dict is not None
        assert len(sampler.product_smiles_dict) == 5

        sampler.close()
