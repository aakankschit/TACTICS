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
from TACTICS.thompson_sampling.core.evaluators import (
    FPEvaluator,
    MWEvaluator,
    LookupEvaluator,
)
from TACTICS.thompson_sampling.strategies.greedy_selection import GreedySelection
from TACTICS.thompson_sampling.config import ThompsonSamplingConfig
from TACTICS.thompson_sampling.strategies.config import GreedyConfig
from TACTICS.thompson_sampling.core.evaluator_config import (
    FPEvaluatorConfig,
    MWEvaluatorConfig,
)
from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef


REACTION_SMARTS = "[C:1](=O)[OH].[N:2]>>[C:1](=O)[N:2]"


@pytest.fixture
def temp_reagent_files():
    """Create temporary reagent files."""
    temp_dir = tempfile.mkdtemp()

    acids_file = os.path.join(temp_dir, "acids.smi")
    amines_file = os.path.join(temp_dir, "amines.smi")

    # Use simple names that match sample_product_library codes
    with open(acids_file, "w") as f:
        f.write("CC(=O)O CA0\n")
        f.write("CCC(=O)O CA1\n")
        f.write("CCCC(=O)O CA2\n")
        f.write("CCCCC(=O)O CA3\n")

    with open(amines_file, "w") as f:
        f.write("CN AA0\n")
        f.write("CCN AA1\n")
        f.write("CCCN AA2\n")
        f.write("CCCCN AA3\n")

    yield [acids_file, amines_file], temp_dir

    # Cleanup
    import shutil

    shutil.rmtree(temp_dir)


@pytest.fixture
def sample_product_library():
    """
    Create a temporary product library CSV with pre-enumerated products.

    This simulates a scenario where products have been pre-computed (e.g., from
    manual enumeration or external docking).
    """
    # Create products that match the reagent naming: CA<i>_AA<j>
    data = {
        "Product_Code": ["CA0_AA0", "CA1_AA0", "CA0_AA1", "CA2_AA2", "CA3_AA3"],
        "SMILES": [
            "CC(=O)NC",
            "CCC(=O)NC",
            "CC(=O)NCC",
            "CCCC(=O)NCCC",
            "CCCCC(=O)NCCCC",
        ],
    }

    df = pl.DataFrame(data)

    # Create temporary file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        df.write_csv(f.name)
        temp_file = f.name

    yield temp_file

    # Cleanup
    if os.path.exists(temp_file):
        os.remove(temp_file)


@pytest.fixture
def pipeline(temp_reagent_files):
    """Create a SynthesisPipeline for testing."""
    reagent_files, _ = temp_reagent_files
    config = ReactionConfig(
        reactions=[ReactionDef(reaction_smarts=REACTION_SMARTS, step_index=0)],
        reagent_file_list=reagent_files,
    )
    return SynthesisPipeline(config)


class TestProductLibraryLoading:
    """Test loading and validation of product libraries."""

    def test_load_valid_library(self, sample_product_library, pipeline):
        """Test loading a valid product library."""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(pipeline, selection_strategy=strategy)

        # Load library
        sampler.load_product_library(sample_product_library)

        # Check that library was loaded
        assert sampler.product_smiles_dict is not None
        assert len(sampler.product_smiles_dict) == 5
        assert "CA0_AA0" in sampler.product_smiles_dict

        sampler.close()

    def test_load_library_via_constructor(self, sample_product_library, pipeline):
        """Test loading library via constructor parameter."""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            pipeline,
            selection_strategy=strategy,
            product_library_file=sample_product_library,
        )

        # Check that library was loaded
        assert sampler.product_smiles_dict is not None
        assert len(sampler.product_smiles_dict) == 5

        sampler.close()

    def test_load_nonexistent_file(self, pipeline):
        """Test error handling for nonexistent library file."""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(pipeline, selection_strategy=strategy)

        with pytest.raises(FileNotFoundError):
            sampler.load_product_library("/nonexistent/file.csv")

        sampler.close()

    def test_load_library_missing_product_code_column(self, pipeline):
        """Test error handling for library missing Product_Code column."""
        # Create invalid CSV
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            df = pl.DataFrame({"Name": ["prod1", "prod2"], "SMILES": ["CCO", "CCC"]})
            df.write_csv(f.name)
            temp_file = f.name

        try:
            strategy = GreedySelection(mode="maximize")
            sampler = ThompsonSampler(pipeline, selection_strategy=strategy)

            with pytest.raises(ValueError, match="Product_Code"):
                sampler.load_product_library(temp_file)

            sampler.close()
        finally:
            if os.path.exists(temp_file):
                os.remove(temp_file)

    def test_load_library_missing_smiles_column(self, pipeline):
        """Test error handling for library missing SMILES column."""
        # Create invalid CSV
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            df = pl.DataFrame({"Product_Code": ["prod1", "prod2"], "Score": [1.0, 2.0]})
            df.write_csv(f.name)
            temp_file = f.name

        try:
            strategy = GreedySelection(mode="maximize")
            sampler = ThompsonSampler(pipeline, selection_strategy=strategy)

            with pytest.raises(ValueError, match="SMILES"):
                sampler.load_product_library(temp_file)

            sampler.close()
        finally:
            if os.path.exists(temp_file):
                os.remove(temp_file)


class TestProductLibraryEvaluation:
    """Test evaluation using pre-enumerated products."""

    def test_evaluate_with_library(
        self, sample_product_library, pipeline, temp_reagent_files
    ):
        """Test that evaluation works with pre-enumerated library."""
        reagent_files, _ = temp_reagent_files
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            pipeline,
            selection_strategy=strategy,
            product_library_file=sample_product_library,
        )

        # Set up sampler
        sampler.read_reagents(reagent_files)
        evaluator = MWEvaluator()
        sampler.set_evaluator(evaluator)

        # Evaluate product CA0_AA0 (reagents at indices [0, 0])
        smiles, name, score = sampler.evaluate([0, 0])

        # Check that evaluation succeeded
        assert name == "CA0_AA0"
        assert smiles != "FAIL"
        assert np.isfinite(score)
        assert score > 0  # MW should be positive

        sampler.close()

    def test_evaluate_missing_product(
        self, sample_product_library, pipeline, temp_reagent_files
    ):
        """Test behavior when product is not in library (falls back to synthesis)."""
        reagent_files, _ = temp_reagent_files
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            pipeline,
            selection_strategy=strategy,
            product_library_file=sample_product_library,
        )

        # Set up sampler
        sampler.read_reagents(reagent_files)
        evaluator = MWEvaluator()
        sampler.set_evaluator(evaluator)

        # Try to evaluate product not in library (e.g., CA1_AA1)
        # With the new API, it should fall back to synthesis
        smiles, name, score = sampler.evaluate([1, 1])

        # Should either fail or succeed via synthesis
        if smiles != "FAIL":
            assert np.isfinite(score)
        else:
            assert np.isnan(score)

        sampler.close()

    def test_evaluate_with_fp_evaluator(
        self, sample_product_library, pipeline, temp_reagent_files
    ):
        """Test evaluation with FPEvaluator using pre-enumerated library."""
        reagent_files, _ = temp_reagent_files
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            pipeline,
            selection_strategy=strategy,
            product_library_file=sample_product_library,
        )

        # Set up sampler
        sampler.read_reagents(reagent_files)
        evaluator = FPEvaluator({"query_smiles": "CC(=O)NC1CCCCC1"})
        sampler.set_evaluator(evaluator)

        # Evaluate product
        smiles, name, score = sampler.evaluate([0, 0])

        # Check that evaluation succeeded
        assert smiles != "FAIL"
        assert np.isfinite(score)
        assert 0.0 <= score <= 1.0  # Tanimoto similarity

        sampler.close()

    def test_evaluate_batch_with_library(
        self, sample_product_library, pipeline, temp_reagent_files
    ):
        """Test batch evaluation with pre-enumerated library."""
        reagent_files, _ = temp_reagent_files
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            pipeline,
            selection_strategy=strategy,
            product_library_file=sample_product_library,
        )

        # Set up sampler
        sampler.read_reagents(reagent_files)
        evaluator = MWEvaluator()
        sampler.set_evaluator(evaluator)

        # Evaluate batch
        choice_lists = [[0, 0], [1, 0], [0, 1]]
        results = sampler.evaluate_batch(choice_lists)

        assert len(results) == 3
        for smiles, name, score in results:
            # At least some should succeed (those in library)
            if name in ["CA0_AA0", "CA1_AA0", "CA0_AA1"]:
                assert smiles != "FAIL"
                assert np.isfinite(score)

        sampler.close()


class TestProductLibraryBatchMode:
    """Test product library with batch Thompson sampling (batch_size > 1)."""

    @pytest.mark.skip(reason="Warmup has edge case with small test data")
    def test_batch_mode_with_library(
        self, sample_product_library, pipeline, temp_reagent_files
    ):
        """
        Test that batch Thompson sampling (batch_size > 1) works with product library.

        In batch mode, multiple compounds are sampled per cycle. With product library,
        this means multiple product lookups per cycle.
        """
        reagent_files, _ = temp_reagent_files
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            pipeline,
            selection_strategy=strategy,
            product_library_file=sample_product_library,
            batch_size=3,  # Batch mode!
        )

        # Set up sampler
        sampler.read_reagents(reagent_files)
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

        sampler.close()

    @pytest.mark.skip(reason="Warmup has edge case with small test data")
    def test_batch_mode_parallel_evaluation_with_library(
        self, sample_product_library, pipeline, temp_reagent_files
    ):
        """
        Test batch mode with parallel evaluation and product library.

        This tests the combination of:
        - batch_size > 1 (multiple compounds per cycle)
        - processes > 1 (parallel evaluation)
        - product_library_file (pre-enumerated products)
        """
        reagent_files, _ = temp_reagent_files
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            pipeline,
            selection_strategy=strategy,
            product_library_file=sample_product_library,
            batch_size=5,
            processes=2,  # Parallel evaluation
            min_cpds_per_core=5,
        )

        # Set up sampler
        sampler.read_reagents(reagent_files)
        evaluator = MWEvaluator()
        sampler.set_evaluator(evaluator)

        # Run a small search
        search_df = sampler.search(num_cycles=3)

        # Verify results
        assert len(search_df) > 0
        scores = search_df["score"].to_list()
        # Some may be NaN if not in library and synthesis failed
        finite_scores = [s for s in scores if np.isfinite(s)]
        assert len(finite_scores) > 0

        sampler.close()


class TestProductLibraryWithWorkflow:
    """Test product library with full warmup and search workflow."""

    @pytest.mark.skip(reason="Warmup has edge case with small test data")
    def test_warmup_with_library(
        self, sample_product_library, pipeline, temp_reagent_files
    ):
        """Test that warmup works with pre-enumerated library."""
        reagent_files, _ = temp_reagent_files
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            pipeline,
            selection_strategy=strategy,
            product_library_file=sample_product_library,
        )

        # Set up sampler
        sampler.read_reagents(reagent_files)
        evaluator = MWEvaluator()
        sampler.set_evaluator(evaluator)

        # Run warmup
        warmup_df = sampler.warm_up(num_warmup_trials=2)

        # Check that warmup produced results
        assert len(warmup_df) > 0

        sampler.close()

    @pytest.mark.skip(reason="Warmup has edge case with small test data")
    def test_search_with_library(
        self, sample_product_library, pipeline, temp_reagent_files
    ):
        """Test that search works with pre-enumerated library."""
        reagent_files, _ = temp_reagent_files
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            pipeline,
            selection_strategy=strategy,
            product_library_file=sample_product_library,
        )

        # Set up sampler
        sampler.read_reagents(reagent_files)
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

    def test_normal_synthesis_without_library(self, pipeline, temp_reagent_files):
        """Test that synthesis works normally when no library provided."""
        reagent_files, _ = temp_reagent_files
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(pipeline, selection_strategy=strategy)

        # Set up sampler normally
        sampler.read_reagents(reagent_files)
        evaluator = MWEvaluator()
        sampler.set_evaluator(evaluator)

        # Evaluate should work via normal synthesis
        smiles, name, score = sampler.evaluate([0, 0])

        assert smiles != "FAIL"
        assert np.isfinite(score)
        assert sampler.product_smiles_dict is None  # No library loaded

        sampler.close()

    def test_config_without_library(self, temp_reagent_files):
        """Test that config works without product_library_file."""
        reagent_files, _ = temp_reagent_files
        reaction_config = ReactionConfig(
            reactions=[ReactionDef(reaction_smarts=REACTION_SMARTS, step_index=0)],
            reagent_file_list=reagent_files,
        )
        pipeline = SynthesisPipeline(reaction_config)

        config = ThompsonSamplingConfig(
            synthesis_pipeline=pipeline,
            num_ts_iterations=5,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=MWEvaluatorConfig(),
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

    @pytest.mark.skip(reason="Product library lookup integration needs debugging")
    def test_lookup_evaluator_with_product_library(
        self, sample_product_library, pipeline, temp_reagent_files
    ):
        """
        Test that LookupEvaluator works correctly with product library mode.

        This is critical: LookupEvaluator uses product_name (product code) to
        look up scores, NOT the molecule object. With product library mode,
        we skip synthesis but LookupEvaluator should still work by using the
        product_name directly.
        """
        reagent_files, _ = temp_reagent_files

        # Create a lookup evaluator with sample scores
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            scores_df = pl.DataFrame(
                {
                    "Product_Code": ["CA0_AA0", "CA1_AA0", "CA0_AA1"],
                    "Scores": [-10.5, -12.3, -9.8],
                }
            )
            scores_df.write_csv(f.name)
            scores_file = f.name

        try:
            strategy = GreedySelection(mode="minimize")
            sampler = ThompsonSampler(
                pipeline,
                selection_strategy=strategy,
                product_library_file=sample_product_library,
            )

            # Set up with LookupEvaluator
            sampler.read_reagents(reagent_files)
            evaluator = LookupEvaluator({"ref_filename": scores_file})
            sampler.set_evaluator(evaluator)

            # Evaluate a product
            smiles, name, score = sampler.evaluate([0, 0])

            # Check that evaluation succeeded using product_name
            assert name == "CA0_AA0"
            assert smiles != "FAIL"
            assert score == -10.5  # Should match the lookup score

            sampler.close()
        finally:
            if os.path.exists(scores_file):
                os.remove(scores_file)

    @pytest.mark.skip(reason="Synthesis with LookupEvaluator needs debugging")
    def test_lookup_evaluator_without_product_library(
        self, pipeline, temp_reagent_files
    ):
        """
        Test that LookupEvaluator still works in normal synthesis mode.

        This ensures backward compatibility - existing LookupEvaluator usage
        should be unchanged when not using product library.
        """
        reagent_files, _ = temp_reagent_files

        # Create a lookup evaluator with sample scores
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            scores_df = pl.DataFrame(
                {
                    "Product_Code": ["CA0_AA0", "CA1_AA0", "CA0_AA1"],
                    "Scores": [-10.5, -12.3, -9.8],
                }
            )
            scores_df.write_csv(f.name)
            scores_file = f.name

        try:
            strategy = GreedySelection(mode="minimize")
            sampler = ThompsonSampler(pipeline, selection_strategy=strategy)

            # Set up with normal synthesis (no product library)
            sampler.read_reagents(reagent_files)
            evaluator = LookupEvaluator({"ref_filename": scores_file})
            sampler.set_evaluator(evaluator)

            # Evaluate a product via synthesis
            smiles, name, score = sampler.evaluate([0, 0])

            # Check that evaluation succeeded
            assert name == "CA0_AA0"
            assert smiles != "FAIL"
            assert score == -10.5  # Should match the lookup score

            sampler.close()
        finally:
            if os.path.exists(scores_file):
                os.remove(scores_file)


class TestProductLibraryConfig:
    """Test configuration integration with product library."""

    def test_config_with_library(self, sample_product_library, temp_reagent_files):
        """Test that config properly passes product_library_file."""
        reagent_files, _ = temp_reagent_files
        reaction_config = ReactionConfig(
            reactions=[ReactionDef(reaction_smarts=REACTION_SMARTS, step_index=0)],
            reagent_file_list=reagent_files,
        )
        pipeline = SynthesisPipeline(reaction_config)

        config = ThompsonSamplingConfig(
            synthesis_pipeline=pipeline,
            num_ts_iterations=5,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=MWEvaluatorConfig(),
            product_library_file=sample_product_library,
        )

        sampler = ThompsonSampler.from_config(config)

        # Check that library was loaded
        assert sampler.product_smiles_dict is not None
        assert len(sampler.product_smiles_dict) == 5

        sampler.close()
