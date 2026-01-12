"""
Tests for Pydantic config model validation.

These tests verify that the configuration models work correctly
with the new synthesis_pipeline-based API.
"""

import pytest
import tempfile
import os
from pydantic import ValidationError

from TACTICS.thompson_sampling import ThompsonSamplingConfig, RandomBaselineConfig
from TACTICS.thompson_sampling.strategies.config import (
    GreedyConfig,
    RouletteWheelConfig,
    UCBConfig,
    EpsilonGreedyConfig,
    BoltzmannConfig,
    BayesUCBConfig,
)
from TACTICS.thompson_sampling.core.evaluator_config import (
    LookupEvaluatorConfig,
    DBEvaluatorConfig,
)
from TACTICS.library_enumeration import SynthesisPipeline
from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef


class TestConfigValidation:
    """
    Tests for Pydantic config model validation.
    These tests verify that the configuration models work correctly without external dependencies.
    """

    def setup_method(self):
        """
        Set up test fixtures before each test method.
        Creates temporary files for testing.
        """
        self.temp_dir = tempfile.mkdtemp()
        self.reagent_file1 = os.path.join(self.temp_dir, "reagents1.smi")
        self.reagent_file2 = os.path.join(self.temp_dir, "reagents2.smi")

        # Write sample reagent data
        with open(self.reagent_file1, "w") as f:
            f.write("CCO\tethanol\nCCCO\tpropanol\n")
        with open(self.reagent_file2, "w") as f:
            f.write("CC(=O)O\tacetic_acid\nCC(C)(C)O\ttert_butanol\n")

        # Create a simple reaction config and pipeline for testing
        self.reaction_config = ReactionConfig(
            reactions=[
                ReactionDef(
                    reaction_smarts="[C:1](=O)[OH].[O:2]>>[C:1](=O)[O:2]",
                    step_index=0,
                )
            ],
            reagent_file_list=[self.reagent_file1, self.reagent_file2],
        )
        self.pipeline = SynthesisPipeline(self.reaction_config)

    def teardown_method(self):
        """
        Clean up test fixtures after each test method.
        """
        import shutil

        shutil.rmtree(self.temp_dir)

    def test_thompson_sampling_config_creation(self):
        """
        Test that ThompsonSamplingConfig can be created with valid data.

        Inputs:
            - Valid configuration parameters for Thompson sampling

        Outputs:
            - ThompsonSamplingConfig instance is created successfully
        """
        config = ThompsonSamplingConfig(
            synthesis_pipeline=self.pipeline,
            num_ts_iterations=10,
            num_warmup_trials=3,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=DBEvaluatorConfig(db_filename="test.db"),
        )

        assert config.num_ts_iterations == 10
        assert config.num_warmup_trials == 3
        assert config.synthesis_pipeline is self.pipeline
        assert config.reagent_file_list == [self.reagent_file1, self.reagent_file2]

    def test_config_with_roulette_wheel_strategy(self):
        """
        Test that config works with roulette wheel selection strategy.

        Inputs:
            - Valid configuration with roulette_wheel strategy and parameters

        Outputs:
            - ThompsonSamplingConfig instance with strategy params
        """
        config = ThompsonSamplingConfig(
            synthesis_pipeline=self.pipeline,
            num_ts_iterations=10,
            num_warmup_trials=3,
            strategy_config=RouletteWheelConfig(
                mode="maximize",
                alpha=0.1,
                beta=0.1,
            ),
            evaluator_config=DBEvaluatorConfig(db_filename="test.db"),
        )

        assert isinstance(config.strategy_config, RouletteWheelConfig)
        assert config.strategy_config.alpha == 0.1
        assert config.strategy_config.beta == 0.1

    def test_config_with_optional_fields(self):
        """
        Test that optional fields work correctly.

        Inputs:
            - Config with optional fields specified

        Outputs:
            - Config is created with optional fields set correctly
        """
        config = ThompsonSamplingConfig(
            synthesis_pipeline=self.pipeline,
            num_ts_iterations=10,
            num_warmup_trials=3,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=DBEvaluatorConfig(db_filename="test.db"),
            results_filename="results.csv",
            log_filename="test.log",
            batch_size=5,
            max_resamples=100,
        )

        assert config.results_filename == "results.csv"
        assert config.log_filename == "test.log"
        assert config.batch_size == 5
        assert config.max_resamples == 100

    def test_random_baseline_config_creation(self):
        """
        Test that RandomBaselineConfig can be created with valid data.

        Inputs:
            - Valid configuration for random baseline

        Outputs:
            - RandomBaselineConfig instance is created successfully
        """
        config = RandomBaselineConfig(
            synthesis_pipeline=self.pipeline,
            evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv"),
            num_trials=100,
            num_to_save=10,
            ascending_output=False,
            outfile_name="baseline_results.csv",
        )

        assert config.num_trials == 100
        assert config.num_to_save == 10
        assert config.ascending_output == False
        assert config.outfile_name == "baseline_results.csv"

    def test_validation_errors_missing_pipeline(self):
        """
        Test that validation errors are raised when synthesis_pipeline is missing.

        Inputs:
            - Config without required synthesis_pipeline

        Outputs:
            - ValidationError is raised
        """
        with pytest.raises(ValidationError):
            ThompsonSamplingConfig(
                num_ts_iterations=10,
                num_warmup_trials=3,
                strategy_config=GreedyConfig(mode="maximize"),
                evaluator_config=DBEvaluatorConfig(db_filename="test.db"),
            )

    def test_validation_errors_invalid_mode(self):
        """
        Test that validation errors are raised for invalid mode.

        Inputs:
            - Invalid mode parameter

        Outputs:
            - ValidationError is raised
        """
        with pytest.raises(ValidationError):
            GreedyConfig(mode="invalid_mode")

    def test_type_validation(self):
        """
        Test that type validation works correctly.

        Inputs:
            - Configuration with wrong types

        Outputs:
            - ValidationError is raised for type mismatches
        """
        # Test string instead of int for num_ts_iterations
        with pytest.raises(ValidationError):
            ThompsonSamplingConfig(
                synthesis_pipeline=self.pipeline,
                num_ts_iterations="not_an_integer",  # Invalid type
                num_warmup_trials=3,
                strategy_config=GreedyConfig(mode="maximize"),
                evaluator_config=DBEvaluatorConfig(db_filename="test.db"),
            )

    def test_file_paths_validation(self):
        """
        Test that config works with actual file paths via pipeline.

        Inputs:
            - Configuration with real file paths in pipeline

        Outputs:
            - Config is created successfully with file paths accessible
        """
        config = ThompsonSamplingConfig(
            synthesis_pipeline=self.pipeline,
            num_ts_iterations=10,
            num_warmup_trials=3,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=DBEvaluatorConfig(db_filename="test.db"),
        )

        # Verify that the file paths are accessible through the pipeline
        assert self.reagent_file1 in config.reagent_file_list
        assert self.reagent_file2 in config.reagent_file_list
        assert len(config.reagent_file_list) == 2

    def test_all_selection_strategies(self):
        """
        Test that all valid selection strategies can be configured.

        Inputs:
            - Configurations with each valid strategy config

        Outputs:
            - All configs created successfully
        """
        strategy_configs = [
            GreedyConfig(mode="maximize"),
            RouletteWheelConfig(mode="maximize"),
            UCBConfig(mode="maximize"),
            EpsilonGreedyConfig(mode="maximize"),
            BoltzmannConfig(mode="maximize_boltzmann"),
            BayesUCBConfig(mode="maximize"),
        ]

        for strategy_config in strategy_configs:
            config = ThompsonSamplingConfig(
                synthesis_pipeline=self.pipeline,
                num_ts_iterations=10,
                num_warmup_trials=3,
                strategy_config=strategy_config,
                evaluator_config=DBEvaluatorConfig(db_filename="test.db"),
            )
            assert config.strategy_config is strategy_config

    def test_lookup_evaluator_config(self):
        """
        Test that LookupEvaluatorConfig works correctly.

        Inputs:
            - Configuration with LookupEvaluatorConfig

        Outputs:
            - Config created with correct evaluator config
        """
        config = ThompsonSamplingConfig(
            synthesis_pipeline=self.pipeline,
            num_ts_iterations=10,
            num_warmup_trials=3,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=LookupEvaluatorConfig(
                ref_filename="scores.csv",
                score_col="Score",
            ),
        )

        assert isinstance(config.evaluator_config, LookupEvaluatorConfig)
        assert config.evaluator_config.ref_filename == "scores.csv"
        assert config.evaluator_config.score_col == "Score"

    def test_pipeline_properties(self):
        """
        Test that pipeline properties are accessible through config.

        Inputs:
            - Valid config with pipeline

        Outputs:
            - Pipeline properties accessible via config
        """
        config = ThompsonSamplingConfig(
            synthesis_pipeline=self.pipeline,
            num_ts_iterations=10,
            num_warmup_trials=3,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=DBEvaluatorConfig(db_filename="test.db"),
        )

        assert config.num_components == 2
        assert config.num_steps == 1
        assert len(config.reagent_file_list) == 2
