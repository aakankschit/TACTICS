import pytest
import tempfile
import os
from TACTICS.thompson_sampling import ThompsonSamplingConfig, RandomBaselineConfig
from pydantic import ValidationError

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
            evaluator_class_name="DBEvaluator",
            evaluator_arg="test_arg",
            reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            num_ts_iterations=10,
            reagent_file_list=[self.reagent_file1, self.reagent_file2],
            num_warmup_trials=3,
            selection_strategy="greedy",
            mode="maximize"
        )

        assert config.selection_strategy == "greedy"
        assert config.mode == "maximize"
        assert config.evaluator_class_name == "DBEvaluator"
        assert config.num_ts_iterations == 10
        assert len(config.reagent_file_list) == 2

    def test_config_with_roulette_wheel_strategy(self):
        """
        Test that config works with roulette wheel selection strategy.

        Inputs:
            - Valid configuration with roulette_wheel strategy and parameters

        Outputs:
            - ThompsonSamplingConfig instance with strategy params
        """
        config = ThompsonSamplingConfig(
            evaluator_class_name="DBEvaluator",
            evaluator_arg="test_arg",
            reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            num_ts_iterations=10,
            reagent_file_list=[self.reagent_file1, self.reagent_file2],
            num_warmup_trials=3,
            selection_strategy="roulette_wheel",
            mode="maximize",
            strategy_params={
                "alpha": 0.1,
                "beta": 0.1,
                "scaling": 1.0
            }
        )

        assert config.selection_strategy == "roulette_wheel"
        assert config.strategy_params["alpha"] == 0.1
        assert config.strategy_params["beta"] == 0.1

    def test_config_with_optional_fields(self):
        """
        Test that optional fields work correctly.

        Inputs:
            - Config with optional fields specified

        Outputs:
            - Config is created with optional fields set correctly
        """
        config = ThompsonSamplingConfig(
            evaluator_class_name="DBEvaluator",
            evaluator_arg="test_arg",
            reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            num_ts_iterations=10,
            reagent_file_list=[self.reagent_file1],
            num_warmup_trials=3,
            selection_strategy="greedy",
            mode="maximize",
            results_filename="results.csv",
            log_filename="test.log",
            batch_size=5,
            max_resamples=100,
            processes=4,
            min_cpds_per_core=20
        )

        assert config.results_filename == "results.csv"
        assert config.log_filename == "test.log"
        assert config.batch_size == 5
        assert config.max_resamples == 100
        assert config.processes == 4
        assert config.min_cpds_per_core == 20

    def test_random_baseline_config_creation(self):
        """
        Test that RandomBaselineConfig can be created with valid data.

        Inputs:
            - Valid configuration for random baseline

        Outputs:
            - RandomBaselineConfig instance is created successfully
        """
        config = RandomBaselineConfig(
            evaluator_class_name="LookupEvaluator",
            evaluator_arg={"ref_filename": "scores.csv"},
            reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            reagent_file_list=[self.reagent_file1, self.reagent_file2],
            num_trials=100,
            num_to_save=10,
            ascending_output=False,
            outfile_name="baseline_results.csv"
        )

        assert config.num_trials == 100
        assert config.num_to_save == 10
        assert config.ascending_output == False
        assert config.outfile_name == "baseline_results.csv"

    def test_validation_errors_invalid_strategy(self):
        """
        Test that validation errors are raised for invalid strategy.

        Inputs:
            - Invalid selection_strategy parameter

        Outputs:
            - ValidationError is raised
        """
        with pytest.raises(ValidationError):
            ThompsonSamplingConfig(
                evaluator_class_name="DBEvaluator",
                evaluator_arg="test_arg",
                reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
                num_ts_iterations=10,
                reagent_file_list=[self.reagent_file1],
                num_warmup_trials=3,
                selection_strategy="invalid_strategy",  # Invalid
                mode="maximize"
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
            ThompsonSamplingConfig(
                evaluator_class_name="DBEvaluator",
                evaluator_arg="test_arg",
                reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
                num_ts_iterations=10,
                reagent_file_list=[self.reagent_file1],
                num_warmup_trials=3,
                selection_strategy="greedy",
                mode="invalid_mode"  # Invalid
            )

    def test_type_validation(self):
        """
        Test that type validation works correctly.

        Inputs:
            - Configuration with wrong types

        Outputs:
            - ValidationError is raised for type mismatches
        """
        # Test string instead of int
        with pytest.raises(ValidationError):
            ThompsonSamplingConfig(
                evaluator_class_name="DBEvaluator",
                evaluator_arg="test_arg",
                reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
                num_ts_iterations="not_an_integer",  # Invalid type
                reagent_file_list=[self.reagent_file1],
                num_warmup_trials=3,
                selection_strategy="greedy",
                mode="maximize"
            )

        # Test invalid list type for reagent_file_list
        with pytest.raises(ValidationError):
            ThompsonSamplingConfig(
                evaluator_class_name="DBEvaluator",
                evaluator_arg="test_arg",
                reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
                num_ts_iterations=10,
                reagent_file_list="not_a_list",  # Invalid type
                num_warmup_trials=3,
                selection_strategy="greedy",
                mode="maximize"
            )

    def test_file_paths_validation(self):
        """
        Test that config works with actual file paths.

        Inputs:
            - Configuration with real file paths

        Outputs:
            - Config is created successfully with file paths
        """
        config = ThompsonSamplingConfig(
            evaluator_class_name="DBEvaluator",
            evaluator_arg="test_arg",
            reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            num_ts_iterations=10,
            reagent_file_list=[self.reagent_file1, self.reagent_file2],
            num_warmup_trials=3,
            selection_strategy="greedy",
            mode="maximize"
        )

        # Verify that the file paths are correctly stored
        assert self.reagent_file1 in config.reagent_file_list
        assert self.reagent_file2 in config.reagent_file_list
        assert len(config.reagent_file_list) == 2

    def test_all_selection_strategies(self):
        """
        Test that all valid selection strategies can be configured.

        Inputs:
            - Configurations with each valid strategy

        Outputs:
            - All configs created successfully
        """
        strategies = ["greedy", "roulette_wheel", "ucb", "epsilon_greedy", "boltzmann"]

        for strategy in strategies:
            config = ThompsonSamplingConfig(
                evaluator_class_name="DBEvaluator",
                evaluator_arg="test_arg",
                reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
                num_ts_iterations=10,
                reagent_file_list=[self.reagent_file1],
                num_warmup_trials=3,
                selection_strategy=strategy,
                mode="maximize"
            )
            assert config.selection_strategy == strategy

    def test_dict_evaluator_arg(self):
        """
        Test that evaluator_arg can be a dictionary.

        Inputs:
            - Configuration with dict evaluator_arg

        Outputs:
            - Config created with dict arg
        """
        config = ThompsonSamplingConfig(
            evaluator_class_name="LookupEvaluator",
            evaluator_arg={"ref_filename": "scores.csv", "ref_colname": "Score"},
            reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            num_ts_iterations=10,
            reagent_file_list=[self.reagent_file1],
            num_warmup_trials=3,
            selection_strategy="greedy",
            mode="maximize"
        )

        assert isinstance(config.evaluator_arg, dict)
        assert config.evaluator_arg["ref_filename"] == "scores.csv"
