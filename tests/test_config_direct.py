import pytest
from TACTICS.thompson_sampling import ThompsonSamplingConfig, RandomBaselineConfig
from pydantic import ValidationError

# Test ThompsonSamplingConfig with valid input
def test_thompson_sampling_config_valid():
    """
    This test checks that a valid ThompsonSamplingConfig can be created without errors.
    Inputs:
        - All required fields with correct types and values.
    Outputs:
        - A ThompsonSamplingConfig instance is created.
    """
    config = ThompsonSamplingConfig(
        evaluator_class_name="DBEvaluator",
        evaluator_arg="some_arg",
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        num_ts_iterations=10,
        reagent_file_list=["file1.smi", "file2.smi"],
        num_warmup_trials=2,
        selection_strategy="greedy",
        mode="maximize"
    )
    assert config.selection_strategy == "greedy"
    assert config.mode == "maximize"

# Test ThompsonSamplingConfig with roulette wheel strategy
def test_thompson_sampling_config_roulette_wheel():
    """
    This test checks that a valid ThompsonSamplingConfig with roulette_wheel strategy can be created.
    Inputs:
        - Valid config with roulette_wheel strategy and parameters.
    Outputs:
        - A ThompsonSamplingConfig instance is created with strategy parameters.
    """
    config = ThompsonSamplingConfig(
        evaluator_class_name="DBEvaluator",
        evaluator_arg="some_arg",
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        num_ts_iterations=10,
        reagent_file_list=["file1.smi", "file2.smi"],
        num_warmup_trials=2,
        selection_strategy="roulette_wheel",
        mode="maximize",
        strategy_params={"alpha": 0.1, "beta": 0.1, "scaling": 1.0},
        processes=2
    )
    assert config.selection_strategy == "roulette_wheel"
    assert config.processes == 2

# Test ThompsonSamplingConfig with missing required field
def test_thompson_sampling_config_missing_field():
    """
    This test checks that missing a required field raises a ValidationError.
    Inputs:
        - Missing 'evaluator_class_name' field.
    Outputs:
        - ValidationError is raised.
    """
    with pytest.raises(ValidationError):
        ThompsonSamplingConfig(
            # evaluator_class_name missing
            evaluator_arg="some_arg",
            reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            num_ts_iterations=10,
            reagent_file_list=["file1.smi", "file2.smi"],
            num_warmup_trials=2,
            selection_strategy="greedy",
            mode="maximize"
        )

# Test ThompsonSamplingConfig with invalid selection_strategy
def test_thompson_sampling_config_invalid_strategy():
    """
    This test checks that an invalid selection_strategy raises a ValidationError.
    Inputs:
        - selection_strategy="invalid_strategy" (invalid)
    Outputs:
        - ValidationError is raised.
    """
    with pytest.raises(ValidationError):
        ThompsonSamplingConfig(
            evaluator_class_name="DBEvaluator",
            evaluator_arg="some_arg",
            reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            num_ts_iterations=10,
            reagent_file_list=["file1.smi", "file2.smi"],
            num_warmup_trials=2,
            selection_strategy="invalid_strategy",  # invalid
            mode="maximize"
        )


# Test config model field validation
def test_config_field_types():
    """
    Test that config models enforce correct field types.
    Inputs:
        - Invalid field types (string instead of int, etc.)
    Outputs:
        - ValidationError is raised for type mismatches
    """
    # Test that string for int field raises error
    with pytest.raises(ValidationError):
        ThompsonSamplingConfig(
            evaluator_class_name="DBEvaluator",
            evaluator_arg="some_arg",
            reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            num_ts_iterations="not_an_integer",  # invalid type
            reagent_file_list=["file1.smi", "file2.smi"],
            num_warmup_trials=2,
            selection_strategy="greedy",
            mode="maximize"
        )

# Test config model with optional fields
def test_config_optional_fields():
    """
    Test that optional fields work correctly.
    Inputs:
        - Config without optional fields
    Outputs:
        - Config is created successfully with default values for optional fields
    """
    config = ThompsonSamplingConfig(
        evaluator_class_name="DBEvaluator",
        evaluator_arg="some_arg",
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        num_ts_iterations=10,
        reagent_file_list=["file1.smi", "file2.smi"],
        num_warmup_trials=2,
        selection_strategy="greedy",
        mode="maximize"
        # results_filename has default value "results.csv"
        # log_filename is optional (None)
    )
    assert config.results_filename == "results.csv"  # default value
    assert config.log_filename is None

# Test RandomBaselineConfig
def test_random_baseline_config_valid():
    """
    Test that RandomBaselineConfig can be created with valid data.
    Inputs:
        - All required fields with correct types and values.
    Outputs:
        - A RandomBaselineConfig instance is created.
    """
    config = RandomBaselineConfig(
        evaluator_class_name="LookupEvaluator",
        evaluator_arg={"ref_filename": "scores.csv"},
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        reagent_file_list=["file1.smi", "file2.smi"],
        num_trials=100,
        num_to_save=10
    )
    assert config.num_trials == 100
    assert config.num_to_save == 10

# Test all valid selection strategies
def test_all_selection_strategies():
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
            reagent_file_list=["file1.smi"],
            num_warmup_trials=3,
            selection_strategy=strategy,
            mode="maximize"
        )
        assert config.selection_strategy == strategy

# Test dict evaluator_arg
def test_dict_evaluator_arg():
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
        reagent_file_list=["file1.smi"],
        num_warmup_trials=3,
        selection_strategy="greedy",
        mode="maximize"
    )
    assert isinstance(config.evaluator_arg, dict)
    assert config.evaluator_arg["ref_filename"] == "scores.csv"
