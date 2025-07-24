import pytest
import sys
import os

# Add the legacy directory to the path so we can import the config module directly
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'PRISMS', 'thompson_sampling', 'legacy'))

# Import the config module directly
from PRISMS.thompson_sampling import StandardSamplerConfig, EnhancedSamplerConfig
from pydantic import ValidationError

# Test StandardSamplerConfig with valid input
def test_standard_sampler_config_valid():
    """
    This test checks that a valid StandardSamplerConfig can be created without errors.
    Inputs:
        - All required fields with correct types and values.
    Outputs:
        - A StandardSamplerConfig instance is created.
    """
    config = StandardSamplerConfig(
        sampler_type="standard",
        ts_mode="maximize",
        evaluator_class_name="DBEvaluator",
        evaluator_arg="some_arg",
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        num_ts_iterations=10,
        reagent_file_list=["file1.smi", "file2.smi"],
        num_warmup_trials=2
    )
    assert config.sampler_type == "standard"
    assert config.ts_mode == "maximize"

# Test EnhancedSamplerConfig with valid input
def test_enhanced_sampler_config_valid():
    """
    This test checks that a valid EnhancedSamplerConfig can be created without errors.
    Inputs:
        - All required fields with correct types and values.
    Outputs:
        - An EnhancedSamplerConfig instance is created.
    """
    config = EnhancedSamplerConfig(
        sampler_type="enhanced",
        processes=2,
        scaling=1.0,
        percent_of_library=0.5,
        minimum_no_of_compounds_per_core=1,
        stopping_criteria=10,
        evaluator_class_name="DBEvaluator",
        evaluator_arg="some_arg",
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        num_ts_iterations=10,
        reagent_file_list=["file1.smi", "file2.smi"],
        num_warmup_trials=2
    )
    assert config.sampler_type == "enhanced"
    assert config.processes == 2

# Test StandardSamplerConfig with missing required field
def test_standard_sampler_config_missing_field():
    """
    This test checks that missing a required field raises a ValidationError.
    Inputs:
        - Missing 'ts_mode' field.
    Outputs:
        - ValidationError is raised.
    """
    with pytest.raises(ValidationError):
        StandardSamplerConfig(
            sampler_type="standard",
            # ts_mode missing
            evaluator_class_name="DBEvaluator",
            evaluator_arg="some_arg",
            reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            num_ts_iterations=10,
            reagent_file_list=["file1.smi", "file2.smi"],
            num_warmup_trials=2
        )

# Test EnhancedSamplerConfig with invalid percent_of_library
def test_enhanced_sampler_config_invalid_percent():
    """
    This test checks that percent_of_library outside (0,1] raises a ValidationError.
    Inputs:
        - percent_of_library=1.5 (invalid)
    Outputs:
        - ValidationError is raised.
    """
    with pytest.raises(ValidationError):
        EnhancedSamplerConfig(
            sampler_type="enhanced",
            processes=2,
            scaling=1.0,
            percent_of_library=1.5,  # invalid
            minimum_no_of_compounds_per_core=1,
            stopping_criteria=10,
            evaluator_class_name="DBEvaluator",
            evaluator_arg="some_arg",
            reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            num_ts_iterations=10,
            reagent_file_list=["file1.smi", "file2.smi"],
            num_warmup_trials=2
        )

# Test EnhancedSamplerConfig with negative processes
def test_enhanced_sampler_config_negative_processes():
    """
    This test checks that a negative value for 'processes' raises a ValidationError.
    Inputs:
        - processes=-1 (invalid)
    Outputs:
        - ValidationError is raised.
    """
    with pytest.raises(ValidationError):
        EnhancedSamplerConfig(
            sampler_type="enhanced",
            processes=-1,  # invalid
            scaling=1.0,
            percent_of_library=0.5,
            minimum_no_of_compounds_per_core=1,
            stopping_criteria=10,
            evaluator_class_name="DBEvaluator",
            evaluator_arg="some_arg",
            reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            num_ts_iterations=10,
            reagent_file_list=["file1.smi", "file2.smi"],
            num_warmup_trials=2
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
        StandardSamplerConfig(
            sampler_type="standard",
            ts_mode="maximize",
            evaluator_class_name="DBEvaluator",
            evaluator_arg="some_arg",
            reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            num_ts_iterations="not_an_integer",  # invalid type
            reagent_file_list=["file1.smi", "file2.smi"],
            num_warmup_trials=2
        )

# Test config model with optional fields
def test_config_optional_fields():
    """
    Test that optional fields work correctly.
    Inputs:
        - Config without optional fields
    Outputs:
        - Config is created successfully with None values for optional fields
    """
    config = StandardSamplerConfig(
        sampler_type="standard",
        ts_mode="maximize",
        evaluator_class_name="DBEvaluator",
        evaluator_arg="some_arg",
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        num_ts_iterations=10,
        reagent_file_list=["file1.smi", "file2.smi"],
        num_warmup_trials=2
        # results_filename and log_filename are optional
    )
    assert config.results_filename is None
    assert config.log_filename is None 