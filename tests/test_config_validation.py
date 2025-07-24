import pytest
import tempfile
import os
from PRISMS.thompson_sampling import StandardSamplerConfig, EnhancedSamplerConfig
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
    
    def test_standard_config_creation(self):
        """
        Test that StandardSamplerConfig can be created with valid data.
        
        Inputs:
            - Valid configuration parameters for standard sampler
        
        Outputs:
            - StandardSamplerConfig instance is created successfully
        """
        config = StandardSamplerConfig(
            sampler_type="standard",
            ts_mode="maximize",
            evaluator_class_name="DBEvaluator",
            evaluator_arg="test_arg",
            reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            num_ts_iterations=10,
            reagent_file_list=[self.reagent_file1, self.reagent_file2],
            num_warmup_trials=3
        )
        
        assert config.sampler_type == "standard"
        assert config.ts_mode == "maximize"
        assert config.evaluator_class_name == "DBEvaluator"
        assert config.num_ts_iterations == 10
        assert len(config.reagent_file_list) == 2
    
    def test_enhanced_config_creation(self):
        """
        Test that EnhancedSamplerConfig can be created with valid data.
        
        Inputs:
            - Valid configuration parameters for enhanced sampler
        
        Outputs:
            - EnhancedSamplerConfig instance is created successfully
        """
        config = EnhancedSamplerConfig(
            sampler_type="enhanced",
            processes=4,
            scaling=1.0,
            percent_of_library=0.5,
            minimum_no_of_compounds_per_core=2,
            stopping_criteria=10,
            evaluator_class_name="DBEvaluator",
            evaluator_arg="test_arg",
            reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            num_ts_iterations=10,
            reagent_file_list=[self.reagent_file1, self.reagent_file2],
            num_warmup_trials=3
        )
        
        assert config.sampler_type == "enhanced"
        assert config.processes == 4
        assert config.percent_of_library == 0.5
        assert config.minimum_no_of_compounds_per_core == 2
    
    def test_config_with_optional_fields(self):
        """
        Test that optional fields work correctly in both config types.
        
        Inputs:
            - Config with optional fields specified
        
        Outputs:
            - Config is created with optional fields set correctly
        """
        # Test standard config with optional fields
        standard_config = StandardSamplerConfig(
            sampler_type="standard",
            ts_mode="maximize",
            evaluator_class_name="DBEvaluator",
            evaluator_arg="test_arg",
            reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            num_ts_iterations=10,
            reagent_file_list=[self.reagent_file1],
            num_warmup_trials=3,
            results_filename="results.csv",
            log_filename="test.log"
        )
        
        assert standard_config.results_filename == "results.csv"
        assert standard_config.log_filename == "test.log"
        
        # Test enhanced config with optional fields
        enhanced_config = EnhancedSamplerConfig(
            sampler_type="enhanced",
            processes=2,
            scaling=1.0,
            percent_of_library=0.5,
            minimum_no_of_compounds_per_core=1,
            stopping_criteria=10,
            evaluator_class_name="DBEvaluator",
            evaluator_arg="test_arg",
            reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            num_ts_iterations=10,
            reagent_file_list=[self.reagent_file1],
            num_warmup_trials=3,
            results_filename="enhanced_results.csv",
            log_filename="enhanced_test.log"
        )
        
        assert enhanced_config.results_filename == "enhanced_results.csv"
        assert enhanced_config.log_filename == "enhanced_test.log"
    
    def test_validation_errors(self):
        """
        Test that validation errors are raised for invalid configurations.
        
        Inputs:
            - Invalid configuration parameters
        
        Outputs:
            - ValidationError is raised for each invalid case
        """
        # Test invalid sampler_type
        with pytest.raises(ValidationError):
            StandardSamplerConfig(
                sampler_type="invalid",  # Invalid
                ts_mode="maximize",
                evaluator_class_name="DBEvaluator",
                evaluator_arg="test_arg",
                reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
                num_ts_iterations=10,
                reagent_file_list=[self.reagent_file1],
                num_warmup_trials=3
            )
        
        # Test invalid ts_mode
        with pytest.raises(ValidationError):
            StandardSamplerConfig(
                sampler_type="standard",
                ts_mode="invalid_mode",  # Invalid
                evaluator_class_name="DBEvaluator",
                evaluator_arg="test_arg",
                reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
                num_ts_iterations=10,
                reagent_file_list=[self.reagent_file1],
                num_warmup_trials=3
            )
        
        # Test invalid percent_of_library for enhanced config
        with pytest.raises(ValidationError):
            EnhancedSamplerConfig(
                sampler_type="enhanced",
                processes=2,
                scaling=1.0,
                percent_of_library=1.5,  # Invalid (> 1)
                minimum_no_of_compounds_per_core=1,
                stopping_criteria=10,
                evaluator_class_name="DBEvaluator",
                evaluator_arg="test_arg",
                reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
                num_ts_iterations=10,
                reagent_file_list=[self.reagent_file1],
                num_warmup_trials=3
            )
        
        # Test negative processes
        with pytest.raises(ValidationError):
            EnhancedSamplerConfig(
                sampler_type="enhanced",
                processes=-1,  # Invalid (negative)
                scaling=1.0,
                percent_of_library=0.5,
                minimum_no_of_compounds_per_core=1,
                stopping_criteria=10,
                evaluator_class_name="DBEvaluator",
                evaluator_arg="test_arg",
                reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
                num_ts_iterations=10,
                reagent_file_list=[self.reagent_file1],
                num_warmup_trials=3
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
            StandardSamplerConfig(
                sampler_type="standard",
                ts_mode="maximize",
                evaluator_class_name="DBEvaluator",
                evaluator_arg="test_arg",
                reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
                num_ts_iterations="not_an_integer",  # Invalid type
                reagent_file_list=[self.reagent_file1],
                num_warmup_trials=3
            )
        
        # Test int instead of float
        with pytest.raises(ValidationError):
            EnhancedSamplerConfig(
                sampler_type="enhanced",
                processes=2,
                scaling="not_a_float",  # Invalid type
                percent_of_library=0.5,
                minimum_no_of_compounds_per_core=1,
                stopping_criteria=10,
                evaluator_class_name="DBEvaluator",
                evaluator_arg="test_arg",
                reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
                num_ts_iterations=10,
                reagent_file_list=[self.reagent_file1],
                num_warmup_trials=3
            )
    
    def test_file_existence_validation(self):
        """
        Test that config works with actual file paths.
        
        Inputs:
            - Configuration with real file paths
        
        Outputs:
            - Config is created successfully with file paths
        """
        config = StandardSamplerConfig(
            sampler_type="standard",
            ts_mode="maximize",
            evaluator_class_name="DBEvaluator",
            evaluator_arg="test_arg",
            reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            num_ts_iterations=10,
            reagent_file_list=[self.reagent_file1, self.reagent_file2],
            num_warmup_trials=3
        )
        
        # Verify that the file paths are correctly stored
        assert self.reagent_file1 in config.reagent_file_list
        assert self.reagent_file2 in config.reagent_file_list
        assert len(config.reagent_file_list) == 2 